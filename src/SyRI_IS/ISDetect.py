from re import compile
import sys
import os
import subprocess
import pandas as pd
import logging
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sklearn.cluster import AgglomerativeClustering
from utility import format_time
from alignment import alignment
import pysam
import matplotlib.pyplot as plt
from collections import Counter
import pyarrow
import pyarrow.parquet as pq

from pathlib import Path

###############################
# Based on LoRTEv1.2
##############################


def Fastq2Fasta(fasta, fastq):
    os.system("awk '(NR - 1) % 4 < 2' " + fastq + " | sed 's/@/>/' > ", fasta)
    os.system(r'sed -e -i "s/\(^>\S*\) .*/\1/" ' + fasta) # remove discription 

def FastqDir2Fasta(outFile, dir):
    C1, C2, C3 = '', '', ''
    if os.path.isfile(dir):
        tmp_fq = outFile.replace('.fasta', '.fastq')
        C1 = 'filtlong --min_length 1000 --target_bases ' + str(4000000 * 50) + ' ' + dir + ' > ' + tmp_fq
        C2 = 'seqkit fq2fa ' + tmp_fq + ' -o ' + outFile
    else:
        C2 = 'seqkit fq2fa ' + dir + '/* -o ' + outFile
    # remove discription
    C3 = r'sed -i -e "s/\(^>\S*\) .*/\1/" ' + outFile 
    os.system(C1) if C1 else None
    os.system(C2) if C2 else None
    os.system(C3) if C3 else None
    logging.info(C1) if C1 else None
    logging.info(C2) if C2 else None
    logging.info(C3) if C3 else None
    return outFile

def makedb(fastaFile, logName, tmpdir='tmp'):
    """
    Creates a database file: nucl at the file location. Required for alignments.
    """
    start_time = time.time()
    logging.info("Making database (" + logName + ")")
    db_prefix = os.path.join(tmpdir, 'blastdb', os.path.basename(fastaFile))
    try:
        command_file = (" -in " + fastaFile + " -dbtype nucl")
        command_file += " -out " + db_prefix
        os.system("makeblastdb" + command_file)
    except Exception as e:
        print(e)
        print("Database construction failed (" +
              logName + "). Check input reads, exiting...")
        sys.exit(1)
    proc_time = time.time() - start_time
    logging.info("Constructed database (" + logName + ") " +
                 format_time(proc_time))
    return db_prefix


def BLAST(queryFasta, targetDB, eValue, prefix, threads=1, dbdir = ''):
    """
    Used to detect query in targetDB. Use makedb to create database.
    """
    C1 = "blastn -task blastn -query "+queryFasta
    C2 = " -evalue "+str(eValue)+ " -num_threads " + str(threads) + \
        " -outfmt \"7 sseqid qseqid qstart qend sstart send length evalue\" -out " + \
        prefix+"-blastn"
    C3 = " -db "+dbdir if dbdir else ''
    blast = C1+C2+C3
    logging.info(blast)
    os.system(blast)
    return(prefix+"-blastn")

def genDFOfHowReadsMatchIS_(BlastFile, ISlength = 3000, bp_from_end = 300, min_match = 1000):
    """
    Split the reads by IS matches, and get the flank positions.
    """
    # parse reads and get the flanking matches in fasta
    dat_ = pd.read_table(BlastFile, header=None, comment='#')
    dat_ = dat_.set_axis(['subject_id', 'query_id', 'q_start', 'q_end', 's_start', 's_end', 'alignment_length', 'evalue'],
                         axis=1)

    # take valid reads
    if ISlength < 1000:
        logging.error("ISlength = " + str(ISlength) + " bp. ISs must be larger than 1000 bp.")
    originalN = dat_.shape[0]
    valid_reads = 'q_start < ' + str(bp_from_end) + ' & q_end > ' + str(ISlength - bp_from_end)
    valid_reads += ' & alignment_length > ' + str(min_match)
    dat_ = dat_.query(valid_reads).eval('positive_strand = s_end - s_start > 0')
    filteredN = dat_.shape[0]
    print(str(filteredN) + " reads out of " +
          str(originalN) + " reads passed QC.")

    # get end positions the order is reversed for negative strand
    dat_ = dat_.assign(i1=1, i2=lambda df: df.s_start,
                       i3=lambda df: df.s_end, i4=-1, matchid=0)
    dat_ = dat_.assign(i2=lambda df: df.i2.mask(~df.positive_strand, df.s_end),
                       i3=lambda df: df.i3.mask(~df.positive_strand, df.s_start))

    dat_ = dat_.sort_values('s_start')
    groups_ = dat_.groupby('subject_id')
    # id = dat_.columns.get_loc("i1")

    # split multiple hit reads
    for name, group in groups_:
        if len(group) > 1:
            ids = group.index
            for i, _ in enumerate(group.iterrows()):
                r, _1 = _
                dat_.loc[ids[i], 'matchid'] = i  # matchid
                if i > 0: # if there are other IS hits in the same read save the end of previous is 
                    dat_.loc[ids[i], 'i1'] = group.i3.iloc[i-1]
                if i < len(group)-1: # if this is not the last IS hit in the same read save the start of next is
                    dat_.loc[ids[i], 'i4'] = group.i2.iloc[i+1]
    dat_ = dat_.sort_values(['s_start', 'subject_id'])

    # add column of flank ids
    # f: left on read; r: righ on read
    dat_ = dat_.\
        assign(f=lambda df: df.subject_id + df.matchid.map('_{0:0>2}f'.format)).\
        assign(r=lambda df: df.subject_id + df.matchid.map('_{0:0>2}r'.format)).\
        melt(id_vars=list(dat_.columns), var_name='flank', value_name='flank_id')

    # add column of flank directions.
    # This is to ensure that the matches of flanks to references start from zero,
    # where zero is closer to the IS match.
    dat_ = dat_.\
        assign(inv=lambda df: ((df.flank.map(lambda x: x == 'f') & df.positive_strand) |
               (df.flank.map(lambda x: x == 'r') & ~ df.positive_strand))
               )

    return dat_


def GetReadNameSeqDictOfMatchedReadsWithoutRedundancy(dat_, readFile, tempFolder, prefix):
    """
        Polish the matchlist by removing redundant reads.
        use
    """
    # write match list (list of read names that matched IS) for seqtk etc
    matchList = os.path.join(tempFolder, prefix + ".ISmatch.lst")

    os.makedirs(tempFolder, exist_ok=True)
    with open(matchList, 'w') as f:
        for item in dat_.subject_id:
            f.write("%s\n" % item)

    matchFile = os.path.join(tempFolder, prefix+".ISmatch.fasta")

    try:
        # get and save a list of unique read names
        with open(matchList + "_uniq", "w") as outfile:
            subprocess.run(
                'cat ' + matchList + ' | sort | uniq',
                stdout=outfile,
                shell=True
            )

        # get and save the sequences of the unique reads obtained above
        with open(matchFile, "w") as output:
            subprocess.run(
                'seqtk subseq ' + readFile + ' ' + matchList+"_uniq" +
                ' | seqkit rmdup',
                stdout=output,
                shell=True
            )
    except Exception as e:
        logging.error("Read alignment failed, check input reads, exiting...")
        logging.error(str(e))

    seqdict_ = SeqIO.to_dict(SeqIO.parse(matchFile, "fasta"))

    return seqdict_


def WriteFastaOfSeqFlankingIS_(ISMatchedReadsDF, readFile, tempFolder, prefix):
    """
        Write FlankFile (.fasta) with the flanking sequences.
    """
    # key: queryid, value: sequence of the matched reads
    # retrives the sequences of the matched reads
    seqdict_ = GetReadNameSeqDictOfMatchedReadsWithoutRedundancy(ISMatchedReadsDF, readFile, tempFolder, prefix)

    minSeqLen = 100  # minimum match with reference to accept must be large enough to avoid false positives such as rrnB
    maxSeqLen = 1000  # maximum match with reference to accept to shorten runs
    flankingSequenceList = []

    for i, row in enumerate(ISMatchedReadsDF.iterrows()):
        id, item = row
        if item.inv:
            if item.i2 > item.i1 + minSeqLen:
                flank_length = item.i2 - item.i1 - 1
                flank_length = flank_length if flank_length < maxSeqLen else maxSeqLen
                record = SeqRecord(
                    Seq(seqdict_[
                        item.subject_id].seq[(item.i2-2-flank_length+1):(item.i2-2)
                                            ].reverse_complement()),
                    id=item.flank_id,
                    description=""
                )
                flankingSequenceList.append(record)
        elif ~item.inv:
            i4 = item.i4 if item.i4 > 0 else len(seqdict_[item.subject_id])
            if i4 > item.i3 + minSeqLen:
                flank_length = i4 - item.i3 -1
                flank_length = flank_length if flank_length < maxSeqLen else maxSeqLen
                record = SeqRecord(
                    Seq(seqdict_[item.subject_id].seq[(item.i3):(item.i3+flank_length-1)]),
                    id=item.flank_id,
                    description=""
                )
                flankingSequenceList.append(record)

    logging.info('%d out of %d potential flanks passed QC.' %
          (len(flankingSequenceList), ISMatchedReadsDF.shape[0]))

    FlankFile = os.path.join(tempFolder, prefix+".flank.fasta")
    SeqIO.write(flankingSequenceList, FlankFile, "fasta")

    return FlankFile


def ParseFlankGenomeMappingBLAST_(BLASTFile):
    flankBlast = pd.read_table(BLASTFile, header=None, comment='#')
    flankBlast = flankBlast.set_axis(['subject_id', 'query_id', 'q_start', 'q_end', 's_start', 's_end', 'alignment_length', 'evalue'],
                                     axis=1)
    # get high quality flanks that match well with the reference from its end
    flankBlast = flankBlast.query('q_start < 10 & q_end > 100')
    flankBlast = flankBlast.loc[flankBlast.groupby(
        'query_id').evalue.idxmin()].reset_index(drop=True)

    goodFlankN = flankBlast['query_id'].nunique()
    allFlankN = flankBlast['query_id'].nunique() 
    logging.info('%d out of %d flanks aligned reference sequence well.'
          % (allFlankN, goodFlankN)
          )
    return flankBlast 


def GetClusters(readFlankInsPosDF, key = 's_start_y', distance_threshold=100, minSupportNum=3, ignoreStrand=False):
    """
        returns the insertion positions identified by agglomerative clustering,
        and the assignment of the flank of IS in reads to each cluster.
    """
    if not ignoreStrand:
        assert 'IS_strand' in readFlankInsPosDF.columns, 'IS_strand not in readFlankInsPosDF.columns'
    # Cluster by IS seq genome seq border
    clustering = AgglomerativeClustering(
        n_clusters=None,
        linkage='single',  # closest within groups
        distance_threshold=distance_threshold
    ).fit([[i] for i in readFlankInsPosDF[key]])

    # remove insert_id, cluster_id if exist
    if 'insert_id' in readFlankInsPosDF.columns:
        readFlankInsPosDF = readFlankInsPosDF.drop('insert_id', axis=1)
    if 'cluster_id' in readFlankInsPosDF.columns:
        readFlankInsPosDF = readFlankInsPosDF.drop('cluster_id', axis=1)
    clustersRaw_ = pd.concat([readFlankInsPosDF, pd.Series(clustering.labels_, name='cluster_id')], axis=1)

	# get the median position of each cluster
    groups_ = ['cluster_id'] if ignoreStrand else ['cluster_id', 'IS_strand']
    clusters = clustersRaw_.\
        groupby(groups_).\
        agg(count=(key, 'count'),
            medianpos=(key, 'median')
            ).\
        query('count >= '+str(minSupportNum)).\
        reset_index()

    clusters = clusters.sort_values('medianpos').reset_index()
    clusterN_ = clusters.shape[0]
    clusters['insert_id'] = range(1, 1 + clusterN_)
    clustersRaw = pd.merge(clustersRaw_,
        clusters[['cluster_id', 'insert_id']],
        on = 'cluster_id',
        how = 'left'
    )
    clustersRaw.loc[clustersRaw.insert_id.isna(), 'insert_id'] = -1
    clustersRaw.insert_id = clustersRaw.insert_id.astype(int)

    # Check whether the flank is in 5' from 0 in reference genome
    if 'flank' in clustersRaw.columns:
        clustersRaw = clustersRaw.assign(flankMatchDir=lambda df: ~df.flank.map(lambda x: x=='f') ^ df.IS_strand)

    clusters = clusters.drop(['index', 'cluster_id'], axis=1, errors='ignore')
    return clusters, clustersRaw

def getFastaLength(fastaFile):
    """
        Returns the length of the first sequence in fastaFile
    """
    for seq in SeqIO.parse(fastaFile, "fasta"):
        return len(seq.seq)


"""
q_start_x: 1-based index of the start of the alignment in the IS
s_start_x: 1-based index of the start of the alignment in the fastq read
positive_strand: whether the IS direction is the same as the read direction
i1-4: 1-based index of the start and end of the flanking sequence in the read. -1 means the flank is at the right end of the read.
flank: f or r, indicating whether the flank is upstream or downstream of the IS
inv: whether the flank sequence is inverted before mathching to the reference genome.
For instance, if the flank is upstream of the IS and the read is in the same direction as the IS, inv is True.
This inversion makes all the starting query be the position of IS referece genome border.
"""
def GetInsertions(readFile, refgenome, ISfile, prefix,
                  threads=4, tempFolder="tmp", distance_threashold=20, minSupportNum=3,
                  flank_file =''):
    """
        Get IS insertion position based on reads in `readFile` that have blast hits of
        ISs in `ISfile`, and matching those reads to `refgenome.`
    """
    os.makedirs(tempFolder, exist_ok=True)

    # Blast search of IS within reads
    FlankingSeqPosDF_ = None
    FlankFile =''
    if flank_file == '':
        logName = "nanopore reads"
        dbdir = os.path.join(tempFolder, 'blastdb', os.path.basename(readFile))
        makedb(readFile, logName, tmpdir=tempFolder)
        BLASTprefix_ = os.path.join(tempFolder, prefix+'-ISinReads')
        ISinReadsBLAST = BLAST(ISfile, readFile, 1e-40,
                            BLASTprefix_, threads=threads, dbdir = dbdir)
        logging.info("Ran Blast search of IS in reads.")

        # Get the position of sequences flanking ISs
        fasta_length = getFastaLength(readFile)
        FlankingSeqPosDF_ = genDFOfHowReadsMatchIS_(ISinReadsBLAST, ISlength=fasta_length, bp_from_end = fasta_length, min_match = 300)

        # Write FlankFile (.fasta) with the flanking sequences.
        FlankFile = WriteFastaOfSeqFlankingIS_(FlankingSeqPosDF_, readFile, tempFolder, prefix)
    else:
        FlankFile = flank_file+'.flank.fasta'
        FlankingSeqPosDF_ =  pd.read_csv(flank_file+'.csv')
        assert all([x in FlankingSeqPosDF_.columns for x in ['flank_id', 'flank']]), 'FlankingSeqPosDF_ must have columns: flank_id, flank'
        FlankingSeqPosDF_ = FlankingSeqPosDF_.assign(flank_id = lambda x: f'{x.pos_id}_{x.flank} flank')
        assert os.path.isfile(FlankFile), 'FlankFile does not exist: ' + FlankFile

    # Find the flanking sequences within reference genome
    logName = "reference fasta"
    dbdir = os.path.join(tempFolder, 'blastdb', os.path.basename(refgenome))
    makedb(refgenome, logName, tmpdir=tempFolder)
    BLASTprefix_ = os.path.join(tempFolder, prefix+'-FlkinRef')
    FlkinRefBLAST = BLAST(FlankFile, refgenome, 1e-40, BLASTprefix_, threads=threads, dbdir=dbdir)
    logging.info("Ran Blast search of flanking sequences in references.")

    # Get the matched positions of sequence surrounding ISs within the reference genomes
    # join the BLAST result matching flanking sequences to the genome with readFlankDF having where the flank seqs originates from.
    flankBlast = ParseFlankGenomeMappingBLAST_(FlkinRefBLAST)
    # Detect the strand of IS corresponding to the flank.
    FlankingSeqPosDF_ = pd.merge(FlankingSeqPosDF_, flankBlast, left_on='flank_id', right_on='query_id').\
        assign(IS_strand=lambda df: (df.s_end_y > df.s_start_y)
               ^ df.flank.map(lambda x: x == 'f'))

    # Estimate IS insertion position
    clusters, rawClusterDF_ = GetClusters(FlankingSeqPosDF_, distance_threshold=distance_threashold, minSupportNum=minSupportNum)
    rawClusterDF_.to_csv(os.path.join(tempFolder, prefix +'-FlkPosDF.csv'), index = False)
    logging.info('GetInsertions: Detected %d IS insertions.' % (clusters.shape[0]))

    FlankingSeqPosDF_.to_csv(os.path.join(tempFolder, prefix+'-ISFlankingSeqPos.csv'), index = False)

    return clusters, (FlankingSeqPosDF_, rawClusterDF_)

def parse_flank_is_pair_genome_mapping_blast(BLASTFile, flank_size = 100):
    flankBlast = pd.read_table(BLASTFile, header=None, comment='#')
    flankBlast = flankBlast.set_axis(['subject_id', 'query_id', 'q_start', 'q_end', 's_start', 's_end', 'alignment_length', 'evalue'],
                                     axis=1)
    # get high quality flanks that match well with the reference from its end
    flankBlast = flankBlast.query('q_start < 10 & q_end > ' + str(flank_size - 10))
    flankBlast = flankBlast.loc[flankBlast.groupby(
        'query_id').evalue.idxmin()].reset_index(drop=True)

    goodFlankN = flankBlast['query_id'].nunique()
    allFlankN = flankBlast['query_id'].nunique() 
    logging.info('%d out of %d flanks aligned reference sequence well.'
          % (allFlankN, goodFlankN)
          )
    return flankBlast 

# THe flanks must be: flnak (over 100 bp) - IS (over 1000 bp)  format
# FlankingSeqPosDF_: flank_id, s_start, s_end
def get_insertion_position(FlankingSeqPosDF_, FlankFile, refgenome, prefix, tempFolder,
                           distance_threshold=20, minSupportNum=1, flank_size = 100, savesubdir =''):
    assert all([x in FlankingSeqPosDF_.columns for x in ['flank_id', 's_start', 's_end']]), 'FlankingSeqPosDF_ must have columns: flank_id, s_start, s_end'
    assert all([x not in FlankingSeqPosDF_.columns for x in ['q_start', 'q_end']]), 'FlankingSeqPosDF_ must have columns: flank_id, s_start, s_end'

    BLASTprefix_ = os.path.join(tempFolder, prefix+'-FlkinRef')
    FlkinRefBLAST = BLAST(FlankFile, refgenome, 1e-40, BLASTprefix_, dbdir=os.path.join(tempFolder, 'blastdb', os.path.basename(refgenome)))
    logging.info("Ran Blast search to find flanking sequences in references.")

    # Get the matched positions of sequence surrounding ISs within the reference genomes
    # join the BLAST result matching flanking sequences to the genome with readFlankDF having where the flank seqs originates from.
    flankBlast = parse_flank_is_pair_genome_mapping_blast(FlkinRefBLAST, flank_size=flank_size)
    # Detect the strand of IS corresponding to the flank.
    FlankingSeqPosDF_ = pd.merge(FlankingSeqPosDF_, flankBlast, left_on='flank_id', right_on='query_id').\
        assign(IS_strand=lambda df: (df.s_end_y > df.s_start_y)
               ^ df.flank.map(lambda x: x == 'r'))

    # Estimate IS insertion position
    FlankingSeqPosDF_['is_match_length'] = FlankingSeqPosDF_.q_end - flank_size
    # IS match starts after the flank
    #FlankingSeqPosDF_['is_ins_pos'] = FlankingSeqPosDF_.apply(lambda x: x.s_start_y + (1 if x.s_start_y < x.s_end_y else -1) * (x.is_match_length - x.q_start + 1), axis = 1)
    FlankingSeqPosDF_['is_ins_pos'] = FlankingSeqPosDF_.apply(lambda x: x.s_start_y + (1 if x.s_start_y < x.s_end_y else -1) * (flank_size - x.q_start + 1), axis = 1)
    clusters, rawClusterDF_ = GetClusters(FlankingSeqPosDF_, key ='is_ins_pos',  distance_threshold=distance_threshold, minSupportNum=minSupportNum)
    rawClusterDF_.to_csv(os.path.join(tempFolder, savesubdir, prefix +'-FlkPosDF.csv'), index = False)
    logging.info('GetInsertions: Detected %d IS insertions.' % (clusters.shape[0]))

    FlankingSeqPosDF_.to_csv(os.path.join(tempFolder, savesubdir, prefix+'-ISFlankingSeqPos.csv'), index = False)

    return clusters, (FlankingSeqPosDF_, rawClusterDF_)

def MakeGeneList(genbankFile, savedir='no-name', save=False):
    """
        Returns a dataframe of genes with their position and strand
    """
    records = SeqIO.read(genbankFile, 'genbank')
    DF = pd.DataFrame(columns=['gene', 'synonyms', 'bnumber', 'strand', 'start', 'end'])
    bregex = compile(r'b\d{4}')
    for i, f in enumerate(records.features):
        if f.type == 'gene':
            start_ =  int(f.location.start)
            end_ = int(f.location.end)
            strand_ = start_ < end_
            gene_ = f.qualifiers['gene'][0]
            syns_ = f.qualifiers['gene_synonym'][0]
            bnum_ = bregex.match(f.qualifiers['gene_synonym'][0]).group(0)
            DF.loc[i] = [gene_, syns_, bnum_, strand_, start_, end_]
    if save:
        DF.to_csv(savedir, index=False)
        print('Saved new gene list in ' + savedir)
    return DF

def IdentifyGeneInInsertedPos(insertDF, geneList, savedir='no-name', save=False, originalIS = []):
    s_min, e_max = geneList.iloc[0].start, geneList.iloc[-1].end
    firstGene, lastGene = geneList.gene.values[[0,-1]] 
    insertedGeneList = []
    for i, c in insertDF.iterrows():
        if (c.medianpos < s_min) | (c.medianpos > e_max):
            insertedGeneList.append(lastGene+'-'+firstGene)
        else :
            leftGene = geneList[geneList.start < c.medianpos].gene.iloc[-1]
            rightGene = geneList[geneList.end > c.medianpos].gene.iloc[0]
            if leftGene == rightGene :
                insertedGeneList.append(leftGene)
            else :
                insertedGeneList.append(leftGene + '-' + rightGene)
    insertDF_ = pd.concat([insertDF, pd.Series(insertedGeneList, name = 'gene')], axis=1)
    insertDF_ = pd.concat([insertDF_, originalIS], join = 'inner', axis = 0)
    if save:
        insertDF_.to_csv(savedir, index=False)
        print('Saved list of psesudogenized genes in' + savedir)
    return insertDF_

# Create a new fasta file with reads duplicated n times.
# This is to pass the requirement for Clustering algorithm that there should be more than 3 or so supports
def MakeCopyOfFastaFile(newFastaFile, templateFastaFile, n = 5):
    with open(newFastaFile, 'w') as f:
        for seq in SeqIO.parse(templateFastaFile, "fasta"):
            for i in range(n):
                # add suffix to read names so that there are no duplicated names
                seq.id = seq.id + '_' + str(i)
                SeqIO.write(seq, f, "fasta")

def GetDepth_(bamFile, prefix='', tmpFolder='', readFile = '', refgenome ='',
              skipalignment = True, threads=4, aligner = 'minimap2', readtype = 'ont'):
    if not skipalignment:
        if os.path.isfile(readFile) & os.path.isfile(refgenome):
            # example minimap2 --cs --MD -Y -L -ax map-ont ../../test/refs/MDS42.fasta ../fasta/MDS42_IS1.fa
            C = 'minimap2 --cs --MD -Y -L -ax map-ont  -H -f 100 -r1k,10k --rmq=no -t ' + str(threads)
            C += ' ' + refgenome + ' ' + readFile 
            C += ' | samtools sort -m 200M -O BAM - > ' + bamFile
            C += ' && samtools index ' + bamFile
            logging.info(C)
            os.system(C)
        else :
            raise Exception('nonexisting input was set when skipalignment is set False.')

    #import bam file and parse to pd, output all positions (-a) and ignore reads shorter than 1000bp(-l)
    pysamOut = pysam.depth('-a','-l', '1000', bamFile, split_lines=True)
    # parse pysam output
    depthDF = pd.DataFrame([x.split('\t') for x in pysamOut])
    depthDF = depthDF.drop(depthDF.columns[[0]], axis=1).set_axis(['locus', 'depth'], axis =1)
    depthDF = depthDF.astype(int)
    return depthDF

# depthDF: dataframe of depth of every base pair
# insertPosDF: dataframe of insertions
def SplitIntervalsOfDepth_(depthDF, insertPosDF):
    clusterLst_ = [-1] * depthDF.locus.iloc[-1] # initialize the cluster each position belongs to with -1
    prevPos_ = 0
    insertPosDF_ = insertPosDF.copy()
    insertPosDF_ = insertPosDF_.sort_values('medianpos')
    for i, c in insertPosDF_.iterrows():
        pos = int(c.medianpos) 
        clusterLst_[(prevPos_):(pos)] = [c.insert_id] * len(clusterLst_[(prevPos_):(pos)]) 
        if c.insert_id == insertPosDF_.insert_id.iloc[-1]: # if the interval passes the last position of the genome, assign the last interval to the first interval
            clusterLst_[pos:] = [insertPosDF_.insert_id[0]] * len(clusterLst_[pos:])
        prevPos_ = pos
    depthDF = pd.concat([depthDF, pd.Series(clusterLst_, name = 'interval_id')], axis=1)
    return depthDF

def CalcDepthInterval(insertPosDF, bamDir, prefix, tmpFolder, reads, refs, skipalignment = False, threads=4, refgenomeLength = 4.6e6):
    """
        Returns `clusterDepth` with depth of each interval between ISs
        in `insertPosDF` based on aligment by NGMLR
        saved in `bamDir` (can be skipped with option) and
        `depthDF`, with depth and interval_id of every base pair.
    """
    # get the depth by mapping using aligner
    depthDF = GetDepth_(bamDir, prefix=prefix, tmpFolder=tmpFolder, readFile = reads,
                        refgenome =refs, skipalignment = skipalignment, threads=threads,
                        aligner = 'minimap2', readtype = 'ont')
    depthDF = SplitIntervalsOfDepth_(depthDF, insertPosDF)
    depthTable = pyarrow.Table.from_pandas(depthDF)
    pq.write_table(depthTable, os.path.join(tmpFolder, prefix+'-depthOfEveryBase.parquet'),
                                compression='zstd')
    #depthDF.to_csv(os.path.join(tmpFolder, prefix+'-depthOfEveryLocus.csv'), index=False)

	# get the mean depth of each interval
    intervalDepth = depthDF.groupby('interval_id').\
        agg(count = ('depth', 'count'),
        depth = ('depth', 'mean')).\
            reset_index()

    intervalDepth = getStartEndPositionsOfIntervals(insertPosDF, intervalDepth, refgenomeLength = refgenomeLength)
    return intervalDepth, depthDF

def PlotInsertions(insertPosDF, savePrefix = '', genomeLength = 4.6e6):
    logging.getLogger('matplotlib.font_manager').disabled = True # remove font search error from log.
    plt.figure(figsize=(10, 7))
    plt.scatter(insertPosDF.insert_id[insertPosDF.IS_strand],
        insertPosDF.medianpos[insertPosDF.IS_strand],
        label='True Position', s=30, marker = "^", c= 'black')
    plt.scatter(insertPosDF.insert_id[~insertPosDF.IS_strand],
        insertPosDF.medianpos[~insertPosDF.IS_strand],
        label='True Position', s=30, marker = "v", c='grey')
    plt.ylim(bottom=0)
    plt.ylim(top=genomeLength)
    if savePrefix != '':
        plt.savefig(savePrefix+'-insertPos.jpg')

def AssignIntervalID2GeneList(geneDF, depthOfEveryBaseDF):
    geneDF_ = geneDF.copy()
    geneDF_['interval_id'] = -1

    for i, g in geneDF_.iterrows():
        start, end =  depthOfEveryBaseDF.interval_id[g.start-1], depthOfEveryBaseDF.interval_id[g.end-1]
        if start == end :
            geneDF_.loc[i, 'interval_id'] = start

    return geneDF_ 

# return pandas with start and end positions and length of intervals
def getStartEndPositionsOfIntervals(insertionPosDF, intervalDF, refgenomeLength = 4.6e6):
    insertPosDF_ = insertionPosDF.copy()
    insertPosDF_ = insertPosDF_.sort_values('medianpos')
    insertPosDF_ = insertPosDF_.reset_index(drop=True)
    insertPosDF_['interval_id'] = intervalDF.interval_id
    # to int
    insertPosDF_.medianpos = insertPosDF_.medianpos.astype(int)
    insertPosDF_['start'] = insertPosDF_.medianpos.shift(1)
    insertPosDF_['end'] = insertPosDF_.medianpos
    insertPosDF_.loc[0, 'start'] = insertPosDF_.medianpos.iloc[-1]
    insertPosDF_.start = insertPosDF_.start.astype(int) # some why it turns into float
    insertPosDF_['length'] = insertPosDF_.end - insertPosDF_.start
    insertPosDF_.loc[insertPosDF_.length < 0, 'length'] = insertPosDF_.length + refgenomeLength
    insertPosDF_ = pd.merge(insertPosDF_.loc[:,['interval_id', 'start', 'end', 'length']],
                                 intervalDF,
                                 on = 'interval_id',
                                 how = 'left')
    return insertPosDF_ 

def IdentifyCopyNumberChangedGenes(insertPosDepthDF, depthOfEveryBaseDF, geneDF,
                                   threashold = 0.3, save = False, savePrefix=''):
    geneDF_ = AssignIntervalID2GeneList(geneDF, depthOfEveryBaseDF)
    meanDepth = depthOfEveryBaseDF.depth.mean()

	# defined inside to use inside 'apply'
    def CalcDepth(x):
        normalizedDepth = x/meanDepth
        if (normalizedDepth>threashold) & (normalizedDepth < 2-threashold):
            return 1
        else:
            return round(normalizedDepth)

	# round the depth of each interval
    insertPosDepthDF['copies'] = insertPosDepthDF.depth.apply(CalcDepth)
    CopyNumberChangedInterval = insertPosDepthDF.interval_id[insertPosDepthDF.copies != 1]

	# get the genes that have changed copy number
    geneDF_ =  geneDF_[geneDF_.interval_id.isin(CopyNumberChangedInterval)]
    changedGenes = pd.merge(
        geneDF_.loc[:,['gene','bnumber', 'start', 'end', 'interval_id']],
        insertPosDepthDF.loc[:,['interval_id','copies']],
        on = 'interval_id',
        how = 'left'
        )
    
    #save
    changedGenes.to_csv(savePrefix+'-foldChangeGenes.csv', index = False)
    insertPosDepthDF.to_csv(savePrefix+'-foldChangeIntervals.csv', index = False)
    return changedGenes, insertPosDepthDF

# based on the alignment data of each flanking seq, get the f r pair supporting each insertion
def GetFlankingPairs(rawClusterDF, savePrefix ='', save=False):
    groups_ = rawClusterDF.groupby(['subject_id_x', 'matchid']) # group by reads supporting specific IS

    # for every read, if the f r flanks are both present, add the position data of the pair
    pairs = []
    for name, group in groups_:
        if len(group.i1) == 2:
            id1_, id2_ = group.insert_id
            d1_, d2_ = group.flankMatchDir
            pairs.append({(id1_,d1_), (id2_,d2_)})
    
    count_ = Counter([tuple(u) for u in pairs])

    FlankingPairDF = pd.concat([
        pd.DataFrame.from_records(
            [[*(c[0]), *(c[1])] for  c in list(count_.keys())],
            columns = ['insert1', 'from5prime1', 'insert2', 'from5prime2']
            ),
            pd.Series(list(count_.values()), name = 'count'),
        ],
        axis = 1
    )
    
    if save:
        FlankingPairDF.to_csv(savePrefix+'-flankPairs.csv', index=False)

    return FlankingPairDF


def COGDFofDeleteDuplicatedGenes(cogDF, deletedDF, copyNumberChangeDF, savePrefix = ''):
    interGenicInsertions = deletedDF.gene.apply(lambda x: '-' in x)
    outDF = deletedDF.loc[~ interGenicInsertions, ['gene']]
    outDF['copies'] = 0
    outDF = pd.concat([
        copyNumberChangeDF.loc[:,['gene', 'copies']],
        outDF
    ], axis = 0)
    outDF = outDF.drop_duplicates()
    outDF = outDF.merge(
        cogDF,
        left_on = 'gene',
        right_on = 'gene_x',
    )
    if savePrefix !='':
        outDF.to_csv(savePrefix+'-COGDeleteDuplicateGenes.csv', index = False)
    return outDF
