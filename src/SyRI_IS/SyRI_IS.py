import logging
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import time
import sys
import os
import glob

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from scipy.cluster.hierarchy import fclusterdata


from plotsr_ import plotsr_main

def format_time(time):
    d = datetime(1, 1, 1) + timedelta(seconds=time)
    if d.hour == 0 and d.minute == 0:
        return "%d seconds" % (d.second)
    elif d.hour == 0 and d.minute != 0:
        return "%d minutes %d seconds" % (d.minute, d.second)
    else:
        return "%d hours %d minutes %d seconds" % (d.hour, d.minute, d.second)

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

def get_fasta_length(fastafile, rmN = False):
    record = list(SeqIO.parse(fastafile, 'fasta'))
    if len(record) > 1:
        logging.warning('more than one sequence found in '+fastafile)
        print('more than one sequence found in '+fastafile)
    if rmN:
        return len(record[0].seq.replace('N', ''))
    else:
        return len(record[0].seq)

def get_originPos_and_orientation_of_subjectgenome_(subjectgenome, refgenome, tmpdir, queryLength = 300, skip = 1000, attempts = 10, threads = 1):
	# blastn of subjectgenome to the begining of the refgenome
    # find the first hit to the query, which is the first queryLength bases of refgenome
    # if no good hit, move the queryLength bases of refgenome by skip bases and try again for attempts times
    # if still no good hit, return 0
    new_origin = 0
    orientation = 'forward'
    search_start = 0
    tmp_fasta_base = os.path.join(tmpdir, 'refgenomeOrigin.fasta')
    refgenome_seq = SeqIO.read(refgenome, 'fasta')
    db_dir = makedb(subjectgenome, 'subjectgenome', tmpdir = tmpdir)
    logging.info('searching for origin of '+subjectgenome+' in '+refgenome)
    for i in range(attempts):
        tmp_fasta = tmp_fasta_base + '.' + str(i)
        SeqIO.write(refgenome_seq[search_start:search_start+queryLength], tmp_fasta, 'fasta')
        blastn_prefix = os.path.join(tmpdir, 'subjectgenomeOriginSearch-'+str(i))
        blastn = BLAST(tmp_fasta, subjectgenome, 1e-10, blastn_prefix, threads = threads, dbdir = db_dir)
        blastn_df = pd.read_table(blastn, comment='#', header=None)
        blastn_df = blastn_df.set_axis(['subject_id', 'query_id', 'q_start', 'q_end', 's_start', 's_end', 'alignment_length', 'evalue'], axis=1)
        if blastn_df.shape[0] >= 1:
            new_origin = blastn_df.iloc[0]['s_start']
            logging.info('found origin of '+subjectgenome+' in '+refgenome+' at '+str(new_origin))
            orientation = 'forward' if blastn_df.iloc[0]['s_start'] < blastn_df.iloc[0]['s_end'] else 'reverse'
            return new_origin, orientation
        else:
            search_start += skip
    logging.warning('failed to find origin of '+subjectgenome+' in '+refgenome + ' after '+str(attempts)+' attempts')
    return new_origin, orientation

## create a new genome file with the origin resetted
def reset_origin_of_genome_(genome, pos, filename, tmpdir = 'tmp', reverse = False):
    seqkit = ''
    new_genome = os.path.join(tmpdir, filename)
    if reverse:
        if pos == 0: # seqkit errors if position value is 0
            seqkit = 'seqkit seq -r -p '+genome+' > '
        else:   
            seqkit = 'seqkit restart -i '+str(pos) + ' ' + genome + ' | seqkit seq -r -p > ' 
    else:
        if pos == 0:
            seqkit = 'cat '+genome+' > '
        else:
            seqkit = 'seqkit restart -i '+str(pos) + ' ' + genome + ' > ' 
    seqkit += new_genome
    logging.info(seqkit)
    os.system(seqkit)
    return new_genome

def reset_origin_of_subject_genome(subjectgenome, refgenome, filename, tmpdir='tmp'):
    # get the position of the first hit of the first 300 bases of refgenome in subjectgenome
    originPos, orientation = get_originPos_and_orientation_of_subjectgenome_(subjectgenome, refgenome, tmpdir)
    # reset origin of subjectgenome
    new_genome = reset_origin_of_genome_(subjectgenome, originPos, filename=filename, reverse = (orientation == 'reverse'), tmpdir = tmpdir)
    logging.info('new genome file with reset origin of '+subjectgenome+' saved to '+new_genome)
    return new_genome

# mask the genome with Ns based on the blastn results of clustered IS loci
def create_masked_genome_(genome, clustered_blast_hit_df, masked_genome):
    record = SeqIO.read(genome, 'fasta')
    genome_seq_array = np.array(list(str(record.seq)))
    
    # Mask regions specified in the dataframe
    for row in clustered_blast_hit_df.itertuples(index=False):
        if row.start < row.end:
            genome_seq_array[(row.start-1):(row.end)] = 'N'
        else:
            genome_seq_array[(row.end-1):(row.start)] = 'N'
    
    # Convert back to a string
    genome_seq_str = ''.join(genome_seq_array)

    # Write the modified sequence to a file
    new_seq = Seq(genome_seq_str)
    record = SeqRecord(new_seq, id='contig_1')
    SeqIO.write(record, masked_genome, 'fasta')
    logging.info('masked genome saved to '+masked_genome)    

def get_direction_(start, end, Crossing):
    if Crossing:
        return 'forward' if start > end else 'reverse'
    else:
        return 'forward' if start < end else 'reverse'

def get_middle_pos_(start, end, Crossing, genomesize):
    if Crossing:
        end = end + genomesize
        pos = (start + end) / 2
        return pos % genomesize
    else:
        return (start + end) / 2

def get_IS_positions_in_subjectgenome(blastn_df, minCompleteSize=2900, genomesize = 4e6):
    # Crossing: whether the position crosses the end of sequence
    # IS_strand: the strand of ISs in the subjectgenome
    # IS_pos: the position of ISs in the subjectgenome (middle of the IS)
    # Complete: whether the IS is complete or not > minCompleteSize
    blastn_df = blastn_df.copy(deep=False)
    blastn_df['Crossing'] = blastn_df.apply(lambda row: True if abs(row.s_start - row.s_end) > genomesize/2 else False, axis=1)  # Basically set to false because ori assignment should avoid crossing blastn_df.apply(lambda row: True if abs(abs(row.s_start - row.s_end)+1 - row.alignment_length) < 100 else False, axis=1)
    blastn_df['IS_strand'] = blastn_df.apply(lambda row: get_direction_(row.s_start, row.s_end, row.Crossing), axis=1)
    blastn_df['IS_pos'] = blastn_df.apply(lambda row: get_middle_pos_(row.s_start, row.s_end, row.Crossing, genomesize), axis=1)
    blastn_df['Complete'] = blastn_df.apply(lambda row: True if row.alignment_length >= minCompleteSize else False, axis=1)
    return blastn_df

def overlaps_or_near(s1, e1, s2, e2, max_gap = 10):
    return (s1 <= e2 + max_gap) and (e1 >= s2 - max_gap)

def cluster_close_ones(array, max_gap = 10):
    ser = pd.Series(array)
    # values are the size of the group with same values
    grp_ser = ser.groupby((ser.diff() !=0).cumsum()).transform('size')
    ser_copy = ser.copy()
    ser_copy.loc[(grp_ser < max_gap) & (ser.eq(0))] = 1
    res = ((ser_copy.diff() != 0) & (ser_copy != 0)).cumsum()
    res[ser_copy == 0] = 0
    return res

# cluster extremely close TE hits into one
"""
required columns in df_IS_pos: s_start, s_end, IS_strand, alignment_length
"""
def cluster_TE_hits(df_IS_pos, genome_size, max_gap = 10, seed_match_size = 300):
    # get zero one positions
    vector = create_raw_te_position_vec(df_IS_pos, genome_size)

    # cluster
    cluster_id_array = cluster_close_ones(vector, max_gap = max_gap)

    # create vector where range_start is smaller than range_end from sstart and ssend
    df = pd.DataFrame(columns=['range_start', 'range_end'])
    df['range_start'] = df_IS_pos.apply(lambda row: min(row.s_start, row.s_end), axis=1)
    df['range_end'] = df_IS_pos.apply(lambda row: max(row.s_start, row.s_end), axis=1)
    df = pd.concat([df, df_IS_pos[['IS_strand', 'IS_pos', 'alignment_length', 'Complete']]], axis=1)
    df = df.query('alignment_length >= '+str(seed_match_size)).reset_index(drop=True)
    df = df.sort_values(by = ['IS_pos']).reset_index(drop=True)

    # preserve only ranges that overlap with the seed ranges
    overlap_ranges_raw = pd.DataFrame()
    df['id'] = df.index
    for r in df.itertuples(index=False):
        overlap_ranges_raw = pd.concat([overlap_ranges_raw,
        pd.DataFrame({'id': [r.id, r.id],
         'end': ['start', 'end'],
         'cluster_id': [cluster_id_array[r.range_start], cluster_id_array[r.range_end]]
        })], axis=0)
    overlap_ranges_raw.reset_index(drop=True, inplace=True)

    overlap_range_df = pd.DataFrame()
    # For each new cluster, check if it includes any of the original seed ranges
    for cluster_id in overlap_ranges_raw['cluster_id'].unique():
        df_within_cluter = df.iloc[overlap_ranges_raw[overlap_ranges_raw['cluster_id'] == cluster_id]['id'].unique()]
        max_alignment = df_within_cluter['alignment_length'].max()
        row = df_within_cluter[df_within_cluter['alignment_length'] == max_alignment].iloc[0]
        overlap_range_df = pd.concat([overlap_range_df,
                    pd.DataFrame({'start': cluster_id_array[cluster_id_array == cluster_id].index.min(),
                    'end': cluster_id_array[cluster_id_array == cluster_id].index.max(),
                    'IS_strand': row['IS_strand'],
                    'max_alignment_length': row['alignment_length'],
                    'length': cluster_id_array[cluster_id_array == cluster_id].shape[0],
                    }, index = [0])], axis=0)
    overlap_range_df = overlap_range_df.reset_index(drop=True)
    overlap_range_df['cluster_id'] = overlap_range_df.index

    logging.info('cluters detected: '+str(overlap_range_df.shape[0]))
    
    return overlap_range_df

def overlap_but_not_in(s1, e1, s2, e2, max_gap = 10):
    s1_, e1_ = min(s1, e1), max(s1, e1)
    s2_, e2_ = min(s2, e2), max(s2, e2)
    if(overlaps_or_near(s1_, e1_, s2_, e2_, max_gap=max_gap) and not (s2_ - max_gap < s1_ and e1_ < e2_ + max_gap)):
        print(s1_, e1_, s2_, e2_)
        print(s1, e1, s2, e2)
    return overlaps_or_near(s1_, e1_, s2_, e2_, max_gap=max_gap) and not (s2_ - max_gap < s1_ and e1_ < e2_ + max_gap)

def subset_valid_te_hits(blastn_df, alignment_length_cutoff = 300, max_gap = 10):
    # preserve large is matches
    large_matches = blastn_df.query('alignment_length >= '+str(alignment_length_cutoff))

    # and small matches adjaent to large matches
    fragments = blastn_df.query('alignment_length < '+str(alignment_length_cutoff))
    fragments = fragments[fragments.apply(lambda frag:
        large_matches[large_matches.apply(lambda long:
            overlap_but_not_in(frag.s_start, frag.s_end, long.s_start, long.s_end, max_gap = max_gap)
            , axis=1)].shape[0] != 0, axis=1)]
    if(fragments.shape[0] > 0):
        logging.info('small matches adjacent to large matches detected')

    return pd.concat([large_matches, fragments], ignore_index=True).sort_values(by=['s_start']).reset_index(drop=True)

def create_raw_te_position_vec(df_IS_pos, genome_size):
    pos_vec = np.zeros(genome_size + 1)
    for row in df_IS_pos.itertuples(index=False):
        start, end = min(row.s_start, row.s_end), max(row.s_start, row.s_end)
        pos_vec[(start):(end+1)] = 1
    return pos_vec

# get the sequences matching ISs and convert TE matches longer than alginment_length_cutoff to Ns
def mask_TE_sequences(genome, TEFasta, tmpdir = 'tmp',
                      prefix = 'runSyRI_IS', genome_identifier = '', alginment_length_cutoff = 300,
                      genome_size = 4e6, IS_size = 3000, te_name = 'IS', max_cluster_gap = 20):
    dbdir = makedb(genome, 'origin resetted genome', tmpdir = tmpdir)
    blastn_prefix = os.path.join(tmpdir, prefix + '_genome.masked.blastn.' + (genome_identifier+'.fasta' if genome_identifier else 'fasta'))
    blastn = BLAST(TEFasta, genome, 1e-10, blastn_prefix, dbdir = dbdir)
    blastn_df = pd.read_table(blastn, comment='#', header=None)
    blastn_df =  blastn_df.set_axis(['subject_id', 'query_id', 'q_start', 'q_end', 's_start', 's_end', 'alignment_length', 'evalue'], axis=1)
    #blastn_df = subset_valid_te_hits(blastn_df, alignment_length_cutoff = alginment_length_cutoff, max_gap = max_cluster_gap)
    
    # set file names
    masked_genome = prefix + '_genome.masked.' + te_name + '.' + (genome_identifier+'.fasta' if genome_identifier else 'fasta')
    masked_genome = os.path.join(tmpdir, masked_genome)
    df_IS_pos = os.path.join(tmpdir, prefix + '_' + te_name + '_positions_in_subjectgenome.' + (genome_identifier+'.csv' if genome_identifier else 'csv'))
    df_IS_cluster_pos = os.path.join(tmpdir, prefix + '_clustered_' + te_name + '_positions_in_subjectgenome.' + (genome_identifier+'.csv' if genome_identifier else 'csv'))

    # get IS positions in subjectgenome
    blastn_df = get_IS_positions_in_subjectgenome(blastn_df, genomesize=genome_size, minCompleteSize=IS_size*0.9)
    clustered_is_df = cluster_TE_hits(blastn_df, genome_size, max_gap = max_cluster_gap, seed_match_size=alginment_length_cutoff)

    create_masked_genome_(genome, clustered_is_df, masked_genome) 

    # save
    blastn_df.query(f'alignment_length >= {alginment_length_cutoff}').to_csv(df_IS_pos, index=False)
    clustered_is_df.to_csv(df_IS_cluster_pos, index=False)
    logging.info('blastn results saved to '+ df_IS_pos)
    return masked_genome, df_IS_pos, df_IS_cluster_pos

def join_is_dfs(is_df_dirs, prefix, save_dir):
    is_dfs = []
    for i, is_df_dir in enumerate(is_df_dirs):
        is_df = pd.read_csv(is_df_dir)
        is_df['SyRI_genome_id'] = i
        is_df['SyRI_prefix'] = prefix
        is_dfs.append(is_df)
    is_df = pd.concat(is_dfs, ignore_index=True)
    is_df.to_csv(save_dir, index=False)
    logging.info('joined IS dfs saved to '+save_dir)
    return is_df

def align_two_genomes(target_genome, query_genome, outputdir = 'tmp/tmp.sam', threads = 1, output_type = 'sam', option = '-O4,100 -E2,0', preset = 'map-ont'): 
    C = 'minimap2 -ax ' + preset + ' --eqx ' + option + ' -t ' + str(threads) + ' ' + target_genome + ' ' + query_genome 
    if output_type == 'sam':
        C += ' > ' + outputdir
    elif output_type == 'bam':
        C += ' | samtools sort -O BAM - > ' +outputdir 
        C += ' ; samtools index ' + outputdir
    logging.info(C)
    os.system(C)
    return outputdir 
    
# syri -c bam/L01-2_G8_IS1.sam -r fasta/MDS42_IS1.fa -q fasta/L01-2_G8.fasta.restart -F S -k --prefix L01-2_G8_IS1 --invgaplen 1000 --tdmaxolp 0.2 --nosnp
def run_syri(samfile, refgenome, queryFasta, prefix, tmpdir = 'tmp', mapping_type = 'sam'):
    syri_save_dir = os.path.join(tmpdir, 'syri')
    os.makedirs(syri_save_dir, exist_ok=True)
    C = 'syri -c ' + samfile + ' -r ' + refgenome + ' -q ' + queryFasta
    C +=  ' --prefix ' + prefix + ' --allow-offset 5 --cigar --dir ' + syri_save_dir 
    C += ' -F ' + ('S' if mapping_type == 'sam' else 'B')
    C += ' --maxsize 100 '
    logging.info(C)
    os.system(C)
    return syri_save_dir, os.path.join(syri_save_dir, prefix + 'syri.out')

# sniffles -i R01-4_genome.1.is-inserted.bam -v is-ins-g1-2.vcf --minsupport 1 --no-qc --allow-overwrite
def run_sniffles(bamfile, vcf, minsupport = 1):
    C = 'sniffles -i ' + bamfile + ' -v ' + vcf + ' --minsupport ' + str(minsupport) + ' '
    C += ' --no-qc --allow-overwrite'
    logging.info(C)
    os.system(C)
    return vcf

def get_names_from_fastas(fastafiles):
    names = []
    for fastafile in fastafiles:
        name = os.path.basename(fastafile)
        name = os.path.splitext(name)[0]
        names.append(name)
    return names

def gen_genome_info_file_for_plotsr(genome_dirs, genome_names, formats=[], outputdir='tmp/plotsR_genomes.txt'):
    plotsR_genomes_df = pd.DataFrame(columns=['file', 'name', 'tags'])
    plotsR_genomes_df['file'] = genome_dirs
    plotsR_genomes_df['name'] = genome_names
    if formats:
        plotsR_genomes_df['tags'] = formats
    else:
        plotsR_genomes_df['tags'] = 'lw:1.5'
    plotsR_genomes_df.to_csv(outputdir, index=False, header=False, sep='\t')
    return outputdir

# get the AL sequences only to plot local alignments with plotsr instead of global alignments(default)
def create_local_syri_output(syrioutputfiles):
    local_syri_output = [file.replace('.out', '.local.out') for file in syrioutputfiles]
    for i, file in enumerate(syrioutputfiles):
        syri_output = pd.read_csv(file, sep='\t', header=None)
        syri_output.columns = ['chr1', 'start1', 'end1', 'REF', 'ALT', 'chr2', 'start2', 'end2', 'ID', 'parent', 'SVTYPE', 'MISC']
        # Get the aligned positions
        syri_output = syri_output.query('SVTYPE in ["SYNAL", "INVAL", "TRANSAL", "DUPAL", "INVTRAL", "INVDPAL"]').copy()
        syri_output['SVTYPE'] = syri_output['SVTYPE'].apply(lambda x: x.replace('AL', ''))
        last_r, last_q = syri_output['end1'].max(), syri_output['end2'].max()

        # not required as long as the ori is deleted
        whole_genome_dummy = pd.DataFrame({'chr1': ['chr1'], 'start1': [1], 'end1': [last_r],
        'REF': [''], 'ALT': [''], 'chr2': ['chr1'], 'start2': [1], 'end2': [last_q],
        'ID': ['SYN_dummy'], 'parent': ['-'], 'SVTYPE': ['SYN'], 'MISC': ['-']})
        #syri_output = pd.concat([syri_output, whole_genome_dummy], ignore_index=True)

        syri_output.to_csv(local_syri_output[i], sep='\t', index=False, header=False)
        logging.info('local alignment version of syri output saved to ' + local_syri_output[i])
    return local_syri_output

# change CPG CPL to DUP and DEL
def convert_syri_cpg_cpl(syrioutputfile):
    df = pd.read_csv(syrioutputfile, sep='\t', header=None)
    df.columns = ['chr1', 'start1', 'end1', 'REF', 'ALT', 'chr2', 'start2', 'end2', 'ID', 'parent', 'SVTYPE', 'MISC']
    df2_ = pd.DataFrame()
    for i, row in df.iterrows():
        if row.SVTYPE == 'CPG':
            row.SVTYPE = 'DUP'
            df2_ = pd.concat([df2_, pd.DataFrame([row])], ignore_index=True)
            row.SVTYPE = 'DUPAL'
            row.parent = row.ID
            row.ID = row.ID + '_AL'
            df2_ = pd.concat([df2_, pd.DataFrame([row])], ignore_index=True)
        elif row.SVTYPE == 'CPL':
            row.SVTYPE = 'DEL'
            df2_ = pd.concat([df2_, pd.DataFrame([row])], ignore_index=True)
        else:
            df2_ = pd.concat([df2_, pd.DataFrame([row])], ignore_index=True)
    df = df2_
    df.reset_index(drop=True, inplace=True)
    df.to_csv(syrioutputfile, sep='\t', index=False, header=False)
    logging.info('converted syri output saved to ' + syrioutputfile)
    return syrioutputfile

# using plotsr package
def draw_plot_with_plotsr(syrioutputfiles, plotsr_genome_file, tracks = '', markers = '', export_dir = 'export', prefix='', vertical = False, chrs = [], options = ''):
    savedir = os.path.join(export_dir, prefix + '_plotsr.pdf')
    C = 'plotsr '
    for syrioutputfile in syrioutputfiles:
        C += ' --sr ' + syrioutputfile
    C += ' --genomes ' + plotsr_genome_file 
    C += (' --tracks ' + tracks if tracks else '') + (' --markers ' + markers if markers else ' ')
    C += ' -s 100 ' # minimum size of a SR to be plotted 
    for chr in chrs:
        C += ' --chr ' + chr # Select specific chromosome on reference (first genome) and plots them in the given order. Not compatible with --chrord. Can be used multiple time to select more than one chromosomes.
    if vertical:
        C += ' -v '
        savedir = savedir.replace('.pdf', '.vertical.pdf')
    C += ' -o ' + savedir
    C += f' {options}'
    logging.info(C)
    os.system(C)
    logging.info('plotsr output saved to ' + savedir)
    return savedir

# using plotsr_ forked from plotsr v1.1.0
def draw_plot_with_plotsr_(syrioutputfiles, plotsr_genome_file, tracks = '', markers = '', export_dir = 'export', prefix='', vertical = False, chrs = [], options = ''):
    savedir = os.path.join(export_dir, prefix + '_plotsr.pdf')
    C = ''
    for syrioutputfile in syrioutputfiles:
        C += ' --sr ' + syrioutputfile
    C += ' --genomes ' + plotsr_genome_file 
    C += (' --tracks ' + tracks if tracks else '') + (' --markers ' + markers if markers else ' ')
    C += ' -s 100 ' # minimum size of a SR to be plotted 
    for chr in chrs:
        C += ' --chr ' + chr # Select specific chromosome on reference (first genome) and plots them in the given order. Not compatible with --chrord. Can be used multiple time to select more than one chromosomes.
    if vertical:
        C += ' -v '
        savedir = savedir.replace('.pdf', '.vertical.pdf')
    C += ' -o ' + savedir
    C += f' {options}'
    args_list = C.split()
    print(args_list)
    logging.info('plotsr ' + ' '.join(args_list))
    plotsr_main(args_list)
    logging.info('plotsr output saved to ' + savedir)
    return savedir

# example: mt:v;mc:black;ms:3;tt:Inversion 1;tp:0.02;ts:8;tf:Arial;tc:black
def generate_IS_tag_for_plotsr(strand, size=3, markercolor = 'dimgrey', text='', vertical = False):
    tags = pd.DataFrame(columns=['id','value'])
    if not vertical:
        tags = pd.concat([tags, pd.DataFrame({'id':'mt', 'value':'i2' if strand == 'forward' else 'i3'}, index=[0])], ignore_index=True)
    else:
        tags = pd.concat([tags, pd.DataFrame({'id':'mt', 'value':'i0' if strand == 'forward' else 'i1'}, index=[0])], ignore_index=True)
    tags = pd.concat([tags, pd.DataFrame({'id':'ms', 'value':str(size)}, index=[0])], ignore_index=True)
    if markercolor:
        tags = pd.concat([tags, pd.DataFrame({'id':'mc', 'value':markercolor}, index=[0])], ignore_index=True)
    if text:
        tags = pd.concat([tags, pd.DataFrame({'id':'tt', 'value':text}, index=[0])], ignore_index=True)
    
    # join each row with : and join all rows with ;
    tags = ';'.join(tags.apply(lambda row: ':'.join(row), axis=1))
    return tags

def generate_IS_tag_df(is_dfs, plotsR_genomeFile, tmpdir = 'tmp', prefix = 'runSyRI_IS'):
    is_tag_df = pd.DataFrame(columns=['chr', 'start', 'end', 'genome_id', 'tags'])
    plotsR_df = pd.read_csv(plotsR_genomeFile, sep='\t', header=None)
    genome_ids = plotsR_df[1].tolist()
    for i in range(len(is_dfs)):
        #cols:subject_id,query_id,q_start,q_end,s_start,s_end,alignment_length,evalue
        is_tag_df_ = pd.DataFrame(columns=['chr', 'start', 'end', 'genome_id', 'tags'])
        is_df = pd.read_csv(is_dfs[i])
        is_tag_df_['start'] = is_df['IS_pos'].astype(int)-1
        is_tag_df_['end'] = is_df['IS_pos'].astype(int)
        is_tag_df_['genome_id'] = genome_ids[i]
        is_tag_df_['chr'] = 'contig_1'
        is_tag_df_['tags'] = is_df.apply(lambda row: generate_IS_tag_for_plotsr(row.IS_strand), axis=1)
        is_tag_df = pd.concat([is_tag_df, is_tag_df_], ignore_index=True)
    is_tag_dir = os.path.join(tmpdir, prefix + '_IS_tags_for_plotsr.bed')
    is_tag_df.to_csv(is_tag_dir, index=False, header=False, sep='\t')
    logging.info('IS tags for plotsr saved to ' + is_tag_dir)
    return is_tag_dir

# output: s_start, s_end, flank, flank_id
def write_fasta_of_flank_is_pair(cluster_TE_df, ref_fasta, fasta_out, flank_length=100):
	ref_fasta = [record.seq for record in SeqIO.parse(ref_fasta, 'fasta')][0]
	records = []
	flank_df = pd.DataFrame(columns=['is_id_desc', 'flank_id', 's_start', 's_end', 'flank'])
	for i, row in cluster_TE_df.iterrows():
		start = row['start']
		end = row['end']
		dir = 'f' if row['IS_strand'] == 'forward' else 'r'
		records.append(SeqRecord(
			Seq(ref_fasta[(start-flank_length+1):end]),
			id = str(row['cluster_id']) + '_' + dir,
			description = 'flank'
		))
		flank_df = pd.concat([flank_df, pd.DataFrame({
			'is_id_desc': row['cluster_id'],
			'flank_id': str(row['cluster_id']) + '_' + dir,
			's_start': start-flank_length+1,
			's_end': end,
			'flank': dir,
		}, index=[0])])
		dir = 'r' if row['IS_strand'] == 'forward' else 'f'
		records.append(SeqRecord(
			Seq(ref_fasta[start:(end+flank_length-1)]).reverse_complement(),
			id = str(row['cluster_id']) + '_' + dir,
			description = 'flank'
		))
		flank_df = pd.concat([flank_df, pd.DataFrame({
			'is_id_desc': row['cluster_id'],
			'flank_id': str(row['cluster_id']) + '_' + dir,
			's_start': end + flank_length - 1,
			's_end': start,
			'flank': dir, 
		}, index=[0])])
	SeqIO.write(records, fasta_out, 'fasta')
	return flank_df, fasta_out# output: s_start, s_end, flank, flank_id
