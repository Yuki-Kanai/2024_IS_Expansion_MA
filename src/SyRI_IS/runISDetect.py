import logging
import pandas as pd
import os
import sys
import ISDetect
import argparse

#### Rewrote based on TELR

def get_args():
    parser = argparse.ArgumentParser(description='Transposon insertion locus detector for a single IS.')
    parser.add_argument('fastqDir', help = 'directory of fastq file. If a fasta file is given, use it as fastaFile')
    parser.add_argument('refgenome', help = '(.fasta)')
    parser.add_argument('ISseqfile', help = '(.fasta)')
    parser.add_argument('ISancPos', help = 'pandas df of IS position in ancestral genome. cols = IS_strand, count, medianpos, insert_id, gene, (.csv)')
    parser.add_argument('geneDF', help = 'dataframe of gene essentiallity (.csv)')
    parser.add_argument('cogDF', help = 'dataframe of gene functional annotation (.csv)')

	#optional
    parser.add_argument('-p', '--prefix', type=str, default='ISdetect')
    parser.add_argument('-t', '--threads', type=int, default=4)
    parser.add_argument('-o', '--out', type=str, default='export')
    parser.add_argument('--tmp', type=str, default='tmp')
    parser.add_argument('--log', type=str, default='log', help='log directory')
    parser.add_argument('--flank', type=str, default = '', help='sequences flanking IS to use as query for blastn in multifasta format')

    return parser.parse_args()

def set_logging(logfile):
    formatstr = "%(asctime)s: %(levelname)s: %(message)s"
    datestr = "%m/%d/%Y %H:%M:%S"
    logging.basicConfig(
        level=logging.DEBUG,
        filename=logfile,
        filemode="w",
        format=formatstr,
        datefmt=datestr,
    )
    logging.info("Start running ISDetect")
    logging.info("CMD: " + " ".join(sys.argv))

def main():
    args = get_args()
    os.makedirs(args.log, exist_ok=True)
    set_logging(os.path.join(args.log, args.prefix+".ISdetect.log"))

    prefix = args.prefix
    fastqDir = args.fastqDir
    refgenome = args.refgenome
    ISfile = args.ISseqfile

    tmpDir = args.tmp
    expDir = args.out
    threads = args.threads
    bamdir = os.path.join(tmpDir, prefix+'.sort.bam')
    os.makedirs(tmpDir, exist_ok=True)
    os.makedirs(expDir, exist_ok=True)

    fastaFile = os.path.join(tmpDir, prefix+'_raw.fasta')
    geneDF = pd.read_csv(args.geneDF)
    refgenomeLength = ISDetect.getFastaLength(refgenome)

    if fastqDir.endswith(('.fasta', '.fa', '.fna', '.faa')):
        ISDetect.MakeCopyOfFastaFile(fastaFile, fastqDir, n=3)
        logging.info('Detected fasta file, skip converting fastq to fasta')
    else :
        ISDetect.FastqDir2Fasta(fastaFile, fastqDir)
    
    pseudoGeneFile = os.path.join(expDir, prefix+'-pseudoGenes.csv')
    if not os.path.exists(pseudoGeneFile):
        insertions, _ = ISDetect.GetInsertions(
            fastaFile,
            refgenome,
            ISfile,
            prefix,
            threads = threads,
            tempFolder = tmpDir,
            distance_threashold=20, # much larger than typical duplication
            minSupportNum=3,
            flank_file = args.flank,
        )

        originalISDF = pd.read_csv(args.ISancPos)
        if not all([x in originalISDF.columns for x in ['IS_strand', 'count', 'medianpos', 'insert_id', 'gene']]):
            raise ValueError('ISancPos file must have columns: IS_strand, count, medianpos, insert_id, gene')

        insertions = ISDetect.IdentifyGeneInInsertedPos(
            insertions,
            geneDF,
            savedir = pseudoGeneFile,
            save = True,
            originalIS = originalISDF, # add original gpma galk to output
        )
        logging.info('got insert positions')
    else : 
        logging.info('skip blast searches because pseudogeneFile exists.')

    # from this point onwards, duplication of the fasta file is not required
    if fastqDir.endswith(('.fasta', '.fa', '.fna', '.faa')):
        fastaFile = fastqDir

	# get the depths between IS positions to detect copy number changes
    cogDF = pd.read_csv(args.cogDF)
    insertions = pd.read_csv(pseudoGeneFile)
    ISDetect.PlotInsertions(insertions, savePrefix = os.path.join(tmpDir,prefix), genomeLength=refgenomeLength)
    intervalDepth, depthRaw = ISDetect.CalcDepthInterval(insertions, bamdir,prefix, tmpDir, fastaFile, refgenome,
                                                         skipalignment=os.path.exists(bamdir), threads = threads, refgenomeLength=refgenomeLength)

    logging.info('got depth data')

    # get the genes that are affected by copy number changes
    copyNumberChangedGenes, copyNumberChangedIntervals = ISDetect.IdentifyCopyNumberChangedGenes(
        intervalDepth, depthRaw, geneDF,
        save = True, savePrefix = os.path.join(expDir, prefix)
    )
    logging.info('got tandem deletion and duplications')

    # get the COG assignment of those genes
    CNgeneCOGDF = ISDetect.COGDFofDeleteDuplicatedGenes(
        cogDF,
        insertions,
        copyNumberChangedGenes,
        savePrefix = os.path.join(expDir,prefix),
    )
    logging.info('got COG of deleted and duplicated genes')

	# Detect SVs related to the ISs
    rawFlankingPosDF = pd.read_csv(os.path.join(tmpDir, prefix+'-FlkPosDF.csv'))
    ISDetect.GetFlankingPairs(
        rawFlankingPosDF,
        savePrefix = os.path.join(expDir,prefix),
        save=True
        )
    logging.info('got flanking pairs')

    logging.info('finished running ISDetect')

if __name__ == "__main__":
    main()