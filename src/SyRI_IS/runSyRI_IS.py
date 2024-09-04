import logging
import os
import sys
import argparse
import SyRI_IS
import is_event_classification_func as iecf
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser(description='Transposon insertion locus detector for a single IS.')
    parser.add_argument('-g', '--genome', type=str, help = '<required> The sequences of subject genome. First in would be used as the reference of the second input (.fasta)', action='append', required=True)
    parser.add_argument('-i', '--ISseq', help = '<required> (.fasta)', required=True)

	#optional parameters
    parser.add_argument('-p', '--prefix', type=str, default='SyRI_IS')
    parser.add_argument('-t', '--threads', type=int, default=4)
    parser.add_argument('-o', '--out', type=str, default='export')
    parser.add_argument('--log', type=str, default='log', help='log directory')
    parser.add_argument('--tmp', type=str, default='tmp')
    parser.add_argument('--replicationOrigin', type=int, default=0, help='replication origin position in the reference genome. If not given, 0 is used.')

    # force rerun flags
    # By default, the program will skip the steps that have been done. If you want to rerun the steps, use these flags.
    parser.add_argument('--rerun_reset_origin', action='store_true', help='rerun reset_origin_of_input_genomes') 
    parser.add_argument('--rerun_mask_TE_sequences', action='store_true', help='rerun get_is_pos_and_masked_genomes')
    parser.add_argument('--rerun_align_two_genomes', action='store_true', help='rerun align_pairs_of_genomes')
    parser.add_argument('--rerun_run_syri', action='store_true', help='rerun run_multiple_syri')
    parser.add_argument('--rerun', action='store_true', help='rerun all the steps')

    return parser.parse_args()

def set_logging(logfile):
    formatstr = "%(asctime)s: %(levelname)s: %(message)s"
    datestr = "%m/%d/%Y %H:%M:%S"
    logging.basicConfig(
        level=logging.DEBUG,
        filename=logfile,
        filemode="a",
        format=formatstr,
        datefmt=datestr,
    )
    logging.info("Start running SyRI_IS")
    logging.info("CMD: " + " ".join(sys.argv))

## Check if the input arguments are valid
def check_args(args):
    assert args.prefix == args.prefix.replace('.', ''), 'prefix should not include "."'
    # check if all the genome file names exists and there base names are different
    assert len(args.genome) >= 2, 'At least two genomes should be provided.'
    for genome in args.genome:
        assert os.path.isfile(genome), genome + ' does not exist.'
    base_names = [os.path.basename(genome) for genome in args.genome]
    assert len(base_names) == len(set(base_names)), 'Genome file names should be different.'
    if args.rerun:
        args.rerun_reset_origin = True
        args.rerun_mask_TE_sequences = True
        args.rerun_align_two_genomes = True
        args.rerun_run_syri = True

## Reset the origin of the input genomes so that the replication origin is at the beginning of the genome
def reset_origin_of_input_genomes(genomes, replicationOrigin, tmpdir = 'tmp', prefix = 'SyRI_IS', rerun = False):
    # reset origin of reference
    origin_resetted_genome_names = [prefix + '_genome.'+str(i)+'.fasta' for i in range(len(genomes))]
    origin_resetted_genomes = [os.path.join(tmpdir, name) for name in origin_resetted_genome_names]

    if not os.path.isfile(origin_resetted_genomes[0]) or rerun:
        saved_file = SyRI_IS.reset_origin_of_genome_(genomes[0], replicationOrigin, origin_resetted_genome_names[0], tmpdir = tmpdir)
        assert saved_file == origin_resetted_genomes[0], 'saved_file: '+saved_file+' origin_resetted_genomes[0]: '+origin_resetted_genomes[0]
        logging.info('origin resetted genome saved to '+ origin_resetted_genomes[0])
    else:
        logging.info('origin resetted genome found at '+ origin_resetted_genomes[0])
    for i in range(len(genomes)-1):
        if not os.path.isfile(origin_resetted_genomes[i+1]) or rerun:
            saved_file = SyRI_IS.reset_origin_of_subject_genome(genomes[i+1], origin_resetted_genomes[i], origin_resetted_genome_names[i+1], tmpdir = tmpdir)
            assert saved_file == origin_resetted_genomes[i+1], 'saved_file: '+saved_file+' origin_resetted_genomes[i+1]: '+origin_resetted_genomes[i+1]
            logging.info('origin resetted genome saved to '+ origin_resetted_genomes[i+1])
        else:
            logging.info('origin resetted genome found at '+ origin_resetted_genomes[i+1])
    return origin_resetted_genomes

def get_is_pos_and_masked_genomes(origin_resetted_genomes, ISseq, tmpdir = 'tmp', prefix = 'SyRI_IS', alginment_length_cutoff = 300, rerun = False):
    ngenome = len(origin_resetted_genomes)
    masked_genomes = [prefix + '_genome.masked.IS.'+str(i)+'.fasta' for i in range(ngenome)]
    masked_genomes = [os.path.join(tmpdir, name) for name in masked_genomes]
    is_dfs = [prefix + '_IS_positions_in_subjectgenome.'+str(i)+'.csv' for i in range(ngenome)]
    is_dfs = [os.path.join(tmpdir, name) for name in is_dfs]
    genome_sizes = [SyRI_IS.get_fasta_length(genome) for genome in origin_resetted_genomes]
    IS_size = SyRI_IS.get_fasta_length(ISseq)
    clustered_DF = [prefix + '_clustered_IS_positions_in_subjectgenome.'+str(i)+'.csv' for i in range(ngenome)]
    clustered_DF = [os.path.join(tmpdir, name) for name in clustered_DF]

    logging.warning('genome_sizes: '+str(genome_sizes))

    for i in range(ngenome):
        if not (os.path.isfile(is_dfs[i]) and os.path.isfile(masked_genomes[i]) and os.path.isfile(clustered_DF[i])) or rerun:
            saved_masked_genome, saved_IS_DF, saved_clustered_DF = SyRI_IS.mask_TE_sequences(origin_resetted_genomes[i],
                ISseq, tmpdir = tmpdir, prefix = prefix,
                alginment_length_cutoff = alginment_length_cutoff,
                genome_identifier = str(i), genome_size = genome_sizes[i], IS_size = IS_size)
            assert saved_masked_genome == masked_genomes[i]
            assert saved_IS_DF == is_dfs[i]
            assert saved_clustered_DF == clustered_DF[i]
            logging.info('ISs masked genome '+str(i)+' generated and saved to '+masked_genomes[i])
        else:
            logging.info('ISs masked genome '+str(i)+' found at '+masked_genomes[i])
    return is_dfs, masked_genomes, clustered_DF

def align_pairs_of_genomes(genomes, prefix, threads, rerun = False, folder = ''):
    ngenome = len(genomes)
    bamfiles = [file_name.replace('.fasta', '.prev_align.bam') for file_name in genomes[1:]]
    if folder != '':
        bamfiles = [os.path.join(folder, os.path.basename(bamfile)) for bamfile in bamfiles]
        os.makedirs(folder, exist_ok=True)
    for i in range(ngenome-1):
        if not os.path.isfile(bamfiles[i]) or rerun:
            saved_bamfile = SyRI_IS.align_two_genomes(genomes[i], genomes[i+1], bamfiles[i], threads = threads,
                            preset = 'asm20', option = '-H -f 100 -r1k,10k --rmq=no', output_type = 'bam')
            assert saved_bamfile == bamfiles[i]
            logging.info('bamfile '+str(i)+' generated and saved to '+bamfiles[i])
        else:
            logging.info('bamfile '+str(i)+' found at '+bamfiles[i])
    return bamfiles

def run_multiple_syri(bamfiles, masked_genomes, prefix, tmpdir = 'tmp', rerun = False):
    syri_dir = os.path.join(tmpdir, 'syri')
    os.makedirs(syri_dir, exist_ok=True)
    syri_out_files = [os.path.join(syri_dir, prefix + '.' + str(i) + '.syri.out') for i in range(len(masked_genomes)-1)]
    assert len(bamfiles) == len(masked_genomes) - 1
    for i in range(len(bamfiles)):
        if not os.path.isfile(syri_out_files[i]) or rerun:
            tmp_, saved_syri_out_file = SyRI_IS.run_syri(bamfiles[i], masked_genomes[i], masked_genomes[i+1], prefix+'.'+str(i)+'.', tmpdir = tmpdir, mapping_type = 'bam')
            assert saved_syri_out_file == syri_out_files[i]
            logging.info('syri output '+str(i)+' generated and saved to '+syri_out_files[i])
        else:
            logging.info('syri output '+str(i)+' found at '+syri_out_files[i])
    return syri_dir, syri_out_files

def create_is_deleted_genomes(genome_fastas, is_dfs, outdir, prefix, rerun = True):
    os.makedirs(outdir, exist_ok=True)
    deleted_is_genomes = [prefix + '_genome.deleted_IS.'+str(i)+'.fasta' for i in range(len(genome_fastas))]
    deleted_is_genomes = [os.path.join(outdir, name) for name in deleted_is_genomes]
    coor_trans_fws = [genome_fastas[i].replace('.fasta', '.fw_coord_transformation.arrow') for i in range(len(genome_fastas))]
    coor_trans_rvs = [genome_fastas[i].replace('.fasta', '.rv_coord_transformation.arrow') for i in range(len(genome_fastas))]
    for i in range(len(genome_fastas)):
        if not os.path.isfile(deleted_is_genomes[i]) or rerun:
            is_df = pd.read_csv(is_dfs[i]).rename(columns={'cluster_id': 'insert_id'})
            saved_deleted_is_genome = iecf.delete_is_positions(is_df, genome_fastas[i], deleted_is_genomes[i], rewrite=rerun)
            assert saved_deleted_is_genome == deleted_is_genomes[i]
            logging.info('deleted IS genome '+str(i)+' generated and saved to '+deleted_is_genomes[i])
            saved_vecs = iecf.save_coord_transformation_vecs(is_df, genome_fastas[i], rewrite = rerun)
            assert saved_vecs[0] == coor_trans_fws[i] 
            logging.info('saved coord tran vecs')
        else:
            logging.info('deleted IS genome '+str(i)+' found at '+deleted_is_genomes[i])
    return deleted_is_genomes, coor_trans_fws, coor_trans_rvs

def reverse_lift_up_syri_outs(syri_save_dirs_isrm, coor_trans_rvs):
    assert len(syri_save_dirs_isrm) == len(coor_trans_rvs)-1, 'len(syri_save_dirs_isrm) != len(coor_trans_rvs)-1 : ' + f'{len(syri_save_dirs_isrm)} {len(coor_trans_rvs)-1}'
    new_syri_save_dirs = []
    for rid in range(len(syri_save_dirs_isrm)):
        r_transform = pd.read_feather(coor_trans_rvs[rid])
        q_transform = pd.read_feather(coor_trans_rvs[rid+1])
        new_syri_save_dir = iecf.reverse_lift_up_syri_out(syri_save_dirs_isrm[rid],
                            r_transform, q_transform) 
        SyRI_IS.convert_syri_cpg_cpl(new_syri_save_dir)
        new_syri_save_dirs += [new_syri_save_dir]
    logging.info('syri output lifted up')
    return new_syri_save_dirs

def main():
    args = get_args()
    check_args(args)
    os.makedirs(args.log, exist_ok=True)
    set_logging(os.path.join(args.log, args.prefix+".log"))
    tmpdir = os.path.join(args.tmp, args.prefix)

    os.makedirs(tmpdir, exist_ok=True)
    os.makedirs(args.out, exist_ok=True)
    bamdir = os.path.join(tmpdir, 'bam')
    os.makedirs(bamdir, exist_ok=True)

    origin_resetted_genomes = reset_origin_of_input_genomes(args.genome, args.replicationOrigin, tmpdir = tmpdir, prefix=args.prefix, rerun = args.rerun_reset_origin)
    is_dfs, masked_genomes, clustered_df = get_is_pos_and_masked_genomes(origin_resetted_genomes, args.ISseq, tmpdir = tmpdir, prefix = args.prefix, alginment_length_cutoff = 300, rerun = args.rerun_mask_TE_sequences)

    # output IS positions
    is_position_dir = os.path.join(args.out, 'is_position_files')
    os.makedirs(is_position_dir, exist_ok=True)
    SyRI_IS.join_is_dfs(is_dfs, args.prefix, os.path.join(is_position_dir, args.prefix + '_IS_positions.csv'))

    # run SV detection by SyRI
    bamfiles = align_pairs_of_genomes(masked_genomes, args.prefix, args.threads, rerun = args.rerun_align_two_genomes, folder = os.path.join(tmpdir, 'bam_is_masked'))
    syri_save_dirs, syri_out_files = run_multiple_syri(bamfiles, masked_genomes, args.prefix, tmpdir = tmpdir, rerun = args.rerun_run_syri)
    local_syri_out_files = SyRI_IS.create_local_syri_output(syri_out_files)
    plotsR_genomeFile = SyRI_IS.gen_genome_info_file_for_plotsr(masked_genomes, SyRI_IS.get_names_from_fastas(args.genome), formats=[], outputdir=os.path.join(tmpdir, args.prefix + '_plotsR_genomes.txt'))
    is_tag_dir = SyRI_IS.generate_IS_tag_df(is_dfs, plotsR_genomeFile, tmpdir = tmpdir, prefix=args.prefix)

    # Plot SVs
    plotsr_dir = os.path.join(args.out, 'plotsr')
    os.makedirs(plotsr_dir, exist_ok=True)
    SyRI_IS.draw_plot_with_plotsr(local_syri_out_files, plotsR_genomeFile, prefix=args.prefix, export_dir = plotsr_dir, markers = is_tag_dir)

    # SyRI analysis using IS removed genomes
    del_is_dir = os.path.join(tmpdir, 'is_deleted')
    os.makedirs(del_is_dir, exist_ok=True)
    deleted_is_genomes, coor_trans_fws, coor_trans_rvs = create_is_deleted_genomes(origin_resetted_genomes, clustered_df, del_is_dir, args.prefix, rerun = True)
    bamfiles_isrm = align_pairs_of_genomes(deleted_is_genomes, args.prefix, args.threads, rerun = True, folder = os.path.join(tmpdir, 'bam'))
    syri_save_dirs_isrm, syri_out_files_isrm = run_multiple_syri(bamfiles_isrm, deleted_is_genomes, args.prefix, tmpdir = del_is_dir, rerun = True)
    syri_out_files_isrm = reverse_lift_up_syri_outs(syri_out_files_isrm, coor_trans_rvs)
    local_syri_out_files_isrm = SyRI_IS.create_local_syri_output(syri_out_files_isrm)
    plotsR_genomeFile_isrm = SyRI_IS.gen_genome_info_file_for_plotsr(masked_genomes, SyRI_IS.get_names_from_fastas(args.genome), formats=[], outputdir=os.path.join(tmpdir, args.prefix + '_plotsR_genomes.txt'))
    
    # Plot SVs
    os.makedirs(os.path.join(args.out, 'plotsr-isrm'), exist_ok=True)
    SyRI_IS.draw_plot_with_plotsr(local_syri_out_files_isrm, plotsR_genomeFile, prefix=args.prefix+'_isrm', export_dir = os.path.join(args.out, 'plotsr-isrm'), markers = is_tag_dir)

if __name__ == "__main__":
    main()