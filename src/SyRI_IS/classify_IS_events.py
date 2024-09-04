import logging
import os
import sys
import argparse
import pandas as pd

import is_event_classification_func as iecf
import SyRI_IS

start_message = "Start analyzing using SyRI_IS output"

def get_args():
	parser = argparse.ArgumentParser(description='Transposon insertion locus detector for a single IS.')
	parser.add_argument('-p', '--prefix', type=str, default='SyRI_IS')
	parser.add_argument('--input', type=str, default='tmp', required = True, help='tmp directory of SyRI_IS run to read output from')
	parser.add_argument('-r', '--rid', type=int, required=True, help='genome id of SyRI_IS run to read output from')
	parser.add_argument('-q', '--qid', type=int, required=True, help='genome id of SyRI_IS run to read output from')
	parser.add_argument('--genemap', type=str, help='folder containing gene coordinate in anc', default='')
	parser.add_argument('--calc_depth_interval', action='store_true', help='calculate depth interval from bam data')

	#optional parameters
	parser.add_argument('-t', '--threads', type=int, default=4)
	parser.add_argument('-o', '--out', type=str, default='export')
	parser.add_argument('--log', type=str, default='log', help='log directory')

	# force rerun flags
	return parser.parse_args()


def check_args(args):
	return 0

def main():
	args = get_args()
	check_args(args)
	os.makedirs(os.path.join(args.log, 'classify_is_events'), exist_ok=True)
	set_logging(os.path.join(args.log, 'classify_is_events', args.prefix+".classfy_is_events.log"))
	tmpdir = os.path.join(args.input, 'classify_is_events')
	os.makedirs(tmpdir, exist_ok=True)
	outdir = os.path.join(args.out, args.prefix)
	os.makedirs(outdir, exist_ok=True)

	rid = str(args.rid)
	qid = str(args.qid)
	r_fasta = os.path.join(args.input, args.prefix + '_genome.'+rid+'.fasta')
	q_fasta = os.path.join(args.input, args.prefix + '_genome.'+qid+'.fasta')
	r_is_cluster = pd.read_csv(os.path.join(args.input, args.prefix +'_clustered_IS_positions_in_subjectgenome.'+ rid + '.csv'))
	q_is_cluster = pd.read_csv(os.path.join(args.input, args.prefix +'_clustered_IS_positions_in_subjectgenome.'+ qid + '.csv'))
	is_flank_seq_pos_in_ref_genome_df, is_classified_df, ref_genome_is_pos_status = iecf.get_is_flank_seq_pos_in_ref_genome(r_fasta, q_fasta, r_is_cluster, q_is_cluster, args.input, args.prefix, args.qid)

	#debug
	is_classified_df.to_csv(os.path.join(outdir, args.prefix + '_genome.'+qid+'.is_classified.raw.csv'), index=False)
	is_classified_df = iecf.annotate_simple_insertions_of_q_is(is_classified_df)
	is_flank_seq_pos_in_ref_genome_df, is_classified_df, ref_genome_is_pos_status = iecf.merge_ins_id_info_to_ref_genome_is_pos_status(is_flank_seq_pos_in_ref_genome_df, is_classified_df, ref_genome_is_pos_status)

	# create vectors showing the IS residing in the position
	# L02-2_2_merged_1_2_is_pos_conversion.arrow etc
	is_pos_conv_vec_save_path = os.path.join(tmpdir, args.prefix + '_is_pos_conversion_vecs_' + qid)
	is_pos_conv_vec_dirs = iecf.save_is_pos_conversion_vecs(r_is_cluster,
		q_is_cluster, is_flank_seq_pos_in_ref_genome_df,
		SyRI_IS.get_fasta_length(r_fasta), SyRI_IS.get_fasta_length(q_fasta),
		is_pos_conv_vec_save_path, f'{args.prefix}_{qid}', rid, qid,
		buffer_sizes = [0, 20, 100], skip_wrote=False)

	# annotate simple transpositions
	use_syri_out = False
	syri_out_dir = os.path.join(args.input, 'is_deleted', 'syri', f'{args.prefix}.{rid}.syri.out')
	if int(rid)+1 != int(qid): 
		syri_out_dir = os.path.join(args.input, 'is_deleted', 'syri', f'{args.prefix}.{rid}_{qid}.syri.out')
	if os.path.isfile(syri_out_dir):
		logging.info(f'Found syri file {syri_out_dir}')
		use_syri_out = True
		syri_out = iecf.read_syri_out(syri_out_dir)
		syri_out = iecf.remove_short_sv_syri_out(syri_out, 5)
		r_coor_tran_rv = pd.read_feather(os.path.join(args.input, f'{args.prefix}_genome.{rid}.rv_coord_transformation.arrow'))
		q_coor_tran_rv = pd.read_feather(os.path.join(args.input, f'{args.prefix}_genome.{qid}.rv_coord_transformation.arrow'))
		iecf.reverse_lift_up_syri_out_(syri_out, r_coor_tran_rv, q_coor_tran_rv, suffix = '_is_rm', add_rows=True)
		syri_out = iecf.classify_sv_boundary_is_in_SyRI_out(syri_out, is_pos_conv_vec_dirs[f'merged_{rid}_{qid}'], is_pos_conv_vec_dirs[f'query_{qid}'])
		iecf.annotate_simple_transpositions(is_flank_seq_pos_in_ref_genome_df, is_classified_df, ref_genome_is_pos_status, syri_out)
		syri_out.to_csv(os.path.join(outdir, args.prefix + '_genome.'+qid+'.syri_out.out'), index=False, sep='\t')
		#print(syri_out.loc[:,['SVTYPE', 'start1_is_rm', 'end1_is_rm', 'start2_is_rm', 'end2_is_rm', 'start1', 'end1', 'start2', 'end2']])
		syri_out['SVLEN_ISDEL'] = syri_out.apply(lambda x: int(x.end1_is_rm) - int(x.start1_is_rm) + 1 if x.SVTYPE != 'INS' and str(x.start1_is_rm).isdigit() else abs(int(x.end2_is_rm) - int(x.start2_is_rm)) + 1, axis=1)
		syri_out.to_csv(os.path.join(outdir, args.prefix + '_genome.'+qid+'.syri_out.out'), index=False, sep='\t')
	else :
		logging.info(f'No syri file {syri_out_dir} found. Skipping simple transposition annotation.')
	simple_insertion_distances = iecf.analyze_simple_insertion_distances(ref_genome_is_pos_status)

	# save
	is_flank_seq_pos_in_ref_genome_df.to_csv(os.path.join(outdir, args.prefix + '_genome.'+qid+'.is_flank_seq_pos_in_ref_genome.csv'), index=False)
	is_classified_df.to_csv(os.path.join(outdir, args.prefix + '_genome.'+qid+'.is_classified.csv'), index=False)
	ref_genome_is_pos_status.to_csv(os.path.join(outdir, args.prefix + '_genome.'+qid+'.ref_genome_is_pos_status.csv'), index=False)
	simple_insertion_distances.to_csv(os.path.join(outdir, args.prefix + '_genome.'+qid+'.simple_insertion_distances.csv'), index=False)

	if use_syri_out:
		# annotate other SVs
		sv_df = iecf.fix_syri_sv_boundaries_of_syri_out(syri_out, r_is_cluster, q_is_cluster,
						is_flank_seq_pos_in_ref_genome_df, is_classified_df, ref_genome_is_pos_status)
		# merge with insertions
		sv_df.to_csv(os.path.join(outdir, args.prefix + '_genome.'+qid+'.sv_df.csv'), index=False)

	if args.calc_depth_interval:
		# calculate depth interval
		bam_dir = os.path.join(args.input, 'bam', f'{args.prefix}_genome.deleted_IS.{rid}_{qid}.bam')
		if int(rid)+1 == int(qid):
			bam_dir = os.path.join(args.input, 'bam', f'{args.prefix}_genome.deleted_IS.{qid}.prev_align.bam')
		if os.path.isfile(bam_dir):
			depth_interval_df = iecf.get_depth_intervals_from_bam(bam_dir)
			is_pos_conv_vec_ref = pd.read_feather(is_pos_conv_vec_dirs[f'merged_{rid}_{qid}'])
			r_coor_tran_rv = pd.read_feather(os.path.join(args.input, f'{args.prefix}_genome.{rid}.rv_coord_transformation.arrow'))
			iecf.assign_is_id_to_interval_bound(depth_interval_df, r_coor_tran_rv, is_pos_conv_vec_ref)
			depth_interval_df.to_csv(os.path.join(outdir, args.prefix + '_genome.'+qid+'.depth_interval.csv'), index=False)
			logging.info(f'Finished calculating depth interval from {bam_dir}. Saved to {os.path.join(outdir, args.prefix + "_genome."+qid+".depth_interval.csv")}')
		else :
			logging.warning(f'No bam file {bam_dir} found. Skipping depth interval calculation.')

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
	logging.info("Start analyzing using SyRI_IS output")
	logging.info("CMD: " + " ".join(sys.argv))

if __name__ == '__main__':
	main()