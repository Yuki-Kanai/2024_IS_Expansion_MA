import logging
import os
import sys
import argparse
import glob
import re
import SyRI_IS

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

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
	logging.info("Start running SyRI_IS")
	logging.info("CMD: " + " ".join(sys.argv))

def get_args():
	parser = argparse.ArgumentParser(description='Plot multiple runs of SyRI_IS')
	parser.add_argument('dir', type=str, help='tmp dir of syri runs')
	parser.add_argument('single_out', type=str, help='single syri output file')
	parser.add_argument('-p', '--prefix', type=str, required=True, action='append', help='prefixes of syri runs')
	parser.add_argument('-i', '--id', type=str, required=True, action='append', help='genome identifiers of syri runs')

	#parser.add_argument('--syri_path', type=str, default='syri')
	parser.add_argument('--out', type=str, default='export/merged')
	parser.add_argument('--log', help = 'log directory', default='log')
	parser.add_argument('-t', '--tmp', type=str, default='tmp/merged')
	parser.add_argument('--plotsr', type=str, default='', help='additional options for plotsr')
	#parser.add_argument('--include', type=str, action='append', help='include only these genome identifiers', default=[])
	#parser.add_argument('--exclude', type=str, action='append', help='exclude these genome identifiers', default=[])
	#parser.add_argument('--start_from', type=int, default=0, help='start from this genome identifier')

	return parser.parse_args()

def rewrite_syri_out(syri_out_file, syri_save_dir, prefix_name = 'genome', contig_name = 'contig_1', append = False):
	C = "sed 's/" + contig_name + "/" + prefix_name + "/g' " + syri_out_file
	if append:
		C += " >> " + syri_save_dir
	else:
		C += " > " + syri_save_dir 
	logging.info(C)
	os.system(C)


def rewrite_syri_out_per_run(syri_out_files, syri_save_dir, prefix_name = 'genome', contig_name = 'contig_1'):
	for syri_out_file in syri_out_files:
		rewrite_syri_out(syri_out_file, syri_save_dir, prefix_name = prefix_name, contig_name = contig_name)

def rewrite_syri_outs(base_tmp_dir, prefixes, syri_save_dir, ids):
	saved_files = []
	os.makedirs(syri_save_dir, exist_ok=True)
	max_id = len(ids) - 2
	for prefix in prefixes:
		syri_out_files = glob.glob(os.path.join(base_tmp_dir, prefix, 'is_deleted', 'syri', '*rv_coord.syri.out')) 
		for syri_out_file in syri_out_files:
			gen_identifier = os.path.basename(syri_out_file).split('.')[1] # 0 of 'L01-3.0.syri.out'
			if not gen_identifier.isdigit():  # skip cases like 0_1
				continue
			if int(gen_identifier) > max_id: # don't include unspecified ids
				continue
			save_dir = os.path.join(syri_save_dir, 'merged.' + gen_identifier + '.syri.out')
			append_to_old = save_dir in saved_files
			rewrite_syri_out(syri_out_file, save_dir, prefix_name = prefix, contig_name = 'contig_1', append=append_to_old)
			if not append_to_old:
				saved_files.append(save_dir)
	assert len(saved_files) == len(ids)-1, 'saved_files and ids should have the same length but got ' + str(saved_files) + ' and ' + str(ids)
	logging.info('syri output files merged and saved to ' + syri_save_dir)
	# reorder the syris from order found to order of genome identifiers
	# required because this is the order the syri outputs are pasted to plotsr
	gen_identifiers = [int(os.path.basename(syri_out_file).split('.')[-3]) for syri_out_file in saved_files]
	saved_files = np.array(saved_files)[np.argsort(gen_identifiers)]
	return saved_files

def merge_fasta_for_merged_plotsr_(fastafiles, contig_names, save_dir):
	assert len(fastafiles) == len(contig_names), 'fastafiles and contig_names should have the same length'
	os.makedirs(os.path.dirname(save_dir), exist_ok=True)
	# read fasta files
	records = []
	for i, fastafile in enumerate(fastafiles):
		seq = SeqIO.read(fastafile, 'fasta')
		seq.id = contig_names[i]
		records.append(seq)
	# write to a new fasta file
	SeqIO.write(records, save_dir, 'fasta')
	return save_dir

def merge_fasta_for_merged_plotsr(base_tmp_dir, prefixes, save_dir, fasta_pattern = r".*_genome\.[0-9]+\.fasta", tmpdir = 'tmp/merged'):
	# pd with columns: genome_identifier, fastadir, prefix
	existing_fasta_files = []
	prefix_list = []
	for prefix in prefixes:
		fastafiles_ = []
		for fname in os.listdir(os.path.join(base_tmp_dir, prefix)):
			if re.match(fasta_pattern, fname):
				fastafiles_.append(os.path.join(base_tmp_dir, prefix, fname))
		prefix_list += [prefix] * len(fastafiles_)
		existing_fasta_files += fastafiles_
	assert len(existing_fasta_files) > 1, 'No fasta files found in ' + base_tmp_dir + ' with pattern ' + fasta_pattern
	# get genome_identifier
	df = pd.DataFrame(existing_fasta_files, columns=['fastadir'])
	df['genome_identifier'] = df.fastadir.apply(lambda x: os.path.basename(x).split('.')[-2])
	df['prefix'] = prefix_list
	df.to_csv(os.path.join(tmpdir, 'fasta_files.csv'), index=False)
	logging.info(str(len(df.index)) + ' fasta files found')

	# merge by genome_identifier
	ids = df['genome_identifier'].unique()
	ids.sort()
	merged_fasta_files = []
	for id in ids:
		df_ = df[df['genome_identifier'] == id]
		fastafiles_ = df_['fastadir'].tolist()
		prefixes_ = df_['prefix'].tolist()
		save_dir_ = os.path.join(save_dir, 'merged.' + id + '.fasta')
		merged_fasta_files += [save_dir_]
		merge_fasta_for_merged_plotsr_(fastafiles_, prefixes_, save_dir_)
		logging.info('fasta files merged and saved to ' + save_dir_)
	return merged_fasta_files

def merge_IS_tag_DF_for_merged_plotsr(export_dir, prefixes, genome_ids, save_dir, pattern='_IS_positions.csv'):
	csv_dirs = []
	for prefix in prefixes:
		csv_dirs_ = glob.glob(os.path.join(export_dir, 'is_position_files', prefix+pattern))
		if len(csv_dirs_) != 1:
			logging.warning('csv file not found for prefix ' + prefix)
			continue
		csv_dirs += [csv_dirs_[0]]
	for csv_dir in csv_dirs:
		if not os.path.exists(csv_dir):
			logging.warning('csv file not found for prefix ' + prefix)
			csv_dirs.remove(csv_dir)
	df = pd.concat([pd.read_csv(df_) for df_ in csv_dirs])
	max_id = len(genome_ids) - 1
	df = df.query('SyRI_genome_id <= ' + str(max_id)).reset_index(drop=True)
	logging.info('Detected ' + str(len(df.index)) + ' ISs. IS positions merged and saved to ' + save_dir)
	df.to_csv(save_dir, index=False)

	is_tag_df = pd.DataFrame(columns=['chr', 'start', 'end', 'genome_id', 'tags'])
	is_tag_df['start'] = df['IS_pos'].astype(int)-1
	is_tag_df['end'] = df['IS_pos'].astype(int)
	is_tag_df['genome_id'] = df.apply(lambda x: genome_ids[x.SyRI_genome_id], axis=1)
	is_tag_df['chr'] = df['SyRI_prefix']
	is_tag_df['tags'] = df.apply(lambda row: SyRI_IS.generate_IS_tag_for_plotsr(row.IS_strand), axis=1)
	is_tag_df.to_csv(save_dir, index=False, header=False, sep='\t')
	is_tag_df['tags'] = df.apply(lambda row: SyRI_IS.generate_IS_tag_for_plotsr(row.IS_strand, vertical=True), axis=1)
	is_tag_df.to_csv(save_dir+'.vertical', index=False, header=False, sep='\t')
	logging.info('IS tags merged and saved to ' + save_dir)
	return save_dir

def main():
	args = get_args()
	os.makedirs(args.log, exist_ok=True)
	set_logging(os.path.join(args.log, 'merged.log'))
	os.makedirs(args.tmp, exist_ok=True)
	os.makedirs(args.out, exist_ok=True)

	syri_out_files = rewrite_syri_outs(args.dir, args.prefix, os.path.join(args.tmp, 'syri'), args.id)
	syri_local_out_files = SyRI_IS.create_local_syri_output(syri_out_files)
	fasta_files = merge_fasta_for_merged_plotsr(args.dir, args.prefix, os.path.join(args.tmp, 'fasta'), tmpdir=args.tmp)
	genome_info = SyRI_IS.gen_genome_info_file_for_plotsr(fasta_files[0:(len(args.id))], args.id, formats=[], outputdir=os.path.join(args.tmp, 'genome_info.txt'))
	is_tag_dir = merge_IS_tag_DF_for_merged_plotsr(args.single_out, args.prefix, args.id, os.path.join(args.tmp, 'merged_IS_tag.bed'))
	savedir = SyRI_IS.draw_plot_with_plotsr_(syri_out_files, genome_info, prefix='merged', export_dir = args.out, markers = is_tag_dir, options=args.plotsr)
	savedir = SyRI_IS.draw_plot_with_plotsr_(syri_local_out_files, genome_info, prefix='merged-local', export_dir = args.out, markers = is_tag_dir, options=args.plotsr)
	#savedir = SyRI_IS.draw_plot_with_plotsr(syri_out_files, genome_info, prefix='merged', export_dir = args.out, markers = is_tag_dir+'.vertical', vertical = True)

if __name__ == "__main__":
	main()






