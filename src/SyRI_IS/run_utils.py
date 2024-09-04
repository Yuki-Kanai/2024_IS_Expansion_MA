import os
import logging

from Bio import SeqIO
from Bio.Seq import Seq

# Use ragtag to connect contigs based on reference
# if more than one contig in the input fasta
def get_full_contig(filename, tmpdir='tmpdir', basefasta = '', skip_if_exists = True, print_cmd = False):
	full_contig_dir = filename
	base_name = os.path.basename(filename)
	# count contig number of filename fasta
	contig_cnt = len([record for record in SeqIO.parse(filename, 'fasta')])
	print('contig_cnt', contig_cnt)
	if contig_cnt > 1: # if more than one contig, run ragtag
		base_name = os.path.splitext(base_name)[0]
		full_contig_folder = os.path.join(tmpdir, 'ragtag', base_name)
		full_contig_dir = os.path.join(full_contig_folder, base_name + '.fasta')
		if skip_if_exists and os.path.isdir(full_contig_folder):
			message = 'Skipping ragtag for ' + base_name + ' because ' + full_contig_folder + ' exists.'
			print(message) if print_cmd else logging.info(message)
		else:
			os.makedirs(full_contig_folder, exist_ok=True)
			C = 'ragtag.py scaffold ' + basefasta + ' ' + filename + ' -o ' + full_contig_folder
			print(C) if print_cmd else logging.info(C)
			os.system(C)
			C2 = 'cp ' + os.path.join(full_contig_folder, 'ragtag.scaffold.fasta') + ' ' + full_contig_dir
			print(C2) if print_cmd else logging.info(C2)
			os.system(C2)
	return full_contig_dir