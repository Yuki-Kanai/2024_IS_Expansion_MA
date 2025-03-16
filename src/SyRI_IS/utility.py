import os
from datetime import datetime, timedelta
import subprocess
from Bio import SeqIO


def get_full_contig(filename, tmpdir='tmpdir', basefasta = ''):
	# from prior analysis
	#assert basefasta != '', 'basefasta must be specified'
	#assert os.path.exists(basefasta), 'file does not exist'
	unconnected = ['L03-4_G8.fasta', 'L04-1_G8.fasta', 'L04-3_G8.fasta', 'L05-1_G8.fasta', 'L05-2_G8.fasta', 'L05-4_G8.fasta', 'L05_Anc.fasta', 'L07-4_G8.fasta', 'L10-2_G8.fasta', 'L11-4_G8.fasta']
	full_contig_dir = filename
	base_name = os.path.basename(filename)
	if base_name in unconnected:
		base_name = os.path.splitext(base_name)[0]
		full_contig_folder = os.path.join(tmpdir, 'ragtag', base_name)
		os.makedirs(full_contig_folder, exist_ok=True)
		C = 'ragtag.py scaffold ' + basefasta + ' ' + filename + ' -o ' + full_contig_folder
		os.system(C)
		full_contig_dir = os.path.join(full_contig_folder, base_name + '.fasta')
		C2 = 'cp ' + os.path.join(full_contig_folder, 'ragtag.scaffold.fasta') + ' ' + full_contig_dir
		print(C2)
		os.system(C2)
	return full_contig_dir


def format_time(time):
	d = datetime(1, 1, 1) + timedelta(seconds=time)
	if d.hour == 0 and d.minute == 0:
		return "%d seconds" % (d.second)
	elif d.hour == 0 and d.minute != 0:
		return "%d minutes %d seconds" % (d.minute, d.second)
	else:
		return "%d hours %d minutes %d seconds" % (d.hour, d.minute, d.second)

