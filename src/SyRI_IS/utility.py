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

######### FROM TELR

def rm_file(file):
	if os.path.exists(file):
		os.remove(file)


def mkdir(dir):
	if os.path.isdir(dir):
		print("Directory %s exists" % dir)
		return
	try:
		os.mkdir(dir)
	except OSError:
		print("Creation of the directory %s failed" % dir)
	else:
		print("Successfully created the directory %s " % dir)


def check_exist(file):
	if file:
		if os.path.isfile(file) and os.stat(file).st_size != 0:
			return True
		else:
			return False
	else:
		return False


def format_time(time):
	d = datetime(1, 1, 1) + timedelta(seconds=time)
	if d.hour == 0 and d.minute == 0:
		return "%d seconds" % (d.second)
	elif d.hour == 0 and d.minute != 0:
		return "%d minutes %d seconds" % (d.minute, d.second)
	else:
		return "%d hours %d minutes %d seconds" % (d.hour, d.minute, d.second)


def create_loci_set(vcf_parsed):
	all_loci = set()
	with open(vcf_parsed, "r") as input:
		for line in input:
			entry = line.replace("\n", "").split("\t")
			all_loci.add("_".join(entry[0:3]))
	return all_loci


def report_bad_loci(set_raw, set_filtered, report, message):
	with open(report, "a") as output:
		for locus in set_raw:
			if locus not in set_filtered:
				output.write("\t".join([locus, message]) + "\n")


def get_cmd_output(cmd_list):
	"""get output from subprocess"""
	output = subprocess.check_output(cmd_list)
	output = output.decode("utf-8")
	return output


def get_rev_comp_sequence(fasta_in, fasta_out):
	"""get reverse complement of a sequence"""
	with open(fasta_out, "w") as output:
		for record in SeqIO.parse(fasta_in, "fasta"):
			output.write(
				">" + record.id + "\n" + str(record.seq.reverse_complement()) + "\n"
			)


def export_env(file):
	"""export conda environment"""
	file_tmp = file + ".tmp"
	cmd_list = ["conda", "env", "export", "--name", "TELR", "--file", file_tmp]
	subprocess.call(cmd_list)
	with open(file, "w") as output, open(file_tmp, "r") as input:
		for line in input:
			if (
				not line.startswith("prefix:")
				and "- pip:" not in line
				and "- telr==" not in line
			):
				output.write(line)
	os.remove(file_tmp)
	########  end from telr