import os
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

fasta_files = os.listdir(".")
fasta_files = [f for f in fasta_files if f.endswith(".fasta", ".fa", ".fna", ".faa")]
df = pd.DataFrame()
for fasta_file in fasta_files:
	seqs = [str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]
	df = df.append(pd.DataFrame({"fasta_file"": fasta_file,
                              "length": max([len(seq) for seq in seqs]),
                              "n_contigs": len(seqs)}, index=[0]))
df.to_csv("fasta_stats.csv", index=False)
