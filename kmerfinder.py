#!/usr/bin/env python3
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO

k = 30

reads_dict = {}

print("Generating k-mers from reads FASTA file")
with open("fake_reads.fasta", "r") as read_file, open("kmers.fasta", "w") as kmer_file:
	for read in SeqIO.parse(read_file, "fasta"):
		read_len = len(read.seq)
		read_tag = read.id.split()[0]
		if read_tag not in reads_dict:
			reads_dict[read_tag] = [read_len - k  + 1,0]

		for pos in range(read_len - k):
			kmer_seq = read.seq[pos:pos + k]
			kmer_name = read_tag + "_kmer" + str(pos)
			kmer_file.write(">" + kmer_name + "\n")
			kmer_file.write(str(kmer_seq) + "\n")

print("Generating concatenated database")

concatamer = ""

hits_list = []



with open("sciuridae.fasta") as database:
	for record in SeqIO.parse(database, "fasta"):
		record_seq = str(record.seq)
		first_part = record_seq[0:6000]
		concatamer += str(record.seq) + first_part + "X"*100


print("Searching k-mers against local reference database")


with open("kmers.fasta") as kmer_file:
	for kmer in SeqIO.parse(kmer_file, "fasta"):
		if str(kmer.seq) in concatamer:
			#print("Está :D")
			hits_list.append(kmer.id)
			read_length = "_".join(str(kmer.id).split("_")[0:-1])
			reads_dict[read_length][1] += 1

print("Total of %s" % str(len(hits_list)))
print("Total of reads %s" % str(len(reads_dict)))
#print(reads_dict)

acceptance_threshold = 0.75

candidate_reads = []

for read, hit_list in reads_dict.items():
	total_kmers = hit_list[0]
	total_hits = hit_list[1]
	hit_ratio =  float(total_hits)/total_kmers
	if hit_ratio > acceptance_threshold:
		candidate_reads.append(read)

print(len(candidate_reads))