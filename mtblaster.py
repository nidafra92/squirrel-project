#!/usr/bin/env python
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
import sys

max_length = 17000


readsfq_filename = sys.argv[1]
database_filename = sys.argv[2]


filtreads_list = []

print("Filtering reads")
with open(readsfq_filename, "r") as read_file, open("filtered_reads.fasta", "w") as filtered_file:
	for read in SeqIO.parse(read_file, "fastq"):
		read_len = len(read.seq)
		if read_len <= max_length:
			filtreads_list.append(read)
			# 			read_tag = read.id.split()[0]
# 			kmer_seq = read.seq[pos:pos + k]
# 			kmer_name = read_tag + "_kmer" + str(pos)
# 			filtered_file.write(">" + kmer_name + "\n")
# 			filtered_file.write(str(kmer_seq) + "\n")
	print(len(filtreads_list))
	SeqIO.write(filtreads_list, filtered_file, "fasta")

print("Blasting valid reads against local reference database")

cline = NcbiblastnCommandline(query="filtered_reads.fasta", subject=database_filename, \
	evalue=0.01, out="blast_results.xml", outfmt=5, \
	#word_size=15, gapopen=3, gapextend=2, reward=2, num_threads=4, perc_identity=70, word_size=15)
	task="blastn", num_threads=4)

stdout, stderr = cline()

print("Parsing BLAST results")

reads_w_hits = set([])
total_hits = 0

with open("blast_results.xml") as result_handle:
	blast_records = NCBIXML.parse(result_handle)
	for blast_record in blast_records:
		query_len = blast_record.query_length
		for alignment in blast_record.alignments:
			total_qhitlenght = 0
			for hsp in alignment.hsps:
				total_qhitlenght += len(hsp.query)
				total_hits += 1
					#print('sequence:', alignment.title)
					#print('length:', alignment.length)
					#print('e value:', hsp.expect)
			print(str(total_qhitlenght) + " " + str(query_len))
			#if total_qhitlenght >= query_len*0.5:
			if total_qhitlenght >= 100:
				reads_w_hits.add(str(blast_record.query).split()[0])

print("Total of %s" % str(total_hits))
print("belonging to %s reads" % str(len(reads_w_hits)))


record_list = []

with open(readsfq_filename) as reads:
	for read in SeqIO.parse(reads, "fastq"):
		read_name = read.id.split()[0]
		if read_name in reads_w_hits:
			record_list.append(read)

print("Writing to file")

with open("reads_w_hits.fastq", "w") as output_handle:
	SeqIO.write(record_list, output_handle, "fastq")