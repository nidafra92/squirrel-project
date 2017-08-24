#!/usr/bin/env python
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO

k = 24

with open("fake_reads.fasta", "r") as read_file, open("kmers.fasta", "w") as kmer_file:
	for read in SeqIO.parse(read_file, "fasta"):
		read_len = len(read.seq)
		read_tag = read.id.split()[0]
		for pos in range(read_len - k):
			kmer_seq = read.seq[pos:pos + k]
			kmer_name = read_tag + "_kmer" + str(pos)
			kmer_file.write(">" + kmer_name + "\n")
			kmer_file.write(str(kmer_seq) + "\n")


cline = NcbiblastnCommandline(query="fake_reads.fasta", subject="sciuridae.fasta", \
	evalue=0.001, out="blast_results.xml", outfmt=5, perc_identity=100, ungapped=True)

stdout, stderr = cline()

with open("blast_results.xml") as result_handle:
	blast_records = NCBIXML.parse(result_handle)
	for blast_record in blast_records:
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				if len(hsp.query) == k:
					print(blast_record.query)
					print('sequence:', alignment.title)
					print('length:', alignment.length)
					print('e value:', hsp.expect)