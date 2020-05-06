#!/bin/sh

#This shell script iteratively implements the indels_w_difflib_v2.2.py script across all fastq files in the current directory


# fastqfiles=$(ls *.fastq)

for i in $(ls *.fastq)
do
	python indels_w_difflib_find_ODN.py $1 $i $2 $3 $4 $5 $6 $7 $8
	echo $i

done



#Usage:
#sh /Users/Jack/Desktop/Dropbox/maly_lab/Cas9/Cas9-CRISPR/scripts/implement_indels_w_difflib.sh GGGATCGCGCTGAGTATAAAAGCCGGTTTTCGGGGCTTTATCTAACTCGCTGTAGTAATTCCAGCGAGAGGCAGAGGGAGCGAGCGGGCGGCCGGCTAGGGTGGAAGAGCCGGGCGAGCAGAGCTGCGCTGCGGGCGTCCTGGGAAGGGAGATCCGGAGCGAATAGGGGGCTTCGCCTCTGGCCCAGCCCTCCCGCTGATCCCCCAGCCAGCGGTCCGCAACCCTTGCCGCATCCACGAAACTTTGCCCA 30 GTAATTCCAGCGAGAGGCAG 20 8
#
#sh implement_indels_w_difflib.sh ref_amplicon phred_threshold target_seq flanking_length nmer_length
