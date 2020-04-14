#!/usr/bin/env python 


#This script processes all indel.csvs generated with indels_w_difflib_v2.3.py in the current directory. It outputs a single csv
#with pertinent data for all samples/indel.csvs in directory.

#Usage: python /Users/Jack/Desktop/Dropbox/maly_lab/Cas9/Cas9-CRISPR/scripts/process_indel_data.py

import sys, csv, os, glob

# def getFileList():
# 	l = glob.glob("*.csv")


def cat(filename):
	f = open(filename, 'rU')
	indelinfo = f.read().splitlines()
	indelinfo = [s.split('\t') for s in indelinfo]
	return indelinfo



def main():
	# with open(str('indel_analysis.csv','w') as out:
 #  		csv_out = csv.writer(out)

 	with open('indel_analysis.csv','w') as out:
  		csv_out = csv.writer(out)
		csv_out.writerow(['Sample','Indel (stringent) %','Indel %','Reads filered (low Phred)','Reads lacking nmer(s)','Reads passing filters','Total Indels (stringent)','Total Indels'])
		l = glob.glob("*difflib_indels.csv")
		for i in l:
			print i
			indel_stats = cat(i)

			indel_stringent = str(indel_stats[12])
			indel_stringent = indel_stringent.split(",")
			indel_stringent = str(indel_stringent[1])
			indel_stringent = indel_stringent.split("'")
			indel_stringent = float(indel_stringent[0])

			indel_nonstri = str(indel_stats[10])
			indel_nonstri = indel_nonstri.split(",")
			indel_nonstri = str(indel_nonstri[1])
			indel_nonstri = indel_nonstri.split("'")
			indel_nonstri = float(indel_nonstri[0])

			phred_filtered = str(indel_stats[6])
			phred_filtered = phred_filtered.split(",")
			phred_filtered = str(phred_filtered[1])
			phred_filtered = phred_filtered.split("'")
			phred_filtered = int(phred_filtered[0])

			lack_nmer = str(indel_stats[7])
			lack_nmer = lack_nmer.split(",")
			lack_nmer = str(lack_nmer[1])
			lack_nmer = lack_nmer.split("'")
			lack_nmer = int(lack_nmer[0])

			good_reads = str(indel_stats[8])
			good_reads = good_reads.split(",")
			good_reads = str(good_reads[1])
			good_reads = good_reads.split("'")
			good_reads = int(good_reads[0])

			tot_indels = str(indel_stats[9])
			tot_indels = tot_indels.split(",")
			tot_indels = str(tot_indels[1])
			tot_indels = tot_indels.split("'")
			tot_indels = int(tot_indels[0])

			tot_indels_stringent = str(indel_stats[11])
			tot_indels_stringent = tot_indels_stringent.split(",")
			tot_indels_stringent = str(tot_indels_stringent[1])
			tot_indels_stringent = tot_indels_stringent.split("'")
			tot_indels_stringent = int(tot_indels_stringent[0])




			csv_out.writerow([str(i),indel_stringent,indel_nonstri,phred_filtered,lack_nmer,good_reads,tot_indels_stringent,tot_indels])



if __name__ == '__main__':
  main()
