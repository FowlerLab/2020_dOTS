#!/usr/bin/env python 

#A modified version of the indels_w_difflib_v2.2.py script for identifying, collecting, and exporting trimmed
#reads identied as having indels. This will allow better evaluationg of possible artificats and seq problems

import sys, csv, difflib

# sys.argv[1] = reference seq
# sys.argv[2] = input fastq
# sys.argv[3] = minimum avg phred score
# sys.argv[4] = target
# sys.argv[5] = flanking region length (bases between target and nmers)
# sys.argv[6] = nmer length (nmers are the sequence handles which form the new 5' and 3' ends of the ref and reads)
# sys.argv[7] = PAM position. + (PAM on right, 3' of target seq), - (PAM on left, 5' of target seq)
# Example usage:
# $ python /Users/Jack/Desktop/Dropbox/maly_lab/Cas9/Cas9-CRISPR/workarea_for_new_scripts/indels_w_difflib_v2.2.py GAGAAGGGCAGGGCTTCTCAGAGGCTTGGCGGGAAAAAGAACGGAGGGAGGGATCGCGCTGAGTATAAAAGCCGGTTTTCGGGGCTTTATCTAACTCGCTGTAGTAATTCCAGCGAGAGGCAGAGGGAGCGAGCGGGCGGCCGGCTAGGGTGGAAGAGCCGGGCGAGCAGAGCTGCGCTGCGGGCGTCCTGGGAAGGGAGATCCGGAGCGAATAGGGG 20-reads_5-deletions_5-insertions.fq 20 GTAATTCCAGCGAGAGGCAG 15 7

# Function to read in fastq file reads
def cat(filename,minscore,reference, target, flanklength, nmer_len):
	f = open(filename, 'rU')
	fastq = f.read().splitlines()
	fastq = [s.split('\t') for s in fastq]
	reads = []
	n = 1
	reads_filtered = 0
	for line in fastq:
		#only read in sequence lines, ignore names and quality
		if n % 4 == 2:
			reads.append(line)


		# If average read score is below threshold, remove preceding read from reads to be returned
		if n % 4 == 0:
			totscore = 0
			# print line
			for i in str(line):
				totscore = totscore + ord(i) - 33

			avgscore = float(totscore)/len(str(line))
			# print avgscore

			#remove last read if average score below threshold
			if avgscore < int(minscore):
				reads.pop(-1)
				reads_filtered = reads_filtered+1
				# print "deleted line " + str(n-2)
		n = n + 1

	# Initialize list to hold trimmed reads
	trimmed_reads = []

	# Initialize counter for reads removed because they do not have one or both nmers
	nmer_not_found_count = 0
	for line in reads:
			trimmed_read = trim_to_target(reference,target,flanklength,str(line),nmer_len)

			# Trimmed reads without one or both nmers are set to str()
			# If length of trimmed read < 1, increase nmer_not_found_count and do not add trimmed read to trimmed_reads list
			if len(trimmed_read) < 1:
				nmer_not_found_count = nmer_not_found_count + 1

			# Add trimmed_read to trimmed_reads list as long as length > 1
			else:
				trimmed_reads.append(trimmed_read)
	# print trimmed_reads
	# print "Reads removed which did not pass filtering: " + str(reads_filtered)
	# print "Reads lacking left or right nmer: "+str(nmer_not_found_count)

	return (trimmed_reads, reads_filtered, nmer_not_found_count)

def trim_to_target(reference, target, flanklength, read, nmer_len):


	# Find 5' index of target in reference
	target_index1 = reference.find(target)
	# print "Target_index: " + str(target_index1)
	# Add length of target to 5' index to get 3' index of target
	target_index2 = int(target_index1) + len(target)

	# Determine position of left nmer
	left_index1 = int(target_index1 - flanklength - nmer_len)
	left_index2 = int(target_index1 - flanklength)

	# Determine position of right nmer
	right_index1 = int(target_index2 + flanklength)
	right_index2 = int(target_index2 + flanklength + nmer_len)

	#Get left and right nmer sequences
	left_nmer = reference[left_index1:left_index2]
	# print "Left nmer: " + left_nmer
	right_nmer = reference[right_index1:right_index2]
	# print "Right nmer: " + right_nmer

	left_nmer_read_index = read.find(left_nmer)
	if left_nmer_read_index < 0:
		# print "Left nmer not found"
		trimmed_read = str()

	right_nmer_read_index = read.find(right_nmer)
	if right_nmer_read_index < 0:
		# print "Right nmer not found"
		trimmed_read = str()
	# print left_nmer_read_index
	# print right_nmer_read_index

	trimmed_read = read[left_nmer_read_index:(right_nmer_read_index+nmer_len)]

	return trimmed_read


def main():

	target = str(sys.argv[4])
	flanklength = int(sys.argv[5])
	nmer_len = int(sys.argv[6])

	reads = cat(sys.argv[2],sys.argv[3],sys.argv[1],target,flanklength,nmer_len)[0]
	filtered_reads = cat(sys.argv[2],sys.argv[3],sys.argv[1],target,flanklength,nmer_len)[1]
	lacking_nmer = cat(sys.argv[2],sys.argv[3],sys.argv[1],target,flanklength,nmer_len)[2]
	print "Reads removed which did not pass filtering: " + str(filtered_reads)
	print "Reads lacking left or right nmer: "+str(lacking_nmer)

	# Indel status stores indel status for each read 3-tuples: (indelstatus, # of dels, # of insertions)
	# If indelstatus = 1, at least one indel was detected in that read



	ref = trim_to_target(sys.argv[1],target,flanklength,sys.argv[1],nmer_len)
	print "Trimmed reference: " + ref
	if str(sys.argv[7]) == "-":
		cutsite_adjustment = 3
	elif str(sys.argv[7]) == "+":
		cutsite_adjustment = 17
	else:
		print "error: PAM orientation not specified. Please specify as + for left, - for right. "
	cutsite = int(sys.argv[5])+int(sys.argv[6])+cutsite_adjustment
	print "Cut site: "+ str(cutsite)

	cut_proximity = 1 #max distance from cutsite to be called indel
	indelstatus = []

	# Iteratively compare each read to reference and parse difflib operations to determine indel status
	for j in reads:
		# print j
		# trimmed_read = trim_to_target(sys.argv[1],target,flanklength,str(j),nmer_len)
		# print trimmed_read
		# compare reference and read using difflib
		s = difflib.SequenceMatcher(None, ref, str(j),autojunk=False)
		opers = s.get_opcodes() #get operations to convert ref to read

		# Initialize indel counters
		insert = 0
		dels = 0
		ind = 0 #0 if no indel, 1 if indel present
		ind2 = 0 #0 if no indel, 1 if indel present and len != ref
		indelreads = [] #list of reads with indels


		# Parse sequencematcher operations to count insertions and deletions in read
		# If 1 or more indel operation detected, set ind = 1
		for i in opers:

			if i[0] == "insert":
				#Ignore insertions within 5 bp of 5' end and 5bp of 3' of reference
				if (i[1] > (cutsite - cut_proximity) and (i[1] < (cutsite + cut_proximity))):
				# if (i[1] > 74) and (i[2] < (124)):
					insert = insert + 1
			if i[0] == "delete":

				if (i[1] > (cutsite - cut_proximity) and (i[1] < (cutsite + cut_proximity))):
					# If deletion detected, increase del counter by one
					dels = dels + 1
				elif (i[2] > (cutsite - cut_proximity) and (i[2] < (cutsite + cut_proximity))):
					dels = dels + 1
				elif (i[1] < cutsite < i[2]):
					dels = dels + 1


			if (dels + insert) > 0:
				#If indel detected, set ind to 1
				ind = 1
				if len(str(j)) != len(str(ref)):
					ind2 = 1


		indelstatus.append((ind, dels, insert, len(j), ind2, opers,str(j)))
		# if ind == 1:
		# 	print str(j)
		# 	indelreads.append(str(j))

	indlist = []
	ind2list = []

	# Extract indel status for purposes of counting indels with sum and dividing by total reads
	for k in indelstatus:
		indlist.append(k[0])
	for l in indelstatus:
		ind2list.append(l[4])

	print 'number of indels:'
	print sum(indlist)
	print 'number of indels with length != ref:'
	print sum(ind2list)
	print 'total reads:'
	print len(indelstatus)

	indperc = float(sum(indlist))/len(indelstatus)
	print 'indel %: '
	print indperc*100

	ind2perc = float(sum(ind2list))/len(indelstatus)
	print 'indel % (stringent): '
	print ind2perc*100
	#Write csv output
	with open(sys.argv[2]+'wtf-difflib_indels.csv','w') as out:
  		csv_out = csv.writer(out)
  		csv_out.writerow(['Ref seq:',sys.argv[1]])
  		csv_out.writerow(['Target seq:',sys.argv[4]])
  		csv_out.writerow(['Input:',sys.argv[2]])
  		csv_out.writerow(['Minimum avg Phred:',sys.argv[3],"Max indel distance to cutsite: ",str(cut_proximity)])
  		csv_out.writerow(['Flanking region length:',sys.argv[5]])
  		csv_out.writerow(['nmer length:',sys.argv[6]])
  		csv_out.writerow(['Reads filtered (low Phred): ', filtered_reads])
  		csv_out.writerow(['Reads missing lacking nmer(s): ',lacking_nmer])
  		csv_out.writerow(['Total reads (passed filters): ',len(indelstatus)])
  		csv_out.writerow(['Total indels: ',sum(indlist)])
  		csv_out.writerow(['Indel %: ',indperc*100])
  		csv_out.writerow(['Total indels (stringent): ',sum(ind2list)])
  		csv_out.writerow(['Indel (stringent) %: ',ind2perc*100])
  		csv_out.writerow(['Indel?','deletions','insertions','length','stringent indel?','diff_lib_operations','trimmed sequence'])
  		for row in indelstatus:
  			csv_out.writerow(row)

  	# print indelreads
  	# with open(sys.argv[2]+'indel_reads_list.csv','w') as out:
  	# 	csv2_out = csv.writer(out)
  	#  	for row in indelreads:
  	# 		csv2_out.writerow(row)



if __name__ == '__main__':
  main()
