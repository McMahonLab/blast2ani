#!/usr/bin/python

import sys, os, csv

"""poolBLASTS.py takes the metagenome metadata and a set of metagenome blasts, 
	and pools each month from the same year together."""

__author__ = "Sarah Stevens"
__email__ = "sstevens2@wisc.edu"

def usage():
	print "Usage: poolBLASTS.py blastfile metadatafile"
	sys.exit(2)

if len(sys.argv) != 3:
	usage()
	exit()

# Ins and Outs
blastfile=open(sys.argv[1], 'rU')
blast=blastfile.read()
output=open(sys.argv[1]+'.pooled', 'w')

monthyear=[]
tocat=[]
# open file
with open(sys.argv[2], 'rU') as metafile:
	# Read through each line of metadata
	metadata=csv.reader(metafile, delimiter='\t')
	for row in metadata:
		name=row[0]
		year=row[5]
		month=row[6].zfill(2)
		my=year+"_"+month # new names for pooled
		# Replace sample name with year_month
		blast=blast.replace(name,my)
		if my not in monthyear:
			monthyear.append(my)
			tocat.append([])
		index=monthyear.index(my)
		tocat[index].append(name)

output.write(blast)
output.close()

outfile=open('1metaMonYears.txt', 'w') # pools with only 1 sample in that month
outfile2=open('pools.txt','w') # file to match up the sample with the pool easily
for i,group in enumerate(tocat):
	if len(group) == 1:
		outfile.write(monthyear[i]+'\t'+group[0]+'\n')
	for sample in group:
		outfile2.write(sample+'\t'+monthyear[i]+'\n')
outfile.close()
outfile2.close()