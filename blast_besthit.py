#!/usr/local/bin/python3

import sys, os, argparse
import pandas as pd

"""blast_besthit.py: pulling out only the best hit from blast results (in outfmt 6)"""

__author__ = "Sarah Stevens"
__email__ = "sstevens2@wisc.edu"
	
# Inputs
def parseArgs():
	parser = argparse.ArgumentParser(description='blast_besthit.py: pulls out only the best hit from blast results(outfmt 6), option for pooled file types, also can keep only one or all(if equal bit scores)')
	parser.add_argument('--blast_in','-bin' , action="store", dest='blast_in', type=str, required=True, metavar='catblastoutfm6', help="This should be the blast output file in outformat 6")
	parser.add_argument('--pooled','-p', action="store_true", dest='pool', default=False, help='If you are using a pooled(custom file type) blast output file, add this flag')
	parser.add_argument('--combBBH', '-cBBH', action='store_true', dest='cbbh', default=False, help='If you are using a pooled file that is combined bbh files, you need this flag')
	parser.add_argument('--keepAllBest', '-kb', action='store_true', dest='keepboth', default=False,help='If you want to keep all of the best hits, include this flag.  Otherwise, it only keeps the first one')
	args=parser.parse_args()
	return args.blast_in, args.pool, args.keepboth, args.cbbh

# Read in input
blastfile, pooled, keepBoth, cbbh = parseArgs()
if pooled:
	if not cbbh:
		inblast = pd.read_table(blastfile,delim_whitespace=True, header=None,names=['pool','read_info','subject','PID','align_len','mismatches','gaps','q_start','q_end','s_start','s_end','evalue','bit_score'])
		inblast['read']=inblast['read_info'].str.split('.blast:').str.get(1)
	else:
		inblast = pd.read_table(blastfile,delim_whitespace=True, header=None,names=['pool','read_info','subject','PID','align_len','mismatches','gaps','q_start','q_end','s_start','s_end','evalue','bit_score','read'])
else:
	inblast = pd.read_table(blastfile,delim_whitespace=True,header=None,names=['read','subject','PID','align_len','mismatches','gaps','q_start','q_end','s_start','s_end','evalue','bit_score'])

# Remove duplicates
bs_maxes = inblast.groupby('read').bit_score.transform(max) # finds maxes
bbh_df = inblast[inblast.bit_score == bs_maxes] # keeps only the max values
if not keepBoth:
	bbh_df = bbh_df.drop_duplicates('read',keep='first') # removes any duplicates (if both are max) - only executed if -keepAllBest flag not included
	# Checking that there are no duplicate reads left
	checkingdups = bbh_df[bbh_df.duplicated('read') == True] # creates a new dataframe with the repeated rows
	assert len(checkingdups) == 0 #shouldn't have any duplicates left and pass this
print('{}\tStarted with {} hits, kept {} hits'.format(blastfile,str(len(inblast)), str(len(bbh_df)))) # prints info about number of hits before and after

# Output to file
if keepBoth:
	bbh_df.to_csv(blastfile+'.keepAll.bbh', sep='\t', header=False, index=False)
else:
	bbh_df.to_csv(blastfile+'.bbh', sep='\t', header=False, index=False)
