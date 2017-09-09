#!/usr/bin/env python

## Takes in multiple FASTQ files and returns metrics for the quality score at each bp position. Output is a tab-delim .txt file per FASTQ input

## IMPORTS

import numpy as np
import argparse
import zipfile
import gzip

## BACKGROUND DEFINITIONS

## Definte a function that gets a quality score
def convert_phred(letter):                                   
    """Converts a single ASCII character into a quality score, assuming Phred+33"""
    qscore = ord(letter) - 33   # ord takes the ASCII char and turns it into a Phred score 
    return qscore

## Arg parsing
def qual_score():
	parser = argparse.ArgumentParser(description='Fastq files past through will output a .txt file with metrics for the quality score at each bp position')
	parser.add_argument("-f", help="FASTQ file", required=True, type=str, nargs='*')
	return parser.parse_args()

args = qual_score()

## MAIN 

for fastq in args.f:
    with gzip.open(fastq,"rt") as file:      ## open a fastq file
        NR = 0
        for line in file:            
            NR += 4
            line = line.strip('\n')               # strip off \n            
            if NR == 4:                           # look at the forth line
                sequence_length = len(line)       # how many
                break
    with gzip.open(fastq,"rt") as file:       # open file
        NR = 0
       	reads = 0
        scores = np.zeros(sequence_length)
        for line in file:            
            NR += 1
            line = line.strip('\n')               # strip off \n            
            if NR%4==0:                           # for every 4th line                   
                i = 0
                reads += 1
                for bp in line:
                    scores[i] += convert_phred(line[i]) 
                    i += 1  

        for i in range(len(scores)):
            scores[i] = scores[i]/reads

        with open(fastq+"_q_score_metrics.txt",'w')	as fh: # write out the metrics to a new file
            fh.write("Mean Quality Score"+"\n")	  #Print the header you want
            i = 0
            for x in scores:
                fh.write(str(scores[i])+"\n")
                i += 1

print()
print("Created file(s) with mean quality scores at each bp position")
