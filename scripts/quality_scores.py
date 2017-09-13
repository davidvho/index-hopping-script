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
            NR += 1
            line = line.strip('\n')               # strip off \n            
            if NR == 4:                           # look at the forth line
                sequence_length = len(line)       # how many
                break
    with gzip.open(fastq,"rt") as file:       # open file
        NR = 0
       	reads = 0
        all_scores = np.zeros(sequence_length)    # running total of means
        score_freq = {}    # dictionary to keep avg quality score for each read
        
        for line in file:            
            NR += 1
            line = line.strip('\n')               # strip off \n            
            if NR%4==0:                           # for every 4th line                   
                scores = np.zeros(sequence_length)         # quality scores for that read
                i = 0
                reads += 1                         # add another count to the number of reads
                
                for bp in line:                    # Look at the quality score line and convert it the ascii
                    all_scores[i] += convert_phred(line[i])    # Add to running total
                    scores[i] += convert_phred(line[i])        # What are the scores for THIS read
                    i += 1  
                
                read_mean = sum(scores)/len(scores)        # for THIS read, what is the avg quality score
                
                if read_mean in score_freq:                # add count of avg quality read score to dictionary
                    score_freq[read_mean] += 1
                else:
                    score_freq[read_mean] = 1
#for the running total, what is the mean at each position?
        for i in range(len(all_scores)):
            all_scores[i] = all_scores[i]/reads
#write out two files, (1) avg score at each position, (2) dictionary of avg scores of every read where value = count
        with open(fastq+"_mean_score_bp.txt",'w')	as fh: # write out the metrics to a new file
            fh.write("Mean Quality Score"+"\n")	  #Print the header you want
            i = 0
            for x in all_scores:
                fh.write(str(all_scores[i])+"\n")
                i += 1
        with open(fastq+"_freq.txt", "w") as fh:
            fh.write("Average qual score & Frequency"+"\n")
            for key,value in score_freq.items():
                fh.write(str(key)+"\t"+str(value)+"\n")

print("Created file(s) with mean quality scores at each bp position & frequency of avg quality score/read")
