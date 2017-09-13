##### IMPORTS

import argparse
import gzip

##### BACKGROUND FUNCTIONS

def convert_phred(letter):
    """Converts a single ASCII character into a quality score, assuming Phred+33"""
    qscore = ord(letter) - 33   # ord takes the ASCII char and turns it into a Phred score
    return qscore

def qual_score():
    parser = argparse.ArgumentParser(description='Looking at the index combos from paired-end reads. Input -R3 file barcode will be reverse complemented. Requires an index file. Output are statistics for the files & counts for expected index pairs AND counts for index pairs that are swapped. Input sequencing files need to be GZIPPED.')
    parser.add_argument("-R1", help="Foward read FASTQ file, .gz", required=True, type=str)
    parser.add_argument("-R2", help="1st barcode FASTQ file, .gz", required=True, type=str)
    parser.add_argument("-R3", help="2nd barcode FASTQ file, .gz, note: the barcodes will be reverse complemented.", required=True, type=str)
    parser.add_argument("-R4", help="Reverse read FASTQ file, .gz", required=True, type=str)
    parser.add_argument("-i", help="tab-delim file with indexes in second column", required=True, type=str)
    parser.add_argument("-c", help="If a barcode read contains a quality score below this limit, the whole read will be 'thrown out'", required=True, type=int)
    return parser.parse_args()

args = qual_score()

## In this case, since we are not outputting the files, we just need to look at R2/R3, but... let's just include R1/R4 as well.

def reverse_complement(sequence):     ## function to get the reverse complement of a sequence
    """Returns reverse complement of a nucleotide DNA sequence"""
    complement = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
    return "".join([complement.get(nt.upper(), '') for nt in sequence[::-1]])

##### MAIN

index_combos = {}   #create a dictionary to store expected, good index pairs
index_list = []     #create a list to store indexes
all_pairs = {}      #create a dictionary to store index pairs that would occur if there was any index swapping

with open(args.i,"r") as index1:       #open the index file and put the indexes that are in the 2nd column to a list (index_list)
    for line in index1:
        line = line.split()
        index = line[1]
        index_list.append(line[1])     #from index file, store indices into index_list
    
    possible_index_pairs = [index1+"_"+index2 for index1 in index_list for index2 in index_list if index1 == index2] #what are possible GOOD combos?

all_combos = [index1+"_"+index2 for index1 in index_list for index2 in index_list if index1 != index2]  #what are possible BAD combos aka INDEX SWAPS

for pairs in possible_index_pairs:
    index_combos[pairs] = 0                #put possible GOOD combos into dictionary

for pairs in all_combos:     #put possible BAD combos/INDEX SWAP into dictionary
    all_pairs[pairs] = 0


LIMIT = args.c     #what is the limit/cutoff?
determined = 0     #how many reads are above cutoff, good, expected?
bad_pairs = 0      #how many reads are above cutoff but are index swapped?
undetermined = 0   #how many reads have Ns?
sequencing_error = 0   #how many are sequencing errros?
bad_quality = 0    #how many reads do not meet cutoff?
count = 0          #how many reads total?

def convert_phred(letter):
    """Converts a single ASCII character into a quality score, assuming Phred+33"""
    qscore = ord(letter) - 33   # ord takes the ASCII char and turns it into a Phred score
    return qscore

with gzip.open(args.R1, 'rt') as read1, gzip.open(args.R2, 'rt') as barcode1, gzip.open(args.R3, 'rt') as barcode2, gzip.open(args.R4,'rt') as read2:   #open up the files
    
    for line in read1:           #for all files above, read the first four lines of each and save as variables (16 variables in total)
        line = line.strip()
        read1_header = line.strip()
        read1_seq = read1.readline().strip()
        read1_plus = read1.readline().strip()
        read1_qual = read1.readline().strip()
        
        read2_header = read2.readline().strip()
        read2_seq = read2.readline().strip()
        read2_plus = read2.readline().strip()
        read2_qual = read2.readline().strip()
        
        barcode1_header = barcode1.readline().strip()
        barcode1_seq = barcode1.readline().strip()
        barcode1_plus = barcode1.readline().strip()
        barcode1_qual = barcode1.readline().strip()
        
        barcode2_header = barcode2.readline().strip()
        barcode2_seq = reverse_complement(barcode2.readline().strip())    #take the reverse complement of the sequence
        barcode2_plus = barcode2.readline().strip()
        barcode2_qual = barcode2.readline().strip()
        
        count += 1  #add to the running total of records
        
        index_pair = barcode1_seq+"_"+barcode2_seq   #create an index_pair variable
        #header_parts = read2_header.split()         #we want the header line if we were to write out new files
        #new_header = header_parts[0]+"_"+index_pair
        
        score = []  #list with the quality scores for both barcodes/read (len = 16)
        for bp in range(len(barcode1_qual)):        #look at the two barcode quality scores, convert, add to list score
            score.append(convert_phred(barcode1_qual[bp]))
            score.append(convert_phred(barcode2_qual[bp]))
        
        if all(items >= LIMIT for items in score):   #1. are all quality scores in the scores list >= LIMIT?
            if index_pair in index_combos:           #2. is it an expected index pair?
                index_combos[index_pair] += 1
                determined += 1
            
            if index_pair in all_pairs:    #2. is it an index pair that suggests swapping?
                all_pairs[index_pair] += 1
                bad_pairs += 1
            
            if "N" in index_pair:       #2. does it contain an N?
                undetermined += 1
            
            if index_pair not in index_combos and index_pair not in all_pairs and "N" not in index_pair:
                sequencing_error += 1      #if not good index pair, not bad index pair, and contains no Ns, it is therefore a sequencing error
    
        else:
            bad_quality += 1   #one or more of the quialty score of the barcodes is below the cutoff

#### Statistics to print at end of script
print("Cutoff:","\t",LIMIT)
print("Number of reads:","\t",count)
print("Number of retained, EXPECTED index pairs:","\t",determined)
print("Percentage of retained, EXPECTED index pairs out of ALL reads:","\t",determined/count*100,"%")
print("Number of retained, BAD index pairs (AKA Index-swapping):","\t",bad_pairs)
print("Number of retained reads with Ns in the barcode:","\t",undetermined)
print("Number of retained reads that have index not matching list (sequencing errors):","\t",sequencing_error)
print("Number of retained reads with Ns or sequencing errors:","\t",undetermined+sequencing_error)
print("Number of reads below cutoff:","\t",bad_quality)
print()
for key,value in index_combos.items():    #print out the expected indexes and their counts
    print(key,"\t",value)
print()
for key,value in all_pairs.items():       #print out the swapped indexes and their counts
    print(key,"\t",value)
