# Index hopping and quality check

The contents of this repo was part of an assignment for a course.

Index hopping is a result of barcodes being switched in Illumina sequencing, causing \\
the mis-identification of libraries. This problem is seen more when patterned flow cells are used. \\
The python script in `scripts` takes in the four fastq files from sequencing and a list of unique-dual-indexed used.\\

`demultiplex_counts.py`: The output is the frequencies of sequences with expected index pairs AND index pairs that are swapped.\\

`quality_scores.py`: Takes in multiple FASTQ files and returns metrics for the quality score at each bp position. Output is a tab-delim .txt file per FASTQ input

