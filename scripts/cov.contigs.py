import sys
import re

## Goal 
## Retrieve Ids of contigs that have alignments satisfying user-defined criteria

file=sys.argv[1]
cov=sys.argv[2]
num_reads=sys.argv[3]

cov_rf=float(cov)
num_rf=int(num_reads)

## arguments
## file: coverage output file bbmap
## cov: minimum percentage of contig length covered by reads
## num_reads: minimum number of reads covering the contigs 

## output
## file with contig ID, coverage and number os reads aligned.

with open(file, 'r') as f:
    for line in f:
        line=line.strip()
        string='^#.+'
        if not re.search(string,line) :
            list = line.split('\t')
            cov_per = float(list[4])
            total = int(list[6])+int(list[7])
            if cov_per >= cov_rf and total >= num_rf:
                id=list[0]
                pre = [id,cov_per,total]
                imp="\t".join(str(i) for i in pre)
                print(imp)
