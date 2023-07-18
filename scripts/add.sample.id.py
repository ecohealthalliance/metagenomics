import sys
import re
from Bio import SeqIO

input_file=sys.argv[1]


desc = input_file.split('/')
sample = desc[1]

acc = {}

with open(input_file, 'r') as lines:
     for line in lines:      
        line=line.strip()
        info = line.split('\t')
        new = info[0]+"_"+sample
        info[0] = new
        pre = [ str(i) for i in info ]
        string = "\t".join(pre)
        print(string)


