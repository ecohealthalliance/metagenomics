import sys
import re

## Goals 
## build bed files for the target sequence 

input_file=sys.argv[1]
output_file=sys.argv[2]

## arguments
## input_file: report by contig

control = {}

with open(input_file, 'r') as lines, open(output_file, 'w') as  out_line:
    for line in lines:
     string='^qseqid'
     if not re.search(string,line) :
       line=line.strip()
       info=line.split('\t')
       target = info[2]
       if not target in control:
           control[target] = ""
           length = info[3]
           pos_fn = int(length) - 1
           imp = [ target, 0, pos_fn  ]
           pre= [ str(i) for i in imp ]
           string ="\t".join(pre)
           string+="\n"
           out_line.write(string)
       
