import sys
import re
import numpy as np


input_file=sys.argv[1]
output_file=sys.argv[2]

# qseqid

with open(input_file, "r") as input, open(output_file, "w") as output:
   for line in input:      
      string = '^qseqid'
      if not re.search(string,line) :
        line=line.strip()
        info = line.split("\t")
        if info[-1] == "Viruses":
           new= info[0]+"\n"
           output.write(new)
        
