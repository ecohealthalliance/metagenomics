import sys
import re
import numpy as np


input_file=sys.argv[1]
output_file=sys.argv[2]

# qseqid

with open(input_file, "r") as input, open(output_file, "w") as output:
   for line in input:      
      line=line.strip()
      string = '^Genome'
      if not re.search(string,line) :
        info = line.split("\t")
        if info[-2] == "Viruses":
           new=line+"\n"
           output.write(new)
      else:
        new=line+"\n"
        output.write(new)            
 
