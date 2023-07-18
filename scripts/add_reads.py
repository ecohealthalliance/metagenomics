import sys
import re

## Goal 
## new id: add Id sample to Id sequence

input_file=sys.argv[1]

acc = {}

# KP118894.1	S10	11120	1591	27684	0.057470018783412805

with open(input_file, "r") as input:
   for line in input:      
        line=line.strip()
        info = line.split("\t")
        if not info[2] in acc:
          acc[info[2]]=int(info[3])
        else:
          inter = acc[info[2]]
          acc[info[2]] = inter + int(info[3])  

head = [ "readn", "contign", "taxid" ]
sup = [ str(i) for i in head ]
imp_hed = "\t".join(sup)
print(imp_hed)

for ky in acc.keys():
  pre = [ str(acc[ky]) , str(1), str(ky) ]
  imp = "\t".join(pre)
  print(imp) 
  
