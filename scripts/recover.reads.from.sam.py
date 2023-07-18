import sys
import re

sam=sys.argv[1]

list = {}

with open(sam, 'r') as f:
    for line in f:
        line=line.strip()
        str='^@.+'
        if not re.search(str,line) :
            id = line.split('\t')[0]
            if not id in list:
                list[id]="NA"


for keys in list:
    print(keys)
