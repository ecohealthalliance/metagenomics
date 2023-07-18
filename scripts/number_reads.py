import sys
import re
from Bio import SeqIO

input_file=sys.argv[1]
file_reads=sys.argv[2]


reads = {}

with open(file_reads, 'r') as lines:
     for line in lines:   
         line=line.strip()
         string = '^Genome'
         if not re.search(string,line) :         
           info = line.split('\t')
           nre = info[0].replace(".contigs.supported.csv", "")
           reads[nre] = info[2]
         
#LC316817.1	2	892	7355	0.1212780	83.333-87.978	Enterovirus G genomic RNA, nearly complete genome, strain: EVG/Porcine/JPN/Iba464-3-2/2015/G4	106966	Enterovirus G	Enterovirus	Picornaviridae	Picornavirales	Pisoniviricetes	Pisuviricota	Viruses	k141_11611_9,k141_108464_9	NM


with open(input_file, 'r') as lines:
     for line in lines:
      line=line.strip()
      string = '^Genome'
      if not re.search(string,line) :
        info = line.split('\t')
        length = info[3]
        pre_ct = info[15]
        taxaid = info[7]
        total=0
        contigs = pre_ct.split(',')
        for ct in contigs:
            if ct in reads:
                 predt = ct.split("_")
                 sample= predt[-2]
                 inter = total+int(reads[ct])
                 total = inter
            else:
                 print("error")
                 print(line)
        div=int(total)/int(length)
        imp = [  info[0], sample, taxaid, total, length, div ]
        pre_st = [  str(i) for i in imp  ]
        string =  "\t".join(pre_st)
        print(string)
