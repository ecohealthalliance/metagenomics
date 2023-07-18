import sys
import numpy as np
import re

## Goal 
## Identify all sequences with annotations shorter than normal. Send these sequences to files called constructs and remove these sequences from the files in the taxa directory.


input_file=sys.argv[1]
input_multiple=sys.argv[2]
output_clean_sg=sys.argv[3]
output_clean_mt=sys.argv[4]
output_constructs=sys.argv[5]
output_count=sys.argv[6]


#input_file:file from taxa directory 
#output_clean: file without sequences with  incompleta anotación taxonómica
#output_constructs: file with sequences with  incompleta anotación taxonómica
#output_count: count number of sequences with only best hit or with multiple best hits and number of sequences that are constructs 


Head = [ "qseqid" ,	"sseqid"	,"pident"	,"length"	,"mismatch", "gapopen"	,"qstart"    ,"qend"	,"sstart"	,"send"	,"evalue"	,"bitscore",	"qcovs" , "slen", "sgi"	,"sscinames"	,"stitle" , "staxid", "species","genus","family","order","class","phylum","superkingdom"]

string_head  = "\t".join( str(i) for i in Head )
string_head +="\n" 

constr = {}
multiple = {}
Total= {}

with open(input_file, 'r') as handled:
         for line in handled:
            line=line.strip()
            string = '^qseqid'
            if not re.search(string,line) :
               info = line.split('\t')
               Total[info[0]]= ""
               lg = len(info)
               if lg != 25:
                if info[0] in constr:      
                  constr[info[0]][line] = ""
                else:
                  constr[info[0]] = {}
                  constr[info[0]][line] = "" 


with open(input_multiple, 'r') as handled:
         for line in handled:
            line=line.strip()
            string = '^qseqid'
            if not re.search(string,line) :
             info = line.split('\t')
             Total[info[0]]= ""
             lg = len(info)
             if lg != 25:
               if info[0] in constr:      
                  constr[info[0]][line] = ""
               else:
                  constr[info[0]] = {}
                  constr[info[0]][line] = "" 


with  open(input_file, 'r') as handled  , open(output_clean_sg, 'w') as out_line:
         out_line.write(string_head) 
         for line in  handled:
            line=line.strip()
            string = '^qseqid'
            if not re.search(string,line) :
             info =  line.split('\t')
             if not info[0] in  constr:
                 post = [ 'unknown' if  i == ""  else i for i in info]
                 ps_str = "\t".join( str(i) for i in post )
                 ps_str+="\n" 
                 out_line.write(ps_str)


with open(output_clean_mt, 'w') as out_line, open(input_multiple, 'r') as handled:
         out_line.write(string_head) 
         for line in handled:
            line=line.strip()
            string = '^qseqid'
            if not re.search(string,line) :
             info = line.split('\t')
             if not info[0] in  constr:
                 post = [ 'unknown' if  i == ""  else i for i in info]
                 ps_str = "\t".join( str(i) for i in post )
                 ps_str+="\n" 
                 out_line.write(ps_str)       
                 if not info[0] in multiple:      
                          multiple[info[0]]= {}       
                          multiple[info[0]][info[1]] = ""
                 else:
                          multiple[info[0]][info[1]] = ""      


with open(output_constructs, 'w') as out_const:
         for ct, lines in constr.items():
           for line in lines:
               line+="\n"
               out_const.write(line)


cons_ct = len(constr.keys())            

# head = [ "Total", "Contigs_one_best_hit", "Contigs_multiple_best_hot", "Plasmid_sequence"  ]
# str_he = "\t".join(str(i) for i in head)
# str_he +="\n"
#out_ct.write(str_he) 

name=input_file.replace("taxa/", "")
name=name.replace(".best.taxa.csv", "")

with open(output_count, 'w') as out_ct:
         ct_single = 0
         ct_mul = 0
         total = len(Total.keys())  
         for ct, lines in multiple.items(): 
            lg = len(lines.keys())
            if lg == 1:
               ct_single += 1
            else:
                ct_mul += 1      
         pre = [name, total, ct_single,  ct_mul, cons_ct ]
         string ="\t".join(str(i) for i in pre)
         string +="\n"
         out_ct.write(string) 
