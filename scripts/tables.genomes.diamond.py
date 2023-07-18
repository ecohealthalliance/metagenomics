import sys
import argparse
import re
import numpy as np

## Goal 

# Generate the report of the taxonomic classification of the contigs

parser = argparse.ArgumentParser()
parser.add_argument('--coverage_file', type=str, required=True)
parser.add_argument('--best_file', type=str, required=True)
parser.add_argument('--output_file', type=str, required=True)
args = parser.parse_args()


## arguments
## --coverage_file: output of bed software between the contigs alignments and the functional annotation of the genomes 
## --best_file: best-hit from the blast output file (add contigs id to table_genomes) 
## --output_file: for each target it is retrieved from the number of contigs, the bases covered, the length of the genome, and the percentage of the genome covered, the list of contig identifiers, the taxonomy and the type of sequence (genome, fragment)it from blast output aligments



Tg = {}
Tx = {}
It = {}

#qseqid	length	sseqid	slen	pident	qcovs	qstart	qend	sstart	send	evalue	bitscoresgi	sscinames	stitle	staxid	species	genus	family	order	class	phylum	superkingdom

with open(args.best_file, 'r') as lines:
   for line in lines:
       line=line.strip()
       string = '^qseqid'
       if not re.search(string,line) :
          line=line.split('\t')
          if not line[2] in Tg:
             Tg[line[2]] = {}
             Tg[line[2]][line[0]] = ""
          else:
             Tg[line[2]][line[0]] = ""
          if not line[2] in Tx:
             tx_list= [ line[14], line[15], line[16], line[17] , line[18], line[19], line[20], line[21], line[22]  ]
             tx_str= "\t".join(str(i) for i in tx_list)
             Tx[line[2]] = tx_str
          if not line[2] in It:  
             It[line[2]] = [float(line[4])]
          else:
             It[line[2]].append(float(line[4]))

             
## Take the information of the number of contigs, the bases covered, the length of the genome, and the percentage of the genome covered, add the list of the identifiers of the contigs, the taxonomy and the type of sequence (genome, fragment) 

   
with open(args.coverage_file, 'r') as lines, open(args.output_file, 'w') as line_out:
   hd=["Genome", "Number_of_contigs", "Covered_bases", "Length_genome","%_genome_covered", "Min_max_ident", "Stitle", "Staxid",	"Species"	,"Genus",	"Family"	,"Order"	, "Class",	"Phylum",	"Superkingdom"   ,"Contigs"]
   str_hd=  "\t".join(str(i) for i in hd)
   str_hd+="\n"
   line_out.write(str_hd)
   for line in lines:
      line=line.strip()
      line=line.split('\t')
      if  int(line[3]) > 0:
         if line[0] in Tg:
              contigs = Tg[line[0]]
              cts= ",".join(str(i) for i in contigs)
              taxa= Tx[line[0]]
              mini = min(It[line[0]])
              maxi = max(It[line[0]])
#             print(It[line[0]]) 
              mm = [mini, maxi]
              str_mm = "-".join( str(i) for i in  mm)
              data = [line[0], line[3], line[4], line[5], line[6], str_mm , taxa, cts ]
              str_out="\t".join(str(i) for i in data)
              str_out+="\n"
              line_out.write(str_out)
         else:
              print("Error not contigs for line[0]")


       
