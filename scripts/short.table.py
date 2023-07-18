import sys
import argparse

## Goal 

# Generate the report of the taxonomic classification of the contigs

parser = argparse.ArgumentParser()
parser.add_argument('--input_file', type=str, required=True)
parser.add_argument('--output_file', type=str, required=True)
parser.add_argument('--output_accesions', type=str)
args = parser.parse_args()


   
## arguments
## --input_file: output of adding the taxonomy of the target sequence to the blast results
## --output_file: report with the most important information of the blast alignment and the taxonomic classification for each of the contigs
## --output_accesions: file with the Identifiers of the target sequences, used to retrieve the gff3 files  

## output
## Report with the most relevant information of taxonomic classification and blast aligments 
## output_accesions: id target sequences

if args.output_accesions:
   acc = {}


with open(args.input_file, 'r') as lines, open(args.output_file, 'w') as  out_line:

   for line in lines:
       line=line.strip()
       line=line.split('\t')
       if args.output_accesions:
          if not  line[1] in acc and line[1] !=  "sseqid":
             acc[line[1]] = ""

             

       
# If missing values assign NA
        
       new = [  i  if not i == "" else "NA" for i in line ]


# Choose the most important data
          

       impr = [ new[0], new[3],  new[1],  new[13], new[2], new[12], new[6], new[7], new[8], new[9], new[10], new[11] ,  new[14], new[15], new[16], new[17], new[18], new[19], new[20], new[21], new[22], new[23] ,  new[24]  ]
       str_out="\t".join(str(i) for i in impr)
       str_out += "\n"
       out_line.write(str_out)


if args.output_accesions:
   with open(args.output_accesions, 'w') as line:
      claves =  acc.keys()
      string =  "\n".join(str(i) for i in claves)
      string += "\n"
      line.write(string)
