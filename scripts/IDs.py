import sys
from Bio import SeqIO

## Goal 
## new id: add Id sample to Id sequence

input_file=sys.argv[1]
output_file=sys.argv[2]
sample=sys.argv[3]

## arguments
## input_file: fasta file contigs supported by reads
## output_file: output fasta file with new sequences Ids
## sample: sample id


with open(output_file, 'w') as out_file, open(input_file, "r") as input:
    fasta_sequences = SeqIO.parse(input,'fasta')
    for fasta in fasta_sequences:
        fasta.id +="_"+sample
        SeqIO.write(fasta, out_file, "fasta")
