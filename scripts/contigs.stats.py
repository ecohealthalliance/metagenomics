import sys
from Bio import SeqIO
import numpy as np

## Goal 

## Determine the distribution of the length of the sequences and some basic statistics of this distribution

input_file=sys.argv[1]
output_file=sys.argv[2]
sample=sys.argv[3]

## arguments
## file: fasta file contigs supported by reads
## sample: sample id

## output
## basic stats from the length distribution of the sequences

fasta_seqs = SeqIO.parse(open(input_file),'fasta')
lengths = [ len(seq.seq) for seq in fasta_seqs ]
array = np.array(lengths)
qts = np.quantile(array, (0, 0.25, 0.5, 0.75, 1))
mesures = [ 'Min', 'l_quartile', 'Median', '3_quartile', 'Max']
delimiter = '\n'
num_seqs = array.shape[0]


with open(output_file, 'w') as line:
  line.write('{:s}: {:s}\n'.format("Sample", sample))  
  line.write('{:s}: {:d}\n'.format("Number_of_contigs", num_seqs))
  line.write(delimiter.join(["%s: %.2f" % (v,u) for v,u in zip(mesures, qts)]))
