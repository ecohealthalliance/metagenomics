#!/bin/bash

input=$1
sample=$2
nm_cores=$3

#The goal is to split a large fasta file into files containing 10 sequences.


## get the fasta file name, remove directory
file_input="${input/contigs_support\//}"


if [ ! -d "pre_blastn/$sample" ]
then
    mkdir pre_blastn/$sample
fi


cp $input  pre_blastn/$sample/split.fasta
cd pre_blastn/$sample

#### split fasta file in files with 10 sequences
#### In this version the fasta file is divided between files that have the number of sequences (10)

# awk -v size=10 -v pre=pre.blastn -v pad=12 '/^>/{n++;if(n%size==1){close(f);f=sprintf("%s.%0"pad"d",pre,n)}}{print>>f}' $file_input

####
### According to (https://peerj.com/articles/3486/) it is better to split the fasta file between files of the same size and corresponding to the number of available cores.
## pyfasta can fail when the file is divided between a large number of sub-files, in that case send errorf

pyfasta split -n $nm_cores  split.fasta

rm  split.fasta


