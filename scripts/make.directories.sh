#!/bin/bash

sample=$1

tmp_par=tmp_parallel/temporary_$sample
tmp_pt=tmp_parallel/parallel_temporary_$sample

if [ ! -d tmp_parallel ]
then
    mkdir tmp_parallel
    mkdir $tmp_par
    mkdir $tmp_pt
else
    if  [ ! -d $tmp_par ]
	then
        mkdir $tmp_par
    fi
    if  [ ! -d $tmp_pt ]
	then
        mkdir $tmp_pt
    fi 
fi
