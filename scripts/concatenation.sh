#!/bin/bash

sample=$1

tmpe="diamond/concatenation."$sample
touch $tmpe
out=diamond/$sample.dmnd.tsv

for file in diamond/$sample.dmnd.tsv_*
do
        tmpfile=$(mktemp diamond/tempotxt.XXXXXXX)
	cat $file $tmpe > $tmpfile
 	mv $tmpfile $tmpe
done

if [ -s $tmpe  ]
then
    mv $tmpe $out
fi
