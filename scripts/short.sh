#!/bin/bash

input=$1
output=$2

awk -F "\t" '{for (i=1; i<=NF-2; i++) {printf ("%s", $i); if (i<NF-2) printf "\t"} print ""}' $input > $output

