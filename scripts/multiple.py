import sys
import numpy as np

## Goals 

## Identify a best hit for each of the contigs according to the frequency of the target sequence in the alignments. Among all the best hits of a contig, the one that was most observed as the best alignment for the other contigs is chosen.

input_file=sys.argv[1]
blast_file=sys.argv[2]
out_best=sys.argv[3]
out_multiple_best=sys.argv[4]
out_bed=sys.argv[5]

## arguments
## input_file: file output of best.multiple.pl, fields: id_contig id(s)_best_hit(s)
## blast_file: blastn output file  
## out_best: only the best hit for each sequence
## out_multiple_best: multiple best hits for each contig
## out_bed:  bed file of best hit

## Determines the number of times a target is the best hit of a query in the dataset
## arguments: file handled

def num_qrs_by_target (input_file): 
    with open(input_file, 'r') as handled:
        targets = {}
        for line in handled:
            line=line.strip()
            (ct, tags) = line.split('\t')
            list_tgs = tags.split(',')
            for tg in list_tgs:
                if not tg in targets:
                    targets[tg] = {}
                    targets[tg][ct] = "" 
                else:
                    targets[tg][ct] = "" 
    
        targ_ct = { pr:  len(queries) for pr, queries in targets.items()}
        return targ_ct


## Identify in a list the target that is the most observed in the data
## arguments: 
## list : list best hit(s), dictionary: target => number of queries 

def freq_best (list, dictionary):
    
    sub_dic = { i:dictionary[i]  for i in list if i in dictionary }
    sorted_sub = sorted(sub_dic, key=sub_dic.get, reverse=True)  
    return sorted_sub   


## alignment information for each of the query-target pairs  

blast = {}
bed =  {}


with open(blast_file, 'r') as f:
    for line in f:
        list = line.split('\t')        
        query = list[0]
        target = list[1]
        if not query in blast:
            blast[query] = {}
            blast[query][target] = line
        else:
            blast[query][target] = line
        if list[8] >=  list[9]:
            tmp = list[8]
            list[8] = list[9]
            list[9] = tmp
  ## bed features localization based in 0
        pi = int(list[8])
        pf = int(list[9]) 
        pi -= 1
        pf -= 1        
        if pi > pf:
            temporal = pf
            pf = pi
            pi = temporal 
        string = "\t".join([ list[1],  str(pi),  str(pf),   list[0] ])
        string+="\n"
        if not query in bed:
             bed[query] = {}
             bed[query][target] = string
        else:
             bed[query][target] = string
             
## dictionary where the key is a target id and the value is the number of times that it is the best alignment for a query
    
dictic=num_qrs_by_target(input_file)


## Identify the best alignment for each of the contigs
## two output files, one with the best hit for each contig and one with all the best hits when this is the case

hd =  [ 'qseqid', 'sseqid' , 'pident', 'length', 'mismatch' , 'gapopen' ,  'qstart' ,  'qend'  ,  'sstart'  , 'send' , 'evalue' ,  'bitscore' ,  'qcovs' ,  'slen' ,  'sgi' ,  'sscinames',   'stitle' ,  'staxid' ]

header = "\t".join(hd)+"\n"


with open(input_file, 'r') as f,  open(out_best, 'w') as wr_best, open(out_multiple_best, 'w') as wr_mul, open(out_bed, 'w') as wr_bed:
    wr_best.write(header)
    wr_mul.write(header)
    for line in f:      
        line=line.strip()
        (ct, tags) = line.split('\t')    
        list_tgs = tags.split(',')
        list_tgs_ord= freq_best(list_tgs, dictic)
        best = list_tgs_ord[0]
        ln_best = blast[ct][best]
        wr_best.write(ln_best)
        ln_bed = bed[ct][best]
        wr_bed.write(ln_bed)
        for tg  in  list_tgs:
            ln_mul = blast[ct][tg]
            wr_mul.write(ln_mul)
