from os.path import join
import re
import string
import socket
configfile: "./config.yaml"

DIR = "data"
IDS, = glob_wildcards(join(DIR, "{id}_R1.fastq"))


rule all:
     input:
      expand("table_genomes_blastn/{id}.table.genomes.csv", id=IDS)


## Count the number of raw reads 

rule countRaw:
    input:
      expand("data/{id}_R1.fastq", id=IDS)
    output:
      "count/count.raw.txt"
    threads:
      1
    log:"logs/countRaw.log"
    shell:"""

   scripts/count_fastq.sh {input} {output}  2> {log}

    """

# Trimming reads

rule trimming:
    input:
        fw="data/{id}_R1.fastq",
        rv="data/{id}_R2.fastq",
        lk="count/count.raw.txt"
    output:
        tf="trimming/{id}_R1.clean.fastq",
        tr="trimming/{id}_R2.clean.fastq",
        uf="trimming/{id}_R1.unpaired.fq",
        ur="trimming/{id}_R2.unpaired.fq"
    threads:
        2
    params:
        ad=str(config["adapters"]),
    log:"logs/trimming/{id}.log"
    shell:"""

trimmomatic PE -threads {threads} {input.fw} {input.rv} {output.tf} {output.uf} {output.tr} {output.ur} ILLUMINACLIP:{params.ad}:2:30:10  LEADING:6 TRAILING:6  SLIDINGWINDOW:4:15 MINLEN:40 2> {log}


    """

## Quality control with fastqc

rule fastqc:
    input:
        tf="trimming/{id}_R1.clean.fastq",
        tr="trimming/{id}_R2.clean.fastq",
    output:
        qctf="QC/{id}_R1.clean_fastqc.html",
        qctr="QC/{id}_R2.clean_fastqc.html"
    threads:
        1
    params:
        dir="QC"
    log:"logs/qc/{id}.log"
    shell:"""

     fastqc  --threads {threads} {input.tf} --outdir={params.dir}  2> {log} 
     fastqc  --threads {threads} {input.tr} --outdir={params.dir}  2>> {log}

      
    """

## Assembly of contigs with megahit

rule megahit:
     input:
      fr="trimming/{id}_R1.clean.fastq",
      rv="trimming/{id}_R2.clean.fastq",
      uf="trimming/{id}_R1.unpaired.fq",
      ur="trimming/{id}_R2.unpaired.fq",
      lk="QC/{id}_R1.clean_fastqc.html"
     output:
      "megahit/{id}.megahit.contigs.fasta"
     threads:
      5
     params:
      dr="megahit_{id}",
      lg=int(config["length_megahit"]) 
     log:"logs/megahit/{id}.log"
     shell:"""

   megahit  --min-contig-len {params.lg}   -t {threads} -1  {input.fr}  -2 {input.rv}  -r {input.uf},{input.ur} -o {params.dr} 2> {log}
   mv {params.dr}/final.contigs.fa  {output} 2>> {log}
   rm -r {params.dr} 2>> {log}


"""

## Reduce redundancy of sequences with cd-hit-est

rule   Clustering:
       input:
        "megahit/{id}.megahit.contigs.fasta"
       output:
        cl="NoRedundancy/{id}.cluster_scaffolds.clstr",
        sq="NoRedundancy/{id}.cluster_scaffolds.fasta"
       params:
        ref="NoRedundancy/{id}.cluster_scaffolds"
       threads:
        3
       log:"logs/Clustering/{id}.log"
       shell:"""


  cd-hit-est -M 20000 -T {threads} -i {input} -o {params.ref}  -c 0.99  -G 0  -aS 1 2> {log}
  mv {params.ref} {output.sq} 2>> {log}

"""

## Alignment reads to contigs

rule reads_2_contigs:
     input:
       fa="NoRedundancy/{id}.cluster_scaffolds.fasta",
       fr="trimming/{id}_R1.clean.fastq",
       rv="trimming/{id}_R2.clean.fastq",
       uf="trimming/{id}_R1.unpaired.fq",
       ur="trimming/{id}_R2.unpaired.fq"
     output:
       sm="reads_2_contigs/{id}.reads.to.contigs.sam",
       rp="reads_2_contigs/{id}.bowtie.reads2contigs.report.txt"
     threads:
       6
     log:"logs/reads_2_contigs/{id}.log"
     shell:"""


  bowtie2-build  {input.fa} {input.fa} 2> {log}

  bowtie2 -p  {threads} -q -x  {input.fa}  -1 {input.fr} -2 {input.rv}  -U {input.uf},{input.ur}  -S {output.sm} --no-unal  2>> {log}  1>  {output.rp}

 scripts/metrics.bowtie.pl {log} > {output.rp} 

"""

### From sam to bam file

rule getBam_reads_2_contigs:
       input:
         "reads_2_contigs/{id}.reads.to.contigs.sam"
       output:
         "reads_2_contigs/{id}.reads.to.contigs.bam"
       threads:
         4
       log:"logs/getBam_reads_2_contigs/{id}.log"
       shell:"""

        samtools view --threads {threads}  -F 4  -S -b {input} > {output}  2> {log}

   """

# Sort bam file by coordinates

rule sortBowtie_reads_2_contigs:
       input:
         "reads_2_contigs/{id}.reads.to.contigs.bam"
       output:
         "reads_2_contigs/{id}.reads.to.contigs.sort.bam"
       threads:
         4
       log:"logs/sortBowtie_reads_2_contigs/{id}.log"
       shell:"""

           samtools sort --threads {threads} {input} -o {output} 2> {log}

        """

# Get index from the sorted bam file

rule indexBowtie_reads_2_contigs:
       input:
        "reads_2_contigs/{id}.reads.to.contigs.sort.bam"
       output:
        "reads_2_contigs/{id}.reads.to.contigs.sort.bam.bai"
       threads:
         1
       log:"logs/indexBowtie_reads_2_contigs/{id}.log"
       shell:"""

           samtools index {input} 2> {log}

        """


## Estimate the coverage and the RPKM of the alignment of the reads to the contigs

rule coverage:
        input:
          sm="reads_2_contigs/{id}.reads.to.contigs.sam",
          lk="reads_2_contigs/{id}.reads.to.contigs.sort.bam.bai"  
        output:
          rp="reads_2_contigs/{id}.rpkm.csv",
          cv="reads_2_contigs/{id}.coverage.csv"

        threads:
          2
        log:"logs/coverage/{id}.log"
        shell:"""


  pileup.sh  in={input.sm} out={output.cv}   rpkm={output.rp} 2> {log}


"""

## Identification of contigs supported by reads
## Two parameters are used to filter the sequences according to the read alignment to the contigs.
## The first is the length of the contig covered by the reads coverag ; while the second is the number of aligned reads.


rule contigs_supported:
     input:
      "reads_2_contigs/{id}.coverage.csv"
     output:
      "contigs_support/{id}.contigs.supported.csv"
     threads:
      1
     params:
      cv=str(config["cov_contigs"]),
      mn=str(config["min_reads_cont"])
     log:"logs/contigs_supported/{id}.log"
     shell:"""

     python3  scripts/cov.contigs.py {input} {params.cv} {params.mn}  > {output} 2> {log}

"""

## Sequences supported by reads and add sample id to the sequence id

rule contigs_supported_fasta:
     input:
      fl="contigs_support/{id}.contigs.supported.csv",
      sq="NoRedundancy/{id}.cluster_scaffolds.fasta"
     output:
      fa="contigs_support/{id}.contigs.supp.id.fasta",
      st="contigs_support/{id}.contigs.stats.length.csv",
      ft=temp("contigs_support/{id}.contigs.supported.fasta"),
      lt=temp("contigs_support/{id}.support.list.csv")
     threads:
      4
     log:"logs/contigs_supported_fasta/{id}.log"
     shell:"""

    cut -f1 {input.fl}  > {output.lt}  2> {log}
    seqkit grep --pattern-file {output.lt} {input.sq} > {output.ft} 2>> {log}
    python3 scripts/IDs.py  {output.ft} {output.fa} {wildcards.id} 2>> {log}     
    python3 scripts/contigs.stats.py {output.fa} {output.st} {wildcards.id} 2>> {log}

"""

## Diamond alignment against nr NCBI database    
## parameters: evalue blast threshold and the referece sequence database 


rule diamond_initialization:
	input:
	  "contigs_support/{id}.contigs.supp.id.fasta"
	output:
	  "diamond/{id}.dmnd.tsv"
	params:
	  db=str(config["diamond_ref"]),
	  id=int(config["diamond_id"]),
	  ev=float(config["diamond_evalue"]),
	  mm=int(config["diamond_block"]),
	  ci=int(config["diamond_chunck"]),
	  tm="tmp_parallel/temporary_{id}",
	  pt="tmp_parallel/parallel_temporary_{id}"
	threads:
	  45
	log:"logs/diamond_initialization/{id}.log"
	shell:"""

    scripts/make.directories.sh {wildcards.id}

   diamond blastx    --threads  {threads}  --min-orf 1 --fast  --id  {params.id}  --evalue {params.ev}  -b   {params.mm} -c {params.ci}  -d  {params.db} -q  {input}  -o  {output} --outfmt 6  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore  qcovhsp slen qtitle sscinames stitle  staxids --multiprocessing --mp-init --tmpdir {params.tm}  --parallel-tmpdir {params.pt} 2>>  {log}  


   diamond blastx   --threads  {threads}  --min-orf 1 --fast  --id  {params.id}  --evalue {params.ev}  -b   {params.mm} -c {params.ci}  -d  {params.db} -q  {input}  -o  {output} --outfmt 6  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore  qcovhsp slen qtitle sscinames stitle  staxids --multiprocessing --tmpdir {params.tm}  --parallel-tmpdir {params.pt} 2>>  {log}

 scripts/concatenation.sh {wildcards.id}
  

"""

## Processing of diamond alignment results

rule Multiple:
        input:
          "diamond/{id}.dmnd.tsv"
        output:
          "diamond/{id}.multiple.csv"
        threads:
           1
        log:"logs/Multiple/{id}.log"
        shell:"""

  scripts/best.multiple.pl {input}  > {output} 2> {log}


"""

## Identify the best alignment for each of the contigs    

rule Best:
        input:
          mt="diamond/{id}.multiple.csv",
          bs="diamond/{id}.dmnd.tsv"
        output:
          sg="best/{id}.best.csv",
	  bd="table_genomes_diamond/{id}.best.bed",
          mt="best/{id}.best.multiples.csv"
        threads:
          1
        log:"logs/best/{id}.log"	  
        shell:"""

   python3 scripts/multiple.py {input.mt} {input.bs} {output.sg} {output.mt} {output.bd}  2>  {log}


"""

## add taxonomic information 
## parameters: path to nodes.dmp and names.dmp files

rule link_taxa:
        input:
          sg="best/{id}.best.csv",
          mt="best/{id}.best.multiples.csv"
        output:
          sg="taxa/{id}.best.taxa.csv",
          mt="taxa/{id}.best.multiples.taxa.csv"
        threads:
          1
        params:
          no=str(config["nodes"]),
          na=str(config["names"])
	log:"logs/link_taxa/{id}.log"  
        shell:"""


    scripts/link_taxinfo.pl -h -v  {input.sg} {params.no} {params.na}  species,genus,family,order,class,phylum,superkingdom 1>   {output.sg}  2>  {log}

    scripts/link_taxinfo.pl -h -v  {input.mt} {params.no} {params.na}  species,genus,family,order,class,phylum,superkingdom 1>   {output.mt}  2>>  {log}



"""


## Some of the sequences may be non-biological sequences (Exp: viral vectors). 
## These are characterized by the fact that the taxonomic annotation is shorter. 
## The script searches for lines that do not have complete taxonomic annotations; 
## these sequences are sent to the construct  directory.


rule clean_taxa:
        input:
          sg="taxa/{id}.best.taxa.csv",
          mt="taxa/{id}.best.multiples.taxa.csv"
        output:
          cs="clean/{id}.best.clean.csv",
	  cm="clean/{id}.best.multiples.clean.csv",
	  os="constructs/{id}.csv",
          om="constructs/{id}.count.csv"
        threads:
           1
        log:"logs/clean_taxa/{id}.log"  
        shell:"""


       python3 scripts/constructs.py {input.sg} {input.mt} {output.cs} {output.cm} {output.os} {output.om} 


"""
## Table with the taxonomic classification of the sequences according to the results of diamond.

rule contigs_report:
        input:
          sg="clean/{id}.best.clean.csv",
          mt="clean/{id}.best.multiples.clean.csv"
        output:
          sg="table_contigs_diamond/{id}.best.report.csv",
          mt="table_contigs_diamond/{id}.best.multiples.report.csv",
          am="table_contigs_diamond/{id}.best.acc.csv"
        threads:
          1
        log:"logs/contigs_report/{id}.log"  
        shell:"""

        python3 scripts/short.table.py --input_file {input.sg}   --output_file {output.sg}  --output_accesions  {output.am} 2>  {log}  
        python3 scripts/short.table.py --input_file {input.mt}   --output_file {output.mt}     2>>  {log}
 
    

"""


### Build bed files for the target sequence


rule bed_targets:
	input:
	  "table_contigs_diamond/{id}.best.report.csv"
	output:
	  "bed/{id}.bed.target.csv"
	threads:
	  1
	log:"logs/bed_targets/{id}.log"
	shell:"""

   python3 scripts/make.bed.py  {input} {output}


"""


## Table of target sequences. 
## If there are several contigs aligning to the same target sequence,
## this file describes the number of contigs, the genome coverage, the identity range among others variables.

rule table_genomes:
	input:
          bd="table_genomes_diamond/{id}.best.bed",
	  sg="table_contigs_diamond/{id}.best.report.csv",
	  bt="bed/{id}.bed.target.csv"
	output:
          cv="table_genomes_diamond/{id}.coverage.targets.csv",
	  tb="table_genomes_diamond/{id}.table.genomes.csv"
	threads:
	  1
	log:"logs/table_genomes_diamond/{id}.log"
	shell:"""


   bedtools coverage -a   {input.bt}  -b {input.bd} >  {output.cv}
   
   python3 scripts/tables.genomes.diamond.py  --coverage_file  {output.cv} --best_file  {input.sg} --output_file  {output.tb}



"""

## Fasta file of the sequeces associated to viruses by diamond    
## In the second part of the pipeline are sequences identified as viral, they will be aligned with BLASTn to obtain the final classification.


rule viruses_fasta:
	input:
          sg="clean/{id}.best.clean.csv",
	  fa="contigs_support/{id}.contigs.supp.id.fasta",
	  lk="table_genomes_diamond/{id}.table.genomes.csv"
	output:
          ls=temporary("viruses_fasta/{id}.viruses.csv"),
	  fa="viruses_fasta/{id}.viruses.fasta"
	threads:
	  1
	log:"logs/viruses_fasta/{id}.log"
	shell:"""

       python3 scripts/recover.viruses.py {input.sg} {output.ls} 2> {log} 
       seqkit grep --pattern-file  {output.ls} {input.fa} > {output.fa} 2>> {log}    


"""

## split the fasta file into 15 different files and create a list
## with the identifiers of these files to be executed in the following rule.

checkpoint pre_blastn:
        input:
          "viruses_fasta/{id}.viruses.fasta"
        output:
          blast_dir=directory("pre_blastn/{id}")
        threads:
          1
        log:"logs/pre_blastn/{id}.log"  
        shell:"""                                                                                 

      mkdir pre_blastn/{wildcards.id}  2>  {log}
      scripts/blastn.split.sh {input} {wildcards.id} 15   2>>  {log}

"""



rule blastn_alig:
        input:
          "pre_blastn/{id}/split.{ext}.fasta"
        output:
          "blastn_alignment/{id}/{ext}.out.csv"
        params:
          ref=str(config["ncbi_nt"]),
          ev=float(config["evalue_blast"]) 
        threads:
          5
	log:"logs/blastn_alig/{id}.{ext}.log"
        shell:"""   


 blastn -db  {params.ref}  -perc_identity 40   -evalue {params.ev} -word_size 11 -num_threads {threads}   -outfmt '6 std qcovs slen sgi sscinames stitle  staxids' -query {input} >  {output}



"""

def aggre_blastn(wildcards):
     checkpoint_output = checkpoints.pre_blastn.get(**wildcards).output[0]
     return expand("blastn_alignment/{id}/{ext}.out.csv",
           id=wildcards.id,
           ext=glob_wildcards(os.path.join(checkpoint_output, "split.{ext}.fasta")).ext)

## Concatenate the BLASTn results  

rule  aggregate:
        input:
          aggre_blastn
        output:
          "blastn/{id}.csv"
        threads:
          1  
        log:"logs/aggregate/{id}.log"
        shell:"""                                                                                 

   
    cat {input} >  {output}  2>  {log}
     

"""

## Processing of the BLASTn results to obtain the best hit of each of the contigs.

rule multiple_blastn:
        input:
          "blastn/{id}.csv"
        output:
          mp="blastn/{id}.multiple.csv",
          pb="blastn/{id}.pre.best.csv"
        threads:
          1  
        log:"logs/multiple_blastn/{id}.log"
        shell:"""


   scripts/best.multiple.dev.blastp.II.pl  {input}  {output.mp} {output.pb} 2> {log}


"""

## Identify the best alignment for each of the contigs    

rule best_blastn:
        input:
          mt="blastn/{id}.multiple.csv",
          bs="blastn/{id}.pre.best.csv"
        output:
          sg="best_blastn/{id}.best.csv",
          bd="table_genomes_blastn/{id}.best.bed",
          mt="best_blastn/{id}.best.multiples.csv"
        threads:
          1  
        log:"logs/best/{id}.log"	  
        shell:"""

   python3 scripts/multiple.py {input.mt} {input.bs} {output.sg} {output.mt} {output.bd}  2>  {log}


"""


## add taxonomic information 
## parameters: path to nodes.dmp and names.dmp files

rule link_taxa_blastn:
        input:
          sg="best_blastn/{id}.best.csv",
          mt="best_blastn/{id}.best.multiples.csv"
        output:
          sg="taxa_blastn/{id}.best.taxa.csv",
          mt="taxa_blastn/{id}.best.multiples.taxa.csv"
        params:
          no=str(config["nodes"]),
          na=str(config["names"])
	threads:
          1    
	log:"logs/link_taxa_blastn/{id}.log"  
        shell:"""


    scripts/link_taxinfo.pl -h -v  {input.sg} {params.no} {params.na}  species,genus,family,order,class,phylum,superkingdom 1>   {output.sg}  2>  {log}

    scripts/link_taxinfo.pl -h -v  {input.mt} {params.no} {params.na}  species,genus,family,order,class,phylum,superkingdom 1>   {output.mt}  2>>  {log}



"""

##  Identification and removal of contaminating sequences

rule clean_taxa_blastn:
        input:
          sg="taxa_blastn/{id}.best.taxa.csv",
          mt="taxa_blastn/{id}.best.multiples.taxa.csv"
        output:
          cs="clean_blastn/{id}.best.clean.csv",
	  cm="clean_blastn/{id}.best.multiples.clean.csv",
	  os="constructs_blastn/{id}.csv",
          om="constructs_blastn/{id}.count.csv"
        threads:
          1  
        log:"logs/clean_taxa_blastn/{id}.log"  
        shell:"""


       python3 scripts/constructs.py {input.sg} {input.mt} {output.cs} {output.cm} {output.os} {output.om} 


"""

## The "{id}.best.report.csv" tables describe the taxonomy of each of the
## contigs according to the best hit  identified in the NCBI database.

rule contigs_report_blastn:
        input:
          sg="clean_blastn/{id}.best.clean.csv",
          mt="clean_blastn/{id}.best.multiples.clean.csv"
        output:
          sg="table_contigs_blastn/{id}.best.report.csv",
          mt="table_contigs_blastn/{id}.best.multiples.report.csv",
          am="table_contigs_blastn/{id}.best.acc.csv"
        threads:
          1  
        log:"logs/contigs_report_blastn/{id}.log"  
        shell:"""

        python3 scripts/short.table.py --input_file {input.sg}   --output_file {output.sg}  --output_accesions  {output.am} 2>  {log}  
        python3 scripts/short.table.py --input_file {input.mt}   --output_file {output.mt}     2>>  {log}
 
    

"""

rule bed_targets_blastn:
	input:
	  "table_contigs_blastn/{id}.best.report.csv"
	output:
	 "bed_blastn/{id}.bed.target.csv"
	threads:
	  1  
	log:"logs/bed_targets_blastn/{id}.log"
	shell:"""

   python3 scripts/make.bed.py  {input} {output}


"""

## Tables with the summarized information for the target sequences; 
## this includes the number of sequences aligned, The proportion of the length of the reference sequence covered by reads
## and the taxonomic classification. 

rule table_genomes_blastn:
	input:
          bd="table_genomes_blastn/{id}.best.bed",
	  sg="table_contigs_blastn/{id}.best.report.csv",
	  bt="bed_blastn/{id}.bed.target.csv"
	output:
          cv="table_genomes_blastn/{id}.coverage.targets.csv",
	  tb="table_genomes_blastn/{id}.table.genomes.csv"
	threads:
	  1
	log:"logs/table_genomes_blastn/{id}.log"
	shell:"""


   bedtools coverage -a  {input.bt}  -b {input.bd} >  {output.cv}
   python3 scripts/tables.genomes.diamond.py  --coverage_file  {output.cv} --best_file  {input.sg} --output_file  {output.tb}



"""
