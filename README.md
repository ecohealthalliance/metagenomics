#                                                     Metagenomics pipeline

###### This pipeline performs a taxonomic classification of metagenomic sequences resulting from second-generation technologies.

##  Introduction 
###### Metagenomics is a multi-application tool in virology. This technology has allowed the community to identify new viral sequences and expand our knowledge of the genomic structure, diversity, and ecology of currently known viral families. In this exploratory task, metagenomics has succeeded in identifying pathogenic viruses with zoonotic potential. However, it is possible that in the coming years, the impact of metagenomics will be even greater. New studies are using metagenomics as a tool to test or evaluate ecological and/or evolutionary hypotheses. Several of these projects also integrate metagenomics with other "omic" sciences, including genomics and transcriptomics, among others. It is possible that the empirical data generated from these analyzes could in the future improve and recalibrate epidemiological models that try to predict and anticipate new viral outbreaks. We developed a pipeline that allows processing and obtaining the classification of the sequences obtained in these metagenomic studies. This pipeline has been tested on a wide variety of samples, including environmental samples, feces, and swabs of different tissues from wild and domestic animals. 

## Pipeline

###### The pipeline was designed to process and taxonomically classify second-generation metagenomic sequences.  We briefly describe the main moments of the pipeline (figure below). The raw reads are trimmed and then assembled into _de novo_ contigs. The clean reads are aligned against contigs; The contigs are filtered according to coverage and the number of aligned reads. The contigs supported by the reads are aligned with the diamond software against the non-redundant (NR) protein NCBI database. For each of the contigs, the best hit is chosen, and only the contigs aligning to viruses are retained. These sequences are aligned with the BLASTn software against the non-redundant (NT) nucleotide NCBI database. For each of the sequences that obtains significant alignments, the best hit is identified. The pipeline summarizes the main information of the processing and classification in tables. In the case of the contigs, these tables contain information on the taxonomic classification, the target sequence, the identity and coverage of the alignment, the number of reads aligned to the contigs, the RPM, and RPKM, among others. Tables are also produced that summarize the information from the perspective of the target sequences. This information includes the number and identity of the contigs aligning to a specific target, the coverage of the target sequence by the contigs, and the average identity among others.


![Pipeline](https://github.com/ecohealthalliance/metagenomics/assets/72785049/900cb646-5cfb-4d31-8413-9512ef77afdc)

## Documentation

- Installation
