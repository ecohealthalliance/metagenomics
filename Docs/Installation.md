### Operating system
<p align="justify">
The pipeline was developed on Ubuntu 22.04.2. The software was developed with the following system requirements:
 </p>

- Python 3.7.12
- conda 23.1.0
- git 2.34.1
- perl v5.26.2

### Pipeline software installation

<p align="justify">
1. Clone the GitHub pipeline repository.
</p>

```
git clone https://github.com/ecohealthalliance/metagenomics.git

```
<p align="justify">
Change directory
</p>

```
cd metagenomics

```
<p align="justify">
At this point, you could activate the scripts execution permissions with the following command:
</p>

```
chmod +x scripts/*

```
 
2. Installation of software requires conda. Installation of conda is described [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). Once installed conda, run the following commands: 

- Activate conda base environment

```
conda activate base

```
- Install mamba and then create the iSNVs environment by installing snakemake at the same time. 

```
conda install -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake python=3.7

```

- Activate snakemake environment
```
conda activate snakemake

```

- Install software

```
conda install --yes --file packages.txt

```
