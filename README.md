# rna_seq

## Overview:
RNA-Seq pipelines that uses [HISAT2](http://daehwankimlab.github.io/hisat2/) and [Kallisto](https://pachterlab.github.io/kallisto/) for alignment (pseudo-alignment and abundance calculations in case of the latter).

## Download Data:
The dataset downloaded was [GSE120534](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA493296). Once the accession list is downloaded, download the SRA toolkit (https://github.com/ncbi/sra-tools) and run the following:
```bash
prefetch --option-file ../SRR_Acc_List_1.txt
```
Once the files are download, retrieve the FASTQ files the following way:
```bash
for file in *; do fastq-dump --split-files "$file"; done
```

## Requirements:
1. [Kallisto](https://pachterlab.github.io/kallisto). Alternatively, `conda install kallisto`
2. [HISAT2](http://daehwankimlab.github.io/hisat2/).
3. [Transcriptome index](https://github.com/pachterlab/kallisto-transcriptome-indices/releases) for *Homo sapiens* if Kallisto is used. Alternatively, index file can also be built using `kallisto index`.
4. [Genome Index](http://daehwankimlab.github.io/hisat2/download/#index) for *Homo sapiens* if HISAT2  is used. Alternatively, index file can also be built using `hisat2-build`.

## Arguments:
`-a | --aligner-to-use`: Specify `1` if you want to use Kallisto or `2` if you want to use HISAT2. DEFAULT: Kallisto<br/>
`-i | --input-files-directory`: Enter the path of the directory containing FASTQ files.<br/>
`-r | --reference-index`: Enter the path of the reference index. 1. If Kallisto is selected, then enter the path of indexed reference transcriptome (Ends with .idx); 2. If HISAT2 is selected, then enter the path of indexed reference genome with the prefix<br/>
`-o | --output-directory`: Enter name of directory which will contain output of respective aligner.

## Script execution:

### 1. Kallisto:
- Run `./aligner_wrapper.py -a 1 -i <FASTQ files directory> -r <reference index directory> -o <aligner output directory>` to use Kallisto for pseudo-alginment and generate the abundance files.<br/>
Alternatively, the following bash script can be run in the directory where the input files are present:
```bash
for f in `ls *.fastq | sed 's/_[12].fastq//g' | sort -u`
do
kallisto quant -i ../reference/homo_sapiens/transcriptome.idx -o kallisto_output/${f} ${f}_1.fastq ${f}_2.fastq
done
```

### 2. HISAT2:
- Run `./aligner_wrapper.py -a 2 -i <FASTQ files directory> -r <reference index directory>genome -o <aligner output directory>` to use HISAT2 for alignment.<br/>
Alternatively, the following bash script can be run in the directory where the input files are present:
```bash
for f in `ls *.fastq | sed 's/_[12].fastq//g' | sort -u`
do
hisat2 -x ../reference/grch38/genome -1 ${f}_1.fastq -2 ${f}_2.fastq -S ${f}.sam
done
```
