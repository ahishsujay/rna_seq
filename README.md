# rna_seq

## Overview:
RNA-Seq pipelines that uses [HISAT2](http://daehwankimlab.github.io/hisat2/) for alignment.

## Steps:
### 1. HISAT2:
- Run `hisat2_wrapper.py` to use HISAT2 for alignment.<br/>
Alternatively, the following bash script can be run in the directory where the files are present:
```bash
for f in `ls *.fastq | sed 's/_[12].fastq//g' | sort -u`
do
hisat2 -x ../reference/grch38/genome -1 ${f}_1.fastq -2 ${f}_2.fastq -S ${f}.bam
done
```
