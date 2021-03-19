#!/usr/bin/env Rscript

#List of packages required:
packages = c("tximport", "tximportData", "DESeq2", "argparse")

#If package is not present, then install:
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(new.packages)

#Loading packages:
suppressPackageStartupMessages(library("tximport"))
suppressPackageStartupMessages(library("tximportData"))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()
parser$add_argument("-o", "--kallisto-output-directory", help = "Enter Kallisto's output directory name.")
parser$add_argument("-p", "--padj-cutoff", type="double", default=0.05, help = "Enter cutoff value for padj. [default %(default)s]")
args <- parser$parse_args()

#tximport:
dir <- system.file("extdata", package = "tximportData")
output_dir = args$kallisto_output_directory
samples <- dir(file.path(output_dir))
files <- file.path(output_dir, samples, "abundance.h5")
names(files) <- paste0("sample", 1:length(samples))

#Checking if the files exist, and if they do then use tximport:
if (isTRUE(all(file.exists(files)))) {
  print("Files exist, running tximport...")
  txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
} else {
  print("Error. Check the path of the input files.")
}

#DESeq2:
print("Running DESeq2...")
sampleTable <- data.frame(condition = factor(c(rep("control",5), rep("knockout", 5)))) #Change this to reflect for your experiment
rownames(sampleTable) <- colnames(txi.kallisto$counts)

#Construct DESEQDataSet Object:
dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)

#Run the DESeq function:
print("Running DESeq function...")
dds <- DESeq(dds)

#Saving results of DESeq:
res <- results(dds)

#Order by the adjusted pvalue:
print("Ordering by adjusted pvalue...")
table(res$padj<args$padj_cutoff)
res <- res[order(res$padj), ]

#Merging with normalized counts data:
print("Merging with normalized counts data...")
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "transcriptID"

#Writing to csv:
print("Writing results.csv...")
write.csv(x=resdata, file="kallisto_tximport_deseq2_results.csv")