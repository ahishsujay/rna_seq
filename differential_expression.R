#!/usr/bin/env Rscript

#List of packages required:
packages = c("sleuth", "biomaRt", "tximport", "tximportData", "DESeq2", "argparse")

#If package is not present, then install:
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(new.packages)

#Loading argparse to take in arguments:
suppressPackageStartupMessages(library("argparse"))

#Create parser object
parser <- ArgumentParser()
parser$add_argument("-de", "--differential-expression-analysis", type = "character", default = 'DESeq', help = "Enter Sleuth or DESeq to perform DE analysis. [default %(default)s")
parser$add_argument("-o", "--kallisto-output-directory", help = "Enter Kallisto's output directory name.")
parser$add_argument("-m", "--meta-data", help = "Enter file name of meta data.")
parser$add_argument("-p", "--padj-cutoff", type="double", default=0.05, help = "Enter cutoff value for padj. [default %(default)s]")
args <- parser$parse_args()

if( tolower(args$differential_expression_analysis) == "deseq" ) {

  print("Using DESeq...")

  #Loading packages:
  suppressPackageStartupMessages(library("tximport"))
  suppressPackageStartupMessages(library("tximportData"))
  suppressPackageStartupMessages(library("DESeq2"))

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
  #sampleTable <- data.frame(condition = factor(c(rep("control",5), rep("knockout", 5))))
  data_info <- read.table(file.path("data_info.txt"), header=TRUE)
  sampleTable <- data.frame(condition = factor(dplyr::pull(data_info, sample)))
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

} else if (tolower(args$differential_expression_analysis) == "sleuth" ) {

  print("Using Sleuth...")

  #Loading packages:
  suppressPackageStartupMessages(library("sleuth"))
  suppressPackageStartupMessages(library("biomaRt"))

  #Function to get ENSEMBL geneID and geneName for the transcripts:
  transcript2gene <- function(){
    mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
    t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
    return(t2g)
  }
  t2g <- transcript2gene()

  sample_id <- dir(file.path(args$kallisto_output_directory))
  kal_dirs <- file.path(args$kallisto_output_directory, sample_id, "abundance.h5")
  data_info <- read.table(file.path(args$meta_data), header=TRUE)
  s2c <- dplyr::mutate(data_info, path = kal_dirs)

  #Initializing sleuth object:
  print("Initializing sleuth object...")
  so <- sleuth_prep(s2c, target_mapping=t2g)

  #Fitting the reduced model:
  print("Fitting the reduced model...")
  so <- sleuth_fit(so, ~1, 'reduced')

  #Fitting the full model:
  print("Fitting the full model...")
  so <- sleuth_fit(so, ~condition, 'full')

  #Performing test:
  print("Performing test...")
  so <- sleuth_lrt(so, 'reduced', 'full')

  #Examining models that have been fit:
  models(so)

  #Test results:
  sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt')
  sleuth_significant <- dplyr::filter(sleuth_table, qval <= args$padj_cutoff)

  #Writing to csv:
  print("Writing CSVs...")
  write.csv(x=sleuth_table, file="kallisto_sleuth_tableResults.csv")
  write.csv(x=sleuth_significant, file="kallisto_sleuth_significantResults.csv")

} else {
    print("Error.")
}
