# -----------------------------------------------------------------------------#
# Mediterranean sponge fungi - ocean acidification project 
# Processing cutadapt-trimmed Reads
# Author: Geoffrey Zahn
# Requirements: ITSx v 1.1b1
# -----------------------------------------------------------------------------#

# Prepare demultiplexed reads in Bash... ####
# Run ITSxpress on all fasta files 


for i in ./seqs/fastqs/*.fastq.gz; do itsxpress --fastq $i  --outfile $i.FungalITS1.fastq.gz --region ITS1 --taxa Fungi -s --threads 4 --log $i.ITSx.log; done
# rename the files
# rename -v -e 's/^(.{7}).*/$1.fastq.gz/' *.FungalITS1.fastq.gz
# No ITS1 found in REV reads


# Load packages ####
# library("devtools")
# devtools::install_github("benjjneb/dada2")
# devtools::install_github("joey711/phyloseq")
# devtools::install_github("bryandmartin/corncob")
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")javascript:;
library(DESeq2)
library(DECIPHER)
library(decontam)
library(phangorn)
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(gridExtra)
library(ggpubr)
library(dada2); packageVersion("dada2")
library(stringr)
library(purrr)
library(corncob)
