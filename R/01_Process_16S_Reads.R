# -----------------------------------------------------------------------------#
# Mediterranean sponge bacteria - ocean acidification project 
# Processing cutadapt-trimmed Reads
# Author: Geoffrey Zahn
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     dada2 v 1.14.1
#                     decontam v 1.6.0
#                     phyloseq v 1.30.0
#                     purrr v 0.3.4
#                     Biostrings 2.54.0
# -----------------------------------------------------------------------------#


# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")

source("./R/palettes.R")
source("./R/plot_bar2.R")


#################################################################################
#                               Main workflow                                   #
# Filter and trim, denoise, sample inferrence, chimera and contaminant removal, # 
# taxonomic assignment, phylogenetic tree, combine sequence table and metadata  #
#################################################################################

# PARSE FILE PATHS ####

# File parsing - For this, we will use only the forward illumina reads - make sure to move fwd reads into their own directory for simplest processing
path <- "./seqs/16S/fastqs" # CHANGE to the directory containing your demultiplexed fastq files
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present
fns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "L001_R1_001.fastq.gz"))
rns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "L001_R2_001.fastq.gz"))


A <- unlist(map(strsplit(basename(fns), "-16S"), 1))
B <- unlist(map(strsplit(A, "DNA-"), 2))

sample.names <- B 
rm(list = c("A","B"))

# visualize a couple of fwd read quality profiles to help select reasonable filtration parameters
plotQualityProfile(fns[1:4])
plotQualityProfile(rns[1:4])

# FILTER AND TRIM ####
filts_f <- file.path(path, "filtered", paste0(sample.names, "_FWD_filt.fastq.gz"))
filts_r <- file.path(path, "filtered", paste0(sample.names, "_REV_filt.fastq.gz"))

out <- filterAndTrim(fns, filts_f, rns, filts_r,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, truncLen = c(280,200),
                     compress=TRUE, multithread=4) # On Windows set multithread=FALSE
head(out)

# Some input samples had no reads pass the filter!


# sanity check  comparison of before and after filtration
# plotQualityProfile(c(fns[1:2],filts_f[1:2]))
# plotQualityProfile(c(rns[1:2],filts_r[1:2]))

# LEARN ERROR RATES ####
# Since some samples may have had zero reads pass QC, reassign filts
filts_f <- sort(list.files(filtpath, full.names = TRUE,pattern = "FWD"))
filts_r <- sort(list.files(filtpath, full.names = TRUE,pattern = "REV"))
length(fns);length(filts_f)
length(rns);length(filts_r) 

# learn errors
errF <- learnErrors(filts_f, multithread=TRUE, MAX_CONSIST = 20)
errR <- learnErrors(filts_r, multithread=TRUE, MAX_CONSIST = 20)

# sanity check for error model
# plotErrors(errF, nominalQ=TRUE)
# plotErrors(errR, nominalQ=TRUE)

# DEREPLICATION ####
derepF <- derepFastq(filts_f, verbose=TRUE)
derepR <- derepFastq(filts_r, verbose=TRUE)


# Name the derep-class objects by the sample names
# If some samples were removed (no reads passed QC), reassign sample.names
if(length(derepF) != length(sample.names)){
  sample.names <- unlist(map(strsplit(basename(filts_f), "_filt"), 1))
}


if(identical(unlist(map(strsplit(basename(filts_f), "FWD_filt"), 1)),unlist(map(strsplit(basename(filts_r), "REV_filt"), 1)))){
  names(derepF) <- sample.names
  names(derepR) <- sample.names
} else {
  stop("Make sure fwd and rev files are in same order!")
}  



# SAMPLE INFERRENCE ####
dadaFs <- dada(derepF, err=errF, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo")
dadaRs <- dada(derepR, err=errR, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo")



mergers <- mergePairs(dadaFs, filts_f, dadaRs, filts_r, verbose=TRUE)

# MAKE SEQUENCE TABLE ####
seqtab <- makeSequenceTable(mergers)

# REMOVE CHIMERAS ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# reassign "out" to remove any missing reads
out = out[as.data.frame(out)$reads.out > 0,]

# TRACK READS THROUGH PIPELINE ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track = as.data.frame(track)
track$total.loss.proportion = (track[,1]-track$nonchim)/track[,1]
head(track)

write.csv(track, file = "./output/16S_read_counts_at_each_step.csv", row.names = TRUE)


# Save intermediate seqtable object
saveRDS(seqtab.nochim, "./output/16S_seqtab.nochim.RDS")


# IMPORT METADATA ####
meta = read_delim("./data/Metadata_16S.txt",delim = "\t")
row.names(meta) <- as.character(meta$SampleID) %>% str_remove("DNA-") %>% str_remove("-16S")
# subset to samples that passed with >0 reads


identical(row.names(meta),row.names(seqtab.nochim))


# Remove all seqs with fewer than 100 nucleotides ####
keeper_esvs <- nchar(names(as.data.frame(seqtab.nochim))) > 99
seqtab.nochim <- seqtab.nochim[,keeper_esvs]

# ASSIGN TAXONOMY ####
taxa <- assignTaxonomy(seqtab.nochim, "./taxonomy/rdp_train_set_16.fa.gz", multithread=4)

# Save intermediate taxonomy file
saveRDS(taxa, file = "./output/16S_RDP_Taxonomy_from_dada2.RDS")

# add_species
taxa <- addSpecies(taxa, "./taxonomy/rdp_species_assignment_16.fa.gz")

# Save completed taxonomy file
saveRDS(taxa, file = "./output/16S_RDP_Taxonomy_from_dada2_sp.RDS")


# inspect taxonomy
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# re-load point
# seqtab.nochim <- readRDS("./output/dada2_seqtable.RDS")
# taxa <- readRDS("./output/RDP_Taxonomy_from_dada2.RDS")
# meta <- read_delim("./data/Metadata.csv",delim = ",")
# row.names(meta) <- as.character(meta$SampleID)


# Hand off to Phyloseq ####
otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
tax <- tax_table(taxa)
met <- sample_data(meta)
row.names(met) <- row.names(meta)


ps <- phyloseq(otu,met,tax)



# Find non-bacteria ####
ps_nonbact <- subset_taxa(ps, Kingdom != "Bacteria")

ps %>% transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar2(fill="Kingdom")
ggsave("./output/figs/16S_Kingdom-Level_Taxonomic_Proportions.png",dpi=300)

ps_nonbact %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar2(fill="Kingdom")


# REMOVE NON-BACTERIA, CHLOROPLASTS, MITOCHONDRIA, and empty samples/taxa ####
ps <- subset_taxa(ps, Kingdom == "Bacteria")
ps <- subset_taxa(ps,Class != "Chloroplast")
ps <- subset_taxa(ps, taxa_sums(ps) > 0)
ps <- subset_samples(ps, sample_sums(ps) > 0)

# Save DNA sequences apart from rownames (from subsetted ps object)
seqs <- taxa_names(ps)
seqs <- DNAStringSet(seqs)
saveRDS(seqs,"./output/16S_ASV_reference_sequences.RDS")


pretty_names <- paste("16S_ASV",1:length(taxa_names(ps)),":",
      tax_table(ps)[,2],
      tax_table(ps)[,3],
      tax_table(ps)[,4],
      tax_table(ps)[,5],
      tax_table(ps)[,6],
      tax_table(ps)[,7], sep="_") %>%
  str_remove("k__") %>% str_remove("p__") %>% str_remove("c__") %>% str_remove("o__") %>% str_remove("f__") %>% str_remove("g__") %>% str_remove("s__") %>%
  str_replace(pattern = "_:_",replacement = ": ")

df <- data.frame(TaxaName=pretty_names,Sequence=taxa_names(ps))
saveRDS(df,"./output/16S_SequenceNames_and_Taxonomy.RDS")

# Set Seawater as first level of Sponge_Species
ps@sam_data$Sponge_Species <- factor(ps@sam_data$Sponge_Species, 
                                     levels = c("Seawater","Chondrilla","Chondrosia","Crambe","Petrosia"))


# Save RDS object for Phyloseq
saveRDS(ps, file = "./output/16S_clean_phyloseq_object.RDS")


