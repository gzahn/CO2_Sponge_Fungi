# -----------------------------------------------------------------------------#
# Mediterranean sponge bacteria - ocean acidification project 
# Building and adding a phylogeny to the cleaned phyloseq object
# Author: Geoffrey Zahn
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     vegan v 2.5.6
#                     phangorn v 2.2.5
#                     phyloseq v 1.30.0
#                     msa v 1.18.0
#                     ape v 5.4
#                     seqinr v 3.6.1
# -----------------------------------------------------------------------------#


# Packages and functions ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(phangorn); packageVersion("phangorn")
library(msa); packageVersion("msa")
library(ape); packageVersion("ape")
library(seqinr); packageVersion("seqinr")
theme_set(theme_bw())



# Read in phyloseq object ####
ps <- readRDS("./output/16S_clean_phyloseq_object.RDS")


seqs <- rownames(tax_table(ps))
names(seqs) <- paste0("ASV_",1:length(seqs)) # This propagates to the tip labels of the tree

# Multiple sequence alignment  ####
alignment <- msa(seqs,method = "Muscle", type = "dna",verbose = TRUE,order = "input",maxiters = 10)

# save progress 
saveRDS(alignment,"./output/trees/16S_dna_alignment_muscle.RDS")
# writeXStringSet(DNAStringSet(seqs),filepath = "./output/trees/16S_Seqs.fasta",format = "fasta")

# re-load point
# alignment <- readRDS("./output/trees/Chagos_16S_dna_alignment_muscle.RDS")

# Convert to phangorn format
phang.align = as.phyDat(alignment, type = "DNA")
# write.phyDat(phang.align,"./output/trees/16S_dna_alignment_muscle.nex",format="nexus")

# distance max likelihood
dm <- dist.ml(phang.align)

#save
saveRDS(dm,"./output/trees/16S_ML_Distance.RDS")

# Initial neighbor-joining tree ####
treeNJ <- NJ(dm) # Note, tip order != sequence order
treeNJ$edge

treeNJ$tip.label <- seqs
#save
saveRDS(treeNJ, "./output/trees/16S_treeNJ.RDS")
# re-load
# treeNJ <- readRDS("./output/trees/treeNJ.RDS")

# Estimate model parameters ####
fit = pml(treeNJ, data=phang.align)
#save
saveRDS(fit,"./output/trees/16S_fit_treeNJ.RDS")
# re-load
# fit <- readRDS("./output/trees/fit_treeNJ.RDS")

# Likelihood of tree ####
fitJC <- optim.pml(fit, TRUE)
# save
saveRDS(fitJC, "./output/trees/16S_tree_fitJC.RDS") # This is the new tree using optim.pml
write.tree(fitJC$tree, file = "./output/trees/16S_tree_JC.nwk")
# reload point
# fitJC <- readRDS("./output/trees/16S_tree_fitJC.RDS")

names(fitJC$tree$tip.label) <- NULL
identical(fitJC$tree$tip.label, taxa_names(ps))

# names(taxa_names(ps)) <- names(fitJC$tree$tip.label)




# add tree to phyloseq object ####
ps2 <- phyloseq(tax_table(tax_table(ps)),
               otu_table(otu_table(ps)),
               sample_data(sample_data(ps)),
               phy_tree(fitJC$tree))






ps <- ps2

# tidy sampleID names
ps@sam_data$SampleID <- str_remove(ps@sam_data$SampleID,"-16S")
ps@sam_data$SampleID <- str_replace_all(ps@sam_data$SampleID,pattern = "-",replacement = "_")


# Save updated phyloseq object with tree
saveRDS(ps, "./output/16S_clean_ps_object_w_tree.RDS")
