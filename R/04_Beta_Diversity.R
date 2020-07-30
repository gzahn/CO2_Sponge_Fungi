# -----------------------------------------------------------------------------#
# Mediterranean sponge fungi - ocean acidification project 
# Beta-Diversity measures nd Ordinations
# Author: Geoffrey Zahn
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     patchwork v 1.0.1
#                     phyloseq v 1.30.0
#                     corncob v 0.1.0
#                     broom v 0.7.0
#                     purrr v 0.3.4
#                     vegan v 2.5.6
#                     ade4 v 1.7.15
# -----------------------------------------------------------------------------#

# packages ####
library(tidyverse); packageVersion("tidyverse")
library(patchwork); packageVersion("patchwork")
library(phyloseq); packageVersion("phyloseq")
library(corncob); packageVersion("corncob")
library(vegan); packageVersion("vegan")
library(broom); packageVersion("broom")
library(purrr); packageVersion("purrr")
library(ade4); packageVersion("ade4")

#functions
source("./R/bbdml_helper.R")
source("./R/palettes.R")


# IMPORT DATA ####
ps <- readRDS("./output/clean_phyloseq_object.RDS")

# remove "NA" Phylum taxa
ps <- subset_taxa(ps,!is.na(tax_table(ps)[,2]))

# relative abundance transformation
ps_ra <- ps %>% transform_sample_counts(function(x){x/sum(x)})

# MANTEL TEST AND MULTIPLE REGRESSION ON MATRICES ####
spatial.dist = vegdist(as.matrix(cbind(ps_ra@sam_data$Latitude, ps_ra@sam_data$Longitude)),method = "bray")
comm.dist = vegdist(as.matrix(ps_ra@otu_table),method = "bray")
man <- vegan::mantel(spatial.dist,comm.dist)
man
sink("./output/mantel_test_summary.txt")
print("Spatial vs community distance...")
man
sink(NULL)
# No apparent signifigance to spatial correlation and community distance



# Ordinations ####
ord <- ordinate(ps,method = "NMDS")
plot_ordination(ps,ord,color="Acidified")
