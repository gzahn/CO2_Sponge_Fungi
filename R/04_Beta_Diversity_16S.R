# -----------------------------------------------------------------------------#
# Mediterranean sponge bacteria - ocean acidification project 
# Beta-Diversity measures nd Ordinations
# Author: Geoffrey Zahn
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     patchwork v 1.0.1
#                     phyloseq v 1.30.0
#                     broom v 0.7.0
#                     purrr v 0.3.4
#                     vegan v 2.5.6
#                     ade4 v 1.7.15
#                     doParallel v 1.0.15
#                     microbiome v 1.8.0
#                     foreach v 1.5.0
# -----------------------------------------------------------------------------#

# packages ####
library(tidyverse); packageVersion("tidyverse")
library(patchwork); packageVersion("patchwork")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(broom); packageVersion("broom")
library(purrr); packageVersion("purrr")
library(ade4); packageVersion("ade4")
library(vegan); packageVersion("vegan")
library(doParallel); packageVersion("doParallel")
library(microbiome); packageVersion("microbiome")
library(foreach); packageVersion("foreach")


#functions
source("./R/bbdml_helper.R")
source("./R/palettes.R")

# IMPORT DATA ####
ps <- readRDS("./output/16S_clean_ps_object_w_tree.RDS")

# remove "NA" Phylum taxa
ps <- subset_taxa(ps,!is.na(tax_table(ps)[,2]))

# relative abundance transformation
ps_ra <- ps %>% transform_sample_counts(function(x){x/sum(x)}) %>% subset_samples(Sponge_Species != "Seawater")

# MANTEL TEST AND MULTIPLE REGRESSION ON MATRICES ####
spatial.dist = vegdist(as.matrix(cbind(ps_ra@sam_data$Latitude, ps_ra@sam_data$Longitude)),method = "bray")
comm.dist = vegdist(as.matrix(ps_ra@otu_table),method = "bray")
man <- vegan::mantel(spatial.dist,comm.dist)
man
sink("./output/16S_mantel_test_summary.txt")
print("Spatial vs community distance...")
man
sink(NULL)

# Beta-diversity distances and ordinations ####
unifrac.dist <- UniFrac(ps,weighted = TRUE,normalized = TRUE,parallel = TRUE)

glimpse(sample_data(ps))
ordu = ordinate(ps, "PCoA","unifrac", weighted=TRUE)
plot_ordination(ps, ordu, color="Sponge_Species") +
  geom_point(size=3,alpha=.5) + scale_color_manual(values=pal.discrete) +
  labs(caption = "MDS/PCoA on weighted-UniFrac distance")

ggsave("./output/figs/16S_W-Unifrac_Ordination_Plot_by_Sponge-Species.png",dpi=300)

names(ps@sam_data)
plot_ordination(ps, ordu, color="Sampling_Site",shape="Sponge_Species") +
  geom_point(size=5,alpha=.5) + scale_color_manual(values=pal.discrete) +
  labs(caption = "MDS/PCoA on weighted-UniFrac distance") + theme_minimal()
ggsave("./output/figs/16S_W-Unifrac_Ordination_Plot_by_Sponge-Species_and_Sample_Site.png",dpi=300)

# PERMANOVA ####
set.seed(123)
permanova <- vegan::adonis(otu_table(ps_ra) ~ ps_ra@sam_data$Sampling_Site * ps_ra@sam_data$Sponge_Species)
sink("./output/16S_permanova_comm-distance_vs_Site_and_Sponge-Species.txt")
permanova
sink(NULL)

# Beta-Dispersion
w <- betadiver(otu_table(ps_ra),"w")
w.disper <- betadisper(w,group = meta(ps_ra)$Sponge_Species)

png("./output/figs/16S_Beta_Dispersion_Plot_Sponge-Species.png")
plot(w.disper,main = "Beta-Dispersion")
dev.off()
