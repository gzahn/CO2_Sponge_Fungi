# -----------------------------------------------------------------------------#
# Mediterranean sponge fungi - ocean acidification project 
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
#                     microbiome v 1.8.0
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
# BiocManager::install("microbiome")
library(microbiome); packageVersion("microbiome")


#functions
source("./R/bbdml_helper.R")
source("./R/palettes.R")


# IMPORT DATA ####
ps <- readRDS("./output/clean_phyloseq_object.RDS")

# remove "NA" Phylum taxa
ps <- subset_taxa(ps,!is.na(tax_table(ps)[,2]))

#subset to Petrosia only
ps_pet <- subset_samples(ps, Sponge_Species == "Petrosia")

# relative abundance transformation
ps_ra <- ps %>% transform_sample_counts(function(x){x/sum(x)}) %>% subset_samples(Sponge_Species != "Seawater")
ps_pet_ra <- ps_pet %>% transform_sample_counts(function(x){x/sum(x)})

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
set.seed(123)
ord <- ordinate(ps_ra,method = "DCA",distance = "jaccard")
p1 <- plot_ordination(ps_ra,ord,color="Sampling_Site") + scale_color_manual(values=pal.discrete) + theme(legend.position = "top") +
  labs(color="Site")
p2 <- plot_ordination(ps_ra,ord,color="Sponge_Species") + scale_color_manual(values=pal.discrete) + 
  theme(legend.position = "top",
        legend.text = element_text(face="italic")) + labs(color="Sponge species")
p3 <- plot_ordination(ps_ra,ord,color="Acidified") + scale_color_manual(values=pal.discrete) + theme(legend.position = "top") +
  labs(color="Acidification status")

p1+p2+p3
ggsave("./output/figs/DCA_OrdinationPlots_Colored_Variously.png",height = 8,width = 12,dpi=300)


# Petrosia only
ordu = ordinate(ps_pet, "PCoA","jaccard", weighted=TRUE)
plot_ordination(ps_pet, ordu, color="Sampling_Site") +
  geom_point(size=3,alpha=.5) + scale_color_manual(values=pal.discrete) +
  labs(caption = "MDS/PCoA on weighted-UniFrac distance",
       subtitle = "Petrosia sp.") +
  theme(plot.subtitle = element_text(face = "italic"))
ggsave("./output/figs/ITS_Petrosia_W-Unifrac_Ordination_Plot_by_Sample_Site.png",dpi=300)

plot_ordination(ps_pet, ordu, color="Acidified") +
  geom_point(size=3,alpha=.5) + scale_color_manual(values=pal.discrete) +
  labs(caption = "MDS/PCoA on weighted-UniFrac distance",
       subtitle = "Petrosia sp.") +
  theme(plot.subtitle = element_text(face = "italic"))
ggsave("./output/figs/ITS_Petrosia_W-Unifrac_Ordination_Plot_by_Acidification.png",dpi=300)


# PERMANOVA ####
set.seed(123)
permanova <- vegan::adonis(otu_table(ps_ra) ~ ps_ra@sam_data$Sampling_Site * ps_ra@sam_data$Sponge_Species)
sink("./output/permanova_comm-distance_vs_Site_and_Sponge-Species.txt")
permanova
sink(NULL)

set.seed(123)
permanova <- vegan::adonis(otu_table(ps_pet_ra) ~ ps_pet_ra@sam_data$Sampling_Site * ps_pet_ra@sam_data$Acidified)
sink("./output/ITS_Petrosia_permanova_comm-distance_vs_Site_and_Acidification.txt")
permanova
sink(NULL)

# Beta-dispersion
names(meta(ps_ra))
w <- betadiver(otu_table(ps_ra),"w")
w.disper <- betadisper(w,group = meta(ps_ra)$Sponge_Species)

png("./output/figs/Beta_Dispersion_Plot_Sponge-Species.png")
plot(w.disper)
dev.off()

w <- betadiver(otu_table(ps_pet_ra),"w")
w.disper <- betadisper(w,group = meta(ps_pet_ra)$Sampling_Site)

png("./output/figs/ITS_Petrosia_Beta_Dispersion_Plot_Sampling_Site.png")
plot(w.disper)
dev.off()
