# -----------------------------------------------------------------------------#
# Mediterranean sponge fungi - ocean acidification project 
# Initial exploratory analyses of processed ITS reads
# Author: Geoffrey Zahn
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     phyloseq v 1.30.0
#                     purrr v 0.3.4
#                     vegan v 2.5.6
# -----------------------------------------------------------------------------#


# PACKAGES, SCRIPTS, AND SETUP ####
library(phyloseq); packageVersion("phyloseq")
library(tidyverse); packageVersion("tidyverse")
library(vegan); packageVersion("vegan")
library(purrr); packageVersion("purrr")

source("./R/palettes.R")
source("./R/plot_bar2.R")

#Set ggplot theme
theme_set(theme_bw())

# IMPORT DATA ####
ps <- readRDS("./output/clean_phyloseq_object.RDS")

# remove "NA" Phylum taxa
ps <- subset_taxa(ps,!is.na(tax_table(ps)[,2]))

# Look at available metadata
glimpse(sample_data(ps))

# Clean up taxonomy names
for(i in 1:7){
  tax_table(ps)[,i] <- str_remove(tax_table(ps)[,i],".__")  
}

# relative abundance
ps_ra <- ps %>% transform_sample_counts(function(x){x/sum(x)})

# Quick glance at alpha diversity ####
plot_richness(ps,x="Acidified",measures = "Shannon") + 
  facet_grid(~Sponge_Species) + labs(y="Shannon diversity")
ggsave("./output/figs/Shannon_diversity_dotplot_by_Species_and_Acidification.png",dpi=300)

# Calculate alpha diversity measures and add to metadata
ps@sam_data$Shannon <- estimate_richness(ps, measures="Shannon")$Shannon
ps@sam_data$Richness <- specnumber(otu_table(ps))

# save as dataframe object for easy modeling
meta <- as(sample_data(ps),"data.frame")
glimpse(meta)

# Export crosstab of acidification vs species
sink("./output/Sample_N_for_groups_Acidification_and_SpongeSpecies.txt")
table(ps@sam_data$Acidified,ps@sam_data$Sponge_Species)
sink(NULL)

# Relative abundance Bar Plots ####

# merge based on sponge species and acidification
newvar <- paste(ps@sam_data$Sponge_Species,ps@sam_data$Acidified,sep="_")
ps@sam_data$newvar <- newvar
ps_merged <- merge_samples(ps,newvar)

# repair metadata
ps_merged@sam_data$Sponge_Species <- unlist(map(str_split(sample_names(ps_merged),"_"),1))
ps_merged@sam_data$Acidified <- unlist(map(str_split(sample_names(ps_merged),"_"),2))

# Basic plots of diversity for overall data
# stacked boxplots x-axis=Acidification, y-axis relative-abundance

# Phylum
ps_merged %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar2(fill="Phylum",x="Acidified") + scale_fill_manual(values = pal.discrete) + 
  labs(y="Relative abundance",x="Acidification") +
  theme(axis.text.x = element_text(angle = 60,hjust=1),
        axis.title = element_text(face="bold",size=16),
        axis.text = element_text(face="bold",size=12),
        strip.background = element_blank(),
        strip.text = element_text(size=12,face="bold.italic")) +
    facet_wrap(~Sponge_Species,scales = "free")
ggsave("./output/figs/Phylum_Diversity_BarChart_by_Acidification.png",dpi=300)

# Class
ps_merged %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar2(fill="Class",x="Acidified") + scale_fill_manual(values = pal.discrete) + 
  labs(y="Relative abundance",x="Acidification") +
  theme(axis.text.x = element_text(angle = 60,hjust=1),
        axis.title = element_text(face="bold",size=16),
        axis.text = element_text(face="bold",size=12),
        strip.background = element_blank(),
        strip.text = element_text(size=12,face="bold.italic")) +
  facet_wrap(~Sponge_Species,scales = "free") 
ggsave("./output/figs/Class_Diversity_BarChart_by_Acidification.png",dpi=300,width = 10, height = 8)


# merge based on site only
ps_merged <- merge_samples(ps,"Sampling_Site")

# Phylum
ps_merged %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar2(fill="Phylum") + scale_fill_manual(values = pal.discrete) + 
  labs(y="Relative abundance",x="Site") +
  theme(axis.text.x = element_text(angle = 60,hjust=1),
        axis.title = element_text(face="bold",size=16),
        axis.text = element_text(face="bold",size=12),
        strip.background = element_blank(),
        strip.text = element_text(size=12,face="bold.italic"))
ggsave("./output/figs/Phylum_Diversity_BarChart_by_Site.png",dpi=300)

# Class
ps_merged %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar2(fill="Class") + scale_fill_manual(values = pal.discrete) + 
  labs(y="Relative abundance",x="Site") +
  theme(axis.text.x = element_text(angle = 60,hjust=1),
        axis.title = element_text(face="bold",size=16),
        axis.text = element_text(face="bold",size=12),
        strip.background = element_blank(),
        strip.text = element_text(size=12,face="bold.italic"))
ggsave("./output/figs/Class_Diversity_BarChart_by_Site.png",dpi=300)



# Model alpha diversity ####
# as a function of acidification, species, bleaching status

# Consider, first, the data without Seawater samples...not sure how to treat those...
meta2 <- meta[meta$Sponge_Species != "Seawater",]
mod1 <- glm(data = meta2, 
    Richness ~ Acidified * Sponge_Species)
summary(mod1)

mod2 <- glm(data = meta2, 
            Shannon ~ Acidified * Sponge_Species)
summary(mod2)


# ... Does not seem to be an effect for alpha diversity based on sponge species or acidification




