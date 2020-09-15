# -----------------------------------------------------------------------------#
# Mediterranean sponge bacteria - ocean acidification project 
# Differnetial abundance of taxa in sample groups
# Author: Geoffrey Zahn
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     patchwork v 1.0.1
#                     phyloseq v 1.30.0
#                     corncob v 0.1.0
#                     broom v 0.7.0
#                     purrr v 0.3.4
#                     vegan v 2.5.6
# -----------------------------------------------------------------------------#

# packages ####
library(tidyverse); packageVersion("tidyverse")
library(patchwork); packageVersion("patchwork")
library(phyloseq); packageVersion("phyloseq")
library(corncob); packageVersion("corncob")
library(vegan); packageVersion("vegan")
library(broom); packageVersion("broom")
library(purrr); packageVersion("purrr")
library(dada2)
#functions
source("./R/bbdml_helper.R")
source("./R/palettes.R")


# IMPORT DATA ####
ps <- readRDS("./output/16S_clean_ps_object_w_tree.RDS")
seqs_and_taxa <- readRDS("./output/16S_SequenceNames_and_Taxonomy.RDS")
seqs_and_taxa <- seqs_and_taxa %>% map_df(as.character)

# remove "NA" Phylum taxa
ps <- subset_taxa(ps,!is.na(tax_table(ps)[,2]))


# Check taxa names
identical(taxa_names(ps),seqs_and_taxa$Sequence)


# mutate(tax_table(ps),FullTaxonomy=paste(Phylum,Class,Order,Family,Genus,Species))

# Corncob Differential Abundance ####

# Sponge_Species
set.seed(123)
da_analysis <- differentialTest(formula = ~ Sponge_Species, #abundance
                                phi.formula = ~ Sponge_Species, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps,
                                fdr_cutoff = 0.05)


length(da_analysis$significant_models)

#add taxonomic names to models
# sig_taxa <- as.data.frame(tax_table(ps)[da_analysis$significant_taxa,])
# sig_taxa <- mutate(sig_taxa,FullTaxonomy=paste(Phylum,Class,Order,Family,Genus,Species))
# names(da_analysis$significant_models) <- sig_taxa$FullTaxonomy
# names(da_analysis$significant_taxa) <- sig_taxa$FullTaxonomy
# plot(da_analysis)


if(length(da_analysis$significant_models) > 0 & length(da_analysis$significant_models) < 35){ # spare the pointless computation of a hundred separate plots!
  
  
  bbdml_obj <- multi_bbdml(da_analysis,
                           ps_object = ps,
                           mu_predictor = "Sponge_Species",
                           phi_predictor = "Sponge_Species",
                           taxlevels = 2:7)
  length(bbdml_obj)
  
  plot_multi_bbdml(bbdml_list = bbdml_obj,
                   color = "pH",
                   pointsize = 3)
}

# Genus-level tax_glom to tame things a bit ####
ps_genus <- tax_glom(ps,"Genus")

set.seed(123)
da_analysis <- differentialTest(formula = ~ Sponge_Species, #abundance
                                phi.formula = ~ Sponge_Species, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_genus,
                                fdr_cutoff = 0.05)

plot(da_analysis) + theme(axis.text.x = element_text(face="italic"),
                          strip.background = element_blank()) + labs(caption = "Differential abundance of significant genera; Wald test, FDR < 0.05")
ggsave("./output/figs/16S_DA_Plot_Sponge-Species_Overview.png",dpi = 300,width = 20,height = 8)


if(length(da_analysis$significant_models) > 0){
  
  
  bbdml_obj <- multi_bbdml(da_analysis,
                           ps_object = ps_genus,
                           mu_predictor = "Sponge_Species",
                           phi_predictor = "Sponge_Species",
                           taxlevels = 2:7)
  length(bbdml_obj)
  
  plot_multi_bbdml(bbdml_list = bbdml_obj,
                   color = "Sponge_Species",
                   pointsize = 3)
}

# save model results
sink("./output/16S_DA_genus-level_significant_models.txt")
bbdml_obj
sink(NULL)

# plot ALL the bbdml individual models
bbdml_plot_1/bbdml_plot_2/bbdml_plot_3/bbdml_plot_4/bbdml_plot_5/bbdml_plot_6/bbdml_plot_7
ggsave("./output/figs/16S_DA_Taxa_Part1.png",dpi=300,height = 14,width = 12)
bbdml_plot_8/bbdml_plot_9/bbdml_plot_10/bbdml_plot_11/bbdml_plot_12/bbdml_plot_13/bbdml_plot_14
ggsave("./output/figs/16S_DA_Taxa_Part2.png",dpi=300,height = 14,width = 12)
bbdml_plot_15/bbdml_plot_16/bbdml_plot_17/bbdml_plot_18/bbdml_plot_19/bbdml_plot_20/bbdml_plot_21
ggsave("./output/figs/16S_DA_Taxa_Part3.png",dpi=300,height = 14,width = 12)
bbdml_plot_22/bbdml_plot_23/bbdml_plot_24/bbdml_plot_25/bbdml_plot_26/bbdml_plot_27/bbdml_plot_28
ggsave("./output/figs/16S_DA_Taxa_Part4.png",dpi=300,height = 14,width = 12)
bbdml_plot_29/bbdml_plot_30/bbdml_plot_31/bbdml_plot_32/bbdml_plot_33/bbdml_plot_24/bbdml_plot_25
ggsave("./output/figs/16S_DA_Taxa_Part5.png",dpi=300,height = 14,width = 12)




# Subset to Petrosia and run diffabund based on site and acidification ####


ps_pet <- ps %>% subset_samples(Sponge_Species == "Petrosia") %>% tax_glom("Genus")
ps_pet@sam_data

set.seed(123)
da_analysis <- differentialTest(formula = ~ Sampling_Site, #abundance
                                phi.formula = ~ Sampling_Site, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_pet,
                                fdr_cutoff = 0.05)

length(da_analysis$significant_models)

#add taxonomic names to models
sig_taxa <- as.data.frame(tax_table(ps)[da_analysis$significant_taxa,])
sig_taxa <- mutate(sig_taxa,FullTaxonomy=paste(Phylum,Class,Order,Family,Genus,Species))
names(da_analysis$significant_models) <- sig_taxa$FullTaxonomy
names(da_analysis$significant_taxa) <- sig_taxa$FullTaxonomy


plot(da_analysis) + theme(axis.text.x = element_text(face="italic"),
                          strip.background = element_blank()) + labs(caption = "Differential abundance of significant genera; Wald test, FDR < 0.05")
ggsave("./output/figs/16S_DA_Plot_Petrosia_by_Site_Overview.png",dpi = 300,width = 20,height = 8)

# Acidification

set.seed(123)
da_analysis <- differentialTest(formula = ~ Acidified, #abundance
                                phi.formula = ~ Acidified, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_pet,
                                fdr_cutoff = 0.05)

length(da_analysis$significant_models)

#add taxonomic names to models
# sig_taxa <- as.data.frame(tax_table(ps)[da_analysis$significant_taxa,])
# sig_taxa <- mutate(sig_taxa,FullTaxonomy=paste(Phylum,Class,Order,Family,Genus,Species))
# names(da_analysis$significant_models) <- sig_taxa$FullTaxonomy
# names(da_analysis$significant_taxa) <- sig_taxa$FullTaxonomy


plot(da_analysis) + theme(axis.text.x = element_text(face="italic"),
                          strip.background = element_blank()) + labs(caption = "Differential abundance of significant genera; Wald test, FDR < 0.05")
ggsave("./output/figs/16S_DA_Plot_Petrosia_by_Acidification_Overview.png",dpi = 300,width = 12,height = 8)
