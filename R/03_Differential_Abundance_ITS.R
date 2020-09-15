# -----------------------------------------------------------------------------#
# Mediterranean sponge fungi - ocean acidification project 
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
ps <- readRDS("./output/clean_phyloseq_object.RDS")
seqs_and_taxa <- readRDS("./output/SequenceNames_and_Taxonomy.RDS")

# remove "NA" Phylum taxa
ps <- subset_taxa(ps,!is.na(tax_table(ps)[,2]))

# Clean up taxonomy names
for(i in 1:7){
  tax_table(ps)[,i] <- str_remove(tax_table(ps)[,i],".__")  
}


# glom taxa at various levels
ps_order <- tax_glom(ps,"Order")
ps_family <- tax_glom(ps,"Family")
ps_genus <- tax_glom(ps,"Genus")

glimpse(sample_data(ps))

######################## Acidified vs Control #########################

# ASV-Level

set.seed(123)
da_analysis <- differentialTest(formula = ~ Acidified, #abundance
                                phi.formula = ~ Acidified, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps,
                                fdr_cutoff = 0.05)

plot(da_analysis)

if(length(da_analysis$significant_models) > 0){
  
  
  bbdml_obj <- multi_bbdml(da_analysis,
                           ps_object = ps,
                           mu_predictor = "Acidified",
                           phi_predictor = "Acidified",
                           taxlevels = 2:7)
  length(bbdml_obj)
  
  plot_multi_bbdml(bbdml_list = bbdml_obj,
                   color = "Acidified",
                   pointsize = 3)
}



# Order-level ##

set.seed(123)
da_analysis <- differentialTest(formula = ~ Acidified, #abundance
                                phi.formula = ~ Acidified, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_order,
                                fdr_cutoff = 0.05)

if(length(da_analysis$significant_models) > 0){


  bbdml_obj <- multi_bbdml(da_analysis,
                          ps_object = ps_order,
                          mu_predictor = "Acidified",
                          phi_predictor = "Acidified",
                          taxlevels = 2:7)
length(fire_bbdml)

plot_multi_bbdml(bbdml_list = fire_bbdml,
                 color = "Acidified",
                 pointsize = 3)
}


# family-level ##
set.seed(123)
da_analysis <- differentialTest(formula = ~ Acidified, #abundance
                                phi.formula = ~ Acidified, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_family,
                                fdr_cutoff = 0.05)

plot(da_analysis)
if(length(da_analysis$significant_models) > 0){
  
  
  bbdml_obj <- multi_bbdml(da_analysis,
                            ps_object = ps_order,
                            mu_predictor = "Acidified",
                            phi_predictor = "Acidified",
                            taxlevels = 2:7)
  length(fire_bbdml)
  
  plot_multi_bbdml(bbdml_list = fire_bbdml,
                   color = "Acidified",
                   pointsize = 3)
}

# genus-level ##
set.seed(123)
da_analysis <- differentialTest(formula = ~ Acidified, #abundance
                                phi.formula = ~ Acidified, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_genus,
                                fdr_cutoff = 0.05)

plot(da_analysis)

if(length(da_analysis$significant_models) > 0){
  
  
  bbdml_obj <- multi_bbdml(da_analysis,
                            ps_object = ps_order,
                            mu_predictor = "Acidified",
                            phi_predictor = "Acidified",
                            taxlevels = 2:7)
  length(fire_bbdml)
  
  plot_multi_bbdml(bbdml_list = fire_bbdml,
                   color = "Acidified",
                   pointsize = 3)
}
bbdml_plot_1
ggsave("./output/figs/DA_Plot_Genus-Level_Acidification.png",height = 6,width = 12,dpi=300)


#################### Sponge_Species #######################


# ASV-Level

set.seed(123)
da_analysis <- differentialTest(formula = ~ Sponge_Species, #abundance
                                phi.formula = ~ Sponge_Species, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps,
                                fdr_cutoff = 0.05)

plot(da_analysis)
ggsave("./output/figs/DA_Plot_Sponge-Species_Overview.png",dpi=300,height = 4, width = 22)


da_analysis$significant_taxa

bbdml(formula = FungalASV_9: Ascomycota_Eurotiomycetes_Eurotiales_Aspergillaceae_Aspergillus_baarnensis ~ Sponge_Species,
      data = ps_genus)

if(length(da_analysis$significant_models) > 0){
  
  bbdml_obj <- multi_bbdml(da_analysis,
                           ps_object = ps,
                           mu_predictor = "Sponge_Species",
                           phi_predictor = "Sponge_Species",
                           taxlevels = 2:7)
  length(bbdml_obj)
  
  plot_multi_bbdml(bbdml_list = bbdml_obj,
                   color = "Sponge_Species",
                   pointsize = 3)
}

bbdml_plot_1 <- bbdml_plot_1 + theme(plot.title = element_text(face="italic"),axis.text.x = element_text(size=6))
bbdml_plot_2 <- bbdml_plot_3 + theme(plot.title = element_text(face="italic"),axis.text.x = element_text(size=6))
bbdml_plot_3 <- bbdml_plot_3 + theme(plot.title = element_text(face="italic"),axis.text.x = element_text(size=6))

bbdml_plot_1 / bbdml_plot_2 / bbdml_plot_3
ggsave("./output/figs/DA_Plot_Sponge-Species_Individual_Taxa.png",dpi=300,height = 8,width = 14)
# Order-level ##

set.seed(123)
da_analysis <- differentialTest(formula = ~ Sponge_Species, #abundance
                                phi.formula = ~ Sponge_Species, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_order,
                                fdr_cutoff = 0.05)

plot(da_analysis)
ggsave("./output/figs/DA_Plot_Order-Level_Sponge-Species_Overall_Plot.png",width = 14,height = 8,dpi = 300)


if(length(da_analysis$significant_models) > 0){
  
  
  bbdml_obj <- multi_bbdml(da_analysis,
                           ps_object = ps_order,
                           mu_predictor = "Sponge_Species",
                           phi_predictor = "Sponge_Species",
                           taxlevels = 2:7)
  
  plot_multi_bbdml(bbdml_list = bbdml_obj,
                   color = "Sponge_Species",
                   pointsize = 3)
}

length(bbdml_obj)
bbdml_plot_1 + bbdml_plot_2 + bbdml_plot_3 + bbdml_plot_4 + bbdml_plot_5
ggsave("./output/figs/DA_Plot_Order-Level_Sponge-Species.png",dpi=300,height = 8,width = 14)

# family-level ##
set.seed(123)
da_analysis <- differentialTest(formula = ~ Sponge_Species, #abundance
                                phi.formula = ~ Sponge_Species, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_family,
                                fdr_cutoff = 0.05)
plot(da_analysis)
ggsave("./output/figs/DA_Plot_Family-Level_Sponge-Species_Overall_Plot.png",width = 14,height = 8,dpi = 300)

length(da_analysis$significant_models)

if(length(da_analysis$significant_models) > 0){
  bbdml_obj <- multi_bbdml(da_analysis,
                           ps_object = ps_family,
                           mu_predictor = "Sponge_Species",
                           phi_predictor = "Sponge_Species",
                           taxlevels = 2:7)
  length(fire_bbdml)
  
  plot_multi_bbdml(bbdml_list = bbdml_obj,
                   color = "Sponge_Species",
                   pointsize = 3)
}

# genus-level ##
set.seed(123)
da_analysis <- differentialTest(formula = ~ Sponge_Species, #abundance
                                phi.formula = ~ Sponge_Species, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_genus,
                                fdr_cutoff = 0.05)

plot(da_analysis)
ggsave("./output/figs/DA_Plot_Genus-Level_Sponge-Species_Overall_Plot.png",width = 14,height = 8,dpi = 300)



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

(bbdml_plot_1 + bbdml_plot_2 + bbdml_plot_3 + bbdml_plot_4 + bbdml_plot_5)
ggsave("./output/figs/DA_Plot_Genus-Level_Sponge-Species.png", height = 6, width = 14,dpi = 300)


###################### pH level #############################

# ASV-Level

set.seed(123)
da_analysis <- differentialTest(formula = ~ pH, #abundance
                                phi.formula = ~ pH, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps,
                                fdr_cutoff = 0.05)

plot(da_analysis)

if(length(da_analysis$significant_models) > 0){
  
  
  bbdml_obj <- multi_bbdml(da_analysis,
                           ps_object = ps,
                           mu_predictor = "pH",
                           phi_predictor = "pH",
                           taxlevels = 2:7)
  length(bbdml_obj)
  
  plot_multi_bbdml(bbdml_list = bbdml_obj,
                   color = "pH",
                   pointsize = 3)
}

# genus-level ##


set.seed(123)
da_analysis <- differentialTest(formula = ~ factor(pH), #abundance
                                phi.formula = ~ factor(pH), #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_genus,
                                fdr_cutoff = 0.05)

# plot(da_analysis)
# ggsave("./output/figs/DA_Plot_Genus-Level_pH_Overall_Plot.png",width = 12,height = 6,dpi = 300)



if(length(da_analysis$significant_models) > 0){
  
  
  bbdml_obj <- multi_bbdml(da_analysis,
                           ps_object = ps_genus,
                           mu_predictor = "pH",
                           phi_predictor = "pH",
                           taxlevels = 2:7)
  length(bbdml_obj)
  
  plot_multi_bbdml(bbdml_list = bbdml_obj,
                   color = "pH",
                   pointsize = 3)
}



# Subset to Petrosia ... DA analysis based on site and acidification ####

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
ggsave("./output/figs/ITS_DA_Plot_Petrosia_by_Site_Overview.png",dpi = 300,width = 20,height = 6)

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


plot(da_analysis) + theme(axis.text.x = element_text(face="italic"),
                          strip.background = element_blank()) + labs(caption = "Differential abundance of significant genera; Wald test, FDR < 0.05")
ggsave("./output/figs/ITS_DA_Plot_Petrosia_by_Acidification_Overview.png",dpi = 300,width = 12,height = 4)





