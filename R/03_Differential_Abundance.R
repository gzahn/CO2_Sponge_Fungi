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

#functions
source("./R/bbdml_helper.R")
source("./R/palettes.R")


# IMPORT DATA ####
ps <- readRDS("./output/clean_phyloseq_object.RDS")

# remove "NA" Phylum taxa
ps <- subset_taxa(ps,!is.na(tax_table(ps)[,2]))


# glom taxa at various levels
ps_order <- tax_glom(ps,"Order")
ps_family <- tax_glom(ps,"Family")
ps_genus <- tax_glom(ps,"Genus")

glimpse(sample_data(ps))

######################## Acidified vs Control #########################


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
                 color = "FireTreatment",
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


if(length(da_analysis$significant_models) > 0){
  
  
  bbdml_obj <- multi_bbdml(da_analysis,
                            ps_object = ps_order,
                            mu_predictor = "Acidified",
                            phi_predictor = "Acidified",
                            taxlevels = 2:7)
  length(fire_bbdml)
  
  plot_multi_bbdml(bbdml_list = fire_bbdml,
                   color = "FireTreatment",
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
ggsave("./output/figs/DA_Plot_Genus-Level_Acidification.png")


#################### Sponge_Species #######################

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


#####################################################################







