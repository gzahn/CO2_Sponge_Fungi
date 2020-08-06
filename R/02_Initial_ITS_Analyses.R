# -----------------------------------------------------------------------------#
# Mediterranean sponge fungi - ocean acidification project 
# Initial exploratory analyses of processed ITS reads
# Author: Geoffrey Zahn
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     phyloseq v 1.30.0
#                     purrr v 0.3.4
#                     vegan v 2.5.6
#                     ggpubr v 0.4.0
#                     patchwork v 1.0.1
# -----------------------------------------------------------------------------#


# PACKAGES, SCRIPTS, AND SETUP ####
library(phyloseq); packageVersion("phyloseq")
library(tidyverse); packageVersion("tidyverse")
library(vegan); packageVersion("vegan")
library(purrr); packageVersion("purrr")
library(ggpubr); packageVersion("ggpubr")
library(patchwork); packageVersion("patchwork")

source("./R/palettes.R")
source("./R/plot_bar2.R")

#Set ggplot theme
theme_set(theme_bw())

# IMPORT DATA ####
ps <- readRDS("./output/clean_phyloseq_object.RDS")

# Investigate UNITE assignments at each taxon level
ps_sp <- ps
phy <- !is.na(tax_table(ps_sp)[,2])
cla <- !is.na(tax_table(ps_sp)[,3])
ord <- !is.na(tax_table(ps_sp)[,4])
fam <- !is.na(tax_table(ps_sp)[,5])
gen <- !is.na(tax_table(ps_sp)[,6])
spp <- !is.na(tax_table(ps_sp)[,7])
assignments_sponge <- data.frame(Phylum=phy, Class=cla,Order=ord,Family=fam,Genus=gen,Species=spp)

assignments_sponge %>% pivot_longer(1:6) %>% mutate(name=factor(name,levels = c("Phylum","Class","Order","Family","Genus","Species"))) %>%
ggplot(aes(x=name,fill=value)) + geom_bar() + scale_fill_manual(values=c("Gray","Black")) +
  labs(x="Taxonomic level",y="Count",fill="Unambiguous\nassignment")
ggsave("./output/figs/UNITE_Taxonomic_Assignment_Efficiency_at_Each_Taxonomic_Rank.png",dpi=300)
rm(phy,cla,ord,fam,gen,spp,assignments_sponge,ps_sp)

# remove "NA" Phylum taxa
ps <- subset_taxa(ps,!is.na(tax_table(ps)[,2]))

# Look at available metadata
glimpse(sample_data(ps))

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

# Merge samples for plotting ####
# merge based on sponge species and acidification
newvar <- paste(ps@sam_data$Sponge_Species,ps@sam_data$Acidified,sep="_")
ps@sam_data$newvar <- newvar
ps_merged <- merge_samples(ps,newvar)

# repair metadata
ps_merged@sam_data$Sponge_Species <- unlist(map(str_split(sample_names(ps_merged),"_"),1))
ps_merged@sam_data$Acidified <- unlist(map(str_split(sample_names(ps_merged),"_"),2))

# Quick heatmaps ####
# phylum
p <- ps_merged %>% tax_glom("Phylum") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_heatmap(taxa.label = "Phylum",sample.label = "Sponge_Species") + 
  facet_grid(~Acidified,scales="free") + theme(strip.background = element_blank(), strip.text = element_text(face="bold")) + labs(x="Sponge species")
p$scales$scales[[3]]$name <- "Relative\nabundance"
p
ggsave("./output/figs/Heatmap_Phylum_RelAbund_by_Sponge_Species.png",dpi=300)
# class
p <- ps_merged %>% tax_glom("Class") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_heatmap(taxa.label = "Class",sample.label = "Sponge_Species") + 
  facet_grid(~Acidified,scales="free") + theme(strip.background = element_blank(), strip.text = element_text(face="bold")) + labs(x="Sponge species")
p$scales$scales[[3]]$name <- "Relative\nabundance"
p
ggsave("./output/figs/Heatmap_Class_RelAbund_by_Sponge_Species.png",dpi=300)
# order
p <- ps_merged %>% tax_glom("Order") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_heatmap(taxa.label = "Order",sample.label = "Sponge_Species") + 
  facet_grid(~Acidified,scales="free") + theme(strip.background = element_blank(), strip.text = element_text(face="bold")) + labs(x="Sponge species")
p$scales$scales[[3]]$name <- "Relative\nabundance"
p
ggsave("./output/figs/Heatmap_Order_RelAbund_by_Sponge_Species.png",dpi=300)
# alternate view
p <- ps_merged %>% tax_glom("Phylum") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_heatmap(taxa.label = "Phylum",sample.label = "Acidified") + 
  facet_grid(~Sponge_Species,scales="free") + theme(strip.background = element_blank(), strip.text = element_text(face="bold")) + labs(x="Sponge species")
p$scales$scales[[3]]$name <- "Relative\nabundance"
p
ggsave("./output/figs/Heatmap_Phylum_Alternate_Layout.png",dpi=300)


# Diversity BarPlots ####
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
    Richness ~ pH + Sponge_Species)
summary(mod1)

mod2 <- glm(data = meta2, 
            Shannon ~ pH + Sponge_Species)
summary(mod2)


# ... Does not seem to be an effect for alpha diversity based on sponge species or acidification

# look at nestedness of species within sampling site
table(meta2$Sponge_Species,meta2$Sampling_Site)

# ANOVA of diversity by site
mod3 <- aov(data = meta2,
            Shannon ~ Sampling_Site * Sponge_Species)
summary(mod3)

mod3b <- aov(data = meta2,
             Richness ~ Sampling_Site * Sponge_Species)
summary(mod3b)


sink("./output/ANOVA_models_for_alpha-diversity.txt")
mod3
summary(mod3)
mod3b
summary(mod3b)
sink(NULL)

# look at cave vs non-cave sites
meta2 <- meta2 %>% mutate(LightStatus = case_when(Sampling_Site == "GM" ~ "Cave",
                                  Sampling_Site != "GM" ~ "Open"))
mod4 <- glm(data = meta2,
            Shannon ~ LightStatus * Sponge_Species)
summary(mod4) # not significant
mod5 <- glm(data = meta2,
            Shannon ~ LightStatus * Sponge_Species)
summary(mod5) # not significant


# Alpha diversity comparison tables and more plots ####
p1 <- meta %>% 
  ggdensity(x="Richness", add = "median",
            color = "Sponge_Species",fill="Sponge_Species",
            rug = TRUE, facet.by = "Acidified") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face="bold",size=14),
        legend.text = element_text(face = "italic"),
        axis.title = element_text(face="bold",size=16),
        legend.position = "none") +
  labs(x="Richness",y="Density",fill="Sponge species") +
  scale_colour_manual(guide="none",values=pal.discrete) +
  scale_fill_manual(values = pal.discrete)

p2 <- meta %>% 
  ggdensity(x="Shannon", add = "median",
            color = "Sponge_Species",fill="Sponge_Species",
            rug = TRUE, facet.by = "Acidified") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face="bold",size=14),
        legend.text = element_text(face = "italic"),
        axis.title = element_text(face="bold",size=16),
        legend.position = "bottom") +
  labs(x="Shannon diversity",y="Density",fill="Sponge species") +
  scale_colour_manual(guide="none",values=pal.discrete) +
  scale_fill_manual(values = pal.discrete)
p1/p2
ggsave("./output/figs/Alpha-Diversity_Distributions_by_Sponge_Species.png",dpi=300)  


# pH vs alpha diversity
comparisons <- list(c("6","7"),c("7","8"),c("6","8"))
p1 <- ggboxplot(meta, x = "pH", y = "Richness",
               color = "pH", palette =pal.discrete,
               add = "jitter") +
  stat_compare_means(comparisons = comparisons) +
  theme(axis.title = element_text(face="bold",size=14),
        axis.text = element_text(face="bold")) +
  labs(x="")
  

p2 <- ggboxplot(meta, x = "pH", y = "Shannon",
               color = "pH", palette =pal.discrete,
               add = "jitter") + 
  stat_compare_means(comparisons = comparisons) +
  labs(y="Shannon diversity") +
  theme(axis.title = element_text(face="bold",size=14),
        axis.text = element_text(face="bold"),
        legend.position = "none")

p1/p2
ggsave("./output/figs/Alpha-Diversity_Boxplots_over_pH.png",dpi=300,height = 8,width = 6)





######################### CORE MYCOBIOME ##################################
# find core mycobiome members of each sponge species

core_full <- ps_ra %>% core_members(prevalence = .25,detection = .01)
core_Crambe <- ps_ra %>% subset_samples(Sponge_Species == "Crambe") %>% core_members(prevalence = .25,detection = .01)
core_Chondrosia <- ps_ra %>% subset_samples(Sponge_Species == "Chondrosia") %>% core_members(prevalence = .25,detection = .01)
core_Petrosia <- ps_ra %>% subset_samples(Sponge_Species == "Petrosia") %>% core_members(prevalence = .25,detection = .01)
core_Chondrilla <- ps_ra %>% subset_samples(Sponge_Species == "Chondrilla") %>% core_members(prevalence = .25,detection = .01)

# save list in file
sink("./output/Core_Taxa_Names.txt")
print("Prevalence = 0.25; detection threshold = 0.01")
print("core_full")
core_full
print("core_Crambe")
core_Crambe
print("core_Chondrosia")
core_Chondrosia
print("core_Petrosia")
core_Petrosia
print("core_Chondrilla")
core_Chondrilla
sink(NULL)

# subset to those core taxa
ps_coretaxa <- ps %>% subset_taxa(taxa_names(ps) %in% unique(c(core_Chondrilla,core_Chondrosia,core_Crambe,core_Petrosia)))
ps_ra_coretaxa <- transform_sample_counts(ps_coretaxa,function(x){x/sum(x)})
ps_ra_coretaxa <- ps_ra_coretaxa %>% subset_samples(sample_sums(ps_ra_coretaxa) > 0)

sample_data(ps_ra_coretaxa)
otu_table(ps_ra_coretaxa)

core_adonis <- adonis(otu_table(ps_ra_coretaxa) ~ ps_ra_coretaxa@sam_data$Sampling_Site * ps_ra_coretaxa@sam_data$Sponge_Species,na.rm=TRUE)
sink("./output/permanova_comm-distance_vs_Site_and_Sponge-Species_only-core-taxa.txt")
core_adonis
sink(NULL)


# plot just core taxa
plot_composition(ps_ra_coretaxa,group_by = "Sponge_Species",sample.sort = "pH") + scale_fill_manual(values = pal.discrete) + 
  labs(fill="Core taxa") +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(face="bold.italic",size=12),
        axis.title = element_text(face="bold",size=14),
        legend.title = element_text(face="bold",size = 12)) + labs(y="Relative abundance")
ggsave("./output/Relative_Abundance_of_Only_Core_Taxa_by_SpongeSpecies.png",dpi=300,height = 6,width = 14)  
ps_ra_coretaxa@sam_data




# find core taxa for acidification status
core_Acidified <- ps_ra %>% subset_samples(Acidified == "Acidified") %>% core_members(prevalence = .25,detection = .01)
core_Control <- ps_ra %>% subset_samples(Acidified == "Control") %>% core_members(prevalence = .25,detection = .01)

# subset to those core taxa
ps_coretaxa <- ps %>% subset_taxa(taxa_names(ps) %in% unique(c(core_Acidified,core_Control)))
ps_ra_coretaxa <- transform_sample_counts(ps_coretaxa,function(x){x/sum(x)})

# plot just core taxa
plot_composition(ps_ra_coretaxa,group_by = "Acidified",sample.sort = "pH") + scale_fill_manual(values = pal.discrete) + 
  labs(fill="Core taxa") +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(face="bold.italic",size=12),
        axis.title = element_text(face="bold",size=14),
        legend.title = element_text(face="bold",size = 12)) + labs(y="Relative abundance")
ggsave("./output/Relative_Abundance_of_Only_Core_Taxa_by_Acidification.png",dpi=300,height = 6,width = 10)  


# Permanova based on just these core taxa!!!


# make names pretty: just "genus spp."
coretaxa <- taxa_names(ps_ra)[grep(pattern = "ASV_19:|ASV_23:|ASV_14:|ASV_26:|ASV_6:|ASV_18:|ASV_15:|ASV_4:", x = taxa_names(ps_ra))]
ps_temp <- subset_taxa(ps_ra,taxa_names(ps_ra) %in% coretaxa)

genus <- ps_temp %>%
  taxa_names() %>% str_remove(pattern = "_fam_Incertae_sedis") %>%
  str_split("_") %>% map_chr(6)
species <- ps_temp %>%
  taxa_names() %>%
  str_split("_") %>% map_chr(c(7))
species[6] <- "platani"

ps_ra_genus_spp <- ps_temp
taxa_names(ps_ra_genus_spp) <- paste(genus,species)

# Plot core heatmap
plot_core(ps_ra_genus_spp, plot.type = "heatmap",min.prevalence = .25,horizontal = TRUE) +
  theme(axis.text.x = element_text(face = "bold.italic",angle=60,hjust=1),
        axis.title = element_text(face="bold",size=14),
        legend.title = element_text(face="bold",size=12)) +
  labs(title="") + scale_fill_viridis_c()
ggsave("./output/figs/Core_Taxa_Heatmap.png",dpi=300)
