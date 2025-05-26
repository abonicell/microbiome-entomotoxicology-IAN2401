# load packages
suppressPackageStartupMessages({
  library("qiime2R")
  library("phyloseq")
  library("tidyverse")
  library("ggrepel")
  library("ggpubr")
  library("vegan")
  library("ranacapa")
  library("knitr")
  library("microbiomeMarker")
  library("DT")
})

# reprocess data to keep unique ASV IDs -----------------------------------
theme_set(theme_bw())

# remove contaminants
( physeq <- physeq %>%
    subset_taxa(
      Family!= "Mitochondria" | is.na(Family) #&
      # Order!="Chloroplast" | is.na(Class)
    ) %>%
    subset_taxa(
      # Family!= "Mitochondria" | is.na(Family) &
      Order!="Chloroplast" | is.na(Class)
    ) )# 4173 taxa

# prune sample to remove all with the sum of 0
ps_prune <-  prune_samples(sample_sums(physeq)!=0, physeq)
ps_prune

# check for missing data
any(is.na(ps_prune@otu_table))

# Normalize number of reads in each sample using median sequencing depth
total = median(sample_sums(ps_prune))
standf = function(x, t=total) round(t * (x / sum(x)))
ps = transform_sample_counts(ps_prune, standf)

# check for missing data
any(is.na(ps@otu_table))

##-----------------------------------------------------------------------------
# Liner discriminant analysis (LDA) effect size (LEFSe) analysis
mm_lefse_treat <- run_lefse(
  ps,
  group = "treatment",
  wilcoxon_cutoff = 0.05,
  kw_cutoff = 0.05,
  lda_cutoff = 2
)

# print results
mm_lefse_stage_tab <- marker_table(mm_lefse_treat)
mm_lefse_stage_tab %>% datatable()
# dot plot
plot_ef_dot(mm_lefse_treat) 

##-----------------------------------------------------------------------------
# Liner discriminant analysis (LDA) effect size (LEFSe) analysis
mm_lefse_stage <- run_lefse(
  ps,
  group = "stage",
  wilcoxon_cutoff = 0.05,
  kw_cutoff = 0.05,
  lda_cutoff = 2
)

# print results
mm_lefse_stage_tab <- data.frame(marker_table(mm_lefse_stage))

write.csv(mm_lefse_stage_tab,
          '/Users/andreabonicelli/Documents/GitHub/microbiome-entomotoxicology-IAN2401/tables/supplementary Table 2.csv')
          
# dot plot
plot_ef_dot(mm_lefse_stage) 

ggsave('/Users/andreabonicelli/Documents/GitHub/microbiome-entomotoxicology-IAN2401/figures/Figure 6.pdf',
       width = 7.5, height = 6)
