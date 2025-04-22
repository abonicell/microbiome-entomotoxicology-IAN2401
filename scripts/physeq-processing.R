# qiime2 parameters
# --p-trim-left-f 20 
# --p-trunc-len-f 200
# --p-trunc-len-r 260

# load packages -----------------------------------------------------------
suppressPackageStartupMessages({
library("qiime2R") 
library("phyloseq")
library("tidyverse")
library("ggpubr")
library("vegan")
library("ranacapa")
library("microViz")
library("ggsci")
library("tidyr")
library("DT")
})

theme_set(theme_bw())

# read other files and merge into phyloseq object
taxonomy<-qiime2R::read_qza("/Users/andreabonicelli/Documents/GitHub/microbiome-entomotoxicology-IAN2401/data/taxonomy-silva-v4.qza")
head(taxonomy$data);
qiime2R::parse_taxonomy(taxonomy$data)

( metadata <- readr::read_tsv("/Users/andreabonicelli/Documents/GitHub/microbiome-entomotoxicology-IAN2401/data/metadata.tsv") )

( physeq<-qiime2R::qza_to_phyloseq(
  features="/Users/andreabonicelli/Documents/GitHub/microbiome-entomotoxicology-IAN2401/data/table.qza",
  tree = "/Users/andreabonicelli/Documents/GitHub/microbiome-entomotoxicology-IAN2401/data/rooted-tree-rep-seqs.qza",
  taxonomy="/Users/andreabonicelli/Documents/GitHub/microbiome-entomotoxicology-IAN2401/data/taxonomy-silva-v4.qza",
  metadata = "/Users/andreabonicelli/Documents/GitHub/microbiome-entomotoxicology-IAN2401/data/metadata.tsv"
  ) ) # 141 taxa 

# remove contaminants
( physeq <- physeq %>%
    subset_taxa(
      Family!= "Mitochondria" | is.na(Family) #&
      # Order!="Chloroplast" | is.na(Class)
    ) %>%
    subset_taxa(
      # Family!= "Mitochondria" | is.na(Family) &
      Order!="Chloroplast" | is.na(Class)
    ) ) # 136 taxa

# prune sample to remove all with the sum of 0
ps_prune <-  prune_samples(sample_sums(physeq)!=0, physeq)
ps_prune  # 49 samples of 54

# check for missing data
any(is.na(ps_prune@otu_table))

# Normalize number of reads in each sample using median sequencing depth
total = median(sample_sums(ps_prune))
standf = function(x, t=total) round(t * (x / sum(x)))
ps = transform_sample_counts(ps_prune, standf)

# check for missing data
any(is.na(ps@otu_table))

# replace unique identifier with consecutive ASV name
# comment for lefser abalysis
ps <-microViz::tax_names2rank(ps, colname = "unique")
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# arrange groups in developement order
sample_data(ps)$stage <-
  factor(
    sample_data(ps)$stage,
    levels = c(
      "adult",
      "egg",
      "first instar larvae",
      "second instar larvae",
      "third instar larvae",
      "pupa",
      "teneral"
    )
  )

ps@tax_table %>% write.csv('/Users/andreabonicelli/Documents/GitHub/microbiome-entomotoxicology-IAN2401/tables/supplementary Table 1.csv')

ps@tax_table %>% datatable()

# analysis of weight ------------------------------------------------------
data.frame(sample_data(ps)) %>%
  drop_na(weight) -> metadata_weight

# Wilcoxon rank-sum test
metadata_weight %>%
  group_by(stage) %>%
  rstatix::wilcox_test(data = ., weight ~ treatment)

ggplot(metadata_weight, aes(treatment, weight)) +
  geom_boxplot(aes(group = treatment, fill = treatment), alpha = 0.7) +
  facet_wrap( ~ stage, scales = 'free_y') +
  theme(legend.position = 'top') +
  stat_compare_means(size= 2.5)


ggsave('/Users/andreabonicelli/Documents/GitHub/microbiome-entomotoxicology-IAN2401/figures/Figure 2.pdf', 
       width = 8, height = 7)

# t-test
metadata_weight %>%
  group_by(stage) %>%
  rstatix::t_test(data = ., weight ~ treatment)

ggplot(metadata_weight, aes(treatment, weight)) +
  geom_boxplot(aes(group = treatment, fill = treatment), alpha = 0.7) +
  facet_wrap( ~ stage, scales = 'free_y') +
  theme(legend.position = 'top') +
  stat_compare_means(method = 't.test', size= 2.5) 






