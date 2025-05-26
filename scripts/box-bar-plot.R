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
  library("patchwork")
  library("ggsci")
})

# bar plot prepare data ---------------------------------------------------
prepare_bar_data <- function(physeq_object,
                             grouping_var = "grouping_var",
                             tax_level = NULL) {
  # Optional: Agglomerate to the specified taxonomic level
  if (!is.null(tax_level)) {
    physeq_object <- tax_glom(physeq_object, tax_level)
  }
  # Merge samples by the grouping variable
  merged <- merge_samples(physeq_object, grouping_var)
  # Get top ASV names sorted by abundance
  ASVnames <- names(sort(taxa_sums(merged), decreasing = TRUE))
  # Prune taxa from both original and merged phyloseq objects
  grouped <- prune_taxa(ASVnames, merged)
  # Add grouping info to sample metadata
  grouped@sam_data[[grouping_var]] <- rownames(grouped@sam_data)
  # Transform counts to relative abundance (%)
  out <- transform_sample_counts(grouped, function(x)
    (x / sum(x)) * 100)
  
  return(out)
}

# -------------------------------------------------------------------------
# group taxa by similarity
phy <- phyloseq::tax_glom(ps, "Phylum")
# assign name
(phyloseq::taxa_names(phy) <- phyloseq::tax_table(phy)[, "Phylum"])

# abundance
phyloseq::psmelt(phy) %>%
  group_by(OTU) %>% tally(Abundance) %>% arrange(desc(n)) %>%
  mutate(rel.freq = paste0(round(100 * n / sum(n), 5), "%"))

# boxplot of taxa (phylum level) for stage
phyloseq::psmelt(phy) %>%
  ggplot(data = ., aes(x = stage, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "Development", y = "Abundance\n") +
  facet_wrap(~ OTU + treatment, scales = "free_y") +
  labs(title = NULL, subtitle = NULL) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95)) +
  theme(legend.position = "none")


# full sample barplot -----------------------------------------------------
phy <- phyloseq::tax_glom(ps, "Phylum")
phy_tr <- transform_sample_counts(ps, function(x)
  x / sum(x) * (100))

plot_bar(phy_tr, x = "Sample", fill = "Phylum") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Abundance (%)", title = "Control") +
  theme(axis.text.x = element_text(
    size = rel(0.9),
    angle = 45,
    hjust = 1,
    vjust = 1
  )) + facet_wrap( ~ stage, scales = "free_x", nrow = 1)

# barplot -----------------------------------------------------------------
phy_control <- subset_samples(ps, treatment == "control")
ps_phylum_ctrl <- prepare_bar_data(phy_control, grouping_var = "stage", tax_level =
                                     "Phylum")

# arrange groups in developement order
sample_data(ps_phylum_ctrl)$stage <-
  factor(
    sample_data(ps_phylum_ctrl)$stage,
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

a <- plot_bar(ps_phylum_ctrl, x = "stage", fill = "Phylum") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "PMI (weeks)", y = "Abundance (%)", title = "Control") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      size = rel(0.9),
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )

# barplot -----------------------------------------------------------------
phy_control <- subset_samples(ps, treatment == "flunitrazepam")
ps_phylum_flun <- prepare_bar_data(phy_control, tax_level = "Phylum")

# arrange groups in developement order
sample_data(ps_phylum_flun)$stage <-
  factor(
    sample_data(ps_phylum_flun)$stage,
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

b <- plot_bar(ps_phylum_flun, x = "stage", fill = "Phylum") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "PMI (weeks)", y = "Abundance (%)", title = "Flunitrazepam") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      size = rel(0.9),
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )

ggarrange(
  a,
  b,
  labels = "AUTO",
  nrow = 1,
  ncol = 2,
  common.legend = TRUE
)

ggsave(
  '/Users/andreabonicelli/Documents/GitHub/microbiome-entomotoxicology-IAN2401/figures/Figure 5.pdf',
  width = 8,
  height = 4.5
)

# barplot -----------------------------------------------------------------
cls_control <- subset_samples(ps, treatment == "control")
ps_class_ctrl <- prepare_bar_data(cls_control, grouping_var = "stage", tax_level =
                                    "Class")

# arrange groups in developement order
sample_data(ps_class_ctrl)$stage <-
  factor(
    sample_data(ps_class_ctrl)$stage,
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

c <- plot_bar(ps_class_ctrl, x = "stage", fill = "Class") +
  geom_bar(aes(fill = Class), stat = "identity", position = "stack") +
  labs(x = "PMI (weeks)", y = "Abundance (%)", title = "Control") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      size = rel(0.9),
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )

# barplot -----------------------------------------------------------------
cls_control <- subset_samples(ps, treatment == "flunitrazepam")
ps_class_flun <- prepare_bar_data(cls_control, grouping_var = "stage", tax_level =
                                    "Class")

# arrange groups in developement order
sample_data(ps_class_flun)$stage <-
  factor(
    sample_data(ps_class_flun)$stage,
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

d <- plot_bar(ps_class_flun, x = "stage", fill = "Class") +
  geom_bar(aes(fill = Class), stat = "identity", position = "stack") +
  labs(x = "PMI (weeks)", y = "Abundance (%)", title = "Flunitrazepam") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      size = rel(0.9),
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )

ggarrange(
  c,
  d,
  labels = "AUTO",
  nrow = 1,
  ncol = 2,
  common.legend = TRUE
)

ggsave(
  '/Users/andreabonicelli/Documents/GitHub/microbiome-entomotoxicology-IAN2401/figures/supplementary Figure 1.pdf',
  width = 8,
  height = 4.5
)
