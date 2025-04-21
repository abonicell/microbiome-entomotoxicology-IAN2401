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
})

# set comparisons ---------------------------------------------------------

my_comparisons <- list(
  c("adult", "egg"),
  c("adult", "first instar larvae"),
  c("adult", "second instar larvae"),
  c("adult", "third instar larvae"),
  c("adult", "pupa"),
  c("adult", "teneral"),
  c("egg", "first instar larvae"),
  c("egg", "second instar larvae"),
  c("egg", "third instar larvae"),
  c("egg", "pupa"),
  c("egg", "teneral"),
  c("first instar larvae", "second instar larvae"),
  c("first instar larvae", "third instar larvae"),
  c("first instar larvae", "pupa"),
  c("first instar larvae", "teneral"),
  c("second instar larvae", "third instar larvae"),
  c("second instar larvae", "pupa"),
  c("second instar larvae", "teneral"),
  c("third instar larvae", "pupa"),
  c("third instar larvae", "teneral"),
  c("pupa","teneral")
)

my_comparisons <- list(
  c("adult", "first instar larvae"),
  c("adult", "third instar larvae"),
  c("adult", "second instar larvae"),
  c("adult", "pupa"),
  c("adult", "teneral"),
  c("egg", "teneral"),
  c("first instar larvae", "teneral"),
  c("second instar larvae", "teneral"),
  c("third instar larvae", "teneral"),
  c("pupa","teneral")
)

# alpha diversity with wilcoxon comparison --------------------------------
A <- plot_richness(ps,
              x = "stage",
              measures = c("Observed", "Shannon"),
              color = "treatment") +
  geom_boxplot(alpha = 0.6) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12,
    )
  ) + 
  stat_compare_means(
    method = "wilcox.test",
    comparisons = my_comparisons,
    label = "p.signif", size= 2.5
  ) 

A

# alpha diversity with wilcoxon comparison --------------------------------
B <- plot_richness(ps,
              x = "treatment",
              measures = c("Observed", "Shannon"),
              color = "treatment") +
  geom_boxplot(alpha = 0.6) + 
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12,
    )
  ) + stat_compare_means(size= 2.5)

B

(B | A) + plot_annotation(tag_levels = 'A')

# alpha diversity with wilcoxon comparison --------------------------------
ps_sex <- subset_samples(ps, sex %in% c('male', 'female'))

my_comparisons <- list(
  c("qst", "2nd"))

C <- plot_richness(ps_sex,
                   x = "sex",
                   measures = c("Observed", "Shannon"),
                   color = "generation") +
  geom_boxplot(alpha = 0.6) + 
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12,
    )
  ) + stat_compare_means(method = "wilcox.test", size= 2.5) 

C

(B + C) / A + plot_annotation(tag_levels = 'A') +
  plot_layout(heights =  c(1, 1.5))

ggsave('/Users/andreabonicelli/Documents/GitHub/microbiome-entomotoxicology-UND/figures/Figure 3.pdf', width = 11, height = 12)

# beta diversity ----------------------------------------------------------
set.seed(123)
bray_dist1 = phyloseq::distance(ps, method="bray")
ordination = ordinate(ps, method="PCoA", distance=bray_dist1)
plot_ordination(ps, ordination, color = "stage", shape = "treatment")   +
  geom_point(aes(color=`stage`), alpha=0.5, size = 4)  + 
  facet_wrap(~stage) + 
  theme(strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9))

set.seed(123)
adonis2(bray_dist1 ~ sample_data(ps)$stage, permutations = 10000)
# Df SumOfSqs      R2      F    Pr(>F)    
# Model     6   7.5527 0.36703 4.0589 9.999e-05 ***
#   Residual 42  13.0253 0.63297                     
# Total    48  20.5779 1.00000                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
adonis2(bray_dist1 ~ sample_data(ps)$treatment, permutations = 10000)
# Df SumOfSqs      R2      F Pr(>F)
# Model     1   0.5162 0.02509 1.2094 0.2587
# Residual 47  20.0617 0.97491              
# Total    48  20.5779 1.00000   

ggsave('/Users/andreabonicelli/Documents/GitHub/microbiome-entomotoxicology-UND/figures/Figure 4.pdf', width = 8.5, height = 6)


