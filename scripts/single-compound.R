suppressPackageStartupMessages({
  library("qiime2R") 
  library("phyloseq")
  library("tidyverse")
  library("ggpubr")
  library("vegan")
  library("ranacapa")
  library("microViz")
  library("reshape2")
  library("dplyr")
  library("rstatix")
  library("microbiomeMarker")
})

# -------------------------------------------------------------------------
ps_otu <- as.data.frame(t(ps@otu_table@.Data))
ps_otu$stage <- ps@sam_data$stage
ps_otu$treatment <- ps@sam_data$treatment
ps_otu$stage <- as.factor(ps_otu$stage)
ps_otu$treatment <- as.factor(ps_otu$treatment)

df_melted <- ps_otu %>%
  dplyr::select(stage, treatment, "ASV1":"ASV136") %>%
  gather(key = "variable", value = "value", -stage, -treatment)

df.summary <- df_melted %>%
  group_by(variable, stage, treatment) %>%
  summarise(
    sd = sd(value, na.rm = FALSE),
    mean = mean(value))

ggplot(df.summary,
       aes(
         x = stage,
         y = mean,
         # ymin = mean - (sd/2),
         # ymax = mean + (sd/2),
         color = treatment
       )) +
  geom_line(aes(group = treatment, linetype = treatment), data = df.summary) +
  # geom_errorbar(aes(group = treatment),width = 0.3,size = 0.5) +
  geom_point(size = 1.5)  +
  facet_wrap( ~ variable, scales = "free_y") + theme_bw(14) +
  theme(axis.text.x = element_text(angle = 40, vjust = 1.05, hjust=1),
        axis.text=element_text(size=8), legend.position = 'bottom') +
  xlab('Developmental stage') + ggtitle("All variables") + ylab('Abundance')

ggsave("/Users/andreabonicelli/Documents/GitHub/microbiome-entomotoxicology-IAN2401/figures/supplementary Figure 2.pdf", 
       height = 20, width = 20)

# selected ASV ------------------------------------------------------------
df_melted <- ps_otu %>%
  dplyr::select(stage, treatment, "ASV4","ASV17","ASV27","ASV28","ASV98","ASV99") %>%
  gather(key = "variable", value = "value", -stage, -treatment)

df.summary <- df_melted %>%
  group_by(variable, stage, treatment) %>%
  summarise(
    sd = sd(value, na.rm = FALSE),
    mean = mean(value))

scales::hue_pal()(6)

ggplot(df.summary,
       aes(
         x = stage,
         y = mean,
         # ymin = mean - (sd/2),
         # ymax = mean + (sd/2),
         color = treatment
       )) +
  geom_line(aes(group = treatment, linetype = treatment), data = df.summary) +
  # geom_errorbar(aes(group = treatment), width = 0.3,size = 0.8) +
  geom_point(size = 1.5)  +
  facet_wrap( ~ variable, scales = "free_y") + theme_bw(14) +
  theme(axis.text.x = element_text(angle = 40, vjust = 1.05, hjust=1),
        axis.text=element_text(size=12), legend.position = 'top') +
  xlab('Developmental stage') + ylab('Abundance')

ggsave("/Users/andreabonicelli/Documents/GitHub/microbiome-entomotoxicology-IAN2401/figures/Figure 7.pdf", 
       height = 6, width = 8)

# -------------------------------------------------------------------------
# Pairwise t-test
ps_otu %>%
  group_by(treatment) %>% 
  pairwise_t_test(ASV4 ~ stage) %>%
  adjust_pvalue(method = "fdr") %>% 
  add_significance("p.adj") -> ASV4_pair

ps_otu %>%
  group_by(treatment) %>% 
  pairwise_t_test(ASV17 ~ stage) %>%
  adjust_pvalue(method = "fdr") %>% 
  add_significance("p.adj") -> ASV17_pair

ps_otu %>%
  group_by(treatment) %>% 
  pairwise_t_test(ASV27 ~ stage) %>%
  adjust_pvalue(method = "fdr") %>% 
  add_significance("p.adj") -> ASV27_pair

ps_otu %>%
  group_by(treatment) %>% 
  pairwise_t_test(ASV28 ~ stage) %>%
  adjust_pvalue(method = "fdr") %>% 
  add_significance("p.adj") -> ASV28_pair

ps_otu %>%
  group_by(treatment) %>% 
  pairwise_t_test(ASV98 ~ stage) %>%
  adjust_pvalue(method = "fdr") %>% 
  add_significance("p.adj") -> ASV98_pair

ps_otu %>%
  group_by(treatment) %>% 
  pairwise_t_test(ASV99 ~ stage) %>%
  adjust_pvalue(method = "fdr") %>% 
  add_significance("p.adj") -> ASV99_pair

rbind(ASV4_pair, ASV17_pair, ASV27_pair, ASV28_pair, ASV98_pair, ASV99_pair) %>% 
  write.csv("/Users/andreabonicelli/Documents/GitHub/microbiome-entomotoxicology-IAN2401/tables/supplementary Table 3.csv")
