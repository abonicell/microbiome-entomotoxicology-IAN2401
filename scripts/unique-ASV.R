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
})

# isolate environmental contamination and remove ASV with value 0
ps_adult <- subset_samples(ps, stage == "adult")
ps_adult <-  prune_taxa(taxa_sums(ps_adult) > 0, ps_adult)
ps_adult_list <- rownames(ps_adult@otu_table)

ps_egg <- subset_samples(ps, stage == "egg")
ps_egg <-  prune_taxa(taxa_sums(ps_egg) > 0, ps_egg)
ps_egg_list <- rownames(ps_egg@otu_table)

ps_first <- subset_samples(ps, stage == "first instar larvae")
ps_first <-  prune_taxa(taxa_sums(ps_first) > 0, ps_first)
ps_first_list <- rownames(ps_first@otu_table)

ps_second <- subset_samples(ps, stage == "second instar larvae")
ps_second <-  prune_taxa(taxa_sums(ps_second) > 0, ps_second)
ps_second_list <- rownames(ps_second@otu_table)

ps_third <- subset_samples(ps, stage == "third instar larvae")
ps_third <-  prune_taxa(taxa_sums(ps_third) > 0, ps_third)
ps_third_list <- rownames(ps_third@otu_table)

ps_pupa <- subset_samples(ps, stage == "pupa")
ps_pupa <-  prune_taxa(taxa_sums(ps_pupa) > 0, ps_pupa)
ps_pupa_list <- rownames(ps_pupa@otu_table)

ps_teneral <- subset_samples(ps, stage == "teneral")
ps_teneral <-  prune_taxa(taxa_sums(ps_teneral) > 0, ps_teneral)
ps_teneral_list <- rownames(ps_teneral@otu_table)

# find ASV only in sample adult
adult_list <- ps_adult_list[!(ps_adult_list %in% c(ps_egg_list, 
                                             ps_first_list, 
                                             ps_second_list, 
                                             ps_third_list, 
                                             ps_pupa_list, 
                                             ps_teneral_list))]

r <- rownames(tax_table(ps)) %in% adult_list
DT::datatable(tax_table(ps)[r,]) 

tax_table(ps)[r,] 

# find ASV only in sample egg
egg_list <- ps_egg_list[!(ps_egg_list %in% c(ps_adult_list, 
                                                   ps_first_list, 
                                                   ps_second_list, 
                                                   ps_third_list, 
                                                   ps_pupa_list, 
                                                   ps_teneral_list))]

r <- rownames(tax_table(ps)) %in% egg_list
DT::datatable(tax_table(ps)[r,]) 

tax_table(ps)[r,] 

# find ASV only in sample first instar
first_list <- ps_first_list[!(ps_first_list %in% c(ps_adult_list, 
                                             ps_egg_list, 
                                             ps_second_list, 
                                             ps_third_list, 
                                             ps_pupa_list, 
                                             ps_teneral_list))]

r <- rownames(tax_table(ps)) %in% first_list
DT::datatable(tax_table(ps)[r,]) 

tax_table(ps)[r,]
              
# find ASV only in sample second instar
second_list <- ps_second_list[!(ps_second_list %in% c(ps_adult_list, 
                                                   ps_egg_list, 
                                                   ps_first_list, 
                                                   ps_third_list, 
                                                   ps_pupa_list, 
                                                   ps_teneral_list))]

r <- rownames(tax_table(ps)) %in% second_list
DT::datatable(tax_table(ps)[r,]) 

tax_table(ps)[r,] 

# find ASV only in sample third instar
third_list <- ps_third_list[!(ps_third_list %in% c(ps_adult_list, 
                                                      ps_egg_list, 
                                                      ps_first_list, 
                                                      ps_second_list, 
                                                      ps_pupa_list, 
                                                      ps_teneral_list))]

r <- rownames(tax_table(ps)) %in% third_list
DT::datatable(tax_table(ps)[r,]) 

tax_table(ps)[r,] 

# find ASV only in sample pupa
pupa_list <- ps_pupa_list[!(ps_pupa_list %in% c(ps_adult_list, 
                                                   ps_egg_list, 
                                                   ps_first_list, 
                                                   ps_second_list, 
                                                  ps_third_list, 
                                                   ps_teneral_list))]

r <- rownames(tax_table(ps)) %in% pupa_list
DT::datatable(tax_table(ps)[r,]) 

tax_table(ps)[r,] 

# find ASV only in sample teneral
teneral_list <- ps_teneral_list[!(ps_teneral_list %in% c(ps_adult_list, 
                                                ps_egg_list, 
                                                ps_first_list, 
                                                ps_second_list, 
                                                ps_third_list, 
                                                ps_pupa_list))]

r <- rownames(tax_table(ps)) %in% teneral_list
tax_table(ps)[r,] 

# unique ASV per treatment ------------------------------------------------
ps_control <- subset_samples(ps, treatment == "control")
ps_control <-  prune_taxa(taxa_sums(ps_control) > 0, ps_control)
ps_control_list <- rownames(ps_control@otu_table)

ps_treated <- subset_samples(ps, treatment == "flunitrazepam")
ps_treated <-  prune_taxa(taxa_sums(ps_treated) > 0, ps_treated)
ps_treated_list <- rownames(ps_treated@otu_table)

# find ASV only in sample control
control_list <- ps_control_list[!(ps_control_list %in% ps_treated)]

r <- rownames(tax_table(ps)) %in% control_list
tax_table(ps)[r,] 

# find ASV only in sample adult
treated_list <- ps_treated_list[!(ps_treated_list %in% ps_control)]

r <- rownames(tax_table(ps)) %in% ps_treated_list
tax_table(ps)[r,]