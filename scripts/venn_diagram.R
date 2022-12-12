# Author: Daniel Moser, AG Doehlemann, University of Cologne, Germany
# Purpose: Creation of venn diagram to disply shared/unqiue proteins found in psiBLAST vs BLASTp analysis

# Load packages
library("tidyverse")
library("ggvenn")
library("grid")

# File path
path = "ENTER YOUR FILE PATH HERE"
BLAST_file = "BLASTp_ensemble_propeptides_hits.txt"
psiBLAST_file = "psiBLAST_ensemble_propeptides_hits.txt"

# Load data
BLAST_data = read_delim(paste0(path,BLAST_file))
psiBLAST_data = read_delim(paste0(path,psiBLAST_file))

# Process data
BLAST_data %>% 
  arrange("bit score") %>% 
  distinct(`subject id`,.keep_all=T) %>% 
  group_by(`query id`) %>% 
  distinct(ORG0,.keep_all=T) -> BLAST_unique # select best hit per organism per peptide

shared <- intersect(BLAST_unique$`subject id`, psiBLAST_data$`subject id`)
unique_BLAST <- setdiff(BLAST_unique$`subject id`, psiBLAST_data$`subject id`)
unique_psiBLAST <- setdiff(psiBLAST_data$`subject id`, BLAST_unique$`subject id`)

x <- list(BLAST_unique$`subject id`, psiBLAST_data$`subject id`)
names(x) <- c("BLASTp", "psiBLAST")
ggvenn(x, columns = c("BLASTp", "psiBLAST"), fill_color=c("#CCCCCC", "#333333")) #generation of venn diagram

