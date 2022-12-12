# Author: Daniel Moser, AG Doehlemann, University of Cologne, Germany
# Purpose: Statistical tests on qPCR data and generation of plots

# Library loading
library(tidyverse)
library(scales)

# File paths
path <- "ENTER YOUR FILEPATH HERE/"
file <- "qPCR data/3h+24h_foldchange_R.txt"
filepath <- paste0(path,file)

# Reading in data
data <- read_delim(file = filepath, delim="\t", col_select = c(1:13))

# Processing data
data %>% pivot_longer(cols = !(starts_with("Treatment")|starts_with("Time"))) -> long #tidy data
long <- na.omit(long) #removal of NAs
long$Treatment <- factor(long$Treatment, 
                         levels = c("mock", "flg22", "chitin", "IRP", "PSK1", "PIP1", "Zip1","DAMP1"))
long$Time <- factor(long$Time, levels = c("3hpi", "24hpi"))
long$name <- factor(long$name)
long %>% filter(Treatment != "DAMP1") ->long #not-used genes/treatments

# Statistics on different Time points, Genes and Treatments
long %>% group_by(Time, name, Treatment) %>%
  summarise(avg = mean(value, na.rm= T), #average
            max = max(value, na.rm = T), #maximum value for plotting reasons
            sd = sd(value, na.rm= T), #standard-deviation
            len = length(value),
            deg = ifelse(avg<1,"down", "up"),
            ttest = ifelse(cur_group()$Treatment=="mock",#|len<2, #t-test
                          1, 
                          t.test(value, 
                                 long$value[long$Treatment=="mock"&
                                              long$name==cur_group()$name&
                                              long$Time==cur_group()$Time], 
                                 paired=F,
                                 var.equal=T)$p.value
                          ),
            sign = ifelse(ttest<0.001,"***",ifelse(ttest<0.01,"**", ifelse(ttest<0.05,"*","")))
            ) -> stats

# Creation of one plot per tested gene, two groups: 3hpi and 24hpi
for (i in unique(long$name)) {
  long %>% filter(name == i) -> long_i
  stats %>% filter(name == i) -> stats_i
  print(i)
  z <- max(long$value[long$name == i]*1.2) #upper border of each plot
  p <- ggplot() +
    geom_col(data=stats_i, #bar plot
             aes(x=Treatment, y=avg, fill=Time), 
             position = position_dodge(width = 0.8), 
             width=0.7) +
    geom_text(data=stats_i, #significance stars
              aes(x=Treatment, y=z, label = sign, group=Time),
              size=5,
              position = position_dodge(width = 0.8)) +
    geom_errorbar(data=stats_i, #error bars
                   aes(x=Treatment, ymin=avg-sd, ymax=avg+sd, group=Time), 
                   width=.4,
                   position=position_dodge(0.8)) +
    geom_point(data=long_i, #data points
               aes(x=Treatment, y=value, group=Time),
               position = position_dodge(width = 0.8),
               alpha=0.7) + 
    theme_classic() + 
    labs(x= "Treatment", y="Fold change", title=i) +
    theme(text = element_text(size=12),
          plot.title = element_text(hjust = 0.5, face="bold", size=14)) + 
    scale_y_continuous(n.breaks=10, expand = c(0,0), limits = c(0,z*1.1), oob=rescale_none) +
    scale_fill_manual(name="Timepoint",
                      values=c("#d0d3d4", "#616a6b")
                      )
  print(p)
  # saving files
  file_png = paste0(path,
                    "plots/phytocytokine_3+24h_foldchange_",
                    i,
                    "_new.png")
  ggsave(file_png)
  file_pdf = paste0(path,
                    "plots/phytocytokine_3+24h_foldchange_",
                    i,
                    "_new.pdf")
  ggsave(file_pdf)
} 

# overview stats
stats %>% select(Time, name, Treatment, len, avg, ttest, sign, deg) -> stats_all
stats %>% filter(sign!="") %>% select(Time, name, Treatment, len, avg, ttest, sign, deg) -> stats_sig


write_delim(
  stats_all,
  paste0(path,"stats_3+24h_all.txt"),
  delim = "\t",
  na = "NA"
)

write_delim(
  stats_sig,
  paste0(path,"stats_3+24h_sig.txt"),
  delim = "\t",
  na = "NA"
)

# same plots as before but separately for every time point
for (i in unique(long$name)) {
  for (time in unique(long$Time)) {
  long %>% filter(name == i & Time == time) -> long_i
  stats %>% filter(name == i & Time == time) -> stats_i
  print(i)
  z <- max(long$value[long$name == i& long$Time == time]*1.2) 
  p <- ggplot() +
    geom_col(data=stats_i, 
             aes(x=Treatment, y=avg, fill=Time), 
             position = position_dodge(width = 0.8), 
             width=0.7) +
    geom_text(data=stats_i,
              aes(x=Treatment, y=z, label = sign, group=Time),
              size=5,
              position = position_dodge(width = 0.8)) +
    geom_errorbar(data=stats_i, 
                  aes(x=Treatment, ymin=avg-sd, ymax=avg+sd, group=Time), 
                  width=.4,
                  position=position_dodge(0.8)) +
    geom_point(data=long_i, 
               aes(x=Treatment, y=value, group=Time),
               position = position_dodge(width = 0.8),
               alpha=0.7) + 
    theme_classic() + 
    labs(x= "Treatment", y="Fold change", title=i) +
    theme(text = element_text(size=12),
          plot.title = element_text(hjust = 0.5, face="bold", size=14)) + 
    scale_y_continuous(n.breaks=10, expand = c(0,0), limits = c(0,z*1.1), oob=rescale_none) +
    scale_fill_manual(name="Timepoint", 
                      values=c("#d0d3d4", "#616a6b"))
    
  print(p)
  file_png = paste0(path,
                    "plots/phytocytokine_foldchange_",
                    i,
                    "_",
                    time,
                    "_new.png")
  ggsave(file_png)
  file_pdf = paste0(path,
                    "plots/phytocytokine_foldchange_",
                    i,
                    "_",
                    time,
                    "_new.pdf")
  ggsave(file_pdf)
}} 



