# Daniel Moser 01.12.2022
# Evaluation of PR gene expression after phytocytokine treatment

# Load Library
library(tidyverse)

# File Path
path <- "YOUR FOLDER PATH/"
file <- "qPCR data/PR_gene_data_24h_R.txt"
filepath <- paste0(path,file)

# Load Data
data <- read_delim(file = filepath, delim="\t", col_select = c("replicate",
                                                           "treatment",
                                                           "PR3",
                                                           "PR4",
                                                           "PR5",
                                                           "PRm6b",
                                                           "PR10.1"))
# Process Data
data <- na.omit(data)

data %>% pivot_longer(cols=starts_with("PR")) -> long
long$treatment <- factor(long$treatment, levels = c("mock", "flg22", "chitin", "Zip1", "IRP", "PSK1", "PIP1", "DAMP1"))

# Statistics
long %>% group_by(name, treatment) %>%
  summarise(avg = mean(value),
            sde = sd(value),#/sqrt(length(value)-1),
            len = length(value),
            ttest = t.test(value, 
                           long$value[long$treatment=="Mock"&long$name==cur_group()$name], 
                           paired=FALSE)$p.value,
            sign = ifelse(ttest<0.001,"***",ifelse(ttest<0.01,"**", ifelse(ttest<0.05,"*","")))
  ) -> stats

# Plotting Figures
for (i in unique(long$name)) {
  long %>% filter(name == i) -> long_i
  stats %>% filter(name == i) -> stats_i
  print(i)
  p <- ggplot() +
    geom_col(data=stats_i, aes(x=Treatment, y=avg), fill="grey90") +
    geom_text(data=stats_i, aes(x=Treatment, y=1.3*(avg+sde), label = sign),
              size=9) +
    geom_errorbar(data=stats_i, aes(x=Treatment, ymin=avg-sde, ymax=avg+sde), width=.4,
                  position=position_dodge(.9)) +
    geom_point(data=long_i, aes(x=Treatment, y=value)) + 
    theme_classic() + 
    labs(x= "Treatment", y="Relative Expression", title=i) +
    theme(text = element_text(size=12),
          plot.title = element_text(hjust = 0.5, face="bold", size=14)) + 
    scale_y_continuous(n.breaks=10)
  print(p)
  file_png = paste0(path,"phytocytokine 24h_",i,".png")
  ggsave(file_png)
} 

# End
