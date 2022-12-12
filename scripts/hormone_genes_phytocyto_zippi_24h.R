library(tidyverse)
#path <- "C:/Users/dan2m/OneDrive/Dokumente/Sciebo Neu/" #lab
path <- "C:/Users/dan2m/sciebo/" #home
file <- "qPCR/2022-10 Phytocytokine/qPCR24hrep2to4 zippi marker genes.txt"
filepath <- paste0(path,file)

data <- read_delim(file = filepath, delim="\t", col_select = c(1:5))
data <- na.omit(data)

data %>% pivot_longer(cols = !(starts_with("Treatment")|starts_with("replicate"))) -> long

long$Treatment <- factor(long$Treatment, levels = c("Mock", "Flg22", "Chitin", "IRP", "PSK1", "PIP1"))
#long %>% filter(!(Treatment == "Flg22" & Replicate == "4")) -> long

long %>% group_by(name, Treatment) %>%
  summarise(avg = mean(value),
            sde = sd(value),#/sqrt(length(value)-1),
            len = length(value),
            ttest = t.test(value, 
                           long$value[long$Treatment=="Mock"&long$name==cur_group()$name], 
                           paired=F)$p.value,
            sign = ifelse(ttest<0.001,"***",ifelse(ttest<0.01,"**", ifelse(ttest<0.05,"*","")))
            ) -> stats

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
    labs(x= "Treatment", y="Relative Expression", title=paste0(i," ", "24h")) +
    theme(text = element_text(size=12),
          plot.title = element_text(hjust = 0.5, face="bold", size=14)) + 
    scale_y_continuous(n.breaks=10)
  print(p)
  file_png = paste0(path,"qPCR/2022-10 Phytocytokine/plots/phytocytokine_24h_",i,".png")
  ggsave(file_png)
} 


