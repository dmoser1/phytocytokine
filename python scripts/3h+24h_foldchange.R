library(tidyverse)
path <- "C:/Users/dan2m/OneDrive/Dokumente/Sciebo Neu/" #lab
#path <- "C:/Users/dan2m/sciebo/" #home
file <- "qPCR/2022-10 Phytocytokine/3h+24h_foldchange_R.txt"
filepath <- paste0(path,file)

data <- read_delim(file = filepath, delim="\t", col_select = c(1:13))


data %>% pivot_longer(cols = !(starts_with("Treatment")|starts_with("Time"))) -> long
long <- na.omit(long)
long$Treatment <- factor(long$Treatment, 
                         levels = c("mock", "flg22", "chitin", "IRP", "PSK1", "PIP1", "ZIP1","DAMP1"))
long$Time <- factor(long$Time, levels = c("3hpi", "24hpi"))
long$name <- factor(long$name)
long %>% filter(!name %in% c("ATFP4", "GST42", "lox8", "MPI1") & Treatment != "DAMP1") ->long
long %>% group_by(Time, name, Treatment) %>%
  summarise(avg = mean(value, na.rm= T),
            max = max(value, na.rm = T),
            sd = sd(value, na.rm= T),
            len = length(value),
            deg = ifelse(avg<1,"down", "up"),
            ttest = ifelse(cur_group()$Treatment=="mock",#|len<2, 
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

for (i in unique(long$name)) {
  long %>% filter(name == i) -> long_i
  stats %>% filter(name == i) -> stats_i
  print(i)
  z <- max(long$value[long$name == i]*1.2) 
  p <- ggplot() +
    geom_col(data=stats_i, 
             aes(x=Treatment, y=avg, fill=Time), 
             position = position_dodge(width = 0.8), 
             width=0.7) +
    geom_text(data=stats_i, 
              #aes(x=Treatment, y=z, label = sign, group=Time),
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
    scale_y_continuous(n.breaks=10) +
    scale_fill_manual(name="Timepoint", 
                      #labels=c("3hpi", "24hpi"), 
                      values=c("#d0d3d4", "#616a6b")
                      )
  print(p)
  file_png = paste0(path,
                    "qPCR/2022-10 Phytocytokine/plots/phytocytokine_3+24h_foldchange_",
                    i,
                    "_new.png")
  ggsave(file_png)
  file_pdf = paste0(path,
                    "qPCR/2022-10 Phytocytokine/plots/phytocytokine_3+24h_foldchange_",
                    i,
                    "_new.pdf")
  ggsave(file_pdf)
} 

stats %>% select(Time, name, Treatment, len, avg, ttest, sign, deg) -> stats_all
stats %>% filter(sign!="") %>% select(Time, name, Treatment, len, avg, ttest, sign, deg) -> stats_sig


write_delim(
  stats_all,
  paste0(path,"qPCR/2022-10 Phytocytokine/stats_3+24h_all.txt"),
  delim = "\t",
  na = "NA"
)

write_delim(
  stats_sig,
  paste0(path,"qPCR/2022-10 Phytocytokine/stats_3+24h_sig.txt"),
  delim = "\t",
  na = "NA"
)

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
              #aes(x=Treatment, y=z, label = sign, group=Time),
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
    scale_y_continuous(n.breaks=10) +
    scale_fill_manual(name="Timepoint", 
                      #labels=c("3hpi", "24hpi"), 
                      values=c("#d0d3d4", "#616a6b")
    )
  print(p)
  file_png = paste0(path,
                    "qPCR/2022-10 Phytocytokine/plots/phytocytokine_foldchange_",
                    i,
                    "_",
                    time,
                    "_new.png")
  ggsave(file_png)
  file_pdf = paste0(path,
                    "qPCR/2022-10 Phytocytokine/plots/phytocytokine_foldchange_",
                    i,
                    "_",
                    time,
                    "_new.pdf")
  ggsave(file_pdf)
}} 



