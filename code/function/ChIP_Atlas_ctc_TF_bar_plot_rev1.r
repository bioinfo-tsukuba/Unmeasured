ChIP_Atlas_ctc_TF_bar_plot_rev1 <- function(calc_opt){
  
  # revise1, expressed TF across TF (stepminer)
  library(tidyverse)
  
  mt_unmeasured_tib <- readRDS("/Users/saeko/Unmeasured/data/tib_unmeasured_stepminer_Ctc.rds")
  
  tib_barplot <- mt_unmeasured_tib %>%
    select(Antigen, Cell_type_class, expressed, num_unmeasured) %>%
    group_by(Cell_type_class) %>%
    summarise(
      n = n(),
      total_expressed_TF_Ctc = sum(expressed),
      #total_measured_ChIPseq = sum(num_ChIPseq),
      total_unmeasured_TF_Ctc = sum(num_unmeasured)
    ) %>%
    mutate(unmeasured_ratio = total_unmeasured_TF_Ctc / total_expressed_TF_Ctc,
           Unmeasured_percentage = (unmeasured_ratio*100),
           Measured_percentage = 100 - Unmeasured_percentage) %>%
    pivot_longer(-c(Cell_type_class, n, total_expressed_TF_Ctc, total_unmeasured_TF_Ctc, unmeasured_ratio), 
                 names_to = "measure", values_to = "percentage")
  
  p1 <- tib_barplot %>% 
    ggplot(aes(x = reorder(Cell_type_class, -unmeasured_ratio), y = percentage, fill = measure, label = round(percentage))) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(y = 110)) +
    xlab("Cell type class (tissue)")+
    ylab("Unmeasured(%)")+
    scale_fill_manual(values = c("gray", "blue3")) +
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          legend.position = c(0.5, 1.1),
          legend.direction = "horizontal",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=20, color = "black"),
          axis.text.x =element_text(size=20, color = "black", angle = 90, hjust = 1),
          axis.text.y =element_text(size=20,color = "black"),
          axis.title=element_text(size=20, color = "black"),
          aspect.ratio = 0.5
    )
  
  return(list(p1))
}