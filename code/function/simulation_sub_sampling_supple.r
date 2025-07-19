simulation_sub_sampling <- function(tgt_ratio = tgt_ratio){
  
  # base: /Users/saeko/Unmeasured/code/analysis/Revise/Fig5_simulation_subsampling.Rmd
  
  df_join <- read_tsv(paste0("/Users/saeko/Unmeasured/data/GWAS_simulation_v1/2025_07_03_sub_sampling/auc_list_plot_", tgt_ratio,".tsv"))
  p <- df_join %>% 
    ggplot(aes(x = reorder(simulation_ID, auc), y = auc, color = color)) +
    geom_boxplot() +
    scale_color_manual(values = c("black", "blue", "red")) +
    xlab("") +
    ylab(paste0("AUC (", tgt_ratio*100, "% sub-sampling)")) +
    theme(plot.title = element_text(size = unit(15, "pt"), face="bold"), 
          legend.position = "none",
          legend.title = element_text(size = unit(20, "pt")),
          legend.text = element_text(size = unit(20, "pt")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.ticks = element_line(colour = "black", linewidth = unit(0.5, "pt")),
          axis.line = element_line(linewidth = unit(0.5, "pt")),
          axis.title = element_text(size = unit(15, "pt"), colour = "black"),
          axis.text = element_text(size = unit(15, "pt"), colour = "black"),
          #axis.text.x = element_blank(),
          # axis.line = element_line(colour="black"),
          aspect.ratio = 0.4
    ) 
  
  return(p)
  
  
  
}