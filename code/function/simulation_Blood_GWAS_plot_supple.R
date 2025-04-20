simulation_Blood_GWAS_plot_supple <- function(){
  
  library(tidyverse)
  library(ggrepel)
  
  df_all3 <- read_tsv("/Users/saeko/Unmeasured/data/GWAS_simulation_v1/2024_07_13_Blood/df_all3.tsv")
  
  # Simulation number = 0を"real"に変更
  df_real <- df_all3 %>% filter(simulation_ID == 0) %>% select(-simulation_ID) %>% 
    mutate("simulation_ID" = "real") %>% 
    select(colnames(df_all3))
  df_all4 <- rbind(df_all3, df_real) %>% filter(simulation_ID != 0)
  
  tmp <- df_all4 %>% filter(simulation_ID != "real" & simulation_ID != 0)
  tmp_real <- df_all4 %>% filter(simulation_ID == "real") 
  max_auc <- max(tmp$auc)
  min_auc <- min(tmp$auc)
  real_auc <- tmp_real$auc %>% unique()
  max_ID <- tmp %>% filter(auc == max_auc) %>% .$simulation_ID %>% unique()
  min_ID <- tmp %>% filter(auc == min_auc) %>% .$simulation_ID %>% unique()
  
  
  # Year - GWAS cover ratio
  df_all5 <-  df_all4 %>% mutate(color2 = ifelse(simulation_ID == max_ID, "red", color)) %>%
    mutate(color3 = ifelse(simulation_ID == min_ID, "blue", color2))
  p1 <- df_all5 %>% #filter(simulation_ID %in% c(min_ID, max_ID, "real")) %>%
    ggplot(aes(x = year, y = GWAS_cover_ratio)) +
    geom_point(aes(color = color3, group = simulation_ID, alpha = color3)) +
    geom_line(aes(color = color3, fill = color3, group = simulation_ID, alpha = color3), linemitre = 20) +
    scale_color_manual(values = c("black", "blue", "gray", "red")) +
    scale_alpha_manual(values = c("gray" = 0.3, "blue" = 1, "red" = 1, "black" = 1)) +
    ylab("GWAS cover ratio") +
    labs(color = "Simulation number") +
    annotate("text", x = mean(df_all4$year)+2, y = max(df_all4$GWAS_cover_ratio) -0.2, 
             label = paste0("AUC = ", max_auc), 
             vjust = -1, color = "red", size = unit(7, "pt")) +
    annotate("text", x = mean(df_all4$year)+2, y = min(df_all4$GWAS_cover_ratio)+0.05, 
             label = paste0("AUC = ", min_auc), 
             vjust = -1, color = "blue", size = unit(7, "pt")) +
    annotate("text", x = mean(df_all4$year)+2, y = min(df_all4$GWAS_cover_ratio)+0.02, 
             label = paste0("AUC = ", real_auc), 
             vjust = -1, color = "gray", size = unit(7, "pt")) +
    theme(plot.title = element_text(size = unit(15, "pt"), face="bold"), 
          legend.position = c(0.2,0.8),
          legend.title = element_text(size = unit(20, "pt")),
          legend.text = element_text(size = unit(20, "pt")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.ticks = element_line(colour = "black", linewidth = unit(0.5, "pt")),
          axis.line = element_line(linewidth = unit(0.5, "pt")),
          axis.title = element_text(size = unit(30, "pt"), colour = "black"),
          axis.text = element_text(size = unit(30, "pt"), colour = "black"),
          axis.text.x =element_text(size = unit(30, "pt"), colour = "black", angle = 45, hjust = 1),
          # axis.line = element_line(colour="black"),
          aspect.ratio = 1
    ) 
  
  p2 <- df_all5 %>%
    ggplot(aes(x = cum_uniq_TF, y = GWAS_cover_ratio)) +
    geom_point(aes(color = color3, group = simulation_ID, alpha = color3)) +
    geom_line(aes(color = color3, group = simulation_ID, alpha = color3), linemitre = 20) +
    scale_color_manual(values = c("black", "blue", "gray", "red")) +
    scale_alpha_manual(values = c("gray" = 0.3, "blue" = 1, "red" = 1, "black" = 1)) +
    ylab("GWAS cover ratio") +
    labs(color = "Simulation number") +
    annotate("text", x = mean(df_all4$cum_uniq_TF)+2, y = max(df_all4$GWAS_cover_ratio) -0.2, 
             label = paste0("AUC = ", max_auc), 
             vjust = -1, color = "red", size = unit(7, "pt")) +
    annotate("text", x = mean(df_all4$cum_uniq_TF)+2, y = min(df_all4$GWAS_cover_ratio)+0.05, 
             label = paste0("AUC = ", min_auc), 
             vjust = -1, color = "blue", size = unit(7, "pt")) +
    annotate("text", x = mean(df_all4$cum_uniq_TF)+2, y = min(df_all4$GWAS_cover_ratio)+0.02, 
             label = paste0("AUC = ", real_auc), 
             vjust = -1, color = "gray", size = unit(7, "pt")) +
    theme(plot.title = element_text(size = unit(15, "pt"), face="bold"), 
          legend.position = c(0.2,0.8),
          legend.title = element_text(size = unit(20, "pt")),
          legend.text = element_text(size = unit(20, "pt")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.ticks = element_line(colour = "black", linewidth = unit(0.5, "pt")),
          axis.line = element_line(linewidth = unit(0.5, "pt")),
          axis.title = element_text(size = unit(30, "pt"), colour = "black"),
          axis.text = element_text(size = unit(30, "pt"), colour = "black"),
          axis.text.x =element_text(size = unit(30, "pt"), colour = "black", angle = 45, hjust = 1),
          # axis.line = element_line(colour="black"),
          aspect.ratio = 1
    ) 
  
  return(list(p1, p2))
    
}