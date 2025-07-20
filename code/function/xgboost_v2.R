xgboost_plot_v2 <- function(){
  
  library(tidyverse)
  
  # scatter plot
  scatter_df <- read_tsv("/Users/saeko/Unmeasured/data/xgboost/rev1/scatter_data.tsv")
  cor <- cor(scatter_df$Observed_ChIP_num, scatter_df$Predicted_ChIP_num_median, method = "spearman")
  p1 <- scatter_df %>% ggplot(aes(x = Observed_ChIP_num, y = Predicted_ChIP_num_median)) +
    geom_point(size = 3) +
    #geom_bin2d(bins = 40) +
    #scale_fill_continuous(type = "viridis") +
    xlab("Observed number of ChIP-seq (log10)") +
    ylab("Predicted number of ChIP-seq (log10)") +
    xlim(c(0, 3))+
    ylim(c(0, 3)) +
    annotate("text", x = Inf, y = -Inf, label = paste0("Cor = ", round(cor, 3)), 
             hjust = 1.1, vjust = -2, size = 8, fontface = "bold") +
    theme(plot.title = element_text(size = unit(15, "pt"), face="bold"), 
          legend.position = "none",
          legend.title = element_text(size = unit(15, "pt")),
          legend.text = element_text(size = unit(15, "pt")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.ticks = element_line(colour = "black", linewidth = unit(0.5, "pt")),
          axis.line = element_line(linewidth = unit(0.5, "pt")),
          axis.title = element_text(size = unit(30, "pt"), colour = "black", face = "bold"),
          axis.text = element_text(size = unit(30, "pt"), colour = "black", face = "bold"),
          axis.text.x =element_text(size = unit(30, "pt"), colour = "black", face = "bold"),
          # axis.line = element_line(colour="black"),
          aspect.ratio = 1
    )
  
  # result (feature importance)
  res_df <- read_tsv("/Users/saeko/Unmeasured/data/xgboost/rev1/feature_importance_summary.tsv")
  colnames(res_df) <- c("number", "Importance", "Std_Importance")
  var_num_df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/ignored_paper/exp_var_number_table.tsv")
  res_df2 <- res_df %>% left_join(var_num_df, by = "number")
  res_df3 <- res_df2 %>% mutate(Feature = ifelse(number == "pub_num", "pub_num", var)) %>% select(-var) %>% arrange(-Importance)
  top_feature <- res_df3[1:30,]$Feature
  top_feature_name <- c("helix_affine_aa", "hydrophobic_aa", "pub_num", "codon_bias_wright", "RNA_expr_reh", "codon_freq_GGT", 
                        "start_codon_ATG", "RNA_expr_hel", "codon_bias_cds", "exome_var_prec", "RNA_expr_urinarybladder", 
                        "RNA_expr_skeletalmuscle", "codon_freq_CAT", "low_complexity_frac", "codon_freq_AGA", "codon_freq_ACT", 
                        "arginine_content", "RNA_expr_k562", "RNA_expr_endometrium", "exome_var_pnull", "codon_freq_ATT", 
                        "codon_freq_AAG", "RNA_expr_u138mg", "codon_freq_CCT", "RNA_expr_smoothmuscle", "RNA_expr_tonsil", 
                        "codon_bias_novembre", "uracil_content", "protein_length", "RNA_expr_karpas707"
                        )
  length(top_feature_name)
  
  p2 <- res_df3 %>% filter(Feature %in% top_feature) %>% mutate(Feature_name = top_feature_name) %>%
    ggplot(aes(x = reorder(Feature_name, Importance), y = Importance)) +
    geom_segment( aes(xend=Feature_name, yend=0)) +
    geom_point( size=4, color="orange") +
    coord_flip() +
    theme_bw() +
    xlab("Feature") +
    ylab("Feature Importance") +
    theme(plot.title = element_text(size = unit(15, "pt"), face="bold"),
          legend.position = "none",
          legend.title = element_text(size = unit(15, "pt")),
          legend.text = element_text(size = unit(15, "pt")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_line(colour = "black", linewidth = unit(0.5, "pt")),
          axis.line = element_line(linewidth = unit(0.5, "pt")),
          axis.title = element_text(size = unit(20, "pt"), colour = "black"),
          axis.text = element_text(size = unit(20, "pt"), colour = "black"),
          axis.text.x =element_text(size = unit(20, "pt"), colour = "black"),
          # axis.line = element_line(colour="black"),
          aspect.ratio = 2.0
    )
  
  return(list(p1, p2))
  
}