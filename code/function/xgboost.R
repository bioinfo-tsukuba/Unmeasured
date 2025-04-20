xgboost_plot <- function(){
  
  
  library(tidyverse)
  
  # data prepareation; ignored_pred_v1.R
  #df_reg = read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/ignored_paper/df_reg_902TF_427var.tsv")
  
  # scatter plot
  scatter_df <- read_tsv("/Users/saeko/Unmeasured/data/xgboost/scatter_df.tsv")
  cor <- cor(scatter_df$y_train, scatter_df$preds, method = "spearman")
  p1 <- scatter_df %>% ggplot(aes(x = y_train, y = preds)) +
    geom_point(size = 3) +
    xlab("Number of ChIP-seq (z-scored)") +
    ylab("Prediction (z-scored)") +
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
  res_df <- read_csv("/Users/saeko/Unmeasured/data/xgboost/feature_importances.csv")
  colnames(res_df) <- c("number", "Importance")
  var_num_df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/ignored_paper/exp_var_number_table.tsv")
  res_df2 <- res_df %>% left_join(var_num_df, by = "number")
  res_df3 <- res_df2 %>% mutate(Feature = ifelse(number == "pub_num", "pub_num", var)) %>% select(-var) %>% arrange(-Importance)
  top_feature <- res_df3[1:30,]$Feature
  top_faeture_name <- c("Pub_num", "low_complex", "exon_var.mis", "glutamin_acid", "exon_var.lof", 
                        "aspartic_acid","SignalP", "exon_var.syn", "antisense", "codon_bias",
                        "gene_essentiality_A375", "code_seq_RNA_AAG", "glutamine", "cystein", "code_seq_RNA_ATC",
                        "arginine", "exome.variation", "total.antisense", "arginine", "total.sense",
                        "hydrophobic.amino.acids", "RNA_CCT", "RNA_CAC", "low_complex_20aa", "RNA_exp",
                        "gene_essentiality_HCT116", "Aminoacids_swiss_or_tremb", "histidine", "RNA_GAG", "isoleucine")
  length(top_faeture_name)
  
  # p2 <- res_df3 %>% filter(Feature %in% top_feature) %>% mutate(Feature_name = top_faeture_name) %>%
  #   ggplot(aes(x = reorder(Feature_name, -Importance), y = Importance)) +
  #   geom_bar(stat = "identity", width = 0.7) +
  #   xlab("Feature") +
  #   ylab("Feature Importance") +
  #   theme(plot.title = element_text(size = unit(15, "pt"), face="bold"), 
  #         legend.position = "none",
  #         legend.title = element_text(size = unit(15, "pt")),
  #         legend.text = element_text(size = unit(15, "pt")),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         panel.background = element_blank(), 
  #         axis.ticks = element_line(colour = "black", linewidth = unit(0.5, "pt")),
  #         axis.line = element_line(linewidth = unit(0.5, "pt")),
  #         axis.title = element_text(size = unit(30, "pt"), colour = "black", face = "bold"),
  #         axis.text = element_text(size = unit(30, "pt"), colour = "black", face = "bold"),
  #         axis.text.x =element_text(size = unit(30, "pt"), colour = "black", face = "bold", angle = 45, hjust = 1),
  #         # axis.line = element_line(colour="black"),
  #         aspect.ratio = 1
  #   )
  
  p2 <- res_df3 %>% filter(Feature %in% top_feature) %>% mutate(Feature_name = top_faeture_name) %>%
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
