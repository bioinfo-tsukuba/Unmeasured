# Figure2 & 3

# Dotplot Unmeasured(%) (TF family-Cell type), including measure;0 
source("/Users/saeko/Unmeasured/code/function/Dotplot_fam_ctc_unmeasured_v2.R")
calc_opt <- "FALSE"
#tgt_threshold <- 0.75 #上位25%
tgt_threshold <- 0.50 #上位50%
dotplot <- Dotplot_fam_ctc_unmeasured_v2(calc_opt, tgt_threshold)
ggsave(paste0("/Users/saeko/Unmeasured/plot/Dotplot_fam_ct_unmeasured_", tgt_threshold,"_v2.pdf"), dotplot ,width = 15, height = 10)


# ChIP-Atlas existence ratio (bar plot) 
source("/Users/saeko/Unmeasured/code/function/ChIP_Atlas_ctc_TF_bar_plot_v2.R")
calc_opt <- FALSE # updateのときにつかう。
p_list <- ChIP_Atlas_ctc_TF_bar_plot_v2(calc_opt)
ggsave("/Users/saeko/Unmeasured/plot/Unmeasured_number_CTC_bar_v2.pdf", p_list[[1]], width = 9, height = 7)
ggsave("/Users/saeko/Unmeasured/plot/Unmeasured_ratio_CTC_point_v2.pdf", p_list[[2]], width = 9, height = 7)
ggsave("/Users/saeko/Unmeasured/plot/Unmeasured_ratio_CTC_bar_percent.pdf", p_list[[3]], width = 9, height = 7)
ggsave("/Users/saeko/Unmeasured/plot/Unmeasured_ratio_CTC_bar_percent_tate.pdf", p_list[[4]], width = 9, height = 14)


# Number of publications/DEGs (violin plot, scatter plot) -----
source("/Users/saeko/Unmeasured/code/function/DEG_pub_violin.R")
p_list <- DEG_pub_violin()
ggsave("/Users/saeko/Unmeasured/plot/DEG_marker_plot/Number_of_DEG_violin_bwn_measure.pdf", p_list[[1]], width = 6, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/DEG_marker_plot/Number_of_pub_violin_bwn_measure.pdf", p_list[[2]], width = 6, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/DEG_marker_plot/Number_of_DEG_violin_bwn_marker.pdf", p_list[[3]], width = 6, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/DEG_marker_plot/Number_of_pub_violin_bwn_marker.pdf", p_list[[4]], width = 6, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/DEG_marker_plot/Scatte_pub_DEG_bwn_measure.pdf", p_list[[5]], width = 12, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/DEG_marker_plot/Scatte_pub_DEG_bwn_marker.pdf", p_list[[6]], width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/DEG_marker_plot/Scatte_pub_DEG_filter_marker.pdf", p_list[[7]], width = 9, height = 9)
