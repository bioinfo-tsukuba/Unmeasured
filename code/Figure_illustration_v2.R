# ChIP-Atlas existence ratio (bar plot) -----
## Cell type class
### version1
source("/Users/saeko/Unmeasured/code/function/ChIP_Atlas_ctc_TF_bar_plot.R")
calc_opt <- FALSE # updateのときにつかう。
p_list <- ChIP_Atlas_ctc_TF_bar_plot(calc_opt)
ggsave("/Users/saeko/Unmeasured/plot/Unmeasured_number_CTC_bar.pdf", p_list[[1]], width = 9, height = 7)
ggsave("/Users/saeko/Unmeasured/plot/Unmeasured_ratio_CTC_point.pdf", p_list[[2]], width = 9, height = 7)

### version2 (CTCごとにFPKMのthresholdを変えた)
source("/Users/saeko/Unmeasured/code/function/ChIP_Atlas_ctc_TF_bar_plot_v2.R")
calc_opt <- FALSE # updateのときにつかう。
p_list <- ChIP_Atlas_ctc_TF_bar_plot_v2(calc_opt)
ggsave("/Users/saeko/Unmeasured/plot/Unmeasured_number_CTC_bar_v2.pdf", p_list[[1]], width = 9, height = 7)
ggsave("/Users/saeko/Unmeasured/plot/Unmeasured_ratio_CTC_point_v2.pdf", p_list[[2]], width = 9, height = 7)
ggsave("/Users/saeko/Unmeasured/plot/Unmeasured_ratio_CTC_bar_percent.pdf", p_list[[3]], width = 9, height = 7)

## Cell type
# source("/Users/saeko/Unmeasured/code/function/ChIP_Atlas_ct_TF_bar_plot.R")
# calc_opt <- FALSE # updateのときにつかう。
# p_list <- ChIP_Atlas_ct_TF_bar_plot(calc_opt)
# ggsave("/Users/saeko/Unmeasured/plot/Unmeasured_number_CT_bar.pdf", p_list[[1]], width = 9, height = 7)
# ggsave("/Users/saeko/Unmeasured/plot/Unmeasured_ratio_CT_point.pdf", p_list[[2]], width = 9, height = 7)


# Number of publications/DEGs (violin plot, scatter plot) -----
source("/Users/saeko/Unmeasured/code/function/DEG_pub_violin.R")
p_list <- DEG_pub_violin()
ggsave("/Users/saeko/Unmeasured/plot/DEG_marker_plot/Number_of_DEG_violin_bwn_measure.pdf", p_list[[1]], width = 4.5, height = 7)
ggsave("/Users/saeko/Unmeasured/plot/DEG_marker_plot/Number_of_pub_violin_bwn_measure.pdf", p_list[[2]], width = 4.5, height = 7)
ggsave("/Users/saeko/Unmeasured/plot/DEG_marker_plot/Number_of_DEG_violin_bwn_marker.pdf", p_list[[3]], width = 4.5, height = 7)
ggsave("/Users/saeko/Unmeasured/plot/DEG_marker_plot/Number_of_pub_violin_bwn_marker.pdf", p_list[[4]], width = 4.5, height = 7)
ggsave("/Users/saeko/Unmeasured/plot/DEG_marker_plot/Scatte_pub_DEG_bwn_measure.pdf", p_list[[5]], width = 9, height = 7)
ggsave("/Users/saeko/Unmeasured/plot/DEG_marker_plot/Scatte_pub_DEG_bwn_marker.pdf", p_list[[6]], width = 9, height = 7)
ggsave("/Users/saeko/Unmeasured/plot/DEG_marker_plot/Scatte_pub_DEG_filter_marker.pdf", p_list[[7]], width = 9, height = 7)


# RR cover ratio (%) plot (*taking several hours)----
source("/Users/saeko/Unmeasured/code/function/RR_ratio_year.R")
tgt_region <- 100
ct_list <- read_tsv("/Users/saeko/Unmeasured/data/cell_type_list_RNA_ChIPAtlas.txt", col_names = F)
colnames(ct_list) <- c("ct", "ct2")
for (tgt_row in 2:nrow(ct_list)) {
  tgt_ct <- ct_list[tgt_row, "ct"] %>% as.character() #RNA-seq
  tgt_ct2 <- ct_list[tgt_row, "ct2"] %>% as.character() #ChIP-Atlas
  p_list <- RR_cover_ratio(tgt_ct, tgt_ct2, tgt_region)
  ggsave(paste0("/Users/saeko/Unmeasured/plot/RR_ratio_year/", tgt_ct2, "_", tgt_region, "bp.pdf" ), p_list[[1]], width = 9, height = 7)
}

# all TF patchwork
library(patchwork)
source("/Users/saeko/Unmeasured/code/function/RR_ratio_year_allTF.R")
tgt_region <- 500
ct_list <- read_tsv("/Users/saeko/Unmeasured/data/cell_type_list_RNA_ChIPAtlas.txt", col_names = F)
colnames(ct_list) <- c("ct", "ct2")
patch1 <- c()
for (tgt_row in 1:nrow(ct_list)) {
  tgt_ct <- ct_list[tgt_row, "ct"] %>% as.character() #RNA-seq
  tgt_ct2 <- ct_list[tgt_row, "ct2"] %>% as.character() #ChIP-Atlas
  p_list <- RR_cover_ratio_allTF(tgt_ct, tgt_ct2, tgt_region)
  if(tgt_row == 1){
    patch1 <- p_list[[1]] 
  }else{
    patch1 <- patch1 + p_list[[1]]
  }
}

ggsave(paste0("/Users/saeko/Unmeasured/plot/RR_ratio_year/allTF_", tgt_region, "bp.pdf" ), patch1)



# Unique TF plot ----
source("/Users/saeko/Unmeasured/code/function/Uniq_TF_year.R")
ct_list <- read_tsv("/Users/saeko/Unmeasured/data/cell_type_list_RNA_ChIPAtlas.txt", col_names = F)
colnames(ct_list) <- c("ct", "ct2")
for (tgt_row in 1:nrow(ct_list)) {
  tgt_ct <- ct_list[tgt_row, "ct"] %>% as.character() #RNA-seq
  tgt_ct2 <- ct_list[tgt_row, "ct2"] %>% as.character() #ChIP-Atlas
  tgt_region <- 500
  p_list <- Uniq_TF_year(tgt_ct, tgt_ct2, tgt_region)
  ggsave(paste0("/Users/saeko/Unmeasured/plot/Uniq_TF_year/", tgt_ct2, "_", tgt_region, "bp.pdf" ), p_list[[1]], width = 9, height = 7)
  ggsave(paste0("/Users/saeko/Unmeasured/plot/SRX_num_year/", tgt_ct2, "_", tgt_region, "bp.pdf" ), p_list[[2]], width = 9, height = 7)
}


# RR ratioとUniq TFとSRX数をcell lineごとに保存 ------
source("/Users/saeko/Unmeasured/code/function/combine_year.R")
ct_list <- read_tsv("/Users/saeko/Unmeasured/data/cell_type_list_RNA_ChIPAtlas.txt", col_names = F)
colnames(ct_list) <- c("ct", "ct2")
for (tgt_row in 1:nrow(ct_list)) {
  tgt_ct <- ct_list[tgt_row, "ct"] %>% as.character() #RNA-seq
  tgt_ct2 <- ct_list[tgt_row, "ct2"] %>% as.character() #ChIP-Atlas
  tgt_region <- 500
  combine_year(tgt_ct, tgt_ct2, tgt_region)
}

# 1つのthresholdで
source("/Users/saeko/Unmeasured/code/function/combine_year_per_threshold.R")
ct_list <- read_tsv("/Users/saeko/Unmeasured/data/cell_type_list_RNA_ChIPAtlas.txt", col_names = F)
colnames(ct_list) <- c("ct", "ct2")
for (tgt_threshold in 1:10) {
  for (tgt_row in 1:nrow(ct_list)) {
    tgt_ct <- ct_list[tgt_row, "ct"] %>% as.character() #RNA-seq
    tgt_ct2 <- ct_list[tgt_row, "ct2"] %>% as.character() #ChIP-Atlas
    tgt_region <- 500
    combine_year(tgt_ct, tgt_ct2, tgt_region, tgt_threshold)
  }
}


# UniqueTF数, RR ratio, scatter ----
source("/Users/saeko/Unmeasured/code/function/UniqueTF_RRratio_scatter.R")
tgt_region <- 500
p <- UniqueTF_RRratio_scatter(tgt_region)
ggsave("/Users/saeko/Unmeasured/plot/UniqueTF_RRratio_scatter/UniqTF_RRratio_scatter.pdf", p)


# gene expression と RR ratio IC50, scatter ------
source("/Users/saeko/Unmeasured/code/function/.R")
tgt_region <- 500
p <- UniqueTF_RRratio_scatter(tgt_region)

# SRX数, RR ratio, scatter ----
source("/Users/saeko/Unmeasured/code/function/SRX_RRratio_scatter.R")
tgt_region <- 500
p <- SRX_RRratio_scatter(tgt_region)
ggsave("/Users/saeko/Unmeasured/plot/SRX_RRratio_scatter/SRX_RRratio_scatter.pdf", p)

