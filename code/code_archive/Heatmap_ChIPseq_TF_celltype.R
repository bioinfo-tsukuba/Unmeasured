# cell lineごとのヒートマップ (Tissue毎/ENCODE/TF family毎に出す)
# heatmap all ----
library(ComplexHeatmap)
library(grid)
library(circlize)

# TF × Cell type ----
df7 <- readRDS("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/TF_CT_sample_num_new.rds")

# TF family annotation ---
df_fam <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_fam.rds")
fam_list <- df_fam %>% filter(Family != "No_annotation") %>% .$Family %>% unique() %>% unlist() %>% as.character()
Ctc_list <- df_fam$Cell_type_class %>% unique() %>% unlist() %>% as.character()
Ctc_list <- setdiff(Ctc_list, c("Unclassified"))
for (tgt_fam in fam_list[31:42]) {
  print(tgt_fam)
  tgt_tf_list <- df_fam %>% filter(Family == tgt_fam) %>% .$Antigen %>% unique() %>% unlist() %>% as.character()
  for (tgt_ctc in Ctc_list) {
    print(tgt_ctc)
    tgt_ct_list <- df_fam %>% filter(Cell_type_class == tgt_ctc) %>% .$Cell_type %>% unique() %>% unlist() %>% as.character()
    tgt_ct_list <- intersect(colnames(df7), tgt_ct_list)
    result <- df7[tgt_tf_list, tgt_ct_list]
    
    if(length(tgt_ct_list) > 1 & length(tgt_tf_list) > 1){
      # rowごとにsample数集計
      row_sum <- apply(result, 1, sum)
      # colごとにsample数集計
      col_sum <- apply(result, 2, sum)
      
      col_fun = colorRamp2(c(-1, 0, max(result)), c("gray", "white", "green3"))
      
      #各ヒートマップのセルに値を表示するためのカスタム描画関数
      cell_fun <- function(j, i, x, y, width, height, fill) {
        grid::grid.text(sprintf("%d", floor(result[i, j])), x, y, gp = grid::gpar(fontsize = 8))
      }
      
      heat <- Heatmap(result, col = col_fun, cluster_rows = TRUE, cluster_columns = TRUE,  
                      top_annotation = HeatmapAnnotation(sample_sum = anno_barplot(col_sum)), 
                      right_annotation =rowAnnotation(sample = anno_barplot(row_sum)),
                      name = paste0(tgt_fam, " ", tgt_ctc),
                      row_names_gp = gpar(fontsize = 7), 
                      column_names_gp = gpar(fontsize = 7),
                      width = unit(60, "cm"), # 列数に基づく幅
                      height = unit(60, "cm"),
                      show_column_dend = FALSE,
                      show_row_dend = FALSE,
                      rect_gp = gpar(col = "gray"),
                      show_heatmap_legend = TRUE, 
                      cell_fun = cell_fun)
      tgt_fam_name <- gsub("/", "_", tgt_fam)
      png(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/plot/heatmap_chipseq/", tgt_fam_name, "_", tgt_ctc, ".png"),
          width = 2400, height = 2400)
      plot(heat)
      dev.off()
    }#if(length(tgt_ct_list) > 1 & length(tgt_tf_list) > 1)
  }#for (tgt_ctc in Ctc_list)
}#for (tgt_fam in fam_list)
