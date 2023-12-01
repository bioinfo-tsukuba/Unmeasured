scatter_expression_cover_ratio_IC50 <- function(tgt_region){
  
  library(tidyverse)
  library(RColorBrewer)
  
  ct_list <- read_tsv("/Users/saeko/Unmeasured/data/cell_type_list_RNA_ChIPAtlas.txt", col_names = F)
  colnames(ct_list) <- c("ct", "ct2")
  
  for (tgt_threshold in 1:10) {
    for (tgt_row in 1:nrow(ct_list)) {
      tgt_ct <- ct_list[tgt_row, "ct"] %>% as.character() #RNA-seq
      tgt_ct2 <- ct_list[tgt_row, "ct2"] %>% as.character() #ChIP-Atlas
      
      # RNAseq
      if(file.exists(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/VirtualChIP/RNAseq/", tgt_ct, "_RNA.tsv"))==TRUE){
        RNA_df <- read_tsv(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/VirtualChIP/RNAseq/", tgt_ct, "_RNA.tsv"))
      }else{
        RNA_df <- read_tsv(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/VirtualChIP/RNAseq/CCLE/", tgt_ct, "_RNA.tsv"))
      }
      colnames(RNA_df) <- c("gene", "RPKM_or_TPM") 
      
      # UniqueTF, Ratio, scatter
      RR_cover_ratio_tib <- suppressMessages(read_tsv(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/RR_ratio_table/RR_cover_ratio_", tgt_ct2, "_",tgt_region, ".tsv")))
      uniq_TF_df <- suppressMessages(read_tsv(paste0("/Users/saeko/Unmeasured/data/uniqTF_table/", tgt_ct2, "_", tgt_region, "bp.tsv")))
      tgt_RR_cover_ratio_tib <- RR_cover_ratio_tib %>% filter(threshold == tgt_threshold)
      df_join1 <- tgt_RR_cover_ratio_tib %>% left_join(uniq_TF_df, by = "year") %>%
        mutate(Cell_type = tgt_ct2, tgt_region = tgt_region)
      
    }
  }
 
  
}