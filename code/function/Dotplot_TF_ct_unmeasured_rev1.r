Dotplot_TF_ct_unmeasured_rev1 <- function(calc_opt, tgt_threshold){
  
  # /Users/saeko/Unmeasured/code/analysis/Revise/Unmeasured_cell_line.Rmd
  library(tidyverse)
  mt_stepMiner_unmeasured <- readRDS("/Users/saeko/Unmeasured/data/mt_stepMiner_unmeasured_TF_ct_v1.rds")
  mt_Chip_ct <- readRDS("/Users/saeko/Unmeasured/data/mt_Chip_ct.rds")
  
  rownames_TF_chip <- rownames(mt_Chip_ct)
  mt_Chip_ct_tib <- mt_Chip_ct %>% as_tibble() %>% mutate(Antigen = rownames_TF_chip) %>%
    pivot_longer(-c(Antigen), names_to = "Cell_type", values_to = "num_ChIPseq") 
  
  rownames_stepMiner_TF_unmeasured <- rownames(mt_stepMiner_unmeasured)
  mt_stepMiner_unmeasured_tib <- mt_stepMiner_unmeasured %>% as_tibble() %>% mutate(Antigen = rownames_stepMiner_TF_unmeasured) %>%
    pivot_longer(-c(Antigen), names_to = "Cell_type", values_to = "unmeasured_TF") 
  
  tib_plot <-  mt_stepMiner_unmeasured_tib %>% left_join(mt_Chip_ct_tib, by = c("Antigen", "Cell_type")) 
  p <- tib_plot %>%
    ggplot(aes(x = reorder(Cell_type, -unmeasured_TF), y = reorder(Antigen, -unmeasured_TF), fill = factor(unmeasured_TF))) +
    geom_tile() +
    scale_fill_manual(values = c("white", "blue")) +
    labs(
      title = "Unmeasured TFs across cell types",
      x = "Cell Type",
      y = "Antigen (TF)",
      fill = "# Unmeasured TFs"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    )
  
  return(p)
}