ChIP_distribution <- function(){
  
  library(tidyverse)
  library(ggrepel)
  library(RColorBrewer)
  library(ggh4x)
  
  # 20231004 version
  df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/SRX_date_20231004.tsv")
  df <- df %>% arrange(year, month, day, time) %>% mutate(order = 1:nrow(df)) %>% 
    filter(Cell_type_class != "Unclassified" & Cell_type_class != "No description") %>%
    filter(Antigen != "GFP" & Antigen != "Epitope tags")
  
  # join TF family
  df_fam <- readRDS("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS_paper_public/data/Fig2/obj/df_fam.rds")
  tgt_IDs <- df$SRX %>% unique() %>% as.character()
  tgt_df_fam <- df_fam %>% filter(ID %in% tgt_IDs) %>% select(Antigen, Family) %>% distinct()
  
  df2 <- df %>% left_join(tgt_df_fam, by = "Antigen") %>% arrange(Family)
  
  #color_list <- c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG"))
  color_list <- c(
    brewer.pal(11, "Spectral"),   # 11色
    brewer.pal(11, "BrBG"),       # 11色
    brewer.pal(11, "Set3"),       # 11色（パステルカラー）
    brewer.pal(11, "Paired")      # 11色（コントラストが強い）
  )
  
  # per TF ------
  p1 <- df2 %>% group_by(Antigen) %>% summarise(n_sample = n()) %>%
    ggplot(aes(x = reorder(Antigen, -n_sample), y = n_sample, label = Antigen)) +
    geom_bar(stat = "identity") +
    ggrepel::geom_label_repel() +
    labs(
      x = "TF",
      y = "Number of TF ChIP-seq"
    )+
    theme(plot.title = element_text(size = unit(12, "pt")), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.ticks = element_line(colour = "black", linewidth = unit(0.25, "pt")),
          axis.ticks.x = element_blank(),
          axis.line = element_line(linewidth = unit(0.25, "pt")),
          axis.title = element_text(size = unit(12, "pt"), colour = "black"),
          axis.text = element_text(size = unit(12, "pt"), colour = "black"),
          axis.text.x = element_blank(),
          # axis.line = element_line(colour="black"),
          aspect.ratio = 0.1
    )
  
  # curcular
  df_pie1 <- df2 %>%
    filter(is.na(Family) != TRUE & Family != "No_annotation") %>%
    group_by(Family) %>%
    summarise(n_sample = n()) %>%
    arrange(desc(row_number())) %>%
    mutate(cumulative = cumsum(n_sample),
           mid_point = cumulative - 0.5 * n_sample) %>%
    arrange(desc(row_number())) %>%
    mutate(Family = factor(Family, levels = Family))
  
  na_no_annotation_num <- df2 %>% filter(is.na(Family) == TRUE | Family == "No_annotation") %>% nrow()
  p1_circ <- df_pie1 %>% 
    ggplot(aes(x = "", y = n_sample, label = Family, fill = Family)) +
    geom_bar(stat = "identity") +
    coord_polar(theta = "y", direction = -1) + #時計回り
    geom_label_repel(
      aes(y = mid_point, label = paste(Family, n_sample)), 
      fill = "white", color = "black", size = 5, 
      label.size = 1, force = 1, nudge_x = 0.5
    ) +
    scale_fill_manual(values = color_list) +
    labs(caption = paste0("No Family name in CIS-BP; ", na_no_annotation_num)) +
    theme(plot.title = element_text(size = unit(12, "pt")), 
          #legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.ticks = element_line(colour = "black", linewidth = unit(0.25, "pt")),
          axis.line = element_line(linewidth = unit(0.25, "pt")),
          axis.title = element_text(size = unit(12, "pt"), colour = "black"),
          axis.text = element_text(size = unit(12, "pt"), colour = "black"),
          axis.text.x = element_blank(),
          # axis.line = element_line(colour="black"),
          aspect.ratio = 1
    )
  
  
  # per Cell type ------
  p2 <- df %>% group_by(Cell_type) %>% summarise(n_sample = n()) %>%
    ggplot(aes(x = reorder(Cell_type, -n_sample), y = n_sample, label = Cell_type)) +
    geom_bar(stat = "identity") +
    ggrepel::geom_label_repel() +
    labs(
      x = "Cell type",
      y = "Number of TF ChIP-seq"
    )+
    scale_fill_manual(values = color_list) +
    theme(plot.title = element_text(size = unit(12, "pt")), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.ticks = element_line(colour = "black", linewidth = unit(0.25, "pt")),
          axis.ticks.x = element_blank(),
          axis.line = element_line(linewidth = unit(0.25, "pt")),
          axis.title = element_text(size = unit(12, "pt"), colour = "black"),
          axis.text = element_text(size = unit(12, "pt"), colour = "black"),
          axis.text.x = element_blank(),
          # axis.line = element_line(colour="black"),
          aspect.ratio = 0.1
    )
  
  # p2_label <- df %>% group_by(Cell_type, Cell_type_class) %>% summarise(n_sample = n()) %>%
  #   ggplot(aes(x = reorder(Cell_type, -n_sample), y = n_sample, label = Cell_type, fill = Cell_type_class)) +
  #   geom_bar(stat = "identity") +
  #   ggrepel::geom_label_repel() +
  #   labs(
  #     x = "Cell type",
  #     y = "Number of TF ChIP-seq"
  #   )+
  #   scale_fill_manual(values = color_list) +
  #   theme(plot.title = element_text(size = unit(12, "pt")), 
  #         #legend.position = "none",
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         panel.background = element_blank(), 
  #         axis.ticks = element_line(colour = "black", linewidth = unit(0.25, "pt")),
  #         axis.line = element_line(linewidth = unit(0.25, "pt")),
  #         axis.title = element_text(size = unit(12, "pt"), colour = "black"),
  #         axis.text = element_text(size = unit(12, "pt"), colour = "black"),
  #         axis.text.x = element_blank(),
  #         # axis.line = element_line(colour="black"),
  #         aspect.ratio = 1
  #   )
  
  #circular plot
  df_pie2 <- df %>%
    group_by(Cell_type_class) %>%
    summarise(n_sample = n()) %>%
    mutate(cumulative = cumsum(n_sample)) %>%  # 累積合計
    arrange(-cumulative) %>%
    mutate(cumulative = cumsum(n_sample)) %>%
    mutate(mid_point = cumulative - 0.5 * n_sample)  # ラベルを配置するy座標
  
  p2_circ <- df_pie2 %>%
    ggplot(aes(x = "", y = n_sample, label = Cell_type_class, fill = Cell_type_class)) +
    geom_bar(stat = "identity") +
    coord_polar(theta = "y", direction = -1) +
    geom_label_repel(
      aes(y = mid_point, label = paste(Cell_type_class, n_sample)), 
      fill = "white", color = "black", size = 5, 
      label.size = 1, force = 1, nudge_x = 0.5
    ) +
    scale_fill_manual(values = color_list) +
    theme(plot.title = element_text(size = unit(12, "pt")), 
          #legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.ticks = element_line(colour = "black", linewidth = unit(0.25, "pt")),
          axis.line = element_line(linewidth = unit(0.25, "pt")),
          axis.title = element_text(size = unit(12, "pt"), colour = "black"),
          axis.text = element_text(size = unit(12, "pt"), colour = "black"),
          axis.text.x = element_blank(),
          # axis.line = element_line(colour="black"),
          aspect.ratio = 1
    )
  
  return(list(p1, p1_circ, p2, p2_circ))
}