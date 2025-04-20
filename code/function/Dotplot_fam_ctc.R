Dotplot_fam_ctc <- function(){
  
  df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/SRX_date_20231004.tsv")
  df <- df %>% arrange(year, month, day, time) %>% mutate(order = 1:nrow(df)) %>% 
    filter(Cell_type_class != "Unclassified" & Cell_type_class != "No description") %>%
    filter(Antigen != "GFP" & Antigen != "Epitope tags")
  total_sample <- length(unique(df$SRX))
  
  df2 <- df %>% select(-c(year, month, day)) %>%
    group_by(Antigen, Cell_type_class) %>%
    summarise(n_sample = n())
  
  library(ggside)
  df_fam <- readRDS("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS_paper_public/data/Fig2/obj/df_fam.rds")
  tgt_IDs <- df$SRX %>% unique() %>% as.character()
  tgt_df_fam <- df_fam %>% filter(ID %in% tgt_IDs) %>% select(Antigen, Family) %>% distinct()
  
  df3 <- df2 %>% left_join(tgt_df_fam, by = "Antigen") %>% 
    group_by(Family, Cell_type_class) %>% 
    summarise(n_sample = sum(n_sample))
  
  tmp <- df3 %>%  filter(Family != "Unknown" & Family != "NA" & Family != "No_annotation") %>%
    mutate(percentage = (n_sample / total_sample)*100) 
  length(unique(tmp$Family))
  
  tmp2 <- tmp %>% group_by(Cell_type_class) %>% summarise(n_sum_sample_ctc = sum(n_sample))
  tmp3 <- tmp %>% left_join(tmp2, by = "Cell_type_class")
  
  p <- tmp3 %>% 
    ggplot(aes(x = reorder(Family, -n_sample), y= reorder(Cell_type_class, n_sum_sample_ctc), size = percentage, label = n_sample)) + 
    geom_point(aes(color = n_sample)) +
    #geom_text(color = "white", size = 2) +
    geom_text(color = "black", size = 2) +
    scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,max(tmp$n_sample)), oob = scales::squish) +
    xlab("TF family") +
    ylab("Cell type class") +
    theme(
      aspect.ratio = 0.7,
      plot.background = element_blank(),
      plot.title = element_text(size = unit(12, "pt")),
      plot.subtitle = element_text(size = unit(10, "pt")),
      plot.caption = element_text(size = unit(8, "pt")),
      plot.tag = element_text(size = unit(18, "pt")),
      axis.ticks = element_line(colour = "black", linewidth = unit(0.5, "pt")),
      axis.line = element_blank(), # element_line(linewidth = unit(0.5, "pt")),
      axis.title = element_text(size = unit(15, "pt")),
      axis.text.x = element_text(size= unit(15, "pt"), color = "black", angle = 90, hjust = 1),
      axis.text.y = element_text(size= unit(15, "pt"), color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = unit(0.25, "pt")),
      panel.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(size = unit(10, "pt")),
      legend.text = element_text(size = unit(10, "pt")),
      strip.text = element_text(colour = "black", size = unit(8, "pt")),
      strip.background = element_rect(fill = NA, color = NA))
  
  p2 <- p + 
    geom_xsidecol(aes(x = reorder(Family, -n_sample), y = n_sample), data = tmp, stat = "identity", width = 0.5, fill = "blue4") +
    geom_ysidecol(aes(x = n_sample, y = Cell_type_class), data = tmp, stat = "identity", width = 0.5, fill = "blue4") +
    theme(
      legend.position = "none"
    ) 
  
  return(p2)
  
}