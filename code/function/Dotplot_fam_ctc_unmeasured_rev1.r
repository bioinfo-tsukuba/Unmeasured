Dotplot_fam_ctc_unmeasured_rev1 <- function(){
  
  # 「１度もChIP-seqされていないTF」も含めたバージョン
  # Revise, TF軸でexpressed TFを決めたver.
  
  mt_unmeasured_tib <- readRDS("/Users/saeko/Unmeasured/data/tib_unmeasured_stepminer_Ctc.rds")
  
  tib_dotplot <- mt_unmeasured_tib %>% 
    filter(is.na(Family) == FALSE & Family != "No_annotation" & Family != "Unknown") %>%
    group_by(Family, Cell_type_class) %>%
    summarise(
      total_expressed_TF = sum(expressed, na.rm = TRUE),
      #total_measured_ChIPseq = sum(num_ChIPseq, na.rm = TRUE),
      total_unmeasured_TF = sum(num_unmeasured, na.rm = TRUE),
      ratio_unmeasured = total_unmeasured_TF/total_expressed_TF
    )
  
  p <- tib_dotplot %>% 
    ggplot(aes(x = reorder(Family, -total_unmeasured_TF), 
               y= reorder(Cell_type_class, total_unmeasured_TF), 
               size = ratio_unmeasured, label = total_unmeasured_TF)) + 
    geom_point(aes(color = ratio_unmeasured)) +
    geom_text(color = "white", size = 2) +
    scale_color_gradientn(
      #colours = c("blue", "white",  "red"),  # お好みの色ベクトル
      colours = c("white",  "blue"),
      limits = c(0, max(tib_dotplot$ratio_unmeasured)), 
      oob = scales::squish
    )+
    xlab("TF family") +
    ylab("Cell type class") +
    theme(
      aspect.ratio = 0.8,
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
  
  return(p)
  
}