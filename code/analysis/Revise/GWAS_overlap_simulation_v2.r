GWAS_overlap_simulation_v2 <- function(tgt_ctc, number, df, GWAS_df, simulation_ID, df_fam){
  
  print(paste0("Target Cell type class : ", tgt_ctc))
  
  df_distinct <- df %>% 
    filter(Cell_type_class == tgt_ctc) %>% 
    select(ID, Antigen, Cell_type_class, Cell_type, year, month, day) %>% 
    distinct()  %>%
    as_tibble()
  print(paste0("Total unique IDs : ", length(unique(df_distinct$ID))))
  print(paste0("Total unique TFs : ", length(unique(df_distinct$Antigen))))
  
  tgt_TFs <- df_distinct$Antigen %>% unique() %>% as.character()
  tgt_df_fam <- df_fam %>% filter(Antigen %in% tgt_TFs) %>% select(Antigen, Family) %>% distinct()
  df_distinct <- df_distinct %>% left_join(tgt_df_fam, by = "Antigen")
  
  date_rand <- df_distinct %>%
    select(year, month, day) %>%
    sample_frac(1) #shuffled date 
  df_distinct_IDs <- df_distinct %>% select(ID, Family, Antigen, Cell_type_class, Cell_type) 
  df_distinct2 <- c(df_distinct_IDs, date_rand) %>% as_tibble()
  
  overlap_df2_selected <- df %>% filter(Cell_type_class == tgt_ctc) %>% 
    select(-c(Cell_type_class, Cell_type, Antigen,  year, month, day))
  overlap_df2_selected_rand <- overlap_df2_selected %>% left_join(df_distinct2, by = "ID")
  unique(overlap_df2_selected_rand$Cell_type_class)
  overlap_df2_selected_rand$year <- as.numeric(overlap_df2_selected_rand$year)
  overlap_df2_selected_rand$month <- as.numeric(overlap_df2_selected_rand$month)
  overlap_df2_selected_rand$day <- as.numeric(overlap_df2_selected_rand$day)
  year_list <- overlap_df2_selected_rand %>% drop_na(year) %>% arrange(year) %>% .$year %>% unique() %>% as.numeric()
  
  rs_list <- c()
  TF_list <- c()
  Familiy_list <- c()
  
  ct_plot_df <- data.frame(matrix(nrow = length(year_list), ncol = 5))
  num <- 1
  for (tgt_year in year_list) {
    
    tgt_rs <- overlap_df2_selected_rand %>% filter(year == tgt_year) %>% .$rs %>% as.character()
    rs_list <- c(rs_list, tgt_rs)
    cum_uniq_rs <- length(unique(rs_list))
    
    tgt_TF <- overlap_df2_selected_rand %>% filter(year == tgt_year) %>% .$Antigen %>% as.character()
    TF_list <- c(TF_list, tgt_TF)
    cum_uniq_TF <- length(unique(TF_list))
    
    tgt_Familiy <- overlap_df2_selected_rand %>% filter(year == tgt_year) %>% .$Family %>% as.character()
    Familiy_list <- c(Familiy_list, tgt_Familiy)
    cum_uniq_Familiy <- length(unique(Familiy_list))
    
    add_row <- tibble(Cell_type_class = tgt_ctc, 
                      year = tgt_year, 
                      cum_uniq_rs = cum_uniq_rs, 
                      cum_uniq_TF = cum_uniq_TF, 
                      cum_uniq_Fam = cum_uniq_Familiy)
    ct_plot_df[num,] <- add_row
    num <- num + 1
  }
  
  total_rs_num <- length(unique(GWAS_df$SNPS))
  
  colnames(ct_plot_df) <- c("Cell_type_class", "year", "cum_uniq_rs", "cum_uniq_TF",  "cum_uniq_Fam")
  ct_plot_df2 <- ct_plot_df %>% mutate(GWAS_cover_ratio = cum_uniq_rs/total_rs_num) %>%
    mutate(simulation_ID = simulation_ID)
  
  p <- ct_plot_df2 %>% 
    ggplot(aes(x = year, y = GWAS_cover_ratio))  +
    geom_point() +
    geom_line()+
    ggtitle(paste0(tgt_ctc, " simu", simulation_ID)) +
    labs(
      x = "Year",
      y = "GWAS-SNP cover ratio"
    )+
    #geom_hline(yintercept = 0.1, linetype="dashed", color = "red") +
    ggtitle(tgt_ctc) +
    theme(
      aspect.ratio = 1,
      plot.background = element_blank(),
      plot.title = element_text(size = unit(12, "pt")),
      plot.subtitle = element_text(size = unit(10, "pt")),
      plot.caption = element_text(size = unit(8, "pt")),
      plot.tag = element_text(size = unit(18, "pt")),
      axis.ticks = element_line(colour = "black", linewidth = unit(0.25, "pt")),
      axis.line = element_line(linewidth = unit(0.5, "pt")),
      axis.title = element_text(size = unit(10, "pt")),
      axis.text = element_text(size = unit(8, "pt"), colour = "black"),
      axis.text.x = element_text(size= unit(8, "pt"), color = "black", angle = 45, hjust = 1),
      panel.grid.major = element_line(linewidth = unit(0.5, "pt")), 
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = unit(0.25, "pt")),
      # panel.background = element_rect(fill = NA, colour = "black", linewidth = unit(0.5, "pt")),
      panel.background = element_blank(),
      # panel.ontop = TRUE,
      legend.key = element_blank(),
      legend.title = element_text(size = unit(10, "pt")),
      legend.text = element_text(size = unit(10, "pt")),
      strip.text = element_text(colour = "black", size = unit(8, "pt")),
      strip.background = element_rect(fill = NA, color = NA)
    )
  
  return(list(overlap_df2_selected_rand, ct_plot_df2, p))
}