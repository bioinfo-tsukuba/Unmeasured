Year_SRX <- function(){
  
  library(tidyverse)
  library(ggrepel)
  library(RColorBrewer)
  library(ggh4x)
  
  # 20231004 version
  df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/SRX_date_20231004.tsv")
  df <- df %>% arrange(year, month, day, time) %>% mutate(order = 1:nrow(df)) %>% 
    filter(Cell_type_class != "Unclassified" & Cell_type_class != "No description") %>%
    filter(Antigen != "GFP" & Antigen != "Epitope tags")
  
  # Year - Number of SRX ----
  p1 <- df %>% #filter(year != 2023) %>% 
    group_by(year) %>% summarise(SRX_num =n()) %>%
    ggplot(aes(x = factor(year), y = SRX_num))+
    geom_point(size = 2)+
    labs(
      x = "Year",
      y = "Total number of SRX (TF ChIP-seq)"
    )+
    theme(plot.title = element_text(size = unit(12, "pt"), face="bold"), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.ticks = element_line(colour = "black", linewidth = unit(0.25, "pt")),
          axis.line = element_line(linewidth = unit(0.25, "pt")),
          axis.title = element_text(size = unit(12, "pt"), colour = "black", face = "bold"),
          axis.text = element_text(size = unit(12, "pt"), colour = "black"),
          axis.text.x =element_text(size = unit(12, "pt"), colour = "black", angle = 45, hjust = 1),
          # axis.line = element_line(colour="black"),
          aspect.ratio = 1
    )
  
  # month-Number of SRX ----
  p2 <- df %>% group_by(year, month) %>% summarise(SRX_num =n()) %>% 
    unite("year_month", c(year, month), sep = "_") %>%
    ggplot(aes(x = factor(year_month), y = SRX_num)) +
    geom_point(size = 2) +
    labs(
      x = "Time",
      y = "Total number of SRX (TF ChIP-seq)"
    )+
    theme(plot.title = element_text(size = unit(12, "pt"), face="bold"), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.ticks = element_line(colour = "black", linewidth = unit(0.25, "pt")),
          axis.line = element_line(linewidth = unit(0.25, "pt")),
          axis.title = element_text(size = unit(12, "pt"), colour = "black", face = "bold"),
          axis.text = element_text(size = unit(12, "pt"), colour = "black"),
          axis.text.x =element_text(size = unit(0, "pt"), colour = "black", angle = 45, hjust = 1),
          # axis.line = element_line(colour="black"),
          aspect.ratio = 1
    )
  
  # yearごとの累積SRX数 ----
  df2 <- df %>% arrange(year, month, day, time) %>% 
    group_by(year) %>% summarise(SRX_num =n()) 
  cum_sum <- cumsum(df2$SRX_num) #累積和
  df2 <- df2 %>% mutate(cum = cum_sum)
  
  p3 <- df2 %>%
    ggplot(aes(x = factor(year), y = cum)) +
    geom_point(size = 2) +
    xlab("Year")+
    ylab("Total number of SRX (TF ChIP-seq)")+
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
          axis.text.x =element_text(size = unit(30, "pt"), colour = "black", face = "bold", angle = 45, hjust = 1),
          # axis.line = element_line(colour="black"),
          aspect.ratio = 1
    ) 
  
  # 新規TFか否か ----
  num_list <- c()
  label_list <- c()
  for (i in 1:nrow(df)) {
    print(paste0(i, "/", nrow(df)))
    target_df <- df[1:i,]
    tmp_num <- length(unique(target_df$Antigen))
    
    if(i == 1){
      num_list <- tmp_num
    }else{
      num_list <- c(num_list, tmp_num)
    }
    
    # 新規TFかどうかの判定
    tmp <- intersect(unique(target_df[1:(i-1),]$Antigen), target_df[i,]$Antigen)
    if(length(tmp) == 0){
      if(i == 1){
        label_list <- target_df[i,]$Antigen
      }else{
        label_list <- c(label_list, target_df[i,]$Antigen)
      }
    }else{
      if(i == 1){
        label_list <- target_df[i,]$Antigen
      }else{
        label_list <- c(label_list, "")
      }
    }
  }
  df3 <- df %>% mutate(TF_uniq_num = num_list, label = label_list)
  
  # Time- 新規TF数 -----
  # CIS-BPで登録されているようなDBDをもつTF以外にもChIP-Atlasには含まれているため、TF数が1800を超えていると考えられる
  p4 <- df3 %>% group_by(year) %>% summarise(TF_uniq_num_max = max(TF_uniq_num)) %>%
    ggplot(aes(x = factor(year), y = TF_uniq_num_max)) +
    geom_point(size = 2) +
    #geom_hline(yintercept=1639, linetype="dashed",size=1) +
    xlab("Year")+
    ylab("Total number of unique TFs")+
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
          axis.text.x =element_text(size = unit(30, "pt"), colour = "black", face = "bold", angle = 45, hjust = 1),
          # axis.line = element_line(colour="black"),
          aspect.ratio = 1
    ) 
    
  p5 <- df3 %>%
    ggplot(aes(x = order, y = TF_uniq_num, label = label)) +
    geom_point(size = 1) +
    geom_label_repel(max.overlaps = 20) +
    #geom_hline(yintercept=1639, linetype="dashed",size=1) +
    xlab("Time")+
    ylab("Total number of unique TFs")+
    theme(plot.title = element_text(size = unit(12, "pt"), face="bold"), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.ticks = element_line(colour = "black", linewidth = unit(0.25, "pt")),
          axis.line = element_line(linewidth = unit(0.25, "pt")),
          axis.title = element_text(size = unit(12, "pt"), colour = "black", face = "bold"),
          axis.text = element_text(size = unit(12, "pt"), colour = "black"),
          axis.text.x =element_text(size = unit(0, "pt"), colour = "black", angle = 45, hjust = 1),
          # axis.line = element_line(colour="black"),
          aspect.ratio = 1
    )
  
  # TFごとにplot -------
  TF_list_order <- df3 %>% filter(Antigen != "Epitope tags" & Antigen != "GFP") %>% 
    group_by(Antigen) %>% summarise(n = n()) %>% 
    arrange(desc(n)) %>% .$Antigen
  top_TF_list <- TF_list_order[1:30]
  
  df4 <- c()
  for(tgt_tf in top_TF_list){
    print(tgt_tf)
    target_df <- df3 %>% filter(Antigen == tgt_tf) %>% group_by(year) %>% summarise(SRX_num =n()) 
    cum_sum <- cumsum(target_df$SRX_num) #累積和
    target_df2 <- target_df %>% mutate(cum_per_TF = cum_sum, Antigen = tgt_tf)
    
    if(tgt_tf == top_TF_list[1]){
      df4 <- target_df2
    }else{
      df4 <- rbind(df4, target_df2)
    }
  }
    
  color_list <- c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG"))
  p6 <- df4 %>% filter(Antigen != "Epitope tags" & Antigen %in% top_TF_list[1:20]) %>%
    ggplot(aes(x = year, y = cum_per_TF, color = Antigen)) +
    geom_point(size = 2) +
    geom_line() +
    scale_color_manual(values = color_list) +
    xlab("Year")+
    ylab("Total number of SRX (TF ChIP-seq)")+
    theme(plot.title = element_text(size = unit(15, "pt"), face="bold"), 
          legend.position = c(0.3,0.75),
          legend.text = element_text(size = unit(15, "pt")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.ticks = element_line(colour = "black", linewidth = unit(0.5, "pt")),
          axis.line = element_line(linewidth = unit(0.5, "pt")),
          axis.title = element_text(size = unit(30, "pt"), colour = "black", face = "bold"),
          axis.text = element_text(size = unit(30, "pt"), colour = "black", face = "bold"),
          axis.text.x =element_text(size = unit(30, "pt"), colour = "black", face = "bold", angle = 45, hjust = 1),
          # axis.line = element_line(colour="black"),
          aspect.ratio = 1
    ) +
    guides(colour = guide_stringlegend(ncol = 2)) # legendを2列にする
  
  
  # Cell type classごとに、新規TFのplot ----
  df5 <- c()
  for(tgt_ctc in unique(df3$Cell_type_class)){
    print(tgt_ctc)
    target_df <- df3 %>% filter(Cell_type_class == tgt_ctc) 
    
    num_list <- c()
    label_list <- c()
    for (i in 1:nrow(target_df)) {
      sub_target_df <- target_df[1:i,]
      tmp_num <- length(unique(sub_target_df$Antigen))
      
      if(i == 1){
        num_list <- tmp_num
      }else{
        num_list <- c(num_list, tmp_num)
      }
      
      # 新規TFかどうかの判定
      tmp <- intersect(unique(sub_target_df[1:(i-1),]$Antigen), sub_target_df[i,]$Antigen)
      if(length(tmp) == 0){
        if(i == 1){
          label_list <- sub_target_df[i,]$Antigen
        }else{
          label_list <- c(label_list, sub_target_df[i,]$Antigen)
        }
      }else{
        if(i == 1){
          label_list <- sub_target_df[i,]$Antigen
        }else{
          label_list <- c(label_list, "")
        }
      }
    }
    target_df2 <- target_df %>% mutate(TF_uniq_num_per_ctc = num_list, label_list_per_ctc = label_list)
    if(tgt_ctc == unique(df3$Cell_type_class)[1]){
      df5 <- target_df2
    }else{
      df5 <- rbind(df5, target_df2)
    }
  }
  
  color_list <- c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG"))
  p7<- df5 %>% group_by(year, Cell_type_class) %>% summarise(TF_uniq_num_max = max(TF_uniq_num_per_ctc)) %>%
    ggplot(aes(x = year, y = TF_uniq_num_max, color = Cell_type_class)) +
    geom_point(size = 3) +
    geom_line() +
    xlab("Year")+
    ylab("Total number of unique TFs")+
    scale_color_manual(values = color_list) +
    theme(plot.title = element_text(size = unit(15, "pt"), face="bold"), 
          legend.position = c(0.3,0.75),
          legend.title = element_text(size = unit(15, "pt")),
          legend.text = element_text(size = unit(15, "pt")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.ticks = element_line(colour = "black", linewidth = unit(0.5, "pt")),
          axis.line = element_line(linewidth = unit(0.5, "pt")),
          axis.title = element_text(size = unit(30, "pt"), colour = "black", face = "bold"),
          axis.text = element_text(size = unit(30, "pt"), colour = "black", face = "bold"),
          axis.text.x =element_text(size = unit(30, "pt"), colour = "black", face = "bold", angle = 45, hjust = 1),
          # axis.line = element_line(colour="black"),
          aspect.ratio = 1
    ) +
    guides(colour = guide_stringlegend(ncol = 2)) # legendを2列にする
  
  return(list(p1, p2, p3, p4, p5, p6, p7))
}