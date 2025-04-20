FPKM_threshold_histogram <- function(){
 
  library(tidyverse)
  library(patchwork)
  
  
  mt_RNA2 <- readRDS("/Users/saeko/Unmeasured/data/mt_RNA2.rds")
  
  # tissue毎にthreshold_FPKMを設定してhistogram内に線を入れる
  p_patch <- c()
  for (j in 1:ncol(mt_RNA2)) {
    threshold_025 <- quantile(mt_RNA2[, j], 1-0.25) #第三四分位数を求める
    threshold_05 <- quantile(mt_RNA2[, j], 1-0.5) 
    #tgt_Chip <- mt_Chip6[i,j]
    tgt_RNA <- mt_RNA2[,j]
    tgt_tib <- tgt_RNA %>% as_tibble() %>% mutate(TF = rownames(mt_RNA2)) 
    colnames(tgt_tib) <- c("FPKM", "TF")
    tgt_p <- tgt_tib %>%
      ggplot(aes(x = FPKM)) +
      geom_histogram() +
      geom_vline(xintercept = threshold_025, linetype = "solid", color = "red", size = 1.2) +
      #annotate("text", x = threshold + 1.8, y = 200, label = paste0("Expressed threshold (", tgt_threshold, ")"), color = "red", vjust = -0.5) +
      geom_vline(xintercept = threshold_05, linetype = "solid", color = "blue", size = 1.2) +
      geom_vline(xintercept = 2.6, linetype = "solid", color = "green3", size = 1.2) +
      #annotate("text", x = 2.6 + 1.2, y = 100,  label = "TF marker", color = "green3", vjust = -0.5)+
      ggtitle(colnames(mt_RNA2)[j]) +
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = c(0.5, 1.1),
            legend.direction = "horizontal",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour="black"),
            axis.text=element_text(size=10, color = "black"),
            axis.text.x =element_text(size=10, color = "black"),
            axis.text.y =element_text(size=10, color = "black"),
            axis.title=element_text(size=10, color = "black"),
            aspect.ratio = 0.5
      )
    tgt_p
    if(j == 1){
      p_patch <- tgt_p
    }else{
      p_patch <- p_patch + tgt_p
    }
  }
  return(p_patch) 
}