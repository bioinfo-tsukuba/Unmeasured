rm(list=ls())
# これはもうやらなくて大丈夫なところ ----
stat_hu2 <- read_tsv( "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix_md/DEG_all_human_annotated.tsv")
TF <- stat_hu2$TF.x %>% unique()

# tissue annotation add
tissue_df <- read_csv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/KnockTFv2_tissue_type.csv")
dim(tissue_df)
tissue_df2 <- tissue_df %>% filter(TF %in% TF)  %>% select(`Tissue type`, `Pubmed ID`, `Dataset ID`) %>% distinct()
colnames(tissue_df2) <- c( "Tissue", "Pubmed_ID", "ID")

stat_hu3 <- stat_hu2 %>% left_join(tissue_df2, by = "ID")
dim(stat_hu2)
dim(stat_hu3)

write_tsv(stat_hu3, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix_md/DEG_all_human_annotated_tissue.tsv")




# TissueとCell type classを対応付ける -------
stat_hu3 <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix_md/DEG_all_human_annotated_tissue.tsv")
Knock_tissue <- unique(stat_hu3$Tissue)
write(Knock_tissue, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/tissue_list.txt")

Chip_df1 <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/experimentList.tab", col_names = F)
Chip_df2 <- Chip_df1 %>%  select(X1, X2, X3, X4, X5, X6)
colnames(Chip_df2) <- c("ID", "Genome", "Antigen_class", "Antigen", "Cell_type_class", "Cell_type")
anno_df <- Chip_df2 %>% filter(Genome == "hg38" & Antigen_class == "TFs and others")

Chip_ctc <- anno_df$Cell_type_class %>% unique()
share_ctc <- intersect(Knock_tissue, Chip_ctc)
#length(knock_ctc)
length(Chip_ctc)
length(share_ctc)

tissue_ctc_df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/tissue_Chip_ctc.txt", col_names = F)
View(tissue_ctc_df)
colnames(tissue_ctc_df) <- c("Tissue", "Cell_type_class")

stat_hu4 <- stat_hu3 %>% left_join(tissue_ctc_df, by = "Tissue")
View(stat_hu4[1:10,])
dim(stat_hu4)
length(unique(stat_hu4$TF.x))
write_tsv(stat_hu4, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix_md/DEG_all_human_annotated_tissue_ctc.tsv")


# unmeasured だったものとそうではなかった組み合わせでDEGの数をplot -------
stat_hu4 <- read_tsv( "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix_md/DEG_all_human_annotated_tissue_ctc.tsv")
Chip_df <- readRDS("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/TF_CTC_sample_num_new.rds")
knock_tf <- stat_hu4$TF.x %>% unique()
Chip_tf <- anno_df$Antigen %>% unique()
share_tf <- intersect(knock_tf, Chip_tf)
length(knock_tf)
length(Chip_tf)
length(share_tf)

rowname <-rownames(Chip_df)
Chip_df2 <- Chip_df %>% as_tibble() %>% mutate(tf = rowname)
Chip_df_selected <- Chip_df2 %>% filter(tf %in% share_tf) %>% select(all_of(share_ct), tf)
dim(Chip_df_selected)
rowname_selected <- Chip_df_selected$tf
Chip_df_selected2 <- Chip_df_selected %>% select(-tf) %>% as.matrix()
rownames(Chip_df_selected2) <- rowname_selected
dim(Chip_df_selected2)
View(Chip_df_selected2[1:100,])

measured_tf <- c()
measured_ctc <- c()
unmeasured_tf <- c()
unmeasured_ctc <- c()
measured_DEG_count <- c()
unmeasured_DEG_count <- c()
for(j in 1:nrow(Chip_df_selected2)){
  print(paste(j, "/", nrow(Chip_df_selected2)))
  for (i in 1:ncol(Chip_df_selected2)) {
    tgt_tf <- rownames(Chip_df_selected2)[j]
    tgt_ctc <- colnames(Chip_df_selected2)[i]
    tgt_chip_num <- Chip_df_selected2[j, i]
    tgt_row <- stat_hu4 %>% filter(TF.x == tgt_tf & Cell_type_class == tgt_ctc) 
    
    if(tgt_chip_num > 1){ #measured
      if(nrow(tgt_row)>0){
        tgt_num_DEG_2FC <- tgt_row$`2FC_up` + tgt_row$`2FC_down`
        measured_tf <- c( measured_tf, tgt_row$TF.x)
        measured_ctc <- c( measured_ctc, tgt_row$Cell_type_class)
        measured_DEG_count <- c(measured_DEG_count, tgt_num_DEG_2FC)
      }
    }else{
      if(nrow(tgt_row)>0){ #unmeasured
        tgt_num_DEG_2FC <- tgt_row$`2FC_up` + tgt_row$`2FC_down`
        unmeasured_tf <- c( unmeasured_tf, tgt_row$TF.x)
        unmeasured_ctc <- c( unmeasured_ctc, tgt_row$Cell_type_class)
        unmeasured_DEG_count <- c(unmeasured_DEG_count, tgt_num_DEG_2FC)
      }
    }
    rm(tgt_row)
  }
}

tmp1 <- tibble(tf = measured_tf, ctc = measured_ctc, count = measured_DEG_count, label = rep("measured", length(measured_DEG_count)))
tmp2 <- tibble(tf = unmeasured_tf, ctc = unmeasured_ctc, count = unmeasured_DEG_count, label = rep("unmeasured", length(unmeasured_DEG_count)))
tibble <- rbind(tmp1, tmp2) %>% distinct()
tibble
View(tibble)

tibble %>% ggplot(aes(x = label, y = count, color = label)) +
  geom_violin() +
  #geom_boxplot()+
  geom_point(size = 1)+
  xlab("label")+
  ylab("Number of 2FC DEGs")+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_line(colour="gray"),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=15,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=10,face="bold"),
        aspect.ratio = 1
  )
write_tsv(tibble, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/tib_2FC_DEGs_ctc_tissue_unm.tsv")

# measured
measured_df <- tibble %>% filter(label == "measured")
dim(measured_df)
summary(measured_df$count)
var(measured_df$count)

# unmeasured
unmeasured_df <- tibble %>% filter(label == "unmeasured")
dim(unmeasured_df)
summary(unmeasured_df$count)
var(unmeasured_df$count)

# 検定 (wilcoxonの順位和検定、対応なし(=マンホイットニーのU検定))
result <- wilcox.test(measured_df$count, unmeasured_df$count)
result
