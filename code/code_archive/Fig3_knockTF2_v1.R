rm(list=ls())
stat_hu2 <- read_tsv( "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix_md/DEG_all_human_annotated.tsv")

# correcter_pval < 0.05をDEGとした
df1 <- stat_hu2 %>% filter(corrected_p < 0.05) %>% select(TF.x, Target_Gene) %>% distinct() %>% group_by(TF.x) %>% summarise(DEG_num = n())
df1 %>% ggplot(aes(x = DEG_num, fill = ))   +
  geom_histogram(bins = 50)+
  scale_fill_manual("blue")+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_line(colour="gray"),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        aspect.ratio = 0.7
  )

df1 %>% arrange(desc(DEG_num)) %>%
  ggplot(aes(x = reorder(TF.x, DEG_num), y =DEG_num)) +
  geom_point(size = 0.7)+
  scale_fill_manual("blue")+
  xlab("TF") +
  ylab("Number of DEGs") +
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_line(colour="gray"),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=3,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        aspect.ratio = 0.7
  )

# ChIP-Atlas
# cell typeの対応検討
Chip_df1 <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/experimentList.tab", col_names = F)
Chip_df2 <- Chip_df1 %>%  select(X1, X2, X3, X4, X5, X6)
colnames(Chip_df2) <- c("ID", "Genome", "Antigen_class", "Antigen", "Cell_type_class", "Cell_type")
anno_df <- Chip_df2 %>% filter(Genome == "hg38" & Antigen_class == "TFs and others")

knock_ct <- stat_hu2$Biosample_name %>% unique()
Chip_ct <- anno_df$Cell_type %>% unique()
share_ct <- intersect(knock_ct, Chip_ct)
length(knock_ct)
length(Chip_ct)
length(share_ct)

Chip_df <- readRDS("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/TF_CT_sample_num_new.rds")
knock_tf <- stat_hu2$TF.x %>% unique()
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
measured_ct <- c()
unmeasured_tf <- c()
unmeasured_ct <- c()
measured_DEG_count <- c()
unmeasured_DEG_count <- c()
for(j in 1:nrow(Chip_df_selected2)){
  print(paste(j, "/", nrow(Chip_df_selected2)))
  for (i in 1:ncol(Chip_df_selected2)) {
    tgt_tf <- rownames(Chip_df_selected2)[j]
    tgt_ct <- colnames(Chip_df_selected2)[i]
    tgt_chip_num <- Chip_df_selected2[j, i]
    tgt_row <- stat_hu2 %>% filter(TF.x == tgt_tf & Biosample_name == tgt_ct) 
    
    if(tgt_chip_num > 1){ #measured
      if(nrow(tgt_row)>0){
        tgt_num_DEG_2FC <- tgt_row$`2FC_up` + tgt_row$`2FC_down`
        measured_tf <- c( measured_tf, tgt_row$TF.x)
        measured_ct <- c( measured_ct, tgt_row$Biosample_name)
        measured_DEG_count <- c(measured_DEG_count, tgt_num_DEG_2FC)
      }
    }else{
      if(nrow(tgt_row)>0){ #unmeasured
        tgt_num_DEG_2FC <- tgt_row$`2FC_up` + tgt_row$`2FC_down`
        unmeasured_tf <- c( unmeasured_tf, tgt_row$TF.x)
        unmeasured_ct <- c( unmeasured_ct, tgt_row$Biosample_name)
        unmeasured_DEG_count <- c(unmeasured_DEG_count, tgt_num_DEG_2FC)
      }
    }
    rm(tgt_row)
  }
}

tmp1 <- tibble(tf = measured_tf, ct = measured_ct, count = measured_DEG_count, label = rep("measured", length(measured_DEG_count)))
tmp2 <- tibble(tf = unmeasured_tf, ct = unmeasured_ct, count = unmeasured_DEG_count, label = rep("unmeasured", length(unmeasured_DEG_count)))
tibble <- rbind(tmp1, tmp2) %>% distinct()
tibble
View(tibble)

tibble %>% ggplot(aes(x = label, y = count, fill = label)) +
  geom_violin() +
  geom_point()+
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
write_tsv(tibble, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/tib_2FC_DEGs_unm.tsv")
