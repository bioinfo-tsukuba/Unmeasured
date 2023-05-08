# TF marker ----
TF_marker_all <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/TFmarker/All_TFmarkers.txt")
taiou <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/TFmarker/TFmarker_tissue_ChIP-Atlas_Ctc_taiou.tsv")
colnames(taiou) <- c("Tissue_Type", "ChIPAtlas_ctc", "description")
taiou <- taiou %>% select(-description) %>% drop_na(ChIPAtlas_ctc)
colnames(TF_marker_all) <- c("PMID", "Gene_Name","Gene_Type","Cell_Name","Cell_Type","Tissue_Type","Experiment_Type","Experimental_Method","Title","Description_of_Gene","Interacting_Gene","CellOntologyID")
df_all2 <- TF_marker_all %>% left_join(taiou, by = "Tissue_Type") %>% select(PMID, Gene_Name, Gene_Type, Cell_Name, Cell_Type,Tissue_Type, ChIPAtlas_ctc) %>%
  unite("key", c(Gene_Name, ChIPAtlas_ctc), sep = "_") %>% select(-PMID, -Cell_Name, -Cell_Type, -Tissue_Type, -Gene_Type) %>%
  mutate(TFmarker = "marker") %>% distinct()

# ChIP-Atlas ----
# Chip_df <- readRDS("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/TF_CTC_sample_num_new.rds")
# tmp_TF_list <- rownames(Chip_df)
# Chip_df2 <- Chip_df %>% as_tibble() %>% mutate(TF = tmp_TF_list) %>% pivot_longer(-TF, names_to = "Cell_type_class", values_to = "count_ChIPseq")
# join_key2 <- Chip_df2 %>%  unite("key", c(TF, Cell_type_class), sep = "_") %>% .$key
# Chip_df3 <- Chip_df2 %>% mutate(key = join_key2) %>% select(key, `count_ChIPseq`) %>% drop_na(`count_ChIPseq`) %>%
#   mutate(measure = ifelse(`count_ChIPseq` == 0, "unmeasured", "measured"))

mt_RNA <- readRDS("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/RNAseq/RefEx_expression_RNAseq10_human_PRJEB2445_heat.rds")
mt_Chip <- readRDS("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/TF_CTC_sample_num_new.rds")
share_TF <- intersect(rownames(mt_RNA), rownames(mt_Chip)) 
length(share_TF)

mt_RNA <- mt_RNA[share_TF, ]
mt_Chip <- mt_Chip[share_TF, ]

# Antigenの行を2つの行列で揃える
order_RNA<- order(rownames(mt_RNA))
order_Chip<- order(rownames(mt_Chip))
mt_RNA2 <- mt_RNA[order_RNA,]
mt_Chip2 <-mt_Chip[order_Chip,]

# Ctc, tissueを2つの行列で揃える
colnames(mt_Chip2)
colnames(mt_RNA2)
mt_Chip3 <- mt_Chip2 %>% as_tibble() %>% select(-Bone, -Embryo, -Epidermis, -`No description`, -Others, -Pancreas, -Placenta, -`Pluripotent stem cell`, -Unclassified)

# ChIP-AtlasのCell type classアノテーションを、RefEXに合わせる
mt_Chip4 <- mt_Chip3 %>% mutate(muscular = rowSums(mt_Chip3[, c(4, 10)])) %>% select(-Cardiovascular, -Muscle) #CardiovascularとMuscleはmuscular"へ
mt_Chip5 <- mt_Chip4 %>% mutate(reproductive = rowSums(mt_Chip3[, c(6, 12, 13)])) %>% select(-Gonad, -Prostate, -Uterus) #GonadとProstateとUterusは"reproductive"へ

mt_Chip6 <- mt_Chip5 %>% select("Neural", "Breast", "Blood", "Adipocyte", "reproductive", "muscular", "Digestive tract", "Liver", "Lung", "Kidney") 
#colnames(mt_Chip6) <- colnames(mt_RNA2)
mt_Chip6 <- mt_Chip6 %>% as.matrix()
rownames(mt_Chip6) <- rownames(mt_Chip2)


# FPKM  > 2.77(上位25%) 以上でcountが0ならunmeasured(0)、countが1以上ならmeasuredとする(1)---
dim(mt_Chip6)
dim(mt_RNA2)

mt_mesured_unmeasured <- matrix(nrow = nrow(mt_Chip6), ncol = ncol(mt_Chip6))
for (i in 1:nrow(mt_Chip6)) {
  for (j in 1:ncol(mt_Chip6)) {
    tgt_Chip <- mt_Chip6[i,j]
    tgt_RNA <- mt_RNA2[i,j]
    if(tgt_Chip == 0 & tgt_RNA >= 2.77){
      mt_mesured_unmeasured[i,j] <- 0
    }else if(tgt_Chip >= 1){
      mt_mesured_unmeasured[i,j] <- 1
    }else{
      mt_mesured_unmeasured[i,j] <- NA
    }
  }
}

rownames(mt_mesured_unmeasured) <- rownames(mt_Chip6)
colnames(mt_mesured_unmeasured) <- colnames(mt_Chip6)


mt_mesured_unmeasured2 <- mt_mesured_unmeasured %>% as_tibble() %>% mutate(TF = rownames(mt_Chip6)) %>%
  pivot_longer(-TF, names_to = "Cell_type_class", values_to = "measured_unmeasured") %>% drop_na(measured_unmeasured)%>%
  unite("key", c(TF, Cell_type_class), sep = "_")



# join ---
join1 <- mt_mesured_unmeasured2 %>% left_join(df_all2, by = "key") %>% mutate(TFmarker2 = ifelse(is.na(TFmarker) ==TRUE, "no_marker", "marker"))
join2 <- join1 %>% mutate(measure = ifelse(measured_unmeasured == 1, "measured", "unmeasured")) %>% group_by(measure, TFmarker2) %>% summarise(n = n())

join2 %>% #separate(label, c("label_measured", "label_marker"), sep = "_") %>% 
  ggplot(aes(x = measure, y = n, fill = TFmarker2))+
  geom_bar(stat = "identity", position = "fill", width = 0.7)+
  scale_fill_manual(values = c("orange", "gray"))+
  xlab("")+
  ylab("count")+
  geom_text(aes(label = n), position = position_fill(vjust = 0.5), size = 7) +
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        #legend.position = "none",
        panel.grid.major = element_line(color = "gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold", color = "black"),
        axis.text.x =element_text(size=10,face="bold", color = "black"),
        axis.text.y =element_text(size=10,face="bold", color = "black"),
        axis.title=element_text(size=15,face="bold", color = "black"),
        ggside.panel.scale = 0.2,
        aspect.ratio = 0.9
  )

# markerだけどunmeasuredなTF-Ctcのペアの分布
join3 <- join1 %>% separate("key", into = c("TF", "Cell_type_class"), sep = "_") %>% 
  filter(measured_unmeasured == 0 & TFmarker2 == "marker")

library(colorspace)
join3 %>% ggplot(aes(x = TFmarker2, fill = Cell_type_class)) +
  geom_bar(width = 0.5) +
  scale_fill_manual(values = c())
