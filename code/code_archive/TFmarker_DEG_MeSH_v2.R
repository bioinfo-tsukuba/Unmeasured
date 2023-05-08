rm(list=ls())
install.packages("ggside")
library(ggside)
library(tidyverse)
library(ggrepel)

# 351TF, human, Knock-TF -----
######## ここはやらなくてok #########################################################################################################################
stat_hu4 <- read_tsv( "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix_md/DEG_all_human_annotated_tissue_ctc.tsv")
dim(stat_hu4)
length(unique(stat_hu4$Profile_ID))
length(unique(stat_hu4$TF.x))
length(unique(stat_hu4$Biosample_name))
View(stat_hu4[1:10,])

colnames(stat_hu4)
stat_hu5 <- stat_hu4 %>% select(Target_Gene, TF.x, corrected_p,Cell_type_class) %>% distinct()
dim(stat_hu5)

stat_hu_DEG <- stat_hu5 %>% filter(corrected_p < 0.05) %>% group_by(TF.x, Cell_type_class)  %>% summarise(DEG_num = n())
dim(stat_hu_DEG)
head(stat_hu_DEG)
summary(stat_hu_DEG$DEG_num)
colnames(stat_hu_DEG) <- c("TF", "Cell_type_class", "DEG_num_thre_p")
write_tsv(stat_hu_DEG, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix_md/DEG_all_human_threshold_p.tsv")
####################################################################################################################################################

stat_hu_DEG <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix_md/DEG_all_human_threshold_p.tsv")
TF_list_Knock <- stat_hu_DEG$TF %>% unique()
length(TF_list_Knock)

# publication数 -----
pub1 <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/MeSH/Entrez_PMID_AH97924_symbol_human.tsv")
pub2 <- pub1 %>% filter(Gene_symbol %in% TF_list_Knock) %>% select(-descriptorid, -descriptorterm) %>%
  distinct() %>% group_by(Gene_symbol) %>% summarise(pub_count = n())
colnames(pub2) <- c("TF", "pub_count")
df_join1 <- stat_hu_DEG %>% left_join(pub2, by = "TF")
View(df_join1)
dim(df_join1)
join_key <- df_join1 %>% unite("key", c(TF, Cell_type_class), sep = "_") %>% .$key
df_join1 <- df_join1 %>% mutate(key = join_key)

# Measured/Unmeasured -----
Chip_df <- readRDS("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/TF_CTC_sample_num_new.rds")
tmp_TF_list <- rownames(Chip_df)
Chip_df2 <- Chip_df %>% as_tibble() %>% mutate(TF = tmp_TF_list) %>% pivot_longer(-TF, names_to = "Cell_type_class", values_to = "count_ChIPseq")
join_key2 <- Chip_df2 %>%  unite("key", c(TF, Cell_type_class), sep = "_") %>% .$key
Chip_df3 <- Chip_df2 %>% mutate(key = join_key2) %>% select(key, `count_ChIPseq`)
dim(Chip_df3)

# KnockTFに載っていた351TFのうち、ChIP-AtlasのTFのみに絞ると、152種類になった
ChIP_TF <- Chip_df2$TF %>% unique()
df_join2 <- df_join1 %>% filter(TF %in% ChIP_TF) %>% left_join(Chip_df3, by = "key")
dim(df_join2)
length(unique(df_join2$TF))
View(df_join2)

# TF marker ---
######## ここはやらなくてok #########################################################################################################################
TF_marker_all <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/TFmarker/All_TFmarkers.txt")
taiou <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/TFmarker/TFmarker_tissue_ChIP-Atlas_Ctc_taiou.tsv")
colnames(taiou) <- c("Tissue_Type", "ChIPAtlas_ctc", "description")
taiou <- taiou %>% select(-description) %>% drop_na(ChIPAtlas_ctc)
colnames(TF_marker_all) <- c("PMID", "Gene_Name","Gene_Type","Cell_Name","Cell_Type","Tissue_Type","Experiment_Type","Experimental_Method","Title","Description_of_Gene","Interacting_Gene","CellOntologyID")
df_all2 <- TF_marker_all %>% left_join(taiou, by = "Tissue_Type") %>% select(PMID, Gene_Name, Gene_Type, Cell_Name, Cell_Type,Tissue_Type)
write_tsv(df_all2, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/TFmarker/All_TFmarkers_annotated.tsv")
####################################################################################################################################################

TF_marker <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/TFmarker/All_TFmarkers_annotated.tsv")
length(unique(TF_marker$Gene_Name))
colnames(TF_marker) <- c("PMID", "TF", "Gene_type", "Cell_Name", "Cell_Type", "Cell_type_class")

TF_list_Knock_ChIPAtlas <- df_join2$TF %>% unique()
join_key3 <- TF_marker %>% unite("key", c(TF, Cell_type_class)) %>% .$key
TF_marker2 <- TF_marker %>% mutate(key = join_key3) %>% select(TF, Cell_type_class, key, Gene_type) %>% filter(TF %in% TF_list_Knock_ChIPAtlas) %>% 
  select(key, Gene_type) %>% distinct()

df_join3 <- df_join2 %>% left_join(TF_marker2, by = "key", multiple = "all")
View(df_join3)
length(unique(df_join3$TF))
write_tsv(df_join3, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/DEG_TFmarker_pub/summary_allstats.tsv")

# plot -----_
df_join3 <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/DEG_TFmarker_pub/summary_allstats.tsv")
df_join3 <- df_join3 %>% drop_na(TF, Cell_type_class)
length(unique(df_join3$TF))
length(unique(df_join3$Cell_type_class))
length(unique(df_join3$key))

cor.test(df_join3$pub_count, df_join3$DEG_num_thre_p, method = "spearman")
cor.test(df_join3$pub_count, df_join3$DEG_num_thre_p, method = "pearson")

color_list <- c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG"), "gray")
df_join3 %>% select(-Gene_type) %>% distinct() %>%
  ggplot(aes(x = log(pub_count, base = 10), y = DEG_num_thre_p, 
                        #color = Cell_type_class,
                        label = key
                        )) +
  geom_point() +
  geom_text_repel() +
  xlab("log10(Number of publication)")+
  ylab("Number of DEGs (corrected p value < 0.05)") +
  #scale_color_manual(values = color_list)+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        #legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=12,face="bold"),
        axis.text.y =element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 0.9
  )

df_join3 %>% select(-Gene_type) %>% distinct() %>%
  ggplot(aes(x = log(pub_count, base = 10), y = DEG_num_thre_p, 
             color = Cell_type_class,
             label = TF
  )) +
  geom_point() +
  geom_text_repel() +
  xlab("log10(Number of publication)")+
  ylab("Number of DEGs (corrected p value < 0.05)") +
  scale_color_manual(values = color_list)+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        #legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=12,face="bold"),
        axis.text.y =element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 0.9
  )

df_join3 %>% select(-Gene_type) %>% distinct() %>%
  ggplot(aes(x = log(pub_count, base = 10), 
                        y = log(DEG_num_thre_p, base = 10), 
                        color = Cell_type_class,
                        label = TF
                        )) +
  geom_point() +
  geom_text_repel() +
  xlab("log10(Number of publication)")+
  ylab("log10(Number of DEGs (corrected p value < 0.05))") +
  scale_color_manual(values = color_list)+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        #legend.position = "none",
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_line(colour = "gray"),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=8,face="bold"),
        axis.text.x =element_text(size=7,face="bold"),
        axis.text.y =element_text(size=8,face="bold"),
        axis.title=element_text(size=10,face="bold"),
        ggside.panel.scale = 0.4,
        aspect.ratio = 0.7
  ) +
  geom_xsideboxplot(aes(y = Cell_type_class), orientation = "y") + 
  scale_xsidey_discrete(guide = guide_axis(angle = 0))+
  geom_ysideboxplot(aes(x = Cell_type_class), orientation = "x") +
  scale_ysidex_discrete(guide = guide_axis(angle = 45))
  
df_join3 %>% select(-Gene_type) %>% distinct() %>%
  ggplot(aes(x = log(pub_count, base = 10), 
                        y = log(DEG_num_thre_p, base = 10), 
                        label = key
)) +
  geom_point() +
  geom_text_repel() +
  xlab("log10(Number of publication)")+
  ylab("log10(Number of DEGs (corrected p value < 0.05))") +
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        #legend.position = "none",
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_line(colour = "gray"),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=8,face="bold"),
        axis.text.x =element_text(size=7,face="bold"),
        axis.text.y =element_text(size=8,face="bold"),
        axis.title=element_text(size=10,face="bold"),
        ggside.panel.scale = 0.2,
        aspect.ratio = 0.7
  ) +
  geom_xsidedensity(aes(y = stat(density)), position = "stack") +
  scale_xsidey_discrete(guide = guide_axis(angle = 0))+
  geom_ysidedensity(aes(x = stat(density)), position = "stack") +
  scale_ysidex_discrete(guide = guide_axis(angle = 45))

# 細胞型ごとに別plotを並べる
df_join3 %>% select(-Gene_type) %>% distinct() %>%
  ggplot(aes(x = log(pub_count, base = 10), 
                        y = log(DEG_num_thre_p, base = 10), 
                        color = Cell_type_class,
                        label = TF
)) +
  geom_point() +
  geom_text_repel() +
  xlab("log10(Number of publication)")+
  ylab("log10(Number of DEGs (corrected p value < 0.05))") +
  scale_color_manual(values = color_list)+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        legend.position = "none",
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_line(colour = "gray"),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=8,face="bold"),
        axis.text.x =element_text(size=7,face="bold"),
        axis.text.y =element_text(size=8,face="bold"),
        axis.title=element_text(size=10,face="bold"),
        ggside.panel.scale = 0.2,
        aspect.ratio = 1
  ) +
  facet_wrap(. ~ Cell_type_class, nrow = 4)

# 点の色を、TF markerごとに分ける。NAならgrayにする
df_join4 <- df_join3 %>% mutate(label_marker = ifelse(is.na(Gene_type) == FALSE, "Marker", "No marker")) %>%
  mutate(label_chip_measure = ifelse(count_ChIPseq == 0, "Unmeasured", "Measured")) 

df_join4 %>% ggplot(aes(x = log(pub_count, base = 10), y = DEG_num_thre_p, 
                        color = Gene_type,
                        label = TF
)) +
  geom_point() +
  geom_text_repel() +
  xlab("log10(Number of publication)")+
  ylab("Number of DEGs (corrected p value < 0.05)") +
  scale_color_manual(values = c( "#5E4FA2", "#2B83BA","green4", "#FDAE61", "#D7191C", "#D3D3D3"))+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        #legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=12,face="bold"),
        axis.text.y =element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        ggside.panel.scale = 0.35,
        aspect.ratio = 0.9
  )+
  geom_xsideviolin(aes(y = label_marker), orientation = "y") + 
  scale_xsidey_discrete(guide = guide_axis(angle = 0))+
  geom_ysideviolin(aes(x = label_marker), orientation = "x") +
  scale_ysidex_discrete(guide = guide_axis(angle = 45))


df_join4 %>% select(-Gene_type) %>% distinct() %>%
  ggplot(aes(x = log(pub_count, base = 10), y = DEG_num_thre_p, 
                        color = label_marker,
                        label = key
)) +
  geom_point() +
  geom_text_repel() +
  xlab("log10(Number of publication)")+
  ylab("Number of DEGs (corrected p value < 0.05)") +
  scale_color_manual(values = c("#D7191C", "gray"))+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        #legend.position = "none",
        panel.grid.major = element_line(color = "gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=10,face="bold"),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        ggside.panel.scale = 0.2,
        aspect.ratio = 0.9
  )+
  geom_xsideviolin(aes(y = label_marker), orientation = "y") + 
  geom_xsideboxplot(aes(y = label_marker), orientation = "y") + 
  scale_xsidey_discrete(guide = guide_axis(angle = 0))+
  geom_ysideviolin(aes(x = label_marker), orientation = "x") +
  geom_ysideboxplot(aes(x = label_marker), orientation = "x") +
  scale_ysidex_discrete(guide = guide_axis(angle = 45))

marker_DEG <- df_join4 %>% filter(label_marker == "Marker") %>% .$DEG_num_thre_p
no_marker_DEG <- df_join4 %>% filter(label_marker == "No marker") %>% .$DEG_num_thre_p
wilcox.test(marker_DEG, no_marker_DEG, paired = FALSE)

marker_pub <- df_join4 %>% filter(label_marker == "Marker") %>% .$pub_count
no_marker_pub <- df_join4 %>% filter(label_marker == "No marker") %>% .$pub_count
wilcox.test(marker_pub, no_marker_pub, paired = FALSE)

# 点の色をChIP-seqサンプルごとに分ける
df_join4 %>% select(-Gene_type) %>% distinct() %>%
  ggplot(aes(x = log(pub_count, base = 10), y = DEG_num_thre_p, 
                        color = log(count_ChIPseq, base = 10),
                        label = key
)) +
  geom_point() +
  geom_text_repel() +
  xlab("log10(Number of publication)")+
  ylab("Number of DEGs (corrected p value < 0.05)") +
  scale_colour_gradientn(colours = c("blue", "red"))+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        #legend.position = "none",
        panel.grid.major = element_line(color = "gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=10,face="bold"),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 0.9
  )

df_join4 %>% select(-Gene_type) %>% distinct() %>%
  ggplot(aes(x = log(pub_count, base = 10), y = DEG_num_thre_p, 
                        color = label_chip_measure,
                        label = key
)) +
  geom_point() +
  geom_text_repel() +
  xlab("log10(Number of publication)")+
  ylab("Number of DEGs (corrected p value < 0.05)") +
  scale_color_manual(values = c("#D7191C", "blue"))+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        #legend.position = "none",
        panel.grid.major = element_line(color = "gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=10,face="bold"),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        ggside.panel.scale = 0.2,
        aspect.ratio = 0.9
  )+
  geom_xsideviolin(aes(y = label_chip_measure), orientation = "y") + 
  geom_xsideboxplot(aes(y = label_chip_measure), orientation = "y") + 
  scale_xsidey_discrete(guide = guide_axis(angle = 0))+
  geom_ysideviolin(aes(x = label_chip_measure), orientation = "x") +
  geom_ysideboxplot(aes(x = label_chip_measure), orientation = "x") +
  scale_ysidex_discrete(guide = guide_axis(angle = 45))

df_join4 %>% select(-Gene_type) %>% distinct() %>%
  ggplot(aes(x = log(pub_count, base = 10), y = DEG_num_thre_p, 
             color = label_chip_measure,
             label = TF
  )) +
  geom_point() +
  geom_text_repel() +
  xlab("log10(Number of publication)")+
  ylab("Number of DEGs (corrected p value < 0.05)") +
  scale_color_manual(values = c("#D7191C", "blue"))+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        #legend.position = "none",
        panel.grid.major = element_line(color = "gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=10,face="bold"),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        ggside.panel.scale = 0.2,
        aspect.ratio = 0.9
  )

# 検定
measured_DEG <- df_join4 %>% filter(label_chip_measure == "Measured") %>% .$DEG_num_thre_p
unmeasured_DEG <- df_join4 %>% filter(label_chip_measure == "Unmeasured") %>% .$DEG_num_thre_p
wilcox.test(measured_DEG, unmeasured_DEG, paired = FALSE)

measured_pub <- df_join4 %>% filter(label_chip_measure == "Measured") %>% .$pub_count
unmeasured_pub <- df_join4 %>% filter(label_chip_measure == "Unmeasured") %>% .$pub_count
wilcox.test(measured_pub, unmeasured_pub, paired = FALSE)


# TF markerのアノテーションがあったものだけでChIP-seqサンプル数をプロットする
df_join5 <- df_join4 %>% drop_na(Gene_type)

df_join5 %>% select(-Gene_type) %>% distinct() %>%
  ggplot(aes(x = log(pub_count, base = 10), y = DEG_num_thre_p, 
                          color = label_chip_measure,
                          label = key
)) +
  geom_point() +
  geom_text_repel() +
  xlab("log10(Number of publication)")+
  ylab("Number of DEGs (corrected p value < 0.05)") +
  scale_colour_manual(values = c( "red", "blue"))+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        #legend.position = "none",
        panel.grid.major = element_line(color = "gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=10,face="bold"),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        ggside.panel.scale = 0.2,
        aspect.ratio = 0.9
  )

df_join5 %>% select(-Gene_type) %>% distinct() %>%
  ggplot(aes(x = log(pub_count, base = 10), y = DEG_num_thre_p, 
             color = log(count_ChIPseq),
             label = key
  )) +
  geom_point() +
  geom_text_repel() +
  xlab("log10(Number of publication)")+
  ylab("Number of DEGs (corrected p value < 0.05)") +
  scale_colour_gradientn(colours = c("blue", "red"))+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        #legend.position = "none",
        panel.grid.major = element_line(color = "gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=10,face="bold"),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        ggside.panel.scale = 0.2,
        aspect.ratio = 0.9
  )
