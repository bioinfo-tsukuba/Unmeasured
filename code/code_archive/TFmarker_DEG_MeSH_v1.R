# TF marker --------
# TF-markerにはあるのにUnmeasuredだったもの 
unmeasured_df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/TFmarker/unmeasured_tf_tissue_TFmarker_selected.tsv")
unmeasured_df <- unmeasured_df %>% drop_na(ChIPAtlas_ctc)
dim(unmeasured_df)
head(unmeasured_df)
unmeasured_df2 <- unmeasured_df %>% select(Gene_Name, Gene_Type, ChIPAtlas_ctc, Cell_Name) %>% distinct()
colnames(unmeasured_df2) <- c("tf", "Gene_Type", "ctc", "Cell_Name")
unmeasured_df3 <- unmeasured_df2 %>% mutate(label_tfmarker = "unmeasured") 
key <- unmeasured_df3 %>% unite("key", c("tf", "ctc"), sep = "_") %>% .$key
unmeasured_df4 <- unmeasured_df3 %>% mutate(key = key)

# DEGs --------
df_DEGs <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/tib_2FC_DEGs_ctc_tissue_unm.tsv")
head(df_DEGs)
colnames(df_DEGs) <- c("tf", "ctc", "count_2FC_DEGs", "label")
key2 <- df_DEGs %>% unite("key", c("tf", "ctc"), sep = "_") %>% .$key
df_DEGs2 <- df_DEGs %>% mutate(key = key2) %>% select(-tf, -ctc)

df_join1 <- unmeasured_df4 %>% left_join(df_DEGs2, by = "key", multiple = "all")
View(df_join1)
df_join2 <- df_join1 %>% drop_na(count_2FC_DEGs) %>%distinct()
dim(df_join2)
View(df_join2)

left_tf <- df_join2$tf %>% unique()
left_ctc <- df_join2$ctc %>% unique()


# publication ------
pub1 <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/MeSH/Entrez_PMID_AH97924_symbol_human.tsv")
pub2 <- pub1 %>% filter(Gene_symbol %in% left_tf) %>% distinct() %>% group_by(Gene_symbol) %>% summarise(pub_count = n())
colnames(pub2) <- c("tf", "pub_count")
df_join3 <- df_join2 %>% left_join(pub2, by = "tf")
View(df_join3)

df_join3 %>% select(tf, count_2FC_DEGs, pub_count, key, Gene_Type) %>% distinct() %>%
  ggplot(aes(x = count_2FC_DEGs, y = pub_count, label = key, color = Gene_Type)) +
  geom_point() +
  geom_text_repel() +
  xlab("Number of 2FC-DEGs")+
  ylab("Number of publication")+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        #legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=0,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 1
  )
