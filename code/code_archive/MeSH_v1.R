df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/MeSH/gene2pubmed_human.tsv", col_names = F)
colnames(df) <- c("tax_id", "Entrez_ID", "PMID")
head(df)
dim(df)
length(unique(df$Entrez_ID))
length(unique(df$PMID))

# ENTREZIDとEnsemblIDやSymbolと紐づけ
# 「Gene ID と MeSH term の対応表」は AnnotationHub の一部として公開されている
BiocManager::install("AnnotationHub")
library("AnnotationHub")
ah <- AnnotationHub()
head(mcols(ah)) # # mcolsでデータフレーム状に一覧
length(ah) # how many resources? -> 67944

# PMIDとMeSH termの対応表
qr <- query(ah, c("PubMedDb"))
qr$title
qr$tags
table(qr$description)
qr2 <- query(ah, c("PubMed", "PMID-DescriptorID-DescriptorTerm","v002"))
qr2 #目的の項目が入っているAH番号(=データベース)を見つける
hoge <- ah[["AH97924"]] # データベースをダウンロードする

head(hoge, n=20)
length(unique(hoge$pmid))

# GENE-PMID-MESH ID -----
tgt_PMID <- df$PMID %>% unique()
hoge2 <- hoge %>% filter(pmid %in% tgt_PMID)
colnames(hoge2) <- c("PMID", "descriptorid", "descriptorterm")
dim(hoge2)
length(unique(hoge2$PMID))
write_tsv(hoge2, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/MeSH/PubMedDb_AH97924.tsv")

df2 <- df %>% left_join(hoge2, by = "PMID", multiple = "all")
dim(df2)
length(unique(df2$PMID))
write_tsv(df2, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/MeSH/Entrez_PMID_AH97924_human.tsv")

# check ----
df2 <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/MeSH/Entrez_PMID_AH97924_human.tsv")
head(df2)
dim(df2)

print("Entrez_ID")
length(unique(df2$Entrez_ID))
print("PMID")
length(unique(df2$PMID))
print("descriptorid")
length(unique(df2$descriptorid))
print("descriptorterm")
length(unique(df2$descriptorterm))

# Entrex -> symbol ----
if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

library(biomaRt)
db <- useMart("ensembl")
listDatasets(db)
hd <- useDataset("hsapiens_gene_ensembl", mart = db)
View(listFilters(hd)) #属性チェック

entrez_list <- df2$Entrez_ID %>% unique()
length(entrez_list)
res <- getBM(attributes = c('entrezgene_id','hgnc_symbol'), filters = 'entrezgene_id', values = entrez_list, mart = hd, useCache = FALSE)
head(res)
dim
View(res)
colnames(res) <- c("Entrez_ID", "Gene_symbol")

df3 <- df2 %>% left_join(res, "Entrez_ID", multiple = "all")
dim(df3)
write_tsv(df3, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/MeSH/Entrez_PMID_AH97924_symbol_human.tsv")


# statstics ----
df3 <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/MeSH/Entrez_PMID_AH97924_symbol_human.tsv")
length(unique(df3$Gene_symbol))

df3 %>% group_by(Gene_symbol)  %>% summarise(num = n()) %>%
  ggplot(aes(x = reorder(Gene_symbol, -num), y = num)) +
  geom_point() +
  xlab("Gene symbol") +
  ylab("Number of publications") +
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=0,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 0.7
  )
  

Chip_df <- readRDS("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/TF_CTC_sample_num_new.rds")
TF_list <- rownames(Chip_df)
tmp  <- df3 %>% filter(Gene_symbol %in% TF_list) %>%
  drop_na(Gene_symbol) %>% select(Gene_symbol, PMID) %>% distinct() %>%
  group_by(Gene_symbol)  %>% summarise(num = n()) 

library(ggrepel)
tmp %>%
  ggplot(aes(x = reorder(Gene_symbol, -num), y = num, label = Gene_symbol)) +
  geom_point() +
  geom_text_repel() +
  xlab("TFs") +
  ylab("Number of publications") +
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=0,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 0.7
  )
summary(tmp$num)
dim(tmp)
View(tmp)

count <- 0
tmp <- tmp %>% arrange(desc(num))
count_list <- c()
ratio_list <- c()
for (i in 1:nrow(tmp)) {
  tgt_tf <- tmp[i,1]
  count <- count + tmp[i,2]
  if(i == 1){
    count_list <- count
    ratio_list <- count / sum(tmp$num)
  }else{
    count_list <- c(count_list, count) 
    ratio_list <- c(ratio_list, count/sum(tmp$num)) 
  }
  
}

tmp2 <- tmp %>% mutate(count = as.numeric(count_list), ratio = as.numeric(ratio_list))
View(tmp2)

tmp2 %>% ggplot(aes(reorder(Gene_symbol, -num), y = ratio)) +
  geom_point(size = 1) +
  geom_hline(yintercept=0.25,linetype=2,alpha=0.7,size=0.5,color='black') +
  geom_hline(yintercept=0.5,linetype=2,alpha=0.7,size=0.5,color='black') +
  geom_hline(yintercept=0.75,linetype=2,alpha=0.7,size=0.5,color='black') +
  xlab("TF") +
  ylab("Ratio") +
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=1,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 0.7
  )
