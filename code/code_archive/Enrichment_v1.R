library(tidyverse)
library(biomaRt)
db <- useMart("ensembl")
listDatasets(db)
hd <- useDataset("hsapiens_gene_ensembl", mart = db)
View(listFilters(hd)) #属性チェック

df2 <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/MeSH/Entrez_PMID_AH97924_human.tsv")
entrez_list <- df2$Entrez_ID %>% unique()
length(entrez_list)
res <- getBM(attributes = c('entrezgene_id','hgnc_symbol'), filters = 'entrezgene_id', values = entrez_list, mart = hd, useCache = FALSE)
head(res)
dim(res)
View(res)
colnames(res) <- c("Entrez_ID", "Gene_symbol")

write_tsv(res, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/MeSH/human_gene_symbol.tsv")

res2 <- res$Gene_symbol %>% as.character()
length(res2)
write(res2, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/MeSH/human_gene_symbol2.txt")

# ChIP-Atlas all genesでEnrichmentを行った結果をSRX番号ごとにplot
result_enrich <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/all_human_gene_enruchment_result_TFs.tsv", col_names = F)
colnames(result_enrich) <- c("ID", "Antigen class", "TF", "Cell type class", "Cell type", "Number of peaks", "Overlaps/Mydata", "Overlaps/Control",
                             "log_pval", "log_qval", "Fold_enrichment", "FE>1?")
View(result_enrich[1:10,])


result_enrich2 <- result_enrich %>% arrange(ID)
result_enrich_SRX <- result_enrich2 %>% filter(str_detect(ID, "SRX")) 
SRX_num <- gsub("SRX", "", result_enrich_SRX$ID) %>% as.numeric()
result_enrich_SRX2 <- result_enrich_SRX %>% mutate(ID_num = SRX_num) %>% arrange(ID_num) 
View(result_enrich_SRX2)
length(unique(result_enrich$ID))
length(unique(result_enrich_SRX2$ID))

summary(result_enrich_SRX2$ID_num)
hist(result_enrich_SRX2$ID_num, breaks = 100)

res_df <- c()
for (i in 1:180) {
  tgt_end <- i*100000
  tgt_start <- tgt_end-100000+1
  #print(paste(tgt_start, tgt_end))
  tgt_df <- result_enrich_SRX2 %>% filter(ID_num >= tgt_start & ID_num < tgt_end) 
  tgt_num <- length(tgt_df$ID)  
  
  if(i == 1){
    tgt_sum <- tgt_num
    res_row <- c(paste0(tgt_start, "_", tgt_end), tgt_num, tgt_sum) %>% t()
    res_df <- res_row
  }else{
    tgt_sum <- sum(as.numeric(res_df[,2])) + tgt_num
    res_row <- c(paste0(tgt_start, "_", tgt_end), tgt_num, tgt_sum)
    res_df <- rbind(res_df, res_row)
  }
}

res_df <- res_df %>% as_tibble()
res_df$V1 <- as.character(res_df$V1)
res_df$V2 <- as.numeric(res_df$V2)
res_df$V3 <- as.numeric(res_df$V3)
colnames(res_df) <- c("ID", "Number_of_SRX", "culumtive_num")
res_df2 <- res_df %>% separate(ID, c("start", "end"), sep = "_")
res_df2$start <- as.numeric(res_df2$start)
res_df2$end <- as.numeric(res_df2$end)
res_df2 %>% ggplot(aes(x = start, y = culumtive_num)) +
  geom_point()+
  xlab("ID number")+
  ylab("Number of SRX")+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        legend.position = "none",
        panel.grid.major = element_line(colour="gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 0.7
  )
write_tsv(res_df2, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/res2.tsv")




#SRX_list <- result_enric2$ID
#write(SRX_list, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/SRX_list.txt")
