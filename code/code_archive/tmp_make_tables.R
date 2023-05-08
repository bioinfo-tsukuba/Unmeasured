Chip_df <- readRDS("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/TF_CTC_sample_num_new.rds")
View(Chip_df )

dim(Chip_df)
TFs <- rownames(Chip_df)
Chip_df2 <- Chip_df %>% as_tibble()  %>% mutate(TF = TFs) %>% pivot_longer(cols = -TF, names_to = "Cell_type_class", values_to = "count")
write_tsv(Chip_df2, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/summary_meta_ChIP_Atlas/summary_ChIP_count.tsv")

library(tidyverse)
df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/experimentList.tab", col_names = F)
head(df)
df2 <- df %>%  select(X1, X2, X3, X4, X5, X6)
head(df2)
colnames(df2) <- c("ID", "Genome", "Antigen_class", "Antigen", "Cell_type_class", "Cell_type")
head(df2)
dim(df2)
anno_df <- df2 %>% filter(Genome == "hg38" & Antigen_class == "TFs and others")
dim(anno_df)
length(unique(anno_df$ID))
TF_summary <- anno_df %>% group_by(Antigen) %>% summarise(n = n())

write_tsv(anno_df, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/summary_meta_ChIP_Atlas/anno_df.tsv")


# publication - release year -------
pub1 <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/MeSH/Entrez_PMID_AH97924_symbol_human.tsv")
