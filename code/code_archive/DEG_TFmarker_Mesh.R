# DEGの数, TFmarker, MeSHを組み合わせる
# DEG---
df_DEGs <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/tib_2FC_DEGs_ctc_tissue_unm.tsv")
head(df_DEGs)

# TFmarker -----
unmeasured_df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/TFmarker/unmeasured_tf_tissue_TFmarker_selected.tsv")
unmeasured_df <- unmeasured_df %>% drop_na(ChIPAtlas_ctc)
dim(unmeasured_df)


unmeasured_df2 <- unmeasured_df %>% unite("combi", c(Gene_Name, ChIPAtlas_ctc), sep = "_")
df_DEGs2 <- df_DEGs %>% unite("combi", c(TF, Cell_type_class), sep = "_")
df_join <- unmeasured_df %>% left_join(df_DEGs, by = TF)


# MeSH -----