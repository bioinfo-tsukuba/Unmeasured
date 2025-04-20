# Supplementary table
S2_join_TF_3 <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/ignored_paper/S2_join_902TF_num.tsv")
S2_join_TF_3$pub_num <- (S2_join_TF_3$pub_num - mean(S2_join_TF_3$pub_num)) / sd(S2_join_TF_3$pub_num)
S2_join_TF_3$ChIPseq_num <- (S2_join_TF_3$ChIPseq_num - mean(S2_join_TF_3$ChIPseq_num)) / sd(S2_join_TF_3$ChIPseq_num)
dim(S2_join_TF_3)
write_tsv(S2_join_TF_3, "/Users/saeko/Unmeasured/data/supple_table/xgboost_table.tsv")

# feature対応表
exp_var_table <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/ignored_paper/exp_var_number_table.tsv")
write_tsv(exp_var_table, "/Users/saeko/Unmeasured/data/supple_table/TF_feature.tsv")

# Unmeasured TF, Cell type class list (shinyのdata_prep_v1.Rで作成)
tib_unmeasured <- read_tsv("/Users/saeko/Unmeasured/shiny/Unmeasured_shiny_v1/data/tib_unmeasured.tsv")
write_tsv(tib_unmeasured, "/Users/saeko/Unmeasured/data/supple_table/tib_unmeasured.tsv")
