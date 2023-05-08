# ファイル名を全て.csvに直す -----
rm(list=ls())
file_list <- list.files(path = paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix"))
for (tgt in file_list) {
  if(str_detect(tgt, pattern = ".txt")==TRUE){
    tgt_new <- gsub(".txt", "", tgt)
    tgt_new2 <- paste0(tgt_new, ".csv")
    system(paste0("mv /Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix/", tgt, " /Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix/", tgt_new2))
  }
}

# bind all dataset ----
count <- 0
count_no_row_sample <- 0
for (tgt in file_list) {
  if(file.exists(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix/", tgt)) == TRUE){
    target_df <- c()
    target_df <- suppressMessages(read_csv(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix/", tgt), col_names = T))
    print(dim(target_df))
    if(nrow(target_df)!=0){
      count <- count + as.numeric(nrow(target_df))
    }else{ #if(is.null(nrow(target_df))==FALSE)
      count_no_row_sample  <- count_no_row_sample  + 1
    }
  } #if(file.exist)
} #for

print(count_no_row_sample)
df_binded <- data.frame(matrix(nrow = count, ncol = 15))

num <- 0 
for (tgt in file_list) {
  if(file.exists(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix/", tgt)) == TRUE){
    target_df <- suppressMessages(read_csv(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix/", tgt), col_names = T))
    target_df <- target_df %>% mutate(file_name = tgt)
    if(is.null(nrow(target_df)) == FALSE & nrow(target_df) != 0){
      df_binded[(num+1):(num+nrow(target_df)),] <- target_df
      num <- num + nrow(target_df)
    }
  } #if(file.exist)
} #for
colnames(df_binded) <- colnames(target_df)
dim(df_binded)
head(df_binded)
df_binded2 <- df_binded[,1:13] %>% as_tibble()
colnames(df_binded2) <- c("Target_Gene", "TF", "Mean_Expr_ctrl", "MeanExpr_Treat", "FC", "logFC", 
                      "Rank", "pval", "corrected_p", "Promoter_TF", "SE_TF", "TE_TF", "file_name")
#write_tsv(df_binded2, paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix_md/DEG_all.tsv"))


# statistics ----
df_all <- read_tsv(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix_md/DEG_all.tsv"))
df_all

# humanとmouseにわける
mouse <- df_all %>% filter(str_detect(TF, regex("[A-Z]{1}[a-z]{1}")) | str_detect(TF, regex("[a-z]{1}[0-9]{1}")) | str_detect(TF, regex("[A-Z]{1}[0-9]{1}[a-z]{1}")))
TF_mouse <- mouse$TF %>% unique()
human <- df_all%>% filter(!TF %in% TF_mouse)
#write_tsv(mouse, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix_md/DEG_all_mouse.tsv")
#write_tsv(human, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix_md/DEG_all_human.tsv")

human <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix_md/DEG_all_human.tsv")
mouse <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix_md/DEG_all_mouse.tsv")

length(unique(df_all$file_name))
length(unique(df_all$TF))
unique(df_all$TF)

length(unique(human$file_name))
length(unique(human$TF))
unique(human$TF)

length(unique(mouse$file_name))
length(unique(mouse$TF))
unique(mouse$TF)

# Add annotation ----
statistics <- read_csv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/KnockTF_Statistics.csv")
colnames(statistics) <- c("ID", "TF", "Species", "Knock_method", "Biosample_name", "Profile_ID", "Platform" , "2FC_up", "2FC_down", "1.5FC_up", "1.5FC_down")
View(statistics)

hu_file_name <- human$file_name 
hu_ID <- gsub(".csv", "", hu_file_name)
human2 <- human %>% mutate(ID = hu_ID)
stat_hu <- statistics %>% filter(ID %in% hu_ID) 
stat_hu2 <- human2 %>% left_join(stat_hu, by = "ID")
write_tsv(stat_hu2, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix_md/DEG_all_human_annotated.tsv")

mo_file_name <- mouse$file_name 
mo_ID <- gsub(".csv", "", mo_file_name)
mouse2 <- mouse %>% mutate(ID = mo_ID)
stat_mo <- statistics %>% filter(ID %in% mo_ID)
stat_mo2 <- mouse2 %>% left_join(stat_mo, by = "ID")
write_tsv(stat_mo2, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix_md/DEG_all_mouse_annotated.tsv")

# statistics_Biosamplename
stat_hu2 <- read_tsv( "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix_md/DEG_all_human_annotated.tsv")
stat_mo2 <- read_tsv( "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/KnockTF_v2/DEG_matrix_md/DEG_all_mouse_annotated.tsv")

length(unique(stat_hu$Biosample_name))
length(unique(stat_mo$Biosample_name))

summary(as.numeric(stat_hu2$Mean_Expr_ctrl))
hist(as.numeric(stat_hu2$Mean_Expr_ctrl), breaks = 100)
hist(log(stat_hu2$Mean_Expr_ctrl, base = 10), breaks = 100)
summary(stat_hu2$MeanExpr_Treat)
hist(log(stat_hu2$MeanExpr_Treat, base = 10), breaks = 100)
summary(stat_hu2$FC)
hist(log(stat_hu2$FC, base = 10), breaks = 100)
summary(stat_hu2$logFC)
summary(stat_hu2$pval)
summary(stat_hu2$corrected_p)
