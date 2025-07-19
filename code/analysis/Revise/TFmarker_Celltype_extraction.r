# PMIDに関連するCell lineをとってくる
TF_marker <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/TFmarker/All_TFmarkers_annotated.tsv")
length(unique(TF_marker$Gene_Name))
colnames(TF_marker) <- c("PMID", "TF", "Gene_type", "Cell_Name", "Cell_Type", "Cell_type_class")

marker_TF_list <- TF_marker %>% filter(Gene_type == "TFMarker") %>% .$TF %>% unique()
length(marker_TF_list)

df_chip <- read_tsv("/Users/saeko/Unmeasured/data/SRX_date_20231004_unmeasured_added.tsv")

# PMIDのAbstract抽出
library(rentrez)
results <- data.frame(PMID = character(), Contains_CellLine = logical(), Abstract = character())
target_cell_lines_all <- df_chip$Cell_type %>% unique()
regex_cell <- paste0("\\b(", paste(target_cell_lines, collapse = "|"), ")\\b")

# 抽出と判定
for (i in 1:nrow(TF_marker)) {
#for (i in 1:10) {
  tgt_row <- TF_marker[i,]
  pmid <- tgt_row$PMID
  
  abstract <- tryCatch({
    entrez_summary(db = "pubmed", id = pmid)$title
    rec <- entrez_fetch(db = "pubmed", id = pmid, rettype = "abstract", parsed = FALSE)
    rec
  }, error = function(e) NA)
  
  if (is.na(abstract)) next
  
  # Cell line名を抽出
  matched_cells <- str_extract_all(abstract, regex_cell, simplify = FALSE)[[1]]
  matched_cells <- unique(toupper(trimws(matched_cells)))  # 重複削除・大文字化
  
  if (length(matched_cells) > 0) {
    results <- rbind(
      results,
      data.frame(
        PMID = pmid,
        Matched_Cell_Lines = paste(matched_cells, collapse = ", "),
        stringsAsFactors = FALSE
      )
    )
  }
}

write_tsv(results, "/Users/saeko/Unmeasured/data/TFmarker/PMID_Cell_type.tsv")

# TF markerのdataをjoin
tgt_PMIDs <- results$PMID %>% unique()
results_joined <- results %>% filter(PMID %in% tgt_PMIDs)%>% left_join(TF_marker, by = "PMID") #%>% 
  #select(PMID, Matched_Cell_Lines, Cell_Name, Cell_Type, Cell_type_class) %>% distinct()
View(results_joined)
write_tsv(results_joined, "/Users/saeko/Unmeasured/data/TFmarker/PMID_Cell_type_joined.tsv")
