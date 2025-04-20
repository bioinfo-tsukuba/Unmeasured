library(tidyverse)
library(ggrepel)
library(ggside)

# /Users/saeko/Documents/MOCCS/important_chipseq_prediction/code/TFmarker_DEG_MeSH_v3_20231125.R より
df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/DEG_TFmarker_pub/summary_allstats_v3.tsv")
df <- df %>% drop_na(TF, Cell_type_class)
df2 <- df %>% 
  mutate(label_marker = ifelse(is.na(Gene_type) == FALSE, "Marker", "No marker")) %>%
  mutate(label_chip_measure = ifelse(count_ChIPseq == 0, "Unmeasured", "Measured")) %>%
  select(-Gene_type) %>% distinct() 
tgt_TFs <- df2$TF %>% unique()

# RefEx
mt_RNA <- readRDS("/Users/saeko/Unmeasured/data/mt_RNA2.rds")
colnames(mt_RNA)
unique(df2$Cell_type_class)
colnames(mt_RNA) <- c("Neural", "Blood", "Adipocyte", "Gonad", "Cardiovascular", "Digestive tract",
                      "Liver", "Lung", "Kidney", "Breast") #RefExデータの組織型をTFmarkerの表記と揃える
genes <- rownames(mt_RNA)

# TF marker, KnockTFのTF -----
mt_RNA2 <- mt_RNA %>% as_tibble() %>% 
  mutate(TF = genes) %>% 
  filter(TF %in% tgt_TFs) %>% 
  pivot_longer(-TF, names_to = "Cell_type_class", values_to = "FPKM")
  
df3 <- mt_RNA2 %>% left_join(df2, by = c("TF", "Cell_type_class")) %>% drop_na(label_chip_measure)
df3 %>% ggplot(aes(x = label_marker, y = FPKM, fill = label_marker)) +
  geom_violin(alpha = 0.7) +
  ggbeeswarm::geom_beeswarm() +
  geom_point() +
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        legend.position = "none",
        panel.grid.major = element_line(color = "gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold", color = "black"),
        axis.text.x =element_text(size=10,face="bold", color = "black"),
        axis.text.y =element_text(size=10,face="bold", color = "black"),
        axis.title=element_text(size=15,face="bold", color = "black"),
        aspect.ratio = 1.5
  )

print("statistical test: FPKM, TF marker vs others")
marker_FPKM <- df3 %>% filter(label_marker == "Marker") %>% .$FPKM
no_marker_FPKM <- df3 %>% filter(label_marker == "No marker") %>% .$FPKM
print(wilcox.test(marker_FPKM, no_marker_FPKM, paired = FALSE))


df3 %>% ggplot(aes(x = label_chip_measure, y = FPKM, fill = label_chip_measure)) +
  geom_violin(alpha = 0.7) +
  geom_point() +
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        legend.position = "none",
        panel.grid.major = element_line(color = "gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold", color = "black"),
        axis.text.x =element_text(size=10,face="bold", color = "black"),
        axis.text.y =element_text(size=10,face="bold", color = "black"),
        axis.title=element_text(size=15,face="bold", color = "black"),
        aspect.ratio = 1.5
  )

print("statistical test: FPKM, Unmeasured vs Measured")
measured_FPKM <- df3 %>% filter(label_chip_measure == "Measured") %>% .$FPKM
unmeasured_FPKM <- df3 %>% filter(label_chip_measure == "Unmeasured") %>% .$FPKM
print(wilcox.test(measured_FPKM, unmeasured_FPKM, paired = FALSE))


# measured / unmeasured -----
tgt_TFs <- mt_RNA2$TF %>% unique()
df_chipAtlas <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/SRX_date_20231004.tsv")
df_chipAtlas <- df_chipAtlas %>%
  filter(Cell_type_class != "Unclassified" & Cell_type_class != "No description") %>%
  filter(Antigen != "GFP" & Antigen != "Epitope tags") %>% 
  mutate(label_chip_measured = "Measured") %>% 
  select(Antigen, Cell_type_class,label_chip_measured) %>%
  filter(Antigen %in% tgt_TFs) %>%
  distinct()
colnames(df_chipAtlas) <- c("TF", "Cell_type_class", "label_chip_measured")

mt_RNA3 <- mt_RNA2 %>% left_join(df_chipAtlas, by = c("TF", "Cell_type_class"))
mt_RNA3$label_chip_measured[is.na(mt_RNA3$label_chip_measured)] <- "Unmeasured"

mt_RNA3 %>% ggplot(aes(x = label_chip_measured, y = FPKM, fill = label_chip_measured))  +
  geom_violin() +
  geom_point() +
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        legend.position = "none",
        panel.grid.major = element_line(color = "gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold", color = "black"),
        axis.text.x =element_text(size=10,face="bold", color = "black"),
        axis.text.y =element_text(size=10,face="bold", color = "black"),
        axis.title=element_text(size=15,face="bold", color = "black"),
        aspect.ratio = 1.5
  )

print("statistical test: FPKM, Unmeasured vs Measured")
measured_FPKM <- mt_RNA3 %>% filter(label_chip_measured == "Measured") %>% .$FPKM
unmeasured_FPKM <- mt_RNA3 %>% filter(label_chip_measured == "Unmeasured") %>% .$FPKM
print(wilcox.test(measured_FPKM, unmeasured_FPKM, paired = FALSE))

summary(measured_FPKM)
summary(unmeasured_FPKM)
