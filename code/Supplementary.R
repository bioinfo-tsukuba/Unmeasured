# Supplementary


# Year_SRX_perTF ------
source("/Users/saeko/Unmeasured/code/function/Year_SRX.R")
p_list <- Year_SRX()
ggsave("/Users/saeko/Unmeasured/plot/Year_chip.pdf", p_list[[1]], width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/Month_chip.pdf",p_list[[2]],width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/Year_accum_chip.pdf",p_list[[3]],width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/Year_accum_uniqTF.pdf",p_list[[4]],width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/Month_accum_uniqTF.pdf",p_list[[5]],width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/Year_chip_perTF.pdf",p_list[[6]],width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/Year_uniqTF_perCtc.pdf",p_list[[7]],width = 12, height = 9)


# Number of Cell type per Cell type class(Supplementary)
source("/Users/saeko/Unmeasured/code/function/Ct_per_Ctc.R")
p <- Ct_per_Ctc()
ggsave("/Users/saeko/Unmeasured/plot/Number_of_Ct_per_Ctc.pdf", p, width = 9, height = 9)

# Number of TFs per TF families(Supplementary)
source("/Users/saeko/Unmeasured/code/function/TFs_per_fam.R")
p <- TFs_per_fam()
ggsave("/Users/saeko/Unmeasured/plot/Number_of_TFs_per_TFfam.pdf", p, width = 9, height = 9)


# Dotplot (TF family-Cell type class, Unique TF)
source("/Users/saeko/Unmeasured/code/function/Dotplot_fam_ctc_uniqTF.R")
dotplot <- Dotplot_fam_ctc_uniqTF()
ggsave("/Users/saeko/Unmeasured/plot/Dotplot_uniqTF.pdf", dotplot ,width = 15, height = 10)

# Dotplot (TF family-Cell type)
source("/Users/saeko/Unmeasured/code/function/Dotplot_fam_ct.R")
dotplot <- Dotplot_fam_ct()
ggsave("/Users/saeko/Unmeasured/plot/Dotplot_S1_fam_ct.pdf", dotplot ,width = 15, height = 10)


# Supplementary FPKM threshold
source("/Users/saeko/Unmeasured/code/function/FPKM_threshold_histogream.R")
p_FPKM <- FPKM_threshold_histogram()
ggsave("/Users/saeko/Unmeasured/plot/FPKM_threshold_histogram.pdf", p_FPKM, width = 9, height = 9)


# Dotplot Unmeasured(%) (TF family-Cell type), including measure;0 
source("/Users/saeko/Unmeasured/code/function/Dotplot_fam_ctc_unmeasured_v2.R")
calc_opt <- "FALSE"
tgt_threshold <- 0.75 #上位25%
#tgt_threshold <- 0.50 #上位50%
dotplot <- Dotplot_fam_ctc_unmeasured_v2(calc_opt, tgt_threshold)
ggsave(paste0("/Users/saeko/Unmeasured/plot/Dotplot_fam_ct_unmeasured_", tgt_threshold,"_v2.pdf"), dotplot ,width = 15, height = 10)
