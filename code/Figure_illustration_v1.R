# TF-related paper ------
source("/Users/saeko/Unmeasured/code/function/TF_pub_scatter.R")
p_list <- TF_pub_scatter()
ggsave("/Users/saeko/Unmeasured/plot/TF_pub_scatter1.pdf", p_list[[1]], width = 9, height = 7)
ggsave("/Users/saeko/Unmeasured/plot/TF_pub_scatter2.pdf", p_list[[2]], width = 9, height = 7)


# Year_SRX_perTF ------
source("/Users/saeko/Unmeasured/code/function/Year_SRX.R")
p_list <- Year_SRX()
ggsave("/Users/saeko/Unmeasured/plot/Year_chip.pdf", p_list[[1]], width = 9, height = 7)
ggsave("/Users/saeko/Unmeasured/plot/Month_chip.pdf",p_list[[2]],width = 9, height = 7)
ggsave("/Users/saeko/Unmeasured/plot/Year_accum_chip.pdf",p_list[[3]],width = 9, height = 7)
ggsave("/Users/saeko/Unmeasured/plot/Year_accum_uniqTF.pdf",p_list[[4]],width = 9, height = 7)
ggsave("/Users/saeko/Unmeasured/plot/Month_accum_uniqTF.pdf",p_list[[5]],width = 9, height = 7)
ggsave("/Users/saeko/Unmeasured/plot/Year_chip_perTF.pdf",p_list[[6]],width = 9, height = 7)
ggsave("/Users/saeko/Unmeasured/plot/Year_uniqTF_perCtc.pdf",p_list[[7]],width = 9, height = 7)



# publication-ChIP-seq number plot ------
## all TFs and others


## only Cis-BP TFs (including DBDs)
source("/Users/saeko/Unmeasured/code/function/pub_chip_plot_CISBP.R")
p <- pub_chip_plot_CISBP()
ggsave("/Users/saeko/Unmeasured/plot/pub_chip_scatter_CISBP.pdf", p,  width = 9, height = 7)


# comparing No. of ChIP-seq (ignored paper, Fig, 2A) ------
source("/Users/saeko/Unmeasured/code/function/ignored_Fig2A_comp_year.R")
comp_year()

source("/Users/saeko/Unmeasured/code/function/ignored_Fig2A_comp_year_all.R")
comp_year_all()
