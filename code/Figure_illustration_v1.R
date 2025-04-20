# TF-related paper ------
source("/Users/saeko/Unmeasured/code/function/TF_pub_scatter.R")
p_list <- TF_pub_scatter()
ggsave("/Users/saeko/Unmeasured/plot/TF_pub_scatter1.pdf", p_list[[1]], width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/TF_pub_scatter2.pdf", p_list[[2]], width = 9, height = 9)


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


# ChIP-seq distribution ------
source("/Users/saeko/Unmeasured/code/function/ChIP_distribution.R")
p_list <- ChIP_distribution()
ggsave("/Users/saeko/Unmeasured/plot/ChIP_distribution_TF.pdf", p_list[[1]], width = 20, height = 5)
ggsave("/Users/saeko/Unmeasured/plot/ChIP_distribution_TF_circ.pdf", p_list[[2]], width = 20, height = 20)
ggsave("/Users/saeko/Unmeasured/plot/ChIP_distribution_Ct.pdf", p_list[[3]], width = 20, height = 5)
ggsave("/Users/saeko/Unmeasured/plot/ChIP_distribution_Ct_circ.pdf", p_list[[4]], width = 20, height = 20)

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


# Number of Cell type per Cell type class(Supplementary)
source("/Users/saeko/Unmeasured/code/function/Ct_per_Ctc.R")
p <- Ct_per_Ctc()
ggsave("/Users/saeko/Unmeasured/plot/Number_of_Ct_per_Ctc.pdf", p, width = 9, height = 9)

# Number of TFs per TF families(Supplementary)
source("/Users/saeko/Unmeasured/code/function/TFs_per_fam.R")
p <- TFs_per_fam()
ggsave("/Users/saeko/Unmeasured/plot/Number_of_TFs_per_TFfam.pdf", p, width = 9, height = 9)

# Dotplot (TF family-Cell type class, Number of ChIP-seq)
source("/Users/saeko/Unmeasured/code/function/Dotplot_fam_ctc.R")
dotplot <- Dotplot_fam_ctc()
ggsave("/Users/saeko/Unmeasured/plot/Dotplot.pdf", dotplot ,width = 15, height = 10)

# Dotplot (TF family-Cell type class, Unique TF)
source("/Users/saeko/Unmeasured/code/function/Dotplot_fam_ctc_uniqTF.R")
dotplot <- Dotplot_fam_ctc_uniqTF()
ggsave("/Users/saeko/Unmeasured/plot/Dotplot_uniqTF.pdf", dotplot ,width = 15, height = 10)

# Dotplot (TF family-Cell type)
source("/Users/saeko/Unmeasured/code/function/Dotplot_fam_ct.R")
dotplot <- Dotplot_fam_ct()
ggsave("/Users/saeko/Unmeasured/plot/Dotplot_S1_fam_ct.pdf", dotplot ,width = 15, height = 10)

# xgboost
source("/Users/saeko/Unmeasured/code/function/xgboost.R")
p_list <- xgboost_plot()
ggsave("/Users/saeko/Unmeasured/plot/xgboost_scatter.pdf", p_list[[1]], width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/xgboost_feature_importance.pdf", p_list[[2]], width = 9, height = 18)


# Lorenz curve
source("/Users/saeko/Unmeasured/code/function/Lorenz_curve.R")
p_list <- Lorenz_curve()
ggsave("/Users/saeko/Unmeasured/plot/Lorenz_curve_TF_remap2.pdf", p_list[[1]], width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/Lorenz_curve_TF_GTRD.pdf", p_list[[2]], width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/Lorenz_curve_TF_ChIPatlas.pdf", p_list[[3]], width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/Lorenz_curve_TF_all.pdf", p_list[[4]], width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/Lorenz_curve_Ct_remap2.pdf", p_list[[5]], width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/Lorenz_curve_Ct_GTRD.pdf", p_list[[6]], width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/Lorenz_curve_Ct_ChIPatlas.pdf", p_list[[7]], width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/Lorenz_curve_Ct_all.pdf", p_list[[8]], width = 9, height = 9)
