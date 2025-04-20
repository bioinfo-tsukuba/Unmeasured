# Figure1

# Dotplot (TF family-Cell type class, Number of ChIP-seq)
source("/Users/saeko/Unmeasured/code/function/Dotplot_fam_ctc.R")
dotplot <- Dotplot_fam_ctc()
ggsave("/Users/saeko/Unmeasured/plot/Dotplot.pdf", dotplot ,width = 15, height = 10)

# ChIP-seq distribution 
source("/Users/saeko/Unmeasured/code/function/ChIP_distribution.R")
p_list <- ChIP_distribution()
ggsave("/Users/saeko/Unmeasured/plot/ChIP_distribution_TF.pdf", p_list[[1]], width = 20, height = 5)
ggsave("/Users/saeko/Unmeasured/plot/ChIP_distribution_TF_circ.pdf", p_list[[2]], width = 20, height = 20)
ggsave("/Users/saeko/Unmeasured/plot/ChIP_distribution_Ct.pdf", p_list[[3]], width = 20, height = 5)
ggsave("/Users/saeko/Unmeasured/plot/ChIP_distribution_Ct_circ.pdf", p_list[[4]], width = 20, height = 20)



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

# xgboost
source("/Users/saeko/Unmeasured/code/function/xgboost.R")
p_list <- xgboost_plot()
ggsave("/Users/saeko/Unmeasured/plot/xgboost_scatter.pdf", p_list[[1]], width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/xgboost_feature_importance.pdf", p_list[[2]], width = 9, height = 18)


# comparing No. of ChIP-seq 
source("/Users/saeko/Unmeasured/code/function/ignored_Fig2A_comp_year.R")
comp_year()

source("/Users/saeko/Unmeasured/code/function/ignored_Fig2A_comp_year_all.R")
comp_year_all()