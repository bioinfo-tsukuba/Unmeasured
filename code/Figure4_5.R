# Figure4&5

# Reg-TF cover ratio plot ----
# Scatter plot (Reg-TF cover ratio - Unique TFs) ----
# calc_Reg_TF_cover_ratioで事前に計算
source("/Users/saeko/Unmeasured/code/function/Reg_TF_cover_ratio_year_plot.R")
calc_opt <- FALSE
p_list <- Reg_TF_cover_ratio_year_plot(calc_opt)
ggsave("/Users/saeko/Unmeasured/plot/Reg_TF_ratio_per_ctc.pdf", p_list[[1]], width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/sctter_uniqTF_Reg_TF_ratio.pdf", p_list[[2]], width = 9, height = 9)




# GWAS-SNP cover ratio plot ----
# Scatter plot (GWAS-SNP cover ratio - Unique TFs) ----
source("/Users/saeko/Unmeasured/code/function/GWAS_cover_ratio_year_plot.R")
p_list <- GWAS_cover_ratio_year_plot()
ggsave("/Users/saeko/Unmeasured/plot/GWAS_ratio_per_ctc.pdf", p_list[[1]], width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/sctter_uniqTF_GWAS_ratio.pdf", p_list[[2]], width = 9, height = 9)



# Simulation ----
## year plot & ratio vs unique TFs 
source("/Users/saeko/Unmeasured/code/function/simulation_Blood_GWAS_plot.R")
p_list <- simulation_Blood_GWAS_plot()
ggsave("/Users/saeko/Unmeasured/plot/simulation_Blood_GWAS_ratio.pdf", p_list[[1]], width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/simulation_Blood_cum_TF.pdf", p_list[[2]], width = 9, height = 9)

# Simulation, supplementary -----
source("/Users/saeko/Unmeasured/code/function/simulation_Blood_GWAS_plot_supple.R")
p_list <- simulation_Blood_GWAS_plot_supple()
ggsave("/Users/saeko/Unmeasured/plot/simulation_Blood_GWAS_ratio_supple1.pdf", p_list[[1]], width = 9, height = 9)
ggsave("/Users/saeko/Unmeasured/plot/simulation_Blood_GWAS_ratio_supple2.pdf", p_list[[2]], width = 9, height = 9)

