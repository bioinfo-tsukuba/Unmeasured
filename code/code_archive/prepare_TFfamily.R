# Fig. 3D and 3E: Cell type comparison
source("~/MOCCS_paper_public/function/sub_function/Read_df.R")
df_raw <- Read_df()

source("~/MOCCS_paper_public/function/sub_function/Annot_DBDs_of_CIS_BP.R")
df_fam <- Annot_DBDs_of_CIS_BP(df_raw, TF_info_dir ="~/MOCCS_paper_public/data/Fig3/TF_Information.txt")

source("~/MOCCS_paper_public/function/sub_function/Annot_pairs.R")
df_p_1 <- Annot_pairs(df_fam)

source("~/MOCCS_paper_public/function/sub_function/Calc_pairs.R")
df_p_2 <- Calc_pairs(df_raw)


# Join
df_p_1 %>%
  mutate(ID_pair = paste0(ID1, ID2)) -> df_p_1_2
df_p_2 %>%
  mutate(ID_pair = paste0(ID1, ID2)) %>%
  dplyr::select(ID_pair, k_sim_1, k_sim_2) -> df_p_2_2
inner_join(df_p_1_2, df_p_2_2, by = "ID_pair") -> df_p_3

source("~/MOCCS_paper_public/function/sub_function/Group_df.R")
df_p_3_gp <- Group_df(df_p_3)

source("~/MOCCS_paper_public/function/sub_function/Compare_cell_type.R")
flag_1 <- !(df_p_3_gp$ID1 %in% all_ns | df_p_3_gp$ID2 %in% all_ns)
Compare_cell_type(df_p_3_gp[flag_1, ],
                  plot_ant = c("FOS", "JUN", "GATA2", "MYC"))
