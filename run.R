args <- (commandArgs(trailingOnly = TRUE))

data_f <- args[1]
data_truth <- args[2]

HuyMarkerBenchmark::run_TPR_real(data_f, data_truth, F)

#
# # TENX molecule coverage plot ---------------------------------------------
# f <- "/Users/bioturing/Data/Dotplot_Results/molecule_info/molecule.tsv"
# TENX_mol <- "~/Data/Dotplot_Results/molecule_info/molecule.tsv"
# TENX_mol_df <- readr::read_tsv(TENX_mol)
# l <- TENX_mol_df[["read_count"]]*150/(TENX_mol_df[["end"]] - TENX_mol_df[["start"]])
# ll <- sample(l , 100000)
# library(plotly)
# library(ggplot2)
# s <- sd(ll)
# mu <- mean(ll)
#
# ll_df <- data.frame(l = ll)
# ggplot(ll_df, aes(x = l)) + geom_histogram(binwidth = 0.01) +xlim(0,0.5)
#
# # UST molecule coverage plot
# f <- "/Users/bioturing/Data/Dotplot_Results/molecule_info/run180522.molecule.tsv"
# UST_mol_df <- readr::read_tsv(f, col_names = F)
# l <- UST_mol_df[[7]]
# library(plotly)
# library(ggplot2)
# ll_df <- data.frame(l = l)
# ggplot(ll_df, aes(x = l)) + geom_histogram(binwidth = 0.01) +xlim(0,0.5)
# s <- sd(l)
# mu <- mean(l)
#
# # SIM2 molecule coverage plot
# f <- "/Users/bioturing/Data/Dotplot_Results/molecule_info/molecule_sim5.inf"
# SIM2_mol_df <- readr::read_tsv(f)
# l <- SIM2_mol_df[[5]]*144/SIM2_mol_df[[4]]
# plot_ly(x = SIM2_mol_df[[5]], type = "histogram")
#
# library(plotly)
# library(ggplot2)
# ll_df <- data.frame(l = l)
# p <- ggplot(ll_df, aes(x = l)) + geom_histogram(binwidth = 0.01) +xlim(0,0.4)
# ggplotly(p)
#
# d <- data.frame(x = rgamma(1000, shape = 4, rate = 0.7) * 0.01)
# p <- ggplot(d, aes(x = x)) + geom_histogram(binwidth = 0.01) +xlim(0,0.5)
# p
#
# d <- data.frame(x = rgamma(1000, shape = 0.8, rate = 9))
# plot_ly(x = d$x, type = "histogram")
# p <- ggplot(d, aes(x = x)) + geom_histogram(binwidth = 0.005) +xlim(0,0.1)
# p
#
# f <- "/Users/bioturing/Data/Assembly/Ecoli/mdup_out/molecule.tsv"
# mol_df <- readr::read_tsv(f, col_names = F)
