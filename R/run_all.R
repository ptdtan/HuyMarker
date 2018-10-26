run_all <- function(obj, type = "real")
{
  source("R/run_Huy.R")
  source("R/run_Seurat.R")
  source("R/run_edgeR.R")
  source("R/run_MAST.R")
  source("R/run_limma.R")
  source("R/run_ttest.R")
  mat.raw <<- obj[["mat"]]
  cells.1 <- obj[[type]][[1]]
  cells.2 <- obj[[type]][[2]]
  cells.1 <- sample(cells.1, min(length(cells.1), 100), replace = F)
  cells.2 <- sample(cells.2, min(length(cells.2), 100), replace = F)

  res <- list(res_Huy = run_Huy(cells.1, cells.2),
              res_SeuratBimod = run_Seurat(cells.1, cells.2, method = "bimod"),
              res_SeuratT = run_Seurat(cells.1, cells.2, method = "t"),
              res_SeuratPoisson = run_Seurat(cells.1, cells.2, method = "poisson"),
              res_Seuratnegbinom = run_Seurat(cells.1, cells.2, method = "negbinom"),
              res_edgeRLRT = run_edgeRLRT(cells.1, cells.2),
              res_edgeRQLFDetRate = run_edgeRQLFDetRate(cells.1, cells.2),
              res_edgeRQLF = run_edgeRQLF(cells.1, cells.2),
              #res_ttest = run_ttest(cells.1, cells.2),
              res_MASTcpmDetRate = run_MASTcpmDetRate(cells.1, cells.2),
              res_limmatrend = run_limmatrend(cells.1, cells.2)
              )

  return(res)
}

run_FDR_sim <- function()
{
  list_files <- paste0("data/sim_", seq(1, 10),".rds")
  obj <- readRDS(list_files[1])
}

run_FDR_real <- function()
{
  files <- list(
    GSE62270 = "data/GSE62270_data.rds",
    GSE81076_GPL16791 = "data/GSE81076_GPL16791_data.rds",
    GSE81076_GPL18573 = "data/GSE81076_GPL18573_data.rds",
    pbmc4k_1 = "data/PBMC4k_1_data.rds",
    pbmc4k_2 = "data/PBMC4k_2_data.rds",
    pbmc4k_3 = "data/PBMC4k_3_data.rds",
    zeisel2015_7 = "data/zeisel2015_7_data.rds",
    zeisel2015_8 = "data/zeisel2015_8_data.rds",
    zeisel2015_9 = "data/zeisel2015_9_data.rds"
    )
  stats <- parallel::mclapply(files, function(file){
                  obj = readRDS(file)
                  timming <- system.time({
                    res = run_all(obj, type = "null")
                    stats = get_FDR_onedata(res)
                  })
                  message(timming)
                }, mc.cores = 12)
  saveRDS(stats, file = "data/stats_results.rds")
}

get_FDR_onedata <- function(res)
{
  lapply(res, function(obj){
    tab <- table(obj$df$padj < 0.05)
    tab[[2]]/sum(tab[[2]] + tab[[1]])
  })
}
