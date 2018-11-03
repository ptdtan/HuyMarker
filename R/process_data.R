suppressPackageStartupMessages(library(NoraSC))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(scDD))

process_pbmc4k <- function()
{
  master_dir <- "data.1/PBMC4k"
  mat.sp <- NoraSC::Read10XData(file.path(master_dir, "main"), type = "hdf5", return.raw = T)
  mat.raw <- as(mat.sp, "matrix")
  graph <- jsonlite::fromJSON(file.path(master_dir, "main/cluster_result.json"))
  clusters <- graph$kmeans$clusters[3, ]
  cluster <- c(1,2,3)

  for(i in cluster){
    message("Cluster", i)
    null.data <- which(clusters == i)
    null.data.1 <- sample(x = null.data, min(100, length(null.data)), replace = F)
    d <- setdiff(null.data, null.data.1)
    null.data.2 <- sample(x = d ,
                          min(100, length(d)), replace = F)
    real.data.1 <- which(clusters == i)
    real.data.2 <- which(clusters != i)
    PBMC4k <- list(real = list(real.data.1, real.data.2),
                     null = list(null.data.1, null.data.2),
                     mat = assay(mat.raw))
    data.file <- file.path(paste0("data.1/PBMC4k_", i, "_data.rds"))
    saveRDS(PBMC4k, data.file)
  }
}

process_zeisel2015 <- function()
{
  master_dir <- "data.1/zeisel2015"
  mat.sp <- NoraSC::Read10XData(file.path(master_dir, "main"), type = "hdf5", return.raw = T)
  mat.raw <- as(mat.sp, "matrix")
  graph <- jsonlite::fromJSON(file.path(master_dir, "main/cluster_result.json"))
  clusters <- graph$graph$clusters
  cluster <- c(7,8,9)

  for(i in cluster){
    message("Cluster", i)
    null.data <- which(clusters == i)
    null.data.1 <- sample(x = null.data, min(100, sum(clusters == i)), replace = F)
    d <- setdiff(null.data, null.data.1)
    null.data.2 <- sample(x = d ,
                          min(100, length(d)), replace = F)
    real.data.1 <- which(clusters == i)
    real.data.2 <- which(clusters != i)
    data <- list(real = list(real.data.1, real.data.2),
                   null = list(null.data.1, null.data.2),
                   mat = assay(mat.raw))
    data.file <- file.path(paste0("data.1/zeisel2015_", i, "_data.rds"))
    saveRDS(data, data.file)
  }
}

process_GSE62270 <- function()
{
  file <- "data.1/GSE62270-GPL17021.rds"
  data <- readRDS(file)
  config <- list(
    groupid = "source_name_ch1",
    keepgroups = c("Randomly extracted cells from whole intestinal organoids",
                   "Randomly extracted ex vivo isolated 5 day YFP positive cells")
  )
  clusters <- colData(data)[[config[["groupid"]]]]
  cond1 <- config$keepgroups[1]
  real.data.1 <- which(clusters == cond1)
  real.data.2 <- which(clusters != cond1)
  cond2 <- "Randomly extracted Lgr5-positive intestinal cells"
  null.data <- which(clusters == cond2)
  null.data.1 <- sample(x = null.data, 100, replace = F)
  d <- setdiff(null.data, null.data.1)
  null.data.2 <- sample(x = d ,
                        min(100, length(d)), replace = F)
  data.file <- file.path("data.1/GSE62270_data.rds")

  GSE62270 <- list(real = list(real.data.1, real.data.2),
                   null = list(null.data.1, null.data.2),
                   mat = assay(data[[1]]))
  saveRDS(GSE62270, data.file)
}

process_GSE81076 <- function()
{

  # GSE81076_GPL18573 -------------------------------------------------------
  file <- "data.1/GSE81076-GPL18573.rds"
  data <- readRDS(file)
  config <- list(
    groupid = "characteristics_ch1",
    keepgroups = c("cell type: CD13+ sorted cells",
                   "cell type: CD24+ CD44+ live sorted cells")
  )
  clusters <- colData(data)[[config[["groupid"]]]]
  cond1 <- config$keepgroups[1]
  real.data.1 <- which(clusters == cond1)
  real.data.2 <- which(clusters != cond1)
  cond2 <- "cell type: live sorted cells"
  null.data <- which(clusters == cond2)
  null.data.1 <- sample(x = null.data, 100, replace = F)
  d <- setdiff(null.data, null.data.1)
  null.data.2 <- sample(x = d ,
                        min(100, length(d)), replace = F)
  GSE81076_GPL18573 <- list(real = list(real.data.1, real.data.2),
                  null = list(null.data.1, null.data.2),
                  mat = assay(data[[1]]))
  data.file <- file.path("data.1/GSE81076_GPL18573_data.rds")
  saveRDS(GSE81076_GPL18573, data.file)

  # GSE81076-GPL16791.rds ---------------------------------------------------
  file <- "data.1/GSE81076-GPL16791.rds"
  data <- readRDS(file)
  config <- list(
    groupid = "characteristics_ch1",
    keepgroups = c("cell type: exocrine fraction, live sorted cells",
                   "cell type: live sorted cells")
  )
  clusters <- colData(data)[[config[["groupid"]]]]
  cond1 <- config$keepgroups[1]
  real.data.1 <- which(clusters == cond1)
  real.data.2 <- which(clusters != cond1)
  cond2 <- "cell type: live sorted cells"
  null.data <- which(clusters == cond2)
  null.data.1 <- sample(x = null.data, 100, replace = F)
  d <- setdiff(null.data, null.data.1)
  null.data.2 <- sample(x = d ,
                        min(100, length(d)), replace = F)
  GSE81076_GPL16791 <- list(real = list(real.data.1, real.data.2),
                            null = list(null.data.1, null.data.2),
                            mat = assay(data[[1]]))
  data.file <- file.path("data.1/GSE81076_GPL16791_data.rds")
  saveRDS(GSE81076_GPL16791, data.file)

}

process_GSE78779 <- function()
{
  file <- "data.1/GSE78779-GPL17021.rds"
  data <- readRDS(file)
  config <- list(
    groupid = "characteristics_ch1.1",
    keepgroups = c("tissue: Ear",
                   "tissue: femora and tibiae")
  )
  clusters <- colData(data)[[config[["groupid"]]]]
  cond1 <- config$keepgroups[1]
  real.data.1 <- which(clusters == cond1)
  real.data.2 <- which(clusters != cond1)
  cond2 <- "tissue: Ear"
  null.data <- which(clusters == cond2)
  null.data.1 <- sample(x = null.data, 100, replace = F)
  d <- setdiff(null.data, null.data.1)
  null.data.2 <- sample(x = d ,
                        min(100, length(d)), replace = F)
  GSE78779 <- list(real = list(real.data.1, real.data.2),
                   null = list(null.data.1, null.data.2),
                   mat = assay(data[[1]]))
  data.file <- file.path("data.1/GSE78779_data.rds")
  saveRDS(GSE78779, data.file)
}

#' Title
#'
#' @param nDE
#' @param nDP
#' @param nDM
#' @param nDB
#' @param nEE
#' @param nEP
#' @param numSamples
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
simulate <-  function(scDatEx, numSamples=100,
                      nDE=250, nDP=250, nDM=250, nDB=250,
                      nEE=5000, nEP=4000,
                      seed = 816)
{
  SD <- simulateSet(scDatEx, numSamples=numSamples, nDE=nDE, nDP=nDP, nDM=nDM,
                   nDB=nDB, nEE=nEE, nEP=nEP, sd.range=c(2,2), modeFC=c(2,3,4),
                   plots=FALSE,
                   random.seed=seed)
  return(SD)
}

process_simulate <- function()
{
  library(scDD)
  data(scDatEx)
  for(i in seq(1, 10)){
    file <- file.path(paste0("data.1/sim_", i, ".rds"))
    message(paste("Doing", file))
    SD <- simulateSet(scDatEx)
    real.data <- colData(SD)[["condition"]]
    real.data.1 <- which(real.data == 1)
    real.data.2 <- which(real.data != 1)
    sim <- list(real = list(real.data.1, real.data.2),
                null = NULL,
                mat = assay(SD))
    saveRDS(sim, file = file)
  }
}

prepare_all <- function()
{
  message("PBMC4k")
  process_pbmc4k()

  message("process_zeisel2015")
  process_zeisel2015()

  message("process_GSE62270")
  process_GSE62270()

  message("process_GSE81076")
  process_GSE81076()

  message("process_GSE78779")
  process_GSE78779()
}


install_packages <- function()
{
  packages <- c("edgeR", "MultiAssayExperiment", "scDD", "genefilter")
  BiocInstaller::biocLite(packages)
}
