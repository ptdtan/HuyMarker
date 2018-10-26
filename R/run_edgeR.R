suppressPackageStartupMessages(library(edgeR))

# run_edgeRLRT ------------------------------------------------------------

run_edgeRLRT <- function(cells.1, cells.2) {
  message("edgeRLRT")
  session_info <- sessionInfo()

  mat <- mat.raw[, c(cells.1, cells.2)]
  clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

  timing <- system.time({
    dge <- DGEList(mat, group = clusters)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~clusters)
    dge <- estimateDisp(dge, design = design)
    fit <- glmFit(dge, design = design)
    lrt <- glmLRT(fit)
    tt <- topTags(lrt, n = Inf)
  })
  list(session_info = session_info,
       timing = timing,
       tt = tt,
       df = data.frame(pval = tt$table$PValue,
                       padj = tt$table$FDR,
                       row.names = rownames(tt$table)))
}

# run_edgeRQLFDetRate -----------------------------------------------------

run_edgeRQLFDetRate <- function(cells.1, cells.2) {
  message("edgeRQLFDetRate")
  session_info <- sessionInfo()

  mat <- mat.raw[, c(cells.1, cells.2)]
  clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

  timing <- system.time({
    dge <- DGEList(mat, group = clusters)
    dge <- calcNormFactors(dge)
    cdr <- scale(colMeans(mat > 0))
    design <- model.matrix(~ cdr + clusters)
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design)
    qlf <- glmQLFTest(fit)
    tt <- topTags(qlf, n = Inf)
  })
  list(session_info = session_info,
       timing = timing,
       tt = tt,
       df = data.frame(pval = tt$table$PValue,
                       padj = tt$table$FDR,
                       row.names = rownames(tt$table)))
}

# run_edgeRQLF ------------------------------------------------------------

run_edgeRQLF <- function(cells.1, cells.2) {
  message("edgeRQLF")
  session_info <- sessionInfo()

  mat <- mat.raw[, c(cells.1, cells.2)]
  clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

  timing <- system.time({
    dge <- DGEList(mat, group = clusters)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~clusters)
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design)
    qlf <- glmQLFTest(fit)
    tt <- topTags(qlf, n = Inf)
  })
  list(session_info = session_info,
       timing = timing,
       tt = tt,
       df = data.frame(pval = tt$table$PValue,
                       padj = tt$table$FDR,
                       row.names = rownames(tt$table)))
}

