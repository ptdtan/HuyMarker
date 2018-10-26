run_Huy <- function(cells.1, cells.2)
{
  message("Huy Marker")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time(
      {
        mat <- mat.raw[, c(cells.1, cells.2)]
        clusters <- c(rep("A", length(cells.1)), rep("B", length(cells.2)))

        scores <- NoraSC::constructGeneMarkersData(expr.mat = mat, celltypes = clusters, cell = "A",
                                                   backend_Sum = DelayedMatrixStats::rowSums2)
        res <- scoringGeneMarkers(scores)
        rate <- log2(res$rate + 1)
        m <- mean(rate)
        s <- sd(rate)
        res$pval <- 1 - pnorm(rate, mean = m, sd = s)
      })
    list(session_info = session_info,
         timing = timing,
         df = data.frame(pval = res$pval, padj = res$pval, row.names = res$name),
         res = res
         )
  }, error = function(e) {
    "Huy results could not be calculated"
    list(session_info = session_info)
  })
}
