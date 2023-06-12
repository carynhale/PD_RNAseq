runLimma <- function(meta, #metadata must have a "Sample" that exact matches the rownames of the counts object.
                     counts,
                     design, #design object created with the model.matrix() function
                     group_variable, #one of the variables used for the model (for gene filtering)
                     contrasts, #contrasts matrix made with the makeContrasts() function
                     experiment_name,
                     date,
                     path,
                     pathways,
                     min.pathways) {
  require(tidyverse)
  require(fgsea)
  require(GSVA)
  require(limma)
  require(Biobase)
  require(edgeR)
  stopifnot(colnames(counts)==meta$Sample) 
  dir.create(path = file.path(path,experiment_name))
  setwd(file.path(path,experiment_name))
  n_samples <- dim(meta)[1]
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, group = group_variable)
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  
  pdf("density.plot.pdf")
  plotDensities(cpm(dge, log = TRUE), legend = FALSE)
  dev.off()
  
  dge <- calcNormFactors(dge, method = "TMM")
  lcpm <- cpm(dge, log = TRUE, prior.count = 3)
  
  vv <- voomWithQualityWeights(dge, design = design, plot = TRUE, normalize.method = "none")
  
  vfit = lmFit(vv,design)
  vfit = contrasts.fit(vfit,contrasts)
  vfit= eBayes(vfit)

  summa.fit <- decideTests(vfit)
  summary(summa.fit)
  
  for (i in seq_along(colnames(contrasts))) {
    results <- topTable(vfit,coef = colnames(contrasts)[i], n = "Inf") %>%
      as_tibble()
    write.table(results,
                file = paste(experiment_name, "_limma_", colnames(contrasts)[i], "_", date, ".txt", sep = ""),
                sep = "\t",
                row.names = FALSE,
                quote = FALSE)
  }
  gene_exp <- lcpm %>%
    as_tibble(rownames = "ID") %>%
    pivot_longer(cols = 2:(n_samples+1), names_to = "Sample", values_to = "lcpm") %>%
    left_join(meta,
              by = "Sample")
  save(gene_exp, file = paste(experiment_name, "lcpm", date, "RData", sep = "."))
  save(meta, file = paste(experiment_name, "samples", date, "RData", sep = "."))
  
  for (i in seq_along(colnames(contrasts))) {
    results <- topTable(vfit,coef = colnames(contrasts)[i], n = "Inf") %>%
      as_tibble()
    rnk <- results$t
    names(rnk) <- results$ID
    gsea_res <- fgsea(pathways, 
                      rnk,
                      minSize = min.pathways,
                      maxSize = Inf,
                      eps = 0)
    write.table(gsea_res %>%
                  dplyr::select(-leadingEdge),
                file = paste(experiment_name, ".gsea.", colnames(contrasts)[i], date, ".txt", sep = ""),
                sep = "\t",
                row.names = FALSE,
                quote = FALSE)
  }
  
  lcpm <- lcpm[is.na(rownames(lcpm))==FALSE,]
  eset <- ExpressionSet(assayData = lcpm)
  gsva.go <- gsva(eset, pathways, min.sz = min.pathways)
  arrayw<-arrayWeights(gsva.go,design=design)
  vfit = lmFit(gsva.go,design,weights = arrayw)
  vfit = contrasts.fit(vfit,contrasts)
  vfit= eBayes(vfit)
  summa.fit <- decideTests(vfit)
  summary(summa.fit)
  
  for (i in seq_along(colnames(contrasts))) {
    results <- topTable(vfit,coef = colnames(contrasts)[i], n = "Inf") %>%
      as_tibble(rownames = "pathway")
    write.table(results,
                file = paste(experiment_name, ".gsva.", colnames(contrasts)[i], date, ".txt", sep = ""),
                sep = "\t",
                row.names = FALSE,
                quote = FALSE)
  }
  
  exprs <- exprs(gsva.go)
  
  toPlot_pathways <- exprs %>%
    as_tibble(rownames = "pathway") %>%
    pivot_longer(cols = 2:(n_samples+1), names_to = "Sample", values_to = "score") %>%
    left_join(meta,
              by = "Sample")
  save(toPlot_pathways, file = paste(experiment_name, "gsva_exprs", date, "RData", sep = "."))
}
