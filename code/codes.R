write_fit <- function (fit, results = NULL, digits = 3, adjust = "none",
                       method = "separate", F.adjust = "none", sep = "\t", ...)
{
  if (!is(fit, "MArrayLM"))
    stop("fit should be an MArrayLM object")
  if (!is.null(results) && !is(results, "TestResults"))
    stop("results should be a TestResults object")
  if (is.null(fit$t) || is.null(fit$p.value))
    fit <- eBayes(fit)
  method <- match.arg(method, c("separate", "global"))
  p.value <- as.matrix(fit$p.value)
  if (adjust == "none") {
    p.value.adj <- NULL
  }
  else {
    p.value.adj <- p.value
    if (method == "separate")
      for (j in 1:ncol(p.value)) p.value.adj[, j] <- p.adjust(p.value[,
                                                                      j], method = adjust)
    if (method == "global")
      p.value.adj[] <- p.adjust(p.value, method = adjust)
  }
  if (F.adjust == "none" || is.null(fit$F.p.value))
    F.p.value.adj <- NULL
  else F.p.value.adj <- p.adjust(fit$F.p.value, method = F.adjust)
  rn <- function(x, digits = digits) if (is.null(x))
    NULL
  else {
    if (is.matrix(x) && ncol(x) == 1)
      x <- x[, 1]
    round(x, digits = digits)
  }
  tab <- list()
  tab$A <- rn(fit$Amean, digits = digits - 1)
  tab$Coef <- rn(fit$coef, digits = digits)
  tab$t <- rn(fit$t, digits = digits - 1)
  tab$p.value <- rn(p.value, digits = digits + 2)
  tab$p.value.adj <- rn(p.value.adj, digits = digits + 3)
  tab$F <- rn(fit$F, digits = digits - 1)
  tab$F.p.value <- rn(fit$F.p.value, digits = digits + 2)
  tab$F.p.value.adj <- tab$F.p.value.adj <- rn(F.p.value.adj,
                                               digits = digits + 3)
  tab$Res <- unclass(results)
  tab$Genes <- fit$genes
  tab <- data.frame(tab, check.names = FALSE)
  #write.table(tab, file = file, quote = FALSE, row.names = FALSE,
  #sep = sep, ...)
  return(tab)
}


## Enrichment analysis

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

goSummaries <- url("https://uofabioinformaticshub.github.io/summaries2GO/data/goSummaries.RDS") %>%
  readRDS() %>%
  mutate(ontology = as.character(ontology))


getGeneLists <- function(pwf, goterms, genome, ids){
  gene2cat <- getgo(rownames(pwf), genome, ids)
  cat2gene <- split(rep(names(gene2cat), sapply(gene2cat, length)),
                    unlist(gene2cat, use.names = FALSE))
  out <- list()
  for(term in goterms){
    tmp <- pwf[cat2gene[[term]],]
    tmp <- rownames(tmp[tmp$DEgenes > 0, ])
    out[[term]] <- tmp
  }
  out


}

getLength_data = function(genome){
  ensembl = useEnsembl(biomart = "genes")
  genome_db = makeTxDbFromBiomart(dataset=genome)
  txsByGene=transcriptsBy(genome_db,"gene")
  lengthData=median(width(txsByGene))
  return(lengthData)}


get_ontologies = function(genes, pathway){

  ## Overlap the DE genes with the genes in the length data vector
  select_genes <- as.vector(names(lengthData)%in%names(genes))
  new_lengthData <- lengthData[select_genes]
  select_genes <- as.vector(names(genes) %in% names(new_lengthData))
  genes = genes[select_genes]

  ## Calculate a Probability Weighting Function
  ## This files a monotonic spline
  pwf=nullp(genes,bias.data=new_lengthData, plot.fit = FALSE)


  if (missing(pathway)){
    message("Using GO from GOseq")

    Go.wall=goseq(pwf, "mm39", "ensGene", use_genes_without_cat = TRUE)
    goterms <- Go.wall$category
    goList <- getGeneLists(pwf, goterms, "mm39", "ensGene")
    Go.wall$EnsemblID <- sapply(Go.wall$category, function(x) goList[[x]])
    Go.wall$EnsemblID <- vapply(Go.wall$EnsemblID, paste, collapse = ", ", character(1L))
    goRes <- Go.wall %>%
      dplyr::select(-under_represented_pvalue) %>%
      dplyr::filter(numDEInCat > 0) %>%
      dplyr::rename(id = category) %>%
      left_join(goSummaries)
    nGenes <- c(BG = sum(genes == 0), DE = sum(genes))
    goFinal <- goRes %>%
      dplyr::filter(numDEInCat > 2) %>%
      arrange(desc(shortest_path), EnsemblID) %>% distinct(EnsemblID, .keep_all = TRUE) %>%
      arrange(over_represented_pvalue) %>% mutate(Expected = round(nGenes[["DE"]] * numInCat / nGenes[["BG"]], 0)) %>%
      mutate(adjP = p.adjust(over_represented_pvalue, "bonferroni"),FDR = p.adjust(over_represented_pvalue, "fdr"),
             Sig_adjP = adjP < 0.05, Sig_FDR = FDR < 0.05) %>%
      dplyr::filter(Sig_adjP, numDEInCat > Expected)


  } else {
    Go.wall=goseq(pwf, gene2cat = pathway, use_genes_without_cat = TRUE)
    goFinal = Go.wall %>%
      dplyr::filter(numDEInCat > 0) %>%
      dplyr::rename(id = category) %>%
      mutate(adjP = p.adjust(over_represented_pvalue, "bonferroni"),FDR = p.adjust(over_represented_pvalue, "fdr"),
             Sig_adjP = adjP < 0.05, Sig_FDR = FDR < 0.05) %>%
      dplyr::filter(Sig_adjP) %>%
      dplyr::arrange(Sig_adjP)

  }

  return(goFinal)

}

calcDiff <- function(x) {
  x %>% mutate(diff = PropUp-PropDown)
}

ens2entrez = function(v){
  c = v$E %>% as.data.frame()
  c %<>%
    mutate(entrez = mapIds(org.Mm.eg.db, keys=rownames(.) , column="ENTREZID",keytype="ENSEMBL")) %>%
    rownames_to_column("EnsID")
  rownames(c) = c(1:nrow(c))

  c %<>% rownames_to_column("rowNum")
  c %<>% drop_na(entrez) %>%
    distinct(entrez, .keep_all = TRUE) %>%
    set_rownames(.$entrez) %>%
    dplyr::select(-entrez)
  row_num = c$rowNum
  c %<>% dplyr::select(-c(EnsID, rowNum))
  w = v$weights %>%
    as.data.frame()
  w = w[rownames(w) %in% row_num,]
  v$weights=data.matrix(w)
  v$E = c

  return(v)
}



combinedGSEA <- function(v, idx, design, contrasts){
  # Performs the following gene set enrichment testing methods: fry, fgsea, camera,
  # and combines the raw p-values obtained with each method to obtain an overall
  # "ensemble" p-value.
  #
  # Args:
  #   v: voom object obtained from using limma::voom() or limma::voomWithQualityWeights().
  #   idx: Named list of gene sets, with numeric values corresponding to rows in v.
  #        Can be prepared using limma::ids2indices().
  #   design: design matrix prepared during limma workflow.
  #   contrasts: contrasts matrix prepared during limma workflow.
  #
  # Returns:
  #   A nested list of the results for each contrast. The indivTest
  #   slot gives the results returned by each gene set enrichment testing method.
  #   The combTest slot gives the results from combining p-values from the individual
  #   methods, along with FDR and Bonferroni-adjusted p-values.


  # Unlike fry and camera, fgsea needs a ranked list of genes, which we will
  # prepare from t-statistics obtained from fitting a linear model.
  fit <- lmFit(v, design) %>%
    contrasts.fit(contrasts) %>%
    eBayes(robust = TRUE)

  # fgsea also requires the index to be gene names rather than indices like
  # fry and camera.
  fgsea_idx <- lapply(idx, function(x){
    v %>% magrittr::extract(x,) %>% rownames()
  })


  # Run the gene set enrichment testing methods for each contrast.
  res <- colnames(contrasts) %>%
    lapply(function(x){

      # For each contrast, we need to extract out a ranked list of genes for
      # fgsea, ranked by t-statistic for DE in that particular comparison.
      fgsea_ranks <- topTable(fit, coef = x, number = Inf) %>%
        rownames_to_column("entrez") %>%
        dplyr::arrange(t)
      fgsea_ranks <- fgsea_ranks$t %>% set_names(fgsea_ranks$entrez)

      # The methods used are limma::fry, limma::camera, and fgsea.
      res2 <- list(
        fry = fry(
          v,
          idx,
          design = design,
          contrast = contrasts[, x]
        ) %>%
          rownames_to_column("Geneset") %>%
          dplyr::mutate(pval = PValue.Mixed),
        mroast = mroast(
          v,
          idx,
          design = design,
          contrast = contrasts[, x]
        ) %>%
          rownames_to_column("Geneset") %>%
          dplyr::mutate(pval = PValue.Mixed),
        camera = camera(
          v,
          idx,
          design = design,
          contrast = contrasts[, x],
          allow.neg.cor = TRUE,
          inter.gene.cor = NA
        ) %>%
          rownames_to_column("Geneset") %>%
          dplyr::mutate(pval = PValue),
        fgsea = fgsea(
          fgsea_idx,
          fgsea_ranks,
          nperm = 1e5,
          nproc = 6
        ) %>%
          dplyr::rename(Geneset = pathway)
      )

    }) %>% set_names(colnames(contrasts))

  # Extract out only the relevant Geneset and p-value columns so that
  # the rows of each data.frame can be bound together.
  res3 <- res %>% lapply(function(x){
    x %<>% lapply(function(y){
      y %<>% dplyr::select(Geneset, pval)
      return(y)
    })
    return(x)
  }) %>%
    magrittr::extract(. != "mroast") %>%
    lapply(function(x){
      x %>% do.call("rbind", .) %>%
        dplyr::group_by(Geneset) %>%
        dplyr::summarise(wilkinsonp = metap::wilkinsonp(pval, r = 1)$p) %>%
        dplyr::mutate(fdr = p.adjust(wilkinsonp, method = "fdr"),
                      bonferroni = p.adjust(wilkinsonp, method = "bonferroni"))%>%
        dplyr::arrange(wilkinsonp)
    })

  # Perform adjustment for multiple testing
  res4 <- res3 %>% bind_rows(.id = "id") %>%
    mutate(fdr = p.adjust(wilkinsonp, "BH"),
           bonferroni = p.adjust(wilkinsonp, "holm")) %>%
    split(f=.$id)

  # Return relevant results
  results <- list(
    indivTest = res,
    combTest = res4
  )

  return(results)
}


run_fgsea = function(v, idx, design, contrasts){

  fgsea_idx = idx %>%
    dplyr::rename(ensemblID = ensembl_gene) %>%
    #dplyr::inner_join(ens2Entrez) %>%
    dplyr::distinct(gs_name, ensemblID, .keep_all = TRUE) %>%
    split(f = .$gs_name) %>%
    lapply(extract2, "ensemblID")

  fit <- lmFit(v, design) %>%
    contrasts.fit(contrasts) %>%
    eBayes(robust = TRUE)

  res <- colnames(contrasts) %>%
    lapply(function(x){

      # For each contrast, we need to extract out a ranked list of genes for
      # fgsea, ranked by t-statistic for DE in that particular comparison.
      fgsea_ranks <- topTable(fit, coef = x, number = Inf) %>%
        rownames_to_column("ensembID") %>%
        mutate(rankstat = sign(logFC)*-log10(P.Value)) %>%
        arrange(desc(rankstat))

      fgsea_ranks <-  fgsea_ranks$rankstat %>% set_names(fgsea_ranks$ensembID)

      # The methods used are limma::fry, limma::camera, and fgsea.
      res2 <- list(fgsea = fgsea(
        fgsea_idx,
        fgsea_ranks,
        nperm = 1e5,
        nproc = 6
      ) %>%
        dplyr::rename(Geneset = pathway)
      )}) %>% set_names(colnames(contrasts))
  return(res)
}


### Transcript expression analysis



ma_plot = function(res, lfc_thres){
  DEColours <- c("-1" = "#2F124B","1" = "#EA2A5F",  "0" = "#D3D5E3")
  temp = res
  temp$diffexpressed <- 0
  temp$diffexpressed[temp$logFC > lfc_thres & temp$FDR < 0.05] <- "1"
  temp$diffexpressed[temp$logFC < lfc_thres & temp$FDR < 0.05] <- "-1"
  plot = temp %>% as.data.frame() %>%
    ggplot(aes(y= logFC,
               x = logCPM,
               colour = as.factor(as.character(diffexpressed)))) +
    geom_point(alpha = 0.5, size =1) +
    scale_colour_manual(values = DEColours) + theme_classic() +
    theme(legend.position = "none",
          axis.title.y = element_text(size = 12))
  return(plot)

}


