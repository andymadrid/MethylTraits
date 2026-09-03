#' Annotate stochastic epigenetic mutations (SEMs) to their nearest genes
#'
#' Takes the long-format \code{cpgSEMs} table produced by \code{findSEMs(..., saveIndSEMs = TRUE)}
#' (one row per hyper/hypo SEM event: which CpG, which sample, which direction) and assigns each
#' unique SEM-carrying CpG to its nearest gene / genomic feature via ChIPseeker + a UCSC knownGene
#' TxDb (the same annotation machinery used by \code{annotateICRs}). SEM events are then aggregated
#' into a gene-by-sample count matrix suitable for a burden-style association test (see
#' \code{testSEMGenes}).
#'
#' @param cpgSEMs The long-format data.frame returned by \code{findSEMs(..., saveIndSEMs = TRUE)}
#'   (columns: cpg (format "chr:pos"), sample, type ("hyper"/"hypo")), or the path to the saved
#'   sampleSEMs.rda file containing it.
#' @param sampleNames Character vector of ALL sample names in the analysis (e.g.
#'   \code{sampleNames(bs)} or \code{rownames(sampleSEMs)} from \code{findSEMs}). Needed so that
#'   samples with zero SEMs in a gene are correctly recorded as 0 rather than silently dropped.
#' @param genomeBuild One of "hg38" (default) or "hg19". Must match the genome build \code{bs} was
#'   aligned to.
#' @param collapseType Logical. If TRUE (default), hyper- and hypo-SEMs are pooled into a single
#'   total SEM count per gene per sample. If FALSE, separate hyper/hypo count matrices are returned.
#' @param universeMaxGap Integer. Before computing the background gene universe, tested CpGs
#'   supplied via \code{universeCpGs} are merged into contiguous blocks (gap <= this many bp) so
#'   that annotation runs on a manageable number of genomic blocks rather than millions of
#'   individual CpGs. Default is 1000.
#' @param universeCpGs Optional GRanges (or the filtered \code{bs} object) of ALL CpGs that were
#'   tested by \code{findSEMs} (i.e. passed the coverage filter), used to build a background gene
#'   universe for later GO enrichment (\code{goSEMs}) that isn't biased toward whatever the genome
#'   happens to contain overall. If NULL, no universe is computed (\code{goSEMs} will fall back to
#'   the whole-genome background, which is less correct but much faster).
#' @param verbose Logical. Print progress messages. Default TRUE.
#' @return A list with elements: \code{geneMat} (genes x samples integer count matrix; or, if
#'   \code{collapseType = FALSE}, a list of two such matrices named \code{hyper}/\code{hypo}),
#'   \code{cpgGeneMap} (data.frame: cpg, gene, featureType, distanceToTSS for every SEM-carrying
#'   CpG), and \code{universeGenes} (character vector of background gene symbols; NULL if
#'   \code{universeCpGs} was not supplied).
#' @export
annotateSEMs <- function(cpgSEMs, sampleNames, genomeBuild = c("hg38", "hg19"),
                          collapseType = TRUE, universeMaxGap = 1000,
                          universeCpGs = NULL, verbose = TRUE) {

    genomeBuild <- match.arg(genomeBuild)

    if (is.character(cpgSEMs) && length(cpgSEMs) == 1) {
        e <- new.env()
        load(cpgSEMs, envir = e)
        cpgSEMs <- get("cpgSEMs", envir = e)
    }
    if (!all(c("cpg", "sample", "type") %in% colnames(cpgSEMs))) {
        stop("cpgSEMs must have columns cpg, sample, type (as returned by findSEMs(..., saveIndSEMs = TRUE))")
    }

    txdb <- .semTxDb(genomeBuild)
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        stop("Package org.Hs.eg.db is required but not installed")
    }

    if (verbose) message("[annotateSEMs]: annotating ", length(unique(cpgSEMs$cpg)), " SEM-carrying CpGs to nearest gene")
    uCpg <- unique(cpgSEMs$cpg)
    cpg.gr <- .parseCpgIDs(uCpg)

    peakAnno <- ChIPseeker::annotatePeak(cpg.gr, TxDb = txdb, annoDb = "org.Hs.eg.db",
                                          verbose = FALSE, sameStrand = FALSE)
    anno.df <- as.data.frame(peakAnno)

    cpgGeneMap <- data.frame(
        cpg = uCpg,
        gene = if ("SYMBOL" %in% colnames(anno.df)) anno.df$SYMBOL else NA_character_,
        featureType = sub(" \\(.*", "", anno.df$annotation),
        distanceToTSS = anno.df$distanceToTSS,
        stringsAsFactors = FALSE
    )

    cpgSEMs <- merge(cpgSEMs, cpgGeneMap[, c("cpg", "gene")], by = "cpg")
    cpgSEMs <- cpgSEMs[!is.na(cpgSEMs$gene), ]

    sampleNames <- unique(sampleNames)
    buildMat <- function(df) {
        if (nrow(df) == 0) return(matrix(0L, nrow = 0, ncol = length(sampleNames), dimnames = list(NULL, sampleNames)))
        genes <- sort(unique(df$gene))
        tab <- table(factor(df$gene, levels = genes), factor(df$sample, levels = sampleNames))
        matrix(tab, nrow = length(genes), ncol = length(sampleNames), dimnames = list(genes, sampleNames))
    }

    if (collapseType) {
        geneMat <- buildMat(cpgSEMs)
    } else {
        geneMat <- list(hyper = buildMat(cpgSEMs[cpgSEMs$type == "hyper", ]),
                         hypo  = buildMat(cpgSEMs[cpgSEMs$type == "hypo", ]))
    }

    universeGenes <- NULL
    if (!is.null(universeCpGs)) {
        if (is(universeCpGs, "BSseq")) universeCpGs <- GenomicRanges::granges(universeCpGs)
        if (verbose) message("[annotateSEMs]: building background gene universe from ", length(universeCpGs), " tested CpGs")
        blocks <- GenomicRanges::reduce(universeCpGs, min.gapwidth = universeMaxGap)
        if (verbose) message("[annotateSEMs]: merged into ", length(blocks), " blocks (gap <= ", universeMaxGap, "bp) for annotation")
        uAnno <- ChIPseeker::annotatePeak(blocks, TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE, sameStrand = FALSE)
        uAnno.df <- as.data.frame(uAnno)
        universeGenes <- unique(stats::na.omit(uAnno.df$SYMBOL))
        if (verbose) message("[annotateSEMs]: background universe contains ", length(universeGenes), " genes")
    }

    list(geneMat = geneMat, cpgGeneMap = cpgGeneMap, universeGenes = universeGenes)
}

#' Test whether genome-wide (global) SEM burden differs across groups
#'
#' A simple, assumption-light omnibus test analogous to the "Global-EML" metric: are per-sample
#' total SEM counts (from \code{findSEMs}) associated with a grouping variable (e.g. diagnosis)?
#' Uses Kruskal-Wallis by default (rank-based, robust to the right-skew typical of SEM-count
#' distributions and to unequal group sizes); optionally a covariate-adjusted Poisson/negative-
#' binomial GLM if \code{covariates} are supplied.
#'
#' @param sampleSEMs The data.frame returned by \code{findSEMs} (needs a count column, default
#'   \code{Total_SEMs}).
#' @param group Factor or character vector, one entry per row of \code{sampleSEMs}, giving group
#'   membership (e.g. \code{pData(bs)$Diagnosis}).
#' @param countCol Name of the column in \code{sampleSEMs} to test. Default "Total_SEMs".
#' @param covariates Optional data.frame of additional covariates (e.g. age, sex, cell
#'   composition), one row per sample, same order as \code{sampleSEMs}. If supplied, fits
#'   \code{MASS::glm.nb(count ~ group + covariates)} and returns a likelihood-ratio test instead
#'   of Kruskal-Wallis.
#' @param verbose Logical. Default TRUE.
#' @return A list with \code{test} (htest object or anova table), \code{method}, and a per-group
#'   summary data.frame (n, median, mean SEM burden).
#' @export
testSEMBurden <- function(sampleSEMs, group, countCol = "Total_SEMs", covariates = NULL, verbose = TRUE) {

    if (!countCol %in% colnames(sampleSEMs)) stop("countCol '", countCol, "' not found in sampleSEMs")
    y <- sampleSEMs[[countCol]]
    group <- factor(group)

    summ <- do.call(rbind, lapply(split(y, group), function(v) {
        data.frame(n = length(v), median = stats::median(v), mean = mean(v), sd = stats::sd(v))
    }))
    summ <- cbind(group = rownames(summ), summ)
    rownames(summ) <- NULL

    if (is.null(covariates)) {
        if (verbose) message("[testSEMBurden]: running Kruskal-Wallis test of ", countCol, " ~ group")
        test <- stats::kruskal.test(y ~ group)
        method <- "kruskal.test"
    } else {
        if (!requireNamespace("MASS", quietly = TRUE)) stop("Package MASS is required for covariate-adjusted testing but not installed")
        if (verbose) message("[testSEMBurden]: running covariate-adjusted negative-binomial GLM of ", countCol, " ~ group + covariates")
        df <- cbind(y = y, group = group, covariates)
        full <- MASS::glm.nb(y ~ ., data = df)
        null <- MASS::glm.nb(y ~ . - group, data = df)
        test <- stats::anova(null, full, test = "Chisq")
        method <- "glm.nb LRT"
    }

    list(test = test, method = method, summary = summ)
}

#' Test gene-level SEM burden for association with a grouping variable
#'
#' Treats the gene x sample SEM count matrix (from \code{annotateSEMs}) exactly like an RNA-seq
#' count matrix and runs a standard edgeR quasi-likelihood negative-binomial pipeline: this is a
#' simpler, fully-in-R alternative to the SKAT-O/RVTESTS rare-variant burden test used for
#' "Gene-EML" association testing, while still properly modeling the sparse, over-dispersed count
#' data (most genes have 0-2 SEMs per sample) and supporting >2 groups and covariates natively.
#'
#' @param geneMat Genes x samples integer count matrix (the \code{geneMat} element of
#'   \code{annotateSEMs}'s return value).
#' @param group Factor or character vector, one entry per column of \code{geneMat}, giving group
#'   membership. Must be in the same sample order as \code{colnames(geneMat)}.
#' @param covariates Optional data.frame of additional covariates (e.g. age, sex), one row per
#'   sample in the same order as \code{colnames(geneMat)}, added to the design matrix alongside
#'   \code{group}.
#' @param minSamplesWithSEM Integer. Genes with fewer than this many samples carrying >=1 SEM are
#'   dropped before testing (too sparse to say anything about). Default 3.
#' @param verbose Logical. Default TRUE.
#' @return A list with \code{results} (data.frame, one row per tested gene, sorted by FDR: gene,
#'   logCPM, F, PValue, FDR, plus per-group mean counts), \code{fit} (the edgeR glmQLFit object,
#'   for custom contrasts), and \code{dge} (the edgeR DGEList).
#' @export
testSEMGenes <- function(geneMat, group, covariates = NULL, minSamplesWithSEM = 3, verbose = TRUE) {

    if (!requireNamespace("edgeR", quietly = TRUE)) stop("Package edgeR is required but not installed")

    group <- factor(group)
    if (length(group) != ncol(geneMat)) stop("length(group) must equal ncol(geneMat)")
    if (!is.null(covariates) && nrow(covariates) != ncol(geneMat)) stop("nrow(covariates) must equal ncol(geneMat)")

    keep <- rowSums(geneMat > 0) >= minSamplesWithSEM
    if (verbose) message("[testSEMGenes]: testing ", sum(keep), " of ", nrow(geneMat), " genes with SEMs in >= ", minSamplesWithSEM, " samples")
    geneMat <- geneMat[keep, , drop = FALSE]
    if (nrow(geneMat) == 0) stop("No genes pass minSamplesWithSEM -- nothing to test")

    designDF <- if (is.null(covariates)) data.frame(group = group) else data.frame(group = group, covariates)
    design <- stats::model.matrix(~ ., data = designDF)

    dge <- edgeR::DGEList(counts = geneMat, group = group)
    dge <- edgeR::calcNormFactors(dge)
    dge <- edgeR::estimateDisp(dge, design)
    fit <- edgeR::glmQLFit(dge, design)

    groupCols <- grep("^group", colnames(design))
    if (verbose) message("[testSEMGenes]: running omnibus QL F-test across ", nlevels(group), " groups")
    qlf <- edgeR::glmQLFTest(fit, coef = groupCols)

    tt <- edgeR::topTags(qlf, n = Inf)$table
    tt$gene <- rownames(tt)

    cpm <- edgeR::cpm(dge)
    groupMeans <- sapply(levels(group), function(g) rowMeans(geneMat[rownames(tt), group == g, drop = FALSE]))
    colnames(groupMeans) <- paste0("meanSEM_", levels(group))
    results <- cbind(tt, groupMeans[tt$gene, , drop = FALSE])
    rownames(results) <- NULL
    results <- results[order(results$FDR), ]

    list(results = results, fit = fit, dge = dge)
}

#' Run GO / KEGG over-representation analysis on genes with significant SEM burden
#'
#' Thin, opinionated wrapper around \code{clusterProfiler::enrichGO}/\code{enrichKEGG}, using a
#' user-supplied background gene universe (strongly recommended -- see \code{annotateSEMs}'s
#' \code{universeCpGs} argument -- so enrichment reflects genes that could plausibly have been
#' tested, not the whole genome).
#'
#' @param genes Character vector of significant gene symbols (e.g.
#'   \code{testSEMGenes(...)$results$gene[testSEMGenes(...)$results$FDR < 0.05]}).
#' @param universe Optional character vector of background gene symbols (from
#'   \code{annotateSEMs(...)$universeGenes}). If NULL, clusterProfiler uses all genes in
#'   org.Hs.eg.db as background, which is faster but less appropriate for a targeted assay like
#'   WGBS-derived SEMs.
#' @param ont One of "BP", "MF", "CC", or "ALL" (default) for \code{enrichGO}.
#' @param runKEGG Logical. Also run \code{enrichKEGG}. Default TRUE.
#' @param pvalueCutoff,qvalueCutoff Significance thresholds passed to clusterProfiler. Defaults
#'   0.05 and 0.2 (clusterProfiler's own defaults).
#' @param verbose Logical. Default TRUE.
#' @return A list with \code{GO} (enrichResult) and, if \code{runKEGG}, \code{KEGG} (enrichResult).
#' @export
goSEMs <- function(genes, universe = NULL, ont = c("ALL", "BP", "MF", "CC"), runKEGG = TRUE,
                    pvalueCutoff = 0.05, qvalueCutoff = 0.2, verbose = TRUE) {

    ont <- match.arg(ont)
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("Package clusterProfiler is required but not installed")
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) stop("Package org.Hs.eg.db is required but not installed")

    genes <- unique(stats::na.omit(genes))
    if (length(genes) == 0) stop("genes is empty -- nothing to test")

    if (verbose) message("[goSEMs]: running GO (", ont, ") over-representation on ", length(genes), " genes",
                          if (!is.null(universe)) paste0(" against a ", length(universe), "-gene background universe") else " against the whole-genome background")
    goRes <- clusterProfiler::enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                        keyType = "SYMBOL", ont = ont, pAdjustMethod = "BH",
                                        pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)

    out <- list(GO = goRes)

    if (runKEGG) {
        if (verbose) message("[goSEMs]: mapping SYMBOL -> ENTREZID for KEGG")
        idMap <- clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db)
        uniMap <- if (!is.null(universe)) {
            clusterProfiler::bitr(universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db)$ENTREZID
        } else NULL
        if (verbose) message("[goSEMs]: running KEGG over-representation on ", nrow(idMap), " mapped genes")
        keggRes <- clusterProfiler::enrichKEGG(gene = idMap$ENTREZID, universe = uniMap, organism = "hsa",
                                                pAdjustMethod = "BH", pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)
        out$KEGG <- keggRes
    }

    out
}

# ---- internal helpers, shared with annotateICRs.R's conventions ----

.semTxDb <- function(genomeBuild) {
    if (genomeBuild == "hg38") {
        if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
            stop("Package TxDb.Hsapiens.UCSC.hg38.knownGene is required but not installed")
        }
        TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    } else {
        if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE)) {
            stop("Package TxDb.Hsapiens.UCSC.hg19.knownGene is required but not installed")
        }
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    }
}

.parseCpgIDs <- function(cpgIDs) {
    parts <- strsplit(cpgIDs, ":", fixed = TRUE)
    chr <- vapply(parts, `[`, character(1), 1)
    pos <- as.integer(vapply(parts, `[`, character(1), 2))
    gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos, width = 1))
    names(gr) <- cpgIDs
    gr
}
