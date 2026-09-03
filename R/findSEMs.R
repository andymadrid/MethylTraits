#' Calculates CpG-level statistics and identifies stochastic epigenetic mutations (SEMs) from sequencing-based data
#' @param bs A bsseq object
#' @param type Type of methylation values to be extracted from the bsseq object. Either raw or smooth values. For smoothed values, bsseq object needs to have been smoothed prior to running this option. Default is "raw".
#' @param minCov Numeric value. The minimum number of reads to overlap a CpG sites by all samples for it to be used for downstream analysis. Note, lowering the number will likely increase the number of CpGs that are tested for SEM identification but at two big costs: 1) increase the computation time as it will be required to analyze more CpGs and 2) likely introduce false positives as less coverage generally equals less certainty of the estimation of methylation at a given CpG site. Default is 10.
#' @param minSamples Numeric value. The minimum number of samples required to perform the SEM analysis. A larger sample size (e.g., >50) is generally required for an analysis such as this. Default is 50.
#' @param runStats Logical. Whether to run CpG-level statistics (and save them) or not. Default is FALSE.
#' @param nCG Numeric. The number of CpGs used to split data into chunks for easier processing. Default is 100000 CpGs per chunk.
#' @param saveIndSEMs Logical. Whether sample-specific SEMs should be saved as a long-format data.frame (cpg, sample, type) of every hyper/hypo SEM event. Note, this file may be large in size, depending on how many SEMs are identified. Default is FALSE.
#' @param saveDir Directory where outputs should be saved. Default is current working directory.
#' @param ncores Number of chunks to process in parallel via \code{parallel::mclapply}. Default is 1. Increase this on multi-core machines for a large speedup; each worker's memory footprint is roughly proportional to nCG x ncol(bs).
#' @param verbose Logical. Whether output of functions (e.g., number of chunks processed) should be verbose or not. Default is TRUE.
#' @examples
#' \dontrun{
#' findSEMs(bs, type = "raw", minCov = 10, minSamples = 50, runStats = FALSE, nCG = 100000, saveIndSEMs = FALSE, saveDir = getwd(), ncores = 1, verbose = TRUE)
#' }
#' @import bsseq
#' @export

findSEMs <- function(bs, type = c("raw", "smooth"), minCov = 10, minSamples = 50, runStats = FALSE, nCG = 100000, saveIndSEMs = FALSE, saveDir = getwd(), ncores = 1, verbose = TRUE) {

    type <- match.arg(type)

    # check if sample size is large enough
    if (ncol(bs) < minSamples) {
        stop("Number of samples is less than ", minSamples, "! Cannot calculate SEMs with such few samples.")
    }

    # filter CpGs based on coverage
    message("[findSEMs]: filtering low-coverage CpGs from further analysis")
    loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov")<minCov) == 0)
    bs <- bs[loci.idx,]

    # get matrix of methylation values
    message("[findSEMs]: getting methylation estimates from CpGs")
    meth.mat <- bsseq::getMeth(bs, type = type)
    rownames(meth.mat) <- paste0(seqnames(bs), ":", start(bs))
    # some BSseq objects carry no sample names at all (colnames(meth.mat) is then NULL);
    # guarantee a usable identifier so downstream chunking/saveIndSEMs never breaks on this
    if (is.null(colnames(meth.mat))) {
        colnames(meth.mat) <- paste0("Sample", seq_len(ncol(meth.mat)))
    }

    nCpG <- nrow(meth.mat)
    nSamp <- ncol(meth.mat)
    sampNames <- colnames(meth.mat)
    chunks <- unname(split(seq_len(nCpG), ceiling(seq_len(nCpG) / nCG)))
    nChunks <- length(chunks)

    if (runStats == TRUE) {

        message("[findSEMs]: getting CpG-level statistics now")
        message("[findSEMs]: partitioning data for easier processing")
        message(paste0("[findSEMs]: splitting methylation data into chunks with up to ", nCG, " CpGs each"))
        message(paste0("[findSEMs]: a total of ", nChunks, " chunks will be processed"))

        # vectorized replacement for the old per-CpG apply(cgStats): computes all
        # row-wise moments/quantiles for a whole chunk in a handful of C-level calls
        # instead of one R function call per CpG.
        statChunk <- function(idx) {
            m <- meth.mat[idx, , drop = FALSE]
            q <- matrixStats::rowQuantiles(m, probs = c(.25, .75), na.rm = TRUE)
            mu <- matrixStats::rowMeans2(m, na.rm = TRUE)
            sdv <- matrixStats::rowSds(m, na.rm = TRUE)
            mn <- matrixStats::rowMins(m, na.rm = TRUE)
            mx <- matrixStats::rowMaxs(m, na.rm = TRUE)
            centered <- m - mu
            m2 <- matrixStats::rowMeans2(centered^2, na.rm = TRUE)
            m3 <- matrixStats::rowMeans2(centered^3, na.rm = TRUE)
            m4 <- matrixStats::rowMeans2(centered^4, na.rm = TRUE)
            data.frame(
                q1 = q[, 1], q3 = q[, 2],
                skewness = m3 / m2^1.5,
                sd = sdv, min = mn, mean = mu, max = mx,
                kurtosis = m4 / m2^2,
                iqr = q[, 2] - q[, 1],
                range = mx - mn,
                row.names = rownames(m)
            )
        }

        if (ncores > 1) {
            statList <- parallel::mclapply(chunks, statChunk, mc.cores = ncores)
        } else {
            statList <- vector("list", nChunks)
            for (ii in seq_len(nChunks)) {
                statList[[ii]] <- statChunk(chunks[[ii]])
                if (verbose == TRUE) {
                    message(paste0("\tDone with chunk ", ii," of ", nChunks))
                }
            }
        }
        cpgStats <- do.call(rbind, statList)

        message("[findSEMs]: identifying CpG clusters")
        complete <- stats::complete.cases(cpgStats)
        kk <- kmeans(scale(cpgStats[complete, ]), centers = 3)
        cpgStats$cluster <- NA_integer_
        cpgStats$cluster[complete] <- kk$cluster

        message(paste0("[findSEMs]: will save CpG-level statistics at ", saveDir, "/cpgStats.rda"))
        save(cpgStats, file = paste0(saveDir, "/cpgStats.rda"))
        rm(statList)
        gc()

    }

    # workhorse . . . time to find SEMs
    message("[findSEMs]: finding SEMs now . . .")
    message(paste0("[findSEMs]: splitting methylation data into chunks with up to ", nCG, " CpGs each"))
    message(paste0("[findSEMs]: a total of ", nChunks, " chunks will be processed"))

    # vectorized replacement for the old per-CpG apply(SEMFinder): computes the
    # Tukey-fence outlier test for every CpG in a chunk with a few matrix-wide
    # comparisons instead of returning (and accumulating!) an nSamples-by-2
    # matrix per individual CpG.
    semChunk <- function(idx) {
        m <- meth.mat[idx, , drop = FALSE]
        q <- matrixStats::rowQuantiles(m, probs = c(.25, .75), na.rm = TRUE)
        iqr <- q[, 2] - q[, 1]
        hiThresh <- q[, 2] + 3 * iqr
        loThresh <- q[, 1] - 3 * iqr
        # rows with iqr == 0 (all samples identical) can have no outliers, by definition
        hiThresh[iqr == 0] <- Inf
        loThresh[iqr == 0] <- -Inf
        hyperMat <- m > hiThresh
        hypoMat  <- m < loThresh
        hyperMat[is.na(hyperMat)] <- FALSE
        hypoMat[is.na(hypoMat)] <- FALSE

        hits <- NULL
        if (saveIndSEMs) {
            hi <- which(hyperMat, arr.ind = TRUE)
            lo <- which(hypoMat, arr.ind = TRUE)
            hits <- rbind(
                if (nrow(hi) > 0) data.frame(cpg = rownames(m)[hi[, 1]], sample = colnames(m)[hi[, 2]], type = "hyper", stringsAsFactors = FALSE) else NULL,
                if (nrow(lo) > 0) data.frame(cpg = rownames(m)[lo[, 1]], sample = colnames(m)[lo[, 2]], type = "hypo",  stringsAsFactors = FALSE) else NULL
            )
        }

        list(hyperCount = matrixStats::colSums2(hyperMat), hypoCount = matrixStats::colSums2(hypoMat), hits = hits)
    }

    if (ncores > 1) {
        results <- parallel::mclapply(seq_along(chunks), function(ii) semChunk(chunks[[ii]]), mc.cores = ncores)
    } else {
        results <- vector("list", nChunks)
        for (ii in seq_len(nChunks)) {
            results[[ii]] <- semChunk(chunks[[ii]])
            if (verbose == TRUE) {
                message(paste0("\tDone with chunk ", ii," of ", nChunks))
            }
        }
    }

    hyperTotal <- integer(nSamp)
    hypoTotal <- integer(nSamp)
    for (r in results) {
        hyperTotal <- hyperTotal + r$hyperCount
        hypoTotal <- hypoTotal + r$hypoCount
    }

    if (saveIndSEMs == TRUE) {
        message("[findSEMs]: Saving individual SEMs")
        cpgSEMs <- do.call(rbind, lapply(results, `[[`, "hits"))
        rownames(cpgSEMs) <- NULL
        save(cpgSEMs, file = paste0(saveDir, "/sampleSEMs.rda"))
    }

    # clean ‘er up
    message("[findSEMs]: bringing it all back home")
    sampleSEMs <- data.frame(Hyper_SEMs = hyperTotal, Hypo_SEMs = hypoTotal, row.names = sampNames)
    sampleSEMs$Total_SEMs <- sampleSEMs$Hyper_SEMs + sampleSEMs$Hypo_SEMs
    sampleSEMs$log10SEMs <- log10(sampleSEMs$Total_SEMs)
    return(sampleSEMs)
}
