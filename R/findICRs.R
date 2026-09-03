#' Identify imprinting control regions (ICRs) from a bsseq object
#'
#' Scans a BSseq object for genomic windows of `numCG` consecutive CpGs,
#' spanning no more than `maxWidth` bp, whose coverage-weighted methylation
#' level is consistently hemimethylated (somatic tissue) or consistently
#' monoallelic (gamete tissue) across a required fraction of samples. This
#' mirrors the approach used by Jima et al. 2022 (Epigenetics 17:1920-1943,
#' "Genomic map of candidate human imprint control regions: the imprintome"),
#' who used >=5 consecutive CpGs within a 300bp window, 50%+-15% methylation
#' in soma, and <10%/>90% methylation in gametes.
#'
#' @param bs A bsseq object
#' @param minCov Minimum coverage of a CpG from a sample to be considered for further analysis. Default is 4.
#' @param minSample Minimum proportion of samples that must meet minCov (at the CpG level) and must have
#'   usable read coverage in a candidate window (at the window level) for that CpG/window to be retained.
#'   Default is 0.5 (50%).
#' @param numCG Number of adjacent CpGs used in the sliding window to identify ICRs. Default is 5.
#' @param maxWidth Maximum genomic span (bp) allowed for a window of numCG CpGs. Windows spanning more
#'   than this are rejected, since real ICRs are compact, densely-CpG-rich elements, not just any numCG
#'   CpGs that happened to survive coverage filtering. Default is 300 (matching Jima et al. 2022).
#' @param tissueType Type of tissue of samples. One of "somatic" or "gamete". Default is "somatic".
#' @param lowerLimitSomatic Lower bound of methylation fraction consistent with hemimethylation in somatic
#'   tissue. Default is 0.35 (35%), i.e. 50%-15%.
#' @param upperLimitSomatic Upper bound of methylation fraction consistent with hemimethylation in somatic
#'   tissue. Default is 0.65 (65%), i.e. 50%+15%.
#' @param lowerLimitGamete Methylation fraction below which a window is considered fully unmethylated in
#'   gametic tissue. Default is 0.1 (10%).
#' @param upperLimitGamete Methylation fraction above which a window is considered fully methylated in
#'   gametic tissue. Default is 0.9 (90%).
#' @param minSampleFrac Minimum proportion of informative samples (i.e. samples with any read coverage in
#'   the window) that must satisfy the hemimethylation/monoallelic criterion for a window to be called a
#'   candidate ICR. Default is 0.75. True ICRs are expected to look hemimethylated in essentially every
#'   sample, but in a large, heterogeneous cohort (different tissues, ages, disease states, individual
#'   SNPs disrupting a CpG, sampling noise) a strict 100%-of-samples requirement is unrealistic: at 11/11
#'   well-established human ICRs (PEG3, PEG10, MEST, GRB10, KCNQ1OT1, H19, MEG3/DLK1, SNRPN, GNAS,
#'   L3MBTL1, BLCAP/NNAT) tested against a 306-sample, 5-diagnosis-group whole-blood WGBS cohort, only
#'   76-98% of samples fell in the hemimethylated band at the true ICR window; requiring 100% recovered
#'   none of them, while 0.75 recovered all of them.
#' @param ncores Number of chromosomes to process in parallel via \code{parallel::mclapply}. Default is 1.
#'   Increase this on multi-core machines for a large speedup; each worker's memory footprint is roughly
#'   proportional to (CpGs on its chromosome) x (number of samples), so scale with available RAM.
#' @param verbose Logical. Print per-chromosome progress messages. Default TRUE.
#' @return A data.frame of candidate ICRs (seqnames, start, end, width, strand), or a zero-row data.frame
#'   with those columns if none were found.
#' @examples
#' \dontrun{
#' findICRs(bs, minCov = 4, minSample = 0.5, numCG = 5, maxWidth = 300, tissueType = "somatic")
#' }
#' @import bsseq
#' @import GenomicRanges
#' @export

findICRs <- function(bs, minCov = 4, minSample = 0.5, numCG = 5, maxWidth = 300,
                      tissueType = NULL,
                      lowerLimitSomatic = 0.35, upperLimitSomatic = 0.65,
                      lowerLimitGamete = 0.1, upperLimitGamete = 0.9,
                      minSampleFrac = 0.75, ncores = 1, verbose = TRUE) {

    if (!is(bs, "BSseq")) {
        stop("bs must be a BSseq object")
    }
    if (is.null(tissueType)) {
        tissueType <- "somatic"
    }
    if ((tissueType != "somatic") & (tissueType != "gamete")) {
        stop("tissueType must be one of gamete or somatic")
    }
    if (tissueType == "somatic") {
        lowerLimit <- lowerLimitSomatic
        upperLimit <- upperLimitSomatic
    }
    if (tissueType == "gamete") {
        lowerLimit <- lowerLimitGamete
        upperLimit <- upperLimitGamete
    }
    if (numCG < 2) {
        stop("numCG must be >= 2")
    }

    message(paste0("[findICRs]: will look for ICRs from ", tissueType, " using limits of [", lowerLimit, ", ", upperLimit, "]"))
    message(paste0("[findICRs]: requiring ", numCG, " consecutive CpGs spanning <= ", maxWidth, "bp"))

    nSamples <- ncol(bs)
    chrs <- as.character(runValue(seqnames(bs)))
    seqlevels(bs) <- seqlevelsInUse(bs)

    processChr <- function(chr) {

        bs.chr <- bs[seqnames(bs) == chr, ]
        n0 <- nrow(bs.chr)
        if (n0 < numCG) return(NULL)

        Cov <- as.matrix(getCoverage(bs.chr, type = "Cov"))

        # CpG-level coverage filter: keep CpGs where enough samples meet minCov
        keep <- DelayedMatrixStats::rowSums2(Cov < minCov) <= (nSamples * (1 - minSample))
        if (!any(keep)) return(NULL)

        bs.chr <- bs.chr[keep, ]
        Cov <- Cov[keep, , drop = FALSE]
        M <- as.matrix(getCoverage(bs.chr, type = "M"))
        pos <- start(bs.chr)
        n <- nrow(Cov)
        if (n < numCG) return(NULL)

        if (verbose) message(paste0("\tWorking on ", chr, " (", n, " CpGs pass coverage filter)"))

        # rolling window sums of methylated- and total-read counts via cumulative sums;
        # this pools reads across the numCG CpGs per sample (coverage-weighted region
        # methylation), which is far more robust to per-CpG sampling noise than averaging
        # (or worse, taking the raw min/max of) individual per-CpG beta values.
        csM <- rbind(0, apply(M, 2, cumsum))
        csCov <- rbind(0, apply(Cov, 2, cumsum))
        idxHi <- (numCG + 1):(n + 1)
        idxLo <- 1:(n - numCG + 1)
        winM <- csM[idxHi, , drop = FALSE] - csM[idxLo, , drop = FALSE]
        winCov <- csCov[idxHi, , drop = FALSE] - csCov[idxLo, , drop = FALSE]
        winMeth <- winM / winCov  # NaN where a sample has zero reads across the whole window

        # genomic span of each window: reject windows where the numCG CpGs are not
        # actually clustered together (e.g. because low-coverage CpGs in between were
        # dropped by the coverage filter, making "adjacent rows" genomically distant)
        span <- pos[numCG:n] - pos[1:(n - numCG + 1)]
        spanOK <- span <= maxWidth

        if (tissueType == "somatic") {
            passMat <- (winMeth >= lowerLimit) & (winMeth <= upperLimit)
        } else {
            passMat <- (winMeth <= lowerLimit) | (winMeth >= upperLimit)
        }

        nInformative <- rowSums(!is.na(passMat))
        fracPass <- rowMeans(passMat, na.rm = TRUE)
        fracPass[nInformative == 0] <- 0

        windowOK <- spanOK & ((nInformative / nSamples) >= minSample) & (fracPass >= minSampleFrac)
        hits <- which(windowOK)
        if (length(hits) == 0) return(NULL)

        gr <- GRanges(chr, IRanges(start = pos[hits], end = pos[hits + numCG - 1]))
        seqlevels(gr) <- chrs
        gr
    }

    message("[findICRs]: identifying ICRs now")
    if (ncores > 1) {
        resList <- parallel::mclapply(chrs, processChr, mc.cores = ncores)
    } else {
        resList <- lapply(chrs, processChr)
    }
    resList <- resList[!vapply(resList, is.null, logical(1))]

    message("[findICRs]: cleaning up now")
    emptyResult <- data.frame(seqnames = character(0), start = integer(0), end = integer(0),
                               width = integer(0), strand = character(0))
    if (length(resList) == 0) {
        message("[findICRs]: no ICRs were identified. Bummer . . .")
        return(emptyResult)
    }

    candidate_icrs.gr <- do.call(c, resList)
    candidate_icrs <- as.data.frame(GenomicRanges::reduce(candidate_icrs.gr))
    if (nrow(candidate_icrs) == 0) {
        message("[findICRs]: no ICRs were identified. Bummer . . .")
        return(emptyResult)
    }
    message(paste0("[findICRs]: found ", nrow(candidate_icrs), " candidate ICRs"))
    return(candidate_icrs)
}
