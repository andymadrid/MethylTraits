.karyotypeRefsDefault <- data.frame(
    karyotype = c("XX", "XY", "X0", "XXY", "XXX", "XYY"),
    sex       = c("F",  "M",  "F",  "M",   "F",   "M"),
    Xcopy     = c(2,    1,    1,    2,     3,     1),
    Ycopy     = c(0,    1,    0,    1,     0,     2),
    stringsAsFactors = FALSE
)

.parRegions <- function(genomeBuild) {
    if (genomeBuild == "hg38") {
        list(par1 = c(10001, 2781479),
             par2X = c(155701383, 156030895),
             par2Y = c(56887903, 57217415))
    } else {
        list(par1 = c(60001, 2699520),
             par2X = c(154931044, 155260560),
             par2Y = c(59034050, 59363566))
    }
}

#' Predict genetic sex (and sex-chromosome aneuploidy) from a bsseq object
#'
#' WGBS read coverage at CpGs is, in effect, free copy-number data: for a diploid
#' autosome the number of reads piling up per CpG is proportional to 2 copies of
#' the genome; on chrX/chrY it is proportional to however many copies of X or Y
#' the sample actually carries. This is the same read-depth logic used by
#' whole-genome sex/relatedness QC tools (e.g. Peddy, somalier, XYalign) -- it does
#' not rely on training a classifier on labeled data, and does not use methylation
#' level at all (X-inactivation makes methylation-based sex-typing noisier and
#' harder to generalize across tissues than a simple depth ratio). Per sample, the
#' function computes the median CpG coverage on the autosomes, on chrX (excluding
#' the pseudoautosomal regions, which are diploid regardless of karyotype and would
#' dilute the X-specific signal), and on chrY (excluding PAR), converts each to a
#' copy-number estimate relative to the diploid autosomal baseline (copy number =
#' 2 x chrCoverage / autosomalCoverage), and assigns each sample to the sex
#' karyotype (XX, XY, X0/Turner, XXY/Klinefelter, XXX/triple X, XYY) whose
#' theoretical (Xcopy, Ycopy) point is nearest in Euclidean distance. Because true
#' aneuploidies of the sex chromosomes are rare, most samples in most cohorts will
#' land squarely on XX or XY; the aneuploidy calls are offered as a bonus for the
#' rare sample where they are informative, flagged as such rather than silently
#' trusted, since the reference centroids below are never fit to more than a
#' handful (if any) of real aneuploid samples.
#'
#' @param bs A bsseq object.
#' @param genomeBuild One of "hg38" (default) or "hg19". Selects which pseudoautosomal
#'   region (PAR) coordinates are used to exclude PAR CpGs from the chrX/chrY signal.
#'   Must match the genome build \code{bs} was aligned to.
#' @param statFUN One of "median" (default) or "mean". Per-chromosome-group summary
#'   statistic used to combine CpG-level coverage into one depth value per sample.
#'   Median is the default since it is robust to the handful of extreme-coverage CpGs
#'   (e.g. repetitive/low-mappability loci) that any large WGBS panel contains.
#' @param karyotypeRefs Optional data.frame with columns \code{karyotype}, \code{sex}
#'   ("M"/"F"), \code{Xcopy}, \code{Ycopy} giving the reference centroids to classify
#'   against. Default NULL uses the bundled 6-karyotype table (XX, XY, X0, XXY, XXX, XYY)
#'   at their theoretical copy numbers. Supply your own if you want to add/remove
#'   karyotypes (e.g. XXYY, XXXY) or use empirically-derived centroids directly.
#' @param knownSex Optional character vector, one entry per sample in \code{bs}
#'   (same order as \code{sampleNames(bs)}/columns of \code{bs}), giving "M"/"F"
#'   (NA where unknown) for samples with independently known sex. If supplied and
#'   \code{calibrate = TRUE}, used to rescale the reference centroids to this assay's
#'   own XX/XY depth ratios before classification (see \code{calibrate}); either way,
#'   it is also used to report predicted-vs-known concordance, which is a useful way
#'   to catch sample-swap/mislabeling QC problems.
#' @param calibrate Logical. If TRUE (default) and \code{knownSex} has at least
#'   \code{minKnownPerGroup} labeled M and F samples, fits a 2-point linear rescaling
#'   per axis (using the theoretical XX/XY centroids as anchors) so the reference
#'   table matches this assay's empirically observed XX/XY depth ratios, then applies
#'   that same rescaling to all karyotypes in the table. Real WGBS coverage tends to
#'   run several percent high on chrX relative to the naive 2x/1x theoretical
#'   expectation (higher CpG-island density on chrX than the autosomal average), so
#'   calibration meaningfully sharpens classification when labels are available; with
#'   no (or too few) labels, the uncalibrated theoretical centroids are used instead.
#' @param minKnownPerGroup Minimum number of labeled M and labeled F samples (each)
#'   required in \code{knownSex} for calibration to run. Default 3.
#' @param ambiguousGapThresh Numeric. A sample is flagged \code{ambiguous} if the gap
#'   between its distance to the nearest and second-nearest reference centroid is
#'   smaller than this (in copy-number units), i.e. it sits roughly between two
#'   karyotypes rather than clearly at one -- the signature of true mosaicism, or of
#'   a low-coverage sample too noisy to call confidently. Default 0.3.
#' @param outlierDistThresh Numeric. A sample is flagged \code{outlier} if its distance
#'   to even the nearest reference centroid exceeds this -- it does not look like any
#'   of the karyotypes in \code{karyotypeRefs} (e.g. a rarer aneuploidy such as XXYY,
#'   sample contamination, or a technical artifact) and warrants manual review.
#'   Default 0.75 (about half the spacing between adjacent reference centroids).
#' @param verbose Logical. Print progress/summary messages. Default TRUE.
#' @return A data.frame, one row per sample, with columns: \code{sample},
#'   \code{medAutoCov}, \code{medXcov}, \code{medYcov} (raw depth summaries),
#'   \code{Xcopy}, \code{Ycopy} (copy-number estimates), \code{karyotype},
#'   \code{predictedSex} ("M"/"F"), \code{distToNearest}, \code{distToSecondNearest},
#'   \code{ambiguous}, \code{outlier}. If chrY has no usable (non-PAR) CpGs in \code{bs}
#'   (e.g. a Y-masked reference), \code{Ycopy} is NA throughout and karyotype/aneuploidy
#'   calls are skipped in favor of a lower-confidence X-copy-only sex call, with a warning.
#' @examples
#' \dontrun{
#' predictSex(bs)
#' predictSex(bs, knownSex = pData(bs)$Sex, verbose = TRUE)
#' }
#' @import bsseq
#' @import GenomicRanges
#' @export

predictSex <- function(bs, genomeBuild = c("hg38", "hg19"), statFUN = c("median", "mean"),
                        karyotypeRefs = NULL, knownSex = NULL, calibrate = TRUE,
                        minKnownPerGroup = 3, ambiguousGapThresh = 0.3, outlierDistThresh = 0.75,
                        verbose = TRUE) {

    if (!is(bs, "BSseq")) {
        stop("bs must be a BSseq object")
    }
    genomeBuild <- match.arg(genomeBuild)
    statFUN <- match.arg(statFUN)
    sumFUN <- if (statFUN == "median") matrixStats::colMedians else matrixStats::colMeans2

    if (is.null(karyotypeRefs)) {
        karyotypeRefs <- .karyotypeRefsDefault
    }
    if (!all(c("karyotype", "sex", "Xcopy", "Ycopy") %in% colnames(karyotypeRefs))) {
        stop("karyotypeRefs must have columns: karyotype, sex, Xcopy, Ycopy")
    }

    nSamples <- ncol(bs)
    samples <- sampleNames(bs)
    if (is.null(samples)) samples <- colnames(bs)
    if (is.null(samples)) samples <- paste0("Sample", seq_len(nSamples))

    if (!is.null(knownSex) && length(knownSex) != nSamples) {
        stop("knownSex must have one entry per sample in bs (length(knownSex) == ncol(bs))")
    }

    par <- .parRegions(genomeBuild)
    chrs <- as.character(runValue(seqnames(bs)))

    message(paste0("[predictSex]: computing autosomal baseline coverage (", statFUN, ") for ", nSamples, " samples"))
    autosomes <- paste0("chr", 1:22)
    bs.auto <- bs[as.character(seqnames(bs)) %in% autosomes, ]
    medAutoCov <- sumFUN(as.matrix(getCoverage(bs.auto, type = "Cov")))
    rm(bs.auto)

    message("[predictSex]: computing chrX (non-PAR) coverage")
    hasX <- "chrX" %in% chrs
    if (hasX) {
        bs.X <- bs[seqnames(bs) == "chrX", ]
        posX <- start(bs.X)
        keepX <- !((posX >= par$par1[1] & posX <= par$par1[2]) |
                    (posX >= par$par2X[1] & posX <= par$par2X[2]))
        bs.X <- bs.X[keepX, ]
        medXcov <- if (nrow(bs.X) > 0) sumFUN(as.matrix(getCoverage(bs.X, type = "Cov"))) else rep(NA_real_, nSamples)
        rm(bs.X)
    } else {
        warning("[predictSex]: bs has no chrX -- cannot predict sex")
        medXcov <- rep(NA_real_, nSamples)
    }

    message("[predictSex]: computing chrY (non-PAR) coverage")
    hasY <- "chrY" %in% chrs
    if (hasY) {
        bs.Y <- bs[seqnames(bs) == "chrY", ]
        posY <- start(bs.Y)
        keepY <- !((posY >= par$par1[1] & posY <= par$par1[2]) |
                    (posY >= par$par2Y[1] & posY <= par$par2Y[2]))
        bs.Y <- bs.Y[keepY, ]
        medYcov <- if (nrow(bs.Y) > 0) sumFUN(as.matrix(getCoverage(bs.Y, type = "Cov"))) else rep(NA_real_, nSamples)
        rm(bs.Y)
    } else {
        warning("[predictSex]: bs has no chrY -- karyotype/aneuploidy calls will be skipped, ",
                "falling back to a lower-confidence chrX-only sex call")
        medYcov <- rep(NA_real_, nSamples)
    }

    Xcopy <- 2 * medXcov / medAutoCov
    Ycopy <- 2 * medYcov / medAutoCov

    # ---- optional empirical calibration against knownSex ----
    calibApplied <- FALSE
    refs <- karyotypeRefs
    if (!is.null(knownSex) && isTRUE(calibrate)) {
        isM <- !is.na(knownSex) & knownSex == "M"
        isF <- !is.na(knownSex) & knownSex == "F"
        if (sum(isM) >= minKnownPerGroup && sum(isF) >= minKnownPerGroup) {
            obsXX_X <- stats::median(Xcopy[isF], na.rm = TRUE)
            obsXY_X <- stats::median(Xcopy[isM], na.rm = TRUE)
            obsXX_Y <- stats::median(Ycopy[isF], na.rm = TRUE)
            obsXY_Y <- stats::median(Ycopy[isM], na.rm = TRUE)
            slopeX <- (obsXX_X - obsXY_X) / (2 - 1)
            interceptX <- obsXY_X - slopeX * 1
            slopeY <- (obsXY_Y - obsXX_Y) / (1 - 0)
            interceptY <- obsXX_Y
            refs$Xcopy <- slopeX * karyotypeRefs$Xcopy + interceptX
            refs$Ycopy <- slopeY * karyotypeRefs$Ycopy + interceptY
            calibApplied <- TRUE
            message(paste0("[predictSex]: calibrated reference centroids from ", sum(isF), " labeled F and ",
                            sum(isM), " labeled M samples (observed XX: Xcopy=", round(obsXX_X, 2),
                            ", Ycopy=", round(obsXX_Y, 2), "; observed XY: Xcopy=", round(obsXY_X, 2),
                            ", Ycopy=", round(obsXY_Y, 2), ")"))
        } else {
            message(paste0("[predictSex]: knownSex has fewer than minKnownPerGroup=", minKnownPerGroup,
                            " labeled samples in one or both groups -- using uncalibrated theoretical centroids"))
        }
    }

    message("[predictSex]: classifying samples against reference karyotypes")
    if (hasY) {
        distMat <- sapply(seq_len(nrow(refs)), function(i) {
            sqrt((Xcopy - refs$Xcopy[i])^2 + (Ycopy - refs$Ycopy[i])^2)
        })
        ord <- t(apply(distMat, 1, order))
        nearestIdx <- ord[, 1]
        secondIdx <- ord[, 2]
        distToNearest <- distMat[cbind(seq_len(nSamples), nearestIdx)]
        distToSecondNearest <- distMat[cbind(seq_len(nSamples), secondIdx)]
        karyotype <- refs$karyotype[nearestIdx]
        predictedSex <- refs$sex[nearestIdx]
        ambiguous <- (distToSecondNearest - distToNearest) < ambiguousGapThresh
        outlier <- distToNearest > outlierDistThresh
    } else {
        # chrY unavailable: fall back to a simple midpoint threshold on Xcopy alone
        # (cannot separate XY from X0, or XX from XXY, without Y coverage)
        threshX <- mean(c(1, 2))
        predictedSex <- ifelse(Xcopy < threshX, "M", "F")
        karyotype <- ifelse(Xcopy < threshX, "XY", "XX")
        distToNearest <- abs(Xcopy - ifelse(Xcopy < threshX, 1, 2))
        distToSecondNearest <- rep(NA_real_, nSamples)
        ambiguous <- distToNearest > (threshX - 1) * 0.7
        outlier <- rep(FALSE, nSamples)
    }
    predictedSex[is.na(Xcopy)] <- NA_character_
    karyotype[is.na(Xcopy)] <- NA_character_

    result <- data.frame(
        sample = samples,
        medAutoCov = medAutoCov,
        medXcov = medXcov,
        medYcov = medYcov,
        Xcopy = Xcopy,
        Ycopy = Ycopy,
        karyotype = karyotype,
        predictedSex = predictedSex,
        distToNearest = distToNearest,
        distToSecondNearest = distToSecondNearest,
        ambiguous = ambiguous,
        outlier = outlier,
        stringsAsFactors = FALSE
    )
    attr(result, "calibrated") <- calibApplied
    attr(result, "referenceCentroids") <- refs

    if (verbose) {
        message("[predictSex]: predicted sex breakdown:")
        print(table(result$predictedSex, useNA = "ifany"))
        if (hasY) {
            aneuploid <- result$karyotype %in% c("X0", "XXY", "XXX", "XYY")
            if (any(aneuploid, na.rm = TRUE)) {
                message(paste0("[predictSex]: ", sum(aneuploid, na.rm = TRUE),
                                " sample(s) called as a sex-chromosome aneuploidy karyotype -- ",
                                "these are rare, so review Xcopy/Ycopy and ambiguous/outlier flags before trusting them:"))
                print(result[which(aneuploid), c("sample", "karyotype", "Xcopy", "Ycopy", "ambiguous", "outlier")])
            }
        }
        nFlagged <- sum(result$ambiguous | result$outlier, na.rm = TRUE)
        if (nFlagged > 0) {
            message(paste0("[predictSex]: ", nFlagged, " sample(s) flagged ambiguous and/or outlier -- recommend manual review"))
        }
    }

    if (!is.null(knownSex)) {
        known <- !is.na(knownSex)
        if (any(known)) {
            match <- result$predictedSex[known] == knownSex[known]
            acc <- mean(match, na.rm = TRUE)
            message(paste0("[predictSex]: concordance with knownSex: ", sum(match, na.rm = TRUE), "/", sum(known),
                            " (", round(100 * acc, 1), "%)"))
            if (any(!match, na.rm = TRUE)) {
                mismatches <- result[known, ][!match, c("sample", "karyotype", "predictedSex", "Xcopy", "Ycopy")]
                mismatches$knownSex <- knownSex[known][!match]
                message("[predictSex]: samples where predictedSex != knownSex (possible sample swap/mislabeling -- worth investigating):")
                print(mismatches)
            }
        }
    }

    return(result)
}
