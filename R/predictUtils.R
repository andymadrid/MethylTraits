#' Internal helpers shared by predictRace(), predictPASC(), and predictAge()
#'
#' Not exported. Two small pieces of logic are common to more than one of this
#' package's methylation-based predictors, so they live here once rather than
#' being copy-pasted into each predictor's file: (1) resolving the directory that
#' holds a given predictor's bundled model artifacts, and (2) auto-detecting
#' whether an input bsseq object's coordinates are 0-based or 1-based by checking
#' which convention reads "CG" at the reported position in the hg38 reference
#' sequence (predictAge() has its own, GRanges-overlap-based, offset-detection
#' logic instead -- see predict_age.R's original design note -- since it is
#' matching against a whole reference GRanges panel rather than a coordinate
#' table read from CSV).
#' @keywords internal
#' @noRd

.mtModelDir <- function(sub) {
    pkgDir <- tryCatch(system.file("models", sub, package = "MethylTraits"), error = function(e) "")
    if (nzchar(pkgDir) && dir.exists(pkgDir)) return(pkgDir)
    # dev-mode fallback: package not installed, source tree is being used directly
    file.path("/media/data/WGBS/Datasets/misc/oliveR2/inst/models", sub)
}

.autoDetectCGOffset <- function(chrs, pos, verbose = TRUE) {
    if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE) ||
        !requireNamespace("Biostrings", quietly = TRUE)) {
        if (verbose) message("[autoDetectCGOffset]: BSgenome.Hsapiens.UCSC.hg38 not available -- ",
                              "defaulting to offset=0 (already 1-based). Pass offset explicitly to override.")
        return(0L)
    }
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    set.seed(1)
    n_sample <- min(2000, length(pos))
    idx <- sample(length(pos), n_sample)
    s_chrs <- chrs[idx]; s_pos <- pos[idx]
    seqlen <- GenomeInfoDb::seqlengths(genome)[s_chrs]
    valid <- s_chrs %in% GenomeInfoDb::seqnames(genome) & !is.na(seqlen) & s_pos >= 2 & (s_pos + 2) <= seqlen
    s_chrs <- s_chrs[valid]; s_pos <- s_pos[valid]

    get_dint <- function(p) {
        tryCatch(as.character(Biostrings::getSeq(genome, s_chrs, start = p, width = 2)),
                 error = function(e) rep(NA_character_, length(p)))
    }
    match0 <- if (length(s_pos) > 0) mean(get_dint(s_pos) == "CG", na.rm = TRUE) else NaN
    match1 <- if (length(s_pos) > 0) mean(get_dint(s_pos + 1) == "CG", na.rm = TRUE) else NaN
    if (verbose) message(sprintf("[autoDetectCGOffset]: CG-match rate offset=0: %.1f%%, offset=1: %.1f%%",
                                  100 * match0, 100 * match1))
    if (max(match0, match1, na.rm = TRUE) < 0.5 || abs(match0 - match1) < 0.3) {
        if (verbose) message("[autoDetectCGOffset]: ambiguous -- defaulting to offset=0")
        0L
    } else if (match1 > match0) {
        1L
    } else {
        0L
    }
}
