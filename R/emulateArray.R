#' Emulate a 450K / EPICv1 / EPICv2 microarray from a WGBS bsseq object
#'
#' Filters the CpGs in a hg38 bsseq object down to just the CpGs profiled by the
#' chosen Illumina array, producing a beta-value matrix directly comparable to
#' real array data (e.g. for applying array-trained models, or for
#' cross-platform QC/harmonization work such as \code{predictPASC}'s training
#' pipeline, which was built this way).
#'
#' Coordinate handling: array probe positions are taken from hg38 Illumina
#' manifests. For 450K/EPICv1 the package's default annotation columns are legacy
#' hg19, so the embedded \code{CHR_hg38}/\code{Start_hg38} columns are used
#' instead (\code{Start_hg38} is 0-based, so +1 gives the 1-based CpG position);
#' EPICv2's hg38 manifest is hg38-native in its default chr/pos columns already.
#' The bsseq object's own coordinate convention (0- vs 1-based) is auto-detected
#' by comparing overlap counts against the array probes as-is vs. shifted by +1,
#' and the better-matching convention is used.
#'
#' @param bs A bsseq object (hg38).
#' @param array One of "EPICv1" (default), "450K", or "EPICv2" -- also accepts
#'   "EPIC"/"EPIC1"/"EPIC2" and "450" as synonyms.
#' @param covMin Minimum read coverage to call a beta value at a CpG; below this
#'   the value is NA. Default 5.
#' @param runClocks Logical. Also run DunedinPACE and every methylCIPHER
#'   epigenetic clock on the emulated array data. Default FALSE. Requires the
#'   DunedinPACE and methylCIPHER packages.
#' @param methylOut,ageOut Optional CSV paths to also write the emulated beta
#'   matrix and (if \code{runClocks}) the age/pace estimates to. Default NULL
#'   (return only).
#' @param verbose Logical. Print progress messages. Default TRUE.
#' @return A list with \code{beta} (CpGs x samples emulated array beta matrix,
#'   row-named by Illumina probe ID) and, if \code{runClocks = TRUE},
#'   \code{ageEstimates} (a data.frame of every methylCIPHER clock plus
#'   DunedinPACE, one row per sample, with \code{pData(bs)} appended if present).
#' @examples
#' \dontrun{
#' emulateArray(bs, array = "EPICv1")
#' emulateArray(bs, array = "EPICv2", runClocks = TRUE)
#' }
#' @import bsseq
#' @import GenomicRanges
#' @export

emulateArray <- function(bs, array = c("EPICv1", "450K", "EPICv2"), covMin = 5, runClocks = FALSE,
                          methylOut = NULL, ageOut = NULL, verbose = TRUE) {

    if (!is(bs, "BSseq")) stop("bs must be a BSseq object")
    if (is.character(array) && length(array) > 1) array <- array[1]
    array <- switch(toupper(array),
                     "450K" = "450K", "450" = "450K",
                     "EPICV1" = "EPICv1", "EPIC" = "EPICv1", "EPIC1" = "EPICv1",
                     "EPICV2" = "EPICv2", "EPIC2" = "EPICv2",
                     stop("Unknown array '", array, "'; must be one of 450K, EPICv1, EPICv2"))
    if (!requireNamespace("minfi", quietly = TRUE)) stop("Package minfi is required but not installed")
    suppressPackageStartupMessages(library(minfi))

    if (verbose) message("[emulateArray]: building ", array, " hg38 probe panel")
    if (array %in% c("450K", "EPICv1")) {
        if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b5.hg38", quietly = TRUE)) {
            stop("Package IlluminaHumanMethylationEPICanno.ilm10b5.hg38 is required but not installed")
        }
        # getAnnotation()'s underlying S4 dispatch needs the annotation package
        # attached (not just namespace-loaded) -- same class of issue as
        # coef.multnet() in predictRace(), just surfacing as an S4 "get" failure here.
        suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38))
        anno <- as.data.frame(minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38))
        if (array == "450K") anno <- anno[!is.na(anno$Methyl450_Loci) & anno$Methyl450_Loci == "TRUE", ]
        anno$Start_hg38 <- as.numeric(anno$Start_hg38)
        anno <- anno[!is.na(anno$Start_hg38) & !is.na(anno$CHR_hg38) & anno$CHR_hg38 != "", ]
        probe_gr <- GRanges(anno$CHR_hg38, IRanges(anno$Start_hg38 + 1, width = 1))
        names(probe_gr) <- rownames(anno)
    } else {
        if (!requireNamespace("IlluminaHumanMethylationEPICv2anno.20a1.hg38", quietly = TRUE)) {
            stop("Package IlluminaHumanMethylationEPICv2anno.20a1.hg38 is required but not installed")
        }
        suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
        anno <- as.data.frame(minfi::getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
        anno$pos <- as.numeric(anno$pos)
        anno <- anno[!is.na(anno$pos) & !is.na(anno$chr) & anno$chr != "", ]
        probe_gr <- GRanges(anno$chr, IRanges(anno$pos, width = 1))
        names(probe_gr) <- rownames(anno)
    }
    if (verbose) message(array, " probe panel size: ", length(probe_gr))

    bs_gr <- granges(bs)
    if (length(intersect(seqlevels(bs_gr), seqlevels(probe_gr))) == 0) {
        if (all(grepl("^chr", seqlevels(probe_gr))) && !any(grepl("^chr", seqlevels(bs_gr)))) {
            seqlevels(bs_gr) <- paste0("chr", seqlevels(bs_gr))
        } else if (!any(grepl("^chr", seqlevels(probe_gr))) && all(grepl("^chr", seqlevels(bs_gr)))) {
            seqlevels(bs_gr) <- sub("^chr", "", seqlevels(bs_gr))
        }
    }

    ov_asis <- findOverlaps(probe_gr, bs_gr, type = "equal")
    shift_gr <- GRanges(seqnames(bs_gr), IRanges(start(bs_gr) + 1, width = 1))
    ov_shift <- findOverlaps(probe_gr, shift_gr, type = "equal")
    n_asis <- length(unique(queryHits(ov_asis)))
    n_shift <- length(unique(queryHits(ov_shift)))
    if (verbose) message(sprintf("[emulateArray]: coordinate check: %d/%d probes match as-is; %d/%d match with a +1 shift",
                                  n_asis, length(probe_gr), n_shift, length(probe_gr)))

    if (n_shift > n_asis) {
        if (verbose) message("[emulateArray]: bsseq object appears 0-based; shifting +1 to align with 1-based array probe positions")
        ov <- ov_shift
        coord_mode <- "0-based (shifted +1)"
    } else {
        if (verbose) message("[emulateArray]: bsseq object appears standard 1-based; using positions as-is")
        ov <- ov_asis
        coord_mode <- "1-based"
    }

    ov_df <- as.data.frame(ov)
    ov_df <- ov_df[!duplicated(ov_df$queryHits), ]
    if (verbose) message(sprintf("[emulateArray]: matched %d / %d %s probes to CpGs in the bsseq object (%s)",
                                  nrow(ov_df), length(probe_gr), array, coord_mode))
    if (nrow(ov_df) == 0) {
        stop("No probes matched. Check that the bsseq object is hg38 and that chromosome naming ",
             "(e.g. 'chr1' vs '1') is compatible with the array manifest.")
    }

    if (verbose) message("[emulateArray]: extracting beta values (min coverage = ", covMin, ")")
    M <- as.matrix(getCoverage(bs, type = "M")[ov_df$subjectHits, , drop = FALSE])
    Cov <- as.matrix(getCoverage(bs, type = "Cov")[ov_df$subjectHits, , drop = FALSE])
    beta <- M / Cov
    beta[Cov < covMin] <- NA
    rownames(beta) <- names(probe_gr)[ov_df$queryHits]

    sample_ids <- sampleNames(bs)
    if (is.null(sample_ids) || anyNA(sample_ids) || any(sample_ids == "") || any(duplicated(sample_ids))) {
        if (verbose) message("[emulateArray]: bsseq object has no valid/unique sampleNames() -- using generic sample IDs")
        sample_ids <- sprintf("sample%d", seq_len(ncol(beta)))
    }
    colnames(beta) <- sample_ids
    beta <- beta[order(rownames(beta)), , drop = FALSE]

    if (verbose) message(sprintf("[emulateArray]: emulated %s array: %d CpGs x %d samples. Fraction NA (Cov<%d): %.2f%%",
                                  array, nrow(beta), ncol(beta), covMin, 100 * mean(is.na(beta))))
    if (!is.null(methylOut)) {
        write.csv(beta, methylOut, row.names = TRUE)
        if (verbose) message("[emulateArray]: wrote ", methylOut)
    }

    result <- list(beta = beta)

    if (runClocks) {
        if (!requireNamespace("DunedinPACE", quietly = TRUE)) stop("Package DunedinPACE is required but not installed")
        if (!requireNamespace("methylCIPHER", quietly = TRUE)) stop("Package methylCIPHER is required but not installed")
        # attached, not just namespace-loaded -- see predictRace()'s note on why
        suppressPackageStartupMessages({ library(DunedinPACE); library(methylCIPHER) })

        if (verbose) message("[emulateArray]: running DunedinPACE")
        pace_res <- DunedinPACE::PACEProjector(beta, proportionOfProbesRequired = 0.8)
        pace <- pace_res$DunedinPACE
        if (verbose) message(sprintf("[emulateArray]: DunedinPACE computed for %d/%d samples", sum(!is.na(pace)), ncol(beta)))

        if (verbose) message("[emulateArray]: running all methylCIPHER clocks (please remember to cite the ",
                              "individual clocks you use; see methylCIPHER's README.md)")
        DNAm <- as.data.frame(t(beta))
        clock_fns <- c(
            "calcHannum", "calcHorvath1", "calcHorvath2", "calcPhenoAge", "calcEpiTOC",
            "calcEpiTOC2", "calcDunedinPoAm38", "calcAlcoholMcCartney", "calcBMIMcCartney",
            "calcBocklandt", "calcBohlin", "calcDNAmClockCortical", "calcDNAmTL",
            "calcGaragnani", "calcHRSInChPhenoAge", "calcHypoClock", "calcKnight",
            "calcLeeControl", "calcLeeRefinedRobust", "calcLeeRobust", "calcLin",
            "calcMayne", "calcMiAge", "calcPEDBE", "calcSmokingMcCartney",
            "calcVidalBralo", "calcWeidner", "calcZhang", "calcZhang2019"
        )
        pheno <- data.frame(row.names = rownames(DNAm))
        for (fn in clock_fns) {
            pheno <- tryCatch(get(fn)(DNAm, pheno, imputation = FALSE),
                               error = function(e) {
                                   if (verbose) message("  [skipped] ", fn, " - ", conditionMessage(e))
                                   pheno
                               })
        }

        age_est <- pheno
        age_est$DunedinPACE <- pace[rownames(age_est)]
        age_est <- data.frame(sample = rownames(age_est), age_est, row.names = NULL,
                               check.names = FALSE, stringsAsFactors = FALSE)

        pd <- tryCatch(pData(bs), error = function(e) NULL)
        if (!is.null(pd)) {
            pd_df <- as.data.frame(pd)
            if (ncol(pd_df) > 0) {
                if (verbose) message("[emulateArray]: appending pData() (", ncol(pd_df), " columns) to age estimates")
                rownames(pd_df) <- sample_ids
                age_est <- cbind(age_est, pd_df[age_est$sample, , drop = FALSE])
                rownames(age_est) <- NULL  # cbind() otherwise reinstates pd_df's sample-ID rownames
            }
        }

        if (!is.null(ageOut)) {
            write.csv(age_est, ageOut, row.names = FALSE)
            if (verbose) message("[emulateArray]: wrote ", ageOut)
        }
        result$ageEstimates <- age_est
    }

    return(result)
}
