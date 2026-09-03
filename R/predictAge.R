#' Estimate epigenetic age from an arbitrary, unfiltered bsseq object (blood WGBS)
#'
#' Reports two predictions per sample, both from the accompanying methodology
#' notes (\code{new_age_clock/METHODS.md} alongside the original model artifacts):
#' \itemize{
#'   \item \code{wgbs_only_stack} -- Tier 1 (PRIMARY): a pure-WGBS ensemble stack
#'     (glmnet + random forest + SVM + xgboost base learners), zero published-clock
#'     inputs; the independently novel contribution.
#'   \item \code{hybrid_stack} -- Tier 2 (secondary/applied): adds Zhang2019 and
#'     Horvath2 (via methylCIPHER) as 2 extra stacking inputs; best raw accuracy,
#'     but not an independently-trained comparison against those published clocks.
#' }
#' No coverage/quality pre-filtering is required on the input -- this applies the
#' project's own >= \code{covMin} threshold internally and imputes anything below
#' it (or simply missing/uncovered) from the training-set per-CpG median.
#' Coordinate handling: a systematic 0-based vs 1-based offset in the input is
#' auto-detected (by comparing overlap counts against the model's own CpG panel at
#' offsets -1/0/+1) and corrected before matching -- this matters because some
#' upstream pipelines (e.g. BED-derived CpG reports) are 0-based, while every
#' model here was built against 1-based Bioconductor GRanges coordinates; a silent
#' mismatch would otherwise cause near-total CpG-matching failure.
#'
#' Uses \code{.mtModelDir()} from predictUtils.R -- source that file
#' alongside this one (or simply \code{library(MethylTraits)} once packaged). Also
#' requires the methylCIPHER package (for calcZhang2019/calcHorvath2) to be
#' installed.
#'
#' @param bs A bsseq object (blood WGBS, hg38, any coverage -- unfiltered is fine).
#' @param modelDir Directory containing the model artifacts (\code{final_models.rds},
#'   \code{xgb_final.model}, \code{cpg_train_median.rds}, \code{probe_gr_hg38.rds},
#'   \code{stack_best.rds}, \code{stack_published_best.rds}). Default NULL resolves
#'   to the bundled model.
#' @param covMin Minimum read coverage to call a beta value at a CpG; below this
#'   (or if absent) the value is imputed from the training-set per-CpG median.
#'   Default 5.
#' @param outputFile Optional CSV path to also write the full result to. Default NULL.
#' @param verbose Logical. Print progress messages. Default TRUE.
#' @return A data.frame, one row per sample: \code{sample}, the 4 base-learner ages
#'   (\code{glmnet}, \code{rf}, \code{svm}, \code{xgb}), \code{wgbs_only_stack}
#'   (Tier 1, PRIMARY), \code{zhang2019}, \code{horvath2} (published-clock inputs
#'   used only by Tier 2), \code{hybrid_stack} (Tier 2), and diagnostics
#'   \code{coordOffsetApplied}, \code{ownPanelPctObserved}, \code{clockProbesPctObserved}.
#' @examples
#' \dontrun{
#' predictAge(bs)
#' }
#' @import bsseq
#' @import GenomicRanges
#' @export

predictAge <- function(bs, modelDir = NULL, covMin = 5, outputFile = NULL, verbose = TRUE) {

    if (!is(bs, "BSseq")) stop("bs must be a BSseq object")
    if (!requireNamespace("methylCIPHER", quietly = TRUE)) stop("Package methylCIPHER is required but not installed")
    # attached (not just namespace-loaded), matching predictRace()/predictPASC() --
    # see predictRace()'s note on coef.multnet() needing glmnet attached.
    suppressPackageStartupMessages({
        library(glmnet); library(ranger); library(e1071); library(xgboost); library(methylCIPHER)
    })

    if (is.null(modelDir)) modelDir <- .mtModelDir("age")

    if (verbose) message("[predictAge]: loading model artifacts")
    final_models <- readRDS(file.path(modelDir, "final_models.rds"))
    xgb_final <- xgboost::xgb.load(file.path(modelDir, "xgb_final.model"))
    cpg_gr <- final_models$glmnet$cpg_gr                              # own Tier-1 panel
    cpg_median <- readRDS(file.path(modelDir, "cpg_train_median.rds"))
    probe_gr <- readRDS(file.path(modelDir, "probe_gr_hg38.rds"))     # published-clock probe panel
    stack4 <- readRDS(file.path(modelDir, "stack_best.rds"))          # Tier 1: glmnet/rf/svm/xgb
    stack_pub <- readRDS(file.path(modelDir, "stack_published_best.rds"))  # Tier 2
    stopifnot(identical(stack4$base_learners, c("glmnet", "rf", "svm", "xgb")))
    stopifnot(identical(stack_pub$base_learners, c("glmnet", "rf", "svm", "xgb", "zhang2019", "horvath2")))
    if (verbose) message(sprintf("[predictAge]: OK -- %d own-panel CpGs, %d clock probes", length(cpg_gr), length(probe_gr)))

    sn <- sampleNames(bs)
    if (is.null(sn) || length(sn) == 0 || all(sn == "")) sn <- rownames(pData(bs))
    if (is.null(sn) || length(sn) == 0 || all(sn == "")) sn <- as.character(seq_len(ncol(bs)))
    n_val <- ncol(bs)

    if (verbose) message("[predictAge]: detecting coordinate convention (0-based vs 1-based)")
    gr_val <- granges(bs); strand(gr_val) <- "*"
    ref_test <- if (length(cpg_gr) > 300000) cpg_gr[sample(length(cpg_gr), 300000)] else cpg_gr
    shift_gr <- function(gr, off) { start(gr) <- start(gr) + off; end(gr) <- end(gr) + off; gr }
    offsets <- c(0L, -1L, 1L)
    match_counts <- sapply(offsets, function(off) {
        sum(!is.na(findOverlaps(shift_gr(ref_test, off), gr_val, type = "equal", ignore.strand = TRUE, select = "first")))
    })
    names(match_counts) <- offsets
    if (verbose) for (i in seq_along(offsets)) message(sprintf("  offset %+d : %d / %d matched (%.1f%%)",
        offsets[i], match_counts[i], length(ref_test), 100 * match_counts[i] / length(ref_test)))
    best_offset <- offsets[which.max(match_counts)]
    best_rate <- max(match_counts) / length(ref_test)
    if (verbose) {
        if (best_offset != 0L) {
            message(sprintf("[predictAge]: applying a %+d bp coordinate correction (best match rate %.1f%% vs %.1f%% at offset 0)",
                             best_offset, 100 * best_rate, 100 * match_counts["0"] / length(ref_test)))
        } else {
            message("[predictAge]: no coordinate correction needed (offset 0 already matches best)")
        }
    }
    if (best_rate < 0.5) {
        warning("[predictAge]: best match rate is only ", round(100 * best_rate, 1), "% even after offset correction -- ",
                "this may indicate a genome build mismatch (e.g. hg19 vs hg38) or a different chromosome-naming ",
                "convention, not just a 0/1-based issue. Predictions below should be treated with caution.")
    }

    extract_beta_at <- function(bs, target_gr, offset) {
        query_gr <- shift_gr(target_gr, offset)
        match_idx <- findOverlaps(query_gr, gr_val, type = "equal", ignore.strand = TRUE, select = "first")
        matched_which <- which(!is.na(match_idx))
        beta_val <- matrix(NA_real_, nrow = n_val, ncol = length(target_gr))
        if (length(matched_which) > 0) {
            val_rows <- match_idx[matched_which]
            M_sub <- as.matrix(getCoverage(bs, type = "M")[val_rows, , drop = FALSE])
            Cov_sub <- as.matrix(getCoverage(bs, type = "Cov")[val_rows, , drop = FALSE])
            beta_sub <- M_sub / Cov_sub
            beta_sub[Cov_sub < covMin] <- NA
            beta_val[, matched_which] <- t(beta_sub)
        }
        list(beta = beta_val, pct_matched = 100 * length(matched_which) / length(target_gr),
             pct_observed = 100 * mean(!is.na(beta_val)))
    }

    if (verbose) message(sprintf("[predictAge]: extracting beta values (>=%dx coverage)", covMin))
    own <- extract_beta_at(bs, cpg_gr, best_offset)
    clk <- extract_beta_at(bs, probe_gr, best_offset)
    if (verbose) {
        message(sprintf("  own panel   (%d CpGs):   matched %.1f%%, observed %.1f%%", length(cpg_gr), own$pct_matched, own$pct_observed))
        message(sprintf("  clock probes (%d probes): matched %.1f%%, observed %.1f%%", length(probe_gr), clk$pct_matched, clk$pct_observed))
    }

    beta_own <- own$beta
    na_mask <- is.na(beta_own)
    med_mat <- matrix(rep(cpg_median, each = n_val), nrow = n_val)
    beta_own[na_mask] <- med_mat[na_mask]
    colnames(beta_own) <- paste0("cpg_", seq_len(ncol(beta_own)))
    beta_clk <- clk$beta
    colnames(beta_clk) <- names(probe_gr)

    if (verbose) message("[predictAge]: predicting age")
    preds <- data.frame(sample = sn)
    preds$glmnet <- as.numeric(predict(final_models$glmnet$fit, newx = beta_own, s = "lambda.min"))
    preds$rf <- predict(final_models$rf$fit, data = beta_own[, final_models$rf$cpg_idx, drop = FALSE])$predictions
    preds$svm <- as.numeric(predict(final_models$svm$fit, beta_own[, final_models$svm$cpg_idx, drop = FALSE]))
    preds$xgb <- predict(xgb_final, xgboost::xgb.DMatrix(beta_own[, final_models$xgb$cpg_idx, drop = FALSE]))
    preds$wgbs_only_stack <- as.numeric(predict(stack4$fit,
        newx = as.matrix(preds[, c("glmnet", "rf", "svm", "xgb")]), s = "lambda.min"))

    preds$zhang2019 <- as.numeric(methylCIPHER::calcZhang2019(beta_clk))
    preds$horvath2 <- as.numeric(methylCIPHER::calcHorvath2(beta_clk))
    preds$hybrid_stack <- as.numeric(predict(stack_pub$fit,
        newx = as.matrix(preds[, c("glmnet", "rf", "svm", "xgb", "zhang2019", "horvath2")]), s = "lambda.min"))

    preds$coordOffsetApplied <- best_offset
    preds$ownPanelPctObserved <- round(own$pct_observed, 2)
    preds$clockProbesPctObserved <- round(clk$pct_observed, 2)

    if (verbose) {
        message("[predictAge]: summary (Tier 1 PRIMARY = wgbs_only_stack; Tier 2 secondary/applied = hybrid_stack)")
        print(preds[, c("sample", "wgbs_only_stack", "hybrid_stack")], row.names = FALSE, digits = 4)
    }
    if (!is.null(outputFile)) {
        write.csv(preds, outputFile, row.names = FALSE)
        if (verbose) message("[predictAge]: wrote ", outputFile)
    }

    return(preds)
}
