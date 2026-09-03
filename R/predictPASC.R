#' Predict PASC (long-COVID) status from DNA methylation
#'
#' Applies a 5-model ensemble (elastic net, linear SVM, radial SVM, GBM, random
#' forest) plus a logistic-regression stacking meta-learner, trained on
#' WGBS-emulated EPIC array data, to estimate PASC probability per sample.
#'
#' \strong{IMPORTANT CALIBRATION CAVEAT}: this model's probability threshold (0.5)
#' is validated ONLY on the platform it was trained on (WGBS-emulated EPIC data).
#' Applied to real microarray or WGBS/bsseq data from a different platform,
#' false-positive rates were high even after ComBat harmonization, and a decisive
#' diagnostic (mixed AUROC of true positives vs. harmonized real-array negatives =
#' 0.11-0.22) showed the absolute score does not reliably transfer across
#' platforms. Treat \code{predictedStatus} on non-emulated-platform input as
#' EXPLORATORY, not a validated clinical call. Relative ranking of samples WITHIN
#' one batch/dataset (via \code{prob_enet}) is more trustworthy than any single
#' sample's absolute probability. See the pipeline's own analysis report
#' (\code{pipeline/report/pasc_report.html} alongside the model bundle) for the
#' full evidence behind this caveat.
#'
#' Uses \code{.mtModelDir()} and \code{.autoDetectCGOffset()} from
#' predictUtils.R -- source that file alongside this one (or simply
#' \code{library(MethylTraits)} once packaged).
#'
#' @param bs A bsseq object aligned to hg38, OR a numeric matrix of beta values
#'   with CpGs in rows (Illumina probe IDs, e.g. "cg04206906") and samples in
#'   columns.
#' @param offset Integer coordinate offset applied to \code{bs} positions before
#'   matching the model's hg38 CpG panel (ignored for matrix input). Default NA
#'   auto-detects it via \code{.autoDetectCGOffset()} (see \code{predictRace}).
#' @param modelDir Directory containing \code{pasc_model_bundle.rds} and
#'   \code{panel_cpgs_hg38.csv}. Default NULL resolves to the bundled model.
#' @param outputFile Optional CSV path to also write the result to. Default NULL.
#' @param verbose Logical. Print progress messages. Default TRUE.
#' @return A data.frame, one row per sample: \code{sample}, \code{prob_enet},
#'   \code{prob_svmlin}, \code{prob_gbm}, \code{prob_rf}, \code{prob_svmrad}
#'   (per-model PASC probability), \code{prob_stack_ensemble} (meta-learner
#'   combination), \code{panelCoverageSelected}, \code{panelCoverageTopK}
#'   (fraction of each model's CpG panel actually observed, vs. filled from
#'   training-set means), and \code{predictedStatus} ("PASC"/"non-PASC" at the
#'   0.5 threshold on \code{prob_enet} -- see the calibration caveat above).
#' @examples
#' \dontrun{
#' predictPASC(bs)
#' }
#' @import bsseq
#' @import GenomicRanges
#' @export

predictPASC <- function(bs, offset = NA_integer_, modelDir = NULL, outputFile = NULL, verbose = TRUE) {

    # these must be attached (not just namespace-loaded) for predict.train() to
    # correctly dispatch to each base learner's underlying prediction method --
    # requireNamespace() alone was found to silently break glmnet's coef.multnet()
    # (see predictRace()), so the same safer pattern is used for every base learner here.
    suppressPackageStartupMessages({
        library(caret); library(glmnet); library(randomForest); library(e1071); library(gbm); library(kernlab)
    })

    if (is.null(modelDir)) modelDir <- .mtModelDir("pasc")
    bundle_file <- file.path(modelDir, "pasc_model_bundle.rds")
    if (!file.exists(bundle_file)) stop("Model bundle not found: ", bundle_file)

    if (verbose) message("[predictPASC]: loading PASC model bundle")
    bundle <- readRDS(bundle_file)
    if (verbose) message(sprintf("[predictPASC]: selected CpG panel: %d (top-K subset used by RF/SVM/GBM: %d)",
                                  length(bundle$selected_cpgs), length(bundle$top_cpgs)))

    build_panel_matrix <- function(beta_mat, panel, means) {
        samples <- colnames(beta_mat)
        avail <- intersect(panel, rownames(beta_mat))
        X <- matrix(means[panel], nrow = length(samples), ncol = length(panel),
                    byrow = TRUE, dimnames = list(samples, panel))
        if (length(avail) > 0) X[, avail] <- t(beta_mat[avail, , drop = FALSE])
        if (anyNA(X)) {
            fallback <- matrix(means[panel], nrow = nrow(X), ncol = length(panel), byrow = TRUE)
            X[is.na(X)] <- fallback[is.na(X)]
        }
        list(X = X, coverage = length(avail) / length(panel))
    }

    predict_all_models <- function(bundle, beta_mat) {
        samples <- colnames(beta_mat)
        if (length(samples) == 0) stop("no sample columns to predict on (colnames(beta_mat) is empty)")
        Xf <- build_panel_matrix(beta_mat, bundle$selected_cpgs, bundle$training_cpg_means)
        Xt <- build_panel_matrix(beta_mat, bundle$top_cpgs, bundle$training_cpg_means)
        if (verbose) {
            message(sprintf("[predictPASC]: selected-panel coverage: %d/%d (%.1f%% filled from training means)",
                             length(intersect(bundle$selected_cpgs, rownames(beta_mat))), length(bundle$selected_cpgs),
                             100 * (1 - Xf$coverage)))
            message(sprintf("[predictPASC]: top-K panel coverage: %d/%d (%.1f%% filled from training means)",
                             length(intersect(bundle$top_cpgs, rownames(beta_mat))), length(bundle$top_cpgs),
                             100 * (1 - Xt$coverage)))
            if (Xf$coverage < 0.9) {
                message("[predictPASC]: WARNING >10% of the selected CpG panel is missing from this input and was ",
                        "filled with training-set means. Expected for a different array generation/platform, but ",
                        "reduces how much this prediction reflects your actual data.")
            }
        }

        out <- list()
        for (mname in names(bundle$final_models)) {
            X <- if (mname == "enet") Xf$X else Xt$X
            p <- tryCatch(predict(bundle$final_models[[mname]], newdata = X, type = "prob")[, "PASC"],
                          error = function(e) { warning("prediction failed for ", mname, ": ", conditionMessage(e)); rep(NA_real_, length(samples)) })
            out[[mname]] <- p
        }
        wide <- as.data.frame(out)
        colnames(wide) <- names(bundle$final_models)
        stack_p <- tryCatch(predict(bundle$stack_meta, newdata = wide, type = "response"),
                             error = function(e) rep(NA_real_, length(samples)))

        result <- data.frame(sample = samples, prob_enet = wide$enet, prob_svmlin = wide$svmlin,
                              prob_gbm = wide$gbm, prob_rf = wide$rf, prob_svmrad = wide$svmrad,
                              prob_stack_ensemble = stack_p, panelCoverageSelected = Xf$coverage,
                              panelCoverageTopK = Xt$coverage, stringsAsFactors = FALSE)
        result$predictedStatus <- ifelse(result$prob_enet > 0.5, "PASC", "non-PASC")
        result
    }

    if (is(bs, "BSseq")) {
        panel_coords_file <- file.path(modelDir, "panel_cpgs_hg38.csv")
        if (!file.exists(panel_coords_file)) stop("Panel coordinate file not found: ", panel_coords_file)

        bs_sample_names <- sampleNames(bs)
        if (is.null(bs_sample_names) || length(bs_sample_names) != ncol(bs) ||
            any(is.na(bs_sample_names)) || any(duplicated(bs_sample_names))) {
            if (verbose) message("[predictPASC]: bsseq object has no usable sampleNames -- using generic Sample1..N labels")
            bs_sample_names <- paste0("Sample", seq_len(ncol(bs)))
        }

        cov <- getCoverage(bs)
        meth <- getCoverage(bs, type = "M")
        beta_vals <- as.matrix(meth / pmax(cov, 1))
        beta_vals[cov == 0] <- NA

        gr <- granges(bs)
        chrs <- as.character(seqnames(gr))
        if (!any(grepl("^chr", chrs))) chrs <- paste0("chr", chrs)
        raw_pos <- start(gr)

        if (is.na(offset)) {
            if (verbose) message("[predictPASC]: auto-detecting coordinate convention (0-based vs 1-based) against hg38")
            offset <- .autoDetectCGOffset(chrs, raw_pos, verbose = verbose)
        } else if (verbose) {
            message("[predictPASC]: using user-specified offset ", offset)
        }

        pos <- raw_pos + offset
        site_id <- paste0(chrs, ":", pos)

        panel_anno <- read.csv(panel_coords_file, stringsAsFactors = FALSE)
        panel_site_id <- paste0(panel_anno$chr, ":", panel_anno$pos)
        m <- match(panel_site_id, site_id)
        n_matched <- sum(!is.na(m))
        if (verbose) message(sprintf("[predictPASC]: matched %d/%d panel CpGs by hg38 coordinate (%.1f%%)",
                                      n_matched, length(panel_site_id), 100 * n_matched / length(panel_site_id)))
        if (n_matched / length(panel_site_id) < 0.5) {
            warning("[predictPASC]: less than half the panel matched by coordinate -- check genome build (hg38 expected) and offset")
        }

        beta_mat <- matrix(NA_real_, nrow = length(panel_site_id), ncol = ncol(bs),
                            dimnames = list(panel_anno$cpg, bs_sample_names))
        ok <- !is.na(m)
        beta_mat[ok, ] <- beta_vals[m[ok], , drop = FALSE]

    } else if (is.matrix(bs) || is.data.frame(bs)) {
        beta_mat <- as.matrix(bs)
        mode(beta_mat) <- "numeric"
        if (verbose) message(sprintf("[predictPASC]: input matrix %d CpGs x %d samples", nrow(beta_mat), ncol(beta_mat)))
    } else {
        stop("bs must be a BSseq object or a numeric matrix (CpGs x samples)")
    }

    if (verbose) message("[predictPASC]: predicting PASC status")
    result <- predict_all_models(bundle, beta_mat)
    rownames(result) <- NULL

    if (verbose) {
        message("[predictPASC]: prediction distribution (elastic net, primary model, threshold=0.5):")
        print(table(result$predictedStatus))
        message("[predictPASC]: *** CALIBRATION CAVEAT: the 0.5 threshold is validated only within the ",
                "training platform (WGBS-emulated EPIC data). On real microarray/bsseq input, treat ",
                "predictedStatus as exploratory -- rely on relative ranking of prob_enet across your ",
                "samples rather than the absolute cutoff. ***")
    }
    if (!is.null(outputFile)) {
        write.csv(result, outputFile, row.names = FALSE)
        if (verbose) message("[predictPASC]: wrote ", outputFile)
    }

    return(result)
}
