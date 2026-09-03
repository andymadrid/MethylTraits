#' Predict self-reported race/ethnicity category from DNA methylation
#'
#' Applies a pre-trained elastic-net multinomial model (glmnet, trained on a large
#' cross-GEO compendium of blood methylation data) to predict the self-reported
#' race/ethnicity category recorded in each training sample's metadata --
#' African_American, Asian, Hispanic, or White. This is a DNAm-based proxy for
#' genetic ancestry and correlated social/environmental exposure, not a biological
#' definition of race; treat predictions accordingly, especially for individuals of
#' mixed or otherwise unrepresented ancestry.
#'
#' Only the model's 1,954 nonzero-coefficient CpGs are needed as input (out of the
#' 100,000 candidate CpGs it was originally selected from) -- every coefficient
#' outside that set is exactly zero at the model's chosen lambda, so it can only
#' ever contribute 0 to the class score regardless of its input value. This mirrors
#' what the original CLI version of this model already did for bsseq/WGBS input
#' (it only differed, needlessly, for plain-matrix input -- see \code{betaMatrix}).
#'
#' Uses \code{.mtModelDir()} and \code{.autoDetectCGOffset()} from
#' predictUtils.R -- source that file alongside this one (or simply
#' \code{library(MethylTraits)} once packaged).
#'
#' @param bs A bsseq object aligned to hg38, OR a numeric matrix of beta values
#'   (methylation fractions in 0 to 1) with CpGs in rows (named with Illumina probe
#'   IDs, e.g. "cg17992056") and samples in columns -- e.g. from an EPIC/450K array,
#'   where genomic coordinates aren't needed since probe IDs already identify sites.
#' @param offset Integer coordinate offset applied to \code{bs} positions before
#'   matching against the model's hg38 CpG coordinates (ignored if \code{bs} is a
#'   plain matrix). Default NA auto-detects it: a sample of candidate positions is
#'   checked against the hg38 reference sequence for whether it reads "CG" at the
#'   reported (0) or +1-shifted (1) start, the same auto-detection approach used by
#'   \code{predictPASC} and \code{predictAge} in this package. Requires
#'   BSgenome.Hsapiens.UCSC.hg38 and Biostrings; if unavailable, falls back to 0
#'   (already 1-based, standard Bioconductor convention) with a warning.
#' @param modelDir Directory containing \code{race_model.rds} and
#'   \code{predictive_cpgs_hg38.csv}. Default NULL resolves to the bundled model
#'   (via \code{system.file} once installed as a package, else the sibling
#'   \code{models/race/} directory next to this script).
#' @param outputFile Optional CSV path to also write the result to. Default NULL
#'   (return only).
#' @param verbose Logical. Print progress messages. Default TRUE.
#' @return A data.frame, one row per sample: \code{sample}, \code{predictedRace},
#'   \code{prob_<class>} for each of the 4 classes, \code{maxProb} (confidence of
#'   the winning class), and (for \code{bs} input) \code{cpgCoveragePct} (percent of
#'   the model's 1,954 CpGs actually observed in this sample, vs. imputed).
#' @examples
#' \dontrun{
#' predictRace(bs)
#' predictRace(betaMatrix)  # e.g. array beta matrix, CpGs x samples
#' }
#' @import bsseq
#' @import GenomicRanges
#' @export

predictRace <- function(bs, offset = NA_integer_, modelDir = NULL, outputFile = NULL, verbose = TRUE) {

    # glmnet must be attached (not just namespace-loaded) for coef.multnet() to
    # dispatch correctly on the multinomial model below -- requireNamespace() alone
    # silently returns NULL coefficients.
    suppressPackageStartupMessages({ library(glmnet); library(Matrix) })

    if (is.null(modelDir)) modelDir <- .mtModelDir("race")
    model_file <- file.path(modelDir, "race_model.rds")
    cpg_file <- file.path(modelDir, "predictive_cpgs_hg38.csv")
    if (!file.exists(model_file)) stop("Model file not found: ", model_file)
    if (!file.exists(cpg_file)) stop("Predictive CpG coordinate table not found: ", cpg_file)

    if (verbose) message("[predictRace]: loading race prediction model")
    model_obj <- readRDS(model_file)
    nonzero_cpgs <- model_obj$nonzero_cpgs
    training_cpgs <- model_obj$top_cpgs_final
    classes <- model_obj$classes
    if (verbose) message(paste0("[predictRace]: predictive CpGs: ", length(nonzero_cpgs),
                                 "; predicting: ", paste(classes, collapse = ", ")))

    if (is(bs, "BSseq")) {
        cpg_coords <- read.csv(cpg_file, stringsAsFactors = FALSE)
        pred_cpg_ids <- paste0(cpg_coords$chr_hg38, ":", cpg_coords$pos_hg38)
        cpg_names <- cpg_coords$cpg

        if (is.na(offset)) {
            if (verbose) message("[predictRace]: auto-detecting coordinate convention (0-based vs 1-based)")
            offset <- .autoDetectCGOffset(as.character(seqnames(bs)), start(bs), verbose = verbose)
        }

        cov <- getCoverage(bs)
        meth <- getCoverage(bs, type = "M")
        beta_vals <- as.matrix(meth / pmax(cov, 1))
        beta_vals[cov == 0] <- NA

        chrs <- as.character(seqnames(bs))
        if (!any(grepl("^chr", chrs))) chrs <- paste0("chr", chrs)
        pos <- start(bs) + offset
        site_ids <- paste0(chrs, ":", pos)

        coord_match <- match(pred_cpg_ids, site_ids)
        n_matched <- sum(!is.na(coord_match))
        pct_missing <- 100 * (length(nonzero_cpgs) - n_matched) / length(nonzero_cpgs)
        if (verbose) message(sprintf("[predictRace]: matched %d/%d predictive CpGs by coordinate (%.1f%% missing)",
                                      n_matched, length(nonzero_cpgs), pct_missing))
        if (pct_missing > 50) warning("[predictRace]: >50% of predictive CpGs could not be matched by coordinate -- ",
                                       "check genome build (hg38 expected) and offset")

        sample_names <- sampleNames(bs)
        if (is.null(sample_names) || length(sample_names) != ncol(bs)) sample_names <- paste0("Sample", seq_len(ncol(bs)))
        X_pred <- matrix(0.5, nrow = length(nonzero_cpgs), ncol = ncol(bs), dimnames = list(cpg_names, sample_names))
        matched_which <- which(!is.na(coord_match))
        if (length(matched_which) > 0) X_pred[matched_which, ] <- beta_vals[coord_match[matched_which], , drop = FALSE]
        cpgCoveragePct <- rep(100 * n_matched / length(nonzero_cpgs), ncol(bs))

    } else if (is.matrix(bs) || is.data.frame(bs)) {
        betas <- as.matrix(bs)
        mode(betas) <- "numeric"
        sample_names <- colnames(betas)
        if (is.null(sample_names)) sample_names <- paste0("Sample", seq_len(ncol(betas)))
        if (verbose) message(sprintf("[predictRace]: input matrix %d CpGs x %d samples", nrow(betas), ncol(betas)))

        available_train <- intersect(training_cpgs, rownames(betas))
        if (verbose) message(sprintf("[predictRace]: training CpG coverage: %d/%d", length(available_train), length(training_cpgs)))
        X_pred_full <- matrix(0.5, nrow = length(training_cpgs), ncol = ncol(betas),
                               dimnames = list(training_cpgs, sample_names))
        X_pred_full[available_train, ] <- betas[available_train, , drop = FALSE]
        X_pred <- X_pred_full[intersect(nonzero_cpgs, training_cpgs), , drop = FALSE]
        cpgCoveragePct <- rep(100 * length(intersect(nonzero_cpgs, rownames(betas))) / length(nonzero_cpgs), ncol(betas))

    } else {
        stop("bs must be a BSseq object or a numeric matrix (CpGs x samples)")
    }

    # impute any remaining NAs (row mean; fall back to neutral 0.5 if a row is all-NA)
    row_means <- rowMeans(X_pred, na.rm = TRUE)
    row_means[is.na(row_means)] <- 0.5
    for (k in which(rowSums(is.na(X_pred)) > 0)) X_pred[k, is.na(X_pred[k, ])] <- row_means[k]

    if (verbose) message("[predictRace]: predicting")
    coef_list <- coef(model_obj$final_model, s = model_obj$best_lambda)
    X_t <- t(X_pred)
    n_samples <- nrow(X_t)
    K <- length(classes)
    scores <- matrix(NA_real_, n_samples, K, dimnames = list(sample_names, classes))
    for (k in seq_along(classes)) {
        cls <- classes[k]
        beta <- as.matrix(coef_list[[cls]])
        intercept_k <- beta[1, 1]
        feat_betas <- beta[-1, 1, drop = TRUE]
        nz_feats <- names(feat_betas)[feat_betas != 0]
        in_both <- intersect(nz_feats, colnames(X_t))
        scores[, k] <- intercept_k + if (length(in_both) > 0) X_t[, in_both, drop = FALSE] %*% feat_betas[in_both] else 0
    }

    predicted_prob <- t(apply(scores, 1, function(s) { e <- exp(s - max(s)); e / sum(e) }))
    colnames(predicted_prob) <- classes
    predicted_class <- classes[apply(scores, 1, which.max)]

    result <- data.frame(sample = sample_names, predictedRace = predicted_class, stringsAsFactors = FALSE)
    prob_df <- as.data.frame(predicted_prob)
    colnames(prob_df) <- paste0("prob_", colnames(prob_df))
    result <- cbind(result, prob_df)
    result$maxProb <- apply(prob_df, 1, max)
    result$cpgCoveragePct <- round(cpgCoveragePct, 1)
    rownames(result) <- NULL

    if (verbose) {
        message("[predictRace]: prediction distribution:")
        print(table(result$predictedRace))
    }
    if (!is.null(outputFile)) {
        write.csv(result, outputFile, row.names = FALSE)
        if (verbose) message("[predictRace]: wrote ", outputFile)
    }

    return(result)
}
