.knownICRsJima2022 <- data.frame(
    name   = c("PEG3/ZIM2", "PEG10", "MEST/MESTIT1", "GRB10", "KCNQ1/KCNQ1OT1",
               "H19/MRPL23", "IGF2R/AIRN", "PLAGL1/HYMAI", "MEG3/DLK1",
               "SNRPN/SNHG14", "BLCAP/NNAT", "GNAS", "L3MBTL1", "INPP5F"),
    chr    = c("chr19", "chr7", "chr7", "chr7", "chr11",
               "chr11", "chr6", "chr6", "chr14",
               "chr15", "chr20", "chr20", "chr20", "chr10"),
    start  = c(56837320, 94656360, 130490640, 50781638, 2698157,
               1997886, 160005401, 144006941, 100824556,
               24954592, 37520202, 58839107, 43513725, 119817943),
    end    = c(56841439, 94658647, 130494200, 50783354, 2701210,
               2003118, 160006584, 144008825, 100828242,
               24956828, 37522993, 58856828, 43515256, 119819030),
    source = "Jima2022_Table1",
    stringsAsFactors = FALSE
)

#' Annotate candidate ICRs with nearest gene and known-ICR overlap
#'
#' Takes the data.frame (or GRanges) of candidate ICRs returned by
#' \code{findICRs} and adds: (1) the nearest gene / genomic feature
#' (promoter, exon, intron, intergenic, etc.) and distance-to-TSS, via
#' ChIPseeker + a UCSC knownGene TxDb, and (2) whether each candidate
#' overlaps a reference set of previously-published ICRs.
#'
#' The bundled reference set (used when \code{knownICRs = NULL}) is the 14
#' classical, literature-validated human ICRs reported as previously-known
#' loci recovered in Table 1 of Jima et al. 2022 (Epigenetics 17:1920-1943,
#' "Genomic map of candidate human imprint control regions: the imprintome"):
#' PEG3/ZIM2, PEG10, MEST/MESTIT1, GRB10, KCNQ1/KCNQ1OT1, H19/MRPL23,
#' IGF2R/AIRN, PLAGL1/HYMAI, MEG3/DLK1, SNRPN/SNHG14, BLCAP/NNAT, GNAS,
#' L3MBTL1, INPP5F (GRCh38 coordinates, transcribed directly from the paper).
#'
#' IMPORTANT CAVEAT: Jima et al. 2022 report 1,488 candidate ICRs genome-wide
#' in total, but the full coordinate table for all 1,488 (and the complete
#' 25-locus classical-ICR benchmark table) lives only in a supplementary
#' file that is bot-gated (proof-of-work JS challenge) behind PMC/Tandfonline
#' and in an interactive-only genome browser at https://humanicr.org/ --
#' neither could be retrieved programmatically. The bundled 14-locus table
#' below is therefore a verified but PARTIAL reference, not the full 1,488-
#' region catalog. If you can obtain the complete list (e.g. via
#' institutional journal access to the supplementary .docx, or by manually
#' exporting features from the humanicr.org browser), pass it via the
#' \code{knownICRs} argument to overlap your candidates against the full set.
#'
#' @param icrs A data.frame as returned by \code{findICRs} (needs seqnames,
#'   start, end columns), or a GRanges of candidate ICRs.
#' @param genomeBuild One of "hg38" (default) or "hg19". Selects which UCSC
#'   knownGene TxDb is used for gene annotation. Must match the genome build
#'   the bsseq object (and knownICRs, if supplied) were aligned to.
#' @param knownICRs Optional data.frame with columns chr, start, end, name
#'   (and optionally source) giving a reference set of known/published ICRs
#'   to overlap candidates against, in the same genome build as \code{icrs}.
#'   Default NULL uses the bundled 14-locus Jima et al. 2022 Table 1 set.
#' @param verbose Logical, print progress messages. Default TRUE.
#' @return The input data.frame with additional columns: \code{featureType}
#'   (collapsed genomic feature category: Promoter, Exon, Intron, 5'/3' UTR,
#'   Downstream, or Distal Intergenic), \code{annotation} (the full-detail
#'   ChIPseeker string, e.g. "Promoter (<=1kb)" or the specific exon/intron),
#'   \code{nearestGene} (gene symbol), \code{distanceToTSS},
#'   \code{overlapsKnownICR} (logical), \code{knownICR_name},
#'   \code{knownICR_source}.
#' @examples
#' \dontrun{
#' icrs <- findICRs(bs, ncores = 20)
#' icrs.annotated <- annotateICRs(icrs)
#' }
#' @export

annotateICRs <- function(icrs, genomeBuild = c("hg38", "hg19"), knownICRs = NULL, verbose = TRUE) {

    genomeBuild <- match.arg(genomeBuild)

    if (is.data.frame(icrs)) {
        if (!all(c("seqnames", "start", "end") %in% colnames(icrs))) {
            stop("icrs data.frame must have seqnames, start, end columns (as returned by findICRs)")
        }
        icrs.gr <- GenomicRanges::GRanges(icrs$seqnames, IRanges::IRanges(icrs$start, icrs$end))
    } else if (is(icrs, "GRanges")) {
        icrs.gr <- icrs
    } else {
        stop("icrs must be a data.frame (from findICRs) or a GRanges")
    }
    if (length(icrs.gr) == 0) {
        stop("icrs is empty -- nothing to annotate")
    }

    if (verbose) message(paste0("[annotateICRs]: loading ", genomeBuild, " gene annotation"))
    if (genomeBuild == "hg38") {
        if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
            stop("Package TxDb.Hsapiens.UCSC.hg38.knownGene is required but not installed")
        }
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    } else {
        if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE)) {
            stop("Package TxDb.Hsapiens.UCSC.hg19.knownGene is required but not installed")
        }
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    }
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        stop("Package org.Hs.eg.db is required but not installed")
    }

    if (verbose) message("[annotateICRs]: annotating nearest gene / genomic feature for each candidate")
    peakAnno <- ChIPseeker::annotatePeak(icrs.gr, TxDb = txdb, annoDb = "org.Hs.eg.db",
                                          verbose = FALSE, sameStrand = FALSE)
    anno.df <- as.data.frame(peakAnno)

    result <- data.frame(
        seqnames = anno.df$seqnames,
        start = anno.df$start,
        end = anno.df$end,
        width = anno.df$width,
        featureType = sub(" \\(.*", "", anno.df$annotation),
        annotation = anno.df$annotation,
        nearestGene = if ("SYMBOL" %in% colnames(anno.df)) anno.df$SYMBOL else NA_character_,
        distanceToTSS = anno.df$distanceToTSS,
        stringsAsFactors = FALSE
    )

    if (verbose) message("[annotateICRs]: checking overlap with known ICRs")
    if (is.null(knownICRs)) {
        knownICRs <- .knownICRsJima2022
        if (verbose) message(paste0("[annotateICRs]: using bundled 14-locus Jima et al. 2022 Table 1 reference ",
                                     "(not the full 1,488-region catalog -- see ?annotateICRs)"))
    }
    if (!all(c("chr", "start", "end", "name") %in% colnames(knownICRs))) {
        stop("knownICRs must have columns: chr, start, end, name (and optionally source)")
    }
    known.gr <- GenomicRanges::GRanges(knownICRs$chr, IRanges::IRanges(knownICRs$start, knownICRs$end))

    ov <- GenomicRanges::findOverlaps(icrs.gr, known.gr)
    result$overlapsKnownICR <- FALSE
    result$knownICR_name <- NA_character_
    result$knownICR_source <- NA_character_
    if (length(ov) > 0) {
        hitTab <- data.frame(query = S4Vectors::queryHits(ov), subject = S4Vectors::subjectHits(ov))
        agg <- do.call(rbind, lapply(split(hitTab$subject, hitTab$query), function(s) {
            data.frame(name = paste(unique(knownICRs$name[s]), collapse = ";"),
                       source = paste(unique(knownICRs$source[s]), collapse = ";"))
        }))
        qIdx <- as.integer(names(split(hitTab$subject, hitTab$query)))
        result$overlapsKnownICR[qIdx] <- TRUE
        result$knownICR_name[qIdx] <- agg$name
        result$knownICR_source[qIdx] <- agg$source
    }

    nOverlap <- sum(result$overlapsKnownICR)
    if (verbose) message(paste0("[annotateICRs]: ", nOverlap, " of ", nrow(result),
                                 " candidate ICRs overlap a known/reference ICR"))

    return(result)
}
