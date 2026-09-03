.liftChainID <- function(current, new) {
    if (current == "hg38" && new == "hg19") {
        "AH14108"
    } else if (current == "hg19" && new == "hg38") {
        "AH14150"
    } else {
        stop("current and new must differ, and both must be one of 'hg19'/'hg38' -- ",
             "currently only hg19<->hg38 lift-over is supported")
    }
}

#' Lift a bsseq object over from one human genome build to another
#'
#' Converts every CpG's coordinates in a \code{BSseq} object from one genome build to
#' another (hg19<->hg38) via a UCSC chain file, fetched through
#' \code{AnnotationHub}. CpGs that fail to lift over (no equivalent locus in the
#' target build) are dropped; the object's \code{genome()} is updated to the new
#' build. Use this to bring externally-aligned hg19 data onto hg38 (the build every
#' other function in this package -- \code{predictRace}, \code{predictPASC},
#' \code{predictAge}, \code{emulateArray}, and the default for \code{predictSex},
#' \code{annotateICRs}, \code{annotateSEMs} -- expects), or the reverse.
#'
#' @param bs A bsseq object, aligned to \code{current}.
#' @param current The genome build \code{bs} is currently aligned to. One of "hg19"
#'   or "hg38".
#' @param new The genome build to lift \code{bs} to. One of "hg19" or "hg38"; must
#'   differ from \code{current}.
#' @param verbose Logical. Print progress messages, including how many CpGs were
#'   lost in the lift-over. Default TRUE.
#' @return A bsseq object with coordinates on the \code{new} build, containing only
#'   the CpGs that successfully lifted over (in the same relative genomic order as
#'   the input).
#' @examples
#' \dontrun{
#' liftBuild(bs, current = "hg38", new = "hg19")
#' }
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom SummarizedExperiment rowRanges<-
#' @export
liftBuild <- function(bs, current = c("hg19", "hg38"), new = c("hg19", "hg38"), verbose = TRUE) {

    current <- match.arg(current)
    new <- match.arg(new)
    direction <- .liftChainID(current, new)

    if (!requireNamespace("AnnotationHub", quietly = TRUE) || !requireNamespace("rtracklayer", quietly = TRUE)) {
        stop("Packages AnnotationHub and rtracklayer are required but not installed")
    }

    if (verbose) message("[liftBuild]: will lift your bsseq object from ", current, " to ", new)
    # AH14108 = hg38ToHg19.over.chain.gz ; AH14150 = hg19ToHg38.over.chain.gz
    mcols(bs)$cpgs <- 1:length(bs)
    hgCurrent <- rowRanges(bs)
    hgCurrent$cpgs <- 1:length(hgCurrent)
    hub <- AnnotationHub::AnnotationHub()
    chain <- hub[[direction]]
    hgNew <- unlist(rtracklayer::liftOver(hgCurrent, chain))
    bsLifted <- bs[which(hgCurrent$cpgs %in% hgNew$cpgs)]
    rowRanges(bsLifted) <- hgNew
    genome(bsLifted) <- new
    lost <- length(bs) - length(bsLifted)
    if (verbose) message("[liftBuild]: lifting over resulted in a loss of ", lost, " CpGs")
    bsLifted
}
