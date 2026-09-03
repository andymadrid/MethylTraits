#' MethylTraits: Discover and Predict Traits from WGBS Methylation Data
#'
#' See the function index for available tools: [findICRs()], [annotateICRs()],
#' [findSEMs()], [annotateSEMs()], [testSEMBurden()], [testSEMGenes()],
#' [goSEMs()], [predictSex()], [predictRace()], [predictPASC()],
#' [predictAge()], [emulateArray()], [liftBuild()].
#'
#' Centralizes the core Bioconductor infrastructure imports (GenomicRanges and
#' its tightly-coupled dependencies IRanges/S4Vectors/GenomeInfoDb, which
#' between them supply constructors and accessors like `IRanges()`,
#' `queryHits()`, `runValue()`, and `seqlevelsInUse()` used bare throughout
#' this package) plus a handful of base-package functions (`methods::is`,
#' `stats::coef`/`kmeans`/`predict`, `utils::read.csv`/`write.csv`) used
#' unqualified in more than one file, rather than repeating `@importFrom` on
#' every individual function.
#'
#' @keywords internal
#' @import bsseq
#' @import GenomicRanges
#' @import IRanges
#' @import S4Vectors
#' @import GenomeInfoDb
#' @importFrom methods is
#' @importFrom stats coef kmeans predict
#' @importFrom utils read.csv write.csv
"_PACKAGE"

# emulateArray() references these two annotation-package names as bare data objects
# (the standard idiom once library(pkg) has attached them at runtime, done immediately
# before each reference) -- static analysis can't see that, so silence the resulting
# "no visible binding for global variable" NOTE rather than restructure working code.
utils::globalVariables(c(
    "IlluminaHumanMethylationEPICanno.ilm10b5.hg38",
    "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
))
