## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
data_path <- "/media/data/dormant/Epigenetic_Age/bs.pasc3.rda"
data_available <- file.exists(data_path)
knitr::opts_chunk$set(eval = data_available)

## ----echo = FALSE, results = "asis"-------------------------------------------
if (!data_available) {
  cat(
    "**Note:** this vignette was built without access to the internal validation ",
    "dataset (`bs.pasc3.rda`) used to demonstrate every function below, so the ",
    "code chunks were not evaluated. Every example runs the same way on any ",
    "`BSseq` object -- just substitute your own."
  )
}

## ----load-data----------------------------------------------------------------
suppressPackageStartupMessages(library(bsseq))
library(MethylTraits)

load(data_path)
bs <- bs.pasc3
rm(bs.pasc3)

# a 3-sample subset keeps the *prediction* examples fast; predictSex/predictRace/
# predictPASC/predictAge/emulateArray all work identically on cohorts of any size
bs_small <- bs[, 1:3]
sampleNames(bs_small) <- pData(bs_small)$"Sample ID"
bs_small

## ----lift-build---------------------------------------------------------------
bs_small_hg38 <- bs_small[as.character(seqnames(bs_small)) %in% c("chr21", "chr22"), ]
bs_hg19 <- liftBuild(bs_small_hg38, current = "hg38", new = "hg19")
dim(bs_small_hg38)
dim(bs_hg19)
unique(genome(bs_hg19))

## ----discovery-subset---------------------------------------------------------
bs_demo <- bs[as.character(seqnames(bs)) %in% c("chr21", "chr22"), ]
sampleNames(bs_demo) <- pData(bs_demo)$"Sample ID" # guarantee unique sample IDs
bs_demo

## ----find-icrs----------------------------------------------------------------
icrs <- findICRs(bs_demo, minCov = 4, minSample = 0.5, ncores = 4)
head(icrs)

## ----annotate-icrs------------------------------------------------------------
icrs_annotated <- annotateICRs(icrs, genomeBuild = "hg38")
head(icrs_annotated)

## ----find-sems----------------------------------------------------------------
sample_sems <- findSEMs(
  bs_demo, minCov = 10, minSamples = 50,
  saveIndSEMs = TRUE, saveDir = tempdir(), ncores = 4
)
head(sample_sems)

## ----annotate-sems------------------------------------------------------------
load(file.path(tempdir(), "sampleSEMs.rda")) # loads cpgSEMs
sem_annotated <- annotateSEMs(
  cpgSEMs = cpgSEMs, sampleNames = sampleNames(bs_demo), genomeBuild = "hg38",
  universeCpGs = bs_demo
)
str(sem_annotated, max.level = 1)
head(sem_annotated$cpgGeneMap)

## ----test-sem-burden----------------------------------------------------------
group <- pData(bs_demo)$Sex
burden <- testSEMBurden(sample_sems, group = group)
burden$test
burden$summary

## ----test-sem-genes-----------------------------------------------------------
gene_test <- testSEMGenes(sem_annotated$geneMat, group = group)
head(gene_test$results)

## ----go-sems------------------------------------------------------------------
sig_genes <- gene_test$results$gene[gene_test$results$FDR < 0.25]
if (length(sig_genes) > 0) {
  go_res <- goSEMs(sig_genes, universe = sem_annotated$universeGenes, runKEGG = FALSE)
  head(as.data.frame(go_res$GO))
} else {
  message("No genes passed the demo FDR cutoff on this small chr21+chr22 subset -- ",
          "expected on a 2-chromosome slice; try the full genome on real data.")
}

## ----predict-sex--------------------------------------------------------------
sex_res <- predictSex(bs_small)
sex_res[, c("sample", "karyotype", "predictedSex")]

## ----predict-sex-calibrated---------------------------------------------------
sex_res_cal <- predictSex(bs_small, knownSex = pData(bs_small)$Sex)
attr(sex_res_cal, "calibrated")

## ----predict-race-------------------------------------------------------------
race_res <- predictRace(bs_small)
race_res[, c("sample", "predictedRace", "maxProb")]

## ----predict-pasc-------------------------------------------------------------
pasc_res <- predictPASC(bs_small)
pasc_res[, c("sample", "prob_enet", "predictedStatus")]

## ----predict-age--------------------------------------------------------------
age_res <- predictAge(bs_small)
age_res[, c("sample", "wgbs_only_stack", "hybrid_stack")]
pData(bs_small)$Age # known chronological ages, for comparison

## ----emulate-array, warning = FALSE-------------------------------------------
# warning = FALSE: methylCIPHER's own clock formulas emit benign vector-
# recycling warnings on small inputs; harmless, just noisy in a vignette
arr <- emulateArray(bs_small, array = "EPICv1", runClocks = TRUE)
dim(arr$beta)
arr$ageEstimates[, 1:6]

## ----session-info-------------------------------------------------------------
sessionInfo()

