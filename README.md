# MethylTraits
MethylTraits: an R package to characterize and predict biological traits and features from whole genome bisulfite/methylation sequencing data

## Formerly known as oliveR (we're sad to see you go).
New recipe, same great taste!!!

## Welcome
Hey, thanks for stopping on by!

The package is pretty straightforward; all you really need to start is a bsseq object.

With this pacakge, you can do a few different things for which there are microarray functionalities for. For example, estimating age, sex, and stochastic epigenetic mutations (SEMs). Microarrays have been around for awhile and sample sizes and datasets are plentiful. The same cannot be said about WGMS/WGBS datasets. But, as costs are coming down to sequence DNA, I would like to believe that WGMS/WGBS datasets will start growing in number and size.

As such, this R packages seeks to bridge that gap by allowing those array-based utilities to be functional with sequencing data. Currently, starting with just a bsseq object, you can estimate sample ages, identify stochastic epigenetic mutations (SEMs), find putative imprinted control regions (ICRs), guess the sex of samples, guess race/ethnicity (is that bioethical? I don't know), estimate telomere length, and more!

## Vignette

Check out the tutorial of the general functionality of the package at:
[Here](https://htmlpreview.github.io/?https://raw.githubusercontent.com/andymadrid/MethylTraits/refs/heads/main/inst/doc/MethylTraits.html)

## Stochastic epigenetic mutations (SEMs)

While the main function of this package is to estimate age from sequencing-based data, others have recently been moving to detecting stochastic epigenetic mutations (SEMs) in DNA methylation data. So, I've implemented a function to do just that here. It calculates the interquartile range (IQR) for each CpG, then goes sample by sample, CpG by CpG to identify SEMs. For more information on SEMs, please refer to Gentilini et al (2015) Aging and Markov et al (2024) GeroScience. I recommend larger sample sizes (>50) and higher coverage (>10) to reduce the number of CpGs analyzed and to only keep those with substantial evidence of their estimated methylation levels.

## Word of warning

### Coverage
The biggest issue that I have run into while developing this package has been coverage. For reliable, robust results, higher coverage from sequencing data is needed. While testing, I found several datasets that had ~1-4x genome-wide coverage. As such, estimates for age, sex, SEMs either would be highly inaccurate, or were unable to be computed. This is one of the major drawbacks of sequencing data, relative to the array where most samples have useable data on just about every CpG tested on the array. Several imputation methods are utilized in the functions throuhgout the package, but those can only do so much. So, just be wary of that fact. As with most sequencing-based analyses: higher coverage gives better, more robust, and more accurate/reliable results.

### Genome reference
Samples used to test and validate this clock were aligned to the human genome (hg38) using the UCSC chromosome naming scheme (e.g,. chr1, chr2, chr3, etc). If you happen to have aligned your samples to, say, hg19, those coordinates will have to be lifted over in order to properly get filtered/selected for the clock to work. I added in a function liftBuild() that does just that.  

## Questions, concerns, collaborations
Do you have your own sequencing-based methylation data that you want to build your own clock with? Something else? If so, feel free to email me at at madrid2[at]wisc.edu and we can work it out, together. I am always happy to collaborate and-or help! Also, if there's some functionality that you think would be of great value to add, you can also let me know and I can work on implementing that, as well.

See you space cowboy...
- AM
