The purpose of "BEASTmasteR" is to convert NEXUS data file(s) (DNA, amino acids,
discrete morphological characters, and/or continuous traits), plus an Excel
settings file, into Beast2 XML format.

The author is Nicholas J. Matzke, ORCID: http://orcid.org/0000-0002-8698-7656 .

BEASTmasteR is primarily aimed at enabling tip-dating analyses using fossils as
dated terminal taxa. The Birth-Death-Skyline-Serial-Sampling tree model is used
(or the Sampled-Ancestors variant), enabling different birth, death, and sampling
rates through time. However, it may be useful for setting up molecular-only
analyses (particularly automated variants), for learning Beast2 XML (since the
BEASTmasteR XML output is annotated and much easier for human reading than the
BEAUTi output), or for automating such analyses. BEASTmasteR also contains
functions for parsing Beast2 output, e.g. plotting a dated tree in R, with
posterior probabilities, 95% HPDs on node dates, and also 95% HPDs on tip dates,
if available.

2026 additions: The BEASTmasteR R code was updated to work on DNA and morphology datasets in the current Beast 2.7+, and to fix various problems due to the change in how the scan() function works in newer versions of R. (scan() now seems to hang on any moderately large text file; this is a memory and/or typing problem. The solution is to add e.g. nmax=1000000 as an argument. This was done throughout.). Note: The BDSS skyline tree prior works, but the Sampled-Ancestor prior needs more work at present.

2016-2017 additions: Added ascertainment bias corrections, various bug checks, 
gene-tree/species tree analysis setup.

<b>Citation</b>

BEASTmasteR will eventually become an R package and have a publication associated
with it. Until then, please cite something like:

Matzke, Nicholas J. (2014). "BEASTmasteR: R tools for automated conversion of
NEXUS data to BEAST2 XML format, for fossil tip-dating and other uses."
<i>PhyloWiki</i>, <a href="http://phylo.wikidot.com/beastmaster">http://phylo.wikidot.com/beastmaster</a> . Accessed (access
date).

Matzke, Nicholas J. (2026). BEASTmasteR code archive. <i>Github</i>:
<a href="https://github.com/nmatzke/BEASTmasteR">https://github.com/nmatzke/BEASTmasteR</a> . Accessed (access date). Release: 0.21. 
DOI: <a href="http://dx.doi.org/10.5281/zenodo.594056">http://dx.doi.org/10.5281/zenodo.594056</a>. 

Ideally, there will be a release with a DOI, but it may be you just use the 
most up-to-date commit. Find the most recent release at: 
https://github.com/nmatzke/BEASTmasteR/releases , and/or DOI 
http://dx.doi.org/10.5281/zenodo.594056 -- and a button: <a href="https://zenodo.org/badge/latestdoi/18687/nmatzke/BEASTmasteR"><img src="https://zenodo.org/badge/18687/nmatzke/BEASTmasteR.svg" alt="10.5281/zenodo.594056"></a>

<b>Acknowledgements</b>

Some of the functions in "tree_utils_v1.R", e.g. read_beast_prt, used for
extracting the bracketed statistics from BEAST NEXUS tree files (MCC files) are
copied/lightly modified/heavily modified from the R package " phyloch", by
Christoph Heibl, licensed under GPL (>=2), available at:
http://www.christophheibl.de/Rpackages.html. Please cite phyloch also, also if you
use this feature.
