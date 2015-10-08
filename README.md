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

Citation

BEASTmasteR will eventually become an R package and have a publication associated
with it. Until then, please cite:

Matzke, Nicholas J. (2015). "BEASTmasteR: R tools for automated conversion of
NEXUS data to BEAST2 XML format, for fossil tip-dating and other uses."
Instructions at PhyloWiki, http://phylo.wikidot.com/beastmaster. Accessed (access
date).

Matzke, Nicholas J. (2015). BEASTmasteR code archive. Github:
https://github.com/nmatzke/BEASTmasteR . Accessed (access date). Release:
(Ideally, there will be a release with a DOI, but it may be you just use the most
up-to-date commit.)

Some of the functions in "tree_utils_v1.R", e.g. read_beast_prt, used for
extracting the bracketed statistics from BEAST NEXUS tree files (MCC files) are
copied/lightly modified/heavily modified from the R package " phyloch", by
Christoph Heibl, licensed under GPL (>=2), available at:
http://www.christophheibl.de/Rpackages.html. Please cite phyloch also, also if you
use this feature.
