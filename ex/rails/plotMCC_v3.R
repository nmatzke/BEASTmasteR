
#######################################################
# Plot figures, and posterior probability of direct ancestry
#######################################################

library(ape)	# for read/write NEXUS
library(BioGeoBEARS)	# for list2str, moref
library(gdata)	# for trim

source('~/GitHub/bioinfRhints/Rtnt_R_utils_v1.R')
source('/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_generics_v1.R')
sourceall('~/Documents/GitHub/BEASTmasteR/R/')

# wd = "/drives/GDrive/__GDrive_projects/2014-11-21_Randy_Irmis_autapomorphies/_06_TreeSim/_03_BEAST/01_Mk_on_strict_clock_alldata/"
# setwd(wd)



#######################################################
# Plot the FigTree-type tree
#######################################################
defaults = '
nexfn
pdffn=TRUE
plot_node_heights=TRUE
plotPP=TRUE
minage=0
ladderize_tree=TRUE
fliptree=FALSE
tips_to_rotate=NULL
tipcex=1
digits=2
xmin=-5
xmax_mult=1.3
pdfheight=11
pdfwidth=9
newick=FALSE
tipdates_table=NULL
xtext="millions of years ago"
vlines=NULL
vline_col="black"
vline_type="dotted"
space_tipnames=NULL
space_spaces=NULL
fonts=NULL
use_substitute=FALSE


show.tip.label=TRUE
default_tiplabels=TRUE
space_tipnames=NULL
space_spaces=NULL
tiplabel_underscores=FALSE
tipnames_right_justified=FALSE
italics_tiplabels=TRUE
tips_with_italics=NULL
fonts=NULL
use_substitute=FALSE

'

###########################################
# Inputs for the plot
###########################################
# CHANGE: change nexfn to match your MCC tree's filename. Also, check your
# CHANGE: working directory with getwd(), and set it with setwd("whatever/directory/you/like")
nexfn = "treeLog.mcc"    # MCC tree from TreeAnnotator (NEXUS format)

# CHANGE: Any of these options, as you see fit.
titletxt = "Canidae morph+DNA data+constraints (ALL); BEASTmasteR setup"
pdffn=TRUE              # If TRUE, output will be nexfn.pdf; 
                        # if FALSE, plot to screen; 
                        # anything else interpreted as output filename
tipnames_right_justified=TRUE   # Plot tipnames at 0 mya instead of at branch tips;
                                 # connected by dotted grey lines
plot_node_heights=TRUE  # Plot blue bars for node height variation
plotPP=TRUE       	    # Print posterior probabilities at middle of branches
minage=0	            # Should the right edge of the plot be 0 mya or some other value?
                        # NOTE: For calendar years, e.g. "2013 2014 2015 2016", 
                        # set minage=-2016.
ladderize_tree=TRUE     # ladderize the tree, e.g. as in FigTree "increasing"
fliptree=FALSE          # flip the whole tree? (e.g. FigTree "decreasing")
tips_to_rotate = NULL   # If you want to manually rotate some nodes, put the 
                        # tipnames here. The first 2 will specify the first node
                        # to rotate, the next 2 the next node, etc. 
                        # Obviously, you must have an even number of tipnames.
                        # e.g.: tips_to_rotate = c("sp1", "sp2")
tipcex = 0.75           # cex (character expansion) on the tipnames
italics_tiplabels=FALSE # When plotting tipnames, plot in italics?
digits=2                # Number of digits for e.g. Bayesian posterior probabilities
xmin=-0                 # Minimum of the x-axis (typically, xmin=0=root node age, which 
                        # cuts off the bar on the bottom/root node. If so, make 
                        # xmin more negative
xmax_mult=1.4           # Extend the x-axis maximum by multiplying the tree height by
                        # this value (helps avoid getting the tipnames cut off)
pdfheight=30             # PDF height (in inches) ('Merica!!)
pdfwidth=12              # PDF width (in inches) ('Merica!!)
newick=FALSE            # Set newick=TRUE if nexfn is Newick instead of NEXUS. Node bars
                        # not available for Newick unless you write your own function
tipdates_table=NULL     # If you want bars on the tipdates also, give the data.frame here.
                        # See formatting below. Use get_tipdates_from_logfile() to get it.
show.tip.label=TRUE     # Use standard APE labels for starting phylo plot, and xlimits (recommended). Default TRUE
default_tiplabels=FALSE	# If TRUE (default), default APE plot behavior for tiplabels
tiplabel_underscores=FALSE # If FALSE (default), underscores converted to spaces for non-default plots
tips_with_italics=NULL  # A list of tips that *should* be converted to italics
vlines=NULL				# Vertical lines, if desired
vline_col="grey50"      # Color of the vertical lines
vline_type="dashed"     # lty (line type) of the vertical lines
xtext="Millions of years ago"       # If you want something other than the default x-axis
                        # label. E.g., "mya", "years ago", "mega-annum", etc.
space_tipnames=NULL     # tipnames to which to add extra spaces
space_spaces=NULL       # number of spaces for each tipname

# Manual setup of tipdates_table:
# 
# row1 = c("outgroup", 35, 45)
# row2 = c("spA", 32, 36)
# row3 = c("spB", 2, 5)
# tipdates_table = as.data.frame(rbind(row1, row2, row3), stringsAsFactors=FALSE, row.names=FALSE)
# names(tipdates_table) = c("tipname", "min_tipage", "max_tipage")
# tipdates_table

# Automatic extraction of the tipdates from the logfile:
logfn = "traceLog.txt"    # Filename of the Beast2 log file
tr = read.nexus(nexfn)    # Tree object required (for the tipnames)
burnin_skipnum=1000        # Number of trees to skip as burnin
sample_every=1            # After excluding burnin, subsample if desired

# Run it:
res = get_tipdates_from_logfile(logfn, tr, burnin_skipnum, sample_every)

# Check out the summary!
res$summary_tipdates

# Extract from results (res)
sampled_tipdates = res$sampled_tipdates
summary_tipdates = res$summary_tipdates
tipdates_table = res$tipdates_table

write.table(summary_tipdates, file="summary_tipdates.txt", sep="	", quote=FALSE)

plotMCC(nexfn=nexfn, 
titletxt=titletxt, pdffn=pdffn, 
tipnames_right_justified=tipnames_right_justified, 
plot_node_heights=plot_node_heights, 
plotPP=plotPP, 
minage=minage, 
ladderize_tree=ladderize_tree, 
fliptree=fliptree, 
tips_to_rotate=tips_to_rotate, 
tipcex=tipcex, 
italics_tiplabels=italics_tiplabels, 
digits=digits,
xmin=xmin,
xmax_mult=xmax_mult,
pdfheight=pdfheight,
pdfwidth=pdfwidth,
newick=newick,
tipdates_table=tipdates_table,
show.tip.label=show.tip.label,
default_tiplabels=default_tiplabels,
tiplabel_underscores=tiplabel_underscores,
tips_with_italics=tips_with_italics,
vlines=vlines,
vline_col=vline_col,
vline_type=vline_type,
xtext=xtext,
space_tipnames=space_tipnames,
space_spaces=space_spaces)