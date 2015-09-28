
#######################################################
# BEASTmasteR: P(direct ancestry)
#######################################################

# Go through a list of post-burnin trees, asses
# probability of direct ancestry

prob_direct_ancestry <- function(trs, burnin=0, fraction=TRUE)
	{
	defaults='
	trs = read.nexus("/drives/GDrive/__GDrive_projects/2015-08-01_AFA_cladogram/_03_BEAST/2015-08-30/treeLog2.txt")
	burnin=0.5
	fraction=TRUE
	'
	
	# Fraction or number burnin
	if (fraction == TRUE)
		{
		burnin_number = 1+floor( (length(trs) * burnin) )
		} else {
		burnin_number = burnin
		}
	if (burnin_number < 1)
		{
		burnin_number = 1
		}
	
	# Do burnin
	trs2 = trs[burnin_number:length(trs)]
	direct_ancestor_TF = sapply(X=trs2, FUN=tips_w_edgelength_zero)
	direct_ancestor_TF = unlist_df(as.data.frame(direct_ancestor_TF, stringsAsFactors=FALSE))
	head(direct_ancestor_TF)
	
	# Sums_direct
	num_direct_ancestry = apply(X=direct_ancestor_TF, MARGIN=1, FUN=sum)
	fraction_direct_ancestry = num_direct_ancestry / ncol(direct_ancestor_TF)
	
	fraction_direct_ancestry
	
	return(fraction_direct_ancestry)
	} # END prob_direct_ancestry <- function(trs)


# Get the tips with edgelength == 0
tips_w_edgelength_zero <- function(tr)
	{
	defaults='
	tr = trs[[1000]]
	tr
	'
	ntips = length(tr$tip.label)
	tipnums = 1:ntips
	edge_below_tip_TF = tr$edge[, 2] <= ntips
	edge_below_tip_TF

	
	edgenums_below_tips = (1:nrow(tr$edge))[edge_below_tip_TF]
	branchlengths_below_tips = tr$edge.length[edgenums_below_tips]
	direct_ancestor = (branchlengths_below_tips == 0)
	nodenums_at_tips = tr$edge[, 2][edgenums_below_tips]
	
	# Put TF in tip.label order:
	order_tips = order(nodenums_at_tips)
	nodenums_at_tips[order_tips]
	direct_ancestor = as.data.frame(matrix(data=direct_ancestor[order_tips], nrow=1), row.names=NULL, stringsAsFactors=FALSE)
	names(direct_ancestor) = tr$tip.label
	direct_ancestor
	
	return(direct_ancestor)
	}


# Barplot for direct ancestry
barplot_direct_ancestry <- function(fraction_direct_ancestry, tipname_cex=0.45, tipname_line=0.15, xaxis_line=2.0, titletxt=NULL)
	{
	defaults = '
	tipname_cex = 0.45
	tipname_line = 0.15
	xaxis_line = 2.0
	titletxt = NULL
	'
	
	if (is.null(titletxt))
		{
		titletxt = "Probability of direct ancestry for each OTU"
		}
	
	par(mar=c(5,10,4,2))
	par(yaxs = "i")
	B = barplot(rev(fraction_direct_ancestry), names.arg=names(rev(fraction_direct_ancestry)), horiz=TRUE, xlim=c(0,1), axisnames=FALSE, col="grey90")
	label_vals = names(rev(fraction_direct_ancestry))
	mtext(text=label_vals, side=2, at=B, las=1, cex=tipname_cex, line=tipname_line)
	mtext(text="posterior probability", side=1, line=xaxis_line)
	abline(v=1.0, lty="dashed", col="grey50")
	title(titletxt)
	
	return(B)
	}


