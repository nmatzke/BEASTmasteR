
#######################################################
# Plot BEAST and MrBayes trees, auto-consensus, etc.
#######################################################


# Get plot coordinates from the last plotted phylogeny
get_phylo3_plotcoords <- function()
	{
	zz = get(x="last_plot.phylo", envir=.PlotPhyloEnv)
	
	x = zz$xx
	y = zz$yy

	xy = as.data.frame(cbind(x,y))
	return(xy)
	}



#######################################################
# Do an automatic FigTree-like PDF
# Read a Beast2 MCC tree, and add node bars, tip bars,
# posterior probs, etc.
#######################################################

plotMCC <- function(nexfn, titletxt="", pdffn=TRUE, tipnames_right_justified=TRUE, plot_node_heights=TRUE, plotPP=TRUE, minage=0, ladderize_tree=TRUE, fliptree=FALSE, tips_to_rotate=NULL, tipcex=1, italics_tiplabels=TRUE, digits=2, xmin=-5, xmax_mult=1.3, pdfheight=11, pdfwidth=9, newick=FALSE, tipdates_table=NULL, plot_tiplabels=TRUE, plot_axisPhylo=TRUE, tipbarcol=rgb(red=255,green=0,blue=0,alpha=100,max=255), nodebarcol=rgb(red=0,green=0,blue=255,alpha=100,max=255), bar_width=5, vlines=NULL, vline_col="black", vline_type="dotted", xtext="millions of years ago", space_tipnames=NULL, space_spaces=NULL)
	{
	defaults='
	library(BioGeoBEARS)	# for axisPhylo2()
	wd = "/drives/SkyDrive/NIMBioS_projects/2015-03-18_Tumamoc/doggies/r2/"
	setwd(wd)
	nexfn = "r2_treeLog.mcc"
	
	plotMCC(nexfn)
	
	titletxt=""
	pdffn=TRUE
	tipnames_right_justified=TRUE
	plot_node_heights = TRUE
	plotPP = TRUE
	minage = 0
	
	ladderize_tree = TRUE
	fliptree=FALSE
	tips_to_rotate = c("Canis_adustus", "Borophagus_diversidens")
	tipcex = 0.55
	italics_tiplabels=TRUE
	
	xmin=-5
	xmax_mult=1.3
	pdfheight=16
	pdfwidth=8.5
	newick=FALSE
	tipdates_table=NULL

	# This is a table with tipname, min_tipage, and max_tipage
	tipdates_table=NULL
	row1 = c("outgroup", 35, 45)
	row2 = c("Prohesperocyon_wilsoni", 32, 36)
	row3 = c("Borophagus_diversidens", 2, 5)
	tipdates_table = as.data.frame(rbind(row1, row2, row3), stringsAsFactors=FALSE, row.names=FALSE)
	names(tipdates_table) = c("tipname", "min_tipage", "max_tipage")
	tipdates_table
	
	
	# DEFAULTS
	plot_tiplabels=TRUE
	plot_axisPhylo=TRUE
	tipbarcol = rgb(red=255, green=0, blue=0, alpha=100, max=255)
	nodebarcol=rgb(red=150, green=150, blue=150, alpha=100, max=255)
	bar_width = 5
	vlines=NULL
	vline_col="black"
	vline_type="dotted"
	space_tipnames=NULL
	space_spaces=NULL
	
	' # END defaults

	#######################################################
	# Plot a Beast2 MCC tree with nice node bars, etc.
	#######################################################

	#pdffn = "treeLog.mcc.pdf"
	if (pdffn == TRUE)
		{
		pdffn = paste0(nexfn, ".pdf")
		}
	if (pdffn != FALSE)
		{
		pdf(pdffn, height=pdfheight, width=pdfwidth)
		}

	# Move the inner margins in some on the right (space for labels)
	# c(bottom, left, top, right) = c(5, 4, 4, 2) + 0.1.
	#par(mar=c(7,4,6,2))
	if (plot_axisPhylo == TRUE)
		{
		par(xaxs = "i")
		par(yaxs = "i") 
		}

	# Allow Newick if desired
	#nexfn = "treeLog.mcc"
	if (newick == FALSE)
		{
		beastcon = read_beast_prt(file=nexfn, digits=digits, get_tipnames=TRUE, printflag=FALSE)
		names(beastcon)
		dim(beastcon$prt_beast_nodestats)
		tr = beastcon$tr
		} else {
		tr = read.tree(nexfn)
		} # END if (newick == FALSE)
	
	
	if (ladderize_tree == TRUE)
		{
		ltr = ladderize(tr, right=FALSE)
		} else {
		ltr = tr
		} # END if (ladderize_tree == TRUE)

	# Get the MRCA node to rotate this one group
	if (!is.null(tips_to_rotate))
		{
		# Check if even
		if (length(tips_to_rotate)/2 != floor(length(tips_to_rotate)/2))
			{
			txt = "STOP ERROR in plotMCC(): 'tips_to_rotate' must have an EVEN number of tips, since they are processed in pairs."
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			} # END if (length(tips_to_rotate)/2 != floor(length(tips_to_rotate)))
		
		is_by_2 = seq(1, length(tips_to_rotate), by=2)
		for (i in is_by_2)
			{
			j = i+1
			node_to_rotate = getMRCA(phy=ltr, tip=c(tips_to_rotate[i],tips_to_rotate[j]))
			ltr = rotate(phy=ltr, node=node_to_rotate)
			} # END for (i in 1:(length(tips_to_rotate)/2))
		} # END if (!is.null(tips_to_rotate))
	
	if (fliptree == TRUE)
		{
		ltr = phytools::rotateNodes(tree=ltr, nodes="all")
		} # END if (fliptree == TRUE)

	# Relabel nodes etc.
	ltr = read.tree(file="", text=write.tree(phy=ltr, file=""))
	ltr_table = prt(ltr, printflag=FALSE, get_tipnames=TRUE)

	# Add to the xlims on the right
	tr_height = max(ltr_table$node_ht)
	xlims = c(xmin, tr_height*xmax_mult)
	ylims = c(0, length(ltr$tip.label)+1)

	# Needed, nodes not in correct order; also the edge numbers that change!
	tr_table = prt(tr, printflag=FALSE, get_tipnames=TRUE)
	convert_trtable_to_ltr_table_order = match(ltr_table$tipnames, tr_table$tipnames)
	
	if (newick == FALSE)
		{
		beastcon$prt_beast_nodestats = beastcon$prt_beast_nodestats[convert_trtable_to_ltr_table_order,]
		#head(cbind(tr_table$tipnames[convert_trtable_to_ltr_table_order], ltr_table$tipnames))
		}
	
	# Get node numbers
	tipnums = 1:length(tr$tip.label)
	intnums = (length(tr$tip.label)+1) : (length(tr$tip.label)+tr$Nnode)
	

	# Blank plot
	# Plot, no borders (bty="n"), no labels (xlab, ylab), no tick marks (xaxt, yaxt)
	#plot(1, 1, xlim=xlims, ylim=ylims, pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")

	# Tree plot on top
	#par(new=FALSE)
	#par(usr=c(xlims[1], xlims[2], ylims[1], ylims[2]))
	#tipcex = 0.55
	
	# Plot the tree, and tiplabels if desired
	if (tipnames_right_justified == TRUE)
		{
		# Tipnames plotted in clear (can't see them)
		tipcol = rgb(red=0, green=0, blue=0, alpha=0, maxColorValue=255)
		} else {
		# Tipnames plotted in black
		tipcol = "black"
		} # END if (tipnames_right_justified == TRUE)
		
	if (!is.null(space_tipnames))
		{
		for (i in 1:length(space_tipnames))
			{
			space_tipname = space_tipnames[i]
			TF = ltr$tip.label == space_tipname
			spaces = paste0(rep("_", times=space_spaces[i]), collapse="")
			newname = paste0(spaces, ltr$tip.label[TF])
			ltr$tip.label[TF] = newname
			}
		}		
	plotvals = plot(ltr, show.tip.label=plot_tiplabels, cex=tipcex, tip.color=tipcol, x.lim=xlims, y.lim=ylims)
	title(titletxt)
	
	
	

	
	# Get the node coordinates
	xy = get_phylo3_plotcoords()
	tr_height = get_max_height_tree(ltr)

	# Plot the tiplabels
	if (italics_tiplabels == TRUE)
		{
		tiplabels_to_plot = gsub(pattern="_", replacement=" ", x=ltr_table$label[tipnums])
		} else {
		tiplabels_to_plot = ltr_table$label[tipnums]
		} # END if (italics_tiplabels == TRUE)
	
	# Plot the tip labels at the right edge
	if (tipnames_right_justified == TRUE)
		{
		for (i in 1:length(tiplabels_to_plot))
			{
			if (italics_tiplabels == TRUE)
				{
				label_txt = paste0("substitute(italic('", tiplabels_to_plot[i], "'))")
				} else {
				label_txt = tiplabels_to_plot[i]
				} # END if (italics_tiplabels == TRUE)
	
			#text(x=tr_height, y=tipnums, labels=tiplabels_to_plot, pos=4, cex=tipcex)
			cmdstr = paste0("text(x=tr_height, y=tipnums[i], labels=", label_txt, ", pos=4, cex=tipcex)")
			eval(parse(text=cmdstr))
			} # END for (i in 1:length(tiplabels_to_plot))
		} # END if (tipnames_right_justified == TRUE)

	# Plot dotted lines up to present
	# tips_in_plot_order = ltr$tip.label
	#rows_in_ltr_table = match(x=tips_in_plot_order, table=ltr_table$label)
	if (tipnames_right_justified == TRUE)
		{
		tip_TF = ltr_table$node.type == "tip"
		x0 = tr_height-ltr_table$time_bp[tip_TF]
		x1 = rep(tr_height, length(x0))
		y0 = ltr_table$node[tip_TF]
		y1 = y0
		segments(x0=x0, x1=x1, y0=y0, y1=y1, col="gray60", lwd=0.5, lty="dotted")
		} # END if (tipnames_right_justified == TRUE)

	# Plot bars at nodes
	if ( (plot_node_heights==TRUE) && (newick==FALSE))
		{
		x0 = tr_height-beastcon$prt_beast_nodestats$"height_95%_HPD_MIN"
		x1 = tr_height-beastcon$prt_beast_nodestats$"height_95%_HPD_MAX"
		nodebarcol = nodebarcol
		segments(x0=x0, x1=x1, y0=xy$y, y1=xy$y, col=nodebarcol, lwd=bar_width)
		} # END if ( (plot_node_heights==TRUE) && (newick==FALSE))

	# Plot bars at tips
	if ( !is.null(tipdates_table) )
		{
		# Go through the rows in the tipdates table
		for (i in 1:nrow(tipdates_table))
			{
			tipname = tipdates_table$tipname[i]
			min_tipage = as.numeric(tipdates_table$min_tipage[i])
			max_tipage = as.numeric(tipdates_table$max_tipage[i])
			
			# Find the tipnum that matches, if one does
			matchval = match(x=tipname, table=ltr_table$label)
			nodenum = ltr_table$node[matchval]
			
			# skip if not found
			if (is.na(matchval))
				{
				next()
				} # END if (is.na(matchval))
			
			x0 = tr_height-min_tipage
			x1 = tr_height-max_tipage
			y0 = nodenum
			y1 = nodenum
			tipbarcol = tipbarcol
			segments(x0=x0, x1=x1, y0=y0, y1=y1, col=tipbarcol, lwd=bar_width)
			} # END for (i in 1:nrow(tipdates_table))
		} # END if ( !isnull(tipdates_table) )



	# Plot posterior probabilities, if desired
	if (plotPP == TRUE) 
		{
		if (newick == FALSE)
			{
			# Beast2 NEXUS MCC files
			nontip_TF = beastcon$prt_beast_nodestats$node.type != "tip"
			tip_TF = nontip_TF == FALSE
			PPtxt = beastcon$prt_beast_nodestats$posterior[nontip_TF]
			edges = ltr_table$parent_br[nontip_TF]
			edgelabels(text=PPtxt, edge=edges, adj=c(0.5,-0.35), frame="none", bg="none", cex=0.5)
			#tiplabels(text=tipnames_to_plot, tip=tipnums, adj=c(0.0,0.5), cex=0.7, frame="none", bg="none", font=3)
			} else {
			# Newick files that have node labels
			nontip_TF = ltr_table$node.type != "tip"
			tip_TF = nontip_TF == FALSE
			PPtxt = ltr_table$label[nontip_TF]
			edges = ltr_table$parent_br[nontip_TF]
			edgelabels(text=PPtxt, edge=edges, adj=c(0.5,-0.35), frame="none", bg="none", cex=0.5)
			} # END if (newick == FALSE)
		} # END if (plotPP == TRUE) 

	

	# x-axis has time
	# Use axisPhylo2() to round it:
	if (plot_axisPhylo == TRUE)
		{
		axis_at = axisPhylo2(side=1, minage=minage)
		mtext(text=xtext, side=1, line=2)

		if (!is.null(vlines))
			{
			vlines_xvals = tr_height - vlines
			print(tr_height)
			print(axis_at)
			print(vlines_xvals)
			abline(v=vlines_xvals, col=vline_col, lty=vline_type)
			} # END if (!is.null(vlines))
		} else {
		if (!is.null(vlines))
			{
			txt = "WARNING: vlines (vertical lines) have been specified, but this requires plot_axisPhylo=TRUE. You can plot lines manually with abline()."
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			warning(txt)
			} # END if (!is.null(vlines))
		} # END if (plot_axisPhylo == TRUE)
	
	if (pdffn != FALSE)
		{
		dev.off()
		cmdstr = paste0("open ", pdffn)
		cat("\n\n")
		cat("PDF generated.\n\nOpening with:\n\n'")
		cat(cmdstr)
		cat("'\n\n")
		system(cmdstr)
		} # END if (pdffn != FALSE)
	
	
	# Write the trees to Newick and NEXUS
	out_trfn = paste0(nexfn, "_asDisplayed.newick")
	write.tree(ltr, file=out_trfn)
	out_trfn = paste0(nexfn, "_asDisplayed.nexus")
	write.nexus(ltr, file=out_trfn, translate=TRUE)

	# Write the trees to Newick and NEXUS
	out_trfn = paste0(nexfn, "_asDisplayed2.newick")
	write.tree(phytools::rotateNodes(tree=ltr, nodes="all"), file=out_trfn)
	out_trfn = paste0(nexfn, "_asDisplayed2.nexus")
	write.nexus(phytools::rotateNodes(tree=ltr, nodes="all"), file=out_trfn, translate=TRUE)

	
	
	return(ltr)
	} # END plotMCC


#######################################################
# Flip the tree at all nodes
# (useful, since R displays trees backwards from
#  e.g. FigTree)
#######################################################
fliptree <- function(tr)
	{
	tr = phytools::rotateNodes(tree=tr, nodes="all")
	return(tr)
	}




#######################################################
# Get the tip heights from the logfile
# (this assumes you have logged the tip heights, and 
# they are called "ehight". This is the BEASTmasteR
# default.
#######################################################
get_tipdates_from_logfile <- function(logfn, tr, burnin_skipnum=500, sample_every=1)
	{
	defaults='
	wd = "/drives/GDrive/__github/BEASTmasteR/inst/extdata/LW12_plot/"
	setwd(wd)
	logfn = "traceLog.txt"
	tr = read.nexus("treeLog.mcc")
	burnin_skipnum=500
	sample_every=1
	'
	
	tdf = read.table(logfn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
	head(tdf)
	names(tdf)
	
	col_names = names(tdf)
	height_TF = substr(x=col_names, start=1, stop=6) == "height"
	sum(height_TF)
	colnums = 1:ncol(tdf)
	height_colnums = colnums[height_TF]
	
	# Row numbers to extract (post-burnin)
	rownums = seq(burnin_skipnum, nrow(tdf), by=sample_every)
	
	# Go through the height columns. If they correspond to tree tip
	# labels, calculate stats
	sampled_tipdates = NULL
	for (i in height_colnums)
		{
		# R converts
		# height(sp1) --> height.sp1.
		# extract just "sp1":
		
		col_name = col_names[i]
		startpos = nchar("height.") + 1
		endpos = nchar(col_name) - 1
		col_tipname = substr(x=col_name, start=startpos, stop=endpos)
		
		# Is this in the list of tips?
		# If so, get extract numbers
		TF = col_tipname %in% tr$tip.label
		if (TF == FALSE)
			{
			next()
			} # END if (TF == FALSE)

		# Extract the column
		cmdstr1 = paste0(col_tipname, " = tdf$", col_name, "[rownums]")
		eval(parse(text=cmdstr1))
		
		# Add to big table
		cmdstr2 = paste0("sampled_tipdates = cbind(sampled_tipdates, ", col_tipname, ")")
		eval(parse(text=cmdstr2))
		} # END for (i in height_colnums)
	dims = dim(sampled_tipdates)
	dims
	
	# If there is NO variation in tip dates, or no tip dates recorded, just return NULL
	if (is.null(dims))
		{
		res = NULL
		txt = "WARNING in get_tipdates_from_logfile(): No tipdates were found, the logfile needs to have columns with name 'height.[tipname]'. NOTE: BEASTmasteR does NOT log tip ages, unless they are free to vary. Returning 'res = NULL'."
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		warning(txt)
		return(res)
		}
	
	cat("\n\n")
	txt = paste0("Success: get_tipdates_from_logfile() extracted ", dims[2], " logged tipdates, with ", dims[1], " samples each.")
	cat(txt)
	cat("\n\n")
	
	# 
	
	# Summary states
	means = apply(X=sampled_tipdates, MARGIN=2, FUN=mean)
	sds = apply(X=sampled_tipdates, MARGIN=2, FUN=sd)
	medians = apply(X=sampled_tipdates, MARGIN=2, FUN=median)
	mins = apply(X=sampled_tipdates, MARGIN=2, FUN=min)
	maxs = apply(X=sampled_tipdates, MARGIN=2, FUN=max)
	lower95 = apply(X=sampled_tipdates, MARGIN=2, FUN=quantile, probs=0.025)
	upper95 = apply(X=sampled_tipdates, MARGIN=2, FUN=quantile, probs=0.975)
	
	summary_tipdates = cbind(means, sds, medians, mins, maxs, lower95, upper95)
	summary_tipdates = as.data.frame(summary_tipdates, stringsAsFactors=FALSE)
	summary_tipdates
	
	# Input for plotMCC
	tipdates_table = cbind(row.names(summary_tipdates), summary_tipdates$lower95, summary_tipdates$upper95)
	tipdates_table = as.data.frame(tipdates_table, stringsAsFactors=FALSE, row.names=NULL)
	names(tipdates_table) = c("tipname", "min_tipage", "max_tipage")
	tipdates_table$min_tipage = as.numeric(tipdates_table$min_tipage)
	tipdates_table$max_tipage = as.numeric(tipdates_table$max_tipage)
	
	res = NULL
	res$sampled_tipdates = sampled_tipdates
	res$summary_tipdates = summary_tipdates
	res$tipdates_table = tipdates_table

	tipdates_table

	extract='
	sampled_tipdates = res$sampled_tipdates
	summary_tipdates = res$summary_tipdates
	tipdates_table = res$tipdates_table
	' 
	return(res)
	} # END summarize_tipdates_from_logfile <- function(logfn, tr)




plotTrace <- function(logfn, pdffn=TRUE, trace_title="Trace plots", colnames_to_use=NULL, colnums_to_use=NULL, pdfheight=11, pdfwidth=8.5, subplot_basenum=4, subplot_bignum=subplot_basenum+2, numsamps=500)
	{
	defaults='
	library(BioGeoBEARS)	# for axisPhylo2()
	wd = "/drives/SkyDrive/NIMBioS_projects/2015-03-18_Tumamoc/doggies/mb1/"
	setwd(wd)
	logfn = "traceLog.txt"
	pdffn=TRUE
	
	trace_title = "Trace plots"
	colnums_to_use=NULL
	pdfheight=11
	pdfwidth=8.5
	subplot_basenum=4
	subplot_bignum=subplot_basenum + 2
	numsamps=500
	' # END defaults

	#######################################################
	# Plot a Beast2 MCC tree with nice node bars, etc.
	#######################################################

	#pdffn = "treeLog.mcc.pdf"
	if (pdffn == TRUE)
		{
		pdffn = paste0(logfn, ".pdf")
		}
	if (pdffn != FALSE)
		{
		pdf(pdffn, height=pdfheight, width=pdfwidth)
		}
	
	# Load and subset table
	logtable = read.table(file=logfn, header=TRUE, stringsAsFactors=FALSE, skip=1)
	if (is.null(colnums_to_use))
		{
		if (is.null(colnames_to_use))
			{
			colnums_to_use = 1:ncol(logtable)
			} else {
			colnums_to_use = match(x=names(logtable), table=colnames_to_use)
			colnums_to_use = colnums_to_use[is.na(colnums_to_use) == FALSE]
			}
		}
	logtable = logtable[, colnums_to_use]
	head(logtable)

	numcols = ncol(logtable) - 1
	#par(mfrow=c(numcols, 1))

	# Format the plot
	#subplot_basenum = 4
	#subplot_bignum = subplot_basenum + 2
	subplot_nums = NULL
	for (i in 2:ncol(logtable))
		{
		if (i == 2)
			{
			subplot_nums = c(subplot_nums, rep(i-1, subplot_bignum))
			}

		if ( (i != ncol(logtable)) && (i != 2) )
			{
			subplot_nums = c(subplot_nums, rep(i-1, subplot_basenum))
			}

		if (i == ncol(logtable))
			{
			subplot_nums = c(subplot_nums, rep(i-1, subplot_bignum))
			}
		}
	layout(matrix(subplot_nums, ncol=1))

	#trace_title = "mb3.2.5_mb2orig"
	#numsamps = 500
	rows_to_use = seq(1, nrow(logtable), length.out=numsamps)

	for (i in 2:ncol(logtable))
		{
		titletxt = names(logtable)[i]
		xvals = logtable[rows_to_use,1]
		yvals = logtable[rows_to_use,i] 

		# Remove 4% space
		par(xaxs = "i")
		par(yaxs = "i")
		
		#######################################################
		# Adjust axes to show everything better
		#######################################################
		# If the mean of the post-burnin is closer to max
		postburnin_samples = yvals[floor(length(yvals)/2) : (length(yvals))]
		dist_mean_to_max = max(yvals) - mean(postburnin_samples)
		dist_mean_to_min = mean(postburnin_samples) - min(yvals)
		
		if (dist_mean_to_max <= dist_mean_to_min)
			{
			# mean closer to the max
			# E.g. the first one is LnL which displays poorly
			ymax = max(yvals)
			ymin_diff = abs((ymax - min(postburnin_samples)))
			if (mean(postburnin_samples) > 0)
				{
				ymin = max(c(0, ymax - 7*ymin_diff))
				} else {
				ymin = max(c(min(yvals), ymax - 7*ymin_diff))
				}
			} else {
			# mean closer to the min
			ymin = min(yvals)
			ymax_diff = abs((max(postburnin_samples) - ymin))
			if (mean(postburnin_samples) > 0)
				{
				ymax = min(c(max(postburnin_samples), ymin + 7*ymax_diff))
				} else {
				ymax = min(c(max(postburnin_samples), ymin + 7*ymax_diff))
				} # END if (mean(postburnin_samples) > 0)
			} # END if (dist_mean_to_max >= dist_mean_to_min)

	
		# c(bottom, left, top, right)
		if (i == 2)
			{
			# c(bottom, left, top, right)
			par(mar=c(0,4,4,2))
		
			plot(xvals, yvals, pch=".", xlab="", xaxt="n", ylab=titletxt, ylim=c(ymin, ymax))
			title(trace_title)
			}

		if ( (i != ncol(logtable)) && (i != 2) )
			{
			par(mar=c(0,4,0,2))
			plot(xvals, yvals, pch=".", xlab="", xaxt="n", ylab=titletxt, ylim=c(ymin, ymax))
			}

		if (i == ncol(logtable))
			{
			par(mar=c(5,4,0,2))
			plot(xvals, yvals, pch=".", xlab="generation", ylab=titletxt, ylim=c(ymin, ymax))
			}
	
		# Plot the lines
		lines(logtable[,1], logtable[,i])
		} # END for (i in 2:ncol(logtable))


	

	if (pdffn != FALSE)
		{
		dev.off()
		cmdstr = paste0("open ", pdffn)
		cat("\n\n")
		cat("PDF generated.\n\nOpening with:\n\n'")
		cat(cmdstr)
		cat("'\n\n")
		system(cmdstr)
		} # END if (pdffn != FALSE)
	
	return(logtable)
	} # END plotTrace



#######################################################
# Merge two Tracer log files
#######################################################

merge_traceLogs <- function(fns, outfn=NULL, sample_num_col="Sample", pburnin=10)
	{
	defaults='

	fns = c(
	"traceLog_pt1.txt",
	"traceLog_pt2.txt"
	)

	outfn = "traceLog_pt12_last90.txt"
	sample_num_col = "Sample"
	pburnin = 10
	'

	runit = '
	#######################################################
	# Merge two Tracer log files
	#######################################################
	source("/drives/GDrive/__github/BEASTmasteR/R/plotBEAST_v1.R")
	wd = "/drives/GDrive/__GDrive_projects/2015-08-01_AFA_cladogram/_03_BEAST/2015-08-31/all_tips_fixed_WORKED/"
	setwd(wd)
	fns = c(
	"traceLog_pt1.txt",
	"traceLog_pt2.txt"
	)
	outfn = "traceLog_pt12_last90.txt"
	sample_num_col = "Sample"
	pburnin = 10
	newlog = merge_traceLogs(fns, outfn=outfn, sample_num_col="Sample", pburnin=10)
	head(newlog)
	'

	
	newlog = NULL

	cat("\nBEASTmasteR's merge_traceLogs() is reading files...\n")

	for (i in 1:length(fns))
		{
		fn = fns[i]
		txt = paste("\nReading file #", i, "/", length(fns), ": '", fn, "'", sep="")
		cat(txt)

		tmp = read.table(fn, sep="\t", header=TRUE, strip.white=TRUE, fill=FALSE)
		dim(tmp)
		head(tmp)
		
		# Strip all NA column
		if (all(is.na(tmp[,ncol(tmp)])) == TRUE)
			{
			tmp = tmp[,-ncol(tmp)]
			}
		
		# Subset
		fburnin = pburnin / 100
		startrow = ceiling(nrow(tmp) * fburnin)
		tmp2 = tmp[startrow:nrow(tmp),]
	
		newlog = rbind(newlog, tmp2)
		}

	# Re-label the sample rows
	# Get the column
	sample_col = newlog[,sample_num_col]
	increment = sample_col[3] - sample_col[2]
	minval = 0
	maxval = (length(sample_col) * increment)-1

	newlog[,sample_num_col] = as.numeric(seq(from=minval, to=maxval, by=increment))
	newlog[,sample_num_col]= as.numeric(newlog[,sample_num_col])

	# Write to outfn
	if (!is.null(outfn))
		{
		cat("\n\n...merge_traceLogs() is done reading files. Writing combined file to: '", outfn, "'...\n\n")

		# Turn off stupid scientific notation (crashes Tracer)
		# http://stackoverflow.com/questions/3978266/number-format-writing-1e-5-instead-of-0-00001
		default_option = getOption(x="scipen")
		options(scipen=10)
		write.table(x=newlog, file=outfn, sep="\t", quote=FALSE, row.names=FALSE)
		options(scipen=default_option)
		}
	
	return(newlog)
	}


