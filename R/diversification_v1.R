

#######################################################
# Functions for diversification analysis on non-ultrametric trees
#######################################################
# Most algorithms used to study diversification, speciation/extinction rates, etc., 
# are designed for ultrametric trees, where all of the tips live in the Recent.

# If one has a tree with fossils included, however, these methods do not apply
# (e.g., most of the standard likelihood calculations condition on 
# the tips being living).
#
# Fixing this in a comprehensive way would require taking BDSS-type models, and
# modifying them to include the kinds of things that people want to study (SSE models, 
# diversity- and environment-dependent models, etc.)
#
# I do not attempt to solve those problems here, just do a few simple diversification
# analyses that can be done with paleo-trees, given existing methods.


#######################################################
# Method 1: symmetree-type analysis, of tree balance, 
# on a tree progressively "grown" 
# (or, in reverse, cut) from the root.
# 
# (We can use the BioGeoBEARS 'chainsaw' function to do
# the cutting)
# 
# Significant shifts in diversification (as measured by 
# imbalance in the tree topology) are recorded at each 
# time-step.
#
# This helps distinguish the times when the 
# diversification is actually happening, vs. later 
# diversification within the clade
#
# References:
# Moore (2004), SymmeTREE papers
# Ruta (2006), paleo method
# Lloyd (2008), paleo method applied to dinosaur supertree
# apTreeShape, which implements the basic symmetree statistics
#
#######################################################

example1 ='
## Detecting diversification rate variation in bird families (135 tips)
data(bird.families)
#bird.families = multi2di(bird.families)
tree.birds <- as.treeshape(bird.families, model = "yule")
class(tree.birds) <- "treeshape"


trshp = tree.birds
pv <- sapply(1:135, FUN=testshift2, trshp=trshp)

## Significant shifts detected at nodes = 67 and 78
pv[c(67,78)]
shift.test(tree.birds, node = 67, lambda1 = 1, lambda2 = 100, nrep = 10000, silent = TRUE)
shift.test(tree.birds, node = 78, lambda1 = 1, lambda2 = 100, nrep = 10000, silent = TRUE)
## visualize the shifts
par(mfrow=c(2,1))
plot(cutreeshape(tree.birds, ancestor(tree.birds, 67) , "bottom"))
plot(cutreeshape(tree.birds, 78 , "bottom"))
'

testshift2 <- function(i, trshp)
	{
	shift.test(trshp, i, lambda1 = 1, lambda2 = 100, nrep = 1000, silent = TRUE)
	}


#######################################################
# Modify apTreeshape::shift.test to provide more
# output
#######################################################
shift_test2 <- function (tree, node, lambda1=1, lambda2=100, nrep=1000, silent=FALSE, returnall=TRUE) 
	{
    if (class(tree)[1] != "treeshape") 
        stop("Error: object 'tree' not of class 'treeshape'")
    if ((node + 2) > length(tree$names)) 
        stop("Error: node > number of taxa - 2 ")
    obj <- Delta(tree, node, lambda1, lambda2)
    
    
    if (obj$clade.s[1, 1]%%2 == 0)
    	{
        L <- obj$clade.s[1, 1]/2
        l0 <- sapply(sample(1:L, nrep, replace = TRUE), FUN = function(l) {
            logratio(obj$clade.s[1, 1], l, lambda1, lambda2)
        	})
    	} else {
        L <- (obj$clade.s[1, 1] + 1)/2
        l0 <- sapply(sample(1:L, nrep, replace = TRUE, prob = c(rep(2, 
            L - 1), 1)), FUN = function(l) {
            logratio(obj$clade.s[1, 1], l, lambda1, lambda2)
        	})
    	} # END if (obj$clade.s[1, 1]%%2 == 0)
    
    # Probably: if more than 3 tips
    if (obj$clade.s[2, 1] > 3)
    	{
        if (obj$clade.s[2, 1]%%2 == 0)
        	{
            L <- obj$clade.s[2, 1]/2
            l1 <- sapply(sample(1:L, nrep, replace = TRUE), FUN = function(l) {
                logratio(obj$clade.s[2, 1], l, lambda1, lambda2)
            	})
        	} else {
            L <- (obj$clade.s[2, 1] + 1)/2
            l1 <- sapply(sample(1:L, nrep, replace = TRUE, prob = c(rep(2, 
                L - 1), 1)), FUN = function(l) {
                logratio(obj$clade.s[2, 1], l, lambda1, lambda2)
            	})
        	}
    	} else {
        l1 <- logratio(obj$clade.s[2, 1], 1, lambda1, lambda2)
        l1 <- rep(l1, nrep)
	    } # END if (obj$clade.s[2, 1] > 3)
    Pleft <- mean((l0 - l1) <= obj$delta)
    Pright <- mean((l0 - l1) >= obj$delta)
    Delta <- obj$delta
    P <- min(Pleft, Pright)

	if (returnall == FALSE)
		{
		res = P
		} else {
		res = c(node, lambda1, lambda2, nrep, L, mean(l0), mean(l1), mean(l0 - l1), Delta, Pleft, Pright, P)
		res = t(as.data.frame(res, stringsAsFactors=FALSE))
		names(res) = c("node", "lambda1", "lambda2", "nrep", "L", "mean(l0)", "mean(l1)", "mean(l0 - l1)", "Delta", "Pleft", "Pright", "P")
		}
    
    # Return results
    if (silent) 
    	{
    	return(res)
	    } else {
		cat("Test of diversification rate shift at node", node, 
			"\n")
		cat("Delta1 statistic = ", Delta, "\n")
		cat("P-value = ", P, "\n")
		cat("")
		cat("alternative hypothesis: the diversification rate shift is less than ", 
			lambda2/lambda1, " fold \nthe ancestral rate \n")
		cat("Note: The P-value was computed using ", nrep, " Monte-Carlo replicates. \n")
		return(res)
		} # END if (silent) 
	} # END shift_test2



#######################################################
# Calculate Moore (2004) Delta1 time statistic by timeslice
# After Ruta 2006, Lloyd 2008
#######################################################
# Returns stats for each node, and the relative (sub)
# and corresponding absolute (trfn) node numbers
# lambda2 amount of increase over lambda1 (which equals 1)
delta1_as_tree_grows <- function(trfn, timeperiods, lambda2=2)
	{
	defaults='
	trfn = "master_tree.newick"
	timeperiods = c(1,2,3,4,5,6,7,8,9,10,20)
	'
	
	master_tree = read.tree(trfn)
	master_trtable = prt(master_tree, printflag=FALSE, get_tipnames=TRUE)

	inputs = NULL
	inputs$timeperiods = timeperiods
	inputs$trfn = trfn
	inputs = section_the_tree(inputs, make_master_table=TRUE, plot_pieces=TRUE, cut_fossils=FALSE, fossils_older_than=0.1, save_phys_before_they_are_chopped=TRUE) 
	inputs$master_table
	inputs$phys_before_they_are_chopped
	
	
	# Correct the fossil branch lengths on phys_before_they_are_chopped
	# (since they have been artificially extended, temporarily)
	trs = inputs$phys_before_they_are_chopped
	for (i in 2:(length(trs)))
		{
		tmptr = trs[[i]]
		timeval = timeperiods[(i-1)]
		tmptable = prt(tmptr, printflag=FALSE, get_tipnames=TRUE)
		
		tips_are_OTUs_TF = grepl(pattern=",", x=tmptable$tipnames) == FALSE
		tipnames_to_edit = tmptable$label[tips_are_OTUs_TF]
		
		for (j in 1:length(tipnames_to_edit))
			{
			tipname_to_edit = tipnames_to_edit[j]
			TF = tmptable$label == tipname_to_edit
			
			# Find it on the Master Table
			mTF = master_trtable$label == tipname_to_edit
			# Find its time_bp
			time_bp_original = master_trtable$time_bp[mTF]
			# Modify by timeperiod
			time_bp_modified = time_bp_original - timeperiods[i-1]
			if (time_bp_modified < 0)
				{
				time_bp_modified = 0
				}
			
			# Tip time_bp on the chainsawed tmptable tree
			chainsawed_time_bp = tmptable$time_bp[TF]
			
			dif = chainsawed_time_bp - time_bp_modified
			if (dif < 0)
				{
				TF = tmptr$tip.label == tipname_to_edit
				nodenum = (1:length(tmptr$tip.label))[TF]
				edgeTF = tmptr$edge[,2] == nodenum
				edgenum = (1:nrow(tmptr$edge))[edgeTF]
				tmptr$edge.length[edgenum] = tmptr$edge.length[edgenum] + dif
				}
			}
		trs[[i]] = tmptr
		}
	# Back into object
	inputs$phys_before_they_are_chopped = trs

	# WHOLE TREE
	tr = master_tree
	timeperiod = 0

	# Convert tree to a apTreeshape object
	trtable = prt(tr, printflag=FALSE, get_tipnames=TRUE)

	# Convert tree to a apTreeshape object
	trshp = as.treeshape(tr, model = "yule")
	class(trshp)
	class(trshp) = "treeshape"
	nodes_to_test = 1:nrow(trshp$merge)
	# Exclude the root node
	nodes_to_test = nodes_to_test[-length(nodes_to_test)]

	# Get Delta1 stats for each node
	delta_table = NULL
	for (i in nodes_to_test)
		{
		#cat(i, "\n")
		res = shift_test2(trshp, nodes_to_test[i], lambda1=1, lambda2=lambda2, nrep=1000, silent=TRUE, returnall=TRUE)
		delta_table = rbind(delta_table, res)
		}

	delta_table = adf2(delta_table)
	names(delta_table) = c("node", "lambda1", "lambda2", "nrep", "L", "mean_l0", "mean_l1", "meanl0ml1", "Delta1", "Pleft", "Pright", "P")
	cls.df(delta_table)
	#head(delta_table)

	# Most significant shifts
	head(delta_table[order(delta_table$P),])

	# Find Delta1s < 0.1
	TF = delta_table$P <= 1

	# Get the node num that is significant
	subtable = delta_table[TF, ]

	if (sum(TF) > 0)
		{
		timeperiod_i = NULL
		subnodenums = NULL
		nodenums = NULL
		timeperiod = 0
		for (i in 1:nrow(subtable))
			{
			# Get the subtree above this node
			subtree = apTreeshape::cutreeshape(tree=trshp, apTreeshape::ancestor(tree=trshp, i=subtable$node[i]), type="bottom")
		
			# Unique identifier of clade
			subtree_tips_txt = paste0(sort(subtree$names), collapse=",")

			# Get the subnodenums
			node_match_TF = trtable$tipnames == subtree_tips_txt
			subnodenum = trtable$node[node_match_TF]
		
		
			# Split and re-sort to match with master tree
			words = strsplit(subtree_tips_txt, split=",")[[1]]
			subtree_tips_txt = paste0(sort(words), collapse=",")

			# Get the master tree nodenums
			node_match_TF = master_trtable$tipnames == subtree_tips_txt
			nodenum = master_trtable$node[node_match_TF]
		
			timeperiod_i = c(timeperiod_i, timeperiod)
			subnodenums = c(subnodenums, subnodenum)
			nodenums = c(nodenums, nodenum)
			} # END for (i in 1:nrow(subtable))

		subtable = cbind(subtable, subnodenums, nodenums, timeperiod_i)	
		} # END if (sum(TF) > 0)
	subtable


	subtables = NULL
	subtables = rbind(subtables, subtable)
	cat("\n\nCalculating Delta1 (Moore et al. 2004) by node as the tree grows (Ruta 2006, Lloyd 2008) for:\n")
	for (ti in 2:length(inputs$timeperiods))
		{
		cat("timeperiod ", ti-1, "/", length(inputs$timeperiods), "\n")
		tr = inputs$phys_before_they_are_chopped[[ti]]
		timeperiod = inputs$timeperiods[(i-1)]

		# Make the tree table
		trtable = prt(tr, printflag=FALSE, get_tipnames=TRUE)

		# Convert tree to a apTreeshape object
		trshp = as.treeshape(tr, model = "yule")
		class(trshp)
		class(trshp) = "treeshape"
		nodes_to_test = 1:nrow(trshp$merge)
		# Exclude the root node
		nodes_to_test = nodes_to_test[-length(nodes_to_test)]

		# Get Delta1 stats for each node
		delta_table = NULL
		for (i in nodes_to_test)
			{
			#cat(i, "\n")
			res = shift_test2(trshp, nodes_to_test[i], lambda1=1, lambda2=lambda2, nrep=1000, silent=TRUE, returnall=TRUE)
			delta_table = rbind(delta_table, res)
			}

		delta_table = adf2(delta_table)
		names(delta_table) = c("node", "lambda1", "lambda2", "nrep", "L", "mean_l0", "mean_l1", "meanl0ml1", "Delta1", "Pleft", "Pright", "P")
		cls.df(delta_table)
		#head(delta_table)

		# Most significant shifts
		head(delta_table[order(delta_table$P),])

		# Find Delta1s < 0.1
		TF = delta_table$P <= 1

		# Get the node num that is significant
		subtable = delta_table[TF, ]

		if (sum(TF) > 0)
			{
			timeperiod_i = NULL
			subnodenums = NULL
			nodenums = NULL
			timeperiod = inputs$timeperiods[(ti-1)]
			for (i in 1:nrow(subtable))
				{
				# Get the subtree above this node
				subtree = apTreeshape::cutreeshape(tree=trshp, apTreeshape::ancestor(tree=trshp, i=subtable$node[i]), type="bottom")

				# Unique identifier of clade
				subtree_tips_txt = paste0(sort(subtree$names), collapse=",")
		
				# Get the subnodenums
				node_match_TF = trtable$tipnames == subtree_tips_txt
				subnodenum = trtable$node[node_match_TF]
			
			
				# Split and re-sort to match with master tree
				words = strsplit(subtree_tips_txt, split=",")[[1]]
				subtree_tips_txt = paste0(sort(words), collapse=",")

				# Get the master tree nodenums
				node_match_TF = master_trtable$tipnames == subtree_tips_txt
				nodenum = master_trtable$node[node_match_TF]
			
				timeperiod_i = c(timeperiod_i, timeperiod)
				subnodenums = c(subnodenums, subnodenum)
				nodenums = c(nodenums, nodenum)
				} # END for (i in 1:nrow(subtable))

			subtable = cbind(subtable, subnodenums, nodenums, timeperiod_i)	
			} # END if (sum(TF) > 0)
		subtables = rbind(subtables, subtable)
		} # END for (ti in 1:length(inputs$timeperiods))
	
	res = NULL
	res$subtables = subtables
	res$phys_before_they_are_chopped = inputs$phys_before_they_are_chopped
	
	extract='
	subtables = res$subtables
	phys_before_they_are_chopped = res$phys_before_they_are_chopped
	'
	
	return(res)
	} # END delta1_as_tree_grows <- function(trfn, timeperiods)



#######################################################
# Get Delta1 stats for each node
# (following Moore 2004, Ruta et al. 2007, Lloyd 2008)
#######################################################

#######################################################
# Plot possible diversification shifts
#######################################################
plot_Delta1_diversification_shifts <- function(trfn, timeperiods=c(1,2), lambda2=10, minage=0, pdffn="Delta1s_by_time_period.pdf", width=7, height=12, mtext_txt="Ma")
	{
	defaults='
	trfn = "master_tree.newick"
	timeperiods = c(2,4,6,8,10,20)
	lambda2=10
	minage=-2016
	pdffn = "Delta1s_by_time_period.pdf"
	width=7
	height=12
	mtext_txt="Ma"
	'

	#######################################################
	# Run the diversification analysis
	#######################################################
	res = delta1_as_tree_grows(trfn, timeperiods, lambda2=lambda2)
	names(res)

	subtables = res$subtables
	phys_before_they_are_chopped = res$phys_before_they_are_chopped
	
	# Most significant shifts
	head(subtables[order(subtables$P),])


	
	#######################################################
	# Make the PDF
	#######################################################	
	pdf(file=pdffn, width=width, height=height)

	if (minage < 0)
		{
		tree_top_date = abs(minage)
		axisPhylo2_minage = minage
		} else {
		tree_top_date = minage
		axisPhylo2_minage = minage
		}


	master_tree = read.tree(trfn)
	orig_master_tree = master_tree
	orig_master_trtable = prt(orig_master_tree, printflag=FALSE, get_tipnames=TRUE)

	i = 0

	plot(master_tree, cex=0.75)
	axisPhylo2(minage=axisPhylo2_minage)
	vert_lines = get_max_height_tree(master_tree) - timeperiods
	abline(v=vert_lines, lty="dotted", col="grey60")
	abline(v=max(vert_lines), lty="dotted", col="black")
	legend_txt = paste0("Tree top at ", axisPhylo2_minage)
	legend(x="topleft", legend="Delta1 p-values", fill="yellow")
	title("Delta1 analysis for diversification shifts: whole tree")
	mtext(text=mtext_txt, side=1, line=2.5)

	TF1 = subtables$timeperiod_i == i
	TF2 = subtables$P <= 0.1
	TF = (TF1 + TF2) == 2
	sum(TF)

	nodes_to_plot = subtables$nodenums[TF]
	master_trtable = prt(master_tree, printflag=FALSE, get_tipnames=FALSE)

	edges_to_plot = master_trtable$parent_br[nodes_to_plot]
	if (sum(TF) > 0)
		{
		#edgelabels(text=subtables$P[TF], edge=edges_to_plot, adj=c(0.5,-0.35), frame="none", bg="none", cex=0.5)
		nodelabels(text=subtables$P[TF], node=nodes_to_plot, bg="yellow")#, frame="none", bg="grey80", cex=0.5)
		} # END if (sum(TF) > 0)



	# Plot by tree section
	par(mfrow=c(2,1))
	master_tree = orig_master_tree
	for (i in 1:(length(timeperiods)-1))
		{
		master_tree = orig_master_tree
		if (minage < 0)
			{
			tree_top_date = abs(minage) - timeperiods[i]
			axisPhylo2_minage = minage + timeperiods[i]
			} else {
			tree_top_date = minage + timeperiods[i]
			axisPhylo2_minage = minage + timeperiods[i]
			}
		tmp_timeperiods = c(0, timeperiods)
		vert_line_max = get_max_height_tree(master_tree) - timeperiods[i]
		
		
		# Plot the full tree
		plot(master_tree, show.tip.label=TRUE, cex=0.5)
		axisPhylo2(minage=minage)
		vert_lines = get_max_height_tree(master_tree) - tmp_timeperiods
		abline(v=vert_lines, lty="dotted", col="grey60")
		abline(v=vert_line_max, lty="dashed", col="black")
		legend_txt = paste0("Tree top pruned at: ", tree_top_date)
		legend(x="topleft", legend=legend_txt, border="white", col="black", lty="dashed")
		title("Diversification shift p-values (plotted on full tree)")
		mtext(text=mtext_txt, side=1, line=2.5)
	
		TF1 = subtables$timeperiod_i == timeperiods[i]
		TF2 = subtables$P <= 0.1
		TF = (TF1 + TF2) == 2
		sum(TF)
	
		#nodes_to_plot = subtables$subnodenums[TF]
		nodes_to_plot = subtables$nodenums[TF]
		master_trtable = prt(master_tree, printflag=FALSE, get_tipnames=FALSE)
		edges_to_plot = master_trtable$parent_br[nodes_to_plot]
		if (sum(TF) > 0)
			{
			#edgelabels(text=subtables$P[TF], edge=edges_to_plot, adj=c(0.5,-0.35), frame="none", bg="none", cex=0.5)
			nodelabels(text=subtables$P[TF], node=nodes_to_plot, bg="yellow")
			} # END if (sum(TF) > 0)


		# Plot the cut tree (subtree)
		master_tree = phys_before_they_are_chopped[[i+1]]
		master_tree_full_labels = master_tree
		#master_tree = master_tree
		master_tree$tip.label = substr(master_tree$tip.label, start=1, stop=15)
		plot(master_tree, show.tip.label=TRUE, cex=0.5)
		axisPhylo2(minage=axisPhylo2_minage)
		times_cut_from_pruned_tree = (timeperiods-timeperiods[i])
		times_cut_from_pruned_tree = times_cut_from_pruned_tree[times_cut_from_pruned_tree >= 0]
		vert_lines = get_max_height_tree(master_tree) - times_cut_from_pruned_tree
		abline(v=vert_lines, lty="dotted", col="grey60")
		abline(v=max(vert_lines), lty="dashed", col="black")
		legend_txt = paste0("Tree top pruned at: ", tree_top_date)
		legend(x="topleft", legend=legend_txt, border="white", col="black", lty="dashed")
		title("Diversification shift p-values (plotted on pruned tree)")
		mtext(text=mtext_txt, side=1, line=2.5)

		TF1 = subtables$timeperiod_i == timeperiods[i]
		TF2 = subtables$P <= 0.1
		TF = (TF1 + TF2) == 2
		sum(TF)
	
		nodes_to_plot = subtables$subnodenums[TF]
		#nodes_to_plot = subtables$nodenums[TF]
		master_trtable = prt(master_tree, printflag=FALSE, get_tipnames=FALSE)
		edges_to_plot = master_trtable$parent_br[nodes_to_plot]
		if (sum(TF) > 0)
			{
			#edgelabels(text=subtables$P[TF], edge=edges_to_plot)
			nodelabels(text=subtables$P[TF], node=nodes_to_plot, bg="yellow")
			} # END if (sum(TF) > 0)

		} # END for (i in 1:length(timeperiods[-length(timeperiods)]))
	dev.off()
	cmdstr = paste0("open ", pdffn)
	system(cmdstr)
	
	return(res)
	} # END plot_Delta1_diversification_shifts <- function(trfn, timeperiods, lambda2, pdffn, width=7, height=12)

