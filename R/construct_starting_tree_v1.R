
# Write a tree, without the scientific notation
write_tree_noSci <- function (phy, file = "", append = FALSE, digits = 10, tree.names = FALSE, sprintf_option="f") 
{
    if (!(inherits(phy, c("phylo", "multiPhylo")))) 
        stop("object \"phy\" has no trees")
    if (inherits(phy, "phylo")) 
        phy <- c(phy)
    N <- length(phy)
    res <- character(N)
    if (is.logical(tree.names)) {
        if (tree.names) {
            tree.names <- if (is.null(names(phy))) 
                character(N)
            else names(phy)
        }
        else tree.names <- character(N)
    }
    for (i in 1:N) res[i] <- write_tree3(phy[[i]], digits = digits, 
        tree.prefix = tree.names[i], sprintf_option=sprintf_option)
    if (file == "") 
        return(res)
    else cat(res, file = file, append = append, sep = "\n")
}

# ape:::.write.tree2 forces a "g" in the sprintf statement; I like "f"
# (and so does e.g. Beast2)
write_tree3 <- function (phy, digits = 10, tree.prefix = "", sprintf_option="g") 
{
    brl <- !is.null(phy$edge.length)
    nodelab <- !is.null(phy$node.label)
    phy$tip.label <- checkLabel(phy$tip.label)
    if (nodelab) 
        phy$node.label <- checkLabel(phy$node.label)
    f.d <- paste("%.", digits, sprintf_option, sep = "")
    cp <- function(x) {
        STRING[k] <<- x
        k <<- k + 1
    }
    add.internal <- function(i) {
        cp("(")
        desc <- kids[[i]]
        for (j in desc) {
            if (j > n) 
                add.internal(j)
            else add.terminal(ind[j])
            if (j != desc[length(desc)]) 
                cp(",")
        }
        cp(")")
        if (nodelab && i > n) 
            cp(phy$node.label[i - n])
        if (brl) {
            cp(":")
            cp(sprintf(f.d, phy$edge.length[ind[i]]))
        }
    }
    add.terminal <- function(i) {
        cp(phy$tip.label[phy$edge[i, 2]])
        if (brl) {
            cp(":")
            cp(sprintf(f.d, phy$edge.length[i]))
        }
    }
    n <- length(phy$tip.label)
    parent <- phy$edge[, 1]
    children <- phy$edge[, 2]
    kids <- vector("list", n + phy$Nnode)
    for (i in 1:length(parent)) kids[[parent[i]]] <- c(kids[[parent[i]]], 
        children[i])
    ind <- match(1:max(phy$edge), phy$edge[, 2])
    LS <- 4 * n + 5
    if (brl) 
        LS <- LS + 4 * n
    if (nodelab) 
        LS <- LS + n
    STRING <- character(LS)
    k <- 1
    cp(tree.prefix)
    cp("(")
    getRoot <- function(phy) phy$edge[, 1][!match(phy$edge[, 
        1], phy$edge[, 2], 0)][1]
    root <- getRoot(phy)
    desc <- kids[[root]]
    for (j in desc) {
        if (j > n) 
            add.internal(j)
        else add.terminal(ind[j])
        if (j != desc[length(desc)]) 
            cp(",")
    }
    if (is.null(phy$root.edge)) {
        cp(")")
        if (nodelab) 
            cp(phy$node.label[1])
        cp(";")
    }
    else {
        cp(")")
        if (nodelab) 
            cp(phy$node.label[1])
        cp(":")
        cp(sprintf(f.d, phy$root.edge))
        cp(";")
    }
    paste(STRING, collapse = "")
}


# Fix a tree where tips don't come to zero
#######################################################
# level_tree_tips
#######################################################
level_tree_tips <- function(tr, method="mean", printflag=TRUE, fossils_older_than=0.6)
	{
	defaults='
	method="mean"
	printflag=TRUE
	fossils_older_than=0.6
	'
	# Look at the tree table:
	trtable = prt(tr, printflag=printflag, fossils_older_than=fossils_older_than)
	trtable$time_bp

	ntips = length(tr$tip.label)
	fossils_TF = trtable$fossils[1:ntips]

	tipnums_to_change = (1:ntips)[fossils_TF == FALSE]
	tip_ages = trtable$time_bp[tipnums_to_change]
	
	# Get the mean of the tipages, add difference
	if (method == "mean")
		{
		mean_tipage = mean(tip_ages)
		}
	# Adjust everything to match the highest tip
	if (method == "highest")
		{
		mean_tipage = 0
		}
	# Adjust everything to match the lowest tip
	# (could introduce errors)
	if (method == "lowest")
		{
		mean_tipage = max(tip_ages)
		}


	change_branchlength_by = tip_ages - mean_tipage

	# Edit the tree:
	edgenums = trtable$parent_br[tipnums_to_change]
	tr$edge.length[edgenums] = tr$edge.length[edgenums] + change_branchlength_by
	
	# Check for negative branchlengths
	brlen_negative_count = sum( tr$edge.length[edgenums] < 0)
	
	if (brlen_negative_count > 0)
		{
		error_txt = paste0("STOP ERROR in level_tree_tips(): correcting uneven tip ages with method='", method, "', and fossils_older_than=", fossils_older_than, " has resulted in a tree with ", brlen_negative_count, " branches with negative branchlengths. This is Very Bad. You should correct your tree with another method or by hand.  See e.g. impose_min_brlen().")
		cat("\n\n")
		cat(error_txt)
		cat("\n\n")
		stop(error_txt)
		} # END if (brlen_negative_count > 0)
	
	return(tr)
	}





impose_min_brlen_OLD <- function(phy, min_brlen=0.01, leave_BL0_terminals=TRUE)
	{
	num_internal_nodes = phy$Nnode
	numtips = length(phy$tip.label)
	rootnodenum = numtips+1
	# likelihoods are computed at all nodes
	# make a list to store the 
	numnodes = numtips + num_internal_nodes
	
	# Reorder the edge matrix into pruningwise order
	# This is CRUCIAL!!
	phy2 <- reorder(phy, "pruningwise")
	phy2_orig = phy2
	
	observed_min_orig = min(phy2_orig$edge.length)
	TF = phy2_orig$edge.length < min_brlen
	num_branches_below_min = sum(TF)
	if (observed_min_orig >= min_brlen)
		{
		txt = paste0("No branches found =< min_brlen (", min_brlen, "). Returning original tree, pruningwise reordered.")
		cat("\n")
		cat(txt)
		cat("\n")
		return(phy2)
		} else {
		
		txt = paste0(num_branches_below_min, " branches found < min_brlen (", min_brlen, "). Running downpass to edit branchlengths by pushing nodes down so that minimum branchlength=", min_brlen, ". Adjusting i:Leftnode,j:Rightnode,anc:Ancnode;...")
		cat("\n")
		cat(txt)
		cat("\n\n")
		} 
		
		# END if (observed_min_orig > min_brlen)



	# DEFINE DOWNPASS THROUGH THE BRANCHES	
	tipnums <- 1:numtips
	i = 1
	edges_to_visit = seq(from=1, by=2, length.out=num_internal_nodes)

	#######################################################
	#######################################################
	# THIS IS A DOWNPASS FROM THE TIPS TO THE ROOT
	#######################################################
	#######################################################
	for (i in edges_to_visit)
		{
		# First edge visited is i
		#print(i)
		
		# Its sister is j 
		j <- i + 1
		#print(j)

		# Get the node numbers at the tips of these two edges		
		left_desc_nodenum <- phy2$edge[i, 2]
		right_desc_nodenum <- phy2$edge[j, 2]
		
		# And for the ancestor edge (i or j shouldn't matter, should produce the same result!!!)
		ancnode <- phy2$edge[i, 1]
		anc_edge_TF = phy2$edge[,2] == ancnode
		anc_edgenum = (1:nrow(phy2$edge))[anc_edge_TF]
		if (length(anc_edgenum) == 0)
			{
			anc_edgenum = NA
			}

		# Get the edge length (left, right, anc)
		Llength = phy2$edge.length[i]
		Rlength = phy2$edge.length[j]
		Alength = phy2$edge.length[anc_edgenum]
		
		txt = paste("ancnode:", ancnode, " left:", left_desc_nodenum, " right:", right_desc_nodenum, sep="")
		#print(txt)
		
		# If 0-length terminal branches are allowed (e.g. a 
		# sampled-ancestor tree), and
		# If a node is a terminal node,
		# and if it is of length 0 and the other branch > 0,
		# then don't change this branch length
		edit_brlength = TRUE
		if (leave_BL0_terminals == TRUE)
			{
			if (i <= numtips)
				{
				if ( (Llength == 0) && (Rlength > 0) )
					{
					edit_brlength = FALSE
					}
				} # END if (i <= numtips)
			if (j <= numtips)
				{
				if ( (Rlength == 0) && (Llength > 0) )
					{
					edit_brlength = FALSE
					}
				} # END if (i <= numtips)
			} # END if (leave_BL0_terminals == TRUE)
		
		
		# If either the left or right branch is < min_brlen, change both
		if (edit_brlength == TRUE) {
		if ( (Llength < min_brlen) || (Rlength < min_brlen) )
			{
			cat(i, ":", left_desc_nodenum, ",", j, ":", right_desc_nodenum, ",anc:", ancnode, "; ", sep="")
			
			amount_to_add_L = min_brlen - Llength
			amount_to_add_R = min_brlen - Rlength
			amount_to_add = max(c(amount_to_add_L, amount_to_add_R))
			
			# Change both branchlengths by the same amount
			phy2$edge.length[i] = phy2$edge.length[i] + amount_to_add
			phy2$edge.length[j] = phy2$edge.length[j] + amount_to_add
			
			# Subtract that amount from the branch below
			# (may result in negative, which will be fixed further in the downpass)
			if (is.na(anc_edgenum) == FALSE)
				{
				phy2$edge.length[anc_edgenum] = phy2$edge.length[anc_edgenum] - amount_to_add
				}
			}
			} # END if (edit_brlength == TRUE)
		
		} # end downpass
	#######################################################
	# END DOWNPASS FROM THE TIPS TO THE ROOT
	#######################################################
	rootnode = ancnode
	phy3 = phy2

	cat("\n")
	cat("\n")
	txt = paste0("Originally, ", num_branches_below_min, " branches were found < min_brlen (", min_brlen, "). Running downpass to edit branchlengths by pushing nodes down so that minimum branchlength=", min_brlen, ". Adjusting i:Leftnode,j:Rightnode,anc:Ancnode;...")
	cat("\n")
	cat(txt)
		
	observed_min_brlen_new = min(phy3$edge.length)
	cat("\n\n")
	cat("New observed_min_brlen_new=", observed_min_brlen_new, sep="")
	cat("\n")
	return(phy3)
	} # END impose_min_brlen <- function(phy, min_brlen=0.1)







# Avoid negative branchlengths

impose_min_brlen <- function(phy, min_brlen=0.01, leave_BL0_terminals=TRUE, direct_ancestor_brlen=1e-07)
	{
	defaults='
	min_brlen = 1e-6
	leave_BL0_terminals=TRUE
	direct_ancestor_brlen=1e-07
	'
	
	num_internal_nodes = phy$Nnode
	numtips = length(phy$tip.label)
	rootnodenum = numtips+1
	# likelihoods are computed at all nodes
	# make a list to store the 
	numnodes = numtips + num_internal_nodes
	
	# Reorder the edge matrix into pruningwise order
	# This is CRUCIAL!!
	phy2 <- reorder(phy, "pruningwise")
	phy2_orig = phy2
	
	observed_min_orig = min(phy2_orig$edge.length)
	TF = phy2_orig$edge.length < min_brlen
	num_branches_below_min = sum(TF)
	if (observed_min_orig >= min_brlen)
		{
		txt = paste0("No branches found =< min_brlen (", min_brlen, "). Returning original tree, pruningwise reordered.")
		cat("\n")
		cat(txt)
		cat("\n")
		return(phy2)
		} else {
		
		txt = paste0(num_branches_below_min, " branches found < min_brlen (", min_brlen, "). Running downpass to edit branchlengths by pushing nodes down so that minimum branchlength=", min_brlen, ". Adjusting i:Leftnode,j:Rightnode,anc:Ancnode;...")
		cat("\n")
		cat(txt)
		cat("\n\n")
		} 
		
		# END if (observed_min_orig > min_brlen)



	# DEFINE DOWNPASS THROUGH THE BRANCHES	
	tipnums <- 1:numtips
	i = 1
	edges_to_visit = seq(from=1, by=2, length.out=num_internal_nodes)

	#######################################################
	#######################################################
	# THIS IS A DOWNPASS FROM THE TIPS TO THE ROOT
	#######################################################
	#######################################################
	direct_ancestor_count = 0
	for (i in edges_to_visit)
		{
		# First edge visited is i
		#print(i)
		
		# Its sister is j 
		j <- i + 1
		#print(j)

		# Get the node numbers at the tips of these two edges		
		left_desc_nodenum <- phy2$edge[i, 2]
		right_desc_nodenum <- phy2$edge[j, 2]
		
		# And for the ancestor edge (i or j shouldn't matter, should produce the same result!!!)
		ancnode <- phy2$edge[i, 1]
		anc_edge_TF = phy2$edge[,2] == ancnode
		anc_edgenum = (1:nrow(phy2$edge))[anc_edge_TF]
		if (length(anc_edgenum) == 0)
			{
			anc_edgenum = NA
			}

		# Get the edge length (left, right, anc)
		Llength = phy2$edge.length[i]
		Rlength = phy2$edge.length[j]
		Alength = phy2$edge.length[anc_edgenum]
		
		txt = paste("ancnode:", ancnode, " left:", left_desc_nodenum, " right:", right_desc_nodenum, sep="")
		#print(txt)
		
		# If 0-length terminal branches are allowed (e.g. a 
		# sampled-ancestor tree), and
		# If a node is a terminal node,
		# and if it is of length 0 and the other branch > 0,
		# then don't change this branch length
		edit_brlength = TRUE

		if (leave_BL0_terminals == TRUE)
			{
			if (i <= numtips)
				{
				TF1  = (Llength >= 0) && (Llength <= direct_ancestor_brlen)
				TF2  = (Rlength >= 0) && (Rlength <= direct_ancestor_brlen)
				#cat(i, Llength, Rlength, sep="\t")
				#cat("\n")
				if ( TF1 && (Rlength > 0) )
					{
					edit_brlength = FALSE
					direct_ancestor_count = direct_ancestor_count + 1
					}
				} # END if (i <= numtips)
			if (j <= numtips)
				{
				if ( TF2 && (Llength > 0) )
					{
					edit_brlength = FALSE
					direct_ancestor_count = direct_ancestor_count + 1
					}
				} # END if (i <= numtips)
			} # END if (leave_BL0_terminals == TRUE)
		
		
		# If either the left or right branch is < min_brlen, change both
		if (edit_brlength == TRUE) {
		if ( (Llength < min_brlen) || (Rlength < min_brlen) )
			{
			cat(i, ":", left_desc_nodenum, ",", j, ":", right_desc_nodenum, ",anc:", ancnode, "; ", sep="")
			
			amount_to_add_L = min_brlen - Llength
			amount_to_add_R = min_brlen - Rlength
			amount_to_add = max(c(amount_to_add_L, amount_to_add_R))
			
			# Change both branchlengths by the same amount
			phy2$edge.length[i] = phy2$edge.length[i] + amount_to_add
			phy2$edge.length[j] = phy2$edge.length[j] + amount_to_add
			
			# Subtract that amount from the branch below
			# (may result in negative, which will be fixed further in the downpass)
			if (is.na(anc_edgenum) == FALSE)
				{
				phy2$edge.length[anc_edgenum] = phy2$edge.length[anc_edgenum] - amount_to_add
				}
			}
			} # END if (edit_brlength == TRUE)
		
		} # end downpass
	#######################################################
	# END DOWNPASS FROM THE TIPS TO THE ROOT
	#######################################################
	rootnode = ancnode
	phy3 = phy2

	cat("\n")
	cat("\n")
	num_branches_below_min2 = num_branches_below_min - direct_ancestor_count
	
	txt = paste0("Originally, ", num_branches_below_min2, " (non-ancestor) branches were found < min_brlen (", min_brlen, "). Running downpass to edit branchlengths by pushing nodes down so that minimum branchlength=", min_brlen, ". Adjusting i:Leftnode,j:Rightnode,anc:Ancnode;...")
	cat("\n")
	cat(txt)
		
	observed_min_brlen_new = min(phy3$edge.length)
	cat("\n\n")
	cat("New observed_min_brlen_new=", observed_min_brlen_new, sep="")
	cat("\n")
	return(phy3)
	} # END impose_min_brlen <- function(phy, min_brlen=0.1)








add_in_clades <- function(tmp_tree, names_of_groups, list_of_subtrees_w_good_branchlengths)
	{

	# Do any of the tip labels correspond to groups?
	TF = names_of_groups %in% tmp_tree$tip.label 
	if (sum(TF) == 0)
		{
		# None found, you're done
		final_tree = tmp_tree
		return(tmp_tree)
		} # END if (sum(TF) == 0)

	names_tmp = names_of_groups[TF]
	if (sum(TF) >= 1)
		{
		if (length(names_tmp) > length(unique(names_tmp)))
			{
			txt2 = paste0(names_tmp, sep=",")
			txt = paste0("ERROR in add_in_clade(): the clade names '", txt, "' match the tip labels in the growing master_tree ", sum(TF), " times.  The algorithm isn't working properly, or you have a tip name that is also a clade name (bad idea!).")
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			} # END if (length(names_tmp) > length(unique(names_tmp)))
	
		# Otherwise, go through the clades and add them
		nums = (1:length(names_of_groups))[TF]
		for (i in 1:length(nums))
			{
			tree_to_add_in = list_of_subtrees_w_good_branchlengths[[nums[i]]]
			TF2 = tmp_tree$tip.label == names_of_groups[nums[i]]
			tipnum_to_add_at = (1:length(tmp_tree$tip.label))[TF2]
			tmp_tree2 = bind.tree(x=tmp_tree, y=tree_to_add_in, where=tipnum_to_add_at, position=0)
			tmp_tree3 = read.tree(file="", text=write.tree(phy=tmp_tree2, file=""))
			tmp_tree4 = drop.tip(phy=tmp_tree3, tip=names_of_groups[nums[i]])
			tmp_tree5 = read.tree(file="", text=write.tree(phy=tmp_tree4, file=""))
			} # END for (i in 1:length(nums))
		} # END if (sum(TF) >= 1)

	# Check again; if any cladenames are found, add those in
	# (this will proceed recursively)
	TF = names_of_groups %in% tmp_tree5$tip.label 
	if (sum(TF) == 0)
		{
		# None found, you're done
		final_tree = tmp_tree5
		return(final_tree)
		} else {
		# Dynamic programming (recursive function)
		final_tree = add_in_clades(tmp_tree=tmp_tree5, names_of_groups=names_of_groups, list_of_subtrees_w_good_branchlengths=list_of_subtrees_w_good_branchlengths)
		return(final_tree)		
		} # END if (sum(TF) == 0)
	} # END add_in_clades <- function()



#' @OTUs_df=NULL A data.frame with, at least, columns labeled 
#' "OTUs" (tip names), 
#' "use" (which OTUs will be used, "yes" or "no"; blanks etc. mean "yes"), and 
#' "tipdate" (the desired tipdates; use 0, if all the tips are contemporaneous)
#' If OTUs_df is NULL, it can be read from the Excel file xlsfn; the program 
#' will skip the first 14 lines of the worksheet; the 15th row will have the 
#' column headers.  This is the exact command used by construct_starting_tree():
#' readWorksheetFromFile(xlsfn, sheet="OTUs", startRow=15)
#'
#' @nodes_df=NULL A data.frame with, at least, columns labeled 
#' "Taxon" (clade names, or non-clade groups of OTUs), 
#' "use" (which OTUs will be used, "yes" or "no"; blanks etc. mean "yes") 
#' "mono" Is the group monophyletic? "yes" or "no". Only "yes" will be used
#' for building the starting tree.
#' "make_age_prior" Should be "yes" or "no". Only "yes" will be used
#' in the starting tree
#' "distribution" What is the prior distribution on the node age? Can be
#' one of: normal, uniform, lognormal, exponential
#' meanInRealSpace For the lognormal distribution, is the mean parameter
#' (param1) the mean in real space ("yes"), or not ("no"). The function 
#' will convert between these.
#' offset: Offset to the mean (for lognormal or exponential)
#' param1: the mean, lower bound, mean, and mean, respectively
#' param2: the standard deviation, upper bound, standard deviation, 
#' and ignored, respectively
#' If nodes_df is NULL, it can be read from the Excel file xlsfn; the program 
#' will skip the first 14 lines of the worksheet; the 15th row will have the 
#' column headers.  This is the exact command used by construct_starting_tree():
#' readWorksheetFromFile(xlsfn, sheet="nodes", startRow=15)
#'
#' @taxa_df A data.frame with, at least, columns listing the OTUs that 
#' go into each clade. The column headers are the names of the clades.
#' construct_starting_tree() assumes that one column will have ALL of the 
#' OTUs, and will be named "total_group_LCA". This clade represents the 
#' Last Common Ancestor (LCA) of the total group.
#' 
#' All clade names in taxa_df should match a clade name in nodes_df$Taxon.
#' Those that don't will be cut.
#' 
#' As clades will have OTU lists of different lengths, the extra spaces in the 
#' table should be blanks ("") or NA.
#' OTUs missing from OTUs_df, or marked $use="no", will be cut.
#' If taxa_df is NULL, it can be read from the Excel file xlsfn; the program 
#' will skip the first 14 lines of the worksheet; the 15th row will have the 
#' column headers.  This is the exact command used by construct_starting_tree():
#' readWorksheetFromFile(xlsfn, sheet="taxa", startRow=15)
#' 
#' @xlsfn If the above inputs are NULL, then OTUs_df, nodes_df, and taxa_df
#' can be read from an Excel file formated according to the BEASTmasteR
#' example settings file.
#' 
#' @min_brlen During construction of the starting tree, all internal 
#' (non-tip) branchlengths will be min_brlen=0.001 (default). Changing 
#' this will make these larger, but might cause difficulties if you have 
#' node-date constraints that are close together
#' @outfn Write the tree to this file, if not "" (default is "", i.e. blank)
construct_starting_tree <- function(OTUs_df=NULL, taxa_df=NULL, nodes_df=NULL, xlsfn=NULL, min_brlen=0.0001, outfn="")
	{
	defaults='
	source("/drives/GDrive/__github/BEASTmasteR/R/construct_starting_tree_v1.R")
	construct_starting_tree(xlsfn=xlsfn, outfn="starting_tree.newick")
	
	trstr = construct_starting_tree(OTUs_df, taxa_df, nodes_df)
	
	OTUs_df=NULL; taxa_df=NULL; nodes_df=NULL
	'

	

	# Gather the inputs

	# OTUs
	if (is.null(OTUs_df) && !is.null(xlsfn))
		{
		OTUs_df = readWorksheetFromFile(xlsfn, sheet="OTUs", startRow=15)
		} else {
		if (is.null(OTUs_df) == TRUE)
			{
			stoptxt = "STOP ERROR in construct_starting_tree(): Your inputs to this function have OTUs_df=NULL and xlsfn=NULL. At least one of these must be non-null."
	
			cat("\n\n")
			cat(stoptxt)
			cat("\n\n")
			stop(stoptxt)
			}
		} # END if (is.null(OTUs_df) && !is.null(xlsfn))

	prune_OTUs = TRUE
	if (prune_OTUs == TRUE)
		{
		# Remove OTUs with "no" in "use" column of "OTUs" worksheet
		OTUs_df$use[isblank_TF(OTUs_df$use)] = "yes"
		keepTF = OTUs_df$use != "no"
		OTUs_df = OTUs_df[keepTF,]
		OTUs = trim(OTUs_df$OTUs)
		OTUs
		}


	# Taxa groups (list of OTUs in clades)
	if (is.null(nodes_df) && !is.null(xlsfn))
		{
		nodes_df = readWorksheetFromFile(xlsfn, sheet="nodes", startRow=15)
		} else {
		if (is.null(nodes_df) == TRUE)
			{
			stoptxt = "STOP ERROR in construct_starting_tree(): Your inputs to this function have nodes_df=NULL and xlsfn=NULL. At least one of these must be non-null."
	
			cat("\n\n")
			cat(stoptxt)
			cat("\n\n")
			stop(stoptxt)
			}
		} # END if (is.null(OTUs_df) && !is.null(xlsfn))

	# Filter taxa_df by "use"
	# AND by monophyletic (we will ignore non-monophyletic date constraints for now)
	nodes_df$use[isblank_TF(nodes_df$use)] = "yes"
	keepTF1 = nodes_df$use != "no"
	keepTF2 = nodes_df$use == "yes"
	keepTF = (keepTF1 + keepTF2) == 2
	nodes_df = nodes_df[keepTF,]

	# Error check:
	if (any(isblank_TF(nodes_df$mono)))
		{
		stoptxt = "STOP ERROR in construct_starting_tree(): The column nodes_df$mono has blank cells. All cells must be either 'yes' or 'no', indicating whether or not the taxa form a monophyletic group. Printing nodes_df$mono so you can inspect it."
	
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		cat("nodes_df$mono:\n")
		print(nodes_df$mono)
	
		stop(stoptxt)
		}




	# Taxa groups (list of OTUs in clades)
	if (is.null(taxa_df) && !is.null(xlsfn))
		{
		taxa_df = readWorksheetFromFile(xlsfn, sheet="taxa", startRow=15)
		} else {
		if (is.null(nodes_df) == TRUE)
			{
			stoptxt = "STOP ERROR in construct_starting_tree(): Your inputs to this function have taxa_df=NULL and xlsfn=NULL. At least one of these must be non-null."
	
			cat("\n\n")
			cat(stoptxt)
			cat("\n\n")
			stop(stoptxt)
			}
		} # END if (is.null(OTUs_df) && !is.null(xlsfn))


	# NAs in tipdates?
	if (any(isblank_TF(OTUs_df$tipdate)) == TRUE)
		{
		stoptxt = "STOP ERROR in construct_starting_tree(): There are blanks or NAs in OTUs_df$tipdate. Check your OTUs worksheet. Printing 'OTUs_df$tipdate'..."

		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		print("OTUs_df$tipdate:")
		print(OTUs_df$tipdate)
		cat("\n\n")
		stop(stoptxt)
		}
	

	# Gather clade names, remove empty clades
	names_of_groups = names(taxa_df)
	
	# Remove blanks!
	keepTF = isblank_TF(names(taxa_df)) == FALSE
	
	# Don't change, if sum(keepTF == 1)
	if (sum(keepTF) > 1)
		{
		taxa_df = taxa_df[,keepTF]
		}
	
	# Bug fix for keepTF=1
	if (sum(keepTF) == 1)
		{
		taxa_vector = taxa_df[,keepTF]
		taxa_df2 = as.data.frame(matrix(taxa_vector, ncol=1), stringsAsFactors=FALSE)
		names(taxa_df2) = names(taxa_df)[keepTF]
		taxa_df = taxa_df2
		}
	names_of_groups = names(taxa_df)
	
	
	# Check that taxa_df 
	# Are all names in taxa_df found in nodes_df$Taxon?
	TF1 = (nodes_df$Taxon %in% names_of_groups) == FALSE
	if (sum(TF1) > 0)
		{
		txt = "STOP ERROR in construct_starting_tree(): not all of the group names (with 'yes') in nodes_df$Taxon match names in taxa_df."
		cat("\n\n")
		cat(txt)
		cat("\n\nnodes_df$Taxon:\n\n")
		print(nodes_df$Taxon)
		cat("\n\ntaxa_df, 'names_of_groups':\n\n")
		print(names_of_groups)
		stop(txt)
		}

	
	# Check names
	keepTF = rep(TRUE, times=length(names_of_groups))
	for (i in 1:length(names_of_groups))
		{
		items = trim(taxa_df[,names_of_groups[i]])
		items2 = remove_blanks_NAs_etc(items)

		# Keep if the column has two or more non-blank taxa
		if (length(items2) < 2)
			{
			keepTF[i] = FALSE
			} # END if (length(items2) < 1)
		} # END for (i in 1:length(names_of_groups))

	names_of_groups = names_of_groups[keepTF]


	# Error check for unique groups
	if (length(names_of_groups) != length(unique(names_of_groups)))
		{
		stoptxt = "STOP ERROR in construct_starting_tree(): The names of your clades have to all be unique, but they are not. Printing the names of the clades."
	
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
	
		print(names_of_groups)
		cat("\n")
		stop(stoptxt)
		} # END if (length(names_of_groups) != length(unique(names_of_groups)))



	# Keep a list of the taxa which end up EMPTY -- we will CUT 
	# these from the priors list (or else it causes an error)
	cat("\n\n")
	list_of_empty_taxa = NULL
	list_of_clades = list()
	keepTF = rep(TRUE, length(names_of_groups))
	for (i in 1:length(names_of_groups))
		{
		groupname = names_of_groups[i]
	
		# Check if the groupname is used in the subset nodes_df
		if ((groupname %in% nodes_df$Taxon) == FALSE)
			{
			# Skip it
			keepTF[i] = FALSE
			next()
			}
	
		items = taxa_df[,groupname]
		items2 = remove_blanks_NAs_etc(items)

		# Filter out names not found in OTUs; throw warning for each
		if (is.null(OTUs) == FALSE)
			{
			TF = items2 %in% OTUs
		
			if (sum(TF) == 0)
				{
				txt = paste0("WARNING in make_taxa_groups(): None of the taxa names in '", groupname, "' are found in overall list in 'OTUs'. The group '", groupname, "' is being left out of the XML.")
				cat(txt)
				cat("\n")
				warning(txt)
			
				list_of_empty_taxa = c(list_of_empty_taxa, groupname)
				
				# Skip it
				keepTF[i] = FALSE
				
				# Skip to next
				next()
				} # END if (sum(TF) == 0)

			if (sum(TF) == 1)
				{
				txt = paste0("WARNING in make_taxa_groups(): In taxa '", groupname, "', there is only one species, after filtering according to the overall list in 'OTUs'. The group '", groupname, "' is being left out of the XML.")
				cat(txt)
				cat("\n")
				warning(txt)
			
				list_of_empty_taxa = c(list_of_empty_taxa, groupname)

				# Skip it
				keepTF[i] = FALSE
			
				# Skip to next
				next()
				} # END if (sum(TF) == 0)

		
			items_being_removed = items2[TF == FALSE]
		
			if (length(items_being_removed) > 0)
				{
				for (j in 1:length(items_being_removed))
					{
					txt = paste0("'", items_being_removed[j], "' is missing from overall 'OTUs' so is being cut from '", groupname, "'.")
					cat(txt)
					cat("\n")
					warning(txt)
					} # END for (j in 1:length(items_being_removed))
				} # END if (length(items_being_removed) > 0)
			# Subset the OTUs
			items2 = items2[TF]
			} # END if (is.null(OTUs) == FALSE)
	
		list_of_clades = c(list_of_clades, list(items2))
		} # END for (i in 1:length(names_of_groups))
	list_of_clades
	names_of_groups = names_of_groups[keepTF]

	# Order the clades from smallest to largest
	lengths = sapply(X=list_of_clades, FUN=length)

	#######################################################
	# If any of the clades are identical to total_group_LCA, cut them.
	#######################################################
	total_group_LCA_TF = (names_of_groups == "total_group_LCA")
	total_group_LCA_num = (1:length(total_group_LCA_TF))[total_group_LCA_TF]
	length_total_group_LCA = lengths[total_group_LCA_num]
	length_equals_max_TF = lengths == length_total_group_LCA
	length_equals_max_TF[total_group_LCA_num] = FALSE
	
	# Cut any groups that match total_group_LCA length
	length_equals_max_nums = (1:length(length_equals_max_TF))[length_equals_max_TF==FALSE]
	names_of_groups = names_of_groups[length_equals_max_nums]
	list_of_clades = list_of_clades[length_equals_max_nums]
	
	# Save these as the final (original) list of groups
	names_of_groups_orig = names_of_groups
	
	# Order the clades from smallest to largest
	lengths = sapply(X=list_of_clades, FUN=length)

	
	# Get the master_tree_name (the biggest clade, ie the total group)
	longest_TF = lengths == max(lengths)
	master_tree_name = names_of_groups_orig[longest_TF]
	
	clade_order_by_size = order(lengths)
	list_of_clades = list_of_clades[clade_order_by_size]
	names_of_groups = names_of_groups[clade_order_by_size]
	names(list_of_clades) = names_of_groups

	#######################################################
	# Now, if any of the clades are identical to each other,
	# cut them also
	#######################################################
	clade_txts = NULL
	for (i in 1:length(list_of_clades))
		{
		tmptxt = paste(sort(list_of_clades[[i]]), sep="", collapse=",")
		clade_txts = c(clade_txts, tmptxt)
		}
	# The last one is total_group_LCA
	# We need to check if any are identical to any later ones
	drop_TF = rep(FALSE, length(clade_txts))
	for (i in 1:(length(clade_txts)-1))
		{
		drop_TFs = unname(sapply(X=clade_txts[(i+1):length(clade_txts)], FUN=identical, y=clade_txts[i]))
		drop_TF[i] = sum(drop_TFs) > 0
		}
	keep_TFnums = (1:length(drop_TF))[drop_TF==FALSE]
	list_of_clades = list_of_clades[keep_TFnums]
	names_of_groups = names_of_groups[keep_TFnums]
	
	
	# List of which subclades go in which
	list_of_subclades_of_each_clade = list()
	if (length(names_of_groups) > 0)
		{
		for (i in 1:length(names_of_groups))
			{
			cmdtxt = paste0("list_of_subclades_of_each_clade$", names_of_groups[i], " = '", names_of_groups[i], "'")
			eval(parse(text=cmdtxt))
			}
		}

	# For each clade, check that it doesn't contradict larger clades
	if (length(list_of_clades) > 1)
		{
		for (i in 1:(length(list_of_clades)-1))
			{
			for (j in (i+1):length(list_of_clades))
				{
				smaller_clade_OTUs_in_larger_TF = list_of_clades[[i]] %in% list_of_clades[[j]]
		
				# Record the match, if it occurs
				if (sum(smaller_clade_OTUs_in_larger_TF) == length(smaller_clade_OTUs_in_larger_TF))
					{
					list_of_subclades_of_each_clade[[j]] = c(list_of_subclades_of_each_clade[[j]], names_of_groups[i])
					}
		
				if (sum(smaller_clade_OTUs_in_larger_TF) != length(smaller_clade_OTUs_in_larger_TF))
					{
					if (sum(smaller_clade_OTUs_in_larger_TF) != 0)
						{
						stoptxt = "STOP ERROR in construct_starting_tree(). You have a contradiction between clades that you have specified should be monophyletic -- but the first clade only partially fits inside the second clade. Printing the offending clades for your inspection."

						cat("\n\n")
						cat(stoptxt)
						cat("\n\n")
						cat("Clade #1: ", names_of_groups[i], ", printing list of OTUs:")
						cat("\n")
						print(list_of_clades[[i]])
						cat("\n")
						cat("Clade #2: ", names_of_groups[j], ", printing list of OTUs:")
						cat("\n")
						print(list_of_clades[[j]])
						cat("\n")
						stop(stoptxt)
						} # END if (sum(smaller_clade_OTUs_in_larger_TF) != 0)
					} # END if (sum(smaller_clade_OTUs_in_larger_TF) != length(smaller_clade_OTUs_in_larger_TF))
				} # END for (j in (i+1):length(list_of_clades))
			} # END for (i in 1:(length(list_of_clades)-1))
		list_of_subclades_of_each_clade

		# List of tips in each subclade, after collapsing 1 level down
		list_of_clades_collapsed = list_of_clades

		# Collapse tipnames that go in a subclade
		for (j in (length(list_of_clades)):2)
			{
			for (i in (j-1):1)
				{
				smaller_clade_OTUs_in_larger_TF = list_of_clades_collapsed[[i]] %in% list_of_clades_collapsed[[j]]
		
				# Record the match, if it occurs
				if (sum(smaller_clade_OTUs_in_larger_TF) == length(smaller_clade_OTUs_in_larger_TF))
					{
					bigger_OTUs_list = list_of_clades_collapsed[[j]]
					smaller_list = list_of_clades_collapsed[[i]]
					OTUs_to_collapse_TF = bigger_OTUs_list %in% smaller_list
					bigger_OTUs_list = bigger_OTUs_list[OTUs_to_collapse_TF == FALSE]
					bigger_OTUs_list = c(bigger_OTUs_list, names_of_groups[i])
					list_of_clades_collapsed[[j]] = bigger_OTUs_list
					} # END if (sum(smaller_clade_OTUs_in_larger_TF) != length(smaller_clade_OTUs_in_larger_TF))
				} # END for (j in (i+1):length(list_of_clades))
			} # END for (i in 1:(length(list_of_clades)-1))
		} else {
		list_of_clades_collapsed = list_of_clades
		} # END if (length(list_of_clades) > 1)

	list_of_clades_collapsed

	# Make subtrees, from smallest to biggest
	list_of_subtrees = list()
	if (length(list_of_clades_collapsed) > 0)
		{
		for (j in 1:length(list_of_clades_collapsed))
			{
			tmpOTUs = list_of_clades_collapsed[[j]]
			tmptr = rtree(n=length(tmpOTUs), rooted=TRUE, tip.label=tmpOTUs, br=runif, min=min_brlen, max=min_brlen)
			list_of_subtrees[[j]] = tmptr
			}
		}
	list_of_subtrees




	# Put dates on these subtrees (where available)
	nodes_df$offset[isblank_TF(nodes_df$offset)] = 0

	# Order nodes_df to match order of names_of_groups
	matches = match(nodes_df$Taxon, table=names_of_groups)
	keep_TF = !is.na(matches)
	nodes_df = nodes_df[keep_TF,]
	matches = match(nodes_df$Taxon, table=names_of_groups)
	nodes_df = nodes_df[matches,]

	desired_node_ages = NULL
	names_of_groups_w_node_ages = NULL
	if (length(names_of_groups) > 0)
		{
		for (i in 1:length(names_of_groups))
			{
			# Match OTUs to tree tiplabel order,
			# get the dates of all of the included tips
			OTUs_in_subtree_TF = OTUs %in% list_of_clades[[i]]
			subtree_OTUs = OTUs[OTUs_in_subtree_TF]
			subtree_tipdates = OTUs_df$tipdate[OTUs_in_subtree_TF]
			if (length(subtree_tipdates) == 0)
				{
				subtree_tipdates = 0
				}
			oldest_tip_age = max(subtree_tipdates)
		
			# If a desired age has been specified...
			if (nodes_df$make_age_prior[i] == "yes")
				{
				names_of_groups_w_node_ages = c(names_of_groups_w_node_ages, desired_node_ages[i])
			
				if (tolower(nodes_df$distribution[i]) == "normal")
					{
					desired_node_age = nodes_df$param1[i] + nodes_df$offset[i]
					}
				if (tolower(nodes_df$distribution[i]) == "uniform")
					{
					desired_node_age = (nodes_df$param1[i] + nodes_df$param2[i])/2 + nodes_df$offset[i]
					}
				if (tolower(nodes_df$distribution[i]) == "lognormal")
					{
					if (nodes_df$meanInRealSpace[i] == "yes")
						{
						desired_node_age = nodes_df$param1[i] + nodes_df$offset[i]
						}
					if (nodes_df$meanInRealSpace[i] == "no")
						{
						meanval = exp(nodes_df$param1[i])
						desired_node_age = meanval + nodes_df$offset[i]
						}
					}
				if (tolower(nodes_df$distribution[i]) == "exponential")
					{
					if (nodes_df$meanInRealSpace[i] == "yes")
						{
						desired_node_age = nodes_df$param1[i] + nodes_df$offset[i]
						}
					if (nodes_df$meanInRealSpace[i] == "no")
						{
						meanval = 1/(nodes_df$param1[i])
						desired_node_age = meanval + nodes_df$offset[i]
						}
					}
		
	
				# If the desired node age is less than oldest_tip_age, re-set to tip age + 1%
				if (desired_node_age <= oldest_tip_age)
					{
					desired_node_age = oldest_tip_age * 1.01
					}
				} # END if (nodes_df$make_age_prior[i] == "yes")
		
			# If no desired age specified, just make it a little older than the tips
			if (nodes_df$make_age_prior[i] != "yes")
				{
				desired_node_age = oldest_tip_age * 1.01
				} # END if (nodes_df$make_age_prior[i] == "yes")
		
			# Save the desired node age
			desired_node_ages = c(desired_node_ages, desired_node_age)
			} # END for (i in 1:length(names_of_groups))
		} else {
		# nada
		}


	# Apply names to desired_node_ages
	if (is.null(desired_node_ages) == FALSE)
		{
		names(desired_node_ages) = names_of_groups
		desired_node_ages
		}


	# OK, now you have:
	desired_node_ages # the goal for each node age
	list_of_subtrees  # the random subtrees, with ultrashort branch lengths
	list_of_clades_collapsed	# which OTUs are in each subclade, with OTUs collapse into subclades
	list_of_subclades_of_each_clade	# which subclades are in each clade
	names_of_groups	# name of each monophyletic clade


	# For each subtree:
	# Starting from the root, make subtrees, adding in the sub-subtrees as we go
	j = 2
	list_of_subtrees_w_good_branchlengths = list_of_subtrees
	if (length(list_of_subtrees) > 0)
		{
		for (j in 1:length(list_of_subtrees))
			{
			subtree = list_of_subtrees[[j]]

			# Get the desired tip date for each OTU or clade root
			clade_or_tipnames = list_of_clades_collapsed[[j]]
			tip_ages = rep(0, length(clade_or_tipnames))
			names(tip_ages) = clade_or_tipnames

			for (i in 1:length(clade_or_tipnames))
				{
				tmpname = clade_or_tipnames[i]
				TF = (names_of_groups %in% tmpname) 
				if (sum(TF) == 1)
					{
					# Then it's a node, get the desired node age
					tip_ages[i] = desired_node_ages[TF]
					} else {
					# Find it in the OTUs list
					TF = (OTUs_df$OTUs %in% tmpname)
					if (sum(TF) == 1)
						{
						tip_ages[i] = OTUs_df$tipdate[TF]
						} else {
						txt = paste0("STOP ERROR in construct_starting_tree(): the tip or clade '", tmpname, "' had no date specified that we could find.")
						cat("\n\n")
						cat(txt)
						cat("\n\n")
						stop(txt)
						} # END if (sum(TF) == 1)
					} # END if (sum(TF) == 1)
				} # END for (i in 1:length(clade_or_tipnames))

			tip_ages

			# Get the desired age of the root node
			name_of_group = names_of_groups[j]
			TF = (names_of_groups %in% name_of_group) 
			if (sum(TF) == 1)
				{
				node_age = desired_node_ages[TF]
				} else {
				# If not found, the node age will be 10% of the range of the tip ages
				range_val = max(tip_ages) - min(tip_ages)
				node_age = max(tip_ages) + 0.1 * max(tip_ages)
				}

			# Now, adjust these ages relative to the highest tip
			youngest_tip_age = min(tip_ages)
			node_age = node_age - youngest_tip_age
			tip_ages = tip_ages - youngest_tip_age

			# Get the subtree table
			subtr_table = prt(subtree, printflag=FALSE)
			subtr_table

			# How much to add to tips to get the subtree root the correct relative date?
			root_nodenum = length(subtree$tip.label)
			subtr_current_root_age = subtr_table$time_bp[root_nodenum]
			add_to_tips = node_age - tip_ages
			add_to_tips

			# Add to the tip branchlengths
			for (i in 1:length(add_to_tips))
				{
				tip_to_add_to = names(add_to_tips)[i]
				TFrow = subtr_table$label == tip_to_add_to
				edgenum = subtr_table$parent_br[TFrow]
				subtree$edge.length[edgenum] = subtree$edge.length[edgenum] + add_to_tips[i]
				}
	
			list_of_subtrees_w_good_branchlengths[[j]] = subtree
			} # END for (j in 1:length(list_of_subtrees))
		} else {
		
		} # END if (length(list_of_subtrees) > 0)

	# Let's assume the beginning tree, that we will build from,
	# is called total_group_LCA
	#master_tree_name = "total_group_LCA"
	TF = names_of_groups == master_tree_name
	num = (1:length(names_of_groups))[TF]
	master_tree = list_of_subtrees_w_good_branchlengths[[num]]
	master_tree



	# Add clades
	tmp_tree = master_tree


	final_tree = add_in_clades(tmp_tree=tmp_tree, names_of_groups=names_of_groups, list_of_subtrees_w_good_branchlengths=list_of_subtrees_w_good_branchlengths)


	tbl = prt(final_tree, printflag=FALSE)
	tbl[,c(9,10,11,12)][1:length(final_tree$tip.label),]


	# Average the living tips to have an age of 0
	# (fixes minor rounding errors)
	#final_tree = average_tr_tips(tr=final_tree, fossils_older_than=0.0001)
	#tbl = prt(final_tree, printflag=FALSE)
	#tbl[,c(9,10,11,12)][1:length(final_tree$tip.label),]

	# Run several times for complex trees
	final_tree2 = level_tree_tips(tr=final_tree, method="highest", printflag=FALSE, fossils_older_than=0.0001)

	tbl = prt(final_tree2, printflag=FALSE)
	tbl[,c(9,10,11,12)][1:length(final_tree2$tip.label),]

	# plot(final_tree)
	# axisPhylo()

	# Space the branchlengths out a bit more
	if (any(final_tree$edge.length < 0.0000000000001))
		{
		final_tree = impose_min_brlen(phy=final_tree2, min_brlen=min_brlen, leave_BL0_terminals=FALSE)
		} else {
		final_tree = final_tree2
		}


	# write.tree(final_tree, file="starting_tree.newick")


	# Run some error checks
	ntips = length(final_tree$tip.label)

	# Error check: assembled tree has the same number of tips as OTUs
	problem = FALSE
	if (ntips != length(OTUs))
		{
		problem = 1
		}
	OTUs_in_final_tree_TF = OTUs %in% final_tree$tip.label
	if (sum(OTUs_in_final_tree_TF) != length(OTUs_in_final_tree_TF))
		{
		problem = 2
		}
	final_tree_in_OTUs_TF =  final_tree$tip.label %in% OTUs
	if (sum(final_tree_in_OTUs_TF) != length(final_tree_in_OTUs_TF))
		{
		problem = 3
		}
	if (problem != FALSE)
		{
		stoptxt = "STOP ERROR in construct_starting_tree(): Something went wrong in assembling the topology of the starting tree. The OTUs in 'final_tree$tip.label' and in 'OTUs' do not match. Printing them for inspection."
	
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		cat("(Problem #", problem, ")", sep="")
		cat("\n\n")
	
		cat("final_tree$tip.label:")
		cat("\n")
		cat(sort(final_tree$tip.label), sep=", ")
		cat("\n")

		cat("OTUs:")
		cat("\n")
		cat(sort(OTUs), sep=", ")
		cat("\n")
		
		cat("\nThis OTU might not be in the final tree:\n")
		print(OTUs[OTUs_in_final_tree_TF==FALSE])

		stop(stoptxt)
		}
	
	
	#######################################################
	# Check constructed tree for negative branch lengths
	#######################################################
	if (any(final_tree$edge.length < 0.0000000000001))
		{
		stoptxt = "STOP ERROR in construct_starting_tree(): Some edges in the constructed tree 'final_tree' have negative branchlengths."
		
		cat("\n\nPlotting the tree and printing Newick string...\n\n")
		cat(write.tree(final_tree, file=""))
		plot(final_tree)
		
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		stop(stoptxt)
		}
	
	# Turn off scientific notation
	#options(scipen=999)
	
	cat("\n\n") 
	cat("construct_starting_tree() is returning this starting tree. Newick string is pasted below:\n\n")
	final_tree_to_print = final_tree
	#new_brlengths = format(as.numeric(final_tree$edge.length), scientific=FALSE)
	#new_brlengths = gsub(pattern=" ", replacement="", x=new_brlengths)
	#final_tree_to_print$edge.length = new_brlengths
	trstr = write_tree_noSci(final_tree, file="")
	
	
	cat(trstr)
	cat("\n\n")
		
	if (outfn != "")
		{
		write_tree_noSci(phy=final_tree, file=outfn)
		}

	# Turn scientific notation back on
	#options(scipen=999)
	
	return(final_tree)
	} # END construct_starting_tree <- function(OTUs_df=NULL, taxa_df=NULL, nodes_df=NULL, xlsfn=NULL)
