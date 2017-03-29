# Extract the string (strings?) matching the regexp
extract_regexp <- function(tmpstr, tmpre_string)
	{
	matches = gregexpr(tmpre_string, tmpstr)[[1]]
	matches_end = matches-1+attr(matches,"match.length")
	x = mapply(substr, tmpstr, matches, matches_end)
	return(x)
	}



headtail <- function(fn, startline=1, endline=3, outputtxt="| more", outfn=NULL, printall=FALSE)
	{
	
	
	
	# 'head' gets your ending line
	headstart = endline
	
	# 'tail' gets endline-startline+1
	tailstart = endline - startline + 1
	
	if (is.null(outfn))
		{
		cmdstr = paste("head -", headstart, " ", fn, " | tail -", tailstart, " ", outputtxt, sep="")
		} else {
		cmdstr = paste("head -", headstart, " ", fn, " | tail -", tailstart, " > ", outfn, sep="")
		}
	
	if (printall != FALSE)
		{
		cat("\n")
		cat(paste("Excerpting lines ", startline, "-", endline, " of ", fn, ":\n", sep=""))
		cat(paste("Running: ", cmdstr, sep=""))
		cat("\n")
		}
	#system("pwd")
	
	
	system(cmdstr)
	
	return(outfn)
	}



# Get a few lines from a large file...
# http://compgroups.net/comp.unix.solaris/how-to-extract-a-few-lines-from-a-big-file/73315
# 
tailhead <- function(fn, startline=1, endline=3, outputtxt="| more", outfn=NULL, printall=FALSE)
	{
	
	
	
	# 'head' gets your ending line
	#headstart = endline
	
	# 'tail' gets endline-startline+1
	head_end = endline - startline + 1
	
	if (is.null(outfn))
		{
		cmdstr = paste("tail +", startline, " ", fn, " | head -", head_end, " ", outputtxt, sep="")
		} else {
		cmdstr = paste("tail +", startline, " ", fn, " | head -", head_end, " > ", outfn, sep="")
		}
	
	if (printall != FALSE)
		{
		cat("\n")
		cat(paste("Excerpting lines ", startline, "-", endline, " of ", fn, ":\n", sep=""))
		cat(paste("Running: ", cmdstr, sep=""))
		cat("\n")
		}
	#system("pwd")
	
	
	system(cmdstr)
	
	return(outfn)
	}


# Extract numbers / get numbers / get all numbers from a text string
getnums <- function(tmpstr)
	{
	# Example string
	# tmpstr = "The first number is: 32.  342 342.1   -3234e-10  3234e-1 Another one is: 32.1. Here's a number in scientific format, 0.3523e10, and another, 0.3523e-10, and a negative, -313.1"
	
	# To get gsubfn to work, you need the default package tcltk, but this needs this to work on Mac OS X:
	# http://rug.mnhn.fr/seewave/inst.html
	# 
	
	require("gsubfn")	# for strapply
	
# 	patternslist = NULL
# 	p=0
# 	patternslist[[(p=p+1)]] = "(\\d+)"				# positive integer
# 	patternslist[[(p=p+1)]] = "(-\\d+)"				# negative integer
# 	patternslist[[(p=p+1)]] = "(\\d+\\.\\d+)"		# positive float
# 	patternslist[[(p=p+1)]] = "(\\d+\\.\\d+e\\d+)"	# positive float, scientific w. positive power
# 	patternslist[[(p=p+1)]] = "(\\d+\\.\\d+e-\\d+)" # positive float, scientific w. negative power
# 	patternslist[[(p=p+1)]] = "(-\\d+\\.\\d+)"		# negative float
# 	patternslist[[(p=p+1)]] = "(-\\d+\\.\\d+e\\d+)"	# negative float, scientific w. positive power
# 	patternslist[[(p=p+1)]] = "(-\\d+\\.\\d+e-\\d+)"# negative float, scientific w. negative power
# 	
# 	patternslist[[(p=p+1)]] = "(\\d+e\\d+)"			# positive int, scientific w. positive power
# 	patternslist[[(p=p+1)]] = "(\\d+e-\\d+)" 		# positive int, scientific w. negative power
# 	patternslist[[(p=p+1)]] = "(-\\d+e\\d+)"		# negative int, scientific w. positive power
# 	patternslist[[(p=p+1)]] = "(-\\d+e-\\d+)"		# negative int, scientific w. negative power
# 	
# 	pattern = paste(patternslist, collapse="|", sep="")

	# set up the pattern
	# long pattern: https://stat.ethz.ch/pipermail/r-help/2013-June/355429.html
	# require("stringr")
	# pattern = "(\\d+)|(-\\d+)|(\\d+\\.\\d+)|(\\d+\\.\\d+e\\d+)|(\\d+\\.\\d+e-\\d+)|(-\\d+\\.\\d+)|(-\\d+\\.\\d+e\\d+)|(-\\d+\\.\\d+e-\\d+)|(\\d+e\\d+)|(\\d+e-\\d+)|(-\\d+e\\d+)|(-\\d+e-\\d+)"
	# Get the numbers
	# nums_from_tmpstr = as.numeric(str_extract_all(tmpstr,pattern)[[1]])
	
	# shorter pattern!
	# https://stat.ethz.ch/pipermail/r-help/2013-June/355435.html
	pattern = "[-+.e0-9]*\\d"
	
	# Get the numbers
	nums_from_tmpstr = strapply(tmpstr, pattern, as.numeric)[[1]]

	# Return them
	return(nums_from_tmpstr)
	}




# Extract a sample of trees from BEAST posterior sample, write each to a 
# nice simple labeled NEXUS file
subset_BEAST_trees_to_NEXUS <- function(fn, outfn, numlines=NULL, numtrees=3, burnin_fraction=0.5, header_in_first=500)
	{
	defaults='
	fn = nexfn
	numlines = 323970
	numtrees=1000
	burnin_fraction=0.5
	runthis = TRUE
	'

	# Scan the beginning of the file and figure out where trees start
	# Make sure nlines is big enough to get to the trees (i.e. it has to be at least ntaxa + some)
	# This program assumes the header ends within the first 500 (by default of header_in_first)
	X <- scan(file = fn, what = "", sep = "\n", quiet = TRUE, blank.lines.skip=FALSE, nlines=header_in_first)
	
	# Get line numbers for END; and ENDBLOCK;
	endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
	
	# Get line numbers with semicolons
	semico <- grep(";", X)
		
	# Line # for beginning of TREES block
	i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
	
	# End of everything before trees
	semico_after_begin_trees = semico[semico > i1][1]
	start_of_actual_trees = semico[semico > i1][2]


	
	# Get the overall file length
	if (is.null(numlines))
		{
		numlines = linecount(fn)
		} else {
		numlines = numlines
		}
	# 75091
	numlines
	
	# Set burnin
	burnin = round(burnin_fraction * (numlines - start_of_actual_trees - 1))
	
	# Number of trees to sample
	numtrees = numtrees
	
	#treedir = paste(wd, "sampled_newicks", "/", sep="")
	#treedir = "/Users/nickm/Desktop/__projects/_hominin_phylo_dating/_graphics/sampled_newicks/"
	


	# Extract the header
	headerstartline = 1
	headerendline = semico_after_begin_trees
	
	trees_startline = start_of_actual_trees
	trees_endline = numlines - 2
	trees_startline = trees_startline + burnin
	
	# Extract the header
	tmpfn = "temp_NEXUS_header.txt"
	header_txt_fn = headtail(fn=fn, startline=headerstartline, endline=headerendline, outfn=tmpfn, printall=TRUE)
	header_txt_fn
	moref(header_txt_fn)
	
	# Pick numtrees random trees WITHOUT replacement
	randlinenums = sort(sample(trees_startline:trees_endline, numtrees, replace=FALSE))
	
	# The differences between each randline (so you know how many to skip)
	# http://tolstoy.newcastle.edu.au/R/e2/help/07/02/9709.html
	skips = diff(c(1, randlinenums)) - 1 
	skips
	
	# Open a connection to the HUGE text file
	input <- file(fn, "r") 
	remaining_lines = numlines
	sel = length(randlinenums)
	# Vector to allocate my selected lines
	mysel <- vector('character', sel) 
	
	print(system.time())
	for (i in 1:sel)
		{
		tmpstr = paste("Seeking line ", randlinenums[i], "...", sep="")
		cat(tmpstr)
		mysel[i] <- scan(input, what="", sep="\n", skip=skips[i], n=1, quiet=TRUE)
		tmpstr = paste("extracted ", i, " of ", sel, ".\n", sep="")
		cat(tmpstr)
		}
	print(system.time())	
	
	# Now, cat the extracted lines to the header
	write_lines_good(dtf=mysel, outfn=tmpfn, sepval="\n", tmpappend=TRUE)
	
	# And add "END;"
	write("END;", file=tmpfn, append=TRUE, sep="\n")
	
	
	
	cmdstr = paste("mv ", tmpfn, " ", outfn, sep="")
	system(cmdstr)
	
	txt = paste("Subset trees output to: ", outfn, sep="")
	cat("\n", txt, "\n", sep="")
	
	return(outfn)
	}





# Extract a sample of trees from BEAST posterior sample, write each to a 
# nice simple labeled NEWICK file
sample_BEAST_trees_to_newick_files <- function(fn, treedir, numtrees=100, burnin_fraction=0.5, runthis=TRUE)
	{
	defaults='
	fn = morph_only_trees_fn
	treedir = "/Users/nickm/Desktop/__projects/_hominin_phylo_dating/_graphics/sampled_newicks/"
	numtrees=100
	burnin_fraction=0.5
	runthis = TRUE
	'

	# Scan the beginning of the file and figure out where trees start
	# Make sure nlines is big enough to get to the trees (i.e. it has to be at least ntaxa + some)
	X <- scan(file = fn, what = "", sep = "\n", quiet = TRUE, blank.lines.skip=FALSE, nlines=500)
	
	# Get line numbers for END; and ENDBLOCK;
	endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
	
	# Get line numbers with semicolons
	semico <- grep(";", X)
		
	# Line # for beginning of TREES block
	i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
	
	# End of everything before trees
	semico_after_begin_trees = semico[semico > i1][1]
	start_of_actual_trees = semico[semico > i1][2]


	
	# Get the overall file length
	numlines = linecount(fn)
	# 75091
	numlines
	
	# Set burnin
	burnin = round(burnin_fraction * (numlines - start_of_actual_trees - 1))
	
	# Number of trees to sample
	numtrees = numtrees
	
	#treedir = paste(wd, "sampled_newicks", "/", sep="")
	#treedir = "/Users/nickm/Desktop/__projects/_hominin_phylo_dating/_graphics/sampled_newicks/"
	


	# Extract the header
	headerstartline = 1
	headerendline = semico_after_begin_trees
	
	trees_startline = start_of_actual_trees
	trees_endline = numlines - 2
	trees_startline = trees_startline + burnin
	
	# Extract the header
	outfn = "temp_NEXUS_header.txt"
	header_txt_fn = headtail(fn=fn, startline=headerstartline, endline=headerendline, outfn=outfn, printall=TRUE)
	header_txt_fn
	moref(header_txt_fn)
	
	# Pick numtrees random trees WITHOUT replacement
	randlinenums = sample(trees_startline:trees_endline, numtrees, replace=FALSE)
	
	# Select the random trees
	#runthis = TRUE
	if (runthis)
		{
		for (i in 1:numtrees)
			{
			# Print
			cat("\nProcessing tree #", i, sep="")
			
			# Select a random line number
			randlinenum = randlinenums[i]
			
			# Extract a line from the treefile
			tmpoutfn = gsub("\\.nexus", "\\.nexus_subset", fn)
			tmpoutfn = extract_fn_from_path(tmpoutfn)
			tmpoutfn = headtail(fn, startline=randlinenum, endline=randlinenum, outfn=tmpoutfn)
			
			# cat the header and the tree together
			cmdstr = paste("cat ", header_txt_fn, " ", tmpoutfn, " > temp.nexus", sep="")
			system(cmdstr)
			
			# Append "END;" to end of file
			write("END;", file="temp.nexus", append=TRUE, sep="\n")
		
			# Read the temporary NEXUS file
			tmptr = read.nexus("temp.nexus")
			
			# make e.g. 0001
			tmpi = sprintf("%05.0f", i)
			
			# Write to NEWICK
			outnewick = paste(treedir, tmpoutfn, "_", tmpi, ".newick", sep="")
			cat("\n")
			cat(paste("Writing: ", outnewick, sep=""))
			write.tree(tmptr, file=outnewick)
			}
		
		# Plot the last tree you processed...
		# plot(tmptr2)
		}
	
	# Get the newick files
	newick_fns = slashslash(list.files(path=treedir, pattern=".newick", full.names=TRUE))
	newick_fns
	
	return(newick_fns)
	}






#######################################################
# extractBEASTstats_orig: Getting stats from a BEAST MCC summary tree NEXUS file
#######################################################
# extractBEASTstats_orig was copied/lightly modified from 
# phyloch::extractBEASTstats
# R package "phyloch", by Christoph Heibl <christoph.heibl@gmx.net>, 
# licensed under GPL (>=2), 
# available at: http://www.christophheibl.de/Rpackages.html


# Getting stats from a BEAST consensus tree
# From phyloch
'
file = trfn
file = confn
'
extractBEASTstats_orig <- function (file, digits=4, printflag=FALSE) 
	{
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    X <- X[grep("tree TREE1[[:space:]]+=", X)]
    X <- gsub("tree TREE1[[:space:]]+= \\[&R\\] ", "", X)
    tab <- unlist(strsplit(X, "\\["))[-1]
    tab <- gsub("&|;|\\]", "", tab)
    tab <- gsub(":.+$", "", tab)


    # This extracts ~19 objects delimited by [brackets
    # In the consensus tree output
    foo <- function(x)
    	{
        x <- unlist(strsplit(x, ","))
        x
    	}
    tab <- lapply(tab, foo)
    

	#for (i in seq(along = tab))
	# same as
	for (i in 1:length(tab))
    	{
    	# in each element in a tab(le) item, get the
    	# indices of anything containing "{"
        ind <- grep("[{]", tab[[i]])
        
        # Change the name text of those elements (remove e.g. {)
        names <- gsub("=.+$", "", tab[[i]][ind])
        tab[[i]][ind] <- gsub("[{]", "", tab[[i]][ind])
        
        # Replace = with _MIN=
        tab[[i]][ind] <- gsub("=", "_MIN=", tab[[i]][ind])
        
        # Remove the } bracket also
        tab[[i]][ind + 1] <- gsub("[}]", "", tab[[i]][ind + 1])
        
        # And call this item MAX=
        tab[[i]][ind + 1] <- paste(paste(names, "MAX=", sep = "_"), tab[[i]][ind + 1])
    	}


	# Empty data frame
    ttab <- data.frame()
    
    # Take everything before the = signs, and uniqify it
    stats <- unique(gsub("=.+$", "", unlist(tab)))
    
    # Remove "!rotate" (A FigTree code)
    stats = stats[stats != "!rotate"]
    stats
    
    # Go through all the tabs/nodes
    # Put them in a big table for each node
    for (i in seq(along = tab))
    	{
    	
    	# Go through all the summary statistics by name
        for (j in seq(along = stats))
        	{
        	# Get the index of the summary statistic you are looking for
            ind <- grep(paste("^", stats[j], "=", sep = ""), tab[[i]])
            
            # If the summary statistic is found at this node,
            if (length(ind) > 0)
            	{
            	# Remove the name information, and make the number numeric
                v <- as.numeric(gsub(paste(stats[j], "=", sep = ""), "", tab[[i]][ind]))
                
                if (printflag == TRUE)
                	{
	                cat(i, j, ind, v, "\n", sep="	")
	                }
                
                ttab[i, j] <- v
				}
			}
		}
		
	# Name the columns of the table
    colnames(ttab) <- stats
    
    # Figure out which are the tips (doesn't matter, really)
    tip <- which(is.na(ttab$posterior))
    ttab
    
    
# 	phy <- read.beast(file, digits = digits)
# 	int <- phy$Nnode
#     tips <- length(phy$tip.label)
#     node <- (tips + 1):(tips + int)
#     M <- cbind(node, ttab)
#  	   
#     return(M)
	}



#######################################################
# read.beast.table_original: Getting stats from a BEAST consensus tree
#######################################################
# read.beast.table_original was copied/lightly modified from 
# phyloch::read.beast.table
# R package "phyloch", by Christoph Heibl <christoph.heibl@gmx.net>, 
# licensed under GPL (>=2), 
# available at: http://www.christophheibl.de/Rpackages.html

# Getting stats from a BEAST consensus tree
# From phyloch
read.beast.table_original <- function (file, digits = 2) 
	{
    phy <- read.beast_original(file, digits = digits)
    int <- phy$Nnode
    stats <- phy[-(1:4)]
    M <- matrix(unlist(stats), nrow = int, byrow = FALSE)
    colnames(M) <- names(stats)
    tips <- length(phy$tip.label)
    node <- (tips + 1):(tips + int)
    M <- cbind(node, M)
    M
	} # END read.beast.table_original <- function (file, digits = 2) 


#######################################################
# extractBEASTstats3: Another mod of extractBEASTstats
#######################################################
# extractBEASTstats3 was copied/moderately modified from 
# phyloch::extractBEASTstats
# R package "phyloch", by Christoph Heibl <christoph.heibl@gmx.net>, 
# licensed under GPL (>=2), 
# available at: http://www.christophheibl.de/Rpackages.html
'
file=confn
'

extractBEASTstats3 <- function (file) 
{
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    phy <- read.nexus(file)
    X <- X[grep("tree TREE1[[:space:]]+=", X)]
    X <- gsub("^.*tree TREE1[[:space:]]+= \\[&R\\] ", "", X)
    X <- gsub("[!]rotate=true,*|[!]rotate=false,*", "", X)
    X <- gsub(";$", "", X)
    vals <- unlist(strsplit(X, "\\][[:alnum:]:.)]*\\[*&*"))
    foo <- function(x) {
        x <- gsub(",([[:lower:]])", "xxx\\1", x)
        x <- unlist(strsplit(x, "xxx"))
        names(x) <- sapply(x, function(x) {
            unlist(strsplit(x, split = "="))[1]
        })
        x <- sapply(x, function(x) {
            unlist(strsplit(x, split = "="))[2]
        })
        x <- gsub("[{}]", "", x)
        x <- strsplit(x, ",")
        lapply(x, as.numeric)
    }
    vals <- lapply(vals, foo)
    edges <- gsub("\\[(&[[:alnum:]_=%!.,{}-]+)\\]", "", X)
    edges <- unlist(strsplit(gsub("[()]*", "", edges), ":"))
    tips <- gsub("^[0-9]+.[0-9]+,*", "", edges)


	re1='(,)'								# Any Single Character 1
	#re2='(?:[a-z][a-z]*[0-9]+[a-z0-9]*)'	# Alphanum 1
	re2='(?:[a-zA-Z0-9]*)'					# Alphanum 1
    regstr = paste(re1, re2, sep="")
    
    #gsub(regstr, "", "asd,234g")


    # remove everything after comma
    # original
    # edges <- as.numeric(gsub("\\[,$|[[:alpha:]].+$\\]", "", edges))
    # modified
    edges <- as.numeric(gsub(pattern=regstr, "", edges))
    
    tab <- cbind(edges, nodes = c("root", head(tips, -1)), id = rep(NA, 
        length(edges)))
    tab[tab[, 2] == "", 2] <- "internal"
    internal <- grep("internal|root", tab[, 2])
    tab[-internal, 3] <- seq(along = phy$tip.label)
    intnodes <- match(tab[internal, 1], phy$edge.length[phy$edge[, 
        2] > 50], nomatch = 0) + 51
    tab[internal, 3] <- intnodes
    names(vals)[internal] <- paste("NODE", tab[internal, 3], 
        sep = ":")
    names(vals)[-internal] <- paste("TIP", tab[-internal, 3], 
        sep = ":")
    tips <- vals[-internal]
    nodes <- vals[internal]
}


#######################################################
# extractMrBayesStats4: mod of extractMRBAYESstats
#######################################################
# extractMrBayesStats4 was copied/moderately modified from 
# phyloch::extractMRBAYESstats
# R package "phyloch", by Christoph Heibl <christoph.heibl@gmx.net>, 
# licensed under GPL (>=2), 
# available at: http://www.christophheibl.de/Rpackages.html

extractMrBayesStats4 <- function (file) 
{
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    phy <- read.nexus(file)
    X <- X[grep("tree TREE1[[:space:]]+=", X)]
    X <- gsub("^.*tree TREE1[[:space:]]+= \\[&R\\] ", "", X)
    X <- gsub("[!]rotate=true,*|[!]rotate=false,*", "", X)
    X <- gsub(";$", "", X)
    vals <- unlist(strsplit(X, "\\][[:alnum:]:.)]*\\[*&*"))
    foo <- function(x) {
        x <- gsub(",([[:lower:]])", "xxx\\1", x)
        x <- unlist(strsplit(x, "xxx"))
        names(x) <- sapply(x, function(x) {
            unlist(strsplit(x, split = "="))[1]
        })
        x <- sapply(x, function(x) {
            unlist(strsplit(x, split = "="))[2]
        })
        x <- gsub("[{}]", "", x)
        x <- strsplit(x, ",")
        lapply(x, as.numeric)
    }
    vals <- lapply(vals, foo)
    edges <- gsub("\\[(&[[:alnum:]_=%!.,{}-]+)\\]", "", X)
    edges <- unlist(strsplit(gsub("[()]*", "", edges), ":"))
    tips <- gsub("^[0-9]+.[0-9]+,*", "", edges)


	re1='(,)'								# Any Single Character 1
	#re2='(?:[a-z][a-z]*[0-9]+[a-z0-9]*)'	# Alphanum 1
	re2='(?:[a-zA-Z0-9]*)'					# Alphanum 1
    regstr = paste(re1, re2, sep="")
    
    #gsub(regstr, "", "asd,234g")


    # remove everything after comma
    # original
    # edges <- as.numeric(gsub("\\[,$|[[:alpha:]].+$\\]", "", edges))
    # modified
    edges <- as.numeric(gsub(pattern=regstr, "", edges))
    
    tab <- cbind(edges, nodes = c("root", head(tips, -1)), id = rep(NA, 
        length(edges)))
    tab[tab[, 2] == "", 2] <- "internal"
    internal <- grep("internal|root", tab[, 2])
    tab[-internal, 3] <- seq(along = phy$tip.label)
    intnodes <- match(tab[internal, 1], phy$edge.length[phy$edge[, 
        2] > 50], nomatch = 0) + 51
    tab[internal, 3] <- intnodes
    names(vals)[internal] <- paste("NODE", tab[internal, 3], 
        sep = ":")
    names(vals)[-internal] <- paste("TIP", tab[-internal, 3], 
        sep = ":")
    tips <- vals[-internal]
    nodes <- vals[internal]
}


#######################################################
# read.beast_original: Reading a BEAST tree, with the stats
#######################################################
# read.beast_original was copied/lightly modified from 
# phyloch::read.beast
# R package "phyloch", by Christoph Heibl <christoph.heibl@gmx.net>, 
# licensed under GPL (>=2), 
# available at: http://www.christophheibl.de/Rpackages.html


# Reading a BEAST tree, with the stats
# from phyloch
setup='
file=confn
digits = NULL
'
read.beast_original <- function (file, digits = NULL, printflag=FALSE) 
	{
	# Scan the files in
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    
    # LEFT is the lines containing [
    # -- basically the tree lines
    LEFT <- grep("\\[", X)
    
    # Extract all of the saved information at each node
    tab <- extractBEASTstats_orig(file, printflag=FALSE)
    if (!is.null(digits)) 
        tab <- round(tab, digits = digits)
        
    # Nodes without posterior values are tip nodes (always 100%!)
    interior <- which(!is.na(tab$posterior))

    # Right is the lines containing ]
    # -- basically the tree lines
    RIGHT <- grep("\\]", X)    
		
	# Convert lines with fully annotated trees into 
	# lines with just simple Newick-format trees
    if (length(LEFT))
    	{
    	
    	# Do ANY of the lines with [ correspond to lines with ]?
        w <- LEFT == RIGHT
        if (any(w))
        	{        	
            s <- LEFT[w]
            
            # For any lines with [ and ], 
            # remove everything between square brackets
            X[s] <- gsub("\\[[^]]*\\]", "", X[s])
        	}
        	
        # For any lines that have ONLY [ or ONLY ]
        # remove anything anything after [ or before ] (?)
        w <- !w
        if (any(w))
        	{
            s <- LEFT[w]
            X[s] <- gsub("\\[.*", "", X[s])
            sb <- RIGHT[w]
            X[sb] <- gsub(".*\\]", "", X[sb])
            if (any(s < sb - 1)) 
                X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
        	}
    	}
    
    # Get line numbers for END; and ENDBLOCK;
    endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
    
    # Get line numbers with semicolons
    semico <- grep(";", X)
    
    # Line # for beginning of TREES block
    i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)

	# Find and use the TRANSLATE block, if it exists
    # Line # for beginning of TRANSLATE block
    i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
    
    # No TRANSLATE block
    if (length(i2) != 0)
    	{    
		# Last line with TRANSLATE items
		end <- semico[semico > i2][1]
		
		# Lines with TRANSLATE items
		x <- X[(i2 + 1):end]
		x <- unlist(strsplit(x, "[,; \t]"))
		x <- x[nzchar(x)]
		
		# Get the translation information
		TRANS <- matrix(x, ncol = 2, byrow = TRUE)
		TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
		n <- dim(TRANS)[1]

		# Start & end lines of the trees block
		start <- semico[semico > i2][1] + 1
		end <- endblock[endblock > i1][1] - 1
		tree <- X[start:end]
		} else {
		# WITHOUT TRANSLATE block
		# Start & end lines of the trees block
		start <- i1+1
		end <- endblock[endblock > i1][1] - 1
		tree <- X[start:end]

		}
    
    
    # Remove everything up to the first =
    # (e.g. tree TREE1 =  
    #  in NEXUS files...)
    tree <- gsub("^.*= *", "", tree)
    
    # Branch lengths are everything after a :
    brl <- unlist(strsplit(tree, ":"))[-1]
    
    # Remove the junk after some branchlengths
    brl <- gsub("[( | ) | ;]", "", brl)
    
    # Remove cases where there was a ,
    brl <- strsplit(brl, ",")
    
    # Take the first item (before a comma)
    foo <- function(x) x <- head(x, 1)
    brl <- unlist(lapply(brl, foo))
    
    # Put the colons back into the list : (?)
    brl <- paste("", brl, sep = ":")
    brl <- c(brl, ";")
    
    # Make an empty vector, 1 cell for each stat
    nodestats <- vector(mode = "list", length = dim(tab)[2])

	# Go through each nodestat
    for (i in seq(along = nodestats))
    #for (i in 1:1)
    	{
        newtree <- tree
        
        # Paste together the nodestat & branchlength
        val <- tab[, i]
        ggg <- paste(val, brl, sep = "")
        ggg[length(ggg)] <- paste(tail(val, 1), ";", sep = "")
        
        # Go through the interior nodes, and replace
        # the branchlength with the node statistic in question
        for (j in interior)
        	{
        	newtree <- gsub(brl[j], ggg[j], newtree)
        	}
        dt <- read.tree(text = newtree)
        
        # Then get these as node labels
        z <- dt$node.label
        z[z == "NA"] <- 9999
        z <- as.numeric(z)
        z[z == 9999] <- NA
        
        # Tabulate the nodestats
        nodestats[[i]] <- z
        names(nodestats)[i] <- colnames(tab)[i]
    	}
    
    
    # Read the input file as plain NEXUS
    tr <- read.nexus(file)
    
    # Add the stats to the standard tree object
    tr <- c(tr[1:length(tr)], nodestats[1:length(nodestats)])
    class(tr) <- ("phylo")
    
    # Add the origin attribute
    attr(tr, "origin") <- file
    tr
	}



#######################################################
# read_beast_prt
#######################################################
# Reading a BEAST tree, with the stats put into a prt()-like table
# 
# Some code for read_beast_prt was copied/inspired/heavily modified from phyloch:
# R package "phyloch", by Christoph Heibl <christoph.heibl@gmx.net>, 
# licensed under GPL (>=2), 
# available at: http://www.christophheibl.de/Rpackages.html
#
# New features include: then putting into a table with prt()
# adding printflags, etc.

# read_beast_prt() produces these columns in the prt table:
# names(beast_nodestats_dtf) = c("internal_nodenums_from1", "internal_nodenums", "rate_range_MIN", "rate_range_MAX", "height_95%_HPD_MIN", "height_95%_HPD_MAX", "height_median", "rate", "height", "rate_median", "height_range_MIN", "height_range_MAX", "rate_95%_HPD_MIN", "rate_95%_HPD_MAX", "posterior", "height_HPD_width", "height_range", "rate_HPD_width", "height_CV", "rate_CV")

read_beast_prt <- function (file, digits = 9, get_tipnames=TRUE, printflag=FALSE, include_rates=FALSE) 
	{
	setup='
	file=confn
	digits = NULL
	get_tipnames=TRUE # If get_tipnames==TRUE, add a row containing the tipnames in the clade, in alphabetical order
	'

	defaults='
	library(BioGeoBEARS)
	sourceall("/drives/GDrive/__github/BEASTmasteR/R/")
	
	# MrBayes doggies (no rates!)
	wd = "/drives/SkyDrive/NIMBioS_projects/2015-03-18_Tumamoc/doggies/mb1/"
	setwd(wd)
	nexfn = "treeLog.mcc"
	
	# Matzke & Irmis, Muller & Reisz 2006 early eureptiles
	wd="/drives/GDrive/__GDrive_projects/2014-11-21_Randy_Irmis_autapomorphies/_03_BEAST_IrmisDates/2016-01-12_Mk_unif_stripFalse"
	setwd(wd)
	nexfn = "treeLog2.mcc"
	file=nexfn
	digits=9
	get_tipnames=TRUE
	printflag=FALSE
	include_rates=TRUE
	'
	
	# Scan the files in
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    
    # LEFT is the lines containing [
    # -- basically the tree lines
    LEFT <- grep("\\[", X)
    
    # Extract all of the saved information at each node
    tab <- extractBEASTstats_orig(file, printflag=FALSE)
    if (!is.null(digits)) 
        tab <- round(tab, digits = digits)
        
    # Nodes without posterior values are tip nodes (always 100%!)
    interior <- which(!is.na(tab$posterior))
	
	# The OTHER nodes are external (tip nodes)
	external <- which(is.na(tab$posterior))
	
	# The lowest height (0) is the root node
	# NO!! Height=0 is the top node, i.e. the present
	# Height=age
	#root_node <-  which(tab$height == 0)
	
	# Remove the root from the external nodes
	#external = external[external != root_node]
	
	
    # Right is the lines containing ]
    # -- basically the tree lines
    RIGHT <- grep("\\]", X)    
		
	# Convert lines with fully annotated trees into 
	# lines with just simple Newick-format trees
    if (length(LEFT))
    	{
    	
    	# Do ANY of the lines with [ correspond to lines with ]?
        w <- LEFT == RIGHT
        if (any(w))
        	{        	
            s <- LEFT[w]
            
            # For any lines with [ and ], 
            # remove everything between square brackets
            X[s] <- gsub("\\[[^]]*\\]", "", X[s])
        	}
        	
        # For any lines that have ONLY [ or ONLY ]
        # remove anything anything after [ or before ] (?)
        w <- !w
        if (any(w))
        	{
            s <- LEFT[w]
            X[s] <- gsub("\\[.*", "", X[s])
            sb <- RIGHT[w]
            X[sb] <- gsub(".*\\]", "", X[sb])
            if (any(s < sb - 1)) 
                X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
        	}
    	}
    
    # Get line numbers for END; and ENDBLOCK;
    endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
    
    # Get line numbers with semicolons
    semico <- grep(";", X)
    
    # Line # for beginning of TREES block
    i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)


 
 	# Find and use the TRANSLATE block, if it exists
    # Line # for beginning of TRANSLATE block
    i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
    
    # If a TRANSLATE BLOCK IS FOUND (in i2), use the first one
    if (length(i2) != 0)
    	{    
		# Last line with TRANSLATE items
		end <- semico[semico > i2][1]
		
		# Lines with TRANSLATE items
		x <- X[(i2 + 1):end]
		x <- unlist(strsplit(x, "[,; \t]"))
		x <- x[nzchar(x)]
		
		# Get the translation information
		TRANS <- matrix(x, ncol = 2, byrow = TRUE)
		TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
		n <- dim(TRANS)[1]

		# Start & end lines of the trees block
		start <- semico[semico > i2][1] + 1
		end <- endblock[endblock > i1][1] - 1
		tree <- X[start:end]
		} else {
		# WITHOUT TRANSLATE block
		# Start & end lines of the trees block
		start <- i1+1
		end <- endblock[endblock > i1][1] - 1
		tree <- X[start:end]
		} # END if (length(i2) != 0)
    

     
    # Remove everything up to the first =
    # (e.g. tree TREE1 =  
    #  in NEXUS files...)
    tree <- gsub("^.*= *", "", tree)
    
    # Branch lengths are everything after a :
    brl <- unlist(strsplit(tree, ":"))[-1]
    
    # Remove the junk after some branchlengths
    brl <- gsub("[( | ) | ;]", "", brl)
    
    # Remove cases where there was a ,
    brl <- strsplit(brl, ",")
    
    # Take the first item (before a comma)
    foo <- function(x) x <- head(x, 1)
    brl <- unlist(lapply(brl, foo))
    
    # Put the colons back into the list : (?)
    brl <- paste("", brl, sep = ":")
    brl <- c(brl, ";")
    
    # Add totally unique 10-character code to brl text
    # (for gsub ability later)
    # (don't do the last brl, it's a ;)
    codenums = sprintf("%010.0f", 1:(length(brl)-1))
    codetxt = paste0("remove", codenums, "remove", collapse=NULL)
    numchars_to_remove = nchar(codetxt[1])
    brl[1:(length(brl)-1)] = paste0(codetxt, brl[1:(length(brl)-1)], collapse=NULL)
    brl
    
    # tree with codetxt
    tree_w_codetxt = tree
    match_colon_locations = gregexpr(":", text=tree_w_codetxt)[[1]]
    startpos=1
    
    insert_text_before_colon <- function(colonposition, textval, original_txt)
    	{
    	original_txt_part1 = substr(x=original_txt, start=1, stop=colonposition-1)
    	original_txt_part2 = substr(x=original_txt, start=colonposition+1, stop=nchar(original_txt))
    	newtxt = paste0(original_txt_part1, textval, ":", original_txt_part2)
    	return(newtxt)
    	}
    
    newtxt = tree_w_codetxt
    # You HAVE to insert things backwards for the match position numbers to keep working
    for (pos in length(match_colon_locations):1)
    	{
    	(newtxt = insert_text_before_colon(colonposition=match_colon_locations[pos], textval=codetxt[pos], original_txt=newtxt))
    	}
    tree_w_codetxt = newtxt
    tree_w_codetxt
    # Make an empty list, 1 item for each stat
    # nodestats: for internal nodes
    nodestats <- vector(mode = "list", length = dim(tab)[2])
	tipstats <- vector(mode = "list", length = dim(tab)[2])
	
	# Go through each nodestat
    for (i in seq(along = nodestats))
    #for (i in 1:1)
    	{
        newtree <- tree_w_codetxt
        
        # Paste together the nodestat & branchlength
        val <- tab[, i]
        ggg <- paste(val, brl, sep = "")
        ggg[length(ggg)] <- paste(tail(val, 1), ";", sep = "")
        
        # Edit the ":" text in the newtree to be e.g.
        # :remove0000000001remove
        
        
        
        # Go through the interior nodes, and replace
        # the branchlength with the node statistic in question
        for (j in interior)
        	{
        	newtree <- gsub(brl[j], ggg[j], newtree)
        	}
        # External nodes
        for (j in external)
        	{
        	# Add a delimiter for the labels below the tipnames
        	newtree <- gsub(brl[j], paste0("|", ggg[j]), newtree)
        	}
        # Root node (actually already in the tipnodes
        #for (j in root_node)
        #	{
        #	newtree <- gsub(brl[j], ggg[j], newtree)
        #	}
        # I think you DON'T have to remove
        # numchars_to_remove, since you've
        # already gsubbed it out!
        
        # Now remove that codetext
        for (removei in 1:length(codetxt))
        	{
        	newtree = gsub(pattern=codetxt[removei], replacement="", x=newtree)
        	}
        newtree
        dt <- read.tree(text = newtree)
        dt
        
        
        # Then get these as node labels
        z <- dt$node.label
        z[z == "NA"] <- 9999
        z <- as.numeric(z)
        z[z == 9999] <- NA
        
        # Also get the tiplabels
        tipvals = as.numeric(matrix(data=unlist(strsplit(dt$tip.label, split="\\|")), ncol=2, byrow=TRUE)[,2])
        tipvals
        
        # Tabulate the nodestats
        nodestats[[i]] <- z
        names(nodestats)[i] <- colnames(tab)[i]

        tipstats[[i]] <- tipvals
        names(tipstats)[i] <- colnames(tab)[i]
    	}
    
    tipstats
    
    # Read the input file as plain NEXUS
    tr <- read.nexus(file)
    
    # Add the stats (a bunch of sub-objects) to the standard 
    # tree object (4 sub-objects)
    tr <- c(tr[1:length(tr)], nodestats[1:length(nodestats)])
    class(tr) <- ("phylo")
	
    # Add the origin attribute
    attr(tr, "origin") <- file
    tr
    
    
    #######################################################
	# Get the prt table    
	#######################################################
    tr_table = prt(tr, printflag=FALSE, relabel_nodes = FALSE, time_bp_digits=7, get_tipnames=get_tipnames)
    
    # make the internal node numbers
    tipnums = 1:length(tr$tip.label)
    tipnums_from1 = rep(NA, times=length(tr$tip.label))
    internal_nodenums = seq(from=length(tr$tip.label)+1, to=length(tr$tip.label)+tr$Nnode, by=1)
	internal_nodenums_from1 = seq(from=1, to=tr$Nnode, by=1)

    # Make a table of internal node stats from the BEAST stats
    height_HPD_width = tr$"height_95%_HPD_MAX" - tr$"height_95%_HPD_MIN"
    height_range = tr$"height_range_MAX" - tr$"height_range_MIN"
    rate_HPD_width = tr$"rate_95%_HPD_MAX" - tr$"rate_95%_HPD_MIN"
    
    # Calculate the CVs of height and rate
    height_CV = (height_HPD_width/(1.96*2)) / tr$height
    rate_CV = (rate_HPD_width/(1.96*2)) / tr$rate
    
     # Make a table of tip stats from the BEAST stats
    tipstats$height_HPD_width = tipstats$"height_95%_HPD_MAX" - tipstats$"height_95%_HPD_MIN"
    tipstats$height_range = tipstats$"height_range_MAX" - tipstats$"height_range_MIN"
    tipstats$rate_HPD_width = tipstats$"rate_95%_HPD_MAX" - tipstats$"rate_95%_HPD_MIN"
    
    # Calculate the tip stats CVs of height and rate
    tipstats$height_CV = (tipstats$height_HPD_width/(1.96*2)) / tipstats$height
    tipstats$rate_CV = (tipstats$rate_HPD_width/(1.96*2)) / tipstats$rate
   
    
    if (include_rates == FALSE)
    	{
		beast_nodestats_table_temp = cbind(internal_nodenums_from1, internal_nodenums, tr$"height_95%_HPD_MIN", tr$"height_95%_HPD_MAX", tr$height_median, tr$height, tr$height_range_MIN, tr$height_range_MAX, tr$posterior, height_HPD_width, height_range, height_CV)

		# Add tip rows as NAs -- NO! GET THE DATA IN THERE!!
		#tiprows = matrix(data=NA, nrow=length(tr$tip.label), ncol=ncol(beast_nodestats_table_temp))
		tiprows = cbind(tipnums_from1, tipnums, tipstats$"height_95%_HPD_MIN", tipstats$"height_95%_HPD_MAX", tipstats$height_median, tipstats$height, tipstats$height_range_MIN, tipstats$height_range_MAX, tipstats$posterior, tipstats$height_HPD_width, tipstats$height_range, tipstats$height_CV)
	
		beast_nodestats_table = rbind(tiprows, beast_nodestats_table_temp)
		beast_nodestats_dtf = adf2(beast_nodestats_table)
		names(beast_nodestats_dtf) = c("internal_nodenums_from1", "internal_nodenums", "height_95%_HPD_MIN", "height_95%_HPD_MAX", "height_median", "height", "height_range_MIN", "height_range_MAX", "posterior", "height_HPD_width", "height_range", "height_CV")
		} else {
		beast_nodestats_table_temp = cbind(internal_nodenums_from1, internal_nodenums, tr$rate_range_MIN, tr$rate_range_MAX, tr$"height_95%_HPD_MIN", tr$"height_95%_HPD_MAX", tr$height_median, tr$rate, tr$height, tr$rate_median, tr$height_range_MIN, tr$height_range_MAX, tr$"rate_95%_HPD_MIN", tr$"rate_95%_HPD_MAX", tr$posterior, height_HPD_width, height_range, rate_HPD_width, height_CV, rate_CV)

		# Add tip rows as NAs -- NO! GET THE DATA IN THERE!!
		# 2016-06-02_bug: the RATES WILL NOT be NA!
		#tiprows = matrix(data=NA, nrow=length(tr$tip.label), ncol=ncol(beast_nodestats_table_temp))
		tiprows = cbind(tipnums_from1, tipnums, tipstats$rate_range_MIN, tipstats$rate_range_MAX, tipstats$"height_95%_HPD_MIN", tipstats$"height_95%_HPD_MAX", tipstats$height_median, tipstats$rate, tipstats$height, tipstats$rate_median, tipstats$height_range_MIN, tipstats$height_range_MAX, tipstats$"rate_95%_HPD_MIN", tipstats$"rate_95%_HPD_MAX", tipstats$posterior, tipstats$height_HPD_width, tipstats$height_range, tipstats$rate_HPD_width, tipstats$height_CV, tipstats$rate_CV)
		
		beast_nodestats_table = rbind(tiprows, beast_nodestats_table_temp)
		beast_nodestats_dtf = adf2(beast_nodestats_table)
		names(beast_nodestats_dtf) = c("internal_nodenums_from1", "internal_nodenums", "rate_range_MIN", "rate_range_MAX", "height_95%_HPD_MIN", "height_95%_HPD_MAX", "height_median", "rate", "height", "rate_median", "height_range_MIN", "height_range_MAX", "rate_95%_HPD_MIN", "rate_95%_HPD_MAX", "posterior", "height_HPD_width", "height_range", "rate_HPD_width", "height_CV", "rate_CV")
		}
    
    
	# head(beast_nodestats_dtf)
	# tail(beast_nodestats_dtf)
    
    # NAs in these tips represent tips with 0 variation in height. Convert NA to 0:
    TF = is.na(beast_nodestats_dtf$"height_95%_HPD_MIN"[tipnums])
    beast_nodestats_dtf$"height_95%_HPD_MIN"[tipnums][TF] = 0

    TF = is.na(beast_nodestats_dtf$"height_95%_HPD_MAX"[tipnums])
    beast_nodestats_dtf$"height_95%_HPD_MAX"[tipnums][TF] = 0

    TF = is.na(beast_nodestats_dtf$height_median[tipnums])
    beast_nodestats_dtf$height_median[tipnums][TF] = 0

    TF = is.na(beast_nodestats_dtf$height_range_MIN[tipnums])
    beast_nodestats_dtf$height_range_MIN[tipnums][TF] = 0
    
    TF = is.na(beast_nodestats_dtf$height_range_MAX[tipnums])
    beast_nodestats_dtf$height_range_MAX[tipnums][TF] = 0
    
    TF = is.na(beast_nodestats_dtf$height_HPD_width[tipnums])
    beast_nodestats_dtf$height_HPD_width[tipnums][TF] = 0
    
    TF = is.na(beast_nodestats_dtf$height_range[tipnums])
    beast_nodestats_dtf$height_range[tipnums][TF] = 0
    
    TF = is.na(beast_nodestats_dtf$height_CV[tipnums])
    beast_nodestats_dtf$height_CV[tipnums][TF] = 0
    beast_nodestats_dtf
    
    # Assemble the final table
    prt_beast_nodestats = cbind(tr_table, beast_nodestats_dtf)
    
	# Return a list with the tree and the megatable
	beastcon = NULL
	beastcon$tr = tr
	beastcon$prt_beast_nodestats = prt_beast_nodestats
	
	return(beastcon)
	}








#######################################################
# read_beasttree_rates: Reading a BEAST tree, with the stats
#######################################################

# Some code for read_beasttree_rates was copied/inspired/heavily modified from phyloch:
# R package "phyloch", by Christoph Heibl <christoph.heibl@gmx.net>, 
# licensed under GPL (>=2), 
# available at: http://www.christophheibl.de/Rpackages.html

read_beasttree_rates <- function (trfn, digits = NULL) 
	{
	setup='
	trfn="temp.nexus"
	digits = NULL
	'

	# Scan the trfns in
    X <- scan(file = trfn, what = "", sep = "\n", quiet = TRUE)
    X_orig = X
    
    # LEFT is the lines containing [
    # -- basically the tree lines
    lines_with_brackets = grep("\\[", X)
    LEFT <- grep("\\[", X)
    
    # Extract all of the saved information at each node
    tab <- extractBEASTstats2(trfn)
    if (!is.null(digits)) 
        tab <- round(tab, digits = digits)

	#interior <- which(!is.na(tab$posterior))
    
    # Right is the lines containing ]
    # -- basically the tree lines
	RIGHT <- grep("\\]", X)    
	
	
	# Convert lines with fully annotated trees into 
	# lines with just simple Newick-format trees
    if (length(LEFT) > 0)
    	{
    	
    	# Do ANY of the lines with [ correspond to lines with ]?
        w <- LEFT == RIGHT
        if (any(w))
        	{        	
            s <- LEFT[w]
            
            # For any lines with [ and ], 
            # remove everything between square brackets
            X[s] <- gsub("\\[[^]]*\\]", "", X[s])
        	}
        	
        # For any lines that have ONLY [ or ONLY ]
        # remove anything anything after [ or before ] (?)
        w <- !w
        if (any(w))
        	{
            s <- LEFT[w]
            X[s] <- gsub("\\[.*", "", X[s])
            sb <- RIGHT[w]
            X[sb] <- gsub(".*\\]", "", X[sb])
            if (any(s < sb - 1)) 
                X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
                
                # also remove from X_org
                X_orig <- X_orig[-unlist(mapply(":", (s + 1), (sb - 1)))]
                
        	}
    	}
    
    # Get line numbers for END; and ENDBLOCK;
   	endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
	
	# Get line numbers with semicolons
	semico <- grep(";", X)
	
	# Line # for beginning of TREES block
	i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)

	# Line # for beginning of TRANSLATE block
	i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
	
	# Last line with TRANSLATE items
	end <- semico[semico > i2][1]
	
	# Lines with TRANSLATE items
	x <- X[(i2 + 1):end]
	x <- unlist(strsplit(x, "[,; \t]"))
	x <- x[nzchar(x)]
	
	# Get the translation information
	TRANS <- matrix(x, ncol = 2, byrow = TRUE)
	TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
	n <- dim(TRANS)[1]
	
	# Start & end lines of the trees block
	start <- semico[semico > i2][1] + 1
	end <- endblock[endblock > i1][1] - 1
	tree <- X[start:end]

	
	# Otherwise do this...
	#  If there isn't an endblock, then the last tree is the last line of interest
	#line_of_last_tree = 
	
	endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
	
	# Get line numbers with semicolons
	semico <- grep(";", X)
	
	# Line # for beginning of TREES block
	i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)

	# Line # for beginning of TRANSLATE block
	i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
	
	# Last line with TRANSLATE items
	end <- semico[semico > i2][1]
	
	# Lines with TRANSLATE items
	x <- X[(i2 + 1):end]
	x <- unlist(strsplit(x, "[,; \t]"))
	x <- x[nzchar(x)]
	
	# Get the translation information
	TRANS <- matrix(x, ncol = 2, byrow = TRUE)
	TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
	n <- dim(TRANS)[1]
	
	# Start & end lines of the trees block
	start <- semico[semico > i2][1] + 1
	end <- endblock[endblock > i1][1] - 1
	tree <- X[start:end]
	
    # Remove everything up to the first =
    # (e.g. tree TREE1 =  
    #  in NEXUS trfns...)
    tree <- gsub("^.*= *", "", tree)
    

	#====================================
	# Get branch lengths in text trfn order    
	#====================================
    # Branch lengths are everything after a :
    brl <- unlist(strsplit(tree, ":"))[-1]
    
    # Remove the junk after some branchlengths
    brl <- gsub("[( | ) | ;]", "", brl)
    
    # Remove cases where there was a ,
    brl <- strsplit(brl, ",")
    
    # Take the first item (before a comma)
    foo <- function(x) x <- head(x, 1)
    brl <- unlist(lapply(brl, foo))
    
    # Put the colons back into the list : (?)
    brl <- paste("", brl, sep = ":")
    brl <- c(brl, ";")
    ########################################
    
    ########################################
    # Make an empty vector, 1 cell for each stat
    nodestats <- vector(mode = "list", length = dim(tab)[2])
    ########################################


	# Go through ALL the nodes, and replace
	# the branchlength with the rate statistic
	#
	# This produces a rateogram

	rateogram <- tree
	
	for (j in 1:nrow(tab))
		{		
		new_branch_length = paste(":", as.character(tab[j,]), sep="")
		
		# Replace JUST THE FIRST match
		matches = gregexpr(brl[j], rateogram)[[1]]
		matches_end = matches-1+attr(matches,"match.length")
		
		newtree_chars = strsplit(rateogram, "")[[1]]
		firstchars = newtree_chars[1:(matches[1]-1)]
		lastchars = newtree_chars[(matches_end[1]+1):length(newtree_chars)]
		
		middle_chars = new_branch_length
		
		
		# Merge back together
		rateogram = paste(c(firstchars, middle_chars, lastchars), sep="", collapse="")
		
		#newtree <- gsub(brl[j], new_branch_length, newtree)
		}
	rateogram_tr <- read.tree(text = rateogram)
	orig_tr <- read.nexus(trfn)
	
	# Put the edge lengths into the original tree
	orig_tr$edge.length = rateogram_tr$edge.length
	rateogram_tr = orig_tr
	
	
	return(rateogram_tr)
	}







# Some code for read_beasttree_states was copied/inspired/heavily modified from phyloch:
# R package "phyloch", by Christoph Heibl <christoph.heibl@gmx.net>, 
# licensed under GPL (>=2), 
# available at: http://www.christophheibl.de/Rpackages.html

read_beasttree_states <- function (trfn, digits = NULL) 
	{
	# Scan the trfns in
    X <- scan(file = trfn, what = "", sep = "\n", quiet = TRUE)
    X_orig = X
    
    # LEFT is the lines containing [
    # -- basically the tree lines
    lines_with_brackets = grep("\\[", X)
    LEFT <- grep("\\[", X)
    
    # Extract all of the saved information at each node
    tab <- extractBEASTstats2(trfn)
    if (!is.null(digits)) 
        tab <- round(tab, digits = digits)

	#interior <- which(!is.na(tab$posterior))
    
    # Right is the lines containing ]
    # -- basically the tree lines
	RIGHT <- grep("\\]", X)    
	
	
	# Convert lines with fully annotated trees into 
	# lines with just simple Newick-format trees
    if (length(LEFT) > 0)
    	{
    	
    	# Do ANY of the lines with [ correspond to lines with ]?
        w <- LEFT == RIGHT
        if (any(w))
        	{        	
            s <- LEFT[w]
            
            # For any lines with [ and ], 
            # remove everything between square brackets
            X[s] <- gsub("\\[[^]]*\\]", "", X[s])
        	}
        	
        # For any lines that have ONLY [ or ONLY ]
        # remove anything anything after [ or before ] (?)
        w <- !w
        if (any(w))
        	{
            s <- LEFT[w]
            X[s] <- gsub("\\[.*", "", X[s])
            sb <- RIGHT[w]
            X[sb] <- gsub(".*\\]", "", X[sb])
            if (any(s < sb - 1)) 
                X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
                
                # also remove from X_org
                X_orig <- X_orig[-unlist(mapply(":", (s + 1), (sb - 1)))]
                
        	}
    	}
    
    # Get line numbers for END; and ENDBLOCK;
   	endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
	
	# Get line numbers with semicolons
	semico <- grep(";", X)
	
	# Line # for beginning of TREES block
	i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)

	# Line # for beginning of TRANSLATE block
	i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
	
	# Last line with TRANSLATE items
	end <- semico[semico > i2][1]
	
	# Lines with TRANSLATE items
	x <- X[(i2 + 1):end]
	x <- unlist(strsplit(x, "[,; \t]"))
	x <- x[nzchar(x)]
	
	# Get the translation information
	TRANS <- matrix(x, ncol = 2, byrow = TRUE)
	TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
	n <- dim(TRANS)[1]
	
	# Start & end lines of the trees block
	start <- semico[semico > i2][1] + 1
	end <- endblock[endblock > i1][1] - 1
	tree <- X[start:end]

	
	# Otherwise do this...
	#  If there isn't an endblock, then the last tree is the last line of interest
	#line_of_last_tree = 
	
	endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
	
	# Get line numbers with semicolons
	semico <- grep(";", X)
	
	# Line # for beginning of TREES block
	i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)

	# Line # for beginning of TRANSLATE block
	i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
	
	# Last line with TRANSLATE items
	end <- semico[semico > i2][1]
	
	# Lines with TRANSLATE items
	x <- X[(i2 + 1):end]
	x <- unlist(strsplit(x, "[,; \t]"))
	x <- x[nzchar(x)]
	
	# Get the translation information
	TRANS <- matrix(x, ncol = 2, byrow = TRUE)
	TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
	n <- dim(TRANS)[1]
	
	# Start & end lines of the trees block
	start <- semico[semico > i2][1] + 1
	end <- endblock[endblock > i1][1] - 1
	tree <- X[start:end]
	
    # Remove everything up to the first =
    # (e.g. tree TREE1 =  
    #  in NEXUS trfns...)
    tree <- gsub("^.*= *", "", tree)
    

	#====================================
	# Get branch lengths in text trfn order    
	#====================================
    # Branch lengths are everything after a :
    brl <- unlist(strsplit(tree, ":"))[-1]
    
    # Remove the junk after some branchlengths
    brl <- gsub("[( | ) | ;]", "", brl)
    
    # Remove cases where there was a ,
    brl <- strsplit(brl, ",")
    
    # Take the first item (before a comma)
    foo <- function(x) x <- head(x, 1)
    brl <- unlist(lapply(brl, foo))
    
    # Put the colons back into the list : (?)
    brl <- paste("", brl, sep = ":")
    brl <- c(brl, ";")
    ########################################
    
    ########################################
    # Make an empty vector, 1 cell for each stat
    nodestats <- vector(mode = "list", length = dim(tab)[2])
    ########################################

	return(nodestats)
	}



#######################################################
# extractBEASTstats2: Another mod of phyloch::extractBEASTstats
#######################################################
# Some code for extractBEASTstats2 was copied/inspired/heavily modified from phyloch:
# R package "phyloch", by Christoph Heibl <christoph.heibl@gmx.net>, 
# licensed under GPL (>=2), 
# available at: http://www.christophheibl.de/Rpackages.html

setup='
fn = trfn
fn = confn
regexp = "(tree)( )(STATE)(_)(\\d+)"
'
extractBEASTstats2 <- function (fn, regexp = "(tree)( )(STATE)(_)(\\d+)") 
	{
	
	# Scan in the NEXUS file
    X1 <- scan(file = fn, what = "", sep = "\n", quiet = TRUE)
    
    # phyloch original search string (assumes consensus tree)
    #Y <- X[grep("tree [[:space:]]+=", X)]
	
	# use your desired regular expression for the beginning of each tree line
	# The default gets e.g. "tree STATE_2323092"
	#regexp = "(tree)( )(STATE)(_)(\\d+)"

    X2 <- X1[grep(regexp, X1)]
        
    #Y <- X[grep("tree", X)]
    
    
    tmpre = NULL
	ri = 0
	tmpre[[(ri=ri+1)]] = '(tree)'	# Word 1
	tmpre[[(ri=ri+1)]] = '( )'	# White Space 1
	tmpre[[(ri=ri+1)]] = '(STATE)'	# Word 2
	tmpre[[(ri=ri+1)]] = '(_)'	# Any Single Character 1
	tmpre[[(ri=ri+1)]] = '(\\d+)'	# Integer Number 1
	tmpre[[(ri=ri+1)]] = '( )'	# White Space 2
	tmpre[[(ri=ri+1)]] = '(\\[.*?\\])'	# Square Braces 1
	tmpre[[(ri=ri+1)]] = '( )'	# White Space 3
	tmpre[[(ri=ri+1)]] = '(=)'	# Any Single Character 2
	tmpre[[(ri=ri+1)]] = '( )'	# White Space 4
	tmpre[[(ri=ri+1)]] = '(\\[&R\\])'	# Square Braces 2
	tmpre[[(ri=ri+1)]] = '( )'	# White Space 5
	
	tmpre_string = list2str(tmpre, spacer="")
	(tmpre_string)

    # Remove the header   
    X3 <- gsub(tmpre_string, "", X2)
    
	# Extract the tree header
	treeheader_txt = extract_regexp(tmpstr=X2, tmpre_string)
    
    # Parse the tree header
    words = strsplit(treeheader_txt, split=" ")[[1]]
    
    # get the STATE
    tmpstate = gsub(pattern="STATE_", replacement="", words[2])
    
    # Extract everything within square braces
	regexp_braces = "(\\[.*?\\])"
    brackets_strings = extract_regexp(tmpstr=X3, tmpre_string=regexp_braces)
    names(brackets_strings) = NULL
    brackets_strings = unlist(brackets_strings)
	smm(brackets_strings)
	
	# Get the strings after each "["    
    #tab <- unlist(strsplit(X3, "\\["))[-1]
    brackets_strings1 = gsub("\\[&rate=", "", brackets_strings)
    brackets_strings2 = gsub("\\]", "", brackets_strings1)
    
    
    tab = as.numeric(brackets_strings2)
    
    
    
    # PARSE BRACKETS STRINGS & THIS WILL WORK
    
    
    
    tab = matrix(tab, ncol=1)
    tab = adf(tab)
    names(tab) = c("rates")
    
    #tab2 <- gsub("&|;|\\]", "", tab)
    #tab <- gsub(":.+$", "", tab)
    #foo <- function(x)
    #	{
    #    x <- unlist(strsplit(x, ","))
    #    x
   	# 	}
    #tab <- lapply(tab, foo)
    #for (i in seq(along = tab))
    #	{
    #    ind <- grep("[{]", tab[[i]])
    #    names <- gsub("=.+$", "", tab[[i]][ind])
    #    tab[[i]][ind] <- gsub("[{]", "", tab[[i]][ind])
    #    tab[[i]][ind] <- gsub("=", "_MIN=", tab[[i]][ind])
    #    tab[[i]][ind + 1] <- gsub("[}]", "", tab[[i]][ind + 1])
    #    tab[[i]][ind + 1] <- paste(paste(names, "MAX=", sep = "_"), tab[[i]][ind + 1])
	#    }
    #ttab <- data.frame()
    #stats <- unique(gsub("=.+$", "", unlist(tab)))
    #for (i in seq(along = tab)) 
    #	{
     #   for (j in seq(along = stats))
      #  	{
       #     ind <- grep(paste("^", stats[j], "=", sep = ""), 
        #        tab[[i]])
         #   if (length(ind) > 0)
          #  	{
           #     v <- as.numeric(gsub(paste(stats[j], "=", sep = ""), 
            #      "", tab[[i]][ind]))
             #   ttab[i, j] <- v
            	#}
   #     	}
    #	}
    #colnames(ttab) <- stats
    #tip <- which(is.na(ttab$posterior))
    return(tab)
	}



#######################################################
# From phyloch (released with GPL >2)
# MUCH faster than read.nexus.data or read_nexus_data2
# (but, way less flexible)
#######################################################

read_nex_phyloch <- function(fn){
	
	x <- scan(fn, what = "c", quiet = TRUE)
		
	## eliminate comments
	## ------------------
	left <- grep("\\[", x)
	right <- grep("\\]", x)
	if ( length(left) > 0 ){
	  m <- cbind(left, right)
	  x <- x[-unlist(apply(m, 1, function(x) x[1]:x[2]))]
	}
	
	x <- x[x != ""]
	
  ## getting number of taxa
  ## ----------------------
	ntax <- x[grep("ntax", x, ignore.case = TRUE)]
	ntax <- gsub("[[:alpha:]]|[[:punct:]]", "", ntax )
	nb <- ntax <- as.numeric(unique(ntax))
		
	## getting number of characters	
  ## ----------------------------
	ncha <- x[grep("nchar", x, ignore.case = TRUE)]
	ncha <- gsub("[[:alpha:]]|[[:punct:]]", "", ncha )
	ncha <- as.numeric(unique(ncha))
	
	## get beginning and end of matrix
  ## -------------------------------
	start <- grep("^\t?matrix$", x, ignore.case = TRUE)
	end <- grep(";", x)
	end <- min(end[end > start])
	M <- x[(start + 1):(end - 1)]
	
	# assemble DNAbin object:
	# -----------------------
	nblock <- ceiling(ncha / nchar(M[2]))
	id <- seq(1, 2 * ntax, by = 2)
	nam <- M[id]
	fuse <- function(s, M, nblock, ntax){
	  paste(M[seq(s, length.out = nblock, by = ntax * 2)], collapse = "")
	}
	seq <- lapply(id + 1, fuse, M = M, nblock = nblock, ntax = ntax)
	obj <- list(nb = ntax, seq = seq, nam = nam, com = NA)
	class(obj) <- "alignment"
	as.DNAbin(obj)
} # END read.nex


#######################################################
# Much FASTER read of a large NEXUS DNA file
#######################################################
read_nex_phyloch_to_list <- function(fn)
	{
	obj = read_nex_phyloch(fn)
	charmat = as.character(obj)
	
	charlist = list()
	for (i in 1:nrow(charmat))
		{
		charlist[[i]] = charmat[i,]
		}
	names(charlist) = row.names(charmat)
	return(charlist)
	}



