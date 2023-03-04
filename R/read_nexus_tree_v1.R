#######################################################
# read_nexus2: Better version of read.nexus
#
# ie avoiding a bug about
# "numbers of left and right parentheses in Newick string not equal"
#
# Also, can read in node data
#
#######################################################



# For Debugging NEXUS input -- see TNT R stuff 
# for better functions.
defaults='
file=outfn2
tree.names=NULL
force.multi=FALSE
'

# read_nexus_data2 is based on APE's "read.data.nexus", but has been 
# substantially modified to take into account morphology data, common 
# ambiguities, etc.

# Slight modification of the APE function read.data.nexus, to allow the input of characters
# listed as DATATYPE=STANDARD
# (defaults allow only "DNA" and "PROTEIN")
#
#' Obj The output is a "charslist": a list of species, each with a 
#' list of textual character states. This is the same format as 
#' produced by read_nexus_data2(). A charlist can be converted to a 
#' charsdf (a data.frame) with as.data.frame(x=..., stringsAsFactors=FALSE).
#' The reverse operation can occur with as.list().

read_nexus2 <- function (file, tree.names=NULL, force.multi=FALSE) 
	{
	X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
	LEFT <- grep("\\[", X)
	RIGHT <- grep("\\]", X)
	if (length(LEFT)) {
	    w <- LEFT == RIGHT
	    if (any(w)) {
	        s <- LEFT[w]
	        X[s] <- gsub("\\[[^]]*\\]", "", X[s])
	    }
	    w <- !w
	    if (any(w)) {
	        s <- LEFT[w]
	        X[s] <- gsub("\\[.*", "", X[s])
	        sb <- RIGHT[w]
	        X[sb] <- gsub(".*\\]", "", X[sb])
	        if (any(s < sb - 1)) 
	            X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
	    }
	}
	endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
	semico <- grep(";", X)
	i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
	i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
	if (length(i2) == 1 && i2 > i1)
		{
	  translation = TRUE
	  } else {
	  translation = FALSE
	  }
	if (translation == TRUE)
		{
		end <- semico[semico > i2][1]
		x <- X[(i2 + 1):end]
		x <- gsub("^\\s+", "", x)
		x <- gsub("[,;]", "", x)
		x <- unlist(regmatches(x, regexpr("\\s+", x), invert = TRUE))
		x <- x[nzchar(x)]
		TRANS <- matrix(x, ncol = 2, byrow = TRUE)
		TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
		n <- dim(TRANS)[1]
		}
	if (translation)
		{
		start = semico[semico > i2][1] + 1
		} else {
		start = i1 + 1
		}
	end <- endblock[endblock > i1][1] - 1
	tree <- X[start:end]
	rm(X)
	tree <- tree[tree != ""]
	semico <- grep(";", tree)
	Ntree <- length(semico)
	if (Ntree == 1 && length(tree) > 1)
		{
		STRING <- paste(tree, collapse = "")
		} else {
		if (any(diff(semico) != 1))
			{
			STRING <- character(Ntree)
			s <- c(1, semico[-Ntree] + 1)
			j <- mapply(":", s, semico)
			if (is.list(j))
				{
				for (i in 1:Ntree)
					{
					STRING[i] <- paste(tree[j[[i]]], collapse = "")
					}
				} else {
				for (i in 1:Ntree)
					{
					STRING[i] <- paste(tree[j[, i]], collapse = "")
					}
				} # END if (is.list(j))
	    } else {
	    STRING <- tree
	    } # END if (any(diff(semico) != 1))
		} # END if (Ntree == 1 && length(tree) > 1)
	rm(tree)
	STRING <- STRING[grep("^[[:blank:]]*tree.*= *", STRING, ignore.case = TRUE)]
	Ntree <- length(STRING)
	#  THIS CAN CRASH, IF THERE IS AN "=" IN THE TIP NAMES
	nms.trees <- sub(" *= *.*", "", STRING)
	nms.trees <- sub("^[[:blank:]]*tree[[:blank:]\\*]*", "", 
	    nms.trees, ignore.case = TRUE)
	STRING <- sub("^.*= *", "", STRING)
	STRING <- gsub(" ", "", STRING)
	colon <- grep(":", STRING)
	if (!length(colon))
		{
	  trees <- lapply(STRING, .cladoBuild)
		} else if (length(colon) == Ntree) {
	  if (translation)
	  	{
	  	trees = lapply(STRING, .treeBuildWithTokens)
	  	} else {
	  	trees = lapply(STRING, .treeBuild)
	  	} # END if (translation)
		} else {
	    trees <- vector("list", Ntree)
	    trees[colon] <- lapply(STRING[colon], .treeBuild)
	    nocolon <- (1:Ntree)[!1:Ntree %in% colon]
	    trees[nocolon] <- lapply(STRING[nocolon], .cladoBuild)
	    if (translation) {
	        for (i in 1:Ntree) {
	            tr <- trees[[i]]
	            for (j in 1:n) {
	              ind <- which(tr$tip.label[j] == TRANS[, 1])
	              tr$tip.label[j] <- TRANS[ind, 2]
	            }
	            if (!is.null(tr$node.label)) {
	              for (j in 1:length(tr$node.label)) {
	                ind <- which(tr$node.label[j] == TRANS[, 
	                  1])
	                tr$node.label[j] <- TRANS[ind, 2]
	              }
	            }
	            trees[[i]] <- tr
	        }
	        translation <- FALSE
	    }
	}
	for (i in 1:Ntree) {
	    tr <- trees[[i]]
	    if (!translation) 
	        n <- length(tr$tip.label)
	}
	if (Ntree == 1 && !force.multi) {
	    trees <- trees[[1]]
	    if (translation) {
	        trees$tip.label <- if (length(colon)) 
	            TRANS[, 2]
	        else TRANS[, 2][as.numeric(trees$tip.label)]
	    }
	}
	else {
	    if (!is.null(tree.names)) 
	        names(trees) <- tree.names
	    if (translation) {
	        if (length(colon) == Ntree) 
	            attr(trees, "TipLabel") <- TRANS[, 2]
	        else {
	            for (i in 1:Ntree) trees[[i]]$tip.label <- TRANS[, 
	              2][as.numeric(trees[[i]]$tip.label)]
	            trees <- .compressTipLabel(trees)
	        }
	    }
	    class(trees) <- "multiPhylo"
	    if (!all(nms.trees == "")) 
	        names(trees) <- nms.trees
	}
	trees
}