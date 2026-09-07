# 2026-09-07

# Modified from phytools:::readNexusData to put nmax default into scan(),
# to stop scan() from hanging

# phytools::readNexus
phytools_readNexus2 <- function(file = "", format = c("standard", "raxml"), nmax=1000000) 
{
    format <- format[1]
    if (tolower(format) == "standard") 
        tree <- ape::read.nexus(file)
    else if (tolower(format) == "raxml") {
        XX <- phytools_readNexusData2(file, version = 3.5, nmax=nmax)
        text <- XX$text
        trans <- XX$trans
        Ntree <- XX$Ntree
        tree <- lapply(text, phytools_modified.text_to_tree2, trans = trans)
        if (length(tree) == 1) 
            tree <- tree[[1]]
        else class(tree) <- "multiPhylo"
    }
    else {
        cat("Do not recognize format\n")
        tree <- NULL
    }
    tree
} # END phytools_readNexus2 <- function(file = "", format = c("standard", "raxml"), nmax=1000000) 


# Modified from phytools:::readNexusData to put nmax default into scan(),
# to stop scan() from hanging
phytools_readNexusData2 <- function(file, version, nmax=1000000) 
{
    if (version <= 1) {
        X <- scan(file = file, what = "", sep = "\n", quiet = TRUE, nmax=nmax)
        left <- grep("\\[", X)
        right <- grep("\\]", X)
        if (length(left)) {
            w <- left == right
            if (any(w)) {
                s <- left[w]
                X[s] <- gsub("\\[[^]]*\\]", "", X[s])
            }
            w <- !w
            if (any(w)) {
                s <- left[w]
                X[s] <- gsub("\\[.*", "", X[s])
                sb <- right[w]
                X[sb] <- gsub(".*\\]", "", X[sb])
                if (any(s < sb - 1)) 
                  X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
            }
        }
        endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
        semico <- grep(";", X)
        i1 <- grep("begin smptrees;", X, ignore.case = TRUE)
        i2 <- grep("translate", X, ignore.case = TRUE)
        translation <- if (length(i2) == 1 && i2 > i1) 
            TRUE
        else FALSE
        if (translation) {
            end <- semico[semico > i2][1]
            x <- X[(i2 + 1):end]
            x <- unlist(strsplit(x, "[,;\t]"))
            x <- unlist(strsplit(x, " "))
            x <- x[nzchar(x)]
            trans <- matrix(x, ncol = 2, byrow = TRUE)
            trans[, 2] <- gsub("['\"]", "", trans[, 2])
            n <- dim(trans)[1]
        }
        start <- if (translation) 
            semico[semico > i2][1] + 1
        else semico[semico > i1][1]
        end <- endblock[endblock > i1][1] - 1
        tree <- X[start:end]
        tree <- gsub("^.*= *", "", tree)
        tree <- tree[tree != ""]
        semico <- grep(";", tree)
        Ntree <- length(semico)
        if (Ntree == 1 && length(tree) > 1) {
            STRING <- paste(tree, collapse = "")
        }
        else {
            if (any(diff(semico) != 1)) {
                STRING <- character(Ntree)
                s <- c(1, semico[-Ntree] + 1)
                j <- mapply(":", s, semico)
                if (is.list(j)) {
                  for (i in 1:Ntree) STRING[i] <- paste(tree[j[[i]]], 
                    collapse = "")
                }
                else {
                  for (i in 1:Ntree) STRING[i] <- paste(tree[j[, 
                    i]], collapse = "")
                }
            }
            else STRING <- tree
        }
        text <- STRING
        if (translation == TRUE) {
            rownames(trans) <- trans[, 1]
            trans <- trans[, 2]
            return(list(text = text, trans = trans, Ntree = Ntree))
        }
        else return(list(text = text, Ntree = Ntree))
    }
    else if (version > 1) {
        X <- scan(file = file, what = "", sep = "\n", quiet = TRUE, nmax=nmax)
        left <- grep("\\[", X)
        right <- grep("\\]", X)
        skip <- if (version <= 2) 
            grep("\\&map", X)
        else if (version > 2 && version <= 3) 
            grep("\\&prob", X)
        else if (version > 3) 
            grep("\\&label", X)
        left <- setdiff(left, skip)
        right <- setdiff(right, skip)
        if (length(left)) {
            w <- left == right
            if (any(w)) {
                s <- left[w]
                X[s] <- gsub("\\[[^]]*\\]", "", X[s])
            }
            w <- !w
            if (any(w)) {
                s <- left[w]
                X[s] <- gsub("\\[.*", "", X[s])
                sb <- right[w]
                X[sb] <- gsub(".*\\]", "", X[sb])
                if (any(s < sb - 1)) 
                  X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
            }
        }
        endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
        semico <- grep(";", X)
        i1 <- grep("begin trees;", X, ignore.case = TRUE)
        i2 <- grep("translate", X, ignore.case = TRUE)
        translation <- if (length(i2) == 1 && i2 > i1) 
            TRUE
        else FALSE
        if (translation) {
            end <- semico[semico > i2][1]
            x <- X[(i2 + 1):end]
            x <- unlist(strsplit(x, "[,;\t]"))
            x <- unlist(strsplit(x, " "))
            x <- x[nzchar(x)]
            trans <- matrix(x, ncol = 2, byrow = TRUE)
            trans[, 2] <- gsub("['\"]", "", trans[, 2])
            n <- dim(trans)[1]
        }
        start <- if (translation) 
            semico[semico > i2][1] + 1
        else semico[semico > i1][1]
        end <- endblock[endblock > i1][1] - 1
        tree <- X[start:end]
        tree <- sub(".* = *", "", tree)
        tree <- tree[tree != ""]
        semico <- grep(";", tree)
        Ntree <- length(semico)
        if (Ntree == 1 && length(tree) > 1) {
            STRING <- paste(tree, collapse = "")
        }
        else {
            if (any(diff(semico) != 1)) {
                STRING <- character(Ntree)
                s <- c(1, semico[-Ntree] + 1)
                j <- mapply(":", s, semico)
                if (is.list(j)) {
                  for (i in 1:Ntree) STRING[i] <- paste(tree[j[[i]]], 
                    collapse = "")
                }
                else {
                  for (i in 1:Ntree) STRING[i] <- paste(tree[j[, 
                    i]], collapse = "")
                }
            }
            else STRING <- tree
        }
        text <- STRING
        if (translation == TRUE) {
            rownames(trans) <- trans[, 1]
            trans <- trans[, 2]
            return(list(text = text, trans = trans, Ntree = Ntree))
        }
        else return(list(text = text, Ntree = Ntree))
    }
} # END phytools_readNexusData2 <- function(file, version, nmax=1000000) 


# phytools:::getEdgeLength
phytools_getEdgeLength2 <- function (text, start) 
	{
	i <- start + 1
	l <- 1
	temp <- vector()
	# 2026-09-07: If you have a root branch, there is no concluding , or )
	# So the search goes forever
	# But, this can cause issues for other tables of node/branch values
	# Better to cut the root branch.
	# Adding ";"
	# while (is.na(match(x=text[i], table=c(",", ")")))) {
	while (is.na(match(x=text[i], table=c(",", ")", ";"))))
		{
		if (text[i] != ";")
			{
			temp[l] <- text[i]
			l <- l + 1
			i <- i + 1
			}
		}
		list(edge.length = as.numeric(paste(temp, collapse = "")), end = i)
	} # phytools_getEdgeLength2 <- function (text, start) 


# phytools:::getLabel
phytools_getLabel2 <- function (text, start, stop.char = c(",", ":", ")", ";")) 
{
    i <- 0
    while (is.na(match(text[i + start], stop.char))) i <- i + 
        1
    label <- paste(text[0:(i - 1) + start], collapse = "")
    return(list(label = paste(label, collapse = ""), end = i + 
        start))
} # END phytools_getLabel2 <- function (text, start, stop.char = c(",", ":", ")", ";")) 

#phytools:::getBS
phytools_getBS2 <- function (text, start) 
{
    i <- start
    if (text[i] == "[") {
        j <- 1
        xx <- vector()
        while (text[i] != "]") {
            i <- i + 1
            xx[j] <- text[i]
            j <- j + 1
        }
        label <- paste(xx[1:(length(xx) - 1)], collapse = "")
        label <- sub("&label=", "", label)
    }
    list(label = label, end = i + 1)
} # END phytools_getBS2 <- function (text, start) 

#phytools:::modified.text_to_tree
phytools_modified.text_to_tree2 <- function(text, trans) 
{
    text <- unlist(strsplit(text, NULL))
    tip.label <- vector(mode = "character")
    edge <- matrix(c(1, NA), 1, 2)
    edge.length <- vector()
    node.label <- vector(mode = "character")
    currnode <- 1
    Nnode <- currnode
    i <- j <- k <- 1
    while (text[i] != "(") i <- i + 1
    while (text[i] != ";") {
        cat(i)
        cat(",")
        if (text[i] == "(") {
            cat("here1")
            if (j > nrow(edge)) 
                edge <- rbind(edge, c(NA, NA))
            edge[j, 1] <- currnode
            i <- i + 1
            # Pull out just the core tree elements
            if (is.na(match(text[i], c("(", ")", ",", ":", ";")))) {
                temp <- phytools_getLabel2(text, i)
                tip.label[k] <- temp$label
                i <- temp$end
                edge[j, 2] <- -k
                k <- k + 1
                if (text[i] == ":") {
                  temp <- phytools_getEdgeLength2(text, i)
                  edge.length[j] <- temp$edge.length
                  i <- temp$end
                }
            }
            else if (text[i] == "(") {
                Nnode <- Nnode + 1
                currnode <- Nnode
                edge[j, 2] <- currnode
            }
            j <- j + 1
        }
        else if (text[i] == ")") {
            cat("here2")
            i <- i + 1
            if (text[i] == "[") {
                temp <- phytools_getBS2(text, i)
                node.label[currnode] <- as.character(temp$label)
                i <- temp$end
            }
            if (text[i] == ":") {
                temp <- phytools_getEdgeLength2(text, i)
                ii <- match(currnode, edge[, 2])
                edge.length[ii] <- temp$edge.length
                i <- temp$end
            }
            currnode <- edge[match(currnode, edge[, 2]), 1]
        }
        else if (text[i] == ",") {
            cat("here3")
            if (j > nrow(edge)) 
                edge <- rbind(edge, c(NA, NA))
            edge[j, 1] <- currnode
            i <- i + 1
            if (is.na(match(text[i], c("(", ")", ",", ":", ";")))) {
                temp <- phytools_getLabel2(text, i)
                tip.label[k] <- temp$label
                i <- temp$end
                edge[j, 2] <- -k
                k <- k + 1
                if (text[i] == ":") {
                  temp <- phytools_getEdgeLength2(text, i)
                  edge.length[j] <- temp$edge.length
                  i <- temp$end
                }
            }
            else if (text[i] == "(") {
                Nnode <- Nnode + 1
                currnode <- Nnode
                edge[j, 2] <- currnode
            }
            j <- j + 1
        }
    }
    Ntip <- k - 1
    if (!is.null(trans)) 
        tip.label <- trans[tip.label]
    ntip <- abs(min(edge))
    edge[which(edge > 0)] <- ntip + edge[which(edge > 0)]
    edge[which(edge < 0)] <- abs(edge[which(edge < 0)])
    tree <- list(edge = edge, Nnode = as.integer(Nnode), tip.label = tip.label, 
        edge.length = edge.length, node.label = node.label)
    class(tree) <- "phylo"
    tree
} # END phytools_modified.text_to_tree2 <- function(text, trans)

