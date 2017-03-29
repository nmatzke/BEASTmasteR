#######################################################
# How to continue on errors during Rscript runs
#######################################################
continue_on_error <- function()
	{
	print("NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())'")
	}

testfunc <- function(a,b)
	{
	example='
	options(error=continue_on_error)

	print(1)
	print(2)

	# This should work
	testfunc(a=10, b=10)

	print(3)

	# This should produce an error and stop the Rscript run normally
	testfunc(a=1)

	print(4)
	print(5)
	'

	print(a+b)
	}



#######################################################
# Combine XML lists appropriately in XML
#######################################################
cl <- function(A, B, C=NULL, D=NULL, E=NULL, F=NULL, G=NULL, H=NULL, I=NULL, J=NULL, K=NULL, L=NULL, M=NULL, N=NULL, O=NULL, P=NULL, Q=NULL, R=NULL)
	{
	if (class(A)[1] != "list")
		{
		# Check if first item is a list; if not, and if xmlNode, make
		# part of a bigger list
		classA = class(A)
		if ("XMLNode" %in% classA)
			{
			A = c(list(xmlCommentNode(" blank tag to make a bigger list ")), list(A))
			} else {
			A = list(A)
			}
		} # END if (is.list(A) == FALSE)
	
	# Check if first item is a list; if not, and if xmlNode, make
	# part of a bigger list
	if (class(B)[1] != "list")
		{
		# Check the class; if it's an XMLNode, 
		# make it a part of a bigger node
		classB = class(B)
		if ("XMLNode" %in% classB)
			{
			B = c(list(xmlCommentNode(" blank tag to make a bigger list ")), list(B))
			} else {
			B = list(B)
			}
		} # END if (is.list(B) == FBLSE)


	# Additional...
	if ((is.null(C)==FALSE) && (class(C)[1] != "list"))
		{
		# Check the class; if it's an XMLNode, 
		# make it a part of a bigger node
		classC = class(C)
		if ("XMLNode" %in% classC)
			{
			C = c(list(xmlCommentNode(" blank tag to make a bigger list ")), list(C))
			} else {
			C = list(C)
			}
		} # END 
	if ((is.null(D)==FALSE) && (class(D)[1] != "list"))
		{
		# Check the class; if it's an XMLNode, 
		# make it a part of a bigger node
		classD = class(D)
		if ("XMLNode" %in% classD)
			{
			D = c(list(xmlCommentNode(" blank tag to make a bigger list ")), list(D))
			} else {
			D = list(D)
			}
		} # END
	if ((is.null(E)==FALSE) && (class(E)[1] != "list"))
		{
		# Check the class; if it's an XMLNode, 
		# make it a part of a bigger node
		classE = class(E)
		if ("XMLNode" %in% classE)
			{
			E = c(list(xmlCommentNode(" blank tag to make a bigger list ")), list(E))
			} else {
			E = list(E)
			}
		} # END 
	if ((is.null(F)==FALSE) && (class(F)[1] != "list"))
		{
		# Check the class; if it's an XMLNode, 
		# make it a part of a bigger node
		classF = class(F)
		if ("XMLNode" %in% classF)
			{
			F = c(list(xmlCommentNode(" blank tag to make a bigger list ")), list(F))
			} else {
			F = list(F)
			}
		} # END 
	if ((is.null(G)==FALSE) && (class(G)[1] != "list"))
		{
		# Check the class; if it's an XMLNode, 
		# make it a part of a bigger node
		classval = class(G)
		if ("XMLNode" %in% classval)
			{
			G = c(list(xmlCommentNode(" blank tag to make a bigger list ")), list(G))
			} else {
			G = list(G)
			}
		} # END 
	if ((is.null(H)==FALSE) && (class(H)[1] != "list"))
		{
		# Check the class; if it's an XMLNode, 
		# make it a part of a bigger node
		classval = class(H)
		if ("XMLNode" %in% classval)
			{
			H = c(list(xmlCommentNode(" blank tag to make a bigger list ")), list(H))
			} else {
			H = list(H)
			}
		} # END 
	if ((is.null(I)==FALSE) && (class(I)[1] != "list"))
		{
		# Check the class; if it's an XMLNode, 
		# make it a part of a bigger node
		classval = class(I)
		if ("XMLNode" %in% classval)
			{
			I = c(list(xmlCommentNode(" blank tag to make a bigger list ")), list(I))
			} else {
			I = list(I)
			}
		} # END 
	if ((is.null(J)==FALSE) && (class(J)[1] != "list"))
		{
		# Check the class; if it's an XMLNode, 
		# make it a part of a bigger node
		classval = class(J)
		if ("XMLNode" %in% classval)
			{
			J = c(list(xmlCommentNode(" blank tag to make a bigger list ")), list(J))
			} else {
			J = list(J)
			}
		} # END 
	if ((is.null(K)==FALSE) && (class(K)[1] != "list"))
		{
		# Check the class; if it's an XMLNode, 
		# make it a part of a bigger node
		classval = class(K)
		if ("XMLNode" %in% classval)
			{
			K = c(list(xmlCommentNode(" blank tag to make a bigger list ")), list(K))
			} else {
			K = list(K)
			}
		} # END 
	if ((is.null(L)==FALSE) && (class(L)[1] != "list"))
		{
		# Check the class; if it's an XMLNode, 
		# make it a part of a bigger node
		classval = class(L)
		if ("XMLNode" %in% classval)
			{
			L = c(list(xmlCommentNode(" blank tag to make a bigger list ")), list(L))
			} else {
			L = list(L)
			}
		} # END 
	if ((is.null(M)==FALSE) && (class(M)[1] != "list"))
		{
		# Check the class; if it's an XMLNode, 
		# make it a part of a bigger node
		classval = class(M)
		if ("XMLNode" %in% classval)
			{
			M = c(list(xmlCommentNode(" blank tag to make a bigger list ")), list(M))
			} else {
			M = list(M)
			}
		} # END 
	if ((is.null(N)==FALSE) && (class(N)[1] != "list"))
		{
		# Check the class; if it's an XMLNode, 
		# make it a part of a bigger node
		classval = class(N)
		if ("XMLNode" %in% classval)
			{
			N = c(list(xmlCommentNode(" blank tag to make a bigger list ")), list(N))
			} else {
			N = list(N)
			}
		} # END 
	if ((is.null(O)==FALSE) && (class(O)[1] != "list"))
		{
		# Check the class; if it's an XMLNode, 
		# make it a part of a bigger node
		classval = class(O)
		if ("XMLNode" %in% classval)
			{
			O = c(list(xmlCommentNode(" blank tag to make a bigger list ")), list(O))
			} else {
			O = list(O)
			}
		} # END 
	if ((is.null(P)==FALSE) && (class(P)[1] != "list"))
		{
		# Check the class; if it's an XMLNode, 
		# make it a part of a bigger node
		classval = class(P)
		if ("XMLNode" %in% classval)
			{
			P = c(list(xmlCommentNode(" blank tag to make a bigger list ")), list(P))
			} else {
			P = list(P)
			}
		} # END 
	if ((is.null(Q)==FALSE) && (class(Q)[1] != "list"))
		{
		# Check the class; if it's an XMLNode, 
		# make it a part of a bigger node
		classval = class(Q)
		if ("XMLNode" %in% classval)
			{
			Q = c(list(xmlCommentNode(" blank tag to make a bigger list ")), list(Q))
			} else {
			Q = list(Q)
			}
		} # END 
	if ((is.null(R)==FALSE) && (class(R)[1] != "list"))
		{
		# Check the class; if it's an XMLNode, 
		# make it a part of a bigger node
		classval = class(R)
		if ("XMLNode" %in% classval)
			{
			R = c(list(xmlCommentNode(" blank tag to make a bigger list ")), list(R))
			} else {
			R = list(R)
			}
		} # END 
				
	AB = c(A, B, C, D, E, F)
	return(AB)
	} # END cl <- function(A, B)


get_fn_prefix <- function(fn) 
	{
    defaults = "\n\tfn = \"example-data/aP6.fas\"\n\t"
    words = strsplit(fn, split = "\\.")[[1]]
    words_to_keep = words[1:(length(words) - 1)]
    prefix = paste(words_to_keep, collapse = ".", sep = "")
    return(prefix)
	}


get_fn_suffix <- function(fn) 
	{
    defaults = "\n\tfn = \"example-data/aP6.fas\"\n\t"
    words = strsplit(fn, split = "\\.")[[1]]
    words_to_keep = words[length(words)]
    suffix = paste(words_to_keep, collapse = ".", sep = "")
    return(suffix)
	}



# This should be a faster list2 string
list2str_fast <- function(list1, spacer="")
	{
	# convert to character
	# split into list of lists of characters
	# merge lists with unlist
	# paste with collapse argument
	tmpstr = paste(unlist(strsplit(as.character(list1), split="")), collapse=spacer)
	return(tmpstr)
	}




list2str_fast_nosplit <- function(list1, spacer="")
	{
	# convert to character
	# split into list of lists of characters
	# merge lists with unlist
	# paste with collapse argument
	tmpstr = paste(unlist(as.character(list1)), collapse=spacer)
	return(tmpstr)
	}


#######################################################
# List subsetting
#######################################################
listextract_bynum <- function(listnum, listvals)
	{
	return(listvals[[listnum]])
	} # END listextract_bynum <- function(listnum, listvals)

# Subset a list by a vector of indices (1-based)
list_subset_bynum <- function(listnums, listvals)
	{
	newlist = lapply(X=listnums, FUN=listextract_bynum, listvals=listvals)
	return(newlist)
	} # END list_subset_bynum <- function(listnums, listvals)

# Subset a list by a vector of TF (TRUE/FALSE)
list_subset_byTF <- function(TF, listvals)
	{
	listnums = (1:length(listvals))[TF]
	newlist = lapply(X=listnums, FUN=listextract_bynum, listvals=listvals)
	return(newlist)
	} # END list_subset_byTF <- function(TF, listvals)

# TRUE/FALSE to numbers
TFs_to_nums <- function(TF)
	{
	listnums = (1:length(TF))[TF]
	return(listnums)
	} # END TFs_to_nums <- function(TF)



# return matching TRUE/FALSE values
# list1 (.e.g. a big list) TRUE if it is found in list2 (e.g. a smaller list)

#######################################################
# match_list1_in_list2
#######################################################
#' Return TRUE for list1 items when they occur in list2
#' 
#' Return matching TRUE/FALSE values.  E.g. list1 (e.g. a big list) TRUE if it is found
#' in list2 (e.g. a smaller list)
#'
#' Utility function for %in%, when one's brain gets confused.
#' 
#' @param list1 The list of things you want to check
#' @param list2 The list of things you want to check against
#' @return \code{matchlist} The TRUE/FALSE list for list1
#' @export
#' @seealso \code{\link[base]{match}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
match_list1_in_list2 <- function(list1, list2)
	{
	matchlist = list1 %in% list2
	return(matchlist)
	}



#######################################################
# unlist_dtf_cols
#######################################################
#' Unlist the columns in a data.frame
#' 
#' Utility function. What it says.
#' 
#' @param dtf Input \code{\link[base]{data.frame}}
#' @param printflag Print the results if TRUE.
#' @return \code{dtf} The data.frame, hopefully without lists for columns
#' @export
#' @seealso \code{\link[base]{unlist}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
unlist_dtf_cols <- function(dtf, printflag=FALSE)
	{
	# Sometimes cbind makes each column a list, this can screw up use/searching of
	#  the column later on.  
	# Unlist each column...
	for (i in 1:ncol(dtf))
		{
		tmpstr = paste("unlisting col: ", names(dtf)[i], "...", sep="")
		prflag(tmpstr, printflag=printflag)		
		
		# catch a possible error from unlisting
		# the number of rows needs to stay the same!!
		tmpcol = unlist(dtf[, i])
		if (length(tmpcol) != length(dtf[, i]))
			{
			tmpstr = paste("...failed! unlist(col) length=", length(tmpcol), "; nrow(dtf) = ", nrow(dtf), sep="")
			prflag(tmpstr, printflag=printflag)
			} 
		else
			{
			dtf[, i] = tmpcol
			tmpstr = paste(" ", " ", sep="")
			prflag(tmpstr, printflag=printflag)
			}
		}
	
	#dtf2 = adf(dtf)
	
	return(dtf)
	}


# NOTE!!! THESE MATCH FUNCTIONS JUST RETURN THE *FIRST* MATCH, *NOT* ALL MATCHES
# (argh)
# return indices in 2nd list matching the first list
# It WILL return one match for each item in the list, though...

#######################################################
# get_indices_where_list1_occurs_in_list2
#######################################################
#' Return (first!) indices in second list matching the first list
#' 
#' This function will return one match (the first) for each item in the list; i.e. the second-list
#' index for each item in the first list.  Only the first hit in the second list is returned.
#' 
#' This is used by \code{\link{prt}}.
#'
#' @param list1 The first list. 
#' @param list2 The second list list.
#' @return \code{match_indices} The match indices.
#' @export
#' @seealso \code{\link{prt}}, \code{\link[base]{LETTERS}}, \code{\link{get_indices_where_list1_occurs_in_list2_noNA}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' list1 = c("N", "I", "C", "K")
#' list2 = LETTERS
#' get_indices_where_list1_occurs_in_list2(list1, list2)
get_indices_where_list1_occurs_in_list2 <- function(list1, list2)
	{
	match_indices = match(list1, list2)
	return(match_indices)
	}


# return indices in 2nd list matching the first list
#######################################################
# get_indices_where_list1_occurs_in_list2_noNA
#######################################################
#' Return (first!) indices in second list matching the first list, excluding NAs
#' 
#' This function will return one match (the first) for each item in the list; i.e. the second-list
#' index for each item in the first list.  Only the first hit in the second list is returned.  Unlike 
#' \code{\link{get_indices_where_list1_occurs_in_list2}}, non-hits (NAs) are excluded.
#' 
#' This is used by get_indices_of_branches_under_tips, which is used by \code{\link{extend_tips_to_ultrametricize}}, which can be used by section_the_tree.
#'
#' @param list1 The first list. 
#' @param list2 The second list list.
#' @return \code{match_indices} The match indices.
#' @export
#' @seealso \code{\link{prt}}, \code{\link[base]{LETTERS}}, \code{\link{get_indices_where_list1_occurs_in_list2}}, 
#' \code{\link{extend_tips_to_ultrametricize}}, \code{\link{section_the_tree}}, \code{\link{return_items_not_NA}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' list1 = c("N", "I", "C", "K")
#' list2 = LETTERS
#' get_indices_where_list1_occurs_in_list2_noNA(list1, list2)
#' 
get_indices_where_list1_occurs_in_list2_noNA <- function(list1, list2)
	{
	match_indices = match(list1, list2)
	match_indices = return_items_not_NA(match_indices)
	return(match_indices)
	}


#######################################################
# return_items_not_NA
#######################################################
#' Remove NAs from a vector/list
#' 
#' Utility function. This function returns the non-NA values from a vector.
#' 
#' This is used by \code{\link{get_indices_where_list1_occurs_in_list2_noNA}}, which is used 
#' by \code{\link{get_indices_of_branches_under_tips}}, which is used by 
#' \code{\link{extend_tips_to_ultrametricize}}, which can be used by \code{\link{section_the_tree}}.
#'
#' @param x The vector of items to check for being not NA.
#' @return \code{y} The surviving, non-NA cells of a vector.
#' @export
#' @seealso \code{\link{prt}}, \code{\link[base]{LETTERS}}, \code{\link{get_indices_where_list1_occurs_in_list2_noNA}},  
#' \code{\link{get_indices_where_list1_occurs_in_list2}}, \code{\link{extend_tips_to_ultrametricize}}, \code{\link{section_the_tree}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' list1 = c("N", "I", NA, "C", "K")
#' return_items_not_NA(list1)
#' 
return_items_not_NA <- function(x)
	{
	y = x[!is.na(x)]
	return(y)
	}


#######################################################
# order_tipranges_by_tr
#######################################################
#' Order the tipranges in a tipranges object so they match the order of tips in a tree
#' 
#' Utility function. What it says.  Life can get very confusing if you don't do this before plotting.
#' 
#' @param tipranges A tipranges object.
#' @param tr An ape tree object.
#' @return \code{tipranges} The reordered data.frame
#' @export
#' @seealso \code{\link[base]{unlist}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
order_tipranges_by_tr <- function(tipranges, tr)
	{
	tipranges_names = rownames(tipranges@df)
	tr_names = tr$tip.label
	
	match_indices = get_indices_where_list1_occurs_in_list2(list1=tr_names, list2=tipranges_names)
	
	tmpdf = tipranges@df[match_indices, ]
	tipranges@df = tmpdf
	
	return(tipranges)
	}



# Extract just the numbers from a string
# just the numbers INCLUDING THOSE CONNECTED BY DECIMAL POINTS!!!!
#######################################################
# extract_numbers
#######################################################
#' Extract just the numbers from a string, including decimal points
#' 
#' This function extracts numbers from a string.  Contiguous digits, including
#' decimal points, are made into a single number. A list of numbers is returned.
#' 
#' This saves you having to remember the \code{regexp}/\code{\link[base]{gregexpr}} code for this sort of thing, and
#' makes it much easier to parse numbers out of the text output of various programs.
#' 
#' @param tmpstr An input string.
#' @return \code{x2} The list of numbers
#' @export
#' @seealso \code{\link[base]{gregexpr}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' tmpstr = "190Ma - 65Ma"
#' extract_numbers(tmpstr)
#' 
#' tmpstr = "190.1Ma - 65.5Ma"
#' extract_numbers(tmpstr)
#' 
extract_numbers <- function(tmpstr)
	{
	defaults ='
	tmpstr = "190Ma - 65Ma"
	'
	# pull out the numbers / extract the numbers / extract numbers
	# just the numbers INCLUDING THOSE CONNECTED BY DECIMAL POINTS!!!!
	# (but not negative symbols)
	matches = gregexpr("(?:([0-9\\.]+))+", tmpstr)[[1]]
	
	# Get the ending points of the matches
	matches_end = matches-1+attr(matches,"match.length")
	
	# Extract the numbers from the string
	x = mapply(substr, tmpstr, matches, matches_end)

	# Convert to numeric
	x2 = as.numeric(x)
	return(x2)
	}


#######################################################
# list2str
#######################################################
#' Convert a list of items to a string
#' 
#' This is a shortcut to save time when converting a list of items to a string.
#' 
#' @param list1 The list to convert.
#' @param spacer The space between each item. Default " ".
#' @return \code{tmpstr} The output string.
#' @export
#' @seealso \code{\link[base]{paste}}, \code{\link[base]{as.character}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#'	 @cite FosterIdiots
#' @examples
#' test=1
#' 
list2str <- function(list1, spacer=" ")
	{
	
	for (i in 1:length(list1))
		{
		if (i == 1)
			{
			tmpstr = as.character(list1[1])
			if (length(list1) == 1)
				{
				return(tmpstr)
				}
			next
			}
		addstr = as.character(list1[i])
		tmpstr = paste(tmpstr, addstr, sep=spacer)
		}
	return(tmpstr)
	}





# Get the classes of the columns in a data frame
#######################################################
# cls.df
#######################################################
#' Get the class for each column in a list
#' 
#' This function returns the \code{\link[base]{class}} of each column in a \code{\link[base]{data.frame}}.
#' 
#' R does lots of weird and unpredictable things when you build up tables/matrices/data.frames
#' by e.g. \code{\link[base]{cbind}} and \code{\link[base]{rbind}} on vectors of results.  The major problems 
#' are (1) columns get made into class \code{\link[base]{list}}; (2) \code{\link[base]{numeric}}
#' columns are converted to class \code{\link[base]{factor}}; (3) \code{\link[base]{numeric}} columns
#' are converted to class \code{\link[base]{character}}; (4) you have a \code{\link[base]{matrix}} when
#' you think you have a \code{\link[base]{data.frame}}.
#' 
#' All of this could be taken care of by detailed understanding and tracking of when R recasts values in 
#' vectors, matrices, and data frames...but this is a huge pain, it is easier to just have a function
#' that jams everything back to a \code{\link[base]{data.frame}} with no lists, no factors, and with columns being numeric
#' where possible.  See \code{\link{dfnums_to_numeric}} and \code{\link{unlist_df4}} for these options.
#' 
#' @param dtf Input \code{\link[base]{data.frame}}.
#' @param printout Print the results to screen, if desired.
#' @return \code{dtf_classes} A \code{\link[base]{data.frame}} showing the column, column name, and column class.
#' @export
#' @seealso \code{\link{dfnums_to_numeric}}, \code{\link{unlist_df4}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
#' cls.df(x)
#' dfnums_to_numeric(adf(x))
#' unlist_df4(x)
#'
#' x = matrix(c(1,2,3,4,5,"A"), nrow=3, ncol=2)
#' cls.df(x)
#' dfnums_to_numeric(adf(x))
#' unlist_df4(x)
#'
#' x = adf(matrix(c(1,2,3,4,5,"A"), nrow=3, ncol=2))
#' names(x) = c("A","B")
#' cls.df(x)
#' dfnums_to_numeric(adf(x))
#' unlist_df4(x)
#' 
cls.df <- function(dtf, printout=FALSE)
	{
	# Convert to data.frame if needed
	if (class(dtf) == "matrix")
		{
		dtf = as.data.frame(dtf, stringsAsFactors=FALSE)
		}
	
	dtf_names = names(dtf)
	numcols = ncol(dtf)
	
	cls_col_list = c()
	for (i in 1:numcols)
		{
		# Initialize cls_col
		cls_col = NA
		
		# Get one column:
		cmdstr = paste("cls_col = class(dtf$'", dtf_names[i], "')", sep="")
		eval(parse(text = cmdstr))
		
		#cat(i, ": ", dtf_names[i], "	=	", cls_col, "\n", sep="")
		cls_col_list[i] = cls_col
		}
	
	# row names
	colnum = 1:numcols
	
	dtf_classes = cbind(colnum, dtf_names, cls_col_list)
	dtf_classes = data.frame(dtf_classes, row.names=colnum)
	
	# Print the output if true
	if (printout)
		{
		cat("\n")
		cat("cls.df(dtf) reports: dataframe 'dtf' has ", nrow(dtf), " rows, ", numcols, " columns.\n", sep="")
		cat("...names() and classes() of each column below...\n", sep="")
		cat("\n")
		print(dtf_classes)
		cat("\n")
		}	
	return(dtf_classes)
	}



# Get the classes of the columns in a data frame
#######################################################
# dfnums_to_numeric
#######################################################
#' Get the class for each column in a list
#' 
#' This function converts each column to class \code{\link[base]{numeric}} where possible, and
#' class \code{\link[base]{character}} otherwise.
#' 
#' R does lots of weird and unpredictable things when you build up tables/matrices/data.frames
#' by e.g. \code{\link[base]{cbind}} and \code{\link[base]{rbind}} on vectors of results.  The major problems 
#' are (1) columns get made into class \code{\link[base]{list}}; (2) \code{\link[base]{numeric}}
#' columns are converted to class \code{\link[base]{factor}}; (3) \code{\link[base]{numeric}} columns
#' are converted to class \code{\link[base]{character}}; (4) you have a \code{\link[base]{matrix}} when
#' you think you have a \code{\link[base]{data.frame}}.
#' 
#' All of this could be taken care of by detailed understanding and tracking of when R recasts values in 
#' vectors, matrices, and data frames...but this is a huge pain, it is easier to just have a function
#' that jams everything back to a \code{\link[base]{data.frame}} with no lists, no factors, and with columns being numeric
#' where possible.  See \code{\link{unlist_df4}} for more, and \code{\link{cls.df}} to see the class of each column.
#' 
#' \bold{WARNING: IF A COLUMN IS A MIX OF NUMBERS AND NON-NUMBERS, THE NON-NUMBERS WILL BE CONVERTED TO NA IF 
#' THE COLUMN IS MAJORITY NUMBERS (on default; see \code{max_NAs}).}
#' 
#' @param dtf Input \code{\link[base]{data.frame}}.
#' @param max_NAs Non-numeric cells will get converted to NA, up to the fraction of cells specified by \code{max_NAs}.  Above this
#' fraction, the column is converted to class \code{character}.
#' @param printout Print the results to screen, if desired.
#' @param roundval If not NULL, \code{\link[base]{round}} will be run using this for the number of digits.
#' @return \code{dtf} The output \code{\link[base]{data.frame}}.
#' @export
#' @seealso \code{\link{cls.df}}, \code{\link{unlist_df4}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
#' cls.df(x)
#' dfnums_to_numeric(adf(x))
#' unlist_df4(x)
#'
#' x = matrix(c(1,2,3,4,5,"A"), nrow=3, ncol=2)
#' cls.df(x)
#' dfnums_to_numeric(adf(x))
#' unlist_df4(x)
#'
#' x = adf(matrix(c(1,2,3,4,5,"A"), nrow=3, ncol=2))
#' names(x) = c("A","B")
#' cls.df(x)
#' dfnums_to_numeric(adf(x))
#' unlist_df4(x)
#' 
dfnums_to_numeric <- function(dtf, max_NAs=0.5, printout=FALSE, roundval=NULL)
	{
	dtf_classes = cls.df(dtf, printout=FALSE)
	
	dtf_names = names(dtf)
	numcols = ncol(dtf)
	
	cls_col_list = c()
	for (i in 1:numcols)
		{
		# Initialize cls_col
		cls_col = NA
		
		# Get one column:
		cmdstr = paste("cls_col = class(dtf$'", dtf_names[i], "')", sep="")
		eval(parse(text = cmdstr))
		
		#cat(i, ": ", dtf_names[i], "	=	", cls_col, "\n", sep="")
		cls_col_list[i] = cls_col
		}
	
	for (i in 1:numcols)
		{
		if (cls_col_list[i] == "list")
			{
			next()	# skip the lists
			}
		if (cls_col_list[i] != "numeric")
			{
			# Initialize cls_col
			newcol = NA
	
			# Get one column, convert to numeric:
			cmdstr = paste("newcol = as.numeric(as.character(dtf$'", dtf_names[i], "'))", sep="")
			# Evaluates to newcol
			suppressWarnings(eval(parse(text = cmdstr)))
			
			# If it's less than 50% NAs (or max_NA NAs), then convert to numeric
			#print(newcol)
			#print(max_NAs * length(newcol))
			#print(sum(is.na(newcol)))
			if (sum(is.na(newcol)) < (max_NAs * length(newcol)))
				{
				# Get the column, convert to numeric:
				# (if it's a factor, you have to convert to character, then to numeric)
				# (if it's a character, you can convert to character anyway...)
				#cmdstr = paste("dtf$", dtf_names[i], " = as.numeric(as.character(dtf$", dtf_names[i], "))", sep="")
				cmdstr = paste("dtf$'", dtf_names[i], "' = newcol", sep="")
				suppressWarnings(eval(parse(text = cmdstr)))
				
				
				# If a rounding val is specified, do the rounding
				if (!is.null(roundval))
					{
					cmdstr = paste("dtf$'", dtf_names[i], "' = round(dtf$'", dtf_names[i], "', digits=roundval)", sep="")
					suppressWarnings(eval(parse(text = cmdstr)))
					}
				
				}
			}
		}
	tmp_classes = cls.df(dtf)
	dtf_classes$newclasses = tmp_classes[,ncol(tmp_classes)]
	
	if(printout)
		{
		cat("\n")
		cat("dfnums_to_numeric(dtf, max_NAs=", max_NAs, ") reports: dataframe 'dtf_classes' has ", nrow(dtf_classes), " rows, ", ncol(dtf_classes), " columns.\n", sep="")
		cat("...names() and classes() of each column below...\n", sep="")
		cat("\n")
		print(dtf_classes)
		}
	return(dtf)
	}



#######################################################
# A few generic functions from genericR_v1.R, used in BioGeoBEARS
#######################################################
#' Convert to data.frame, without factors
#' 
#' Shortcut for: \code{as.data.frame(x, row.names=NULL, stringsAsFactors=FALSE)}
#' 
#' This function, and \code{\link{adf2}}, are useful for dealing with errors due to 
#' automatic conversion of some columns to factors.  Another solution may be to prepend
#' \code{options(stringsAsFactors = FALSE)} at the start of one's script, to turn off all default stringsAsFactors silliness.
#' 
#' @param x matrix or other object transformable to data.frame
#' @return data.frame
#' @export
#' @seealso \code{\link{adf2}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
#' adf(x)
#' 
adf <- function(x)
	{
	return(as.data.frame(x, row.names=NULL, stringsAsFactors=FALSE))
	}



#' Convert to data.frame, without factors
#' 
#' Shortcut for: \code{tmp_rownames = 1:nrow(x); as.data.frame(x, row.names=tmp_rownames, stringsAsFactors=FALSE)}
#' 
#' This function, and \code{\link{adf2}}, are useful for dealing with errors due to 
#' automatic conversion of some columns to factors.  Another solution may be to prepend
#' \code{options(stringsAsFactors = FALSE)} at the start of one's script, to turn off all default stringsAsFactors silliness.
#'
#' In adf2, rownames are forced to be numbers; this can prevent errors due to e.g. repeated rownames
#' after an \code{rbind} operation.
#'
#' @param x matrix or other object transformable to data.frame
#' @return data.frame
#' @export
#' @seealso \code{\link{adf}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
#' adf2(x)
#' 
adf2 <- function(x)
	{
	# Deals with the problem of repeated row names
	rownames = 1:nrow(x)
	return(as.data.frame(x, row.names=rownames, stringsAsFactors=FALSE))
	}


#' String splitting shortcut
#' 
#' \code{\link[base]{strsplit}} returns the results inside a list, which is annoying. \code{strsplit2} shortens the process.
#'
#' @param x A string to split
#' @param ... Other arguments to \code{\link[base]{strsplit}}.  The argument \code{split} is \emph{required}.
#' @return \code{out} The output from inside the list.
#' @export
#' @seealso \code{\link[base]{strsplit}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
#' 
#' # strsplit returns the results inside a list element
#' out = strsplit("ABC", split="")
#' out
#' # I.e....
#' out[[1]]
#' 
#' # If this is annoying/ugly in the code, use strsplit2:
#' out = strsplit2("ABC", split="")
#' out
#' 
strsplit2 <- function(x, ...)
	{
	out = strsplit(x, ...)[[1]]
	return(out)
	}


#######################################################
# unlist_df:
#######################################################
#' Unlist the columns in a data.frame
#' 
#' Sometimes, matrices or data.frames will malfunction due to their having lists as columns
#' and other weirdness.  This is a shortcut for \code{data.frame(lapply(df, function(x) unlist(x)))}.
#' 
#' @param df matrix or other object transformable to data.frame
#' @return data.frame
#' @export
#' @seealso \code{\link{unlist_df2}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
#' unlist_df2(x)
#'
unlist_df <- function(df)
	{
	outdf <- data.frame(lapply(df, function(x) unlist(x)))
	}



#######################################################
# unlist_df2:
#######################################################
#' Unlist the columns in a data.frame, with more checks
#' 
#' Sometimes, matrices or data.frames will malfunction due to their having lists as columns
#' and other weirdness. This runs \code{\link{unlist}} and additional checks.
#' 
#' @param df matrix or other object transformable to data.frame
#' @return \code{outdf} A \code{\link[base]{matrix}}.
#' @export
#' @seealso \code{\link{unlist_df}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
#' unlist_df2(x)
#'
unlist_df2 <- function(df)
	{
	store_colnames = names(df)
	
	outdf = NULL
	
	numrows = dim(df)[1]
	
	for (i in 1:ncol(df))
		{
		#print(names(df)[i])
		tmpcol = unlist(df[, i])
		#print(length(tmpcol))
		
		# Error check; e.g. blank cells might screw it up
		if (length(tmpcol) < numrows)
			{
			tmpcol2 = df[,i]
			tmpcol = as.character(tmpcol2)
			}
		
		outdf = cbind(outdf, tmpcol)
		}
	
	#outdf = adf2(outdf)
	
	names(outdf) = store_colnames
	return(outdf)
	}



#######################################################
# unlist_df3:
#######################################################
#' Unlist the columns in a data.frame, with more checks and adf
#' 
#' Sometimes, matrices or data.frames will malfunction due to their having lists as columns
#' and other weirdness. This runs \code{\link[base]{unlist}} and additional checks, and
#' forces conversion to a \code{\link[base]{data.frame}} at the end.
#' 
#' @param df matrix or other object transformable to data.frame
#' @return \code{outdf} data.frame
#' @export
#' @seealso \code{\link{unlist_df}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
#' unlist_df3(x)
#'
unlist_df3 <- function(df)
	{
	store_colnames = names(df)
	store_rownames = rownames(df)
	
	outdf = NULL
	
	numrows = dim(df)[1]
	
	for (i in 1:ncol(df))
		{
		#print(names(df)[i])
		tmpcol = unlist(df[, i])
		#print(length(tmpcol))
		
		# Error check; e.g. blank cells might screw it up
		if (length(tmpcol) < numrows)
			{
			tmpcol2 = df[,i]
			tmpcol = as.character(tmpcol2)
			}
		
		outdf = cbind(outdf, tmpcol)
		}
	
	outdf = adf2(outdf)
	
	names(outdf) = store_colnames
	rownames(outdf) = store_rownames
	return(outdf)
	}


#######################################################
# unlist_df4:
#######################################################
#' Unlist the columns in a data.frame, with more checks, adf, and dfnums_to_numeric
#' 
#' Sometimes, matrices or data.frames will malfunction due to their having lists as columns
#' and other weirdness. This runs \code{\link[base]{unlist}} and additional checks, and
#' forces conversion to a \code{\link[base]{data.frame}} at the end.  It also adds
#' \code{\link{dfnums_to_numeric}} which should remove the problem of numbers columns being of
#' class \code{\link[base]{character}}.
#' 
#' See especially  \code{\link[base]{data.matrix}} for a possibly simpler alternative.
#' 
#' @param df matrix or other object transformable to data.frame
#' @param ... Additional options passed to \code{\link{dfnums_to_numeric}}.
#' @return \code{outdf} data.frame
#' @export
#' @seealso \code{\link{unlist_df}}, \code{\link{dfnums_to_numeric}}, \code{\link{cls.df}}, \code{\link[base]{data.matrix}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
#' cls.df(x)
#' unlist_df4(x)
#'
#' x = matrix(c(1,2,3,4,5,"A"), nrow=3, ncol=2)
#' cls.df(x)
#' unlist_df4(x)
#'
#' x = adf(matrix(c(1,2,3,4,5,"A"), nrow=3, ncol=2))
#' names(x) = c("A","B")
#' cls.df(x)
#' unlist_df4(x)
#' 
unlist_df4 <- function(df, ...)
	{
	store_colnames = names(df)
	store_rownames = rownames(df)
	
	outdf = NULL
	
	numrows = dim(df)[1]
	numcols = dim(df)[2]
	
	for (i in 1:ncol(df))
		{
		#print(names(df)[i])
		tmpcol = unlist(df[, i])
		#print(length(tmpcol))
		
		# Error check; e.g. blank cells might screw it up
		if (length(tmpcol) < numrows)
			{
			tmpcol2 = df[,i]
			tmpcol = as.character(unlist(tmpcol2))
			}
		
		outdf = cbind(outdf, tmpcol)
		}

	# Unlist each row
# 	outdf2 = NULL
# 	for (i in 1:nrow(df))
# 		{
# 		#print(names(df)[i])
# 		tmprow = unlist(df[i, ])
# 		#print(length(tmpcol))
# 		
# 		# Error check; e.g. blank cells might screw it up
# 		if (length(tmprow) < numcols)
# 			{
# 			tmprow2 = df[i,]
# 			tmprow = as.character(unlist(tmprow2))
# 			}
# 		
# 		outdf2 = rbind(outdf2, tmprow)
# 		}
# 	outdf = outdf2
	
	outdf_tmp = adf2(outdf)
	
	# Remove factors and character silliness from numbers
	outdf = dfnums_to_numeric(outdf_tmp, ...)
	
	names(outdf) = store_colnames
	rownames(outdf) = store_rownames
	return(outdf)
	}




#######################################################
# slashslash:
#######################################################
#' Remove double slash (slash a slash)
#' 
#' Shortcut for: \code{gsub(pattern="//", replacement="/", x=tmpstr)}
#' 
#' This function is useful for removing double slashes that can
#' appear in full pathnames due to inconsistencies in trailing
#' slashes in working directories etc.
#'
#' @param tmpstr a path that you want to remove double slashes from
#' @return outstr a string of the fixed path
#' @export
#' @seealso \code{\link[base]{getwd}}, \code{\link[base]{setwd}}, \code{\link[base]{gsub}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' tmpstr = "/Library/Frameworks//R.framework/Versions/"
#'
#' outstr = slashslash(tmpstr)
#' outstr
#'
slashslash <- function(tmpstr)
	{
	outstr = gsub(pattern="//", replacement="/", x=tmpstr)
	
	# Check for single starting slash, and remove that
	chars = strsplit(outstr, split="")[[1]]
	chars
	
	TF = chars == "/"
	if ( (sum(TF) == 1) && (chars[1] == "/") )
		{
		chars2 = chars[2:length(chars)]
		outstr = paste(chars2, sep="", collapse="")
		}

	outstr = np(outstr)
	
	return(outstr)
	}


# Add a slash to a directory name if needed
#######################################################
# addslash:
#######################################################
#' Add a slash to a directory name if needed
#' 
#' This function adds a slash to the end of the string, if one is not present. Handy for standardizing paths.
#'
#' @param tmpstr a path that you want to possibly add a slash to 
#' @return outstr a string of the fixed path
#' @export
#' @seealso \code{\link[base]{getwd}}, \code{\link[base]{setwd}}, \code{\link[base]{gsub}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' tmpstr = "/Dropbox/_njm/__packages"
#' tmpstr
#' outstr = addslash(tmpstr)
#' outstr
#'
#' # Annoying, getwd() often doesn't return the ending slash, which 
#' # can make life hard for paste() later on
#' tmpstr = getwd()
#' tmpstr
#' outstr = addslash(tmpstr)
#' outstr
#'
addslash <- function(tmpstr)
	{
	tmpchars = strsplit(tmpstr, split="")[[1]]
	last_char = tmpchars[length(tmpchars)]
	
	if (nchar(tmpstr) == 0)
		{
		return("/")
		}
	
	if (last_char != "/")
		{
		outstr = paste(tmpstr, "/", sep="")
		} else {
		outstr = tmpstr
		}
		
	return(outstr)
	}



# source all .R files in a directory, except "compile" and "package" files

#######################################################
# sourceall
####################################ex###################
#' Source all .R files in a directory, except "compile" and "package" files
#' 
#' Utility function.
#' 
#' @param path The path to source
#' @param pattern Default is .R
#' @param ... Additional arguments to source
#' @return \code{path} The path that was sourced.
#' @export
#' @seealso \code{\link[base]{source}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
sourceall <- function(path=path, pattern="\\.R", ...)
	{
	tmppath = np(addslash(path))
	Rfiles = list.files(path=tmppath, pattern="\\.R", ...)
	
	# Files to remove
	Rfiles_remove_TF1 = grepl("compile", Rfiles)
	Rfiles_remove_TF2 = grepl("package", Rfiles)
	Rfiles_remove_TF = (Rfiles_remove_TF1 + Rfiles_remove_TF2) >= 1
	
	Rfiles = Rfiles[Rfiles_remove_TF == FALSE]

	cat("\nSourcing Rfiles in ", path, "...\n", sep="")

	
	for (Rfile in Rfiles)
		{
		cat("Sourcing Rfile: ", Rfile, "\n", sep="")
		fullfn = np(slashslash(paste(addslash(path), Rfile, sep="")))
		source(fullfn, chdir=TRUE, ...)
		}

	cat("\nDone sourcing Rfiles in ", path, "...\n", sep="")
	return(path)
	}




#######################################################
# printall
#######################################################
#' Print an entire table to screen
#' 
#' Utility function.  This prints a table to screen in chunks of \code{chunksize_toprint} 
#' (default=40).  This avoids the annoying situation of not being able to see the bottom 
#' of a table. Note that if you print something huge, you will be waiting for awhile (try
#' ESC or CTRL-C to cancel such an operation).
#'
#' Another option is to reset options to something like: \code{options(max.print=99999)}, but this
#' is hard to remember.  Your current setting is \code{getOption("max.print")}.
#' 
#' @param dtf The \code{\link[base]{data.frame}} to \code{\link[base]{print}}.
#' @param chunksize_toprint Number of lines to print. Default 50.
#' @param printflag For optional printing. Passed to \code{\link{prflag}}.
#' @return NULL
#' @export
#' @seealso \code{\link[base]{print}}, \code{\link{prflag}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
printall <- function(dtf, chunksize_toprint = 40, printflag=TRUE)
	{
	# Print everything in a data frame, in chunks of e.g. 50 rows
	if (nrow(dtf) <= chunksize_toprint)
		{
		prflag(dtf, printflag=printflag)
		return(dtf)
		}
	rows_toprint = seq(1, nrow(dtf), chunksize_toprint)
	
	if (printflag == TRUE)
		{
		for (i in 1 : (length(rows_toprint)-1) )
			{
			tmp_rows_toprint_start = rows_toprint[i]
			tmp_rows_toprint_end = rows_toprint[i+1]
			prflag(dtf[tmp_rows_toprint_start:tmp_rows_toprint_end, ])
			}
		
		# Then print the end
		tmp_rows_toprint_start = rows_toprint[length(rows_toprint)]
		tmp_rows_toprint_end = nrow(dtf)
		prflag(dtf[tmp_rows_toprint_start:tmp_rows_toprint_end, ])
		}	
	}





#######################################################
# prflag
#######################################################
#' Utility function to conditionally print intermediate results
#'
#' Just a handy shortcut function, allowing other functions to optionally 
#' print, depending on the value of \code{printflag}.
#' 
#' @param x What to print.
#' @param printflag If TRUE, do the printing
#' @return nothing
#' @export
#' @seealso \code{\link{get_daughters}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
prflag <- function(x, printflag=TRUE)
	{
	# A standard function to print (or not) certain variables,
	#   based on a master printflag
	# This avoids having to comment in/out various code chunks
	#   while debugging.
	if (printflag == TRUE)
		{
		# CAT instead of PRINT if it's a string or numeric
		if (is.character(x))
			{
			cat(x, "\n", sep="")
			}
		if (is.numeric(x))
			{
			cat(x, "\n", sep="")
			} else {
			print(x)
			}
		}
	else
		{
		pass="BLAH"
		}
	}




#######################################################
# np
#######################################################
#' normalizePath shortcut
#' 
#' Utility function that runs \code{\link[base]{normalizePath}}. Useful for
#' running on Mac vs. Windows.
#' 
#' @param path The path to run \code{\link[base]{normalizePath}} on.
#' @param ... Additional arguments to \code{\link[base]{normalizePath}}.
#' @return \code{path} The path that was normalized.
#' @export
#' @seealso \code{\link[base]{normalizePath}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' # Get a path
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' extdata_dir
#'
#' path = paste(extdata_dir, "//", "Psychotria_5.2.newick", sep="")
#' path
#'
#' path = np(path)
#' path
#' 
np <- function(path=path, ...)
	{
	path = suppressWarnings(normalizePath(path, ...))
	return(path)
	}





#######################################################
# strsplit_whitespace
#######################################################
#' Split strings on whitespace
#' 
#' This function splits strings on whitespace (spaces and tabs), so you don't have
#' to remember the \code{regexp}/\code{grep} format codes.
#' 
#' @param tmpline A string containing text.
#' @return \code{list_of_strs} 
#' @export
#' @seealso \code{\link[base]{strsplit}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' tmpline = "Hello world see	my	tabs."
#' strsplit_whitespace(tmpline)
#' 
strsplit_whitespace <- function(tmpline)
	{
	# split on 1 or more whitespaces
	temp = strsplit(tmpline, "[ \t]+")
	
	# get the list
	list_of_strs = temp[[1]]
	
	# remove any leading/trailing ""
	list_of_strs = list_of_strs[list_of_strs != ""]
	
	return(list_of_strs)
	}


#######################################################
# moref
#######################################################
#' print to screen the header of a file
#' 
#' This does the rough equivalent of the \code{UNIX} function \code{more}, but within R.
#' 
#' @param fn A filename.
#' @param printnotcat If \code{TRUE}, use \code{\link[base]{print}} instead of \code{\link[base]{cat}}. Default \code{FALSE}.
#' @return Nothing returned.
#' @export
#' @seealso \code{\link[base]{scan}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
moref <- function(fn, printnotcat = FALSE)
	{
	lines = scan(file=fn, what="character", sep="\n")
	
	if (printnotcat == TRUE)
		{
		for (i in 1:length(lines))
			{
			print(lines[i])
			}
		}
	else
		{
		for (i in 1:length(lines))
			{
			cat(paste(lines[i], "\n", sep=""))
			}
		}
	}









#######################################################
# TREE FUNCTIONS
#######################################################





#######################################################
# get_daughters
#######################################################
#' Get all the direct daughters nodes of a node
#' 
#' @param nodenum The node number to get the daughters of
#' @param t An ape phylo object
#' @return \code{daughter_nodenums} List of the daughter node numbers
#' @export
#' @seealso \code{\link{findall}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_daughters <- function(nodenum, t)
	{
	daughter_edgenums = findall(nodenum, t$edge[,1])
	daughter_nodenums = t$edge[,2][daughter_edgenums]
	return(daughter_nodenums)
	}




# Get indices of all matches to a list
#######################################################
# findall
#######################################################
#' Get indices of all matches to a list
#'
#' Just a handy shortcut function
#' 
#' @param what The item to find
#' @param inlist The list to search in 
#' @return \code{matching_indices} List of the matching indices
#' @export
#' @seealso \code{\link{get_daughters}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
findall <- function(what, inlist)
	{
	TFmatches = inlist == what
	indices = 1:length(inlist)
	matching_indices = indices[TFmatches]
	return(matching_indices)
	}





#######################################################
# get_parent
#######################################################
#' Get the direct parent node of a node
#' 
#' @param nodenum The node number to get the parent of
#' @param t An ape phylo object
#' @return \code{parent_nodenum}The parent node number
#' @export
#' @seealso \code{\link{findall}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_parent <- function(nodenum, t)
	{
	matching_edges = findall(nodenum, t$edge[,2])
	parent_nodenum = t$edge[,1][matching_edges][1]
	return(parent_nodenum)
	}


#######################################################
# get_level
#######################################################
#' Get a node's level in the tree
#'
#' Finds how many nodes deep a node is.
#' 
#' @param nodenum The node number to get the parent of
#' @param t An ape phylo object
#' @param tmplevel A starting level (the function is recursive)
#' @return \code{tmplevel} The level of the node.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_level <- function(nodenum, t, tmplevel=0)
	{
	parent_nodenum = get_parent(nodenum, t)
	if (is.na(parent_nodenum))
		{
		#tmplevel = 0
		return(tmplevel)
		}
	else
		{
		#print(paste("parent_nodenum: ", parent_nodenum, " level: ", tmplevel, sep=""))
		tmplevel = tmplevel + 1
		tmplevel = get_level(parent_nodenum, t, tmplevel)
		return(tmplevel)
		}
	# If an error occurs
	return(NA)
	}


#######################################################
# get_TF_tips
#######################################################
#' Get TRUE/FALSE for nodes being tips
#'
#' A utility function that returns \code{TRUE}/\code{FALSE} for whether or not each node is a tip.
#' 
#' @param obj An ape phylo object
#' @return \code{TF_tips} The \code{TRUE}/\code{FALSE} list for each tip.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}, \code{\link{match_list1_in_list2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_TF_tips <- function(obj)
	{
	# Get TF for nodes being tips
	
	# BIG CHANGE?
	#TF_tips = match_list1_in_list2(1:length(dists_from_root), obj$tip.label)
	TF_tips = match_list1_in_list2(1:length(obj$edge), 1:length(obj$tip.label))
	#TF_tips = obj$tip.label[TF_tips_indices]
	return(TF_tips)
	}



#######################################################
# get_TF_tips
#######################################################
#' Get TRUE/FALSE for nodes being tips
#'
#' A utility function that returns indices (node numbers) of the tips. This mostly saves typing.
#' 
#' @param obj An ape phylo object
#' @return \code{tip_indices} The node numbers of the tips.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}, \code{\link[ape]{phylo}}, \code{\link{get_indices_of_branches_under_tips}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_indices_of_tip_nodes <- function(obj)
	{
	tip_indices = 1:length(obj$tip.label)
	return(tip_indices)
	}

#######################################################
# get_indices_of_branches_under_tips
#######################################################
#' Get the indices of the branches (row number in edge matrix) below each tip
#'
#' A utility function. Gets the indices of the branches (row number in edge matrix) below each tip.
#' 
#' @param obj An \code{\link[ape]{ape}} \code{\link[ape]{phylo}} object
#' @return \code{branchnums_under_tips} The indices of the branches (row number in edge matrix) below each tip.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}, \code{\link{get_indices_of_tip_nodes}}, \code{\link{get_indices_where_list1_occurs_in_list2_noNA}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_indices_of_branches_under_tips <- function(obj)
	{
	tip_indices = get_indices_of_tip_nodes(obj)
	branchnums_under_tips = get_indices_where_list1_occurs_in_list2_noNA(tip_indices, obj$edge[, 2])
	return(branchnums_under_tips)
	}




#######################################################
# get_node_ages_of_tips
#######################################################
#' Get the ages of each tip above the root
#'
#' A utility function.
#' 
#' @param obj An ape phylo object
#' @return \code{TF_tips} The age (from the root) of each tip.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_node_ages_of_tips <- function(obj)
	{
	TF_tips = get_TF_tips(obj)
	root_node_num = get_nodenum_structural_root(obj)
	dists_from_root = dist.nodes(obj)[root_node_num, ]
	node_ages_of_tips = dists_from_root[TF_tips]
	return(node_ages_of_tips)
	}


#######################################################
# get_all_node_ages
#######################################################
#' Get the ages of all the nodes in the tree (above the root)
#'
#' A utility function. Use of \code{\link[ape]{dist.nodes}} may be slow.
#' 
#' @param obj An ape phylo object
#' @return \code{TF_tips} The age (from the root) of each node.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_all_node_ages <- function(obj)
	{
	node_ages = dist.nodes(obj)[get_nodenum_structural_root(obj), ]
	return(node_ages)
	}


#######################################################
# get_max_height_tree
#######################################################
#' Get the maximum age of all the nodes (above the root)
#'
#' I.e., the distance of the highest node above the root.  A utility function. 
#' Use of \code{\link[ape]{dist.nodes}} may be slow.
#' 
#' @param obj An ape phylo object
#' @return \code{max_height} The age (from the root) of the highest node.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_max_height_tree <- function(obj)
	{
	max_height = max(get_node_ages_of_tips(obj))
	return(max_height)
	}



#######################################################
# get_edge_times_before_present
#######################################################
#' Get the times of the top and bottom of each edge
#'
#' A utility function. 
#' 
#' @param t An ape phylo object
#' @return \code{edge_times_bp} A 2-column matrix with the age (from the present) of the top and bottom of each edge.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_edge_times_before_present <- function(t)
	{
	#height above root
	hts_at_end_of_branches_aka_at_nodes = t$edge.length
	hts_at_end_of_branches_aka_at_nodes = get_all_node_ages(t)
	h = hts_at_end_of_branches_aka_at_nodes

	# times before present, below (ultrametric!) tips
	# numbers are positive, i.e. in millions of years before present
	#                       i.e. mybp, Ma
	times_before_present = get_max_height_tree(t) - h

	
	# fill in the ages of each node for the edges
	edge_ages = t$edge
	edge_ages[,1] = h[t$edge[,1]]	# bottom of branch
	edge_ages[,2] = h[t$edge[,2]]	# top of branch

	# fill in the times before present of each node for the edges
	edge_times_bp = t$edge
	edge_times_bp[,1] = times_before_present[t$edge[,1]]	# bottom of branch
	edge_times_bp[,2] = times_before_present[t$edge[,2]]	# top of branch
	
	return(edge_times_bp)
	}








#######################################################
# extend_tips_to_ultrametricize
#######################################################
#' Take a tree, extend all tips (including fossils) up to 0.0 my before present
#' 
#' Makes tree precisely ultrametric by extending the terminal branches up to the highest tip (which is treated as 0 my before present).
#'
#' This function ADDS the time_before_present to everything, including fossils.  You have been warned.
#' 
#' @param obj An \code{\link[ape]{ape}} \code{\link[ape]{phylo}} object.
#' @param age_of_root The length of the branch below the root. Default 0.
#' @param tips_end_at_this_date The tips can be set to something other than 0, if desired.  (This could produce negative branchlengths, however.)
#' @return \code{obj} The corrected phylogeny
#' @export
#' @seealso \code{\link[ape]{read.tree}}, \code{\link{prt}}, \code{\link{average_tr_tips}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
extend_tips_to_ultrametricize <- function(obj, age_of_root=0, tips_end_at_this_date=NA)
	{
	#print("node ages of tips:")
	tip_ages = age_of_root + get_node_ages_of_tips(obj)
	#print(tip_ages)
	
	
	if (is.na(tips_end_at_this_date))
		{
		tips_end_at_this_date = max(tip_ages)
		}
	
	nums_to_add_to_tip_to_ultrametricize = tips_end_at_this_date - tip_ages
	
	indices_of_branches_under_tips = get_indices_of_branches_under_tips(obj)

	obj$edge.length[indices_of_branches_under_tips] = obj$edge.length[indices_of_branches_under_tips] + nums_to_add_to_tip_to_ultrametricize
	
	return(obj)
	}




#######################################################
# average_tr_tips
#######################################################
#' Average the heights of (non-fossil) tips to make ultrametric-ish.
#'
#' When you have a digitized tree, or other slightly uneven source tree, average the tips 
#' to get them all to line up at 0 my before present.  This makes an ultrametric tree
#' if and only if there are no fossil tips in the tree.
#'
#' If the user includes fossils accidentally, this function can easily lead to pathological
#' results (negative branch lengths etc.), so use with care!!
#' 
#' @param tr An ape phylo object
#' @param fossils_older_than Tips that are older than \code{fossils_older_than} will be excluded from the tips that 
#' are going to be averaged. This is not currently set to 0, because Newick files can have slight precision issues etc.
#' that mean not all tips quite come to zero (which is why you need \code{\link{average_tr_tips}} in the first place!).  
#' Obviously you should be cautious about the value of , depending on the absolute timescale of your tree. Make sure you do
#' not inappropriately average in fossils!!
#' @return \code{edge_times_bp} A 2-column matrix with the age (from the present) of the top and bottom of each edge.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}, \code{\link{extend_tips_to_ultrametricize}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
average_tr_tips <- function(tr, fossils_older_than=0.6)
	{
	#require(BioGeoBEARS)	# for prt()

	# Check for negative branchlengths
	blren_equal_below_0_TF = tr$edge.length <= 0
	if (sum(blren_equal_below_0_TF) > 0)
		{
		tmptxt = paste(tr$edge.length[blren_equal_below_0_TF], collapse=", ", sep="")
		stoptxt = paste("\nFATAL ERROR in average_tr_tips(): the INPUT tree has branchlengths <= 0:\n", tmptxt, 
		"\nThis can sometimes happen if you (A) are accidentally including fossil tips (change 'fossils_older_than'), or\n",
		"if your input tree was e.g. an MCC (majority clade consensus) trees output by BEAST's TreeAnnotator.\n", 
		"In that case, you must fix the input Newick file. See ?check_BioGeoBEARS_run for comments.\n", sep="")
		cat(stoptxt)
		stop(stoptxt)
		}



	tr_table = prt(tr, printflag=FALSE)

	tipnums_tmp = 1:length(tr$tip.label)
	
	# Cut out fossils!!
	tips_are_fossils_TF = tr_table$time_bp[tipnums_tmp] > fossils_older_than
	tipnums_tmp = tipnums_tmp[tips_are_fossils_TF == FALSE]
	
	edges_w_tips_TF = tr$edge[,2] %in% tipnums_tmp
	edges_rownums = (1:nrow(tr$edge))
	edges_w_tips_rownums = edges_rownums[edges_w_tips_TF]
	
	tipnums = tr$edge[edges_w_tips_rownums, 2]
	

	#tipnums = tr$edge[edges_w_tips_rownums, 2]
	#tipnums
	meanval = mean(tr_table$node_ht[tipnums])
	meanval
	diffs_from_mean = tr_table$node_ht[tipnums] - meanval
	diffs_from_mean

	tr5 = tr
	tr5$edge.length[edges_w_tips_rownums] = tr5$edge.length[edges_w_tips_rownums] - diffs_from_mean
	
	min(tr5$edge.length)
	min(tr$edge.length)


	# Check the output; if IT has negative branchlengths, return NA!!
	# Check for negative branchlengths
	blren_equal_below_0_TF = tr5$edge.length <= 0
	if (sum(blren_equal_below_0_TF) > 0)
		{
		tmptxt = paste(tr$edge.length[blren_equal_below_0_TF], collapse=", ", sep="")
		stoptxt = paste("\nFATAL ERROR in average_tr_tips(): the OUTPUT tree has branchlengths <= 0:\n", tmptxt, 
		"\nThis can sometimes happen if you (A) are accidentally including fossil tips (change 'fossils_older_than'), or\n",
		"(B) if average_tr_tips() introduced more negative branches (especially can happen with shallow branches).\n", 
		"Returning nada\n", sep="")
		cat(stoptxt)
		stop(stoptxt)
		return(NA)
		}
	
	#tr5_table = prt(tr5)
	#tr5_table
		
	return(tr5)
	}






#######################################################
# is.not.na
#######################################################
#' Check for not NA
#'
#' A utility function. 
#' 
#' @param x Thing to check for NA
#' @return \code{TRUE} or \code{FALSE}
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
is.not.na <- function(x)
	{
	return(is.na(x) == FALSE)
	}







# print tree in hierarchical format
#######################################################
# prt
#######################################################
#' Print tree in table format
#' 
#' Learning and using APE's tree structure can be difficult and confusing because much of the information is
#' implicit.  This function prints the entire
#' tree to a table, and makes much of the implicit information explicit.  It is not particularly fast, but
#' it is useful.
#'
#' See \url{http://ape.mpl.ird.fr/ape_development.html} for the official documentation of R tree objects.
#' 
#' @param t A \code{\link[ape]{phylo}} tree object.
#' @param printflag Should the table be printed to screen?  Default TRUE.
#' @param relabel_nodes Manually renumber the internal nodes, if desired. Default FALSE.
#' @param time_bp_digits The number of digits to print in the time_bp (time before present) column. Default=7.
#' @param add_root_edge Should a root edge be added?  Default \code{TRUE}.
#' @param get_tipnames Should the list of tipnames descending from each node be printed as a string in another column?  
#' This is slow-ish, but useful for matching up nodes between differing trees. Default \code{FALSE}.
#' @param fossils_older_than Tips that are older than \code{fossils_older_than} will be marked as \code{TRUE} in a column called \code{fossil}.
#' @param silence_warnings Suppress warnings about missing branchlengths (prt makes each branchlength equal 1)
#' This is not currently set to 0, because Newick files can have slight precision issues etc. that mean not all tips quite come to zero.  You 
#' can attempt to fix this with \code{\link{average_tr_tips}} (but make sure you do not inappropriately average in fossils!!).
#' @return \code{dtf} A \code{\link[base]{data.frame}} holding the table. (Similar to the printout of a \code{\link[phylobase]{phylo4}} object.)
#' @export
#' @seealso \code{\link[ape]{phylo}}, \code{\link{average_tr_tips}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://ape.mpl.ird.fr/ape_development.html}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
prt <- function(t, printflag=TRUE, relabel_nodes = FALSE, time_bp_digits=7, add_root_edge=TRUE, get_tipnames=FALSE, fossils_older_than=0.6, silence_warnings=FALSE)
	{
	defaults='
	t = tr
	printflag=TRUE;
	relabel_nodes = FALSE;
	time_bp_digits=7;
	add_root_edge=TRUE;
	get_tipnames=FALSE;
	fossils_older_than=0.6
	silence_warnings=FALSE
	'
	
	if (class(t) != "phylo")
		{
		txt = paste0("STOP ERROR IN prt(): the input 't' must be of class 'phylo'. You had class(t)='", class(t), "'.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if (class(t) != "phylo")
	

	# assemble beginning table
	
	# check if internal node labels exist
	if ("node.label" %in% attributes(t)$names == FALSE)
		{
		rootnum = get_nodenum_structural_root(t)
		
		new_node_labels = paste("inNode", rootnum:(rootnum+t$Nnode-1), sep="")
		t$node.label = new_node_labels
		}
	
	# or manually relabel the internal nodes, if desired
	if (relabel_nodes == TRUE)
		{
		rootnum = get_nodenum_structural_root(t)
		
		new_node_labels = paste("inNode", rootnum:(rootnum+t$Nnode-1), sep="")
		t$node.label = new_node_labels
		}
	
	labels = c(t$tip.label, t$node.label)
	ordered_nodenames = get_nodenums(t)
	#nodenums = 1:length(labels)
	node.types1 = rep("tip", length(t$tip.label))
	node.types2 = rep("internal", length(t$node.label))
	node.types2[1] = "root"
	node.types = c(node.types1, node.types2)
	

	# These are the index numbers of the edges below each node
	parent_branches = get_indices_where_list1_occurs_in_list2(ordered_nodenames, t$edge[,2])


	# Error check for branchlengths
	edgelength_NULL_warning = FALSE
	if (is.null(t$edge.length))
		{
		edgelength_NULL_warning = TRUE
		txt = paste0("\n\nWARNING: 'brlen_to_parent = t$edge.length[parent_branches]'...produced NULL as a result. This is probably a parsimony tree/cladogram with no branchlengths.\nprt() is inserting '1' for each branchlength. You may or may not want this. THESE ARE NOT REAL BRANCHLENGTHS AND NEITHER THE BRANCHLENGTHS NOR THE DERIVED 'TIMES' ETC SHOULD BE USED!!! You have been warned!\n\n")
		if (silence_warnings == FALSE)
			{
			cat(txt)
			}
		
		t$edge.length =rep(1, times=length(parent_branches))
		} # END if (is.null(brlen_to_parent))


	#parent_edges = parent_branches
	brlen_to_parent = t$edge.length[parent_branches]
	
	
	
	parent_nodes = t$edge[,1][parent_branches]
	daughter_nodes = lapply(ordered_nodenames, get_daughters, t)
	
	# print out the structural root, if desired
	root_nodenum = get_nodenum_structural_root(t)
	tmpstr = paste("prt(t): root=", root_nodenum, "\n", sep="")
	prflag(tmpstr, printflag=printflag)
	
	levels_for_nodes = unlist(lapply(ordered_nodenames, get_level, t))
	#tmplevel = get_level(23, t)
	#print(tmplevel)
	
	
	#height above root
	hts_at_end_of_branches_aka_at_nodes = t$edge.length
	hts_at_end_of_branches_aka_at_nodes = get_all_node_ages(t)
	h = hts_at_end_of_branches_aka_at_nodes

	# times before present, below (ultrametric!) tips
	# numbers are positive, i.e. in millions of years before present
	#                       i.e. mybp, Ma
	times_before_present = get_max_height_tree(t) - h

	
	# fill in the ages of each node for the edges
	edge_ages = t$edge
	edge_ages[,1] = h[t$edge[,1]]	# bottom of branch
	edge_ages[,2] = h[t$edge[,2]]	# top of branch


	# fill in the times before present of each node for the edges
	edge_times_bp = t$edge
	edge_times_bp[,1] = times_before_present[t$edge[,1]]	# bottom of branch
	edge_times_bp[,2] = times_before_present[t$edge[,2]]	# top of branch
	
	
	# If desired, get the list of all tipnames descended from a node, in alphabetical order
	if (get_tipnames == TRUE)
		{
		# Make the empty list
		list_of_clade_members_lists = rep(list(NA), length(ordered_nodenames))
		
		# Tips have only one descendant
		list_of_clade_members_lists[1:length(t$tip.label)] = t$tip.label
		list_of_clade_members_lists
		
		
		nontip_nodenums = (length(t$tip.label)+1) : length(ordered_nodenames)
		if (length(nontip_nodenums) > 1)
			{
			# More than 1 node
			nontip_nodenames = ordered_nodenames[nontip_nodenums]
			nontip_cladelists = sapply(X=nontip_nodenames, FUN=get_all_daughter_tips_of_a_node, t=t)
			nontip_cladelists
			
			nontip_cladelists_alphabetical = sapply(X=nontip_cladelists, FUN=sort)
			nontip_cladelists_alphabetical
			
			nontip_cladelists_alphabetical_str = sapply(X=nontip_cladelists_alphabetical, FUN=paste, collapse=",")
			nontip_cladelists_alphabetical_str
			
			# Store the results
			list_of_clade_members_lists[nontip_nodenums] = nontip_cladelists_alphabetical_str
			list_of_clade_members_lists
			} else {
			# Just one node
			nontip_nodenames = ordered_nodenames[nontip_nodenums]
			nontip_cladelists = sapply(X=nontip_nodenames, FUN=get_all_daughter_tips_of_a_node, t=t)
			nontip_cladewords = unlist(sapply(X=nontip_cladelists, FUN=strsplit, split=","))
			
			nontip_cladelists_alphabetical = sort(nontip_cladewords)
			nontip_cladelists_alphabetical
			
			nontip_cladelists_alphabetical_str = paste(nontip_cladelists_alphabetical, collapse=",", sep="")
			nontip_cladelists_alphabetical_str
			
			# Store the results
			list_of_clade_members_lists[nontip_nodenums] = nontip_cladelists_alphabetical_str
			list_of_clade_members_lists			
			}
			
		}

	
	# Add fossils TRUE/FALSE column.  You can turn this off with fossils_older_than=NULL.
	fossils = times_before_present > fossils_older_than

	# Obviously, internal nodes are irrelevant and should be NA
	tmpnodenums = (length(t$tip.label)+1) : ( length(t$tip.label) + t$Nnode )
	fossils[tmpnodenums] = NA
	
	if (get_tipnames == FALSE)
		{
		# Don't put in the list of clade names
		tmpdtf = cbind(1:length(ordered_nodenames), ordered_nodenames, levels_for_nodes, node.types, parent_branches, brlen_to_parent, parent_nodes, daughter_nodes, h, times_before_present, fossils, labels)
		
		dtf = as.data.frame(tmpdtf, row.names=NULL)
		# nd = node
		
		# edge.length is the same as brlen_2_parent
		names(dtf) = c("node", "ord_ndname", "node_lvl", "node.type", "parent_br", "edge.length", "ancestor", "daughter_nds", "node_ht", "time_bp", "fossils", "label")
		
		# convert the cols from class "list" to some natural class
		dtf = unlist_dtf_cols(dtf, printflag=FALSE)
		} else {
		# Put in the list of clade names
		tmpdtf = cbind(1:length(ordered_nodenames), ordered_nodenames, levels_for_nodes, node.types, parent_branches, brlen_to_parent, parent_nodes, daughter_nodes, h, round(times_before_present, digits=time_bp_digits), fossils, labels, list_of_clade_members_lists)
		
		dtf = as.data.frame(tmpdtf, row.names=NULL)
		# nd = node
		
		# edge.length is the same as brlen_2_parent
		names(dtf) = c("node", "ord_ndname", "node_lvl", "node.type", "parent_br", "edge.length", "ancestor", "daughter_nds", "node_ht", "time_bp", "fossils", "label", "tipnames")
		
		# convert the cols from class "list" to some natural class
		dtf = unlist_dtf_cols(dtf, printflag=FALSE)		
		}
	

	
	
	
	
	# Add the root edge, if desired
	# (AND, only if t$root.edge exists)
	if ( (add_root_edge == TRUE) && (!is.null(t$root.edge)) )
		{
		root_row_TF = dtf$node.type == "root"
		root_edge_length = t$root.edge
		
		# Stick in this edge length
		dtf$edge.length[root_row_TF] = root_edge_length
		
		# Add the root edge length to all node heights
		dtf$node_ht = dtf$node_ht + root_edge_length
		}
	
	# print if desired
	prflag(dtf, printflag=printflag)

	if ((edgelength_NULL_warning == TRUE) && (printflag == TRUE))
		{
		txt = paste0("\n\nWARNING: 'brlen_to_parent = t$edge.length[parent_branches]'...produced NULL as a result. This is probably a parsimony tree/cladogram with no branchlengths.\nprt() is inserting '1' for each branchlength. You may or may not want this. THESE ARE NOT REAL BRANCHLENGTHS AND NEITHER THE BRANCHLENGTHS NOR THE DERIVED 'TIMES' ETC SHOULD BE USED!!! You have been warned!\n\n")
		if (silence_warnings == FALSE)
			{
			cat(txt)
			}
		} # END if (edgelength_NULL_warning == TRUE)
	
	#tree_strings = c()
	#root_str = get_node_info(root_nodenum, t)
	return(dtf)
	} # END prt()






# Drop tip NAs

# If you remove 4 tips from a monophyletic group with drop.tip, you can get 2 NA tips
# This is annoying!
# To fix:
# 1. identify monophyletic NAs
# 2. drop those tips
# 3. repeat until there is just 1 NA

drop_tips_NA <- function(tmptr, copied_tiplabel_to_drop="NA")
	{
	TF = tmptr$tip.label == "NA"

	if (sum(TF) > 1)
		{
		tipnums_to_drop = (1:length(tmptr$tip.label))[TF]
	
		# Check for monophyly of the NA clade
		nodenum = getMRCA(phy=tmptr, tip=tipnums_to_drop)
		skeleton_subtree = extract.clade(phy=tmptr, node=nodenum)
		ntips = length(skeleton_subtree$tip.label)

		# If the clade isn't actually monophyletic in the sampled
		# skeleton tree, don't try to add it; skip this skeleton tree
		if (ntips != sum(TF))
			{
			stoptxt = paste0("STOP ERROR in drop_tips_NA(): in your input tree 'tmptr', the ", sum(TF), " tips labeled 'NA' did not form a monophyletic group. Monophyly of NA tips is required for this function.")
			cat("\n\n")
			cat(stoptxt)
			cat("\n\n")
			stop(stoptxt)
			} # END if (ntips != sum(TF))

		tmptr = drop.tip(phy=tmptr, tip=tipnums_to_drop, trim.internal=FALSE)
		tmptr = read.tree(file="", text=write.tree(phy=tmptr, file=""))
		} # END if (sum(TF) > 1)
	
	# Nest the function:
	TF = tmptr$tip.label == "NA"
	if (sum(TF) > 1)
		{
		tmptr = drop_tips_NA(tmptr=tmptr, copied_tiplabel_to_drop=copied_tiplabel_to_drop)
		tmptr = read.tree(file="", text=write.tree(phy=tmptr, file=""))
		} # END if (sum(TF) > 1)
	
	return(tmptr)
	} # END drop_tips_NA <- function(tmptr, copied_tiplabel_to_drop="NA")







manual_boxplot <- function(x, y, y_med, y0025, y0975, ymin, ymax, width=1, lwd=1, col="gray", border="black", extreme_width=0.1, logt=TRUE)
	{
	if (logt == TRUE)
		{
		# Absolute minimum
		y_all = c(y, y_med, y0025, y0975, ymin, ymax)
		y_pretty = pretty(y_all)
		y_observed_nonzero_min = min(y_all[y_all > 0])
		y_replacement_for_0 = y_observed_nonzero_min / 10
		y_pretty = c(y_observed_nonzero_min, y_pretty)
		log10_y_pretty = log10(y_pretty)
		
		y[y==0] = y_replacement_for_0
		y_med[y_med==0] = y_replacement_for_0
		y0025[y0025==0] = y_replacement_for_0
		y0975[y0975==0] = y_replacement_for_0
		ymin[ymin==0] = y_replacement_for_0
		ymax[ymax==0] = y_replacement_for_0
		
		y = log10(y)
		y_med = log10(y_med)
		y0025 = log10(y0025)
		y0975 = log10(y0975)
		ymin = log10(ymin)
		ymax = log10(ymax)
		} # END if (logt == TRUE)
	
	# Draw the box
	xl = x - width/2
	xr = x + width/2
	
	# Draw the little extreme line
	ew = extreme_width
	x_extr_l = x - ew
	x_extr_r = x + ew
	
	
	# Draw a rectangle
	rect(xleft=xl, ybottom=y0025, xright=xr, ytop=y0975, col=col, border=border, lwd=lwd)
	
	# Draw the mean line
	segments(x0=xl, x1=xr, y0=y, y1=y, col=border, lwd=lwd*2)
	
	# Median
	segments(x0=xl, x1=xr, y0=y, y1=y, col=border, lwd=lwd*2, lty="dashed")
	
	# Extremes
	segments(x0=x, x1=x, y0=y0975, y1=ymax, col=border, lwd=lwd, lty="solid")
	segments(x0=x_extr_l, x1=x_extr_r, y0=ymax, y1=ymax, col=border, lwd=lwd, lty="solid")
	segments(x0=x, x1=x, y0=y0025, y1=ymin, col=border, lwd=lwd, lty="solid")
	segments(x0=x_extr_l, x1=x_extr_r, y0=ymin, y1=ymin, col=border, lwd=lwd, lty="solid")
	
	return(NULL)
	}









