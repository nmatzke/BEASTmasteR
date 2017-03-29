

# https://gist.github.com/nmatzke/8f80723b6e1fc80ed5ac
#
#
#######################################################
# Calculating the number of unobservable site patterns 
# for the Mk-parsimony-informative ascertainment bias 
# correction
#
# by Nick Matzke, 2015-2016
#
# Free to re-use; we are preparing an ms on this issue,
# please email nick.matzke@anu.edu.au if interested.
# 
#######################################################




# Allman et al. 2010 as Mk pars-inf
# 
# Allman ES, Holder MT, Rhodes JA (2010).
# Estimating trees from filtered data: Identifiability of models 
# for morphological phylogenetics. Journal of Theoretical Biology 
# 263 (1):108-119


#######################################################
# Permutations
# http://stackoverflow.com/questions/7906332/how-to-calculate-combination-and-permutation-in-r
#######################################################
perm <- function(n, k)
	{
	return(factorial(n) / factorial(n-k))
	} # END perm <- function(n, k)

# Number of unobservable patterns in Mk-parsimony-informative 
# ascertainment bias correction (unordered characters)
num_unobservable_patterns_ParsInf <- function(ntaxa, nstates, ordering="unordered", printflag=TRUE)
	{
	defaults='
	ntaxa=5
	nstates=2
	printflag=TRUE
	'
	
	total_number_of_possible_patterns = nstates^ntaxa
	number_of_unobservable_patterns_given_nstates_observed = 0
	number_of_unobservable_patterns_given_j_states_in_pattern = 0
	number_of_unobservable_patterns_given_j_states_in_pattern_AND_patterns_are_uniform = 0
	
	# Initialize table
	pattern_counts_table = NULL
	
	
	
	# First, get the number of states in the patterns
		# NOTE: For an ordered character, the number of states to choose
		# for an unobservable, parsimony-uninformative character 
		# can only be 2, or 1 (anything with 3 observed states
		# will have some parsimony information, for an 
		# ordered character)
	max_nstates_in_unobservable_pattern = "none_error"
	if (ordering == "unordered")
		{
		max_nstates_in_unobservable_pattern = nstates
		}
	if (ordering == "ordered")
		{
		max_nstates_in_unobservable_pattern = min(2, nstates)
		}
	for (i in max_nstates_in_unobservable_pattern:1)
		{
		# Then, make a list of the number of combinations of states that are possible
		number_of_combinations_of_states = choose(n=nstates, k=i)
		
		# For each of these combinations, 1 of the states will be the 
		# dominant state, and (i-1) will be autapomorphic states
		# This can happen nstates i choose 1 ways
		number_of_ways_to_pick_the_dominant_state = choose(n=i, k=1)
		
		# After a dominant state is chosen, there are ntaxa choose (i-1)
		# ways to pick the autapomorphic taxa
		number_of_ways_to_pick_autapomorphic_taxa = choose(n=ntaxa, k=(i-1))
		
		# Given that the autapomorphic taxa have been picked, there are
		# (i-1) PERMUTE (i-1) ways to assign the autapomorphic states
		number_of_ways_to_assign_autapo_states_to_autapo_taxa = perm(n=(i-1), k=(i-1))
		
		product = number_of_combinations_of_states * number_of_ways_to_pick_the_dominant_state * number_of_ways_to_pick_autapomorphic_taxa * number_of_ways_to_assign_autapo_states_to_autapo_taxa
		
		tmprow = c(number_of_combinations_of_states, number_of_ways_to_pick_the_dominant_state, number_of_ways_to_pick_autapomorphic_taxa, number_of_ways_to_assign_autapo_states_to_autapo_taxa, product)
		
		pattern_counts_table = rbind(pattern_counts_table, tmprow)
		}
	
	
	abbreviations = c("ncomb_sts", "npicks_domstate", "npicks_autapo_taxa", "n_assign_autapo_states", "product")
	meaning = c("number of states in the pattern", "number_of_combinations_of_states", "number_of_ways_to_pick_the_dominant_state", "number_of_ways_to_pick_autapomorphic_taxa", "number_of_ways_to_assign_autapo_states_to_autapo_taxa", "number of unobservable patterns given this number of states in the pattern")
	headers_explained = cbind(c("(row names)", abbreviations), meaning)
	headers_explained = as.data.frame(headers_explained, stringsAsFactors=FALSE)
	names(headers_explained) = c("abbreviation", "meaning")
	row.names(headers_explained) = NULL
	
	
	# Label the table
	pattern_counts_table = as.data.frame(pattern_counts_table, stringsAsFactors=FALSE)
	names(pattern_counts_table) = abbreviations
	row.names(pattern_counts_table) = paste0("when_nstates=", max_nstates_in_unobservable_pattern:1)
	
	number_of_unobservable_patterns = sum(pattern_counts_table$product)
	
	percentage_of_patterns_unobservable = number_of_unobservable_patterns / total_number_of_possible_patterns * 100
	
	summary_table = matrix(data=c(ntaxa, nstates, number_of_unobservable_patterns, total_number_of_possible_patterns, percentage_of_patterns_unobservable), nrow=1)
	summary_table = as.data.frame(summary_table, stringsAsFactors=FALSE)
	names(summary_table) = c("ntaxa", "nstates", "number_of_unobservable_patterns", "total_number_of_possible_patterns", "percentage_of_patterns_unobservable")
	
	res = NULL
	res$headers_explained = headers_explained
	res$pattern_counts_table = pattern_counts_table
	res$summary_table = summary_table
	
	if (printflag == TRUE)
		{
		require(BioGeoBEARS)
		cat("\n\nTable counting the number of unobservable patterns:\n\n")
		print(conditional_format_table(pattern_counts_table))
		cat("\n\nExplanation of the table headers:\n\n")
		print(headers_explained)
		cat("\n\nSummary of possible and unobservable patterns:\n\n")
		print(conditional_format_table(summary_table))
		}
	
	
	return(res)
	} # END num_unobservable_patterns_ParsInf <- function(ntaxa, nstates, ordering="unordered", printflag=TRUE)





example_count_patterns = '
res = num_unobservable_patterns_ParsInf(ntaxa=2, nstates=2)
res = num_unobservable_patterns_ParsInf(ntaxa=3, nstates=3)
res = num_unobservable_patterns_ParsInf(ntaxa=4, nstates=4)
res = num_unobservable_patterns_ParsInf(ntaxa=5, nstates=5)


res = num_unobservable_patterns_ParsInf(ntaxa=10, nstates=2)
res = num_unobservable_patterns_ParsInf(ntaxa=10, nstates=3)
res = num_unobservable_patterns_ParsInf(ntaxa=10, nstates=4)
res = num_unobservable_patterns_ParsInf(ntaxa=10, nstates=5)


res = num_unobservable_patterns_ParsInf(ntaxa=100, nstates=2)
res = num_unobservable_patterns_ParsInf(ntaxa=100, nstates=3)
res = num_unobservable_patterns_ParsInf(ntaxa=100, nstates=4)
res = num_unobservable_patterns_ParsInf(ntaxa=100, nstates=5)




res

sum(pattern_counts_table$prod)

pattern_counts_table = num_unobservable_patterns_ParsInf(ntaxa=10, nstates=5)
sum(pattern_counts_table$prod)

pattern_counts_table = num_unobservable_patterns_ParsInf(ntaxa=100, nstates=5)
sum(pattern_counts_table$prod)

pattern_counts_table = num_unobservable_patterns_ParsInf(ntaxa=1000, nstates=5)
sum(pattern_counts_table$prod)

pattern_counts_table = num_unobservable_patterns_ParsInf(ntaxa=1000, nstates=2)
sum(pattern_counts_table$prod)

' # END example_count_patterns


# List unobservable site patterns
# (assumes unordered characters)
#' assumed_nstates MrBayes (as of this writing) appears to only actually have code for
#' MkParsInf for binary characters. What it is doing for non-binary characters (3 states, 4 states, etc.)
#' is unknown, but it may be only correcting for the patterns as if it were a binary character.
#' By setting assumed_nstates=2, only the patterns for a binary character will be generated -- i.e.,
#' nstates will be reset to assumed nstates. Default: assumed_nstates=nstates.
list_unobservable_patterns_ParsInf <- function(ntaxa, nstates, printflag=TRUE, max_num_patterns=2500, ordering="unordered", assumed_nstates=nstates)
	{
	defaults='
	ntaxa=10
	nstates=2
	printflag=TRUE
	max_num_patterns=2500
	newcols = list_unobservable_patterns_ParsInf(ntaxa=10, nstates=2, printflag=TRUE, max_num_patterns=2500)
	newcols
	
	newcols = list_unobservable_patterns_ParsInf(ntaxa=10, nstates=3, printflag=TRUE, max_num_patterns=2500)
	newcols
	' # END defaults
	
	print("Assumed number of states, for calculating ascertainment-bias correction:")
	print(assumed_nstates)
	
	# If settings.xlsx is blank, just use nstates
	if ( (isblank_TF(assumed_nstates)) || (is.na(assumed_nstates)) || (assumed_nstates=="NA")  )
		{
		assumed_nstates = nstates
		} # END if ( (is.na(assumed_nstates)) || (assumed_nstates=="NA")  )
	
	
	# Reset nstates if desired
	if (assumed_nstates != nstates)
		{
		txt = paste0("Warning from list_unobservable_patterns_ParsInf(): because nstates != assumed_nstates, nstates is being reset from '", nstates, "', to '", assumed_nstates, "'.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		warning(txt)
		nstates = assumed_nstates
		nstates = as.numeric(nstates)
		} # END if (assumed_nstates != nstates)
	

	# Error check on ordering
	if ( (ordering != "unordered") && (ordering != "ordered") )
		{
		txt = paste0("\n\nSTOP ERROR in list_unobservable_patterns_ParsInf(): 'ordering' must be either\n'unordered' (default) or 'ordered'. Instead, you have: ", ordering, "\n\n")
		cat(txt)
		stop(txt)
		}
	
	if ((ordering != "ordered") && (ordering != "unordered"))
		{
		#error_txt = "STOP ERROR in list_unobservable_patterns_ParsInf(): MkParsInf for ordered characters will be different than for unordered characters, and I have not implemented it yet.\n"
		error_txt = "STOP ERROR in list_unobservable_patterns_ParsInf: 'ordering' must be 'unordered' or 'ordered'.\n"
		cat("\n\n")
		cat(error_txt)
		cat("\n\n")
		stop(error_txt)
		} # END if ((ordering != "ordered") && (ordering != "unordered"))
	
	
	# Get the numbers
	total_number_of_possible_patterns = nstates^ntaxa
	result = num_unobservable_patterns_ParsInf(ntaxa=ntaxa, nstates=nstates, printflag=printflag)
	num_patterns = result$summary_table$number_of_unobservable_patterns
	
	if (printflag)
		{
		txt = paste0("list_unobservable_patterns_ParsInf() is creating ", num_patterns, " unobservable site patterns. This will slow down the likelihood calculations!\n\nUse the function 'num_unobservable_patterns_ParsInf()' in order to see how many unobservable columns you are adding to your data.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		}

	if (num_patterns > max_num_patterns)
		{
		txt = paste0("STOP ERROR in list_unobservable_patterns_ParsInf(): You are trying to create ", num_patterns, " unobservable patterns, but your 'max_num_patterns' is set to: ", max_num_patterns, ".  You can manually increase max_num_patterns, but note that creating huge numbers of patterns (millions/billions) will get very slow, and eventually crash, either R or BEAST.  Use the function num_unobservable_patterns_ParsInf() to see how many unobservable columns you are adding to your data for a given number of taxa and character states.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		}
	
	# New columns
	if (printflag)
		{
		cat("\n\n")
		cat(paste0("list_unobservable_patterns_ParsInf() is creating an empty matrix for ", num_patterns, " new unobservable columns"))
		cat("\n\n")
		}
	newcols = matrix(data=0, nrow=ntaxa, ncol=num_patterns)
	colnum = 0
	
	# First, get the number of states in the patterns
		# NOTE: For an ordered character, the number of states to choose
		# for an unobservable, parsimony-uninformative character 
		# can only be 2, or 1 (anything with 3 observed states
		# will have some parsimony information, for an 
		# ordered character)
		# Then, add those patterns
	states = seq(0, nstates-1)
	max_nstates_in_unobservable_pattern = "none_error"
	if (ordering == "unordered")
		{
		max_nstates_in_unobservable_pattern = nstates
		}
	if (ordering == "ordered")
		{
		max_nstates_in_unobservable_pattern = min(2, nstates)
		}
	for (i in 1:max_nstates_in_unobservable_pattern)
		{
		# Then, make a list of the number of combinations of states that are possible
		number_of_combinations_of_states = choose(n=nstates, k=i)
		
		combn_of_states = combn(x=states, m=i, simplify=FALSE)
		
		for (c1 in 1:length(combn_of_states))
			{
			states_in_this_pattern = combn_of_states[[c1]]
			
			# If it's invariant, just do that
			if (length(states_in_this_pattern) == 1)
				{
				colnum=colnum+1
				tmpcol = rep(states_in_this_pattern, ntaxa)
				newcols[,colnum] = tmpcol
				next()
				} # END if (length(states_in_this_pattern) == 1)
			
			
			# For each of these combinations, 1 of the states will be the 
			# dominant state, and (i-1) will be autapomorphic states
			# This can happen nstates i choose 1 ways
			number_of_ways_to_pick_the_dominant_state = choose(n=i, k=1)
			combn_of_dominant = combn(x=states_in_this_pattern, m=1, simplify=FALSE)
			
			for (c2 in 1:length(combn_of_dominant))
				{
				dominant_state = combn_of_dominant[[c2]]
				TF = states_in_this_pattern == dominant_state
				nondominant_states = states_in_this_pattern[TF == FALSE]
								
				# After a dominant state is chosen, there are ntaxa choose (i-1)
				# ways to pick the autapomorphic taxa
				number_of_ways_to_pick_autapomorphic_taxa = choose(n=ntaxa, k=(i-1))
				combn_of_autapo_taxa = combn(x=1:ntaxa, m=(i-1), simplify=FALSE)
				
				for (c3 in 1:length(combn_of_autapo_taxa))
					{
					indexes_of_autapomorphic_taxa = combn_of_autapo_taxa[[c3]]
					
					# Given that the autapomorphic taxa have been picked, there are
					# (i-1) PERMUTE (i-1) ways to assign the autapomorphic states
					number_of_ways_to_assign_autapo_states_to_autapo_taxa = perm(n=(i-1), k=(i-1))
					
					matrix_of_permutations = gtools::permutations(n=(i-1), r=(i-1), v=nondominant_states)
					
					#matrix_of_permutations = gtools::permutations(n=(i-1), r=0, v=nondominant_states)
					
					
					for (c4 in 1:nrow(matrix_of_permutations))
						{
						colnum=colnum+1
						tmpcol = rep(dominant_state, ntaxa)
						
						# Add autapomorphies
						tmpcol[indexes_of_autapomorphic_taxa] = matrix_of_permutations[c4,]
						
						newcols[,colnum] = tmpcol
						} # END for (c4 in 1:nrow(matrix_of_permutations))
					} # END for (c3 in 1:length(combn_of_autapo_taxa))
				} # END for (c2 in 1:length(combn_of_dominant))
			} # END for (c1 in 1:length(combn_of_states))
		} # END for (i in max_nstates_in_unobservable_pattern:1)
	
	return(newcols)
	} # END list_unobservable_patterns_ParsInf <- function(ntaxa, nstates, printflag=TRUE)






# Calculate # of unobservable patterns for a bunch of combinations
make_table_num_unobservable_patterns <- function(ntaxa_list=c(4,5,10,20,50,100,200,500,1000), nstates_list=c(2,3,4,5,6), ordering="unordered", savefiles=TRUE, printflag=FALSE)
	{
	defaults='

	source("/drives/GDrive/__github/BEASTmasteR/R/Mk_ParsimonyInformative_v1.R")
	ntaxa_list = c(4,5,10,20,50,100,200,500,1000)
	nstates_list = c(2,3,4,5,6)
	savefiles=TRUE
	printflag=TRUE
	ordering="unordered"

	
	res1 = make_table_num_unobservable_patterns(ntaxa_list=ntaxa_list, nstates_list=nstates_list, ordering=ordering, savefiles=savefiles, printflag=printflag)

	ordering="ordered"
	res2 = make_table_num_unobservable_patterns(ntaxa_list=ntaxa_list, nstates_list=nstates_list, ordering=ordering, savefiles=savefiles, printflag=printflag)
	
	names(res1)
	res1$num_patterns_unobservable_matrix_FMT
	res2$num_patterns_unobservable_matrix_FMT

	' # END defaults

	npatterns_matrix = matrix(data=NA, nrow=length(ntaxa_list), ncol=length(nstates_list))
	ttlpats_matrix = matrix(data=NA, nrow=length(ntaxa_list), ncol=length(nstates_list))

	for (i in 1:length(ntaxa_list))
		{
		for (j in 1:length(nstates_list))
			{
			ntaxa = ntaxa_list[i]
			nstates = nstates_list[j]
			res = num_unobservable_patterns_ParsInf(ntaxa=ntaxa, nstates=nstates, printflag=printflag, ordering=ordering)
			num_unobservable_patterns = res$summary_table$number_of_unobservable_patterns
			#num_unobservable_patterns = num_unobservable_patterns_biogeography(ntaxa=ntaxa, nstates=NULL, numareas=numareas, maxareas=NULL)
		
			ttl_num_patterns = res$summary_table$total_number_of_possible_patterns
		
			npatterns_matrix[i,j] = num_unobservable_patterns
			ttlpats_matrix[i,j] = ttl_num_patterns
			} # END for (j in 1:max_numareas)
		} # END for (i in 1:max_ntaxa)

	# Write to file
	npatterns_matrix = as.data.frame(npatterns_matrix, stringsAsFactors=FALSE)
	names(npatterns_matrix) = nstates_list
	row.names(npatterns_matrix) = ntaxa_list
	npatterns_matrix

	conditional_format_table(npatterns_matrix)
	
	if (savefiles == TRUE)
		{
		xlsfn = paste0("num_unobservable_patterns_MkParsInf_", ordering, "_matrix_v1.txt")
		write.table(x=npatterns_matrix, file=xlsfn, quote=FALSE, sep="\t")
		if (printflag) { cat("\nWrote file: ", xlsfn, "") }

		xlsfn = paste0("num_unobservable_patterns_MkParsInf_", ordering, "_matrix_FORMATTED_v1.txt")
		write.table(x=conditional_format_table(npatterns_matrix), file=xlsfn, quote=FALSE, sep="\t")
		if (printflag) { cat("\nWrote file: ", xlsfn, "") }
		} # END if (savefiles == TRUE)


	########################################
	# Total number of possible patterns
	########################################
	ttlpats_matrix = as.data.frame(ttlpats_matrix, stringsAsFactors=FALSE)
	names(ttlpats_matrix) = nstates_list
	row.names(ttlpats_matrix) = ntaxa_list
	ttlpats_matrix

	conditional_format_table(ttlpats_matrix)

	########################################
	# Write to file, if desired
	########################################
	if (savefiles == TRUE)
		{
		xlsfn = paste0("num_POSSIBLE_patterns_MkParsInf_", ordering, "_matrix_v1.txt")
		write.table(x=ttlpats_matrix, file=xlsfn, quote=FALSE, sep="\t")
		if (printflag) { cat("\nWrote file: ", xlsfn, "") }
		
		xlsfn = paste0("num_POSSIBLE_patterns_MkParsInf_", ordering, "_matrix__FORMATTED_v1.txt")
		write.table(x=conditional_format_table(ttlpats_matrix), file=xlsfn, quote=FALSE, sep="\t")
		if (printflag) { cat("\nWrote file: ", xlsfn, "") }
		} # END if (savefiles == TRUE)


	res = NULL
	res$num_patterns_unobservable_matrix = npatterns_matrix
	res$num_patterns_unobservable_matrix_FMT = conditional_format_table(npatterns_matrix)
	res$num_patterns_total_matrix = ttlpats_matrix
	res$num_patterns_total_matrix_FMT = conditional_format_table(ttlpats_matrix)
	
	extract='
	npatterns_matrix = res$num_patterns_unobservable_matrix
	num_patterns_unobservable_matrix = res$num_patterns_unobservable_matrix

	num_patterns_unobservable_matrix_FMT = res$num_patterns_unobservable_matrix_FMT
	npatterns_matrix_FMT = res$num_patterns_unobservable_matrix_FMT


	ttlpats_matrix = res$num_patterns_total_matrix = 
	ttlpats_matrix_FMT = res$num_patterns_unobservable_matrix_FMT
	num_patterns_unobservable_matrix_FMT = res$num_patterns_unobservable_matrix_FMT
	' # END extract
	
	return(res)	
	} # END make_table_num_unobservable_patterns







