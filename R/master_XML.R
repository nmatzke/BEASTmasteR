
#######################################################
# Set up the overall structure of the Beast2 file
#######################################################

# Blankline
bl <- function()
	{
	txt = " blank line "
	blankline = xmlCommentNode(txt)
	return(blankline)
	}

print_master_XML_section_headers <- function(xml)
	{
	for (i in 1:length(xml))
		{
		# cat just the text contained in the 3rd XMLCommentNode
		cat("\n", xml[[i]][[3]]$value)
		}
	cat("\n\n")
	}

setup_master_XML <- function()
	{
	# Sections of a standard Beast2 XML file for tip-dating:
	# 1. Header
	# 2. Taxa
	# 3. Sequences
	# 4. Partitions
	# 5. Miscellaneous
	# 6. Site models
	# 7. Clock model
	# 8. Starting tree
	# 9. Tree model
	# 10. MCMC states
	# 11. Priors
	# 12. Likelihoods
	# 13. Operators
	# 14. Trace Log
	# 15. Screen Log
	# 16. Tree Log
	# 17. Substitutions Tree Log
	
	# Setup
	master_list_of_sections = NULL
	blankline = bl()
	
	
	# Header
	txt = " XML SECTION 1: Header and Beast2 Java class references / mappings "
	section_description = xmlCommentNode(txt)
	
	header = list(blankline, blankline, section_description, blankline)
	header


	# Taxa section
	txt1 = " XML SECTION 2: Taxa and clade definitions (and tip-dates if desired) "
	txt2 = "                (Anything you want to constrain and/or log)           "
	section_description1 = xmlCommentNode(txt1)
	section_description2 = xmlCommentNode(txt2)
	
	taxa = list(blankline, blankline, section_description1, section_description2, blankline)
	taxa
	
	# Sequences
	txt = " XML SECTION 3: Sequence alignments (e.g. DNA, morphology); filtered in later section to produce partitions "
	section_description = xmlCommentNode(txt)
	
	sequences = list(blankline, blankline, section_description, blankline)
	sequences


	# Partitions
	txt1 = " XML SECTION 4: Partitions (filters applied to the sequence alignments produce partitions) "
	txt2 = "                The tree+model confer likelihood on each partition.                        "
	section_description1 = xmlCommentNode(txt1)
	section_description2 = xmlCommentNode(txt2)
	
	partitions = list(blankline, blankline, section_description1, section_description2, blankline)
	partitions


	# Miscellaneous section
	txt1 = " XML SECTION 5: Miscellaneous                                         "
	txt2 = "                (space for various extra parameters etc.)             "
	section_description1 = xmlCommentNode(txt1)
	section_description2 = xmlCommentNode(txt2)
	
	misc = list(blankline, blankline, section_description1, section_description2, blankline)
	misc


	# Site models section
	txt1 = " XML SECTION 6: Site Models (sequence/morphology evolution models for each partition) "
	section_description1 = xmlCommentNode(txt1)
	
	sitemodels = list(blankline, blankline, section_description1, blankline)
	sitemodels


	# Clock model section
	txt1 = " XML SECTION 7: Shared clock model (strict, relaxed, etc.) "
	section_description1 = xmlCommentNode(txt1)
	
	clock = list(blankline, blankline, section_description1,  blankline)
	clock



	# Starting tree section
	txt1 = " XML SECTION 8: Starting tree (random, user-fixed, etc.) "
	section_description1 = xmlCommentNode(txt1)
	
	starting_tree = list(blankline, blankline, section_description1,  blankline)
	starting_tree


	# Tree model section
	txt1 = " XML SECTION 9: Shared tree model (Yule, Birth-Death (BD), Birth-Death Skyline (BDSKY), etc.) "
	section_description1 = xmlCommentNode(txt1)
	
	tree = list(blankline, blankline, section_description1,  blankline)
	tree
	

	# MCMC state section
	txt1 = " XML SECTION 10: MCMC states: the state of the MCMC chain saved at various timepoints "
	txt2 = "                 (as a .state file, to allow the chain to be re-started later on)          "
	section_description1 = xmlCommentNode(txt1)
	section_description2 = xmlCommentNode(txt2)
	
	state = list(blankline, blankline, section_description1, section_description2, blankline)
	state


	# Priors section
	txt1 = " XML SECTION 11: Priors "
	section_description1 = xmlCommentNode(txt1)
	
	priors = list(blankline, blankline, section_description1,  blankline)
	priors


	# Likelihoods section
	txt1 = " XML SECTION 12: Likelihoods "
	section_description1 = xmlCommentNode(txt1)
	
	likes = list(blankline, blankline, section_description1,  blankline)
	likes
	
	
	# Operators section
	txt1 = " XML SECTION 13: Operators. These specify how the parameter values can change at MCMC time-steps. "
	section_description1 = xmlCommentNode(txt1)
	
	operators = list(blankline, blankline, section_description1)
	operators
	
	
	# Trace log section
	txt1 = " XML SECTION 14: Trace log. This logs the parameters in a .log file, which may be viewed in Tracer. "
	section_description1 = xmlCommentNode(txt1)
	
	tracelog = list(blankline, blankline, section_description1)
	tracelog


	# Screen log section
	txt1 = " XML SECTION 15: Screen log. Which parameters print to screen, and how often? "
	section_description1 = xmlCommentNode(txt1)
	
	screenlog = list(blankline, blankline, section_description1,  blankline)
	screenlog	


	# Tree log section
	txt1 = " XML SECTION 16: Tree log. How often should time trees from the MCMC chain be saved, and what data should they have? "
	txt1a = "       (these are trees where the branch lengths are in time units, e.g. millions of years )                    "
	txt2 = "   NOTE: this file can get HUGE if you sample too frequently. You just want a few thousand samples at the end.       "   
	section_description1 = xmlCommentNode(txt1)
	section_description1a = xmlCommentNode(txt1a)
	section_description2 = xmlCommentNode(txt2)
	
	treelog = list(blankline, blankline, section_description1, section_description1a, section_description2, blankline)
	treelog


	# Substitutions Tree log section
	txt1 = " XML SECTION 17: Substitutions Tree log. How often should these trees be saved, and what data should they have?      "
	txt1a = "         (these are trees where the branch lengths are the expected number of substitutions, i.e. rate * time)      "
	txt2 = "   NOTE: this file can get HUGE if you sample too frequently. You just want a few thousand samples at the end.      "
	section_description1 = xmlCommentNode(txt1)
	section_description1a = xmlCommentNode(txt1a)
	section_description2 = xmlCommentNode(txt2)
	
	subslog = list(blankline, blankline, section_description1, section_description1a, section_description2, blankline)
	subslog

	
	# Make a list of the sections
	master_list_of_sections$header = header
	master_list_of_sections$taxa = taxa
	master_list_of_sections$sequences = sequences
	master_list_of_sections$partitions = partitions
	master_list_of_sections$misc = misc
	master_list_of_sections$sitemodels = sitemodels
	master_list_of_sections$clock = clock
	master_list_of_sections$starting_tree = starting_tree
	master_list_of_sections$tree = tree
	master_list_of_sections$state = state
	master_list_of_sections$priors = priors
	master_list_of_sections$likes = likes
	master_list_of_sections$operators = operators
	master_list_of_sections$tracelog = tracelog
	master_list_of_sections$screenlog = screenlog
	master_list_of_sections$treelog = treelog
	master_list_of_sections$subslog = subslog
	master_list_of_sections
	
	return(master_list_of_sections)
	}
