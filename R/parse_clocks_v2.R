#######################################################
# Clock models
# (2014-10-24: default clock model is:
# uncorrelated, logNormal, relaxed clock (UCLD)
#  shared across branches)
#######################################################

# Get the shared clock name(s)
get_clock_names <- function(seqs_df)
	{
	# Get the lines we are actually using
	TF = seqs_df$use == "yes"
	seqs_df2 = seqs_df[TF,]
	
	clockModel_names = unique(seqs_df2$clockModel_name)
	# Alternate spelling
	if (is.null(clockModel_names))
		{
		clockModel_names = unique(seqs_df2$clockmodel_name)
		}
	
	# Find the first matches
	#rownums = match(x=clockModel_names, table=seqs_df2$clockModel_name)
	
	return(clockModel_names)
	} # END get_clock_names(seqs_df)

# Get the shared clock types(s)
get_clock_types <- function(seqs_df)
	{
	# Get the lines we are actually using
	TF = seqs_df$use == "yes"
	seqs_df2 = seqs_df[TF,]
	
	clockModel_names = unique(seqs_df2$clockModel_name)
	
	# Find the first matches
	rownums = match(x=clockModel_names, table=seqs_df2$clockModel_name)
	
	clock_types = seqs_df2$clockmodel_type[rownums]
	 
	
	# ERROR check
	if (sum(isblank_TF(clock_types)) > 0)
		{
		txt = paste0("STOP ERROR in get_clock_names(): in 'seqs_df', every new unique clockModel_name needs a clockmodel_type on the first line of that clockModel_name.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		cat("Printing clock names and types from seqs_df:")
		cat("\n\n")
		clockModel_name = seqs_df2$clockModel_name
		clockmodel_type = seqs_df2$clockmodel_type
		print(as.data.frame(cbind(clockModel_name, clockmodel_type)))
		cat("\n\n")

		stop(txt)
		} # END if (sum(isblank_TF(clock_types)) > 0)
	
	return(clock_types)
	} # END get_clock_names(seqs_df)






######################################################
# CHECK seqs_df for missing clockModel_name or clock_type
######################################################
check_seqs_df <- function(seqs_df)
	{
	# Get the names of each clock (usually there is just 1)
	clockModel_names = get_clock_names(seqs_df)

	# ERROR CHECK: missing clockModel names
	if (any(isblank_TF(seqs_df$clockModel_name)) == TRUE)
		{
		TF = isblank_TF(seqs_df$clockModel_name)
	
		stoptxt = paste0("STOP ERROR #1 in check_clockmodel_setup(): all of the cells in the column 'clockModel_name' of worksheet 'data' MUST be filled in. This is how BEASTmasteR determines which clock model each partition is following. However, you have ", sum(TF), " blank cells. Please fix and re-run.")
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		stop(stoptxt)
		}
	
	
	
	# ERROR FIX
	# Fill out the clock_type if some left blank
	for (i in 1:length(clockModel_names))
		{
		rowsTF = seqs_df$clockModel_name == clockModel_names[i]
		seqs_df_subset = seqs_df[rowsTF,]
		
		# ERROR CHECK: have relRates for morph partitions
		if (nrow(seqs_df_subset) > 1)
			{
			if (any(isblank_TF(seqs_df_subset$clockModel_relRates))==TRUE)
				{
				txt = paste0("STOP ERROR in check_seqs_df(): in worksheet 'data', in clockModel '", clockModel_names[i], "', you have more than one partition, but some relRates are missing. Please fix to proceed.")
				cat("\n\n")
				cat(txt)
				cat("\n\n")
				stop(txt)
				} # END if (any(isblank_TF(seqs_df_subset$clockModel_relRates))==TRUE)
			} # END if (nrow(seqs_df_subset) > 1)


	
		# If the first one is blank, stop.
		clockmodel_type_to_copy = seqs_df_subset$clockmodel_type[1]
		if (isblank_TF(clockmodel_type_to_copy))
			{
			stoptxt = paste0("STOP ERROR #2 in check_clockmodel_setup(): In your Excel settings file, worksheet 'data', for a particular 'clockModel_name' category, at least the FIRST cell in the column 'clockmodel_type' MUST be filled in. This is how BEASTmasteR determines which clock type that particular clockModel_name is using.  Please fix and re-run.")
			cat("\n\n")
			cat(stoptxt)
			cat("\n\n")
			stop(stoptxt)
			}
	
		# Fill in for the rest
		seqs_df_subset$clockmodel_type = clockmodel_type_to_copy
	
		# Stick back in the main seqs_df
		seqs_df[rowsTF,] = seqs_df_subset
		} # END for (i in 1:length(clockModel_names))
	return(seqs_df)
	}




# Define a relaxed clock, and starting values
define_a_shared_clock <- function(seqs_df, ntaxa=NULL, clockModel_name="shared_clock", species_tree_name="shared_tree", gene_tree_name=NULL, log_rateCategories=FALSE, xml=NULL, first=NULL)
	{
	defaults='
	# The dimension of stateNode is the (number of branches - 1)
	length(tr$tip.label)
	relaxed_clock_dimension = length(tr$tip.label) + tr$Nnode - 1
	
	# Or
	relaxed_clock_dimension = ntaxa + (ntaxa-1) - 1
	'
	
	# Error check
	# (make sure you've removed "no" rows from seqs_df)
	TF = seqs_df$use == "no"
	if (sum(TF) > 0)
		{
		txt = paste0("STOP ERROR in define_a_shared_clock(): in 'seqs_df' input, all rows where seqs_df$use=='no' need to be removed from seqs_df.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if (sum(TF) > 0)
	

	# Find the first match
	rownum = match(x=clockModel_name, table=seqs_df$clockModel_name)
	
	# Error check
	# If NA (clockModel not found), stop
	if (is.na(rownum))
		{
		txt = paste0("STOP ERROR in define_a_shared_clock(): clockModel_name=", clockModel_name, " but no match found in seqs_df$clockModel_name.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if (is.na(rownum))
	
	
	# Get the "clock_df" row
	clock_df = seqs_df[rownum, ]
	clock_type = clock_df$clockmodel_type	# e.g. "ucld"
	
	clock_types_allowed = c("strict", "ucld", "uced", "rlc", "SpeciesTreeUCLN", "SpeciesTreeUCED", "SpeciesTreeRLC")
	if (tolower(clock_type) %in% tolower(clock_types_allowed) == FALSE)
		{
		txt = paste0("STOP ERROR in define_a_shared_clock(): clock_type=", clock_type, " but no match found in clock_types_allowed.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		cat("Printing 'clock_types_allowed':\n\n")
		print(clock_types_allowed)
		cat("\n\n")
		stop(txt)
		} # END if (clock_type %in% clock_types_allowed == FALSE)
	
	if (grepl(pattern=tolower("SpeciesTree"), x=tolower(clock_type)) == TRUE)
		{
		tree_name = species_tree_name
		tree_name_idref = paste0("@", tree_name)
		gene_tree_name_idref = paste0("@", gene_tree_name)
		if (is.null(gene_tree_name) == TRUE)
			{
			txt = "STOP ERROR in define_a_shared_clock(): You have specified a SpeciesTree-level clock model in 'clock_type'. This requires that gene_tree_name be input, but instead gene_tree_name=NULL."
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			}
		} else {
		tree_name = species_tree_name
		tree_name_idref = paste0("@", tree_name)
		}

	# Make the tree model ID
	species_tree_name_idref = paste0("@", species_tree_name)
	
	
	# Intialize mean and SD of the relaxed clock model
	# (SD is not used in some clock models)
	mean_of_shared_clock_id = paste0(clock_type, "_Mean_of_", clockModel_name)
	stdev_of_shared_clock_id = paste0(clock_type, "_Stdev_of_", clockModel_name)
	
	# State & starting state of the relaxed clock model
	shared_clock_title_txt = paste0(" Shared clock model, named: ", clockModel_name, " ")
	shared_clock_title_XML = xmlCommentNode(shared_clock_title_txt)
	
	# Get the starting value of the clock rate and clock SD
	starting_clockrate = clock_df$clockrate_starting_val
	if (isblank_TF(starting_clockrate) == TRUE)
		{
		starting_clockrate = 0.004
		} # END if (isblank_TF(starting_clockrate) == TRUE)

	starting_clockSD = clock_df$clockSD_starting_val
	if (isblank_TF(starting_clockSD) == TRUE)
		{
		starting_clockSD = starting_clockrate * 0.5
		} # END if (isblank_TF(starting_clockrate) == TRUE)
	

	#######################################################
	# Cutting -- redundant
	#######################################################
	# Get param1 and param2 for clock mean
# 	dfline = clock_df
# 	tmpXML = make_generic_XML_prior(dfline, colname_prefix="clockrate", param_name=mean_of_shared_clock_id, header_scheme=1, distrib=NULL, param1=NULL, param2=NULL, tmp_offset=NULL, meanInRealSpace=NULL)
# 	priors_on_mean_of_shared_clock_XML = tmpXML$priors
# 
# 	tmpXML = make_generic_XML_prior(dfline, colname_prefix="clockSD", param_name=stdev_of_shared_clock_id, header_scheme=1, distrib=NULL, param1=NULL, param2=NULL, tmp_offset=NULL, meanInRealSpace=NULL)
# 	priors_on_stdev_of_shared_clock_XML = tmpXML$priors
priors_on_mean_of_shared_clock_XML = NULL
priors_on_stdev_of_shared_clock_XML = NULL


	# Min and max on clock rate - defaults
	if (isblank_TF(clock_df$clockrate_min))
		{
		clock_df$clockrate_min = 0.0000001
		}
	if (isblank_TF(clock_df$clockrate_max))
		{
		clock_df$clockrate_max = 10
		}
	if (isblank_TF(clock_df$clockSD_min))
		{
		clock_df$clockSD_min = 0.0000001
		}
	if (isblank_TF(clock_df$clockSD_max))
		{
		clock_df$clockSD_max = 10
		}
	
	# Define XML for relaxed clock, and starting values -- for STATES
	mean_of_shared_clock_XML = xmlNode(name="parameter", starting_clockrate, attrs=list(id=mean_of_shared_clock_id, lower=clock_df$clockrate_min, upper=clock_df$clockrate_max, name="stateNode") )
	stdev_of_shared_clock_XML = xmlNode(name="parameter", starting_clockSD, attrs=list(id=stdev_of_shared_clock_id, lower=clock_df$clockSD_min, upper=clock_df$clockSD_max, name="stateNode") )




	# Strict clock
	if (tolower(clock_type) == "strict")
		{
		# Setup: stateNodes

		clock_type_txt = paste0(" strict: Strict clock with single branch rate ")
		clock_type_XML = xmlCommentNode(clock_type_txt)
		
		# The dimension of the relaxed clock stateNode is the (number of branches - 1)
		clock_rate_mean_idref = paste0("@", mean_of_shared_clock_id)

		state_XMLs = list(bl(), shared_clock_title_XML, bl(), clock_type_XML, mean_of_shared_clock_XML)


		#######################################################
		# Branch rate model -- used in the likelihood calculation of each partition
		#######################################################
		# XML for branchRateModel
		branchRateModel_XML = xmlNode(name="branchRateModel", attrs=list(id=clockModel_name, clock.rate=clock_rate_mean_idref, spec="beast.evolution.branchratemodel.StrictClockModel") )

		# Make into a list
		branchRateModel_XMLs = list(bl(), bl(), shared_clock_title_XML, bl(), clock_type_XML, branchRateModel_XML)



		#######################################################
		# Priors on clock model parameters
		#######################################################
		clock_priors_XMLcomment = xmlCommentNode(" Priors on clock model parameters ")
	
		# Get the prior distributions, the clock rate operators, and the 
		# logs, for this clock
		tmpXML = clock_df_row_to_XML_distribution(clock_df, tree_name=tree_name)
		clockModel_priors_XML = c(tmpXML$priors, priors_on_mean_of_shared_clock_XML)
		clockModel_screenLog_XML = tmpXML$screenlog
		clockModel_traceLog_XML = tmpXML$tracelog


		#######################################################
		# Operators on clock rate
		#######################################################
		ClockRateScaler_operator_id = paste0("ClockRateScaler_", clockModel_name)
		ClockRateScaler_operator_XML = xmlNode(name="operator", attrs=list(id=ClockRateScaler_operator_id, parameter=clock_rate_mean_idref, spec="ScaleOperator", weight="3.0", scaleFactor="0.5") )
	

		rateCategories_xmlComment = xmlCommentNode(paste0(" Operators on Strict Clock (strict) rate (mean): ", clockModel_name, " "))
		clockModel_operators_XML = c(
			list(bl(), bl(), rateCategories_xmlComment), 
			tmpXML$operators,
			list(ClockRateScaler_operator_XML)
			) # END clockModel_operators_XML
		
		# Log 
		clock_rate_mean_log_txt_XML = xmlNode(name="log", attrs=list(idref=mean_of_shared_clock_id))
		clockModel_traceLog_XML = c(clockModel_traceLog_XML, list(bl(), bl(), clock_rate_mean_log_txt_XML))
		
		} # END if (clock_type == "strict")





	# LogNormal relaxed clock ON THE SPECIES TREE
	if (tolower(clock_type) == tolower("SpeciesTreeUCLN")) 
		{
		# Setup: stateNodes

		clock_type_txt = paste0(" SpeciesTreeUCLN: Uncorrelated relaxed clock with lognormally-distributed branch rates ON THE SPECIES TREE; gene tree rates inherited ")
		clock_type_XML = xmlCommentNode(clock_type_txt)

		clock_rate_mean_idref = paste0("@", mean_of_shared_clock_id)
		
		# The dimension of the relaxed clock stateNode is the (number of branches - 1)
		relaxed_clock_dimension = ntaxa + (ntaxa-1) - 1
		rateCategories_of_shared_clock_id = paste0("branchRates_for_SpeciesTree_", species_tree_name, "_under_", clock_type, "_clock_named_", clockModel_name)
		rateCategories_of_shared_clock_idref = paste0("@", rateCategories_of_shared_clock_id)

		#rateCategories_of_shared_clock_XML = xmlNode(name="stateNode", 1, attrs=c(id=rateCategories_of_shared_clock_id, dimension=relaxed_clock_dimension, spec="parameter.IntegerParameter"))
		rateCategories_of_shared_clock_XML = xmlNode(name="stateNode", 1, attrs=c(id=rateCategories_of_shared_clock_id, spec="parameter.IntegerParameter"))

		state_XMLs = list(bl(), shared_clock_title_XML, bl(), clock_type_XML, mean_of_shared_clock_XML, stdev_of_shared_clock_XML, rateCategories_of_shared_clock_XML)


		#######################################################
		# Branch rate prior -- used in the likelihood calculation of each partition
		# This prior is what is lognormally distributed: it is a prior on the
		# mean clock rate
		#######################################################
		# Uncorrelated relaxed clock with lognormally-distributed branch rates
		# XML for mean inside LogNormal distribution inside branchRateModel
# 		branchRateModel_logNormal_mean_XMLcomment = xmlCommentNode(" I have no idea why BEAUTi puts this fixed mean in, or if it even uses it; presumably overridden by clock.rate. ")
# 	
# 		branchRateModel_logNormal_mean_id = paste0("apparently_fixed_mean_for_logNormal_inside_", clockModel_name)
# 		branchRateModel_logNormal_mean_XML = xmlNode(name="parameter", 1.0, attrs=list(id=branchRateModel_logNormal_mean_id, estimate="false", lower="0.0", upper="1.0", name="M") )
# 	
# 		# XML for LogNormal distribution inside branchRateModel
# 		branchRateModel_logNormal_id = paste0("logNormal_distribution_inside_", clockModel_name)
# 
# 		branchRateModel_logNormal_id_XML = xmlNode(name="LogNormal", attrs=list(id=branchRateModel_logNormal_id, S=branchRateModel_logNormal_stdev_idref, meanInRealSpace="true", name="distr"), .children=list(branchRateModel_logNormal_mean_XMLcomment, branchRateModel_logNormal_mean_XML) )
# 		branchRateModel_logNormal_XMLcomment = xmlCommentNode(" logNormal distribution inside branchRateModel ")

		# XML for branchRateModel
		branchRateModel_logNormal_stdev_idref = paste0("@", stdev_of_shared_clock_id)
		branchRateModel_XML_comment = xmlCommentNode(paste0(" speciesTreeRates are specified in clock model: ", clockModel_name, " "))
		branchRateModel_XML = xmlNode(name="speciesTreeRates", attrs=list(id=clockModel_name, estimateRoot="true", spec="starbeast2.UncorrelatedRates", tree=species_tree_name_idref, rates=rateCategories_of_shared_clock_idref, stdev=branchRateModel_logNormal_stdev_idref))

		# Make into a list
		branchRateModel_XMLs = list(bl(), bl(), shared_clock_title_XML, bl(), clock_type_XML, branchRateModel_XML_comment, branchRateModel_XML)



		#######################################################
		# Priors on speciesTree clock model parameters
		#######################################################
		clock_priors_XMLcomment = xmlCommentNode(" Priors on clock model parameters ")
	
		# Get the prior distributions, the clock rate operators, and the 
		# logs, for this clock
		tmpXML = clock_df_row_to_XML_distribution(clock_df, tree_name=species_tree_name)
		clockModel_priors_XML = c(tmpXML$priors, priors_on_mean_of_shared_clock_XML, priors_on_stdev_of_shared_clock_XML)
		clockModel_screenLog_XML = tmpXML$screenlog
		clockModel_traceLog_XML = tmpXML$tracelog


		#######################################################
		# Operators on clockModel branchRates
		#######################################################
		rateCategories_xmlComment = xmlCommentNode(paste0(" Operators on a Species-Tree Level, Uncorrelated Relaxed Clock, with Lognormally-Distributed branch rates (ucld) model parameters for clock model: ", clockModel_name, " "))
		
		rateCategories_branchRatesCycle_operator_id = paste0("branchRatesCycle_for_SpeciesTree_", species_tree_name, "_under_clockModel_", clockModel_name)
		rateCategories_branchRatesCycle_operator_XML = xmlNode(name="operator", attrs=list(id=rateCategories_branchRatesCycle_operator_id, treeRates=rateCategories_of_shared_clock_idref, spec="starbeast2.DiscreteRateCycle", weight="10.0", k="2", optimise="false") )
	
	
		rateCategories_Uniform_operator_id = paste0("branchRatesUniform_operator_for_", species_tree_name, "_under_clockModel_", clockModel_name)
		rateCategories_Uniform_operator_XML = xmlNode(name="operator", attrs=list(id=rateCategories_Uniform_operator_id, treeRates=rateCategories_of_shared_clock_idref, spec="starbeast2.DiscreteRateUniform", weight="10.0", optimise="false") )
	
		clockModel_operators_XML = cl(bl(), bl(), rateCategories_xmlComment, 			tmpXML$operators, rateCategories_branchRatesCycle_operator_XML, rateCategories_Uniform_operator_XML)

		# Add logs for Indicators_idref and LocalRates_idref, but leave them commented out
		rateCategories_logtxt_XML = xmlCommentNode(" Log rate categories, if desired (may be a ridiculous number of columns). ")
		rateCategories_logtxt_XML2 = xmlCommentNode(" (Note: rates are typically also logged in the treeLog) ")

		if (log_rateCategories == FALSE)
			{
			rateCategories_log_txt = paste0('<log idref="', rateCategories_of_shared_clock_id, '"/>')
			rateCategories_log_txt_XML = xmlCommentNode(rateCategories_log_txt)
			
			# Add to the log file
			} else {
			
			rateCategories_log_txt_XML = xmlNode(name="log", attrs=list(idref=rateCategories_of_shared_clock_id))
			
			} # END if (log_rateCategories == FALSE)
		clockModel_traceLog_XML = c(clockModel_traceLog_XML, list(bl(), bl(), rateCategories_logtxt_XML, rateCategories_logtxt_XML2, rateCategories_log_txt_XML))

		} # END if (tolower(clock_type) == tolower("SpeciesTreeUCLN"))
	


	
	
	# LogNormal relaxed clock
	if (tolower(clock_type) == "ucld")
		{
		# Setup: stateNodes

		clock_type_txt = paste0(" ucld: Uncorrelated relaxed clock with lognormally-distributed branch rates ")
		clock_type_XML = xmlCommentNode(clock_type_txt)
		
		# The dimension of the relaxed clock stateNode is the (number of branches - 1)
		relaxed_clock_dimension = ntaxa + (ntaxa-1) - 1
		rateCategories_of_shared_clock_id = paste0("rateCategories_of_", clockModel_name)
		clock_rate_mean_idref = paste0("@", mean_of_shared_clock_id)
		rateCategories_idref = paste0("@", rateCategories_of_shared_clock_id)

		rateCategories_of_shared_clock_XML = xmlNode(name="stateNode", 1, attrs=c(id=rateCategories_of_shared_clock_id, dimension=relaxed_clock_dimension, spec="parameter.IntegerParameter"))


		state_XMLs = list(bl(), shared_clock_title_XML, bl(), clock_type_XML, mean_of_shared_clock_XML, stdev_of_shared_clock_XML, rateCategories_of_shared_clock_XML)


		#######################################################
		# Branch rate model -- used in the likelihood calculation of each partition
		#######################################################
		# Uncorrelated relaxed clock with lognormally-distributed branch rates
		# XML for mean inside LogNormal distribution inside branchRateModel
		branchRateModel_logNormal_mean_XMLcomment = xmlCommentNode(" I have no idea why BEAUTi puts this fixed mean in, or if it even uses it; presumably overridden by clock.rate. ")
	
		branchRateModel_logNormal_mean_id = paste0("apparently_fixed_mean_for_logNormal_inside_", clockModel_name)
		branchRateModel_logNormal_mean_XML = xmlNode(name="parameter", 1.0, attrs=list(id=branchRateModel_logNormal_mean_id, estimate="false", lower="0.0", upper="1.0", name="M") )
	
		# XML for LogNormal distribution inside branchRateModel
		branchRateModel_logNormal_id = paste0("logNormal_distribution_inside_", clockModel_name)
		branchRateModel_logNormal_stdev_idref = paste0("@", stdev_of_shared_clock_id)
		branchRateModel_logNormal_id_XML = xmlNode(name="LogNormal", attrs=list(id=branchRateModel_logNormal_id, S=branchRateModel_logNormal_stdev_idref, meanInRealSpace="true", name="distr"), .children=list(branchRateModel_logNormal_mean_XMLcomment, branchRateModel_logNormal_mean_XML) )
		branchRateModel_logNormal_XMLcomment = xmlCommentNode(" logNormal distribution inside branchRateModel ")

		# XML for branchRateModel
		branchRateModel_XML = xmlNode(name="branchRateModel", attrs=list(id=clockModel_name, clock.rate=clock_rate_mean_idref, rateCategories=rateCategories_idref, spec="beast.evolution.branchratemodel.UCRelaxedClockModel", tree=tree_name_idref), .children=list(bl(), branchRateModel_logNormal_XMLcomment, branchRateModel_logNormal_id_XML))

		# Make into a list
		branchRateModel_XMLs = list(bl(), bl(), shared_clock_title_XML, bl(), clock_type_XML, branchRateModel_XML)



		#######################################################
		# Priors on clock model parameters
		#######################################################
		clock_priors_XMLcomment = xmlCommentNode(" Priors on clock model parameters ")
	
		# Get the prior distributions, the clock rate operators, and the 
		# logs, for this clock
		tmpXML = clock_df_row_to_XML_distribution(clock_df, tree_name=species_tree_name)
		clockModel_priors_XML = c(tmpXML$priors, priors_on_mean_of_shared_clock_XML, priors_on_stdev_of_shared_clock_XML)
		clockModel_screenLog_XML = tmpXML$screenlog
		clockModel_traceLog_XML = tmpXML$tracelog


		#######################################################
		# Operators on rate categories
		#######################################################
		rateCategories_xmlComment = xmlCommentNode(paste0(" Operators on Uncorrelated Relaxed Clock with Lognormally-Distributed branch rates (ucld) model parameters for clock model: ", clockModel_name, " "))

		rateCategories_RandomWalk_operator_id = paste0("CategoriesRandomWalk_", clockModel_name)
		rateCategories_RandomWalk_operator_XML = xmlNode(name="operator", attrs=list(id=rateCategories_RandomWalk_operator_id, parameter=rateCategories_idref, spec="IntRandomWalkOperator", weight="10.0", windowSize="1") )
	
		rateCategories_Swap_operator_id = paste0("CategoriesSwapOperator_", clockModel_name)
		rateCategories_Swap_operator_XML = xmlNode(name="operator", attrs=list(id=rateCategories_Swap_operator_id, intparameter=rateCategories_idref, spec="SwapOperator", weight="10.0") )
	
		rateCategories_Uniform_operator_id = paste0("CategoriesUniformOperator_", clockModel_name)
		rateCategories_Uniform_operator_XML = xmlNode(name="operator", attrs=list(id=rateCategories_Uniform_operator_id, parameter=rateCategories_idref, spec="UniformOperator", weight="10.0") )
	
		clockModel_operators_XML = c(
			list(bl(), bl(), rateCategories_xmlComment), 
			tmpXML$operators,
			list( rateCategories_RandomWalk_operator_XML, rateCategories_Swap_operator_XML, rateCategories_Uniform_operator_XML)
			) # END clockModel_operators_XML

		# Add logs for Indicators_idref and LocalRates_idref, but leave them commented out
		rateCategories_logtxt_XML = xmlCommentNode(" Log rate categories, if desired (may be a ridiculous number of columns). ")
		rateCategories_logtxt_XML2 = xmlCommentNode(" (Note: rates are typically also logged in the treeLog) ")

		if (log_rateCategories == FALSE)
			{
			rateCategories_log_txt = paste0('<log idref="', rateCategories_of_shared_clock_id, '"/>')
			rateCategories_log_txt_XML = xmlCommentNode(rateCategories_log_txt)
			
			# Add to the log file
			} else {
			
			rateCategories_log_txt_XML = xmlNode(name="log", attrs=list(idref=rateCategories_of_shared_clock_id))
			
			} # END if (log_rateCategories == FALSE)
		clockModel_traceLog_XML = c(clockModel_traceLog_XML, list(bl(), bl(), rateCategories_logtxt_XML, rateCategories_logtxt_XML2, rateCategories_log_txt_XML))

		} # END if (clock_type == "ucld")
	

	# Exponential Relaxed Clock
	if (tolower(clock_type) == "uced")
		{
		# Setup: stateNodes

		clock_type_txt = paste0(" uced: Uncorrelated relaxed clock with exponentially-distributed branch rates ")
		clock_type_XML = xmlCommentNode(clock_type_txt)
	
		# The dimension of the relaxed clock stateNode is the (number of branches - 1)
		relaxed_clock_dimension = ntaxa + (ntaxa-1) - 1
		rateCategories_of_shared_clock_id = paste0("rateCategories_of_", clockModel_name)
		clock_rate_mean_idref = paste0("@", mean_of_shared_clock_id)
		rateCategories_idref = paste0("@", rateCategories_of_shared_clock_id)

		rateCategories_of_shared_clock_XML = xmlNode(name="stateNode", 1, attrs=c(id=rateCategories_of_shared_clock_id, dimension=relaxed_clock_dimension, spec="parameter.IntegerParameter"))

		state_XMLs = list(bl(), shared_clock_title_XML, bl(), clock_type_XML, mean_of_shared_clock_XML, rateCategories_of_shared_clock_XML)


		#######################################################
		# Branch rate model -- used in the likelihood calculation of each partition
		#######################################################
		# Uncorrelated relaxed clock with exponentially-distributed branch rates

		# XML for mean inside LogNormal distribution inside branchRateModel
		branchRateModel_Exponential_mean_XMLcomment = xmlCommentNode(" I have no idea why BEAUTi puts this fixed mean in, or if it even uses it; presumably overridden by clock.rate. ")
	
		branchRateModel_Exponential_mean_id = paste0("apparently_fixed_mean_for_Exponential_inside_", clockModel_name)
		branchRateModel_Exponential_mean_XML = xmlNode(name="parameter", 1.0, attrs=list(id=branchRateModel_Exponential_mean_id, estimate="false", lower="0.0", upper="10.0", name="mean") )
	
		# XML for Exponential distribution of branchRates inside branchRateModel
		branchRateModel_Exponential_id = paste0("Exponential_distribution_inside_", clockModel_name)
		branchRateModel_Exponential_id_XML = xmlNode(name="Exponential", attrs=list(id=branchRateModel_Exponential_id, name="distr"), .children=list(branchRateModel_Exponential_mean_XMLcomment, branchRateModel_Exponential_mean_XML) )
		branchRateModel_Exponential_XMLcomment = xmlCommentNode(" Exponential distribution inside branchRateModel ")

		# XML for branchRateModel
		branchRateModel_XML = xmlNode(name="branchRateModel", attrs=list(id=clockModel_name, clock.rate=clock_rate_mean_idref, rateCategories=rateCategories_idref, spec="beast.evolution.branchratemodel.UCRelaxedClockModel", tree=tree_name_idref), .children=list(bl(), branchRateModel_Exponential_XMLcomment, branchRateModel_Exponential_id_XML))

		# Make into a list
		branchRateModel_XMLs = list(bl(), bl(), shared_clock_title_XML, bl(), clock_type_XML, branchRateModel_XML)


		#######################################################
		# Priors on clock model parameters
		#######################################################
		clock_priors_XMLcomment = xmlCommentNode(" Priors on clock model parameters ")
	
		# Get the prior distributions, the clock rate operators, and the 
		# logs, for this clock
		tmpXML = clock_df_row_to_XML_distribution(clock_df, tree_name=species_tree_name)
		clockModel_priors_XML = tmpXML$priors
		clockModel_screenLog_XML = tmpXML$screenlog
		clockModel_traceLog_XML = tmpXML$tracelog


		#######################################################
		# Operators on rate categories
		#######################################################
		rateCategories_RandomWalk_operator_id = paste0("CategoriesRandomWalk_", clockModel_name)
		rateCategories_RandomWalk_operator_XML = xmlNode(name="operator", attrs=list(id=rateCategories_RandomWalk_operator_id, parameter=rateCategories_idref, spec="IntRandomWalkOperator", weight="10.0", windowSize="1") )
	
		rateCategories_Swap_operator_id = paste0("CategoriesSwapOperator_", clockModel_name)
		rateCategories_Swap_operator_XML = xmlNode(name="operator", attrs=list(id=rateCategories_Swap_operator_id, intparameter=rateCategories_idref, spec="SwapOperator", weight="10.0") )
	
		rateCategories_Uniform_operator_id = paste0("CategoriesUniformOperator_", clockModel_name)
		rateCategories_Uniform_operator_XML = xmlNode(name="operator", attrs=list(id=rateCategories_Uniform_operator_id, parameter=rateCategories_idref, spec="UniformOperator", weight="10.0") )
	
		rateCategories_xmlComment = xmlCommentNode(paste0(" Operators on Uncorrelated Relaxed Clock with Exponentially-Distributed branch rates (uced) model parameters for clock model: ", clockModel_name, " "))
	
		clockModel_operators_XML = c(
			list(bl(), bl(), rateCategories_xmlComment), 
			tmpXML$operators,
			list( rateCategories_RandomWalk_operator_XML, rateCategories_Swap_operator_XML, rateCategories_Uniform_operator_XML)
			) # END clockModel_operators_XML



		# Add logs for Indicators_idref and LocalRates_idref, but leave them commented out
		rateCategories_logtxt_XML = xmlCommentNode(" Log rate categories, if desired (may be a ridiculous number of columns). ")
		rateCategories_logtxt_XML2 = xmlCommentNode(" (Note: rates are typically also logged in the treeLog) ")
		if (log_rateCategories == FALSE)
			{
			rateCategories_log_txt = paste0('<log idref="', rateCategories_of_shared_clock_id, '"/>')
			rateCategories_log_txt_XML = xmlCommentNode(rateCategories_log_txt)
			
			# Add to the log file
			} else {
			
			rateCategories_log_txt_XML = xmlNode(name="log", attrs=list(idref=rateCategories_of_shared_clock_id))
			
			} # END if (log_rateCategories == FALSE)
		clockModel_traceLog_XML = c(clockModel_traceLog_XML, list(bl(), bl(), rateCategories_logtxt_XML, rateCategories_logtxt_XML2, rateCategories_log_txt_XML))
		
		} # END if (clock_type == "uced")


	
	#######################################################
	# rlc = "RandomLocalClock"
	#######################################################
	if (tolower(clock_type) == "rlc")
		{
		# Setup: stateNodes
		
		clock_type_txt = paste0(" rlc: Correlated, 'Random Local Clock' rates on branches ")
		clock_type_XML = xmlCommentNode(clock_type_txt)
	
		# The dimension of the relaxed clock stateNode is the (number of branches - 1)
		relaxed_clock_dimension = ntaxa + (ntaxa-1) - 1
		
		# Indicators of Random Local Clock
		Indicators_of_shared_clock_id = paste0("Indicators_of_", clockModel_name)
		clock_rate_mean_idref = paste0("@", mean_of_shared_clock_id)
		Indicators_idref = paste0("@", Indicators_of_shared_clock_id)

		Indicators_of_shared_clock_XML = xmlNode(name="stateNode", "true", attrs=c(id=Indicators_of_shared_clock_id, dimension=relaxed_clock_dimension, spec="parameter.BooleanParameter"))

		# LocalRates of Random Local Clock
		LocalRates_of_shared_clock_id = paste0("LocalRates_of_", clockModel_name)
		clock_rate_mean_idref = paste0("@", mean_of_shared_clock_id)
		LocalRates_idref = paste0("@", LocalRates_of_shared_clock_id)

		LocalRates_of_shared_clock_XML = xmlNode(name="parameter", 0.1, attrs=c(id=LocalRates_of_shared_clock_id, dimension=relaxed_clock_dimension, name="stateNode"))
		
		state_XMLs = list(bl(), shared_clock_title_XML, bl(), clock_type_XML, mean_of_shared_clock_XML, Indicators_of_shared_clock_XML, LocalRates_of_shared_clock_XML)



		#######################################################
		# Branch rate model -- used in the likelihood calculation of each partition
		#######################################################
		# Correlated, Random Local Clock (rlc) model

		# XML for mean inside LogNormal distribution inside branchRateModel
		branchRateModel_Exponential_mean_XMLcomment = xmlCommentNode(" The random local clock has no internal distribution of branch rates ")
	
		# XML for branchRateModel
		branchRateModel_XML = xmlNode(name="branchRateModel", attrs=list(id=clockModel_name, clock.rate=clock_rate_mean_idref, indicators=Indicators_idref, rates=LocalRates_idref, spec="beast.evolution.branchratemodel.RandomLocalClockModel", tree=tree_name_idref) )

		# Make into a list
		branchRateModel_XMLs = list(bl(), bl(), shared_clock_title_XML, bl(), clock_type_XML, branchRateModel_XML)	



		#######################################################
		# Priors on clock model parameters
		#######################################################
		clock_priors_XMLcomment = xmlCommentNode(" Priors on clock model parameters ")
	
		# Get the prior distributions, the clock rate operators, and the 
		# logs, for this clock
		tmpXML = clock_df_row_to_XML_distribution(clock_df, tree_name=species_tree_name)
		clockModel_priors_XML = tmpXML$priors
		clockModel_screenLog_XML = tmpXML$screenlog
		clockModel_traceLog_XML = tmpXML$tracelog


		#######################################################
		# Operators on rate categories -- Random Local Clock model
		#######################################################
		IndicatorsBitFlip_operator_id = paste0("IndicatorsBitFlip_", clockModel_name)
		IndicatorsBitFlip_operator_XML = xmlNode(name="operator", attrs=list(id=IndicatorsBitFlip_operator_id, parameter=Indicators_idref, spec="BitFlipOperator", weight="15.0") )
	
		ClockRateScaler_operator_id = paste0("ClockRateScaler_", clockModel_name)
		ClockRateScaler_operator_XML = xmlNode(name="operator", attrs=list(id=ClockRateScaler_operator_id, parameter=LocalRates_idref, spec="ScaleOperator", weight="15.0", scaleFactor="0.5") )
	
		rateCategories_xmlComment = xmlCommentNode(paste0(" Operators on Random Local Clock model (rlc) parameters for clock model: ", clockModel_name, " "))
	
		clockModel_operators_XML = c(
			list(bl(), bl(), rateCategories_xmlComment), 
			tmpXML$operators,
			list( IndicatorsBitFlip_operator_XML, ClockRateScaler_operator_XML)
			) # END clockModel_operators_XML
		

		# Add logs for Indicators_idref and LocalRates_idref, but leave them commented out
		rateCategories_logtxt_XML = xmlCommentNode(" Log rate categories, if desired (may be a ridiculous number of columns). ")
		rateCategories_logtxt_XML2 = xmlCommentNode(" (Note: rates are typically also logged in the treeLog) ")
		if (log_rateCategories == FALSE)
			{
			Indicators_log_txt = paste0('<log idref="', Indicators_of_shared_clock_id, '"/>')
			LocalRates_log_txt = paste0('<log idref="', LocalRates_of_shared_clock_id, '"/>')
			Indicators_log_txt_XML = xmlCommentNode(Indicators_log_txt)
			LocalRates_log_txt_XML = xmlCommentNode(LocalRates_log_txt)
			
			# Add to the log file
			} else {
			
			Indicators_log_txt_XML = xmlNode(name="log", attrs=list(idref=Indicators_of_shared_clock_id))
			LocalRates_log_txt_XML = xmlNode(name="log", attrs=list(idref=LocalRates_of_shared_clock_id))
			
			} # END if (log_rateCategories == FALSE)
		clockModel_traceLog_XML = c(clockModel_traceLog_XML, list(bl(), bl(), rateCategories_logtxt_XML, rateCategories_logtxt_XML2, Indicators_log_txt_XML, LocalRates_log_txt_XML))

		} # END if (clock_type == "rlc")


	#######################################################
	# Store for output
	#######################################################
	if (is.null(xml))
		{
		res = NULL
		res$state_XMLs = state_XMLs
		res$branchRateModel_XMLs = branchRateModel_XMLs
		res$clockModel_priors_XML = clockModel_priors_XML
		res$clockModel_operators_XML = clockModel_operators_XML
		res$clockModel_traceLog_XML = clockModel_traceLog_XML
		res$clockModel_screenLog_XML = clockModel_screenLog_XML
		
		extract='
		state_XMLs = res$state_XMLs
		branchRateModel_XMLs = res$branchRateModel_XMLs
		clockModel_priors_XML = res$clockModel_priors_XML
		clockModel_operators_XML = res$clockModel_operators_XML
		clockModel_traceLog_XML = res$clockModel_traceLog_XML
		clockModel_screenLog_XML = res$clockModel_screenLog_XML
		'
		return(res)
		} else {
		xml$state = cl(xml$state, state_XMLs)
		xml$clock = cl(xml$clock, branchRateModel_XMLs)
		xml$priors = cl(xml$priors, clockModel_priors_XML)
		xml$operators = cl(xml$operators, clockModel_operators_XML)
		xml$tracelog = cl(xml$tracelog, clockModel_traceLog_XML)
		xml$screenlog = cl(xml$screenlog, clockModel_screenLog_XML)
		return(xml)
		} # END if (is.null(xml))
	
	stop("ERROR in define_a_shared_clock: shouldn't get here")
	} # END define_logNormal_shared_clock <- function(clockModel_name="shared_clock", ntaxa)









