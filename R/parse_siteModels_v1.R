
parse_DNA_AA_siteModels <- function(seqs_df, xml=NULL, tree_name="shared_tree", StarBeast2_TF=FALSE)
	{
	# Subset to remove user="no" lines
	seqs_df = seqs_df[seqs_df$use != "no", ]
	
	
	# Add the relRates for morphology, briefly
	TF = seqs_df$type == "morph"
	if (sum(TF) > 0)
		{
		tmp_morph_seqs_df = seqs_df[TF,]
		morph_relRates = tmp_morph_seqs_df$clockModel_relRates
		morph_relRates = morph_relRates[isblank_TF(morph_relRates) == FALSE]
		#gammaShape_suffixes = tmp_morph_seqs_df$gammaShape_suffix
		
		morph_mutationRate_params_starting_states_XML = c(list(bl()), list(xmlCommentNode(" Starting states of morphological mutationRate relRate parameters ")))
		
		for (i in 1:length(morph_relRates))
			{
			# Make the XML for the mutationRate starting states
			mutationRate_id = morph_relRates[i]
			morph_mutationRate_params_starting_state_XML = xmlNode(name="parameter", 1.0, attrs=list(id=mutationRate_id, name="stateNode") )
			morph_mutationRate_params_starting_states_XML = cl(morph_mutationRate_params_starting_states_XML, morph_mutationRate_params_starting_state_XML)
			}
		} else {
		morph_mutationRate_params_starting_states_XML = NULL
		}
	
	
	keepTF1 = seqs_df$type != "morph"
	keepTF2 = seqs_df$type != "continuous"
	keepTF_sum = (keepTF1 + keepTF2)
	keepTF = keepTF_sum >= 2
	
	# If morphology-only, return null
	if (max(keepTF_sum, na.rm=TRUE) < 2)
		{
		# You still need to set up these exponential distributions
		
		# Make an exponential distribution, for use in setting up the 
		# priors on gammas
		# Store this in "xml$miscellaneous"
		param_id = "exponential_mean_of_1"
		param_XML = xmlNode(name="parameter", 1.0, attrs=list(id=param_id, lower="0.0", upper="0.0", name="mean"))
		exponential_distrib_for_gammaShape_priors_XML_id = "exponential_distrib_for_gammaShape_priors"
		exponential_distrib_for_gammaShape_priors_XML = xmlNode(name="Exponential", attrs=list(id=exponential_distrib_for_gammaShape_priors_XML_id, name="distr"), .children=list(param_XML))
		exponential_distrib_for_gammaShape_priors_XML = c(list(bl()), list(xmlCommentNode(" This is an exponential distribution used as the prior on each gamma shape parameter ")), list(exponential_distrib_for_gammaShape_priors_XML) )

		if (is.null(xml))
			{
			res = NULL
			
			# Miscellaneous
			res$exponential_distrib_for_gammaShape_priors_XML = exponential_distrib_for_gammaShape_priors_XML
			
			extract='
			exponential_distrib_for_gammaShape_priors_XML = res$exponential_distrib_for_gammaShape_priors_XML
			'
			return(res)
			} else {
			# Miscellaneous
			xml$misc = c(xml$misc, exponential_distrib_for_gammaShape_priors_XML)
			return(xml)
			} # END if (is.null(xml))
		} # END if (length(seqs_df$type != "morph") == 0)
	
	# Exclude morphology -- discrete AND continuous
	# For non-morphology datasets, extract partitions from the alignments
	# Remove morphology from seqs
	seqs_df = seqs_df[keepTF, ]
	
	# Make filteredAlignments for each unique partition
	# (partitions might share siteModels or clockModels)
	uniq_partitions = unique(seqs_df$filteredAlignmentName)
	uniq_partitions

	# Get the (first!) rownums of the data partitions
	rownums = match(x=uniq_partitions, table=seqs_df$filteredAlignmentName)
	
	
	
	#######################################################
	# Initialize XML lists for parameters, operators, priors, logs
	#######################################################
	
	# Lists of the parameters for each section
	DNA_AA_siteModels_XML = NULL
	frequencies_starting_states_XMLs = c(list(bl()), list(xmlCommentNode(" Starting states of base frequencies ")))
	estimated_baseFreqs_XMLs = c(list(bl()), list(xmlCommentNode(" Estimated base frequencies (if they are being estimated) ")))
	siteModel_params_starting_states_XML = c(list(bl()), list(xmlCommentNode(" Starting states of siteModel parameters ")))
	siteModel_params_starting_state_XML = NULL
	mutationRate_params_starting_states_XML = c(list(bl()), list(xmlCommentNode(" Starting states of mutationRate parameters ")))
	gamma_starting_states_XML = c(list(bl()), list(xmlCommentNode(" Starting states of gamma parameters ")))
	kappa_starting_states_XML = c(list(bl()), list(xmlCommentNode(" Starting states of kappa parameters (for HKY only) ")))
	
	# Operators
	frequencies_operators_XML = c(list(bl()), list(xmlCommentNode(" Operators for base frequencies ")))
	siteModel_params_operators_XML = c(list(bl()), list(xmlCommentNode(" Operators for siteModel parameters ")))
	siteModel_params_operator_XML = NULL
	mutationRate_params_operators_XML = c(list(bl()), list(xmlCommentNode(" Operators for mutationRate parameters ")))
	gammaShape_params_operators_XML = c(list(bl()), list(xmlCommentNode(" Operators for gammaShape parameters ")))
	kappa_params_operators_XML = c(list(bl()), list(xmlCommentNode(" Operators for kappa parameters (for HKY only) ")))
	
	
	# Priors
	# Frequency priors seem to be left out of BEAUTi-derived XML, so 
	# leaving them out here.
	frequencies_priors_XML = c(list(bl()), list(xmlCommentNode(" Priors on the base frequency parameters seem to be left out of BEAUTi-derived XML,")), list(xmlCommentNode(" so leaving them out here. ")) )
	frequencies_prior_XML = NULL
	kappa_priors_XML = c(list(bl()), list(xmlCommentNode(" Priors on the kappa parameters (for HKY only) ")) )
	siteModel_priors_XML = c(list(bl()), list(xmlCommentNode(" Priors on the siteModel parameters (e.g. kappa) ")) )
	siteModel_prior_XML = NULL
	gammaShape_priors_XML = c(list(bl()), list(xmlCommentNode(" Priors on the gamma shape parameters ")) )
	
	# Likelihoods
	likelihoods_XML = c(list(bl()), list(xmlCommentNode(" Likelihoods of the sequence data (given tree, siteModel, clockModel) ")) )
	# Likelihoods references for the MCMC
	likelihoods_idrefs_XML = c(list(bl()), list(xmlCommentNode(" Likelihoods of the sequence data (given tree, siteModel, clockModel) ")) )
	# Likelihoods traceLogs
	likelihoods_traceLog_XML = c(list(bl()), list(xmlCommentNode(" traceLog of Likelihoods of the sequence data (given tree, siteModel, clockModel) ")) )
	
	
	# Relative mutation rates list of ids
	relative_mutationRate_ids = NULL

	# Log
	log_freqParams_XML = c(list(bl()), list(xmlCommentNode(" Trace Log of the base frequencies parameters ")))
	log_siteModels_XML = c(list(bl()), list(xmlCommentNode(" Trace Log of the siteModel parameters ")))
	log_siteModel_XML = NULL
	log_mutationRates_XML = c(list(bl()), list(xmlCommentNode(" Trace Log of the mutationRate parameters ")))
	log_kappas_XML = c(list(bl()), list(xmlCommentNode(" Trace Log of the kappa parameters ")))
	gammaShape_logs_XML = c(list(bl()), list(xmlCommentNode(" Logs for gammaShape parameters ")))

	# Make an exponential distribution, for use in setting up the 
	# priors on gammas
	# Store this in "xml$miscellaneous"
	param_id = "exponential_mean_of_1"
	param_XML = xmlNode(name="parameter", 1.0, attrs=list(id=param_id, lower="0.0", upper="0.0", name="mean"))
	exponential_distrib_for_gammaShape_priors_XML_id = "exponential_distrib_for_gammaShape_priors"
	exponential_distrib_for_gammaShape_priors_XML = xmlNode(name="Exponential", attrs=list(id=exponential_distrib_for_gammaShape_priors_XML_id, name="distr"), .children=list(param_XML))
	exponential_distrib_for_gammaShape_priors_XML = c(list(bl()), list(xmlCommentNode(" This is an exponential distribution used as the prior on each gamma shape parameter ")), list(exponential_distrib_for_gammaShape_priors_XML) )
	

	

	
	#######################################################
	# Loop through each sequence partition
	#######################################################	
	siteModel_names_used = NULL
	frequencyModel_names_used = NULL
	kappaModel_names_used = NULL
	gammaShapeModel_names_used = NULL
	clockModel_names_used = NULL

	geneTree_clockModel_rate_XML_states = NULL
	geneTree_clockModel_rate_XML = NULL
	geneTree_clockRate_prior_XML = NULL
	clockModel_priors_XML = NULL
	clockModel_screenLog_XML = NULL
	clockModel_traceLog_XML = NULL

	
	for (i in 1:length(rownums))
		{
		rownum = rownums[i]
		partitionName = uniq_partitions[i]
		datasetName = seqs_df$datasetName[rownum]
		if (StarBeast2_TF == TRUE)
			{
			tree_name = seqs_df$geneTreeName[rownum]
			} else {
			tree_name = tree_name
			}
		
		# Actively load the clockModel_name and siteModel_name
		clockModel_name = seqs_df$clockModel_name[rownum]
		siteModel_name = seqs_df$siteModel_name[rownum]
		
		# Get the corresponding clock type
		# First example of the clockModel_name
		if (is.na(clockModel_name) == TRUE)
			{
			stoptxt = paste0("STOP ERROR: during looping through partitions, at i=", i, ", clockModel_name=NA. Make sure you've filled in a clockModel_name for each partition, in worksheet 'data'.")
			cat("\n\n")
			cat(stoptxt)
			cat("\n\n")
			stop(stoptxt)
			} # END if (is.na(clockModel_name) == TRUE)
		
		first_rownum_to_get_clockModel = match(x=clockModel_name, table=seqs_df$clockModel_name)
		clock_type = seqs_df$clockmodel_type[first_rownum_to_get_clockModel]
		
		datatype = seqs_df$type[rownum]
		
		#print(seqs_df)
		#print(datatype)
		if (datatype == "AA")
			{
			dimension = 20
			# Not used here, used in writing out data
			dataType_txt = "aminoacid"
			#txt = "\n\nERROR in parse_DNA_AA_siteModels(): Amino acid (AA) siteModels not yet implemented in BEASTmasteR\n\n"
			#stop(txt)
			}
		if (datatype == "DNA")
			{
			dimension = 4
			}
		
		
		# Proportion of Invariant Sites
		id = paste0("proportionInvariant_", siteModel_name)
		
		# Set to 0, if blank
		if (isblank_TF(seqs_df$pInv[rownum]))
			{
			seqs_df$pInv[rownum] = 0.0
			}
		
		if (seqs_df$pInv[rownum] == "estimate")
			{
			estimate="true"
			pInv_startval = 0.5
			txt = "\n\nERROR in parse_DNA_AA_siteModels(): estimating pInv (proportion of invariant sites)\nhas not been implemented yet in BEASTmasteR.\nYou probably don't want to do it anyway:\nmodels that include pInv often have identifiability problems.\n\n"
			stop(txt)
			
			} else {
			estimate="false"
			pInv_startval = seqs_df$pInv[rownum]
			}
		pInv_XML = xmlNode(name="parameter", pInv_startval, attrs=list(id=id, estimate=estimate, name="proportionInvariant", upper="1.0", lower="0.0") )


		# Frequences reference
		if (isblank_TF(seqs_df$baseFreqs_params_suffix[rownum]) == TRUE)
			{
			# Default option: Auto-name the frequencies
			id = paste0(seqs_df$baseFreqs[rownum], "Frequencies_", siteModel_name)
			frequencies = paste0("@freqParameter_", siteModel_name)
			baseFreqs_idref_XML = xmlNode(name="frequencies", attrs=list(id=id, frequencies=frequencies, spec="Frequencies") )
			
			# Frequencies starting states
			starting_frequencies_id = paste0("starting_freqParameter_", siteModel_name)
			} else {
			# If there is a user-specified name for the base frequencies, use that
			# Default option: Auto-name the frequencies
			# Frequencies starting states
			starting_frequencies_id = paste0("starting_freqParameter_", seqs_df$baseFreqs_params_suffix[rownum])
			starting_frequencies_idref = paste0("@", starting_frequencies_id)
			estimated_baseFreqs_id = paste0(seqs_df$baseFreqs[rownum], "Frequencies_", seqs_df$baseFreqs_params_suffix[rownum])

			baseFreqs_idref_XML = xmlNode(name="frequencies", attrs=list(idref=estimated_baseFreqs_id) )
			
			
			} # END if (isblank_TF(seqs_df$baseFreqs_params_suffix[rownum]) == TRUE)
		frequencyModel_name = estimated_baseFreqs_id

		# Re-use siteModel (and frequencies), if it has already been done. Otherwise, make a new one!
		if ( (frequencyModel_name %in% frequencyModel_names_used) == FALSE)
			{
			# Starting states of base frequencies
			frequencies_starting_state_XML = xmlNode(name="parameter", 1/dimension, attrs=list(id=starting_frequencies_id, dimension=dimension, lower="0.0", upper="1.0", name="stateNode") )

			estimated_baseFreqs_XML = xmlNode(name="frequencies", attrs=list(id=estimated_baseFreqs_id, frequencies=starting_frequencies_idref, spec="Frequencies") )

			} else{
			frequencies_starting_state_XML = NULL
			estimated_baseFreqs_XML = NULL
			}


		
		# Don't do base frequencies for amino acids
		# Or for Jukes-Cantor, JC69
		if ( (datatype == "AA") || (seqs_df$model[rownum] == "JC69") || (seqs_df$model[rownum] == "JC") )
			{
			baseFreqs_idref_XML = NULL
			frequencies_starting_state_XML = NULL
			} # END if (datatype == "AA")
		
		
		###################################################
		# Make XML for each substitution Model
		###################################################
		
		# Implemented models
		DNAmodels = c("JC", "JC69", "HKY", "GTR")
		AAmodels = c("BLOSUM62", "CPREV", "JTT", "MTREV", "WAG", "ratesQ")
		implemented_models = c(DNAmodels, AAmodels)
		###################################################

		model_found_TF = seqs_df$model[rownum] %in% implemented_models
		if (sum(model_found_TF) == 0)
			{
			txt = paste0("\n\nERROR in parse_DNA_siteModels(): model '", seqs_df$model[rownum], "' is not recognized,\nor has not yet been implemented in BEASTmasteR.\nCurrently implemented models are:\n", paste(implemented_models, sep=", ", collapse=""), "\n\n")
			cat(txt)
			stop(txt)
			}


		if ((seqs_df$model[rownum] == "JC69") || (seqs_df$model[rownum] == "JC"))
			{
			id = paste0("JC69_", siteModel_name)
			substModel_XML = xmlNode(name="substModel", attrs=list(id=id, spec="JukesCantor"))

			# For JC69 (constant rate matrix and base frequencies)
			# set these to NULL
			frequencies_starting_state_XML = NULL		
			siteModel_params_starting_state_XML = NULL		
			frequencies_operator_XML = NULL
			siteModel_params_operator_XML = NULL
			siteModel_prior_XML = NULL
			log_freqParam_XML = NULL
			log_siteModel_XML = NULL

			# Operator on substitution model parameter(s)
			} # END if (seqs_df$model[rownum] == "JC69")
		
		
		
		# All the standard AA substitution models can be done easily
		if (seqs_df$model[rownum] %in% AAmodels)
			{
			AAmodel = seqs_df$model[rownum]
			if (AAmodel != "ratesQ")
				{
				# Normal AA rate matrix
				id = paste0(AAmodel, "_", siteModel_name)
				substModel_XML = xmlNode(name="substModel", attrs=list(id=id, spec=AAmodel))
				} else {
				# Implement manual AA rates matrix via ratesQ column!
				# ratesQ: AA Q matrix, manually specified
				# just lower triangle, space-delimited
				# NOT TESTED
				# Beast2 help on optional "rates" entry of e.g. 
				# WAG: "Rate parameter which defines the transition
				# rate matrix. Only the off-diagonal entries need to
				# be specified (diagonal makes row sum to zero in a
				# rate matrix). Entry i specifies the rate from 
				# floor(i/(n-1)) to i%(n-1)+delta where n is the
				# number of states and delta=1 if 
				# floor(i/(n-1)) >= i%(n-1) and 0 otherwise."
				id = paste0(AAmodel, "_", rownum, "_", siteModel_name)
				ratesQ_txt = trim(seqs_df$ratesQ[rownum])
				substModel_XML = xmlNode(name="substModel", attrs=list(id=id, spec="WAG", rates=ratesQ_txt))
				} # if (AAmodel != "ratesQ")
			
			# For AAs (constant rate matrix and base frequencies)
			# set these to NULL
			frequencies_starting_state_XML = NULL		
			siteModel_params_starting_state_XML = NULL		
			frequencies_operator_XML = NULL
			siteModel_params_operator_XML = NULL
			siteModel_prior_XML = NULL
			log_freqParam_XML = NULL
			log_siteModel_XML = NULL


			# Operator on substitution model parameter(s)
			# None
			} # END if (seqs_df$model[rownum] %in% AAmodels)

		



		# XML for the substitution model, using spec HKY, the specified kappa, and specified baseFreqs
		if (seqs_df$model[rownum] == "HKY")
			{
			id = paste0("HKY_", siteModel_name)

			# Kappa parameter for HKY model; starting state and prior
			if ((seqs_df$model[rownum] == "HKY") || (seqs_df$model[rownum] == "GTR"))
				{
				# Default name for kappa parameter? Or user-specified (for re-use across siteModels?)
				if (isblank_TF(seqs_df$kappa_for_HKY_suffix[rownum]) == TRUE)
					{
					# Traditional/default name for kappa
					kappa_id = paste0("kappa_", siteModel_name)
					kappaParam_idref = paste0("@kappa_", siteModel_name)
					} else {
					# User-specified name for kappa
					kappa_id = paste0("kappa_", seqs_df$kappa_for_HKY_suffix[rownum])
					kappaParam_idref = paste0("@kappa_", seqs_df$kappa_for_HKY_suffix[rownum])
					}
				kappaModel_name = kappa_id
				
				# Create XML for the new kappa
				if ((kappaModel_name %in% kappaModel_names_used) == FALSE)
					{
					# Starting value for kappa
					kappa_starting_state_XML = xmlNode(name="parameter", 2.0, attrs=list(id=kappa_id, lower="0.0", name="stateNode") )
			
					# Operator on substitution model parameter(s)
					kappa_scaler_id = paste0(kappa_id, "_Scaler")
					param_ref = paste0("@", kappa_id)
					kappa_param_operator_XML = xmlNode(name="operator", attrs=list(id=kappa_scaler_id, parameter=param_ref, scaleFactor="0.5", spec="ScaleOperator", weight="3") )
			
			
					# Make an logNormal distribution, for use in setting up the 
					# priors on kappas
					# Store this in "xml$miscellaneous"
					partition_idref = paste0("@", siteModel_name)
					param1_id = paste0("logNormal_mean_for_kappa_prior_partition_", siteModel_name)
					param2_id = paste0("logNormal_stddev_for_kappa_prior_partition_", siteModel_name)
					param1_XML = xmlNode(name="parameter", 1.0, attrs=list(id=param1_id, name="M", estimate="false"))
					param2_XML = xmlNode(name="parameter", 1.25, attrs=list(id=param2_id, name="S", estimate="false"))
					logNormal_distrib_for_kappa_priors_XML_id = paste0("logNormal_distrib_for_kappa_priors_", siteModel_name)
					logNormal_distrib_for_kappa_priors_XML = xmlNode(name="LogNormal", attrs=list(id=logNormal_distrib_for_kappa_priors_XML_id, name="distr"), .children=list(param1_XML, param2_XML))
					kappa_prior_id = paste0("prior_for_kappa_", siteModel_name)

					kappa_prior_XML = xmlNode(name="prior", attrs=list(id=kappa_prior_id, name="distribution", x=kappaParam_idref), .children=list(logNormal_distrib_for_kappa_priors_XML))

					# Store the kappa priors in the general siteModel priors
					#siteModel_prior_XML = kappa_prior_XML
			
					# Trace Log of siteModel parameter(s)
					log_kappa_XML = xmlNode(name="parameter", attrs=list(idref=kappa_id, name="log") )
					} else {
					# Don't need new kappa starter, prior, etc.
					siteModel_params_starting_state_XML = NULL
					kappa_param_operator_XML = NULL
					siteModel_params_operator_XML = NULL
					logNormal_distrib_for_kappa_priors_XML = NULL
					siteModel_prior_XML = NULL 
					kappa_prior_XML = NULL
					} # END if ((kappaModel_name %in% kappaModel_names_used) == FALSE)
				} # END if ((seqs_df$model[rownum] == "HKY") || (seqs_df$model[rownum] == "GTR"))

			# Write the substitution model (might not use it, if siteModel has already been used)
			substModel_XML = xmlNode(name="substModel", attrs=list(id=id, kappa=kappaParam_idref, spec="HKY"), .children=list(baseFreqs_idref_XML))
			} # END if (seqs_df$model[rownum] == "HKY")

	#################################################################
	# 
	#################################################################		
		
		if (seqs_df$model[rownum] == "GTR")
			{
			id = paste0("GTR_", siteModel_name)
			
			rateAC_id = paste0("rateAC_", siteModel_name)
			rateAG_id = paste0("rateAG_", siteModel_name)
			rateAT_id = paste0("rateAT_", siteModel_name)
			rateCG_id = paste0("rateCG_", siteModel_name)
			rateCT_id = paste0("rateCT_", siteModel_name)
			rateGT_id = paste0("rateGT_", siteModel_name)

			rateAC_idref = paste0("@", rateAC_id)
			rateAG_idref = paste0("@", rateAG_id)
			rateAT_idref = paste0("@", rateAT_id)
			rateCG_idref = paste0("@", rateCG_id)
			rateCT_idref = paste0("@", rateCT_id)
			rateGT_idref = paste0("@", rateGT_id)

			
			substModel_XML = xmlNode(name="substModel", attrs=list(id=id, rateAC=rateAC_idref, rateAG=rateAG_idref, rateAT=rateAT_idref, rateCG=rateCG_idref, rateCT=rateCT_idref, rateGT=rateGT_idref, spec="GTR"), .children=list(baseFreqs_idref_XML))

			# Substitution Model parameters, starting states
			siteModel_params_starting_state_XML_rateAC = xmlNode(name="parameter", 1.0, attrs=list(id=rateAC_id, lower="0.0", name="stateNode") )
			siteModel_params_starting_state_XML_rateAG = xmlNode(name="parameter", 1.0, attrs=list(id=rateAG_id, lower="0.0", name="stateNode") )
			siteModel_params_starting_state_XML_rateAT = xmlNode(name="parameter", 1.0, attrs=list(id=rateAT_id, lower="0.0", name="stateNode") )
			siteModel_params_starting_state_XML_rateCG = xmlNode(name="parameter", 1.0, attrs=list(id=rateCG_id, lower="0.0", name="stateNode") )
			siteModel_params_starting_state_XML_rateCT = xmlNode(name="parameter", 1.0, attrs=list(id=rateCT_id, lower="0.0", name="stateNode") )
			siteModel_params_starting_state_XML_rateGT = xmlNode(name="parameter", 1.0, attrs=list(id=rateGT_id, lower="0.0", name="stateNode") )
			
			siteModel_params_starting_state_XML = c(
			list(siteModel_params_starting_state_XML_rateAC),
			list(siteModel_params_starting_state_XML_rateAG),
			list(siteModel_params_starting_state_XML_rateAT),
			list(siteModel_params_starting_state_XML_rateCG),
			list(siteModel_params_starting_state_XML_rateCT),
			list(siteModel_params_starting_state_XML_rateGT)
			)
			
			# Operator on substitution model parameter(s)
			rateAC_scaler_id = paste0(rateAC_id, "_Scaler")
			rateAG_scaler_id = paste0(rateAG_id, "_Scaler")
			rateAT_scaler_id = paste0(rateAT_id, "_Scaler")
			rateCG_scaler_id = paste0(rateCG_id, "_Scaler")
			rateCT_scaler_id = paste0(rateCT_id, "_Scaler")
			#rateGT is fixed to 1
			#rateGT_scaler_id = paste0(rateGT_id, "_Scaler")

			siteModel_params_operator_XML_rateAC = xmlNode(name="operator", attrs=list(id=rateAC_scaler_id, parameter=rateAC_idref, scaleFactor="0.5", spec="ScaleOperator", weight="1") )
			siteModel_params_operator_XML_rateAG = xmlNode(name="operator", attrs=list(id=rateAG_scaler_id, parameter=rateAG_idref, scaleFactor="0.5", spec="ScaleOperator", weight="1") )
			siteModel_params_operator_XML_rateAT = xmlNode(name="operator", attrs=list(id=rateAT_scaler_id, parameter=rateAT_idref, scaleFactor="0.5", spec="ScaleOperator", weight="1") )
			siteModel_params_operator_XML_rateCG = xmlNode(name="operator", attrs=list(id=rateCG_scaler_id, parameter=rateCG_idref, scaleFactor="0.5", spec="ScaleOperator", weight="1") )
			siteModel_params_operator_XML_rateCT = xmlNode(name="operator", attrs=list(id=rateCT_scaler_id, parameter=rateCT_idref, scaleFactor="0.5", spec="ScaleOperator", weight="1") )
			#rateGT is fixed to 1
			#siteModel_params_operator_XML_rateGT = xmlNode(name="operator", attrs=list(id=rateGT_scaler_id, parameter=rateGT_idref, scaleFactor="0.5", spec="ScaleOperator", weight="0.1") )
			
			
			
			siteModel_params_operator_XML = c(
			list(siteModel_params_operator_XML_rateAC),
			list(siteModel_params_operator_XML_rateAG),
			list(siteModel_params_operator_XML_rateAT),
			list(siteModel_params_operator_XML_rateCG),
			list(siteModel_params_operator_XML_rateCT)
			#rateGT is fixed to 1
			#siteModel_params_operator_XML_rateGT
			)
			
			
			# Make a Gamma distribution, for use in setting up the 
			# priors on rateAC, rateAG, etc.

			# Prior on rateAC
			param1_id = paste0("alpha_for_rateAC_gammaPrior_partition_", siteModel_name)
			param2_id = paste0("beta_for_rateAC_gammaPrior_partition_", siteModel_name)
			param1_XML = xmlNode(name="parameter", 0.05, attrs=list(id=param1_id, name="alpha", estimate="false"))
			param2_XML = xmlNode(name="parameter", 10.0, attrs=list(id=param2_id, name="beta", estimate="false"))
			Gamma_distrib_for_rateAC_priors_XML_id = paste0("Gamma_distrib_for_rateAC_prior_", siteModel_name)
			Gamma_distrib_for_rateAC_priors_XML_XML = xmlNode(name="Gamma", attrs=list(id=Gamma_distrib_for_rateAC_priors_XML_id, name="distr"), .children=list(param1_XML, param2_XML))
			Gamma_prior_rateAC_id = paste0("prior_distrib_on_rateAC_", siteModel_name)
			partition_idref = paste0("@", siteModel_name)
			Gamma_prior_rateAC_XML = xmlNode(name="prior", attrs=list(id=Gamma_prior_rateAC_id, name="distribution", x=rateAC_idref), .children=list(Gamma_distrib_for_rateAC_priors_XML_XML))
			
			# Prior on rateAG
			param1_id = paste0("alpha_for_rateAG_gammaPrior_partition_", siteModel_name)
			param2_id = paste0("beta_for_rateAG_gammaPrior_partition_", siteModel_name)
			param1_XML = xmlNode(name="parameter", 0.05, attrs=list(id=param1_id, name="alpha", estimate="false"))
			param2_XML = xmlNode(name="parameter", 10.0, attrs=list(id=param2_id, name="beta", estimate="false"))
			Gamma_distrib_for_rateAG_priors_XML_id = paste0("Gamma_distrib_for_rateAG_prior_", siteModel_name)
			Gamma_distrib_for_rateAG_priors_XML_XML = xmlNode(name="Gamma", attrs=list(id=Gamma_distrib_for_rateAG_priors_XML_id, name="distr"), .children=list(param1_XML, param2_XML))
			Gamma_prior_rateAG_id = paste0("prior_distrib_on_rateAG_", siteModel_name)
			partition_idref = paste0("@", siteModel_name)
			Gamma_prior_rateAG_XML = xmlNode(name="prior", attrs=list(id=Gamma_prior_rateAG_id, name="distribution", x=rateAG_idref), .children=list(Gamma_distrib_for_rateAG_priors_XML_XML))
			
			# Prior on rateAT
			param1_id = paste0("alpha_for_rateAT_gammaPrior_partition_", siteModel_name)
			param2_id = paste0("beta_for_rateAT_gammaPrior_partition_", siteModel_name)
			param1_XML = xmlNode(name="parameter", 0.05, attrs=list(id=param1_id, name="alpha", estimate="false"))
			param2_XML = xmlNode(name="parameter", 10.0, attrs=list(id=param2_id, name="beta", estimate="false"))
			Gamma_distrib_for_rateAT_priors_XML_id = paste0("Gamma_distrib_for_rateAT_prior_", siteModel_name)
			Gamma_distrib_for_rateAT_priors_XML_XML = xmlNode(name="Gamma", attrs=list(id=Gamma_distrib_for_rateAT_priors_XML_id, name="distr"), .children=list(param1_XML, param2_XML))
			Gamma_prior_rateAT_id = paste0("prior_distrib_on_rateAT_", siteModel_name)
			partition_idref = paste0("@", siteModel_name)
			Gamma_prior_rateAT_XML = xmlNode(name="prior", attrs=list(id=Gamma_prior_rateAT_id, name="distribution", x=rateAT_idref), .children=list(Gamma_distrib_for_rateAT_priors_XML_XML))
			
			# Prior on rateCG
			param1_id = paste0("alpha_for_rateCG_gammaPrior_partition_", siteModel_name)
			param2_id = paste0("beta_for_rateCG_gammaPrior_partition_", siteModel_name)
			param1_XML = xmlNode(name="parameter", 0.05, attrs=list(id=param1_id, name="alpha", estimate="false"))
			param2_XML = xmlNode(name="parameter", 10.0, attrs=list(id=param2_id, name="beta", estimate="false"))
			Gamma_distrib_for_rateCG_priors_XML_id = paste0("Gamma_distrib_for_rateCG_prior_", siteModel_name)
			Gamma_distrib_for_rateCG_priors_XML_XML = xmlNode(name="Gamma", attrs=list(id=Gamma_distrib_for_rateCG_priors_XML_id, name="distr"), .children=list(param1_XML, param2_XML))
			Gamma_prior_rateCG_id = paste0("prior_distrib_on_rateCG_", siteModel_name)
			partition_idref = paste0("@", siteModel_name)
			Gamma_prior_rateCG_XML = xmlNode(name="prior", attrs=list(id=Gamma_prior_rateCG_id, name="distribution", x=rateCG_idref), .children=list(Gamma_distrib_for_rateCG_priors_XML_XML))
			
			# Prior on rateCT
			param1_id = paste0("alpha_for_rateCT_gammaPrior_partition_", siteModel_name)
			param2_id = paste0("beta_for_rateCT_gammaPrior_partition_", siteModel_name)
			param1_XML = xmlNode(name="parameter", 0.05, attrs=list(id=param1_id, name="alpha", estimate="false"))
			param2_XML = xmlNode(name="parameter", 10.0, attrs=list(id=param2_id, name="beta", estimate="false"))
			Gamma_distrib_for_rateCT_priors_XML_id = paste0("Gamma_distrib_for_rateCT_prior_", siteModel_name)
			Gamma_distrib_for_rateCT_priors_XML_XML = xmlNode(name="Gamma", attrs=list(id=Gamma_distrib_for_rateCT_priors_XML_id, name="distr"), .children=list(param1_XML, param2_XML))
			Gamma_prior_rateCT_id = paste0("prior_distrib_on_rateCT_", siteModel_name)
			partition_idref = paste0("@", siteModel_name)
			Gamma_prior_rateCT_XML = xmlNode(name="prior", attrs=list(id=Gamma_prior_rateCT_id, name="distribution", x=rateCT_idref), .children=list(Gamma_distrib_for_rateCT_priors_XML_XML))
			
			# Store the Gamma priors in the general siteModel priors
			siteModel_prior_XML = c(
			list(Gamma_prior_rateAC_XML),
			list(Gamma_prior_rateAG_XML),
			list(Gamma_prior_rateAT_XML),
			list(Gamma_prior_rateCG_XML),
			list(Gamma_prior_rateCT_XML)
			)
			
			# Trace Log of siteModel parameter(s)
			log_rateAC_XML = xmlNode(name="parameter", attrs=list(idref=rateAC_id, name="log") )
			log_rateAG_XML = xmlNode(name="parameter", attrs=list(idref=rateAG_id, name="log") )
			log_rateAT_XML = xmlNode(name="parameter", attrs=list(idref=rateAT_id, name="log") )
			log_rateCG_XML = xmlNode(name="parameter", attrs=list(idref=rateCG_id, name="log") )
			log_rateCT_XML = xmlNode(name="parameter", attrs=list(idref=rateCT_id, name="log") )
			log_rateGT_XML = xmlNode(name="parameter", attrs=list(idref=rateGT_id, name="log") )

			log_siteModel_XML = c(
			list(log_rateAC_XML),
			list(log_rateAG_XML),
			list(log_rateAT_XML),
			list(log_rateCG_XML),
			list(log_rateCT_XML),
			list(log_rateGT_XML)
			)

			} # END if (seqs_df$model[rownum] == "GTR")
		
		
		substModel_XML
		
		# Make the XML node for this site model
		siteModel_id = paste0("siteModel_", siteModel_name)
		
		# Number of gammaShape categories for sitewise rate variation
		gammaNum = seqs_df$gammaNum[rownum]
		
		# For a StarBeast2 analysis, use the clockModel name + _rate
	#	if (StarBeast2_TF == TRUE)
	#		{
	#		mutationRate_id = paste0(clockModel_name, "_fixed_to_1")
	#		mutationRate_idref = paste0("@", mutationRate_id)
	#		} else {
			# For a concatenated analysis, use a relRate on partitionName

		if (isblank_TF(seqs_df$clockModel_relRates[rownum]) == TRUE)
			{
			# Default name for this relRate
			mutationRate_id = paste0(partitionName, "_relRate")			
			} else {
			# A relRate name was found -- but DON'T use for morphology
			if (seqs_df$type[rownum] == "morph")
				{
				txt = "STOP ERROR in parse_DNA_AA_siteModels(): You have specified 'clockModel_relRates' cells for data of type 'morph'. However, morphology data are automatically assigned relRates (as each number of states gets a different relRate). Change these cells to be blank.\n\nIf you *really* want more specific naming of relRates by number of states, you could make multiple morph partitions manually, one for each number of character states observed in your dataset, manually specifying which column numbers are 2-state, 3-state, etc. See the columns 'startchar', 'endchar', 'by', and 'list'."
				cat("\n\n")
				cat(txt)
				cat("\n\n")
				
				cat("Printing the offending cell in column 'clockModel_relRates':\n")
				print(seqs_df$clockModel_relRates[rownum])
				cat("\n\n")
				
				stop(txt)
				} else {
				mutationRate_id = seqs_df$clockModel_relRates[rownum]
				} # END if (seqs_df$type[rownum] == "morph")
			}
	#		}
		mutationRate_idref = paste0("@", mutationRate_id)

			
		# Make the XML for the Gamma starting states
		if (isblank_TF(seqs_df$gammaShape_suffix[rownum]) == TRUE)
			{
			gammaShape_id = paste0("gammaShape_", siteModel_name)			
			} else {
			gammaShape_id = paste0("gammaShape_", seqs_df$gammaShape_suffix[rownum])
			}
		gammaShape_idref = paste0("@", gammaShape_id)
		gammaShapeModel_name = gammaShape_id

		siteModel_XML = xmlNode(name="siteModel", attrs=list(id=siteModel_id, mutationRate=mutationRate_idref, shape=gammaShape_idref, gammaCategoryCount=gammaNum, spec="SiteModel"), .children=list(pInv_XML, substModel_XML) )
		
		# Make a new starting state, if needed
		if ( (gammaShapeModel_name %in% gammaShapeModel_names_used) == FALSE)
			{
			gamma_starting_state_XML = xmlNode(name="parameter", 1.0, attrs=list(id=gammaShape_id, name="stateNode") )

			# Make the XML for the GammaShape priors
			gamma_prior_id = paste0(gammaShape_id, "_prior")
			distrib_idref = paste0("@", exponential_distrib_for_gammaShape_priors_XML_id)
			gammaShape_idref = paste0("@", gammaShape_id)
			gammaShape_prior_XML = xmlNode(name="prior", attrs=list(id=gamma_prior_id, distr=distrib_idref, name="distribution", x=gammaShape_idref))

			} else {
			gamma_starting_state_XML = NULL
			gammaShape_prior_XML = NULL
			} # END if ( (gammaShapeModel_name %in% gammaShapeModel_names_used) == FALSE)
		
		
		# Make the XML for the mutationRate starting states
		mutationRate_params_starting_state_XML = xmlNode(name="parameter", 1.0, attrs=list(id=mutationRate_id, name="stateNode") )
		
		# Make a list of the mutationRate ids (and the length of each partition,
		# for the weightvector)
		relative_mutationRate_ids = c(relative_mutationRate_ids, mutationRate_id)
		
		# Do base frequencies for DNA, not AA (but not Jukes-Cantor JC69)
		JC69_TF = FALSE
		if ( (seqs_df$model[rownum] == "JC69") || (seqs_df$model[rownum] == "JC") )
			{
			JC69_TF = TRUE
			} # END if ( (seqs_df$model[rownum] == "JC69") etc...
		
		
		
		if ( (datatype == "DNA") && (JC69_TF == FALSE) )
			{
			# If blank, the default is to estimate
			if (isblank_TF(seqs_df$baseFreqs[rownum]) == TRUE)
				{
				seqs_df$baseFreqs[rownum] = "estimate"
				}
			if (seqs_df$baseFreqs[rownum] == "estimate")
				{
				weightval = 10
				}
			if ((seqs_df$baseFreqs[rownum] == "fixed") || (seqs_df$baseFreqs[rownum] == "equal"))
				{
				weightval = 0
				}
			
			# Operator on base frequency parameter(s)
			frequencies_exchanger_id = paste0(starting_frequencies_id, "_Exchanger")
			param_XML = xmlNode(name="parameter", attrs=list(idref=starting_frequencies_id) )
			frequencies_operator_XML = xmlNode(name="operator", attrs=list(id=frequencies_exchanger_id, delta="0.1", spec="DeltaExchangeOperator", weight=weightval), .children=list(param_XML) )
			} else {
			frequencies_operator_XML = NULL
			} # END if ( (datatype == "DNA") && (JC69_TF == FALSE) )
		
		# Operator on gammaShape parameter(s)
		if (gammaNum == 1)
			{
			gammaOperatorWeight = 0
			} else {
			gammaOperatorWeight = 1
			} # END if (gammaNum == 1)

		gammaShape_operator_id = paste0("gammaShapeScaler_for_", gammaShape_id)
		gammaShape_params_operator_XML = xmlNode(name="operator", attrs=list(id=gammaShape_operator_id, parameter=gammaShape_idref, scaleFactor="0.5", spec="ScaleOperator", weight=gammaOperatorWeight) )
		
		# Log the gamma shape
		gammaShape_log_XML = xmlNode(name="log", attrs=list(idref=gammaShape_id))
		
		# Likelihoods
		likelihood_id = paste0("treelikelihood_", partitionName)
		
		# Date ref (the partition) for likelihood
		data_idref_XML = xmlNode(name="data", attrs=list(idref=partitionName) )
		
		# SiteModel ref for likelihood
		siteModel_idref_XML = xmlNode(name="siteModel", attrs=list(idref=siteModel_id) )
		
		# The likelihood XML
		tree_name_idref = paste0("@", tree_name)
		
		# Default: overwritten if SpeciesTreeRelaxedClock
		# If the gene-tree clock is inherited from the SpeciesTree clock,
		# make a new name etc.
		clockModel_name_idref = paste0("@", clockModel_name)
		
		# Default: overwritten if SpeciesTreeRelaxedClock
		likelihood_XML = xmlNode(name="distribution", attrs=list(id=likelihood_id, branchRateModel=clockModel_name_idref, spec="TreeLikelihood", tree=tree_name_idref), .children=list(data_idref_XML, siteModel_idref_XML) )
		
		mean_of_shared_clock_id = paste0(clock_type, "Mean_of_", clockModel_name)
		mean_of_shared_clock_idref = paste0("@", mean_of_shared_clock_id)
		
		
		
		# Gene trees operating under the coalescent model
		genetree_speciesCoal_id = paste0("genetree_speciesCoal_", tree_name)
		genetree_speciesCoal_idref = paste0("@", genetree_speciesCoal_id)
		

		# LogNormal relaxed clock ON THE SPECIES TREE
		# This means the genetree clockrate is inherited
		print(clock_type)
		
		if ((tolower(clock_type) == tolower("SpeciesTreeUCLN")))
			{
			geneTree_clockModel_name_id = paste0("geneTree_", tree_name, "_inherits_relative_rates_from_SpeciesTreeUCLN_clock_", clockModel_name)
			
			# Parameter for the geneTree clockRate
			geneTree_clockModel_rate_id = paste0(clock_type, "_Mean_of_clockRate_for_geneTree_", tree_name)
			geneTree_clockModel_rate_idref = paste0("@", geneTree_clockModel_rate_id)
			
			#
			geneTree_clockModel_rate_XML = xmlNode(name="parameter", 1.0, attrs=list(id=geneTree_clockModel_rate_id, lower="0.0", name="stateNode"))
			geneTree_clockModel_rate_XML_states = cl(geneTree_clockModel_rate_XML_states, geneTree_clockModel_rate_XML)
			
			# Create prior on geneTree clockRate
			geneTree_clockRate_prior_id = paste0("prior_on_", geneTree_clockModel_rate_id)
			geneTree_clockRate_prior_idref = paste0("@", geneTree_clockRate_prior_id)
			clock_df = seqs_df[rownum,]
			geneTree_clockRate_prior_XML = clock_df_row_to_XML_distribution(clock_df=clock_df, tree_name=tree_name, clock_type=clock_type)
			clockModel_priors_XML = cl(clockModel_priors_XML, geneTree_clockRate_prior_XML$priors)
			clockModel_screenLog_XML = cl(clockModel_screenLog_XML, geneTree_clockRate_prior_XML$screenlog)
			clockModel_traceLog_XML = cl(clockModel_traceLog_XML, geneTree_clockRate_prior_XML$tracelog)

			
			# speciesTreeRates idref child
			speciesTreeRates_child_XML = xmlNode(name="speciesTreeRates", attrs=list(idref=clockModel_name_idref))
			
			# branchRateModel child
			branchRateModel_child_XML = xmlNode(name="branchRateModel", attrs=list(id=geneTree_clockModel_name_id, spec="starbeast2.StarBeastClock", clock.rate=geneTree_clockModel_rate_idref, geneTree=genetree_speciesCoal_idref, speciesTreeRates=clockModel_name_idref) )
			
			likelihood_XMLcomment1 = xmlCommentNode(" Likelihood of this data partition under a SpeciesTreeRelaxedClock with a lognormal distribution on branch rates. ")
			likelihood_XMLcomment2 = xmlCommentNode(" ...which is different than e.g. a lognormal prior distribution on the clock rate mean. ")
			likelihood_XMLcomment3 = xmlCommentNode(" (The gene tree rate is inherited from the species tree rates) ")
			
			likelihood_XML = xmlNode(name="distribution", attrs=list(id=likelihood_id, spec="TreeLikelihood", tree=tree_name_idref), .children=list(data_idref_XML, siteModel_idref_XML, branchRateModel_child_XML) )
			likelihood_XML = cl(likelihood_XMLcomment1, likelihood_XMLcomment2, likelihood_XMLcomment3, likelihood_XML)
			}


		# LogNormal relaxed clock ON THE SPECIES TREE
		# This means the genetree clockrate is inherited
		if ((tolower(clock_type) == tolower("SpeciesTreeUCED")))
			{
			geneTree_clockModel_name_id = paste0("geneTree_", tree_name, "_inherits_relative_rates_from_SpeciesTreeUCED_clock_", clockModel_name)
			
			# speciesTreeRates idref child
			speciesTreeRates_child_XML = xmlNode(name="speciesTreeRates", attrs=list(idref=clockModel_name_idref))
			
			# branchRateModel child
			branchRateModel_child_XML = xmlNode(name="branchRateModel", attrs=list(id=geneTree_clockModel_name_id, spec="starbeast2.StarBeastClock", clock.rate=mean_of_shared_clock_idref, geneTree=genetree_speciesCoal_idref, speciesTreeRates=clockModel_name_idref) )
			
			likelihood_XMLcomment1 = xmlCommentNode(" Likelihood of this data partition under a SpeciesTreeRelaxedClock with an exponential distribution on branch rates. ")
			likelihood_XMLcomment2 = xmlCommentNode(" ...which is different than e.g. a exponential prior distribution on the clock rate mean. ")
			likelihood_XMLcomment3 = xmlCommentNode(" (The gene tree rate is inherited from the species tree rates) ")
			
			likelihood_XML = xmlNode(name="distribution", attrs=list(id=likelihood_id, spec="TreeLikelihood", tree=tree_name_idref), .children=list(data_idref_XML, siteModel_idref_XML, branchRateModel_child_XML) )
			likelihood_XML = cl(likelihood_XMLcomment1, likelihood_XMLcomment2, likelihood_XMLcomment3, likelihood_XML)
			}



		if ((tolower(clock_type) == tolower("SpeciesTreeRLC")))
			{
			geneTree_clockModel_name_id = paste0("geneTree_", tree_name, "_inherits_relative_rates_from_SpeciesTreeRLC_clock_", clockModel_name)
			
			# speciesTreeRates idref child
			speciesTreeRates_child_XML = xmlNode(name="speciesTreeRates", attrs=list(idref=clockModel_name_idref))
			
			# branchRateModel child
			branchRateModel_child_XML = xmlNode(name="branchRateModel", attrs=list(id=geneTree_clockModel_name_id, spec="starbeast2.StarBeastClock", clock.rate=mean_of_shared_clock_idref, geneTree=genetree_speciesCoal_idref, speciesTreeRates=clockModel_name_idref) )
			
			likelihood_XMLcomment1 = xmlCommentNode(" Likelihood of this data partition under a SpeciesTreeRelaxedClock with Relaxed Local Clock (RLC) on branch rates. ")
			likelihood_XMLcomment2 = xmlCommentNode(" ...meaning the clock rate shifts according to a Poisson process on the tree. ")
			likelihood_XMLcomment3 = xmlCommentNode(" (The gene tree rate is inherited from the species tree rates) ")
			
			likelihood_XML = xmlNode(name="distribution", attrs=list(id=likelihood_id, spec="TreeLikelihood", tree=tree_name_idref), .children=list(data_idref_XML, siteModel_idref_XML, branchRateModel_child_XML) )
			likelihood_XML = cl(likelihood_XMLcomment1, likelihood_XMLcomment2, likelihood_XMLcomment3, likelihood_XML)
			}

		
		
		
		# Likelihood reference
		likelihood_idref_XML = xmlNode(name="distribution", attrs=list(idref=likelihood_id) )
		
		# TraceLog of the likelihood XML
		likelihood_traceLog_XML = xmlNode(name="log", attrs=list(idref=likelihood_id) )
		
		# Log mutation rates
		log_mutationRate_XML = xmlNode(name="parameter", attrs=list(idref=mutationRate_id, name="log"))
		

		# Log base frequencies for DNA, not AA (but not Jukes-Cantor JC69)
		if ( (datatype == "DNA") && (JC69_TF == FALSE) )
			{
			# Log base frequency
			log_freqParam_XML = xmlNode(name="parameter", attrs=list(idref=starting_frequencies_id, name="log"))
			} else {
			log_freqParam_XML = NULL
			} # END if ( (datatype == "DNA") && (JC69_TF == FALSE) )
		
		
		#################################################
		#################################################
		# Add XMLs to the list
		# 
		# NOTE: *Only* add the siteModel-defining XMLs the
		#       *first* time the siteModel is encountered
		#################################################		
		#################################################
		# STORE siteModel parameters
		# Re-use siteModel, if it has already been done. Otherwise, make a new one!
		if ( (siteModel_name %in% siteModel_names_used) == FALSE)
			{
			siteModel_params_starting_states_XML = cl(siteModel_params_starting_states_XML, siteModel_params_starting_state_XML)
			DNA_AA_siteModels_XML = cl(DNA_AA_siteModels_XML, siteModel_XML)
			siteModel_params_operators_XML = cl(siteModel_params_operators_XML, siteModel_params_operator_XML)
			siteModel_priors_XML = cl(siteModel_priors_XML, siteModel_prior_XML)
			log_siteModels_XML = cl(log_siteModels_XML, log_siteModel_XML)
			
			# Add siteModel_name to the list of siteModel_names already used.
			siteModel_names_used = c(siteModel_names_used, siteModel_name)
			}


		if ( (frequencyModel_name %in% frequencyModel_names_used) == FALSE)
			{
			frequencies_starting_states_XMLs = cl(frequencies_starting_states_XMLs, frequencies_starting_state_XML)			
			frequencies_operators_XML = cl(frequencies_operators_XML, frequencies_operator_XML)
			frequencies_priors_XML = cl(frequencies_priors_XML, frequencies_prior_XML)
			log_freqParams_XML = cl(log_freqParams_XML, log_freqParam_XML)
			estimated_baseFreqs_XMLs = cl(estimated_baseFreqs_XMLs, estimated_baseFreqs_XML)
			
			# Add siteModel_name to the list of siteModel_names already used.
			frequencyModel_names_used = c(frequencyModel_names_used, frequencyModel_name)
			}


		if ( (gammaShapeModel_name %in% gammaShapeModel_names_used) == FALSE)
			{
			gamma_starting_states_XML = cl(gamma_starting_states_XML, gamma_starting_state_XML)
			gammaShape_params_operators_XML = cl(gammaShape_params_operators_XML, gammaShape_params_operator_XML)
			gammaShape_priors_XML = cl(gammaShape_priors_XML, gammaShape_prior_XML)
			gammaShape_logs_XML = cl(gammaShape_logs_XML, gammaShape_log_XML)
			
			# Add siteModel_name to the list of siteModel_names already used.
			gammaShapeModel_names_used = c(gammaShapeModel_names_used, gammaShapeModel_name)
			}
		
		# Only use kappa, if it's an HKY model
		if (seqs_df$model[rownum] == "HKY")
			{
			if ( (kappaModel_name %in% kappaModel_names_used) == FALSE)
				{
				kappa_starting_states_XML = cl(kappa_starting_states_XML, kappa_starting_state_XML)
				kappa_params_operators_XML = cl(kappa_params_operators_XML, kappa_param_operator_XML)
				kappa_priors_XML = cl(kappa_priors_XML, kappa_prior_XML)
				log_kappas_XML = cl(log_kappas_XML, log_kappa_XML)
			
				# Add siteModel_name to the list of siteModel_names already used.
				kappaModel_names_used = c(kappaModel_names_used, kappaModel_name)
				}
			} # END if (seqs_df$model[rownum] == "HKY")
			
			

		
		# STORE likelihoods (for EVERY line, i.e. EVERY data section)
		likelihoods_XML = cl(likelihoods_XML, likelihood_XML)
		likelihoods_idrefs_XML = cl(likelihoods_idrefs_XML, likelihood_idref_XML)
		likelihoods_traceLog_XML = cl(likelihoods_traceLog_XML, likelihood_traceLog_XML)

		
		
		# Clock model stuff
		
		#if ( (clockModel_name %in% clockModel_names_used) == FALSE)
		if (TRUE)	# always do this
			{
			mutationRate_params_starting_states_XML = cl(mutationRate_params_starting_states_XML, mutationRate_params_starting_state_XML)
			log_mutationRates_XML = cl(log_mutationRates_XML, log_mutationRate_XML)
			}
		clockModel_names_used = c(clockModel_names_used, clockModel_name)
		} # END for (i in 1:length(uniq_partitions))
	#######################################################
	# END Loop through each sequence partition
	#######################################################	




		# Put the kappas in the siteModel parameters, if it's an HKY model
		if (seqs_df$model[rownum] == "HKY")
			{
			siteModel_params_starting_states_XML = cl(siteModel_params_starting_states_XML, kappa_starting_states_XML)
			siteModel_params_operators_XML = cl(siteModel_params_operators_XML, kappa_param_operator_XML)
			siteModel_priors_XML = cl(siteModel_priors_XML, kappa_priors_XML)
			log_siteModels_XML = cl(log_siteModels_XML, log_kappas_XML)
			}
		# Add the base frequency priors (if any)
		siteModel_priors_XML = cl(siteModel_priors_XML, frequencies_priors_XML)
	
	# Add in the morphological relRates
	mutationRate_params_starting_states_XML = cl(mutationRate_params_starting_states_XML, morph_mutationRate_params_starting_states_XML)

	# Return
	if (is.null(xml))
		{
		res = NULL
		res$DNA_AA_siteModels_XML = cl(estimated_baseFreqs_XMLs, DNA_AA_siteModels_XML)
		res$siteModel_params_starting_states_XML = siteModel_params_starting_states_XML
		res$gamma_starting_states_XML = gamma_starting_states_XML
		res$frequencies_starting_states_XMLs = frequencies_starting_states_XMLs
		res$mutationRate_params_starting_states_XML = mutationRate_params_starting_states_XML
		res$geneTree_clockModel_rate_XML_states = geneTree_clockModel_rate_XML_states
		
		res$frequencies_operators_XML = frequencies_operators_XML
		res$siteModel_params_operators_XML = cl(siteModel_params_operators_XML, kappa_params_operators_XML)
		#res$relativeMutationRates_operators_XML = relativeMutationRates_operators_XML
		res$gammaShape_params_operators_XML = gammaShape_params_operators_XML
		
		res$siteModel_priors_XML = siteModel_priors_XML
		res$gammaShape_priors_XML = gammaShape_priors_XML
		res$clockModel_priors_XML = clockModel_priors_XML
		
		res$likelihoods_XML = likelihoods_XML
		res$likelihoods_idrefs_XML = likelihoods_idrefs_XML
		res$likelihoods_traceLog_XML = likelihoods_traceLog_XML
		
		res$log_freqParams_XML = log_freqParams_XML
		res$log_siteModels_XML = log_siteModels_XML
		res$log_mutationRates_XML = log_freqParams_XML
		res$gammaShape_logs_XML = gammaShapes_log_XML
		#res$screenLog_mutationRates_XML = screenLog_mutationRates_XML
		
		res$clockModel_screenLog_XML = clockModel_screenLog_XML
		res$clockModel_traceLog_XML = clockModel_traceLog_XML
		
		# Miscellaneous
		res$exponential_distrib_for_gammaShape_priors_XML = exponential_distrib_for_gammaShape_priors_XML

		extract='
		DNA_AA_siteModels_XML = res$DNA_AA_siteModels_XML
		siteModel_params_starting_states_XML = res$siteModel_params_starting_states_XML
		gamma_starting_states_XML = res$gamma_starting_states_XML
		frequencies_starting_states_XMLs = res$frequencies_starting_states_XMLs
		mutationRate_params_starting_states_XML = res$mutationRate_params_starting_states_XML
		
		frequencies_operators_XML = res$frequencies_operators_XML
		siteModel_params_operators_XML = res$siteModel_params_operators_XML
		#relativeMutationRates_operators_XML = res$relativeMutationRates_operators_XML
		gammaShape_params_operators_XML = res$gammaShape_params_operators_XML
		
		siteModel_priors_XML = res$siteModel_priors_XML
		gammaShape_priors_XML = res$gammaShape_priors_XML

		likelihoods_XML = res$likelihoods_XML
		likelihoods_idrefs_XML = res$likelihoods_idrefs_XML
		likelihoods_traceLog_XML = res$likelihoods_traceLog_XML
		
		log_freqParams_XML = res$log_freqParams_XML
		log_siteModels_XML = res$log_siteModels_XML
		log_mutationRates_XML = res$log_freqParams_XML
		gammaShape_logs_XML = res$gammaShape_logs_XML
		#screenLog_mutationRates_XML = res$screenLog_mutationRates_XML
		
		exponential_distrib_for_gammaShape_priors_XML = res$exponential_distrib_for_gammaShape_priors_XML
		'
		
		return(res)
		} else {
		xml$sitemodels = cl(xml$sitemodels, estimated_baseFreqs_XMLs, DNA_AA_siteModels_XML)
		xml$state = cl(xml$state, frequencies_starting_states_XMLs)
		xml$state = cl(xml$state, siteModel_params_starting_states_XML)
		xml$state = cl(xml$state, gamma_starting_states_XML)
		xml$state = cl(xml$state, mutationRate_params_starting_states_XML)
		xml$state = cl(xml$state, geneTree_clockModel_rate_XML_states)
		
		xml$operators = cl(xml$operators, frequencies_operators_XML)
		xml$operators = cl(xml$operators, siteModel_params_operators_XML, kappa_params_operators_XML)
		#xml$operators = cl(xml$operators, relativeMutationRates_operators_XML)
		xml$operators = cl(xml$operators, gammaShape_params_operators_XML)
		
		xml$priors = cl(xml$priors, siteModel_priors_XML)
		xml$priors = cl(xml$priors, gammaShape_priors_XML)
		xml$priors = cl(xml$priors, clockModel_priors_XML)
		
		# Put the likelihood calculation in sitemodels, just after the 
		# actual sitemodels
		xml$sitemodels = cl(xml$sitemodels, likelihoods_XML)
		# Put the likelihood references in likelihoods
		xml$likes = cl(xml$likes, likelihoods_idrefs_XML)
		# And trace the likelihoods
		xml$tracelog = cl(xml$tracelog, likelihoods_traceLog_XML)
			
		xml$tracelog = cl(xml$tracelog, log_freqParams_XML)
		xml$tracelog = cl(xml$tracelog, log_mutationRates_XML)
		xml$tracelog = cl(xml$tracelog, log_siteModels_XML)
		xml$tracelog = cl(xml$tracelog, gammaShape_logs_XML)
		xml$tracelog = cl(xml$tracelog, clockModel_traceLog_XML)
		
		#xml$screenlog = cl(xml$screenlog, screenLog_mutationRates_XML)
		xml$screenlog = cl(xml$screenlog, likelihoods_traceLog_XML, clockModel_screenLog_XML)
		
		xml$misc = cl(xml$misc, exponential_distrib_for_gammaShape_priors_XML)
		
		return(xml)
		} # END if (is.null(xml))
	
	stop("ERROR in parse_DNA_siteModels(): shouldn't get here.")
	} # END parse_DNA_siteModels <- function(seqs_df, xml=NULL)
