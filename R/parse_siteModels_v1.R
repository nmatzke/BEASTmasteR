
parse_DNA_AA_siteModels <- function(seqs_df, xml=NULL, tree_name="shared_tree")
	{
	# Subset to remove user="no" lines
	seqs_df = seqs_df[seqs_df$use != "no", ]

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
	
	# Make filters
	uniq_partitions = unique(seqs_df$filteredAlignmentName)
	uniq_partitions

	# Get the (first!) rownums
	rownums = match(x=uniq_partitions, table=seqs_df$filteredAlignmentName)
	
	
	#######################################################
	# Initialize XML lists for parameters, operators, priors, logs
	#######################################################
	
	# Lists of the parameters for each section
	DNA_AA_siteModels_XML = NULL
	frequencies_starting_states_XML = c(list(bl()), list(xmlCommentNode(" Starting states of base frequencies ")))
	siteModel_params_starting_states_XML = c(list(bl()), list(xmlCommentNode(" Starting states of siteModel parameters ")))
	mutationRate_params_starting_states_XML = c(list(bl()), list(xmlCommentNode(" Starting states of mutationRate parameters ")))
	gamma_starting_states_XML = c(list(bl()), list(xmlCommentNode(" Starting states of gamma parameters ")))
	
	# Operators
	frequencies_operators_XML = c(list(bl()), list(xmlCommentNode(" Operators for base frequencies ")))
	siteModel_params_operators_XML = c(list(bl()), list(xmlCommentNode(" Operators for siteModel parameters ")))
	mutationRate_params_operators_XML = c(list(bl()), list(xmlCommentNode(" Operators for mutationRate parameters ")))
	gammaShape_params_operators_XML = c(list(bl()), list(xmlCommentNode(" Operators for gammaShape parameters ")))
	gammaShape_logs_XML = c(list(bl()), list(xmlCommentNode(" Logs for gammaShape parameters ")))
	
	# Priors
	# Frequency priors seem to be left out of BEAUTi-derived XML, so 
	# leaving them out here.
	#frequency_priors_XML = c(list(bl()), list(xmlCommentNode(" Priors on the base frequency parameters ")) )
	siteModel_priors_XML = c(list(bl()), list(xmlCommentNode(" Priors on the siteModel parameters (e.g. kappa) ")) )
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
	log_mutationRates_XML = c(list(bl()), list(xmlCommentNode(" Trace Log of the mutationRate parameters ")))
	

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
	for (i in 1:length(rownums))
		{
		rownum = rownums[i]
		partitionName = uniq_partitions[i]
		datasetName = seqs_df$datasetName[rownum]
		
		# Actively load the clockmodel_name
		clockModel_name = seqs_df$clockmodel_name[rownum]
		
		datatype = seqs_df$type[rownum]
		
		print(seqs_df)
		print(datatype)
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
		id = paste0("proportionInvariant_", partitionName)
		
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
		id = paste0(seqs_df$baseFreqs[rownum], "Frequencies_", partitionName)
		frequencies = paste0("@freqParameter_", partitionName)
		baseFreqs_idref_XML = xmlNode(name="frequencies", attrs=list(id=id, frequencies=frequencies, spec="Frequencies") )
		
		# Frequencies starting states
		frequencies_id = paste0("freqParameter_", partitionName)
		frequencies_starting_state_XML = xmlNode(name="parameter", 1/dimension, attrs=list(id=frequencies_id, dimension=dimension, lower="0.0", upper="1.0", name="stateNode") )
		
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
			id = paste0("JC69_", partitionName)
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
				id = paste0(AAmodel, "_", partitionName)
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
				id = paste0(AAmodel, "_", rownum, "_", partitionName)
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

		

		if (seqs_df$model[rownum] == "HKY")
			{
			id = paste0("HKY_", partitionName)
			kappaParam_idref = paste0("@kappa_", partitionName)
			substModel_XML = xmlNode(name="substModel", attrs=list(id=id, kappa=kappaParam_idref, spec="HKY"), .children=list(baseFreqs_idref_XML))

			# Substitution Model parameters, starting states
			kappa_id = paste0("kappa_", partitionName)
			siteModel_params_starting_state_XML = xmlNode(name="parameter", 2.0, attrs=list(id=kappa_id, lower="0.0", name="stateNode") )
			
			# Operator on substitution model parameter(s)
			kappa_scaler_id = paste0(kappa_id, "_Scaler")
			param_ref = paste0("@", kappa_id)
			siteModel_params_operator_XML = xmlNode(name="operator", attrs=list(id=kappa_scaler_id, parameter=param_ref, scaleFactor="0.5", spec="ScaleOperator", weight="0.1") )
			
			
			# Make an logNormal distribution, for use in setting up the 
			# priors on kappas
			# Store this in "xml$miscellaneous"
			param1_id = paste0("logNormal_mean_for_kappa_prior_partition_", partitionName)
			param2_id = paste0("logNormal_stddev_for_kappa_prior_partition_", partitionName)
			param1_XML = xmlNode(name="parameter", 1.0, attrs=list(id=param1_id, name="M", estimate="false"))
			param2_XML = xmlNode(name="parameter", 1.25, attrs=list(id=param2_id, name="S", estimate="false"))
			logNormal_distrib_for_kappa_priors_XML_id = paste0("logNormal_distrib_for_kappa_priors_", partitionName)
			logNormal_distrib_for_kappa_priors_XML = xmlNode(name="LogNormal", attrs=list(id=logNormal_distrib_for_kappa_priors_XML_id, name="distr"), .children=list(param1_XML, param2_XML))
			kappa_prior_id = paste0("prior_for_kappa_", partitionName)
			partition_idref = paste0("@", partitionName)
			kappa_prior_XML = xmlNode(name="prior", attrs=list(id=kappa_prior_id, name="distribution", x=kappaParam_idref), .children=list(logNormal_distrib_for_kappa_priors_XML))
			
			# Store the kappa priors in the general siteModel priors
			siteModel_prior_XML = kappa_prior_XML
			
			# Trace Log of siteModel parameter(s)
			log_siteModel_XML = xmlNode(name="parameter", attrs=list(idref=kappa_id, name="log") )
			} # END if (seqs_df$model[rownum] == "HKY")
		
		
		if (seqs_df$model[rownum] == "GTR")
			{
			id = paste0("GTR_", partitionName)
			
			rateAC_id = paste0("rateAC_", partitionName)
			rateAG_id = paste0("rateAG_", partitionName)
			rateAT_id = paste0("rateAT_", partitionName)
			rateCG_id = paste0("rateCG_", partitionName)
			rateCT_id = paste0("rateCT_", partitionName)
			rateGT_id = paste0("rateGT_", partitionName)

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

			siteModel_params_operator_XML_rateAC = xmlNode(name="operator", attrs=list(id=rateAC_scaler_id, parameter=rateAC_idref, scaleFactor="0.5", spec="ScaleOperator", weight="0.1") )
			siteModel_params_operator_XML_rateAG = xmlNode(name="operator", attrs=list(id=rateAG_scaler_id, parameter=rateAG_idref, scaleFactor="0.5", spec="ScaleOperator", weight="0.1") )
			siteModel_params_operator_XML_rateAT = xmlNode(name="operator", attrs=list(id=rateAT_scaler_id, parameter=rateAT_idref, scaleFactor="0.5", spec="ScaleOperator", weight="0.1") )
			siteModel_params_operator_XML_rateCG = xmlNode(name="operator", attrs=list(id=rateCG_scaler_id, parameter=rateCG_idref, scaleFactor="0.5", spec="ScaleOperator", weight="0.1") )
			siteModel_params_operator_XML_rateCT = xmlNode(name="operator", attrs=list(id=rateCT_scaler_id, parameter=rateCT_idref, scaleFactor="0.5", spec="ScaleOperator", weight="0.1") )
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
			param1_id = paste0("alpha_for_rateAC_gammaPrior_partition_", partitionName)
			param2_id = paste0("beta_for_rateAC_gammaPrior_partition_", partitionName)
			param1_XML = xmlNode(name="parameter", 0.05, attrs=list(id=param1_id, name="alpha", estimate="false"))
			param2_XML = xmlNode(name="parameter", 10.0, attrs=list(id=param2_id, name="beta", estimate="false"))
			Gamma_distrib_for_rateAC_priors_XML_id = paste0("Gamma_distrib_for_rateAC_prior_", partitionName)
			Gamma_distrib_for_rateAC_priors_XML_XML = xmlNode(name="Gamma", attrs=list(id=Gamma_distrib_for_rateAC_priors_XML_id, name="distr"), .children=list(param1_XML, param2_XML))
			Gamma_prior_rateAC_id = paste0("prior_distrib_on_rateAC_", partitionName)
			partition_idref = paste0("@", partitionName)
			Gamma_prior_rateAC_XML = xmlNode(name="prior", attrs=list(id=Gamma_prior_rateAC_id, name="distribution", x=rateAC_idref), .children=list(Gamma_distrib_for_rateAC_priors_XML_XML))
			
			# Prior on rateAG
			param1_id = paste0("alpha_for_rateAG_gammaPrior_partition_", partitionName)
			param2_id = paste0("beta_for_rateAG_gammaPrior_partition_", partitionName)
			param1_XML = xmlNode(name="parameter", 0.05, attrs=list(id=param1_id, name="alpha", estimate="false"))
			param2_XML = xmlNode(name="parameter", 10.0, attrs=list(id=param2_id, name="beta", estimate="false"))
			Gamma_distrib_for_rateAG_priors_XML_id = paste0("Gamma_distrib_for_rateAG_prior_", partitionName)
			Gamma_distrib_for_rateAG_priors_XML_XML = xmlNode(name="Gamma", attrs=list(id=Gamma_distrib_for_rateAG_priors_XML_id, name="distr"), .children=list(param1_XML, param2_XML))
			Gamma_prior_rateAG_id = paste0("prior_distrib_on_rateAG_", partitionName)
			partition_idref = paste0("@", partitionName)
			Gamma_prior_rateAG_XML = xmlNode(name="prior", attrs=list(id=Gamma_prior_rateAG_id, name="distribution", x=rateAG_idref), .children=list(Gamma_distrib_for_rateAG_priors_XML_XML))
			
			# Prior on rateAT
			param1_id = paste0("alpha_for_rateAT_gammaPrior_partition_", partitionName)
			param2_id = paste0("beta_for_rateAT_gammaPrior_partition_", partitionName)
			param1_XML = xmlNode(name="parameter", 0.05, attrs=list(id=param1_id, name="alpha", estimate="false"))
			param2_XML = xmlNode(name="parameter", 10.0, attrs=list(id=param2_id, name="beta", estimate="false"))
			Gamma_distrib_for_rateAT_priors_XML_id = paste0("Gamma_distrib_for_rateAT_prior_", partitionName)
			Gamma_distrib_for_rateAT_priors_XML_XML = xmlNode(name="Gamma", attrs=list(id=Gamma_distrib_for_rateAT_priors_XML_id, name="distr"), .children=list(param1_XML, param2_XML))
			Gamma_prior_rateAT_id = paste0("prior_distrib_on_rateAT_", partitionName)
			partition_idref = paste0("@", partitionName)
			Gamma_prior_rateAT_XML = xmlNode(name="prior", attrs=list(id=Gamma_prior_rateAT_id, name="distribution", x=rateAT_idref), .children=list(Gamma_distrib_for_rateAT_priors_XML_XML))
			
			# Prior on rateCG
			param1_id = paste0("alpha_for_rateCG_gammaPrior_partition_", partitionName)
			param2_id = paste0("beta_for_rateCG_gammaPrior_partition_", partitionName)
			param1_XML = xmlNode(name="parameter", 0.05, attrs=list(id=param1_id, name="alpha", estimate="false"))
			param2_XML = xmlNode(name="parameter", 10.0, attrs=list(id=param2_id, name="beta", estimate="false"))
			Gamma_distrib_for_rateCG_priors_XML_id = paste0("Gamma_distrib_for_rateCG_prior_", partitionName)
			Gamma_distrib_for_rateCG_priors_XML_XML = xmlNode(name="Gamma", attrs=list(id=Gamma_distrib_for_rateCG_priors_XML_id, name="distr"), .children=list(param1_XML, param2_XML))
			Gamma_prior_rateCG_id = paste0("prior_distrib_on_rateCG_", partitionName)
			partition_idref = paste0("@", partitionName)
			Gamma_prior_rateCG_XML = xmlNode(name="prior", attrs=list(id=Gamma_prior_rateCG_id, name="distribution", x=rateCG_idref), .children=list(Gamma_distrib_for_rateCG_priors_XML_XML))
			
			# Prior on rateCT
			param1_id = paste0("alpha_for_rateCT_gammaPrior_partition_", partitionName)
			param2_id = paste0("beta_for_rateCT_gammaPrior_partition_", partitionName)
			param1_XML = xmlNode(name="parameter", 0.05, attrs=list(id=param1_id, name="alpha", estimate="false"))
			param2_XML = xmlNode(name="parameter", 10.0, attrs=list(id=param2_id, name="beta", estimate="false"))
			Gamma_distrib_for_rateCT_priors_XML_id = paste0("Gamma_distrib_for_rateCT_prior_", partitionName)
			Gamma_distrib_for_rateCT_priors_XML_XML = xmlNode(name="Gamma", attrs=list(id=Gamma_distrib_for_rateCT_priors_XML_id, name="distr"), .children=list(param1_XML, param2_XML))
			Gamma_prior_rateCT_id = paste0("prior_distrib_on_rateCT_", partitionName)
			partition_idref = paste0("@", partitionName)
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
		siteModel_id = paste0("siteModel_", partitionName)
		gammaNum = seqs_df$gammaNum[rownum]
		mutationRate_id = paste0(partitionName, "_relRate")
		mutationRate_idref = paste0("@", mutationRate_id)
		shapeParam = paste0("@gammaShape_", partitionName)
		siteModel_XML = xmlNode(name="siteModel", attrs=list(id=siteModel_id, mutationRate=mutationRate_idref, shape=shapeParam, gammaCategoryCount=gammaNum, spec="SiteModel"), .children=list(pInv_XML, substModel_XML) )

		# Make the XML for the Gamma starting states
		gammaShape_id = paste0("gammaShape_", partitionName)
		gammaShape_idref = paste0("@", gammaShape_id)
		gamma_starting_state_XML = xmlNode(name="parameter", 1.0, attrs=list(id=gammaShape_id, name="stateNode") )
		
		# Make the XML for the Gamma priors
		gamma_prior_id = paste0(gammaShape_id, "_prior")
		distrib_idref = paste0("@", exponential_distrib_for_gammaShape_priors_XML_id)
		gammaShape_idref = paste0("@", gammaShape_id)
		gammaShape_prior_XML = xmlNode(name="prior", attrs=list(id=gamma_prior_id, distr=distrib_idref, name="distribution", x=gammaShape_idref))
		
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
			# Operator on base frequency parameter(s)
			frequencies_exchanger_id = paste0(frequencies_id, "_Exchanger")
			param_XML = xmlNode(name="parameter", attrs=list(idref=frequencies_id) )
			frequencies_operator_XML = xmlNode(name="operator", attrs=list(id=frequencies_exchanger_id, delta="0.01", spec="DeltaExchangeOperator", weight="0.1"), .children=list(param_XML) )
			} else {
			frequencies_operator_XML = NULL
			}# END if (datatype == "DNA")
		
		# Operator on gammaShape parameter(s)
		if (gammaNum == 1)
			{
			gammaOperatorWeight = 0
			} else {
			gammaOperatorWeight = 0.1
			} # END if (gammaNum == 1)
		gammaShape_operator_id = paste0("gammaShapeScaler_", partitionName)
		gammaShape_params_operator_XML = xmlNode(name="operator", attrs=list(id=gammaShape_operator_id, parameter=gammaShape_idref, scaleFactor="0.5", spec="ScaleOperator", weight=gammaOperatorWeight) )
		
		# Log the gamma shape
		gammaShape_log_XML = xmlNode(name="log", attrs=list(idref=gammaShape_id))
		gammaShape_logs_XML = c(gammaShape_logs_XML, list(gammaShape_log_XML))
		
		# Likelihoods
		likelihood_id = paste0("treelikelihood_", partitionName)
		
		# Date ref (the partition) for likelihood
		data_idref_XML = xmlNode(name="data", attrs=list(idref=partitionName) )
		
		# SiteModel ref for likelihood
		siteModel_idref_XML = xmlNode(name="siteModel", attrs=list(idref=siteModel_id) )
		
		# The likelihood XML
		tree_name_idref = paste0("@", tree_name)
		clockModel_name_idref = paste0("@", clockModel_name)
		likelihood_XML = xmlNode(name="distribution", attrs=list(id=likelihood_id, branchRateModel=clockModel_name_idref, spec="TreeLikelihood", tree=tree_name_idref), .children=list(data_idref_XML, siteModel_idref_XML) )
		
		# Likelihood reference
		likelihood_idref_XML = xmlNode(name="distribution", attrs=list(idref=likelihood_id) )
		
		# TraceLog of the likelihood XML
		likelihood_traceLog_XML = xmlNode(name="log", attrs=list(idref=likelihood_id) )
		
		

		# Log base frequencies for DNA, not AA (but not Jukes-Cantor JC69)
		if ( (datatype == "DNA") && (JC69_TF == FALSE) )
			{
			# Log base frequency
			log_freqParam_XML = xmlNode(name="parameter", attrs=list(idref=frequencies_id, name="log"))
			} else {
			log_freqParam_XML = NULL
			} # END if ( (datatype == "DNA") && (JC69_TF == FALSE) )
		
		# Log mutation rates
		log_mutationRate_XML = xmlNode(name="parameter", attrs=list(idref=mutationRate_id, name="log"))
		
		
		# Add them to list
		
		# STORE siteModel parameters
		DNA_AA_siteModels_XML = cl(DNA_AA_siteModels_XML, siteModel_XML)
		frequencies_starting_states_XML = cl(frequencies_starting_states_XML, frequencies_starting_state_XML)
		siteModel_params_starting_states_XML = cl(siteModel_params_starting_states_XML, siteModel_params_starting_state_XML)
		mutationRate_params_starting_states_XML = cl(mutationRate_params_starting_states_XML, mutationRate_params_starting_state_XML)
		gamma_starting_states_XML = cl(gamma_starting_states_XML, gamma_starting_state_XML)
		
		# STORE operators
		#operator_mutationRates_XML -- note: dealt with below
		frequencies_operators_XML = cl(frequencies_operators_XML, frequencies_operator_XML)
		siteModel_params_operators_XML = cl(siteModel_params_operators_XML, siteModel_params_operator_XML)
		gammaShape_params_operators_XML = cl(gammaShape_params_operators_XML, gammaShape_params_operator_XML)
		
		# STORE priors
		siteModel_priors_XML = cl(siteModel_priors_XML, siteModel_prior_XML)
		gammaShape_priors_XML = cl(gammaShape_priors_XML, gammaShape_prior_XML)
		
		# STORE likelihoods
		likelihoods_XML = cl(likelihoods_XML, likelihood_XML)
		likelihoods_idrefs_XML = cl(likelihoods_idrefs_XML, likelihood_idref_XML)
		likelihoods_traceLog_XML = cl(likelihoods_traceLog_XML, likelihood_traceLog_XML)
		
		# STORE logs
		log_freqParams_XML = cl(log_freqParams_XML, log_freqParam_XML)
		log_mutationRates_XML = cl(log_mutationRates_XML, log_mutationRate_XML)
		log_siteModels_XML = cl(log_siteModels_XML, log_siteModel_XML)
		
		
		# END the loop through partitions
		} # END for (i in 1:length(uniq_partitions))
	

	# Return
	if (is.null(xml))
		{
		res = NULL
		res$DNA_AA_siteModels_XML = DNA_AA_siteModels_XML
		res$siteModel_params_starting_states_XML = siteModel_params_starting_states_XML
		res$gamma_starting_states_XML = gamma_starting_states_XML
		res$frequencies_starting_states_XML = frequencies_starting_states_XML
		res$mutationRate_params_starting_states_XML = mutationRate_params_starting_states_XML

		res$frequencies_operators_XML = frequencies_operators_XML
		res$siteModel_params_operators_XML = siteModel_params_operators_XML
		#res$relativeMutationRates_operators_XML = relativeMutationRates_operators_XML
		res$gammaShape_params_operators_XML = gammaShape_params_operators_XML
		
		res$siteModel_priors_XML = siteModel_priors_XML
		res$gammaShape_priors_XML = gammaShape_priors_XML
		
		res$likelihoods_XML = likelihoods_XML
		res$likelihoods_idrefs_XML = likelihoods_idrefs_XML
		res$likelihoods_traceLog_XML = likelihoods_traceLog_XML
		
		res$log_freqParams_XML = log_freqParams_XML
		res$log_siteModels_XML = log_siteModels_XML
		res$log_mutationRates_XML = log_freqParams_XML
		res$gammaShape_logs_XML = gammaShapes_log_XML
		#res$screenLog_mutationRates_XML = screenLog_mutationRates_XML
		
		# Miscellaneous
		res$exponential_distrib_for_gammaShape_priors_XML = exponential_distrib_for_gammaShape_priors_XML

		extract='
		DNA_AA_siteModels_XML = res$DNA_AA_siteModels_XML
		siteModel_params_starting_states_XML = res$siteModel_params_starting_states_XML
		gamma_starting_states_XML = res$gamma_starting_states_XML
		frequencies_starting_states_XML = res$frequencies_starting_states_XML
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
		xml$sitemodels = c(xml$sitemodels, DNA_AA_siteModels_XML)
		xml$state = c(xml$state, frequencies_starting_states_XML)
		xml$state = c(xml$state, siteModel_params_starting_states_XML)
		xml$state = c(xml$state, gamma_starting_states_XML)
		xml$state = c(xml$state, mutationRate_params_starting_states_XML)
		
		xml$operators = c(xml$operators, frequencies_operators_XML)
		xml$operators = c(xml$operators, siteModel_params_operators_XML)
		#xml$operators = c(xml$operators, relativeMutationRates_operators_XML)
		xml$operators = c(xml$operators, gammaShape_params_operators_XML)
		
		xml$priors = c(xml$priors, siteModel_priors_XML)
		xml$priors = c(xml$priors, gammaShape_priors_XML)
		
		# Put the likelihood calculation in sitemodels, just after the 
		# actual sitemodels
		xml$sitemodels = c(xml$sitemodels, likelihoods_XML)
		# Put the likelihood references in likelihoods
		xml$likes = c(xml$likes, likelihoods_idrefs_XML)
		# And trace the likelihoods
		xml$tracelog = c(xml$tracelog, likelihoods_traceLog_XML)
			
		xml$tracelog = c(xml$tracelog, log_freqParams_XML)
		xml$tracelog = c(xml$tracelog, log_mutationRates_XML)
		xml$tracelog = c(xml$tracelog, log_siteModels_XML)
		xml$tracelog = c(xml$tracelog, gammaShape_logs_XML)
		
		#xml$screenlog = c(xml$screenlog, screenLog_mutationRates_XML)
		
		xml$misc = c(xml$misc, exponential_distrib_for_gammaShape_priors_XML)
		
		return(xml)
		} # END if (is.null(xml))
	
	stop("ERROR in parse_DNA_siteModels(): shouldn't get here.")
	} # END parse_DNA_siteModels <- function(seqs_df, xml=NULL)
