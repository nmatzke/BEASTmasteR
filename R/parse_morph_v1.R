#######################################################
# Parsing morphology data from NEXUS to Beast2 XML
#######################################################


add_parens_to_multistate_chars <- function(char_dtf)
	{
	defaults='
	char_dtf = t(morph_df2_corrected)
	'
	
	# Which characters are multistate?
	TF = nchar(char_dtf) > 1

	if (sum(TF) > 0)
		{
		for (i in 1:sum(TF))
			{
			charval = char_dtf[TF][i]
			words = strsplit(charval, split="")[[1]]
			txt1 = paste(words, sep="", collapse=" ")
			txt = paste0("(", txt1, ")")
			char_dtf[TF][i] = txt
			} # END for (i in 1:sum(TF))
		} # END if (sum(TF) > 0)

	return(char_dtf)
	} # END add_parens_to_multistate_chars <- function(char_dtf)



#######################################################
# Writing morphology models for different numbers
# of characters, ordered vs. unordered, etc.
#######################################################
make_Beast2_morph_models <- function(morphList, morphRate_name="morph_relRate", morphGamma_name="morph_gammaShape", clockModel_name="shared_clock", tree_name="shared_tree", xml=NULL)
	{
	defaults='
	morphRate_name="morph_relRate"
	morphGamma_name="morph_gammaShape"
	clockModel_name="shared_clock"
	tree_name="shared_tree"
	xml=NULL
	'
	
	# Calculate the data section names
	data_section_name = morphList$name_of_userDataType
data_section_name = unname(mapply(FUN=gsub, pattern="datatype_", x=data_section_name, MoreArgs=list(replacement="")))

	# Set up lists
	morph_rateMatrices = list(bl(), bl(), xmlCommentNode(" CUSTOM MORPHOLOGY RATE MATRICES FOR DISCRETE MORPHOLOGICAL CHARACTER EVOLUTION "), bl() )
	morph_substModels = list(bl(), bl(), xmlCommentNode(" CUSTOM MORPHOLOGY SUBSTITUTION MODELS FOR DISCRETE MORPHOLOGICAL CHARACTER EVOLUTION "), bl() )
	morph_sitemodels = list(bl(), bl(), xmlCommentNode(" CUSTOM MORPHOLOGY SITE MODELS FOR DISCRETE MORPHOLOGICAL CHARACTER EVOLUTION "), bl() )
	morph_priors = list(bl(), xmlCommentNode(" Prior on parameter(s) for morphology evolution model ") )
	morph_likelihoods = list(bl(), bl(), xmlCommentNode(" Likelihood of morphology data given tree, model, morphology model parameters ") )
	morph_likelihoods_calc = list(bl(), xmlCommentNode(" Likelihood of morphology data given tree, model, morphology model parameters "))
	morph_likelihoods_log = list(bl(), xmlCommentNode(" Likelihood of morphology data given tree, model, morphology model parameters "))
	
	
	morph_operators = NULL
	morph_param_states = NULL
	morph_param_logs = NULL

	# States & starting values for mutationRate and gammaShape
	morphRate_state_XML = xmlNode(name="parameter", 1.0, attrs=list(id=morphRate_name, name="stateNode") )
	morphGamma_state_XML = xmlNode(name="parameter", 1.0, attrs=list(id=morphGamma_name, name="stateNode") )
	morph_States_XMLs = list(bl(), xmlCommentNode(" Starting states for morphology model parameters "), morphRate_state_XML, morphGamma_state_XML)
	
	# Prior on gammaShape
	morph_gammaShape_prior_id = paste0(morphGamma_name, "_prior")
	morph_gammaShape_idref = paste0("@", morphGamma_name)
	morph_gammaShape_prior_XML = xmlNode(name="prior", attrs=list(id=morph_gammaShape_prior_id, distr="@exponential_distrib_for_gammaShape_priors", name="distribution", x=morph_gammaShape_idref) )
	morph_gammaShape_prior_XMLs = list(bl(), xmlCommentNode(" Prior density of the gamma shape parameter for morphology "), morph_gammaShape_prior_XML)
	
	# Log gammaShape and prior on gammaShape
	morph_gammaShape_param_log_XML = xmlNode(name="log", attrs=list(idref=morphGamma_name))
	morph_gammaShape_param_log_XMLs = list(bl(), xmlCommentNode(" Log of the gamma shape parameter for morphology "), morph_gammaShape_param_log_XML)
	morph_gammaShape_prior_log_XML = xmlNode(name="log", attrs=list(idref=morph_gammaShape_prior_id))
	morph_gammaShape_prior_log_XMLs = list(bl(), xmlCommentNode(" Log of prior density of the gamma shape parameter for morphology "), morph_gammaShape_prior_log_XML)
	morph_param_logs = c(morph_gammaShape_param_log_XMLs, morph_gammaShape_prior_log_XMLs)
	
	
	# Operator on gammaShape
	morph_gammaShape_operator_id = paste0(morphGamma_name, "_Scaler")
	morph_gammaShape_operator_XML = xmlNode(name="operator", attrs=list(id=morph_gammaShape_operator_id, parameter=morph_gammaShape_idref, scaleFactor="0.5", spec="ScaleOperator", weight="0.1") )
	morph_gammaShape_operator_XMLs = list(bl(), xmlCommentNode(" Operator for scaling the gamma shape parameter for morphology "), morph_gammaShape_operator_XML)
	
	
	# Log on gammaShape
	morph_gammaShape_log_XML = xmlNode(name="parameter", attrs=list(id=morphGamma_name, name="log") )
	morph_gammaShape_log_XMLs = list(bl(), xmlCommentNode(" Log the gamma shape parameter for morphology "), morph_gammaShape_log_XML)
	
	# siteModels for the morphology dataset sections
	# (section by the number of character states, and ordered/unordered)
	
	print("morphList:")
	print(morphList)
	
	for (i in 1:nrow(morphList))
		{
		# "Mutation" rates (rates of discrete character change, really)
		# (all derived from the single morphological relativeRate, but 
		#  expanded to make the rate matrix for a particular dataset)
		numstates = morphList$numstates[i]
		ordering = morphList$ordering[i]
		numGammaCat = morphList$numGammaCat[i]
		
		# Set up the site model for this morphology section
		# Rate matrix
		numrates = (numstates^2) - numstates
		comment_txt = paste0(" rates for a ", numstates, "x", numstates, ", ", ordering, " character matrix; ", numrates, " off-diagonal rates ")
		if (ordering == "unordered")
			{
			ratematrix_txt = paste(rep("1.0", times=numrates), sep=" ", collapse=" ")
			} else {
			tmpmat = matrix(data="0.0", nrow=numstates, ncol=numstates)
			for (ii in 1:nrow(tmpmat))
				{
				for (jj in 1:ncol(tmpmat))
					{
					if (ii == jj)
						{
						tmpmat[[ii,jj]] = NA
						}
					if (abs(ii-jj) == 1)
						{
						tmpmat[[ii,jj]] = "1.0"
						}
					}
				} # END for (ii in 1:nrow(tmpmat))
			rate_vals = c(tmpmat)
			rate_vals = rate_vals[is.na(rate_vals) == FALSE]
			ratematrix_txt = paste(rate_vals, sep=" ", collapse=" ")
			ratematrix_txt
			}
		rate_matrix_id = paste0(data_section_name[i], "_ratematrix")
		rate_matrix_idref = paste0("@", rate_matrix_id)
		rate_matrix_XML = xmlNode(name="parameter", ratematrix_txt, attrs=list(id=rate_matrix_id, dimension=numrates, name="rates") )
		rate_matrix_XMLs = list(bl(), xmlCommentNode(comment_txt), rate_matrix_XML)
		rate_matrix_XMLs
		
		# Equal base frequencies for morphology characters
		# (as 1/0 labeling is arbitrary)
		morph_basefreqs_id = paste0(data_section_name[i], "_basefreqs")
		morph_basefreqs_idref = paste0("@", morph_basefreqs_id)
		freqs_txt = paste(rep(1/numstates, times=numstates), sep=" ", collapse=" ")
		morph_basefreqs_XML = xmlNode(name="frequencies", attrs=list(id=morph_basefreqs_id, frequencies=freqs_txt, spec="Frequencies") )
		
		# Substitution model for this morphology section
		comment_txt = paste0(" substitution model for a ", numstates, "x", numstates, ", ", ordering, " character matrix; ", numrates, " off-diagonal rates ")
		substModel_id = paste0(data_section_name[i], "_substModel")
		substModel_idref = paste0("@", substModel_id)
		substModel_XML = xmlNode(name="substModel", attrs=list(id=substModel_id, rates=rate_matrix_idref, spec="GeneralSubstitutionModel"), .children=list(morph_basefreqs_XML) )
		substModel_XMLs = list(bl(), xmlCommentNode(comment_txt), substModel_XML)
		
		# Site model for this morphology section
		comment_txt1 = paste0(" siteModel for a ", numstates, "x", numstates, ", ", ordering, " character matrix; ", numrates, " off-diagonal rates; ")
		comment_txt2 = paste0(" gamma-distributed among-site rate variation with 4 categories, no invariant sites ")
		siteModel_id = paste0(data_section_name[i], "_sitemodel")
		siteModel_idref = paste0("@", siteModel_id)
		child1 = xmlNode(name="parameter", attrs=list(idref=morphRate_name, name="mutationRate") )
		child2 = xmlNode(name="parameter", attrs=list(idref=morphGamma_name, name="shape") )
		
		pInv_id = paste0(data_section_name[i], "_pInv")
		pInv_XML = xmlNode(name="parameter", 0.0, attrs=list(id=pInv_id, lower="0.0", upper="1.0", name="proportionInvariant") )
		
		child4 = xmlNode(name="substModel", attrs=list(idref=substModel_id) )
		
		siteModel_XML = xmlNode(name="siteModel", attrs=list(id=siteModel_id, gammaCategoryCount=numGammaCat, spec="SiteModel"), .children=list(child1, child2, pInv_XML, child4) )
		siteModel_XMLs = list(bl(), xmlCommentNode(comment_txt1), xmlCommentNode(comment_txt2), siteModel_XML)
		
		
		# Set up likelihood for this morphology section
		txt = paste0(" Likelihood of ", data_section_name[i], " on the tree ")
		morph_likelihood_id = paste0(data_section_name[i], "_treeLikelihood")
		data_section_idref = paste0("@", data_section_name[i])
		clockModel_name_idref = paste0("@", clockModel_name)
		tree_name_idref = paste0("@", tree_name)
		morph_likelihood_XML = xmlNode(name="distribution", attrs=list(id=morph_likelihood_id, spec="TreeLikelihood", data=data_section_idref, tree=tree_name_idref, siteModel=siteModel_idref, branchRateModel=clockModel_name_idref, useAmbiguities="true") )
		morph_likelihood_XMLs = list(bl(), xmlCommentNode(txt), morph_likelihood_XML)

		# Reference the likelihoods in the Likelihood calculation
		morph_likelihood_calc_XML = xmlNode(name="distribution", attrs=list(idref=morph_likelihood_id) )
		morph_likelihood_calc_XMLs = list(bl(), xmlCommentNode(txt), morph_likelihood_calc_XML)
		
		# Log the likelihoods
		morph_likelihood_log_XML = xmlNode(name="log", attrs=list(idref=morph_likelihood_id) )


		
		# Store these
		# Add to siteModel
		morph_rateMatrices = c(morph_rateMatrices, rate_matrix_XMLs)
		morph_substModels = c(morph_substModels, substModel_XMLs)
		morph_sitemodels = c(morph_sitemodels, siteModel_XMLs)
		morph_likelihoods = c(morph_likelihoods, morph_likelihood_XMLs)
		morph_likelihoods_calc = c(morph_likelihoods_calc, morph_likelihood_calc_XMLs)
		morph_likelihoods_log = c(morph_likelihoods_log, list(morph_likelihood_log_XML))
		
		} # END for (i in 1:nrow(morphList))
	
	# Add just the one morphological prior (gamma shape)
	morph_priors = c(morph_priors, morph_gammaShape_prior_XMLs)
	
	if (is.null(xml))
		{
		res = NULL
		res$morph_rateMatrices = morph_rateMatrices
		res$morph_substModels = morph_substModels
		res$morph_sitemodels = morph_sitemodels
		res$morph_priors = morph_priors
		res$morph_likelihoods = morph_likelihoods
		res$morph_likelihoods_calc = morph_likelihoods_calc
		res$morph_likelihoods_log = morph_likelihoods_log
		res$morph_operators = morph_gammaShape_operator_XMLs
		res$morph_param_states = morph_States_XMLs
		res$morph_param_logs = morph_param_logs
		
		extract='
		morph_rateMatrices = res$morph_rateMatrices
		morph_substModels = res$morph_substModels
		morph_sitemodels = res$morph_sitemodels
		morph_priors = res$morph_priors
		morph_likelihoods = res$morph_likelihoods
		morph_likelihoods_calc = res$morph_likelihoods_calc
		morph_likelihoods_log = res$morph_likelihoods_log
		morph_operators = res$morph_gammaShape_operator_XMLs
		morph_param_states = res$morph_States_XMLs
		morph_param_logs = res$morph_param_logs
		'
		
		return(res)
		
		} else {
		xml$sitemodels = c(xml$sitemodels, morph_rateMatrices)
		xml$sitemodels = c(xml$sitemodels, morph_substModels)
		xml$sitemodels = c(xml$sitemodels, morph_sitemodels)
		xml$priors = c(xml$priors, morph_priors)
		xml$sitemodels = c(xml$sitemodels, morph_likelihoods)
		xml$likes = c(xml$likes, morph_likelihoods_calc)
		xml$operators = c(xml$operators, morph_gammaShape_operator_XMLs)
		xml$state = c(xml$state, morph_States_XMLs)

		xml$tracelog = c(xml$tracelog, morph_likelihoods_log)
		xml$screenlog = c(xml$screenlog, morph_likelihoods_log)
		xml$tracelog = c(xml$tracelog, morph_param_logs)

		return(xml)
		
		} # END if (is.null(xml))
	
	cat("\n\n")
	stop("\n\nERROR in make_Beast2_morph_models(): shouldn't get here")
	} # END make_Beast2_morph_models


























