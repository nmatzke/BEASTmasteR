
clock_df_row_to_XML_distribution <- function(clock_df, tree_name="shared_tree")
	{
	defaults='
	tree_name="shared_tree"
	'
	
	# Get the clock type
	clock_type = clock_df$clockmodel_type	# e.g. "ucld"

	clock_types_allowed = c("ucld", "uced", "rlc")
	if (clock_type %in% clock_types_allowed == FALSE)
		{
		txt = paste0("STOP ERROR in clock_df_row_to_XML_distribution(): clock_type=", clock_type, " but no match found in clock_types_allowed.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		cat("Printing 'clock_types_allowed':\n\n")
		print(clock_types_allowed)
		cat("\n\n")
		stop(txt)
		} # END if (clock_type %in% clock_types_allowed == FALSE)

	# LogNormal relaxed clock
	if (clock_type == "ucld")
		{
		tmpXML = clock_df_row_to_XML_distribution_meanSD(clock_df, tree_name="shared_tree")
		}
	# Exponential relaxed clock
	if (clock_type == "uced")
		{
		tmpXML = clock_df_row_to_XML_distribution_justMean(clock_df, tree_name="shared_tree")
		}
	# Local random clock
	if (clock_type == "rlc")
		{
		tmpXML = clock_df_row_to_XML_distribution_justMean(clock_df, tree_name="shared_tree")
		}
	
	return(tmpXML)
	}


clock_df_row_to_XML_distribution_meanSD <- function(clock_df, tree_name="shared_tree")
	{
	defaults='
	clock_df = seqs_df[1,]
	tree_name="shared_tree"
	'
	
	# Get the clock name
	clockModel_name = clock_df$clockmodel_name
	clockModel_idref = paste0("@", clockModel_name)
		
	# Get the tree name
	tree_idref = paste0("@", tree_name)
	
	# Get the ids for relaxed clockrate mean and SD
	clock_type = clock_df$clockmodel_type	# e.g. "ucld"

	if ( (clock_type == "uced") || (clock_type == "rlc") )
		{
		txt = paste0("STOP ERROR in 'clock_df_row_to_XML_distribution_justMean': your clock_df$clockmodel_type=='uced' or 'rlc', but these models need _justMean, not _meanSD")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		}
	

	mean_of_shared_clock_id = paste0(clock_type, "Mean_of_", clockModel_name) 
	mean_of_shared_clock_idref = paste0("@", mean_of_shared_clock_id)
	SD_of_shared_clock_id = paste0(clock_type, "Stdev_of_", clockModel_name) 
	SD_of_shared_clock_idref = paste0("@", SD_of_shared_clock_id)
	
	
	#######################################################
	# Prior on the mean clockrate
	#######################################################
	
	# Get the offset, if it exists
	if (isblank_TF(clock_df$clockrate_offset) == FALSE)
		{
		offset = clock_df$clockrate_offset
		} else {
		offset = 0
		} # END if (isblank_TF(clock_df$offset) == FALSE)

	# Get the meanInRealSpace, if it exists
	if (( isblank_TF(clock_df$clockrate_meanInRealSpace) == FALSE) && (clock_df$clockrate_meanInRealSpace == "no") )
		{
		meanInRealSpace = FALSE
		meanInRealSpace_txt = "false"
		} else {
		meanInRealSpace = TRUE
		meanInRealSpace_txt = "true"
		}
	
	if (clock_df$clockrate_prior_dist == "normal")
		{
		id_of_prior_distrib_on_clock_mean = paste0("prP_NormalDistrib_on_mean_of_", clockModel_name)
		
		meanval = clock_df$clockrate_prior_param1 + offset
		sdval = clock_df$clockrate_prior_param2
		
		distribution_name = paste0("NormalDistrib_on_mean_of_", clockModel_name)
		param1_name = paste0("mean_of_NormalDistrib_on_mean_of_", clockModel_name)
		param2_name = paste0("sd_of_NormalDistrib_on_mean_of_", clockModel_name)
		param1_XML = xmlNode(name="parameter", clock_df$clockrate_prior_param1, attrs=list(id=param1_name, name="mean", estimate="false"))
		param2_XML = xmlNode(name="parameter", clock_df$clockrate_prior_param2, attrs=list(id=param2_name, name="sigma", estimate="false"))
		distrib_XML = xmlNode(name="Normal", attrs=list(id=distribution_name, name="distr", offset="0"), .children=list(param1_XML, param2_XML))
		
		txt = paste0(" Prior probability density of the mean of '", clockModel_name, "', according to a ", clock_df$clockrate_prior_dist, " distribution with mean=", meanval, ", sd=", sdval, ". ")
		} # END if (clock_df$clockrate_prior_dist == "normal")

	if (clock_df$clockrate_prior_dist == "lognormal")
		{
		id_of_prior_distrib_on_clock_mean = paste0("prP_LogNormalDistrib_on_mean_of_", clockModel_name)

		meanval = clock_df$clockrate_prior_param1
		sdval = clock_df$clockrate_prior_param2

		# Convert to real space, if needed
		if (meanInRealSpace == FALSE)
			{
			meanval = exp(meanval)
			sdval = exp(meanval)
			}
		 
		distribution_name = paste0("LogNormalDistrib_on_mean_of_", clockModel_name)
		param1_name = paste0("mean_of_LogNormalDistrib_on_mean_of_", clockModel_name)
		param2_name = paste0("sd_of_LogNormalDistrib_on_mean_of_", clockModel_name)
		param1_XML = xmlNode(name="parameter", clock_df$clockrate_prior_param1, attrs=list(id=param1_name, name="M", estimate="false"))
		param2_XML = xmlNode(name="parameter", clock_df$clockrate_prior_param2, attrs=list(id=param2_name, name="S", estimate="false"))
		distrib_XML = xmlNode(name="LogNormal", attrs=list(id=distribution_name, name="distr", offset=offset, meanInRealSpace="true"), .children=list(param1_XML, param2_XML))

		txt = paste0(" Prior probability density of the mean of '", clockModel_name, "', according to a ", clock_df$clockrate_prior_dist, " distribution with mean=", meanval, ", sd=", sdval, ", offset=", offset, ". ")
		} # END if (clock_df$clockrate_prior_dist == "normal")

	if (clock_df$clockrate_prior_dist == "exponential")
		{
		id_of_prior_distrib_on_clock_mean = paste0("prP_ExponentialDistrib_on_mean_of_", clockModel_name)

		meanval = clock_df$clockrate_prior_param1

		# Convert to real space, if needed
		if (meanInRealSpace == FALSE)
			{
			meanval = 1/(meanval)
			}
		 
		distribution_name = paste0("ExponentialDistrib_on_mean_of_", clockModel_name)
		param1_name = paste0("mean_of_ExponentialDistrib_on_mean_of_", clockModel_name)
		param1_XML = xmlNode(name="parameter", clock_df$clockrate_prior_param1, attrs=list(id=param1_name, name="mean", estimate="false"))
		distrib_XML = xmlNode(name="Exponential", attrs=list(id=distribution_name, name="distr", offset=offset), .children=list(param1_XML))

		txt = paste0(" Prior probability density of the mean of '", clockModel_name, "', according to a ", clock_df$clockrate_prior_dist, " distribution with mean=", meanval, ", offset=", offset, ". ")
		} # END if (clock_df$clockrate_prior_dist == "normal")

	if (clock_df$clockrate_prior_dist == "uniform")
		{
		id_of_prior_distrib_on_clock_mean = paste0("prP_UniformDistrib_on_mean_of_", clockModel_name)

		lower = clock_df$clockrate_prior_param1 + offset
		upper = clock_df$clockrate_prior_param2 + offset

		distribution_name = paste0("UniformDistrib_on_mean_of_", clockModel_name)
		param1_name = paste0("lower_of_UniformDistrib_on_mean_of_", clockModel_name)
		param2_name = paste0("upper_of_UniformDistrib_on_mean_of_", clockModel_name)
		distrib_XML = xmlNode(name="Uniform", attrs=list(id=distribution_name, name="distr", offset="0", lower=lower, upper=upper))

		txt = paste0(" Prior probability density of the mean of '", clockModel_name, "', according to a ", clock_df$clockrate_prior_dist, " distribution with lower=", lower, ", upper=", upper, ". ")
		} # END if (clock_df$clockrate_prior_dist == "uniform")
	
	
	########################################################################
	# Now, make the prior probability distribution for the mean of the clock model
	########################################################################	
	# Overall Prior Distribution
	prior_XML = xmlNode(name="prior", attrs=list(id=id_of_prior_distrib_on_clock_mean, name="distribution", x=mean_of_shared_clock_idref), .children=list(distrib_XML))
	prior_XML
	comment_XML = xmlCommentNode(txt)
	comment_XML
	
	prior_XML_plus_comment = list(bl(), comment_XML, prior_XML)






	#######################################################
	# Prior on the SD of relaxed clock
	#######################################################
	# Get the offset, if it exists
	if (isblank_TF(clock_df$clockSD_offset) == FALSE)
		{
		offset = clock_df$clockSD_offset
		} else {
		offset = 0
		} # END if (isblank_TF(clock_df$offset) == FALSE)

	# Get the meanInRealSpace, if it exists
	if (( isblank_TF(clock_df$clockSD_meanInRealSpace) == FALSE) && (clock_df$clockSD_meanInRealSpace == "no") )
		{
		meanInRealSpace = FALSE
		meanInRealSpace_txt = "false"
		} else {
		meanInRealSpace = TRUE
		meanInRealSpace_txt = "true"
		}
	
	if (clock_df$clockSD_prior_dist == "normal")
		{
		id_of_prior_distrib_on_clock_SD = paste0("prP_NormalDistrib_on_SD_of_", clockModel_name)
		
		meanval = clock_df$clockSD_prior_param1 + offset
		sdval = clock_df$clockSD_prior_param2
		
		distribution_name = paste0("NormalDistrib_on_SD_of_", clockModel_name)
		param1_name = paste0("mean_of_NormalDistrib_on_SD_of_", clockModel_name)
		param2_name = paste0("sd_of_NormalDistrib_on_SD_of_", clockModel_name)
		param1_XML = xmlNode(name="parameter", clock_df$clockSD_prior_param1, attrs=list(id=param1_name, name="mean", estimate="false"))
		param2_XML = xmlNode(name="parameter", clock_df$clockSD_prior_param2, attrs=list(id=param2_name, name="sigma", estimate="false"))
		distrib_XML = xmlNode(name="Normal", attrs=list(id=distribution_name, name="distr", offset="0"), .children=list(param1_XML, param2_XML))
		
		txt = paste0(" Prior probability density of the SD of '", clockModel_name, "', according to a ", clock_df$clockSD_prior_dist, " distribution with mean=", meanval, ", sd=", sdval, ". ")
		} # END if (clock_df$clockSD_prior_dist == "normal")

	if (clock_df$clockSD_prior_dist == "lognormal")
		{
		id_of_prior_distrib_on_clock_SD = paste0("prP_LogNormalDistrib_on_SD_of_", clockModel_name)

		meanval = clock_df$clockSD_prior_param1
		sdval = clock_df$clockSD_prior_param2

		# Convert to real space, if needed
		if (meanInRealSpace == FALSE)
			{
			meanval = exp(meanval)
			sdval = exp(meanval)
			}
		 
		distribution_name = paste0("LogNormalDistrib_on_SD_of_", clockModel_name)
		param1_name = paste0("mean_of_LogNormalDistrib_on_SD_of_", clockModel_name)
		param2_name = paste0("sd_of_LogNormalDistrib_on_SD_of_", clockModel_name)
		param1_XML = xmlNode(name="parameter", clock_df$clockSD_prior_param1, attrs=list(id=param1_name, name="M", estimate="false"))
		param2_XML = xmlNode(name="parameter", clock_df$clockSD_prior_param2, attrs=list(id=param2_name, name="S", estimate="false"))
		distrib_XML = xmlNode(name="LogNormal", attrs=list(id=distribution_name, name="distr", offset=offset, meanInRealSpace="true"), .children=list(param1_XML, param2_XML))

		txt = paste0(" Prior probability density of the SD of '", clockModel_name, "', according to a ", clock_df$clockSD_prior_dist, " distribution with mean=", meanval, ", sd=", sdval, ", offset=", offset, ". ")
		} # END if (clock_df$clockSD_prior_dist == "normal")

	if (clock_df$clockSD_prior_dist == "exponential")
		{
		id_of_prior_distrib_on_clock_SD = paste0("prP_ExponentialDistrib_on_SD_of_", clockModel_name)

		meanval = clock_df$clockSD_prior_param1

		# Convert to real space, if needed
		if (meanInRealSpace == FALSE)
			{
			meanval = 1/(meanval)
			}
		 
		distribution_name = paste0("ExponentialDistrib_on_SD_of_", clockModel_name)
		param1_name = paste0("mean_of_ExponentialDistrib_on_SD_of_", clockModel_name)
		param1_XML = xmlNode(name="parameter", clock_df$clockSD_prior_param1, attrs=list(id=param1_name, name="mean", estimate="false"))
		distrib_XML = xmlNode(name="Exponential", attrs=list(id=distribution_name, name="distr", offset=offset), .children=list(param1_XML))

		txt = paste0(" Prior probability density of the SD of '", clockModel_name, "', according to a ", clock_df$clockSD_prior_dist, " distribution with mean=", meanval, ", offset=", offset, ". ")
		} # END if (clock_df$clockSD_prior_dist == "normal")

	if (clock_df$clockSD_prior_dist == "uniform")
		{
		id_of_prior_distrib_on_clock_SD = paste0("prP_UniformDistrib_on_SD_of_", clockModel_name)

		lower = clock_df$clockSD_prior_param1 + offset
		upper = clock_df$clockSD_prior_param2 + offset

		distribution_name = paste0("UniformDistrib_on_SD_of_", clockModel_name)
		param1_name = paste0("lower_of_UniformDistrib_on_SD_of_", clockModel_name)
		param2_name = paste0("upper_of_UniformDistrib_on_SD_of_", clockModel_name)
		distrib_XML = xmlNode(name="Uniform", attrs=list(id=distribution_name, name="distr", offset="0", lower=lower, upper=upper))

		txt = paste0(" Prior probability density of the SD of '", clockModel_name, "', according to a ", clock_df$clockSD_prior_dist, " distribution with lower=", lower, ", upper=", upper, ". ")
		} # END if (clock_df$clockSD_prior_dist == "uniform")
	
	
	########################################################################
	# Now, make the prior probability distribution for the SD of the clock model
	########################################################################	
	# Overall Prior Distribution
	prior_SD_XML = xmlNode(name="prior", attrs=list(id=id_of_prior_distrib_on_clock_SD, name="distribution", x=SD_of_shared_clock_idref), .children=list(distrib_XML))
	prior_SD_XML
	comment_SD_XML = xmlCommentNode(txt)
	comment_SD_XML
	
	prior_SD_XML_plus_comment = list(bl(), comment_SD_XML, prior_SD_XML)







	########################################################################
	# Log the prior density of the clockrate mean
	########################################################################
	log_XMLtxt = xmlCommentNode(" Log of the prior probability density of the sampled mean of this clock ")
	log_XML = xmlNode(name="log", attrs=list(idref=id_of_prior_distrib_on_clock_mean))


	########################################################################
	# Log the prior density of the clockrate SD
	########################################################################
	log_SD_XMLtxt = xmlCommentNode(" Log of the prior probability density of the sampled SD of this clock ")
	log_SD_XML = xmlNode(name="log", attrs=list(idref=id_of_prior_distrib_on_clock_SD))

	
	#######################################################
	# Log on clock model parameters
	#######################################################
	screenLog_xmlComment = xmlCommentNode(paste0(" log of the clock mean for model '", clockModel_name, "' "))
	clockMean_log_XML = xmlNode(name="parameter", attrs=list(idref=mean_of_shared_clock_id, name="log") )
	clockStdev_log_XML = xmlNode(name="parameter", attrs=list(idref=SD_of_shared_clock_id, name="log") )
	
	traceLog_xmlComment = xmlCommentNode(paste0(" Trace log of the '", clockModel_name, "'clock model "))

	# Log the "rate summary"
	rate_summary_id = paste0("rate_summary_", clockModel_name)
	branchratemodel_idref = paste0("@", clockModel_name)
	tree_name_idref = paste0("@", tree_name)
	rate_summary_xmlComment1 = xmlCommentNode(paste0(" Log of rate statistics: rate.mean, rate.variance, rate.coefficientOfVariation "))
	rate_summary_xmlComment2 = xmlCommentNode(" coefficient of variation = stdev / mean = sqrt(variance)/mean ")
	rate_summary_xml = xmlNode(name="log", attrs=list(id=rate_summary_id, branchratemodel=branchratemodel_idref, tree=tree_name_idref, spec="beast.evolution.branchratemodel.RateStatistic"))
	
	########################################################################
	# Bring all the logs together
	########################################################################
	clock_log_XMLs = list(bl(), traceLog_xmlComment, clockMean_log_XML, clockStdev_log_XML, bl(), rate_summary_xmlComment1, rate_summary_xmlComment2, rate_summary_xml, bl(), log_XMLtxt, log_XML, log_SD_XMLtxt, log_SD_XML)

	
	########################################################################
	# Operators on clockrate and clockSD
	########################################################################
	# Scale the clock rate (ucldMean)
	ucldMeanScaler_id = paste0(mean_of_shared_clock_id, "_scaler")
	mean_of_shared_clock_idref = paste0("@", mean_of_shared_clock_id)
	ucldMeanScaler_XML = xmlNode(name="operator", attrs=list(id=ucldMeanScaler_id, parameter=mean_of_shared_clock_idref, scaleFactor="0.5", spec="ScaleOperator", weight="3.0"))
	
	# Scale the clock stdev (ucldStdev)
	ucldStdevScaler_id = paste0(SD_of_shared_clock_id, "_scaler")
	SD_of_shared_clock_idref = paste0("@", SD_of_shared_clock_id)
	ucldStdevScaler_XML = xmlNode(name="operator", attrs=list(id=ucldStdevScaler_id, parameter=SD_of_shared_clock_idref, scaleFactor="0.5", spec="ScaleOperator", weight="3.0"))
	
	# Move clockrate up and tree size down at the same time (improves mixing)
	relaxedUpDownOperator_ucldMean_id = paste0(mean_of_shared_clock_id, "_relaxedUpDownOperator")
	moveup_XML = xmlNode(name="parameter", attrs=list(idref=mean_of_shared_clock_id, name="up") )
	movedown_XML = xmlNode(name="tree", attrs=list(idref=tree_name, name="down") )
	
	relaxedUpDownOperator_ucldMean_XML = xmlNode(name="operator", attrs=list(id=relaxedUpDownOperator_ucldMean_id, scaleFactor="0.75", spec="UpDownOperator", weight="3.0"), .children=list(moveup_XML, movedown_XML) )

	
	tmpXML = NULL
	tmpXML$priors = c(prior_XML_plus_comment, prior_SD_XML_plus_comment)
	tmpXML$operators = list(bl(), xmlCommentNode(" Operators on the clock mean and SD "), ucldMeanScaler_XML, ucldStdevScaler_XML, relaxedUpDownOperator_ucldMean_XML)
	tmpXML$screenlog = list(bl(), screenLog_xmlComment, clockMean_log_XML)
	tmpXML$tracelog = clock_log_XMLs
	return(tmpXML)
	} # END clock_df_row_to_XML_distribution_meanSD()






clock_df_row_to_XML_distribution_justMean <- function(clock_df, tree_name="shared_tree")
	{
	defaults='
	tree_name="shared_tree"
	'
	
	# Get the clock name
	clockModel_name = clock_df$clockmodel_name
	clockModel_idref = paste0("@", clockModel_name)
		
	# Get the tree name
	tree_idref = paste0("@", tree_name)
	
	# Get the ids for relaxed clockrate mean and SD
	clock_type = clock_df$clockmodel_type	# e.g. "uced" or "rlc"
	
	if (clock_type == "ucld")
		{
		txt = paste0("STOP ERROR in 'clock_df_row_to_XML_distribution_justMean': your clock_df$clockmodel_type=='ucld', but the ucld model needs _meanSD, not _justMean")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		}
	
	mean_of_shared_clock_id = paste0(clock_type, "Mean_of_", clockModel_name) 
	mean_of_shared_clock_idref = paste0("@", mean_of_shared_clock_id)
	
	
	#######################################################
	# Prior on the mean clockrate
	#######################################################
	
	# Get the offset, if it exists
	if (isblank_TF(clock_df$clockrate_offset) == FALSE)
		{
		offset = clock_df$clockrate_offset
		} else {
		offset = 0
		} # END if (isblank_TF(clock_df$offset) == FALSE)

	# Get the meanInRealSpace, if it exists
	if (( isblank_TF(clock_df$clockrate_meanInRealSpace) == FALSE) && (clock_df$clockrate_meanInRealSpace == "no") )
		{
		meanInRealSpace = FALSE
		meanInRealSpace_txt = "false"
		} else {
		meanInRealSpace = TRUE
		meanInRealSpace_txt = "true"
		}
	
	if (clock_df$clockrate_prior_dist == "normal")
		{
		id_of_prior_distrib_on_clock_mean = paste0("prP_NormalDistrib_on_mean_of_", clockModel_name)
		
		meanval = clock_df$clockrate_prior_param1 + offset
		sdval = clock_df$clockrate_prior_param2
		
		distribution_name = paste0("NormalDistrib_on_mean_of_", clockModel_name)
		param1_name = paste0("mean_of_NormalDistrib_on_mean_of_", clockModel_name)
		param2_name = paste0("sd_of_NormalDistrib_on_mean_of_", clockModel_name)
		param1_XML = xmlNode(name="parameter", clock_df$clockrate_prior_param1, attrs=list(id=param1_name, name="mean", estimate="false"))
		param2_XML = xmlNode(name="parameter", clock_df$clockrate_prior_param2, attrs=list(id=param2_name, name="sigma", estimate="false"))
		distrib_XML = xmlNode(name="Normal", attrs=list(id=distribution_name, name="distr", offset="0"), .children=list(param1_XML, param2_XML))
		
		txt = paste0(" Prior probability density of the mean of '", clockModel_name, "', according to a ", clock_df$clockrate_prior_dist, " distribution with mean=", meanval, ", sd=", sdval, ". ")
		} # END if (clock_df$clockrate_prior_dist == "normal")

	if (clock_df$clockrate_prior_dist == "lognormal")
		{
		id_of_prior_distrib_on_clock_mean = paste0("prP_LogNormalDistrib_on_mean_of_", clockModel_name)

		meanval = clock_df$clockrate_prior_param1
		sdval = clock_df$clockrate_prior_param2

		# Convert to real space, if needed
		if (meanInRealSpace == FALSE)
			{
			meanval = exp(meanval)
			sdval = exp(meanval)
			}
		 
		distribution_name = paste0("LogNormalDistrib_on_mean_of_", clockModel_name)
		param1_name = paste0("mean_of_LogNormalDistrib_on_mean_of_", clockModel_name)
		param2_name = paste0("sd_of_LogNormalDistrib_on_mean_of_", clockModel_name)
		param1_XML = xmlNode(name="parameter", clock_df$clockrate_prior_param1, attrs=list(id=param1_name, name="M", estimate="false"))
		param2_XML = xmlNode(name="parameter", clock_df$clockrate_prior_param2, attrs=list(id=param2_name, name="S", estimate="false"))
		distrib_XML = xmlNode(name="LogNormal", attrs=list(id=distribution_name, name="distr", offset=offset, meanInRealSpace="true"), .children=list(param1_XML, param2_XML))

		txt = paste0(" Prior probability density of the mean of '", clockModel_name, "', according to a ", clock_df$clockrate_prior_dist, " distribution with mean=", meanval, ", sd=", sdval, ", offset=", offset, ". ")
		} # END if (clock_df$clockrate_prior_dist == "normal")

	if (clock_df$clockrate_prior_dist == "exponential")
		{
		id_of_prior_distrib_on_clock_mean = paste0("prP_ExponentialDistrib_on_mean_of_", clockModel_name)

		meanval = clock_df$clockrate_prior_param1

		# Convert to real space, if needed
		if (meanInRealSpace == FALSE)
			{
			meanval = 1/(meanval)
			}
		 
		distribution_name = paste0("ExponentialDistrib_on_mean_of_", clockModel_name)
		param1_name = paste0("mean_of_ExponentialDistrib_on_mean_of_", clockModel_name)
		param1_XML = xmlNode(name="parameter", clock_df$clockrate_prior_param1, attrs=list(id=param1_name, name="mean", estimate="false"))
		distrib_XML = xmlNode(name="Exponential", attrs=list(id=distribution_name, name="distr", offset=offset), .children=list(param1_XML))

		txt = paste0(" Prior probability density of the mean of '", clockModel_name, "', according to a ", clock_df$clockrate_prior_dist, " distribution with mean=", meanval, ", offset=", offset, ". ")
		} # END if (clock_df$clockrate_prior_dist == "normal")

	if (clock_df$clockrate_prior_dist == "uniform")
		{
		id_of_prior_distrib_on_clock_mean = paste0("prP_UniformDistrib_on_mean_of_", clockModel_name)

		lower = clock_df$clockrate_prior_param1 + offset
		upper = clock_df$clockrate_prior_param2 + offset

		distribution_name = paste0("UniformDistrib_on_mean_of_", clockModel_name)
		param1_name = paste0("lower_of_UniformDistrib_on_mean_of_", clockModel_name)
		param2_name = paste0("upper_of_UniformDistrib_on_mean_of_", clockModel_name)
		distrib_XML = xmlNode(name="Uniform", attrs=list(id=distribution_name, name="distr", offset="0", lower=lower, upper=upper))

		txt = paste0(" Prior probability density of the mean of '", clockModel_name, "', according to a ", clock_df$clockrate_prior_dist, " distribution with lower=", lower, ", upper=", upper, ". ")
		} # END if (clock_df$clockrate_prior_dist == "uniform")
	
	
	########################################################################
	# Now, make the prior probability distribution for the mean of the clock model
	########################################################################	
	# Overall Prior Distribution
	prior_XML = xmlNode(name="prior", attrs=list(id=id_of_prior_distrib_on_clock_mean, name="distribution", x=mean_of_shared_clock_idref), .children=list(distrib_XML))
	prior_XML
	comment_XML = xmlCommentNode(txt)
	comment_XML
	
	prior_XML_plus_comment = list(bl(), comment_XML, prior_XML)



	########################################################################
	# Log the prior density of the clockrate mean
	########################################################################
	log_XMLtxt = xmlCommentNode(" Log of the prior probability density of the sampled mean of this clock ")
	log_XML = xmlNode(name="log", attrs=list(idref=id_of_prior_distrib_on_clock_mean))


	
	#######################################################
	# Log on clock model parameters
	#######################################################
	screenLog_xmlComment = xmlCommentNode(paste0(" log of the clock mean for model '", clockModel_name, "' "))
	clockMean_log_XML = xmlNode(name="parameter", attrs=list(idref=mean_of_shared_clock_id, name="log") )
	
	traceLog_xmlComment = xmlCommentNode(paste0(" Trace log of the '", clockModel_name, "'clock model "))

	# Log the "rate summary"
	rate_summary_id = paste0("rate_summary_", clockModel_name)
	branchratemodel_idref = paste0("@", clockModel_name)
	tree_name_idref = paste0("@", tree_name)
	rate_summary_xmlComment1 = xmlCommentNode(paste0(" Log of rate statistics: rate.mean, rate.variance, rate.coefficientOfVariation "))
	rate_summary_xmlComment2 = xmlCommentNode(" coefficient of variation = stdev / mean = sqrt(variance)/mean ")
	rate_summary_xml = xmlNode(name="log", attrs=list(id=rate_summary_id, branchratemodel=branchratemodel_idref, tree=tree_name_idref, spec="beast.evolution.branchratemodel.RateStatistic"))
	
	########################################################################
	# Bring all the logs together
	########################################################################
	clock_log_XMLs = list(bl(), traceLog_xmlComment, clockMean_log_XML, bl(), rate_summary_xmlComment1, rate_summary_xmlComment2, rate_summary_xml, bl(), log_XMLtxt, log_XML)

	
	########################################################################
	# Operators on clockrate
	########################################################################
	# Scale the clock rate (ucldMean, ucedMean, or rlcMean)
	ucldMeanScaler_id = paste0(mean_of_shared_clock_id, "_scaler")
	mean_of_shared_clock_idref = paste0("@", mean_of_shared_clock_id)
	ucldMeanScaler_XML = xmlNode(name="operator", attrs=list(id=ucldMeanScaler_id, parameter=mean_of_shared_clock_idref, scaleFactor="0.5", spec="ScaleOperator", weight="3.0"))
	
	# Move clockrate up and tree size down at the same time (improves mixing)
	relaxedUpDownOperator_ucldMean_id = paste0(mean_of_shared_clock_id, "_relaxedUpDownOperator")
	moveup_XML = xmlNode(name="parameter", attrs=list(idref=mean_of_shared_clock_id, name="up") )
	movedown_XML = xmlNode(name="tree", attrs=list(idref=tree_name, name="down") )
	
	relaxedUpDownOperator_ucldMean_XML = xmlNode(name="operator", attrs=list(id=relaxedUpDownOperator_ucldMean_id, scaleFactor="0.75", spec="UpDownOperator", weight="3.0"), .children=list(moveup_XML, movedown_XML) )

	
	tmpXML = NULL
	tmpXML$priors = c(prior_XML_plus_comment)
	tmpXML$operators = list(bl(), xmlCommentNode(" Operators on the clock mean "), ucldMeanScaler_XML, relaxedUpDownOperator_ucldMean_XML)
	tmpXML$screenlog = list(bl(), screenLog_xmlComment, clockMean_log_XML)
	tmpXML$tracelog = clock_log_XMLs
	return(tmpXML)
	} # END clock_df_row_to_XML_distribution_justMean()



