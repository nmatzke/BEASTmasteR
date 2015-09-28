# Input a single line of the Excel file...
# stem=TRUE only specifies that we desire a prior on the stem of a terminal branch;
# don't use it elsewhere!  Stem vs. node calibrations for nodes are 
make_generic_XML_prior <- function(dfline, colname_prefix="birthRate", param_name="birthRate", header_scheme=1, stem=FALSE, distrib=NULL, param1=NULL, param2=NULL, tmp_offset=NULL, meanInRealSpace=NULL)
	{
	defaults='
	param_name="birthRate"
	header_scheme=1
	distrib=NULL
	param1=NULL
	param2=NULL
	tmp_offset=NULL
	meanInRealSpace=NULL
	'
	
	# Refer to the parameter
	param_name_idref = paste0("@", param_name)
	
	# Get names of each header
	if (header_scheme == 1)	# E.g. for birthRate_prior
		{
		dist_hdr = paste0(colname_prefix, "_prior_dist")
		param1_hdr = paste0(colname_prefix, "_prior_param1")
		param2_hdr = paste0(colname_prefix, "_prior_param2")
		offset_hdr = paste0(colname_prefix, "_offset")
		meanInRealSpace_hdr = paste0(colname_prefix, "_meanInRealSpace")
		}

	if (header_scheme == 2)	# More generic
		{
		dist_hdr = paste0(colname_prefix, "_dist")
		param1_hdr = paste0(colname_prefix, "_param1")
		param2_hdr = paste0(colname_prefix, "_param2")
		offset_hdr = paste0(colname_prefix, "_offset")
		meanInRealSpace_hdr = paste0(colname_prefix, "_meanInRealSpace")
		}


	# Get values describing the distribution, if they are not already specified
	if (is.null(distrib))
		{
		text_to_run = paste0("(distrib = dfline$", dist_hdr, ")")
		eval(parse(text=text_to_run))
		} # END if (!is.null(distribution))
	
	if (is.null(param1))
		{
		text_to_run = paste0("param1 = dfline$", param1_hdr)
		eval(parse(text=text_to_run))
		} # END if (!is.null(param1))
	
	if (is.null(param2))
		{
		text_to_run = paste0("param2 = dfline$", param2_hdr)
		eval(parse(text=text_to_run))
		} # END if (!is.null(param2))
	
	if (is.null(tmp_offset))
		{
		text_to_run = paste0("tmp_offset = dfline$", offset_hdr)
		eval(parse(text=text_to_run))
		} # END if (!is.null(tmp_offset))
	
	if (is.null(meanInRealSpace))
		{
		text_to_run = paste0("meanInRealSpace = dfline$", meanInRealSpace_hdr)
		eval(parse(text=text_to_run))
		} # END if (!is.null(meanInRealSpace))
	

	#######################################################
	# Prior on the value of the parameter
	#######################################################
	
	# Get the tmp_offset, if it exists
	if (isblank_TF(tmp_offset) == FALSE)
		{
		tmp_offset = tmp_offset
		} else {
		tmp_offset = 0
		} # END if (isblank_TF(tmp_offset) == FALSE)

	# Get the meanInRealSpace, if it exists
	if (( isblank_TF(meanInRealSpace) == FALSE) && (meanInRealSpace == "no") )
		{
		meanInRealSpace = FALSE
		meanInRealSpace_txt = "false"
		} else {
		meanInRealSpace = TRUE
		meanInRealSpace_txt = "true"
		}
	
	if (distrib == "normal")
		{
		id_of_distribution_on_param = paste0("prP_NormalDistrib_on_param_", param_name)
		
		meanval = param1 + tmp_offset
		sdval = param2
		
		distribution_name = paste0("NormalDistrib_on_param_", param_name)
		param1_name = paste0("mean_of_NormalDistrib_on_param_", param_name)
		param2_name = paste0("sd_of_NormalDistrib_on_param_", param_name)
		param1_XML = xmlNode(name="parameter", param1, attrs=list(id=param1_name, name="mean", estimate="false"))
		param2_XML = xmlNode(name="parameter", param2, attrs=list(id=param2_name, name="sigma", estimate="false"))
		distrib_XML = xmlNode(name="Normal", attrs=list(id=distribution_name, name="distr", offset="0"), .children=list(param1_XML, param2_XML))
		
		txt = paste0(" Prior probability density on the value of the parameter '", param_name, "', according to a ", distrib, " distribution with mean=", meanval, ", sd=", sdval, ". ")
		} # END if (distribution == "normal")

	if (distrib == "lognormal")
		{
		id_of_distribution_on_param = paste0("prP_LogNormalDistrib_on_param_", param_name)

		meanval = param1
		sdval = param2

		# Convert to real space, if needed
		if (meanInRealSpace == FALSE)
			{
			meanval = exp(meanval)
			sdval = exp(meanval)
			}
		 
		distribution_name = paste0("LogNormalDistrib_on_param_", param_name)
		param1_name = paste0("mean_of_LogNormalDistrib_on_param_", param_name)
		param2_name = paste0("sd_of_LogNormalDistrib_on_param_", param_name)
		param1_XML = xmlNode(name="parameter", param1, attrs=list(id=param1_name, name="M", estimate="false"))
		param2_XML = xmlNode(name="parameter", param2, attrs=list(id=param2_name, name="S", estimate="false"))
		distrib_XML = xmlNode(name="LogNormal", attrs=list(id=distribution_name, name="distr", offset=tmp_offset, meanInRealSpace="true"), .children=list(param1_XML, param2_XML))

		txt = paste0(" Prior probability density on the value of the parameter '", param_name, "', according to a ", distrib, " distribution with mean=", meanval, ", sd=", sdval, ", offset=", tmp_offset, ". ")
		} # END if (distribution == "normal")

	if (distrib == "exponential")
		{
		id_of_distribution_on_param = paste0("prP_ExponentialDistrib_on_param_", param_name)

		meanval = param1

		# Convert to real space, if needed
		if (meanInRealSpace == FALSE)
			{
			meanval = 1/(meanval)
			}
		 
		distribution_name = paste0("ExponentialDistrib_on_param_", param_name)
		param1_name = paste0("mean_of_ExponentialDistrib_on_param_", param_name)
		param1_XML = xmlNode(name="parameter", param1, attrs=list(id=param1_name, name="mean", estimate="false"))
		distrib_XML = xmlNode(name="Exponential", attrs=list(id=distribution_name, name="distr", offset=tmp_offset), .children=list(param1_XML))

		txt = paste0(" Prior probability density on the value of the parameter '", param_name, "', according to a ", distrib, " distribution with mean=", meanval, ", offset=", tmp_offset, ". ")
		} # END if (distribution == "normal")

	if (distrib == "uniform")
		{
		id_of_distribution_on_param = paste0("prP_UniformDistrib_on_param_", param_name)

		lower = param1 + tmp_offset
		upper = param2 + tmp_offset

		distribution_name = paste0("UniformDistrib_on_param_", param_name)
		param1_name = paste0("lower_of_UniformDistrib_on_param_", param_name)
		param2_name = paste0("upper_of_UniformDistrib_on_param_", param_name)
		distrib_XML = xmlNode(name="Uniform", attrs=list(id=distribution_name, name="distr", offset="0", lower=lower, upper=upper))

		txt = paste0(" Prior probability density on the value of the parameter '", param_name, "', according to a ", distrib, " distribution with lower=", lower, ", upper=", upper, ". ")
		} # END if (distribution == "uniform")
	



	
	########################################################################
	# Now, make the prior probability distribution for the value of the parameter
	########################################################################	
	# Overall Prior Distribution
	prior_XML = xmlNode(name="prior", attrs=list(id=id_of_distribution_on_param, name="distribution", x=param_name_idref), .children=list(distrib_XML))
	prior_XML
	comment_XML = xmlCommentNode(txt)
	comment_XML
	
	prior_XML_plus_comment = list(bl(), comment_XML, prior_XML)






	########################################################################
	# Log the prior density of the value of the parameter
	########################################################################
	idref_of_distribution_on_param = paste0(id_of_distribution_on_param)
	log_XMLtxt = xmlCommentNode(paste0(" Log of the prior probability density of value of the parameter '", param_name, "' "))
	log_XML = xmlNode(name="log", attrs=list(idref=idref_of_distribution_on_param))
	log_XMLs = list(bl(), log_XMLtxt, log_XML)


	#######################################################
	# Log on clock model parameters
	#######################################################
	#screenLog_xmlComment = xmlCommentNode(paste0(" log of the value of the parameter '", param_name, "' "))
	#param_log_XML = xmlNode(name="parameter", attrs=list(idref=param_name, name="log") )
	
	
	
	#######################################################
	# Output the XML
	#######################################################
	tmpXML = NULL
	tmpXML$priors = c(prior_XML_plus_comment)
	tmpXML$tracelog = log_XMLs
	
	return(tmpXML)
	} # END make_generic_XML_prior <- function(dfline, param_name="birthRate", header_scheme=1, distribution=NULL, param1=NULL, param2=NULL, tmp_offset=NULL, meanInRealSpace=NULL)


# Make a relative node age statistic
# (e.g., a difference between two node ages)
# (or, hopefully, tip ages also)
make_generic_difference_statistic <- function(nodes_df=NULL, OTUs_df=NULL, xml=NULL)
	{
	defaults='
	nodes_df = readWorksheetFromFile(file=xlsfn, sheet="nodes", startRow=15, startCol=1, header=TRUE)
	OTUs_df = readWorksheetFromFile(file=xlsfn, sheet="OTUs", startRow=15, startCol=1, header=TRUE)
	
	# Subset to just the ones we are using
	nodes_df = nodes_df[nodes_df$use != "no",]
	OTUs_df = OTUs_df[OTUs_df$use != "no",]
	'
	
	# Error check
	if ( is.null(nodes_df) && is.null(OTUs_df) )
		{
		stop("STOP ERROR in make_generic_difference_statistic(): at least nodes_df or OTUs_df must be non-null")
		} # END if ( is.null(nodes_df) && is.null(OTUs_df) )
	
	
	# Do the first submitted table df (e.g. nodes_df)
	# Are there any relative priors in nodes?
	groupnames1b = NULL
	defined_nodenames = NULL
	if (!is.null(nodes_df))
		{
		TF = isblank_TF(nodes_df$rel_prior_Name) == FALSE
		reldf = nodes_df[TF,]
		# Subset by "rel_prior_use"
		TF = reldf$rel_prior_use != "no"
		reldf = reldf[TF,]
		
		# Get OTUs/clade names
		groupnames1 = reldf$rel_prior_groupnames[isblank_TF(reldf$rel_prior_groupnames) == FALSE]
		groupnames1a = paste(groupnames1, sep=",", collapse=",")
		groupnames1b = gdata::trim(strsplit(groupnames1a, split=",")[[1]])
		groupnames1b
		
		# Check if these nodes were defined in nodes_df
		defined_nodenames = nodes_df$Taxon[nodes_df$use != "no"]
		}
	
	# Are there any relative priors in OTUs?
	defined_OTUs = NULL
	groupnames2b = NULL
	if (!is.null(OTUs_df))
		{
		# Are there any relative priors in tips?
		TF = isblank_TF(OTUs_df$rel_prior_Name) == FALSE
		reldf2 = OTUs_df[TF,]
		# Subset by "rel_prior_use"
		TF = reldf2$rel_prior_use != "no"
		reldf2 = reldf2[TF,]
		
		# Get OTUs/clade names
		groupnames2 = reldf2$rel_prior_groupnames[isblank_TF(reldf2$rel_prior_groupnames) == FALSE]
		groupnames2a = paste(groupnames2, sep=",", collapse=",")
		groupnames2b = gdata::trim(strsplit(groupnames2a, split=",")[[1]])
		groupnames2b

		# Check if these nodes were defined in nodes_df
		defined_OTUs = OTUs_df$OTUs[OTUs_df$use != "no"]
		} # END if (!is.null(OTUs_df))
	
	# Check to see that all groupnames have been defined
	col_headers = c("rel_prior_Name", "rel_prior_use", "rel_prior_groupnames", "Equation", "rel_prior_distribution", "rel_prior_param1", "rel_prior_param2", "rel_meanInRealSpace", "rel_offset")
	defined_taxa = c(defined_nodenames, defined_OTUs)
	groupnames = c(groupnames1b, groupnames2b)
	
	# Exit if nothing found
	if (length(groupnames) < 1)
		{
		if (is.null(xml))
			{
			res = NULL
			return(res)
			} else {
			return(xml)
			} # END if (is.null(xml))
		} # END if (length(groupnames) < 1)
	
	groupnames = unique(groupnames)
	TF = groupnames %in% defined_taxa
	
	if (sum(TF) != length(groupnames))
		{
		errortxt = paste0("\n\nSTOP ERROR in make_generic_difference_statistic(): all taxa listed in 'rel_prior_groupnames' in the OTUs and/or nodes worksheets must be defined in the 'OTUs' column of the OTUs worksheet, or the 'Taxon' column of the nodes worksheet (and 'use' must not be 'no').\n\n")
		cat(errortxt)
		
		missing_groupnames = groupnames[TF]
		missing_groupnames_txt = paste(missing_groupnames, sep="", collapse=", ")
		
		cat("These are the rel_prior_groupnames not found:\n\n")
		cat(missing_groupnames_txt)
		cat("\n\n")
		stop(errortxt)
		}
	
	# Get the prior names for each tip or node
	# Prior names for OTUs
	priornames_for_groupnames1 = NULL
	fixed_TF1 = NULL
	fixed_ages1 = NULL
	rel_prior_specs_for_groupnames1 = NULL
	for (i in 1:length(groupnames1b))
		{
		# Find the OTU row
		groupname = groupnames1b[i]
		TF = nodes_df$Taxon == groupname
		nodeline_df = nodes_df[TF,]
		
		priorname_for_groupnames1_tmp = paste0(groupname, "_prior")
		priornames_for_groupnames1 = c(priornames_for_groupnames1, priorname_for_groupnames1_tmp)
		
		fixed_TF1 = c(fixed_TF1, FALSE)
		fixed_ages1 = c(fixed_ages1, NA)
		} # END for (i in 1:length(groupnames1b))

	# Extract the relative distribution specifications
	tmp_rel_prior_specs_for_groupnames1 = nodes_df[,col_headers]
	tmp_rel_prior_specs_for_groupnames1 = tmp_rel_prior_specs_for_groupnames1[isblank_TF(tmp_rel_prior_specs_for_groupnames1$Equation) == FALSE,]
	rel_prior_specs_for_groupnames1 = rbind(rel_prior_specs_for_groupnames1, tmp_rel_prior_specs_for_groupnames1)

	priornames_for_groupnames2 = NULL
	fixed_TF2 = NULL
	fixed_ages2 = NULL
	rel_prior_specs_for_groupnames2 = NULL
	for (i in 1:length(groupnames2b))
		{
		# Find the OTU row
		groupname = groupnames2b[i]
		TF = OTUs_df$OTUs == groupname
		OTUline_df = OTUs_df[TF,]
		
		# Capitalize distributions
		distrib = OTUline_df$distribution
		distrib = gsub(pattern="uniform", replacement="Uniform", x=distrib)
		distrib = gsub(pattern="lognormal", replacement="LogNormal", x=distrib)
		distrib = gsub(pattern="normal", replacement="Normal", x=distrib)
		distrib = gsub(pattern="exponential", replacement="Exponential", x=distrib)

		priorname_for_groupnames2_tmp = paste0("prP_", distrib, "Distrib_", groupname)
		priornames_for_groupnames2 = c(priornames_for_groupnames2, priorname_for_groupnames2_tmp)
		
		# If it is fixed, we will just use the number rather than the prior
		if (distrib == "fixed")
			{
			fixed_TF2 = c(fixed_TF2, TRUE)
			fixed_age = OTUline_df$tipdate
			fixed_ages2 = c(fixed_ages2, fixed_age)
			} else {
			fixed_TF2 = c(fixed_TF2, FALSE)
			fixed_age = NA
			fixed_ages2 = c(fixed_ages2, fixed_age)
			}				

		} # END for (i in 1:length(groupnames2b))

	# Extract the relative distribution specifications
	tmp_rel_prior_specs_for_groupnames2 = OTUs_df[,col_headers]
	tmp_rel_prior_specs_for_groupnames2 = tmp_rel_prior_specs_for_groupnames2[isblank_TF(tmp_rel_prior_specs_for_groupnames2$Equation) == FALSE,]
	rel_prior_specs_for_groupnames2 = rbind(rel_prior_specs_for_groupnames2, tmp_rel_prior_specs_for_groupnames2)

	# Make a list of the groupnames, and whether or not they have fixed dates
	priornames_for_groupnames = c(priornames_for_groupnames1, priornames_for_groupnames2)
	fixed_TF = c(fixed_TF1, fixed_TF2)
	fixed_ages = c(fixed_ages1, fixed_ages2)
	groupnames
	
	# Make a data.frame of the rows with relative distributions (via "Equation" column)
	rel_prior_specs_for_groupnames = rbind(rel_prior_specs_for_groupnames1, rel_prior_specs_for_groupnames2)
	rel_prior_specs_for_groupnames
	
	# OK, now write the equations...
	expressions_XMLs = list(bl(), xmlCommentNode(" Calculating equation(s): "))
	relpriors_XMLs = list(bl(), xmlCommentNode(" Relative priors: "))
	logs_of_expressions_XMLs = list(bl(), xmlCommentNode(" Log of calculated equation result(s): "))
	logs_of_relpriors_XMLs = list(bl(), xmlCommentNode(" Log of relative priors: "))
	screenlogs_of_expressions_XMLs = list(bl(), xmlCommentNode(" Log of calculated equation result(s): "))
	
	
	# Go through each equation and write it to XML, along with the priors and logs!
	for (rownum in 1:nrow(rel_prior_specs_for_groupnames))
		{
		specs = rel_prior_specs_for_groupnames[rownum,]
		equation_to_modify = specs$Equation
		rel_prior_Name = specs$rel_prior_Name
		
		if (equation_to_modify == "diff_from_mean")
			{
			# We need a different procedure for this equation
			# Get the groupnames
			tmp_groupnames_txt = specs$rel_prior_groupnames
			tmp_groupnames = gdata::trim(strsplit(tmp_groupnames_txt, split=",")[[1]])
			num_groups = length(tmp_groupnames)
			
			# Which groupnames match
			TF = rep(FALSE, length(groupnames))
			for (g in 1:length(tmp_groupnames))
				{
				if (tmp_groupnames[g] %in% groupnames)
					{
					TF[g] = TRUE
					} # END if (tmp_groupnames[g] %in% groupnames)
				} # for (g in 1:length(tmp_groupnames))
			groupnames_to_use = groupnames[TF]
			priornames_for_groupnames_to_use = priornames_for_groupnames[TF]
			fixed_TF_to_use = fixed_TF[TF]
			fixed_ages_to_use = fixed_ages[TF]
			
			# Make equation for the mean
			priornames_for_groupnames_to_use2 = paste(priornames_for_groupnames_to_use, "[1]", sep="")
			# If any are fixed, just use the numbers
			if (sum(fixed_TF_to_use, na.rm=TRUE) > 0)
				{
				groupnames_that_have_fixed_ages = groupnames_to_use[fixed_TF_to_use]
				priornames_for_groupnames_to_use2[fixed_TF_to_use] = fixed_ages_to_use[fixed_TF_to_use]
				}
			
			priornames_for_groupnames_to_use3 = paste(priornames_for_groupnames_to_use2, sep="", collapse="+")
			mean_equation_txt = paste0("((", priornames_for_groupnames_to_use3, ") / ", num_groups, ")")
			mean_equation_txt
			
			# Make equations for difference from the mean
			equations_txt = paste(priornames_for_groupnames_to_use2, mean_equation_txt, sep=" - ")
			equations_txt
			
			# Names of the expression results
			tmpname = "diff_from_mean"
			equation_names = paste(groupnames_to_use, tmpname, "FOR", rel_prior_Name, sep="_")
			equation_names
			prior_names = paste0("", equation_names, sep="")
			
			# Write x inputs XML
			x_XMLs = list(bl(), xmlCommentNode(" x inputs to equation "))
			for (g in 1:length(groupnames_to_use))
				{
				# Skip putting the x in, if it's actually a constant...
				if (fixed_TF_to_use[g] == TRUE)
					{
					next()
					} # END if (fixed_TF_to_use[g] == TRUE)
				x_XML = xmlNode(name="x", attrs=list(idref=priornames_for_groupnames_to_use[g]))
				x_XMLs = c(x_XMLs, list(x_XML))
				} # END for (g in 1:length(groupnames_to_use))
			
			# Write and store XML for the equation, priors, etc.
			expression_XML = list(bl(), xmlCommentNode(paste0(" Expressions calculating the difference from the mean for each groupname listed in the row named '", rel_prior_Name, "' ")))
			prior_XMLs = list(bl(), xmlCommentNode(paste0(" Expressions calculating the difference from the mean for each groupname listed in the row named '", rel_prior_Name, "' ")))
			for (g in 1:length(groupnames_to_use))
				{
				# Write the equation XML
				expression_XML = xmlNode(name="parameter", attrs=list(id=equation_names[g], spec="beast.util.Script", expression=equations_txt[g]), .children=x_XMLs)
				expressions_XMLs = c(expressions_XMLs, list(expression_XML))
				
				# Write the prior XML
				# Aaaand, write the priors
				dfline = specs
				tmpXML = make_generic_XML_prior(dfline, colname_prefix="rel", param_name=prior_names[g], header_scheme=1, distrib=NULL, param1=NULL, param2=NULL, tmp_offset=NULL, meanInRealSpace=NULL)
				relpriors_XMLs = c(relpriors_XMLs, tmpXML$priors)
				logs_of_expressions_XMLs = c(logs_of_expressions_XMLs, tmpXML$tracelog)
		
				# Log the expressions
				log_of_expression_XML = xmlNode(name="log", attrs=list(idref=equation_names[g]))
				logs_of_expressions_XMLs = c(logs_of_expressions_XMLs, list(log_of_expression_XML))
				screenlogs_of_expressions_XMLs = c(screenlogs_of_expressions_XMLs, list(log_of_expression_XML))
		
				# Log the priors (redundant)
				#log_of_relprior_XML = xmlNode(name="log", attrs=list(idref=prior_names[g]))
				#logs_of_relpriors_XMLs = c(logs_of_relpriors_XMLs, list(log_of_relprior_XML))
				} # END for (g in 1:length(groupnames_to_use))
			
			next()	# skip to next iteration of the main for-loop
			} # END if (equation_to_modify == "diff_from_mean")
		
		
		
		# Find the groupnames in the equation
		TF = rep(FALSE, length(groupnames))
		for (i in 1:length(groupnames))
			{
			tmp = gregexpr(pattern=groupnames[i], text=equation_to_modify)
			search_result = tmp[[1]]
			if (length(search_result) > 1)
				{
				errortxt = paste0("\n\nSTOP ERROR in make_generic_difference_statistic(): While parsing the 'Equation' columns in the OTUs/nodes worksheets, an error occurred because a groupname appeared twice in an equation. Change groupnames to avoid this. Printing equation and groupname searched on:\n\n") 
				cat(errortxt)
				
				cat("Equation:\n")
				cat(equation_to_modify)
				cat("\n\ngroupname:\n")
				cat(groupnames[i])
				cat("\n\n")
				
				stop(errortxt)
				} # END if (length(search_result) > 1)
			if (search_result != -1)
				{
				TF[i] = TRUE
				}
			} # END for (i in groupnames)
		
		# If no groupnames were found, you have a problem!
		if (sum(TF) == 0)
			{
			errortxt = paste0("\n\nSTOP ERROR in make_generic_difference_statistic(): While parsing the 'Equation' columns in the OTUs/nodes worksheets, no groupnames were found. Printing equation and groupnames searched on:\n\n") 
			cat(errortxt)
			cat("Equation:\n")
			cat(equation_to_modify)
			cat("\n\ngroupnames:\n")
			cat(groupnames)
			cat("\n\n")
			stop(errortxt) 
			} # END if (sum(TF) == 0)
		
		# Now, replace the groupnames (since we checked that they only appear once)
		groupnames_to_use = groupnames[TF]
		priornames_for_groupnames_to_use = priornames_for_groupnames[TF]
		fixed_TF_to_use = fixed_TF[TF]
		fixed_ages_to_use = fixed_ages[TF]
		x_XML_list = list(xmlCommentNode(" x's give the inputs to the equation "))
		for (g in 1:length(groupnames_to_use))
			{
			old_groupname = groupnames_to_use[g]
			xmlScript_txt = priornames_for_groupnames_to_use[g]
			xmlScript_txt2 = paste0(xmlScript_txt, "[1]")
			
			# Check for fixed, replace with number if so
			if (sum(fixed_TF_to_use[g], na.rm=TRUE) > 0)
				{
				xmlScript_txt2 = fixed_ages_to_use[g]
				}
			
			equation_to_modify = gsub(pattern=old_groupname, replacement=xmlScript_txt2, x=equation_to_modify)
			
			
			# Only add to x_XMLs if NOT fixed
			if (sum(fixed_TF_to_use[g], na.rm=TRUE) == 0)
				{
				# Make x IDrefs for inside of beast.util.Script expression XML
				x_XML_tmp = xmlNode(name="x", attrs=list(idref=xmlScript_txt))
				x_XML_list = c(x_XML_list, list(x_XML_tmp))
				} # END if (sum(fixed_TF_to_use[g], na.rm=TRUE) == 0)
			
			} # END for (g in 1:length(groupnames_to_use))
		
		# Now, write the XML inputs for doing the equation
		equation_XML = xmlNode(name="parameter", attrs=list(id=rel_prior_Name, expression=equation_to_modify, spec="beast.util.Script"), .children=x_XML_list)
		expressions_XMLs = c(expressions_XMLs, list(equation_XML))
		
		# Aaaand, write the priors
		dfline = specs
		tmpXML = make_generic_XML_prior(dfline, colname_prefix="rel", param_name=rel_prior_Name, header_scheme=1, distrib=NULL, param1=NULL, param2=NULL, tmp_offset=NULL, meanInRealSpace=NULL)
		relpriors_XMLs = c(relpriors_XMLs, tmpXML$priors)
		logs_of_expressions_XMLs = c(logs_of_expressions_XMLs, tmpXML$tracelog)
		
		# Log the expressions
		log_of_expression_XML = xmlNode(name="log", attrs=list(idref=rel_prior_Name))
		logs_of_expressions_XMLs = c(logs_of_expressions_XMLs, list(log_of_expression_XML))
		screenlogs_of_expressions_XMLs = c(screenlogs_of_expressions_XMLs, list(log_of_expression_XML))
		
		# Log the priors
		log_of_relprior_XML = xmlNode(name="log", attrs=list(idref=rel_prior_Name))
		logs_of_relpriors_XMLs = c(logs_of_relpriors_XMLs, list(log_of_relprior_XML))
		} # END for (rownum in 1:nrow(rel_prior_specs_for_groupnames))
	
	# Output
	if (is.null(xml))
		{
		res = NULL
		res$misc = expressions_XMLs
		res$priors = relpriors_XMLs
		res$tracelog = c(logs_of_expressions_XMLs, logs_of_relpriors_XMLs)
		res$screenlog = screenlogs_of_expressions_XMLs
		extract='
		expressions_XMLs = res$misc
		relpriors_XMLs = res$priors
		tracelog_XMLs = res$tracelog
		screenlogs_of_expressions_XMLs = res$screenlog
		'
		return(res)
		} else {
		xml$misc = c(xml$misc, expressions_XMLs)
		xml$priors = c(xml$priors, relpriors_XMLs)
		xml$tracelog = c(xml$tracelog, logs_of_expressions_XMLs, logs_of_relpriors_XMLs)
		xml$screenlog = c(xml$screenlog, screenlogs_of_expressions_XMLs)
		return(xml)
		} # END if (is.null(xml))
	stop("ERROR in make_generic_difference_statistic(): shouldn't get here.")
	} # END make_generic_difference_statistic <- function(nodes_df=NULL, OTUs_df=NULL)


