library(XML)


#######################################################
# Make OTUs list
#######################################################

# For use in sapply by make_XMLs_for_OTUs
make_XML_for_OTU_id <- function(OTU)
	{
	xmlNode(name="taxon", attrs=list(id=OTU, spec="Taxon"))
	}

make_XML_for_OTU_idref <- function(OTU)
	{
	xmlNode(name="taxon", attrs=list(idref=OTU))
	}

# If OTU_idref=TRUE, the taxa will be referred to with "idref" tags.
# If OTU_idref=FALSE, the taxa will be referred to with "id" tags 
#    and specified to be spec="Taxon"
make_XMLs_for_OTUs <- function(OTUs, OTU_idref=TRUE)
	{
	defaults='
	OTU_idref=TRUE
	'
	
	if (OTU_idref == TRUE)
		{
		list_of_OTUnames_XML = sapply(X=OTUs, FUN=make_XML_for_OTU_idref)
		} else {
		list_of_OTUnames_XML = sapply(X=OTUs, FUN=make_XML_for_OTU_id)
		}
	XML_list_of_OTUs = unname(list_of_OTUnames_XML)
	return(XML_list_of_OTUs)
	}

make_XML_taxon_block <- function(taxon_name, XML_list_of_OTUs, xml=NULL)
	{
	txt = paste0(" A taxonset named ", taxon_name, ". Can be used for clade constraints and/or logging. ")
	XML_comment = xmlCommentNode(txt)
	XML_comment
	
	XML_taxonset = xmlNode(name="taxonset", attrs=list(id=taxon_name, spec="TaxonSet"), .children=XML_list_of_OTUs)
	
	XML_taxonset = c(list(XML_comment), list(XML_taxonset))
	XML_taxonset
	
	if (is.null(xml))
		{
		return(XML_taxonset)
		} else {
		# Add to the list of taxa/clades
		xml$taxa = c(xml$taxa, XML_taxonset)
		return(xml)
		} # END if (is.null(xml))
	
	cat("\n\n")
	stop("ERROR in make_XML_taxon_block(): Shouldn't get here.")
	}







#######################################################
# Tip dates
#######################################################


make_txt_tipdate <- function(OTU, tipdate, leading_space="\t\t\t", trailing=",\r")
	{
	txt = paste(leading_space, OTU, "=", tipdate, sep="")
	return(txt)
	}

make_txt_tipdates_list <- function(OTUs, tipdates, leading_space="\t\t\t", trailing=",\n")
	{
	defaults='
	leading_space="\t\t\t"
	trailing=",\n"
	'
	
	list_of_tipdates_txt = unname(mapply(FUN=make_txt_tipdate, OTU=OTUs, tipdate=tipdates, MoreArgs=list(leading_space=leading_space, trailing="")))
	list_of_tipdates_txt
	
	txt = paste0(list_of_tipdates_txt, collapse=trailing)
	cat(txt)

	return(txt)
	}

# Make priors for individual tipdates
make_XML_tipdate_priors <- function(OTUs_df, min_precision=0.0001, tree_name="shared_tree", monophyletic="false", xml=NULL)
	{
	defaults = '
	min_precision=0.01
	tree_name="shared_tree"
	xml=NULL
	'
	
	# Remove OTUs marked with use="no"
	OTUs_df$use[isblank_TF(OTUs_df$use)] = "yes"
	OTUs_df = OTUs_df[OTUs_df$use != "no", ]
	
	# Find the tipdates that are not fixed, that have non-uniform distributions
	nonfixed_TF = OTUs_df$distribution != "fixed"
	# Exit if no priors on tip dates
	if (sum(nonfixed_TF) < 1)
		{
		if (is.null(xml))
			{
			tmpXML = NULL
			return(tmpXML)
			} else {
			return(xml)
			} # END if (is.null(xml))
		} # END if (sum(num_priors) < 1)

	# Error checks
	check_tipdates(OTUs_df)


	fossil_TF = OTUs_df$tipdate > min_precision
	num_fossils = sum(fossil_TF)
	prior_TF = (nonfixed_TF + fossil_TF) == 2
	num_priors = sum(prior_TF)
	
	
	# Write the XML for the priors and logs
	OTU_priors_df = OTUs_df[prior_TF,]

	if (is.null(xml))
		{
		for (i in 1:num_priors)
			{
			rowdf = OTU_priors_df[i,]
			if (i == 1)
				{
				tmpXML = rowdf_to_XML_distribution(rowdf, tipsonly="true", monophyletic=monophyletic, tree_name=tree_name, xml=NULL)
				} else {
				tmpXML = rowdf_to_XML_distribution(rowdf, tipsonly="true", monophyletic=monophyletic, tree_name=tree_name, xml=tmpXML)
				} # END if (i == 1)
			} # END for (i in 1:num_priors)
		

		# Add some handy notation -- log prP
		logtxt = " Log the sampled (individual) tipdates and the prior probability densities (prP) of them "
		logtxt_XML = xmlCommentNode(logtxt)
		tmpXML$tracelog = c(list(bl()), list(logtxt_XML), tmpXML$tracelog)

		# Add some handy notation -- actual tipdates
		#logtxt = " Log the sampled (individual) tipdates "
		#logtxt_XML = xmlCommentNode(logtxt)
		#tmpXML$tracelog = c(tmpXML$tracelog, list(bl()), list(logtxt_XML), tmpXML$tipdatelog)

		operators_txt = " Operators to change the dates of the (individual) fossil tipdates "
		operators_txt_XML = xmlCommentNode(operators_txt)
		tmpXML$operators = c(list(bl()), list(operators_txt_XML), tmpXML$operators)

		return(tmpXML)
		} else {
		for (i in 1:num_priors)
			{
			rowdf = OTU_priors_df[i,]
			xml = rowdf_to_XML_distribution(rowdf, tipsonly="true", monophyletic=monophyletic, tree_name=tree_name, xml=xml)
			} # END for (i in 1:num_priors)


		# Add some handy notation -- log prP
		logtxt = " Log the sampled (individual) tipdates and the prior probability densities (prP) of them "
		logtxt_XML = xmlCommentNode(logtxt)
		xml$tracelog = c(list(bl()), list(logtxt_XML), xml$tracelog)
		
		# Add some handy notation -- actual tipdates
		#logtxt = " Log the sampled (individual) tipdates "
		#logtxt_XML = xmlCommentNode(logtxt)
		#xml$tracelog = c(xml$tracelog, list(bl()), list(logtxt_XML), xml$tipdatelog)
		
		operators_txt = " Operators to change the dates of the (individual) fossil tipdates "
		operators_txt_XML = xmlCommentNode(operators_txt)
		xml$operators = c(list(bl()), list(operators_txt_XML), xml$operators)

		return(xml)
		} # END if (is.null(xml))
	} # END make_XML_tipdate_priors


rowdf_to_XML_distribution <- function(rowdf, tipsonly="true", monophyletic="false", tree_name="shared_tree", xml=NULL)
	{
	defaults='
	tipsonly="true"
	monophyletic="false"
	tree_name="shared_tree"
	xml=NULL
	'
	
	# Get the name
	OTUname = rowdf$OTUs
	OTUtaxonset = paste0("age_of_", OTUname)
	tree_idref = paste0("@", tree_name)

	# Convert to normal: default "no"
	if (isblank_TF(rowdf$convert_to_normal) == TRUE)
		{
		rowdf$convert_to_normal = "no"
		}
	
	# Get the offset, if it exists
	if (isblank_TF(rowdf$offset) == FALSE)
		{
		offset = rowdf$offset
		} else {
		offset = 0
		} # END if (isblank_TF(rowdf$offset) == FALSE)

	# Get the meanInRealSpace, if it exists
	if (( isblank_TF(rowdf$meanInRealSpace) == FALSE) && (rowdf$meanInRealSpace == "no") )
		{
		meanInRealSpace = FALSE
		meanInRealSpace_txt = "false"
		} else {
		meanInRealSpace = TRUE
		meanInRealSpace_txt = "true"
		}
	
	if (rowdf$distribution == "normal")
		{
		OTU_priorname = paste0("prP_NormalDistrib_", OTUname)
		
		if (offset > 0)
			{
			errortxt = paste0("\n\nSTOP ERROR in rowdf_to_XML_distribution(): ", OTUname, " has a Normal distribution but the offset is >0. offset=", offset, "\n\n")
			cat(errortxt)
			stop(errortxt)
			} # END if (offset > 0)
		
		meanval = rowdf$param1 + offset
		sdval = rowdf$param2
		
		distribution_name = paste0("NormalDistrib_", OTUname)
		param1_name = paste0("mean_of_NormalDistrib_", OTUname)
		param2_name = paste0("sd_of_NormalDistrib_", OTUname)
		param1_XML = xmlNode(name="parameter", rowdf$param1, attrs=list(id=param1_name, name="mean", estimate="false"))
		param2_XML = xmlNode(name="parameter", rowdf$param2, attrs=list(id=param2_name, name="sigma", estimate="false"))
		distrib_XML = xmlNode(name="Normal", attrs=list(id=distribution_name, name="distr", offset="0"), .children=list(param1_XML, param2_XML))
		
		txt = paste0(" Prior probability density of the date of '", OTUname, "', according to a ", rowdf$distribution, " distribution with mean=", meanval, ", sd=", sdval, ". ")
		} # END if (rowdf$distribution == "normal")

	if (rowdf$distribution == "lognormal")
		{
		OTU_priorname = paste0("prP_LogNormalDistrib_", OTUname)

		meanval = rowdf$param1
		sdval = rowdf$param2

		# Convert to real space, if needed
		if (meanInRealSpace == FALSE)
			{
			meanval = exp(meanval)
			sdval = exp(sdval)
			}
		
		# Convert to normal, if desired
		if (rowdf$convert_to_normal == "yes")
			{
			normMeanval = meanval + offset
			#normSDval = sdval
			normSDval = 0.01 * normMeanval
			
			distribution_name = paste0("convertedForStartingTree_NormalDistrib_", OTUname)
			param1_name = paste0("convertedForStartingTree_mean_of_NormalDistrib_", OTUname)
			param2_name = paste0("convertedForStartingTree_sd_of_NormalDistrib_", OTUname)
			param1_XML = xmlNode(name="parameter", normMeanval, attrs=list(id=param1_name, name="mean", estimate="false"))
			param2_XML = xmlNode(name="parameter", normSDval, attrs=list(id=param2_name, name="sigma", estimate="false"))
			distrib_XML = xmlNode(name="Normal", attrs=list(id=distribution_name, name="distr", offset="0"), .children=list(param1_XML, param2_XML))

			txt = paste0(" As convert_to_normal==yes, for getting a starting tree, the '", OTUname, "' prior is a normal with mean=", normMeanval, ", sd=", normSDval, ". ")		
			} else {
			# Use the specified distribution, DON'T convert to normal
			distribution_name = paste0("LogNormalDistrib_", OTUname)
			param1_name = paste0("mean_of_LogNormalDistrib_", OTUname)
			param2_name = paste0("sd_of_LogNormalDistrib_", OTUname)
			param1_XML = xmlNode(name="parameter", rowdf$param1, attrs=list(id=param1_name, name="M", estimate="false"))
			param2_XML = xmlNode(name="parameter", rowdf$param2, attrs=list(id=param2_name, name="S", estimate="false"))
			distrib_XML = xmlNode(name="LogNormal", attrs=list(id=distribution_name, name="distr", offset=offset, meanInRealSpace="true"), .children=list(param1_XML, param2_XML))

			txt = paste0(" Prior probability density of the date of '", OTUname, "', according to a ", rowdf$distribution, " distribution with mean=", meanval, ", sd=", sdval, ". ")
			} # END if (rowdf$convert_to_normal == "yes")
		} # END if (rowdf$distribution == "lognormal")

	if (rowdf$distribution == "exponential")
		{
		OTU_priorname = paste0("prP_ExponentialDistrib_", OTUname)

		meanval = rowdf$param1

		# Convert to real space, if needed
		if (meanInRealSpace == FALSE)
			{
			meanval = 1/(meanval)
			}

		# Convert to normal, if desired
		if (rowdf$convert_to_normal == "yes")
			{
			normMeanval = meanval + offset
			# For conversion to normal
			normSDval = 0.01 * meanval
			
			distribution_name = paste0("convertedForStartingTree_NormalDistrib_", OTUname)
			param1_name = paste0("convertedForStartingTree_mean_of_NormalDistrib_", OTUname)
			param2_name = paste0("convertedForStartingTree_sd_of_NormalDistrib_", OTUname)
			param1_XML = xmlNode(name="parameter", normMeanval, attrs=list(id=param1_name, name="mean", estimate="false"))
			param2_XML = xmlNode(name="parameter", normSDval, attrs=list(id=param2_name, name="sigma", estimate="false"))
			distrib_XML = xmlNode(name="Normal", attrs=list(id=distribution_name, name="distr", offset="0"), .children=list(param1_XML, param2_XML))

			txt = paste0(" As convert_to_normal==yes, for getting a starting tree, the '", OTUname, "' prior is a normal with mean=", normMeanval, ", sd=", normSDval, ". ")		
			} else {
			# Use the EXPONENTIAL distribution, DON'T convert to normal
			distribution_name = paste0("ExponentialDistrib_", OTUname)
			param1_name = paste0("mean_of_ExponentialDistrib_", OTUname)
			param1_XML = xmlNode(name="parameter", rowdf$param1, attrs=list(id=param1_name, name="mean", estimate="false"))
			distrib_XML = xmlNode(name="Exponential", attrs=list(id=distribution_name, name="distr", offset=offset), .children=list(param1_XML))

			txt = paste0(" Prior probability density of the date of '", OTUname, "', according to a ", rowdf$distribution, " distribution with mean=", meanval, ". ")
			} # END if (rowdf$convert_to_normal == "yes")
		} # END if (rowdf$distribution == "exponential")

	if (rowdf$distribution == "uniform")
		{
		OTU_priorname = paste0("prP_UniformDistrib_", OTUname)

		if (offset > 0)
			{
			errortxt = paste0("\n\nSTOP ERROR in rowdf_to_XML_distribution(): ", OTUname, " has a Uniform distribution but the offset is >0. offset=", offset, "\n\n")
			cat(errortxt)
			stop(errortxt)
			} # END if (offset > 0)

		lower = rowdf$param1 + offset
		upper = rowdf$param2 + offset

		# Convert to normal, if desired
		if (rowdf$convert_to_normal == "yes")
			{
			normMeanval = (lower + upper) / 2
			# For conversion to normal
			normSDval = 0.01 * normMeanval
			
			distribution_name = paste0("convertedForStartingTree_NormalDistrib_", OTUname)
			param1_name = paste0("convertedForStartingTree_mean_of_NormalDistrib_", OTUname)
			param2_name = paste0("convertedForStartingTree_sd_of_NormalDistrib_", OTUname)
			param1_XML = xmlNode(name="parameter", normMeanval, attrs=list(id=param1_name, name="mean", estimate="false"))
			param2_XML = xmlNode(name="parameter", normSDval, attrs=list(id=param2_name, name="sigma", estimate="false"))
			distrib_XML = xmlNode(name="Normal", attrs=list(id=distribution_name, name="distr", offset="0"), .children=list(param1_XML, param2_XML))

			txt = paste0(" As convert_to_normal==yes, for getting a starting tree, the '", OTUname, "' prior is a normal with mean=", normMeanval, ", sd=", normSDval, ". ")		
			} else {
			# Use the UNIFORM distribution, DON'T convert to normal
			distribution_name = paste0("UniformDistrib_", OTUname)
			param1_name = paste0("lower_of_UniformDistrib_", OTUname)
			param2_name = paste0("upper_of_UniformDistrib_", OTUname)
			distrib_XML = xmlNode(name="Uniform", attrs=list(id=distribution_name, name="distr", offset="0", lower=lower, upper=upper))

			txt = paste0(" Prior probability density of the date of '", OTUname, "', according to a ", rowdf$distribution, " distribution with lower=", lower, ", upper=", upper, ". ")
			} # END if (rowdf$convert_to_normal == "yes")

		} # END if (rowdf$distribution == "uniform")
	
	
	########################################################################
	# Now, make the prior probability distribution for this taxonset
	########################################################################	
	# TaxonSet
	taxon_XML = xmlNode(name="taxon", attrs=list(idref=OTUname, spec="Taxon"))
	taxonset_XML = xmlNode(name="taxonset", attrs=list(id=OTUtaxonset, spec="TaxonSet"), .children=list(taxon_XML))
	
	# Overall Prior Distribution
	prior_XML = xmlNode(name="distribution", attrs=list(id=OTU_priorname, tipsonly=tipsonly, monophyletic=monophyletic, tree=tree_idref, spec="beast.math.distributions.MRCAPrior"), .children=list(taxonset_XML, distrib_XML))
	
	comment_XML = xmlCommentNode(txt)
	
	prior_XML_plus_comment = list(bl(), comment_XML, prior_XML)


	########################################################################
	# Log the prior density of each tipdate
	########################################################################
	log_XML = xmlNode(name="log", attrs=list(idref=OTU_priorname))
	
	# Log the actual age of each tipdate
	# Doesn't work, tip mismatch
	# TIP HEIGHTS ARE LOGGED AS E.G. "height(H_erectus)"
	#agelog_XML = xmlNode(name="log", attrs=list(idref=taxonset_XML))

	########################################################################
	# Log the actual sampled date of each tipdate
	########################################################################
	# This doesn't work, produces:
	# 
	# type mismatch for input 
	# log. beast.core.Loggable.isAssignableFrom(class
	# beast.evolution.alignment.TaxonSet)=false 
	# expected 'Loggable' but got 'TaxonSet'
	# 
	#date_of_tip_log_XML = xmlNode(name="log", attrs=list(idref=OTUtaxonset))


	########################################################################
	# Operators on each individual tipdates
	########################################################################
	taxonset_idref = paste0("@", OTUtaxonset)
	operator_id = paste0("TipDatesRandomWalker_", OTUname)
	
	# Make the window width 1/2 of the 99.9 CI
	rowdf_as_uniform = convert_nonUniform_dates_to_uniform(rowdf, CI=0.999)
	window_width = (rowdf_as_uniform$param2 - rowdf_as_uniform$param1) / 2
	operator_XML = xmlNode(name="operator", attrs=list(id=operator_id, windowSize=window_width, taxonset=taxonset_idref, tree=tree_idref, weight="1.0", spec="beast.evolution.operators.TipDatesRandomWalker"))

	
	
	if (is.null(xml))
		{
		tmpXML = NULL
		tmpXML$priors = prior_XML_plus_comment
		tmpXML$operators = list(operator_XML)
		#tmpXML$tipdatelog = list(date_of_tip_log_XML)
		tmpXML$tracelog = list(log_XML)
		return(tmpXML)
		} else {
		xml$priors = c(xml$priors, prior_XML_plus_comment)
		xml$operators = c(xml$operators, list(operator_XML))
		#xml$tipdatelog = c(xml$tipdatelog, list(date_of_tip_log_XML))
		xml$tracelog = c(xml$tracelog, list(log_XML))
		return(xml)
		} # END if (is.null(xml))
	} # END rowdf_to_XML_distribution <- function(rowdf, tipsonly="true", monophyletic="false", tree_name="shared_tree", xml=NULL)


make_XML_tipdates <- function(name="tipDates", OTUs_df, alignment_name_w_taxa, backward=TRUE, xml=NULL, min_precision=0.01)
	{
	defaults='
	name="tipDates"
	taxon_name = "all_taxa"
	backward=TRUE
	traitname = "date-backward"
	'
	
	# Subset to just OTUs to use != "no"
	OTUs_df = OTUs_df[OTUs_df$use != "no",]
	
	# Get the OTUs
	OTUs = OTUs_df$OTUs
	
	
	tipdates = OTUs_df$tipdate
	tipdates = convert_blanks_NAs_to_0(tipdates)
	
	num_fossils = sum(tipdates > min_precision)
	
	if ( length(OTUs) != length(tipdates) )
		{
		txt = paste0("\n\nERROR in make_XML_tipdates(): length(OTUs) != length(tipdates)\n (length(OTUs) = ", length(OTUs), ", length(tipdates) = ", length(tipdates), "\n\n")
		cat(txt)
		stop(txt)
		}
	
	
	units = "year"
	if (backward == TRUE)
		{
		traitname = "date-backward"
		} else {
		traitname = "date-forward"
		}
	
	
	tipdates_txt = make_txt_tipdates_list(OTUs, tipdates, leading_space="", trailing=",\n")
	tipdates_txt = paste0(tipdates_txt, sep="")
	
	txt1 = " Tipdates for all taxa. Tipdates that will vary need to be   "
	txt2 = " specified in operators, and also have a prior distribution. "
	XML_comment1 = xmlCommentNode(txt1)
	XML_comment2 = xmlCommentNode(txt2)
	
	# Specify the taxonSet which the tipdate names refer to
	txt3 = " 'taxa=' tag specifies the beast.evolution.alignment.TaxonSet which the tipdate names refer to. "
	XML_comment3 = xmlCommentNode(txt3)

	txt4 = paste0(" Number of fossils: you have ", num_fossils, " tips with ages > 0 Ma ")
	XML_comment4 = xmlCommentNode(txt4)

	alignment_name_w_taxa_idref = paste0("@", alignment_name_w_taxa)
	taxa_names_source_XML = xmlNode(name="taxa", attrs=list(alignment=alignment_name_w_taxa_idref, spec="beast.evolution.alignment.TaxonSet"))
	
	XML_tipdates = xmlNode(name="trait", attrs=list(id=name, spec="beast.evolution.tree.TraitSet", traitname=traitname, units=units, value=tipdates_txt))
	XML_tipdates
	
	# Don't add children...
	XML_tipdates = addChildren(XML_tipdates, kids=list(XML_comment3, taxa_names_source_XML))
	XML_tipdates
	
	# Add the comments
	XML_tipdates = c(list(bl()), list(bl()), list(XML_comment1), list(XML_comment2), list(XML_comment4), list(XML_tipdates))
	
	XML_tipdates
	
	if (is.null(xml))
		{
		return(XML_tipdates)
		} else {
		xml$taxa = c(xml$taxa, XML_tipdates)
		return(xml)
		} # END if (is.null(xml))
	
	cat("\n\n")
	stop("ERROR in make_XML_tipdates(): shouldn't get here.")
	}


add_fossils_taxon_to_xml <- function(OTUs, tipdates, xml)
	{
	TF = tipdates > 0
	
	# Don't modify the XML, if no fossils
	if (sum(TF) == 0)
		{
		return(xml)
		}
	
	# Otherwise, process the XML accordingly...
	fossil_OTUs = OTUs[TF]
	XML_list_of_OTUs = make_XMLs_for_OTUs(fossil_OTUs, OTU_idref=TRUE)
	XML_fossil_OTUs = make_XML_taxon_block(taxon_name="fossil_taxa", XML_list_of_OTUs=XML_list_of_OTUs)
	XML_fossil_OTUs = c(list(bl()), list(bl()), XML_fossil_OTUs)
	XML_fossil_OTUs

	xml$taxa = c(xml$taxa, XML_fossil_OTUs)
	xml$taxa
	return(xml)
	} # END add_fossils_taxon_to_xml



# Convert blanks etc. to 0
convert_blanks_NAs_to_0 <- function(tipdates, newval=0)
	{
	defaults='
	newval=0
	'
	
	TF = tipdates == ""
	tipdates[TF] = newval

	TF = tipdates == " "
	tipdates[TF] = newval

	TF = tipdates == "\t"
	tipdates[TF] = newval

	TF = is.na(tipdates)
	tipdates[TF] = newval

	TF = is.nan(tipdates)
	tipdates[TF] = newval
	
	return(tipdates)
	}


isblank_TF <- function(items)
	{
	TF0 = is.null(items)
	if (TF0 == TRUE)
		{
		blank_TF = TRUE
		return(blank_TF)
		} # END if (TF0 == TRUE)
	TF1 = items == ""
	TF2 = items == " "
	TF3 = items == "\t"
	TF4 = is.na(items)
	TF5 = is.nan(items)
	
	# Keep only the items where none of the above occur
	TFall = TF1 + TF2 + TF3 + TF4 + TF5
	
	# Correct for NA, NaNs
	TFall[is.na(TFall)] = 1
	TFall[is.nan(TFall)] = 1

	blank_TF = TFall > 0
	return(blank_TF)
	}


# Remove blanks etc. from a list
remove_blanks_NAs_etc <- function(items)
	{
	# Gotta remove NAs and NANs first
	NAs_TF = is.na(items)
	items[NAs_TF] = ""
	NANs_TF = is.nan(items)
	items[NANs_TF] = ""
	
	TF1 = items == ""
	TF2 = items == " "
	TF3 = items == "\t"
	TF4 = is.na(items)
	TF5 = is.nan(items)
	
	# Keep only the items where none of the above occur
	TF = TF1 + TF2 + TF3 + TF4 + TF5
	items = items[TF == 0]
	
	return(items)
	}

get_firstword <- function(OTU, split="_")
	{
	firstword = strsplit(x=OTU, split=split)[[1]][1]
	}


get_firstword_from_OTUs <- function(OTUs, split="_")
	{
	unname(sapply(X=OTUs, FUN=get_firstword, split=split))
	}


get_genera_from_OTUs <- function(OTUs, mintaxa=2, split="_")
	{
	defaults='
	mintaxa=2
	split="_"
	'
	
	# Assuming that species in the same genus will have the 
	# same first word, get a list of the genera
	firstwords = sort(get_firstword_from_OTUs(OTUs, split=split))
	firstwords_counts = table(firstwords)
	firstwords_keepTF = firstwords_counts >= mintaxa
	
	# Remove genera which have fewer than mintaxa representatives
	genera = names(firstwords_counts)[firstwords_keepTF]

	return(genera)
	}


make_XML_for_genera_from_OTUs <- function(xml, OTUs, mintaxa=2, split="_")
	{
	defaults='
	mintaxa=2
	split="_"
	'
	
	genera = get_genera_from_OTUs(OTUs, mintaxa=mintaxa, split=split)
	
	if (length(genera) < 1)
		{
		txt = paste("\n\nmake_XML_for_genera_from_OTUs: No genera found with\nmintaxa >= ", mintaxa, ", using split=", split, ". Returning input 'xml'.\n")
		cat(txt)
		return(xml)
		}
	
	for (i in 1:length(genera))
		{
		genus_name = genera[i]
		pattern = paste0(genus_name, split)
		TF = grepl(pattern=pattern, x=OTUs)
		tmp_OTUs = OTUs[TF]
		
		XML_list_of_OTUs = make_XMLs_for_OTUs(OTUs=tmp_OTUs, OTU_idref=TRUE)
		taxon_XML = make_XML_taxon_block(taxon_name=genus_name, XML_list_of_OTUs=XML_list_of_OTUs)
		
		xml$taxa = c(xml$taxa, list(bl()), taxon_XML)
		}
	
	return(xml)
	} # END make_XML_for_genera_from_OTUs

