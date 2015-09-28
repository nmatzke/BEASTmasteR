#######################################################
# Add priors for the node constraints
#######################################################

# If OTUs is not null, then groups will be pruned of taxa missing from OTUs
make_taxa_groups <- function(taxa_df, OTUs=NULL, xml=NULL)
	{
	defaults='
	OTUs=NULL
	'
	
	names_of_groups = names(taxa_df)
	# User should manually create an "all_OTUs" group first
	#names_of_groups = names_of_groups[names_of_groups != "all_taxa"]

	if (length(names_of_groups) < 1)
		{
		if (is.null(xml))
			{
			return(NULL)
			} else {
			return(xml)
			} # END if (is.null(xml))
		} # END if (length(names_of_groups) < 1)

	# Check names
	keepTF = rep(TRUE, times=length(names_of_groups))
	for (i in 1:length(names_of_groups))
		{
		items = taxa_df[,names_of_groups[i]]
		items2 = remove_blanks_NAs_etc(items)
	
		# Keep if the column has one or more non-blank taxa
		if (length(items2) < 1)
			{
			keepTF[i] = FALSE
			} # END if (length(items2) < 1)
		} # END for (i in 1:length(names_of_groups))

	names_of_groups = names_of_groups[keepTF]

	if (length(names_of_groups) < 1)
		{
		if (is.null(xml))
			{
			return(NULL)
			} else {
			return(xml)
			} # END if (is.null(xml))
		} # END if (length(names_of_groups) < 1)
	
	# Keep a list of the taxa which end up EMPTY -- we will CUT 
	# these from the priors list (or else it causes an error)
	list_of_empty_taxa = NULL
	
	# Otherwise, go through and make these taxa
	taxon_XMLs = list(bl(), xmlCommentNode(" Taxon groupings (clades or other groups of interest) "))
	for (i in 1:length(names_of_groups))
		{
		groupname = names_of_groups[i]
		items = taxa_df[,groupname]
		items2 = remove_blanks_NAs_etc(items)
		
		# Filter out names not found in OTUs; throw warning for each
		if (is.null(OTUs) == FALSE)
			{
			TF = items2 %in% OTUs
			
			if (sum(TF) == 0)
				{
				txt = paste0("WARNING in make_taxa_groups(): None of the taxa names in '", groupname, "' are found in overall list in 'OTUs'. The group '", groupname, "' is being left out of the XML.")
				cat("\n")
				cat(txt)
				cat("\n")
				warning(txt)
				cat("\n")
				
				list_of_empty_taxa = c(list_of_empty_taxa, groupname)
				
				# Skip to next
				next()
				} # END if (sum(TF) == 0)
			
			items_being_removed = items2[TF == FALSE]
			
			if (length(items_being_removed) > 0)
				{
				for (j in 1:length(items_being_removed))
					{
					txt = paste0("'", items_being_removed[j], "' is missing from overall 'OTUs' so is being cut from '", groupname, "'.")
					cat("\n")
					cat(txt)
					cat("\n")
					warning(txt)
					cat("\n")
					} # END for (j in 1:length(items_being_removed))
				} # END if (length(items_being_removed) > 0)
			# Subset the OTUs
			items2 = items2[TF]
			} # END if (is.null(OTUs) == FALSE)
		
		
		XML_list_of_OTUs = make_XMLs_for_OTUs(OTUs=items2, OTU_idref=TRUE)
		taxon_XML = make_XML_taxon_block(taxon_name=groupname, XML_list_of_OTUs=XML_list_of_OTUs)
		
		taxon_XMLs = c(taxon_XMLs, taxon_XML)
		} # END for (i in 1:length(names_of_groups))
	
	# Return results
	if (is.null(xml))
		{
		res = NULL
		res$taxon_XMLs = taxon_XMLs
		res$list_of_empty_taxa = list_of_empty_taxa
		
		extract='
		taxon_XMLs = res$taxon_XMLs
		list_of_empty_taxa = res$list_of_empty_taxa
		'
		
		return(res)
		} else {
		xml$taxa = c(xml$taxa, taxon_XMLs)
		xml$list_of_empty_taxa = list_of_empty_taxa
		list_of_empty_taxa = xml$list_of_empty_taxa
		return(xml)
		} # END if (is.null(xml))

	cat("\n\n")
	stop("ERROR in make_taxa_groups(): shouldn't get here")
	} # END make_taxa_groups <- function(taxa_df, OTUs=NULL, xml=NULL)

make_cladePrior_XMLs <- function(nodes_df, tree_name="shared_tree", xml=NULL, list_of_empty_taxa=NULL)
	{
	defaults='
	tree_name="shared_tree"
	xml=NULL
	list_of_empty_taxa=NULL
	'
	
	# Subset nodes_df
	nodes_df = nodes_df[nodes_df$use != "no",]
	
	# Remove taxa groups that are empty (have no species in them)
	if (!is.null(list_of_empty_taxa))
		{
		taxa_names = nodes_df$Taxon
		TF = taxa_names %in% list_of_empty_taxa
		nodes_df = nodes_df[TF==FALSE,]
		} # END if (is.null(list_of_empty_taxa))
	
	# Make the clade priors
	priordists_XML = list()
	for (i in 1:nrow(nodes_df))
		{
		nodeline = nodes_df[i,]
		priordist_XML_list = nodeline_to_cladePrior(nodeline, tree_name=tree_name)
		priordists_XML = c(priordists_XML, priordist_XML_list)
		priordists_XML
		}
	priordists_XML

	# Otherwise...
	txt = paste(" Priors on node dates and/or clade monophyly ", sep="")
	XML_comment = xmlCommentNode(txt)
	txt2 = paste(" (some might be unconstrained, just for logging purposes) ", sep="")
	XML_comment2 = xmlCommentNode(txt2)
	priordists_XML = c(list(bl()), list(XML_comment), list(XML_comment2), priordists_XML)
	
	# Just return the priordists_XML, if no master xml given
	if (is.null(xml))
		{
		return(priordists_XML)
		}
	
	xml$priors = c(xml$priors, priordists_XML)
	
	return(xml)
	} # END make_cladePrior_XMLs <- function(nodes_df, tree_name="shared_tree", xml=NULL)


make_cladePrior_logs <- function(nodes_df, xml=NULL, tree_name="shared_tree", trace_or_screen="both", list_of_empty_taxa=NULL)
	{
	defaults='
	xml = NULL
	trace_or_screen = "both"
	'
	
	# Subset nodes_df
	nodes_df = nodes_df[nodes_df$use != "no",]

	# Remove ones that are empty
	if (!is.null(list_of_empty_taxa))
		{
		taxa_names = nodes_df$Taxon
		TF = taxa_names %in% list_of_empty_taxa
		nodes_df = nodes_df[TF==FALSE,]
		} # END if (is.null(list_of_empty_taxa))
	
	tree_name_idref = paste0("@", tree_name)
	
	traceLogs_XML = list()
	screenLogs_XML = list()
	for (i in 1:nrow(nodes_df))
		{
		nodeline = nodes_df[i,]
		taxon_name = nodeline$Taxon
		prior_id = paste(taxon_name, "prior", sep="_")
		
		log_XML = xmlNode(name="log", attrs=list(idref=prior_id))
		
		if ( (trace_or_screen == "both") || (trace_or_screen == "trace") )
			{
			traceLogs_XML = c(traceLogs_XML, list(log_XML))
			}
		if ( (trace_or_screen == "both") || (trace_or_screen == "trace") )
			{
			screenLogs_XML = c(screenLogs_XML, list(log_XML))
			}		
		} # END for (i in 1:nrow(nodes_df))

	txt = " Reporting node dates/monophyly to trace log "
	XML_comment = xmlCommentNode(txt)
	traceLogs_XML = c(list(XML_comment), traceLogs_XML)

	txt = " Reporting node dates/monophyly to screen log "
	XML_comment = xmlCommentNode(txt)
	screenLogs_XML = c(list(XML_comment), screenLogs_XML)
	
	if (is.null(xml) == TRUE)
		{
		if (trace_or_screen == "both")
			{
			logs_XML = c(traceLogs_XML, screenLogs_XML)
			return(logs_XML)
			}
		if (trace_or_screen == "trace")
			{
			return(traceLogs_XML)
			}
		if (trace_or_screen == "screen")
			{
			return(screenLogs_XML)
			}
		} # END if (is.null(xml) == TRUE)


	if (is.null(xml) == FALSE)
		{
		if (trace_or_screen == "both")
			{
			xml$tracelog = c(xml$tracelog, traceLogs_XML)
			xml$screenlog = c(xml$screenlog, screenLogs_XML)
			return(xml)
			}
		if (trace_or_screen == "trace")
			{
			xml$tracelog = c(xml$tracelog, traceLogs_XML)
			return(xml)
			}
		if (trace_or_screen == "screen")
			{
			xml$screenlog = c(xml$screenlog, screenLogs_XML)
			return(xml)
			}
		} # END if (is.null(xml) == FALSE)
	} # END make_cladePrior_logs <- function(nodes_df, xml=NULL, tree_name="shared_tree", trace_or_screen="both")

# Take a line from the "nodes" worksheet, and
# make an XML prior
nodeline_to_cladePrior <- function(nodeline, tree_name="shared_tree")
	{
	param2_XML = NULL
	
	tree_name_idref = paste0("@", tree_name)
	
	taxon_name = nodeline$Taxon
	taxonset_node = xmlNode(name="taxonset", attrs=list(idref=taxon_name))
	
	# Break if you don't want to use this nodeline
	if (tolower(nodeline$use) == "no")
		{
		return(NULL)
		}
	
	# Convert to normal: default "no"
	if (isblank_TF(nodeline$convert_to_normal) == TRUE)
		{
		nodeline$convert_to_normal = "no"
		}

	# Default tipsOnly is "no"
	tipsOnly = nodeline$tipsOnly
	tipsOnly[isblank_TF(tipsOnly)] = "no"
	if (tipsOnly == "yes")
		{
		tipsOnly_TF = "true"
		} else {
		tipsOnly_TF = "false"
		} # END if (tipsOnly == "yes")
			
	# Check if node age prior is desired & filled in
	distrib_XML = NULL
	TF1 = tolower(nodeline$make_age_prior) == "yes"
	if (TF1 == TRUE)
		{
		items = unname(unlist(c(nodeline[, c("distribution","param1","param2")])))
		itemtxt = c("distribution","param1","param2")
		if (items[1] == "exponential")
			{
			items[3] = "not needed"
			} # END if (items[1] == "exponential")
		
		TF2 = isblank_TF(items)
		
		
		if (sum(TF2) > 0)
			{
			not_specified_txt = paste(itemtxt[TF2], sep=", ")
		
			errortxt1 = paste0("\n\nERROR in nodeline_to_cladePrior(): for taxon '", taxon_name, "', make_age_prior is 'yes',\nbut the following are not specified:\n", not_specified_txt, "\n")
			cat(errortxt1)
		
			stop(errortxt1)
			} # END if (sum(TF2) > 0)
		distrib = nodeline$distribution
		param1 = nodeline$param1
		param2 = nodeline$param2
		meanInRealSpace = nodeline$meanInRealSpace
		offset = nodeline$offset
		node_or_stem = nodeline$node_or_stem
		tipsOnly = nodeline$tipsOnly
		tipsOnly[isblank_TF(tipsOnly)] = "no"
		
		# Fill in blanks with defaults
		if (isblank_TF(meanInRealSpace))
			{
			meanInRealSpace = "no"
			}
		if (isblank_TF(offset))
			{
			offset = 0.0
			}
		if (isblank_TF(node_or_stem))
			{
			node_or_stem = "node"
			}
		if (isblank_TF(tipsOnly))
			{
			tipsOnly = "no"
			}
		if (tipsOnly == "yes")
			{
			tipsOnly_TF = "true"
			} else {
			tipsOnly_TF = "false"
			} # END if (tipsOnly == "yes")
			
		# Specify the distribution
		# Error check
		allowed_distributions = c("uniform", "normal", "lognormal", "exponential")
		if ( (tolower(distrib) %in% allowed_distributions ) == FALSE )
			{
			errortxt1 = paste0("\n\nERROR in nodeline_to_cladePrior(): prior distribution '", distrib, "' has not yet been implemented here.\n")
			cat(errortxt1)
			cat("\nAllowed distributions:\n\n")
			print(allowed_distributions)
			cat("\n\n")
			stop(errortxt1)
			}
		
		if (tolower(distrib) == "uniform")
			{
			distrib = "Uniform"
			param1_txt = "lower"
			param2_txt = "upper"
			normMeanval = ((param1) + (param2)) / 2
			normSDval = 0.01 * normMeanval

			if (offset > 0)
				{
				errortxt = paste0("\n\nSTOP ERROR in nodeline_to_cladePrior(): ", taxon_name, " has a Uniform distribution but the offset is >0. offset=", offset, "\n\n")
				cat(errortxt)
				stop(errortxt)
				} # END if (offset > 0)
			} # END if (tolower(distrib) == "uniform")
		if (tolower(distrib) == "normal")
			{
			distrib = "Normal"
			param1_txt = "mean"
			param2_txt = "stdev"
			if (offset > 0)
				{
				errortxt = paste0("\n\nSTOP ERROR in nodeline_to_cladePrior(): ", taxon_name, " has a Normal distribution but the offset is >0. offset=", offset, "\n\n")
				cat(errortxt)
				stop(errortxt)
				} # END if (offset > 0)
			} # END if (tolower(distrib) == "normal")
		if (tolower(distrib) == "lognormal")
			{
			distrib = "LogNormal"
			param1_txt = "mean"
			param2_txt = "stdev"
			if (meanInRealSpace == "yes")
				{
				meanInRealSpace = "true"
				} else {
				meanInRealSpace = "false"
				}
			if (meanInRealSpace == "true")
				{
				normMeanval = param1 + offset
				normSDval = 0.01 * normMeanval
				} else {
				normMeanval = exp(param1) + offset
				normSDval = 0.01 * normMeanval
				}
			} # END if (tolower(distrib) == "lognormal")
		if (tolower(distrib) == "exponential")
			{
			distrib = "Exponential"
			param1_txt = "mean"
			param2_txt = "stdev"
			if (meanInRealSpace == "yes")
				{
				meanInRealSpace = "true"
				} else {
				meanInRealSpace = "false"
				}
			if (meanInRealSpace == "true")
				{
				normMeanval = param1 + offset
				normSDval = 0.01 * normMeanval
				} else {
				normMeanval = exp(param1) + offset
				normSDval = 0.01 * normMeanval
				}
			} # END if (tolower(distrib) == "exponential")
			
		
		# xmlNodes for the parameters
		if ((distrib == "Uniform") && (nodeline$convert_to_normal != "yes"))
			{
			node_id = paste(taxon_name, "param", distrib, sep="_")
			param1_XML = xmlNode(name=distrib, attrs=list(id=node_id, name="distr", lower=param1, upper=param2))
			param2_XML = NULL
			}
		if ((distrib == "Normal") && (nodeline$convert_to_normal != "yes"))
			{
			node_id = paste(taxon_name, "param", distrib, param1_txt, sep="_")
			param1_XML = xmlNode(name="parameter", param1, attrs=list(id=node_id, name="mean", estimate="false"))
			node_id = paste(taxon_name, "param", distrib, param2_txt, sep="_")
			param2_XML = xmlNode(name="parameter", param2, attrs=list(id=node_id, name="sigma", estimate="false"))
			}
		if ((distrib == "LogNormal") && (nodeline$convert_to_normal != "yes"))
			{
			node_id = paste(taxon_name, "param", distrib, param1_txt, sep="_")
			param1_XML = xmlNode(name="parameter", param1, attrs=list(id=node_id, name="M", estimate="false"))
			node_id = paste(taxon_name, "param", distrib, param2_txt, sep="_")
			param2_XML = xmlNode(name="parameter", param2, attrs=list(id=node_id, name="S", estimate="false"))
			}
		if ((distrib == "Exponential") && (nodeline$convert_to_normal != "yes"))
			{
			node_id = paste(taxon_name, "param", distrib, param1_txt, sep="_")
			param1_XML = xmlNode(name="parameter", param1, attrs=list(id=node_id, name="mean", estimate="false"))
			}
		# If we are converting to normal
		if (nodeline$convert_to_normal == "yes")
			{
			distrib = "Normal"
			node_id = paste("ConvertedForStartTree", taxon_name, "param", distrib, param1_txt, sep="_")
			param1_XML = xmlNode(name="parameter", normMeanval, attrs=list(id=node_id, name="mean", estimate="false"))
			node_id = paste("ConvertedForStartTree", taxon_name, "param", distrib, param2_txt, sep="_")
			param2_XML = xmlNode(name="parameter", normSDval, attrs=list(id=node_id, name="sigma", estimate="false"))
			}


		# xmlNode for the distribution
		node_id = paste(taxon_name, "prior", distrib, sep="_")
		if (!is.null(param2_XML))
			{
			children = list(param1_XML, param2_XML)
			} else {
			children = list(param1_XML)
			}
		#distrib_XML = xmlNode(name=distrib, attrs=list(id=node_id, name="distribution", x="@shared_tree"), .children=children)
		if (distrib == "LogNormal") 
			{
			distrib_XML = xmlNode(name=distrib, attrs=list(id=node_id, name="distr", meanInRealSpace=meanInRealSpace, offset=offset), .children=children)
			}
		if (distrib == "Exponential")
			{
			distrib_XML = xmlNode(name=distrib, attrs=list(id=node_id, name="distr", offset=offset), .children=children)
			}
		if (distrib == "Normal")
			{
			distrib_XML = xmlNode(name=distrib, attrs=list(id=node_id, name="distr"), .children=children)
			}
		if (distrib == "Uniform")
			{
			distrib_XML = param1_XML
			}			
		distrib_XML
		} # END if (TF1 == TRUE)
	
	# xmlNode for the prior
	if (nodeline$convert_to_normal == "yes")
		{
		#prior_id = paste("ConvertedForStartTree", taxon_name, "prior", sep="_")
		prior_id = paste(taxon_name, "prior", sep="_")
		} else {
		prior_id = paste(taxon_name, "prior", sep="_")
		}
	mono_txt = "false"
	if (nodeline$mono == "yes")
		{
		mono_txt = "true"
		}
	useStem_txt = "false"
	if (nodeline$node_or_stem == "stem")
		{
		useStem_txt = "true"
		}
	
	# Date and monophyly constrained
	if ((is.null(distrib_XML) == FALSE) && (nodeline$mono == "yes") ) 
		{
		txt = paste(" NOTE: ConvertedForStartTree...Prior distribution on the date of monophyletic taxon: ", taxon_name, " ", sep="")
		XML_comment = xmlCommentNode(txt)

		priordist_XML = xmlNode(name="distribution", attrs=list(id=prior_id, monophyletic=mono_txt, spec="beast.math.distributions.MRCAPrior", tree=tree_name_idref, useOriginate=useStem_txt, tipsonly=tipsOnly_TF), .children=c(list(taxonset_node), list(distrib_XML)))
	
		priordist_XML_list = list(bl(), XML_comment, priordist_XML)
		}
	
	# Date-only constrained, monophyly NOT required
	if ((is.null(distrib_XML) == FALSE) && (nodeline$mono == "no") ) 
		{
		txt = paste(" Prior distribution on the date of taxon: ", taxon_name, " (not constrained to be monophyletic) ", sep="")
		XML_comment = xmlCommentNode(txt)

		priordist_XML = xmlNode(name="distribution", attrs=list(id=prior_id, monophyletic=mono_txt, spec="beast.math.distributions.MRCAPrior", tree=tree_name_idref, useOriginate=useStem_txt, tipsonly=tipsOnly_TF), .children=c(list(taxonset_node), list(distrib_XML)))
	
		priordist_XML_list = list(bl(), XML_comment, priordist_XML)
		}
	
	
	# Distribution on the date, monophyly required
	if ((is.null(distrib_XML) == TRUE) && (nodeline$mono == "yes") ) 
		{
		txt = paste(" Distribution on the date (no prior) of monophyletic taxon: ", taxon_name, " ", sep="")
		XML_comment = xmlCommentNode(txt)

		priordist_XML = xmlNode(name="distribution", attrs=list(id=prior_id, monophyletic=mono_txt, spec="beast.math.distributions.MRCAPrior", tree=tree_name_idref, useOriginate=useStem_txt, tipsonly=tipsOnly_TF), .children=c(list(taxonset_node)))
	
		priordist_XML_list = list(bl(), XML_comment, priordist_XML)
		}			
		
	# Distribution on the date, monophyly NOT required
	if ((is.null(distrib_XML) == TRUE) && (nodeline$mono == "no") ) 
		{
		txt = paste(" Distribution on the date (no prior) of taxon: ", taxon_name, " (not constrained to be monophyletic) ", sep="")
		XML_comment = xmlCommentNode(txt)

		priordist_XML = xmlNode(name="distribution", attrs=list(id=prior_id, monophyletic=mono_txt, spec="beast.math.distributions.MRCAPrior", tree=tree_name_idref, useOriginate=useStem_txt, tipsonly=tipsOnly_TF), .children=c(list(taxonset_node)))
	
		priordist_XML_list = list(bl(), XML_comment, priordist_XML)		
		}
	
	return(priordist_XML_list)
	} # END nodeline_to_cladePrior <- function(nodeline, tree_name="shared_tree")