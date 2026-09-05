
#######################################################
# List of taxa, for a starBEAST analysis
#######################################################
taxa_list_for_starBEAST <- function(xlsfn, sheet="taxonsets")
	{
	taxonsets_df = readWorksheetFromFile(xlsfn, sheet="taxonsets", startRow=15)
	taxonsets_df$use[isblank_TF(taxonsets_df$use)] = "yes"
	keepTF = (taxonsets_df$use != "no")
	taxonsets_df = taxonsets_df[keepTF, ]
	taxa_list = unique(taxonsets_df$speciesName)
	return(taxa_list)
	}


#######################################################
# speciesTreeLogger
#######################################################
make_speciesTreeLogger <- function(fileName="species.trees.txt", logEvery=50000, speciesTree_default="speciesTree", xml=NULL, treemodel_df=NULL)
	{
	# Check for BEASTmasteR-programmed population models
	if (is.null(treemodel_df$popModel))
		{
		treemodel_df$popModel = "constant"
		} # END if (is.null(treemodel_df$popModel))
	if (isblank_TF(treemodel_df$popModel))
		{
		treemodel_df$popModel = "constant"
		} # END if (is.null(treemodel_df$popModel))
	popModel = treemodel_df$popModel

	if (popModel == "constant")
		{
		popModel_prefix = "constantPopModel"
		}
	if (popModel == "analytic")
		{
		popModel_prefix = "constantPopIOModel"
		}
	if ((treemodel_df$popModel == "linear"))
		{
		popModel_prefix = "lwcrPopModel"
		}
	
	popModel_id = paste0(popModel_prefix, "_for_", speciesTree_default)
	popModel_idref = paste0("@", popModel_id)
	popSizes_id = paste0(popModel_prefix, "_popSizes.Species")
	popSizes_idref = paste0("@", popSizes_id)
	popSizesMean_id = paste0(popModel_prefix, "_popSizesMean.Species")
	popSizesShape_id = paste0(popModel_prefix, "_popSizesShape.Species")
	popSizes_idref = paste0("@", popModel_prefix, "_popSizes.Species")
	popSizesMean_idref = paste0("@", popModel_prefix, "_popSizesMean.Species")
	popSizesShape_idref = paste0("@", popModel_prefix, "_popSizesShape.Species")
	popModelBridge_id = "popModelBridge.Species"
	
	
	speciesTree_default_idref = paste0("@", speciesTree_default)
	
	child = xmlNode(name="log", attrs=list(id="SpeciesTreeLoggerX", spec="starbeast2.SpeciesTreeLogger", populationmodel="@constantPopModel.Species", speciesTree=speciesTree_default_idref))
	
	# Convert logEvery to text
	logEvery = as.numeric(logEvery)
	logEvery = sprintf("%0.0f", as.numeric(logEvery))
	
	
	logger_XML = xmlNode(name="logger", attrs=list(id="speciesTreeLogger", fileName=fileName, logEvery=logEvery, mode="tree"), .children=list(child))
	
	if (is.null(xml))
		{
		return(logger_XML)
		} else {
		xml$species
		}
	return(stop("make_speciesTreeLogger() ERROR: Shouldn't get here!"))
	} # END make_speciesTreeLogger <- function(fileName="species.trees.txt", logEvery=50000, speciesTree_default="speciesTree", xml=NULL)

#make_geneTreeLoggers <- function()



#######################################################
# Multispecies coalescent for StarBeast2
#######################################################
make_speciescoalescent_for_starBeast2 <- function(OTUs_df, seqs_df, treemodel_df, tree_name="shared_tree", speciesTree_default="speciesTree", xml=NULL)
	{
	# Check for ploidy column in the sequences column
	if (is.null(seqs_df$ploidy))
		{
		seqs_df$ploidy = rep(2.0, times=nrow(seqs_df))
		} # END if (is.null(seqs_df$ploidy))
	
	
	# Check for BEASTmasteR-programmed population models
	if (is.null(treemodel_df$popModel))
		{
		treemodel_df$popModel = "constant"
		} # END if (is.null(treemodel_df$popModel))
	if (isblank_TF(treemodel_df$popModel))
		{
		treemodel_df$popModel = "constant"
		} # END if (is.null(treemodel_df$popModel))
	popModel = treemodel_df$popModel

	if (popModel == "constant")
		{
		popModel_prefix = "constantPopModel"
		}
	if (popModel == "analytic")
		{
		popModel_prefix = "constantPopIOModel"
		}
	if ((treemodel_df$popModel == "linear"))
		{
		popModel_prefix = "lwcrPopModel"
		}
	
	popModel_id = paste0(popModel_prefix, "_for_", speciesTree_default)
	popModel_idref = paste0("@", popModel_id)
	popSizes_id = paste0(popModel_prefix, "_popSizes.Species")
	popSizes_idref = paste0("@", popSizes_id)
	popSizesMean_id = paste0(popModel_prefix, "_popSizesMean.Species")
	popSizesShape_id = paste0(popModel_prefix, "_popSizesShape.Species")
	popSizes_idref = paste0("@", popModel_prefix, "_popSizes.Species")
	popSizesMean_idref = paste0("@", popModel_prefix, "_popSizesMean.Species")
	popSizesShape_idref = paste0("@", popModel_prefix, "_popSizesShape.Species")
	popModelBridge_id = "popModelBridge.Species"

	tip_popSizes_id = paste0("tip_", popModel_prefix, "_popSizes.Species")
	tip_popSizes_idref = paste0("@", popSizes_id)
	tip_popShape_id = paste0("Shape_for_", popSizes_id)
	
	top_popSizes_id = paste0("top_", popModel_prefix, "_popSizes.Species")
	top_popSizes_idref = paste0("@", popSizes_id)
	lwcrPopScale_id = paste0(popModel_prefix, "_popScale.Species")
	lwcrPopScale_idref = paste0("@", lwcrPopScale_id)
		
	
	OTUs_df = OTUs_df[OTUs_df$use != "no"]
	numspecies = nrow(OTUs_df)
	popSizes_dimension = numspecies + (numspecies-2)

	
	# Define parameters, put in stateNodes
	if (popModel == "constant")
		{
		# No dimension needed (should just be 1 popSize for everybody)
		stateNode1_XML = xmlNode(name="parameter", 1.0, attrs=list(id=popSizes_id, lower="0.0", name="stateNode"))
		} else {
		stateNode1_XML = xmlNode(name="parameter", 1.0, attrs=list(id=popSizes_id, dimension=numspecies, lower="0.0", name="stateNode"))
		} # END if (popModel == "constant")
		
	# Mean on population size
	stateNode2_XML = xmlNode(name="parameter", 1.0, attrs=list(id=popSizesMean_id, lower="0.0", name="stateNode"))
	stateNode3_XML = NULL
	
	if ((treemodel_df$popModel == "linear"))
		{
		stateNode3_XML = xmlNode(name="parameter", 1.0, attrs=list(id=lwcrPopScale_id, lower="0.0", name="stateNode"))
		}
	if ((treemodel_df$popModel == "analytic"))
		{
		# populationShape
		#stateNode3_XML = xmlNode(name="parameter", 3.0, attrs=list(id=popSizesShape_id, lower="0.0", estimate="false", name="stateNode"))
		}
	
	# Define priors
	if ((treemodel_df$popModel == "old2.4.2"))
		{
		child1 = xmlNode(name="parameter", 3.0, attrs=list(id=popSizesShape_id, estimate="false", lower="0.0", name="alpha"))
		child2 = xmlNode(name="distr", attrs=list(id="AltInverseGamma.0", spec="starbeast2.AltInverseGamma", mean=popSizesMean_idref), .children=list(child1))
		popSizesShape_prior_id = paste0("prior_on_", popSizesShape_id)
		prior1_XML = xmlNode(name="prior", attrs=list(id=popSizesShape_prior_id, name="distribution", x=popSizes_idref), .children=list(child2))
	
		popSizesMeanPrior_id = paste0("OneOnX_prior_for_", popSizesMean_id)
		child2 = xmlNode(name="OneOnX", attrs=list(id="OneOnX_distrib_for_popSizesMeanPrior", name="distr"))
		prior2_XML = xmlNode(name="prior", attrs=list(id=popSizesMeanPrior_id, name="distribution", x=popSizesMean_idref), .children=list(child2))
		
		prior3_XML = NULL
		} # END if ((treemodel_df$popModel == "constant"))


	if ((treemodel_df$popModel == "constant"))
		{
		child1 = xmlNode(name="parameter", 2.0, attrs=list(id=popSizesShape_id, estimate="false", lower="0.0", name="alpha"))
		child2 = xmlNode(name="Gamma", attrs=list(id="Gamma.0", beta=popSizesMean_idref, mode="ShapeMean", name="distr", spec="beast.math.distributions.Gamma"), .children=list(child1))
		popSizesShape_prior_id = paste0("prior_on_", popSizesShape_id)
		prior1_XML = xmlNode(name="prior", attrs=list(id=popSizesShape_prior_id, name="distribution", x=popSizes_idref), .children=list(child2))

		popSizesMeanPrior_id = paste0("OneOnX_prior_for_", popSizesMean_id)
		child2 = xmlNode(name="OneOnX", attrs=list(id="OneOnX_distrib_for_popSizesMeanPrior", name="distr"))
		prior2_XML = xmlNode(name="prior", attrs=list(id=popSizesMeanPrior_id, name="distribution", x=popSizesMean_idref), .children=list(child2))
			
# 		popSizesMeanPrior_id = paste0("GammaShapePrior_for_", popSizesShape_id)
# 		child1 = xmlNode(name="parameter", 1.0, attrs=list(id="FixedMean_for_GammaShapePrior", estimate="false", lower="0.0", name="mean"))
# 		child2 = xmlNode(name="Exponential", attrs=list(id="Exponential_distrib_for_Gamma_distrib_popSizesShape", name="distr"), .children=list(child1))
# 		prior2_XML = xmlNode(name="prior", attrs=list(id=popSizesMeanPrior_id, name="distribution", x=popSizesMean_idref), .children=list(child2))
		
		prior3_XML = NULL
		} # END if ((treemodel_df$popModel == "constant"))



	if ((treemodel_df$popModel == "analytic"))
		{
		child1 = xmlNode(name="parameter", 3.0, attrs=list(id=popSizesShape_id, estimate="false", lower="0.0", name="populationShape"))

		popSizesShapePrior_id = paste0("Uniform_prior_for_", popSizesShape_id)
		child2 = xmlNode(name="Uniform", attrs=list(id="Uniform_distrib_for_popSizesShapePrior", name="distr", lower=0, upper=10000))

		prior1_XML = xmlNode(name="prior", attrs=list(id=popSizesShape_prior_id, name="distribution"), .children=list(child1, child2))

		
		OneOnX_child_XML = xmlNode(name="OneOnX", attrs=list(id="OneOnX_prior_for_popSizesMean", name="distr"))
		prior2_XML = xmlNode(name="prior", attrs=list(id=popSizesMeanPrior_id, name="distribution", x=popSizesMean_idref), .children=list(OneOnX_child_XML))
		
		prior3_XML = NULL
		} # END if ((treemodel_df$popModel == "constant"))


	if ((treemodel_df$popModel == "linear"))
		{
		# Tip Population sizes
		child1 = xmlNode(name="parameter", 4.0, attrs=list(id=tip_popShape_id, estimate="false", lower="0.0", name="alpha"))
		child2 = xmlNode(name="Gamma", attrs=list(id="Gamma.01", beta=lwcrPopScale_idref, name="distr"), .children=list(child1))
		prior_on_tip_popSizes_id = paste0("prior_on_", tip_popSizes_id)
		prior1_XML = xmlNode(name="prior", attrs=list(id=prior_on_tip_popSizes_id, name="distribution", x=tip_popSizes_idref), .children=list(child2))

		# Top Population sizes
		child1 = xmlNode(name="parameter", 2.0, attrs=list(id=top_popShape_id, estimate="false", lower="0.0", name="alpha"))
		child2 = xmlNode(name="Gamma", attrs=list(id="Gamma.02", beta=lwcrPopScale_idref, name="distr"), .children=list(child1))
		prior_on_top_popSizes_id = paste0("prior_on_", top_popSizes_id)
		prior2_XML = xmlNode(name="prior", attrs=list(id=prior_on_top_popSizes_id, name="distribution", x=top_popSizes_idref), .children=list(child2))
	
		popScalePrior_id = paste0("OneOnX_prior_for_", lwcrPopScale_id)
		child2 = xmlNode(name="OneOnX", attrs=list(id="OneOnX_distrib_for_popScalePrior", name="distr"))
		prior3_XML = xmlNode(name="prior", attrs=list(id=popScalePrior_id, name="distribution", x=lwcrPopScale_idref), .children=list(child2))
		} # END if ((treemodel_df$popModel == "linear"))	


	# Define operators
	if ((treemodel_df$popModel == "constant"))
		{
		op1_XML = xmlNode(name="operator", attrs=list(id="popSizesSwap.Species", spec="SwapOperator", parameter=popSizes_idref, weight="3.0"))

		op2_XML = xmlNode(name="operator", attrs=list(id="popSizesScale.Species", spec="ScaleOperator", parameter=popSizes_idref, scaleFactor="0.5", weight="3.0"))

		op3_XML = xmlNode(name="operator", attrs=list(id="popSizesMeanScale.Species", spec="ScaleOperator", parameter=popSizesMean_idref, scaleFactor="0.75", weight="1.0"))
		}

	# Define operators
	if ((treemodel_df$popModel == "analytic"))
		{
		op1_XML = NULL

		op2_XML = NULL

		op3_XML = xmlNode(name="operator", attrs=list(id="popSizesMeanScale.Species", spec="ScaleOperator", parameter=popSizesMean_idref, scaleFactor="0.75", weight="1.0"))
		}

	if ((treemodel_df$popModel == "linear"))
		{
		op1_XML = xmlNode(name="operator", attrs=list(id="tipPopSizesScale.Species", spec="ScaleOperator", parameter=tip_popSizes_idref, scaleFactor="0.5", weight="3.0"))

		op2_XML = xmlNode(name="operator", attrs=list(id="topPopSizesScale.Species", spec="ScaleOperator", parameter=top_popSizes_idref, scaleFactor="0.5", weight="3.0"))

		op3_XML = xmlNode(name="operator", attrs=list(id="lwcrPopScaleScale.Species", spec="ScaleOperator", parameter=lwcrPopScale_idref, scaleFactor="0.75", weight="1.0"))
		}

	# Define logs
	if ((treemodel_df$popModel == "constant"))
		{
		log1_xml = xmlNode(name="log", attrs=list(idref=popSizes_id))
		log2_xml = xmlNode(name="log", attrs=list(idref=popSizesMean_id))
		log3_xml = NULL
		}
	if ((treemodel_df$popModel == "analytic"))
		{
		log1_xml = xmlNode(name="log", attrs=list(idref=popSizesMean_id))
		log2_xml = xmlNode(name="log", attrs=list(idref=popSizesShape_id))
		log3_xml = NULL
		}
	if ((treemodel_df$popModel == "linear"))
		{
		log1_xml = xmlNode(name="log", attrs=list(idref=tip_popSizes_id))
		log2_xml = xmlNode(name="log", attrs=list(idref=top_popSizes_id))
		log3_xml = xmlNode(name="log", attrs=list(idref=lwcrPopScale_id))
		}


	# Likelihood of the multispecies coalescence process
	speciesTree_default_idref = paste0("@", speciesTree_default)
	tree_idref = paste0("@", tree_name)

	# Get the geneTreeNames	
	geneTreeNames = unique(seqs_df$geneTreeName)
	rownums = match(x=geneTreeNames, table=seqs_df$geneTreeName)

	# Get the ploidy:
	# 2 for an autosomal locus
	# 0.5 for an mtDNA or Y-chromosome locus
	# For an X-chromosome?
	ploidy_vals = seqs_df$ploidy[rownums]
	ploidy_vals[isblank_TF(ploidy_vals)] = 2.0
	
	geneTree_XML_list = list()
	speciesTree_XML = xmlNode(name="speciesTree", attrs=list(id=speciesTree_default, spec="starbeast2.SpeciesTree", tree=tree_idref))
	geneTree_XML_list = NULL


	
if (popModel == "old2.4.2")
	{
	geneTree_XML_list = c(geneTree_XML_list, list(speciesTree_XML))
	for (i in 1:length(geneTreeNames))
		{
		genetree_idref = paste0("@", geneTreeNames[i])
		# Beast2, 2.4.2
		genetree_speciesCoal_id = paste0("genetree_speciesCoal_", geneTreeNames[i])
		tmpXML = xmlNode(name="geneTree", attrs=list(id=genetree_speciesCoal_id, spec="starbeast2.GeneTree", speciesTree=speciesTree_default_idref, tree=genetree_idref, ploidy=ploidy_vals[i]))
		geneTree_XML_list = c(geneTree_XML_list, list(tmpXML))
		}

	speciesCoal_XML = xmlNode(name="distribution", attrs=list(id="speciescoalescent", spec="starbeast2.MultispeciesCoalescent", populationModel="@constantPopModel.Species"), .children=geneTree_XML_list)
	} # END if (popModel == "old2.4.2")


if (popModel == "analytic")
	{
	# Add population shape to the genetree distribution list
	param_XML = xmlNode(name="distribution", attrs=list(idref=popSizesShape_id))
	geneTree_XML_list = c(param_XML, geneTree_XML_list)
	
	for (i in 1:length(geneTreeNames))
		{
		genetree_idref = paste0("@", geneTreeNames[i])
		# Beast2, 2.4.3
		genetree_speciesCoal_id = paste0("genetree_speciesCoal_", geneTreeNames[i])
		popModel_idref = paste0("@", popModel_id)
		
		popModelBridge_idref = paste0("@", popModelBridge_id)
		tmpXML = xmlNode(name="distribution", attrs=list(id=genetree_speciesCoal_id, spec="starbeast2.GeneTree", speciesTree=speciesTree_default_idref, tree=genetree_idref, populationModel=popModelBridge_idref, ploidy=ploidy_vals[i]))

		geneTree_XML_list = c(geneTree_XML_list, list(tmpXML))
		}

	speciesCoal_XML = xmlNode(name="distribution", attrs=list(id="speciescoalescent", spec="starbeast2.MultispeciesCoalescent"), .children=geneTree_XML_list, populationMean=popSizesMean_idref)
	} # END if (popModel == "analytic")


if ( (popModel == "constant") || (popModel == "linear"))
	{
	for (i in 1:length(geneTreeNames))
		{
		genetree_idref = paste0("@", geneTreeNames[i])
		# Beast2, 2.4.3
		genetree_speciesCoal_id = paste0("genetree_speciesCoal_", geneTreeNames[i])
		popModel_idref = paste0("@", popModel_id)
		
		popModelBridge_idref = paste0("@", popModelBridge_id)
		tmpXML = xmlNode(name="distribution", attrs=list(id=genetree_speciesCoal_id, spec="starbeast2.GeneTree", speciesTree=speciesTree_default_idref, tree=genetree_idref, populationModel=popModelBridge_idref, ploidy=ploidy_vals[i]))

		geneTree_XML_list = c(geneTree_XML_list, list(tmpXML))
		}

	speciesCoal_XML = xmlNode(name="distribution", attrs=list(id="speciescoalescent", spec="starbeast2.MultispeciesCoalescent"), .children=geneTree_XML_list)
	} # END if (popModel == "constant") || (popModel == "linear"))




	speciesCoal_XMLs = list(bl(), xmlCommentNode(" Likelihood of the multispecies coalescent "), speciesCoal_XML)
	
	# Log the speciescoalescent
	xml1 = xmlCommentNode(" Likelihood of the multispecies coalescent ")
	xml2 = xmlNode(name="log", attrs=list(idref="speciescoalescent"))
	log3_xml = cl(xml1, xml2)
	
	# Return
	if (is.null(xml))
		{
		res = NULL
		res$state = cl(res$state, stateNode1_XML, stateNode2_XML, stateNode3_XML)
		res$priors = cl(res$priors, prior1_XML, prior2_XML, prior3_XML)
		res$likes = cl(res$likes, speciesCoal_XMLs)
		res$operators = cl(res$operators, op1_XML, op2_XML, op3_XML)
		res$tracelog = cl(res$tracelog, log1_xml, log2_xml, log3_xml)
		return(res)
		} else {
		xml$state = cl(xml$state, stateNode1_XML, stateNode2_XML, stateNode3_XML)
		xml$priors = cl(xml$priors, prior1_XML, prior2_XML, prior3_XML)
		xml$likes = cl(xml$likes, speciesCoal_XMLs)
		xml$operators = cl(xml$operators, op1_XML, op2_XML, op3_XML)
		xml$tracelog = cl(xml$tracelog, log1_xml, log2_xml, log3_xml)
		return(xml)
		}
	return(stop("make_popSizes_for_starBeast2() ERROR: Shouldn't get here!"))
	} # END make_popSizes_for_starBeast2 <- function(OTUs_df, tree_name="shared_tree", speciesTree_default="speciesTree", xml=NULL)


#######################################################
# Tree model and starting tree
#######################################################
make_starting_tree_for_starBeast2 <- function(treemodel_df, taxonset_XML, xlsfn=xlsfn, speciesTree_default="shared_tree", min_brlen=0.001, xml=NULL)
	{
	defaults='
	
	'
	
	# Setup / error fixes
	if (is.null(treemodel_df$use))
		{
		treemodel_df$use = rep("yes", times=nrow(treemodel_df))
		}
	treemodel_df$use[isblank_TF(treemodel_df$use) == TRUE] = "yes"
	treemodel_df$use = tolower(treemodel_df$use)
	treemodel_df = treemodel_df[treemodel_df$use == "yes"]
	
	tree_name = treemodel_df$speciesTreeName
	tree_name_idref = paste0("@", tree_name)
	
	# Check for BEASTmasteR-programmed population models
	if (is.null(treemodel_df$popModel))
		{
		treemodel_df$popModel = "constant"
		} # END if (is.null(treemodel_df$popModel))
	if (isblank_TF(treemodel_df$popModel))
		{
		treemodel_df$popModel = "constant"
		} # END if (is.null(treemodel_df$popModel))
	popModel = treemodel_df$popModel


	# file_type: Newick or NEXUS
	if (is.null(treemodel_df$file_type))
		{
		treemodel_df$file_type = rep("newick", times=nrow(treemodel_df))
		}
	treemodel_df$file_type[isblank_TF(treemodel_df$file_type) == TRUE] = "newick"
	treemodel_df$file_type = tolower(treemodel_df$file_type)
	file_type = treemodel_df$file_type[1]
	
	if ((file_type %in% c("newick", "nexus")) == FALSE)
		{
		txt = paste0("STOP ERROR in make_starting_tree_for_starBeast2(): The starting tree type must be 'newick' or 'nexus', in treemodel_df$file_type.  Blank means 'newick'. However, you have: '", file_type, "'. Fix this in the Excel settings file, and re-run.")
		
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		}
	
	
	
	###################################################
	# Make the tree stateNode for the species tree
	###################################################

	# Random starting tree
	random = FALSE
	if (isblank_TF(treemodel_df$starting_tree[1]) || treemodel_df$starting_tree[1]=="random" || treemodel_df$starting_tree[1]=="upgma")
		{
		# Make the tree stateNode for the species tree
		if ((popModel == "constant") || (popModel == "linear")  || (popModel == "analytic"))
			{
			speciesTree_statenode_XML = make_speciesTree_statenode2(tree_name=tree_name, taxonset_XML)
			} else {
			speciesTree_statenode_XML = make_speciesTree_statenode(tree_name=tree_name, taxonset_XML)
			}

		random = TRUE
		} # END random startBeast2 starting tree
	
	# Constructed starting tree
	if (treemodel_df$starting_tree[1]=="construct_starting_tree")
		{
		random = FALSE
		constructed_starting_tree_fn = "constructed_starting_tree.newick"
		starting_tree = construct_starting_tree(xlsfn=xlsfn, min_brlen=min_brlen, outfn=constructed_starting_tree_fn)
		#trstr = write_tree_noSci(starting_tree, file="")
		treemodel_df$starting_tree[1] = constructed_starting_tree_fn
		} # END random startBeast2 starting tree
	
	# Read in the newick tree
	if ((random == FALSE) && (treemodel_df$starting_tree[1] != "construct_starting_tree"))
		{
		if (file_type == "newick")
			{
			trstr = write_tree_noSci(read.tree(file=treemodel_df$starting_tree[1]), file="")
			} 
		if (file_type == "nexus")
			{
			trstr = write_tree_noSci(read.nexus(file=treemodel_df$starting_tree[1]), file="")
			} 
		
		# Make the tree stateNode for the species tree
		if ((popModel == "constant") || (popModel == "linear")  || (popModel == "analytic"))
			{
			speciesTree_statenode_XML = make_speciesTree_statenode2(tree_name=tree_name, taxonset_XML)
			} else {
			speciesTree_statenode_XML = make_speciesTree_statenode(tree_name=tree_name, taxonset_XML, starting_tree_string=trstr)
			}
		}

	
	
	if (is.null(xml) == TRUE)
		{
		res = NULL
		res$state = speciesTree_statenode_XML
		res$starting_tree = speciesTree_starting_XML
		return(res)
		} else {
		xml$state = cl(xml$state, speciesTree_statenode_XML)
		xml$starting_tree = cl(xml$starting_tree, speciesTree_starting_XML)
		return(xml)
		}
	
	return(stop("ERROR in make_starting_tree_for_starBeast2(): shouldn't get here."))
	} # END make_starting_tree_for_starBeast2 <- function(treemodel_df, taxonset_XML, xlsfn=xlsfn, min_brlen=0.001)


make_BD_model_for_starBeast2 <- function(treemodel_df, seqs_df, taxonset_XML, tree_name="shared_tree", speciesTree_default="speciesTree", birthRate_id="netDiversificationRate", deathRate_id="ExtinctionFraction", logEvery=50000, xml=NULL, nodes_df=NULL, min_brlen=0.001)
	{
	defaults='
	tree_name="shared_tree"
	birthRate_id="netDiversificationRate"
	xml=NULL
	'

	# Check for BEASTmasteR-programmed population models
	if (is.null(treemodel_df$popModel))
		{
		treemodel_df$popModel = "constant"
		} # END if (is.null(treemodel_df$popModel))
	if (isblank_TF(treemodel_df$popModel))
		{
		treemodel_df$popModel = "constant"
		} # END if (is.null(treemodel_df$popModel))
	popModel = treemodel_df$popModel

	if (popModel == "constant")
		{
		popModel_prefix = "constantPopModel"
		}
	if (popModel == "analytic")
		{
		popModel_prefix = "constantPopIOModel"
		}
	if ((treemodel_df$popModel == "linear"))
		{
		popModel_prefix = "lwcrPopModel"
		}
	
	speciesTree_default_idref = paste0("@", speciesTree_default)
	
	popModel_id = paste0(popModel_prefix, "_for_", speciesTree_default)
	popModel_idref = paste0("@", popModel_id)
	popSizes_id = paste0(popModel_prefix, "_popSizes.Species")
	popSizes_idref = paste0("@", popSizes_id)
	tip_popSizes_id = paste0("tip_", popModel_prefix, "_popSizes.Species")
	tip_popSizes_idref = paste0("@", popSizes_id)
	top_popSizes_id = paste0("top_", popModel_prefix, "_popSizes.Species")
	top_popSizes_idref = paste0("@", popSizes_id)
	lwcrPopScale_id = paste0(popModel_prefix, "_popScale.Species")
	lwcrPopScale_idref = paste0("@", lwcrPopScale_id)
	
	popSizesMean_id = paste0(popModel_prefix, "_popSizesMean.Species")
	popSizesShape_id = paste0(popModel_prefix, "_popSizesShape.Species")
	popSizes_idref = paste0("@", popModel_prefix, "_popSizes.Species")
	popSizesMean_idref = paste0("@", popModel_prefix, "_popSizesMean.Species")
	popSizesShape_idref = paste0("@", popModel_prefix, "_popSizesShape.Species")
	popModelBridge_id = "popModelBridge.Species"



	# file_type: Newick or NEXUS
	if (is.null(treemodel_df$file_type))
		{
		treemodel_df$file_type = rep("newick", times=nrow(treemodel_df))
		}
	treemodel_df$file_type[isblank_TF(treemodel_df$file_type) == TRUE] = "newick"
	treemodel_df$file_type = tolower(treemodel_df$file_type)
	file_type = treemodel_df$file_type[1]
	
	if ((file_type %in% c("newick", "nexus")) == FALSE)
		{
		txt = paste0("STOP ERROR in make_starting_tree_for_starBeast2(): The starting tree type must be 'newick' or 'nexus', in treemodel_df$file_type.  Blank means 'newick'. However, you have: '", file_type, "'. Fix this in the Excel settings file, and re-run.")
		
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		}


	
	tree_name_idref = paste0("@", tree_name)
	
	# Convert logEvery to text
	logEvery = as.numeric(logEvery)
	logEvery = sprintf("%0.0f", as.numeric(logEvery))

	
	# Error check
	treemodel_df$use[isblank_TF(treemodel_df$use)==TRUE] = "yes"
	useTF = treemodel_df$use == "yes"
	treemodel_df = treemodel_df[useTF,]
	if (nrow(treemodel_df) > 1)
		{
		txt = paste0("STOP ERROR in make_BD_model_for_starBeast2(): this function assumes you are doing a starBeast2 analysis. However, starBeast2 assumes a simple birth-death model (constant rate, not skyline models), so you should only have one row in treemodel. However, you have ", nrow(treemodel_df), " rows with 'use' as 'yes' or blank.  Printing treemodel_df:")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		}





	
	###################################################
	# Make the tree stateNode for the species tree
	###################################################

	# Random starting tree
	random = FALSE
	if (isblank_TF(treemodel_df$starting_tree[1]) || treemodel_df$starting_tree[1]=="random" || treemodel_df$starting_tree[1]=="upgma")
		{
		# Make the tree stateNode for the species tree
		if ((popModel == "constant") || (popModel == "linear")  || (popModel == "analytic"))
			{
			speciesTree_statenode_XML = make_speciesTree_statenode2(tree_name=tree_name, taxonset_XML)
			} else {
			speciesTree_statenode_XML = make_speciesTree_statenode(tree_name=tree_name, taxonset_XML)
			}

		random = TRUE
		} # END random startBeast2 starting tree
	
	# Constructed starting tree
	if (treemodel_df$starting_tree[1]=="construct_starting_tree")
		{
		random = FALSE
		constructed_starting_tree_fn = "constructed_starting_tree.newick"
		starting_tree = construct_starting_tree(xlsfn=xlsfn, min_brlen=min_brlen, outfn=constructed_starting_tree_fn)
		#trstr = write_tree_noSci(starting_tree, file="")
		treemodel_df$starting_tree[1] = constructed_starting_tree_fn
		} # END random startBeast2 starting tree
	
	# Read in the newick tree
	if ((random == FALSE) && (treemodel_df$starting_tree[1] != "construct_starting_tree"))
		{
		if (file_type == "newick")
			{
			trstr = write_tree_noSci(read.tree(file=treemodel_df$starting_tree[1]), file="")
			} 
		if (file_type == "nexus")
			{
			trstr = write_tree_noSci(read.nexus(file=treemodel_df$starting_tree[1]), file="")
			} 
		
		# Make the tree stateNode for the species tree
		if ((popModel == "constant") || (popModel == "linear")  || (popModel == "analytic"))
			{
			speciesTree_statenode_XML = make_speciesTree_statenode2(tree_name=tree_name, taxonset_XML)
			} else {
			speciesTree_statenode_XML = make_speciesTree_statenode(tree_name=tree_name, taxonset_XML, starting_tree_string=trstr)
			}
		}






	
	# Check for BEASTmasterR-programmed treemodels in starBeast2 
	treeprior = treemodel_df$treeModel
	allowed_models = c("starBeast2_BD", "starBeast2_CalibratedYule")
	if (treeprior %in% allowed_models)
		{
		pass = 1
		} else {
		txt = paste0("STOP ERROR in make_BD_model_for_starBeast2(): Your choice of tree model is not programmed for make_BD_model_for_starBeast2(). In worksheet 'treemodel', column 'treeModel', you have '", treeprior, "', but the programmed models are just: ", paste0(allowed_models, collapse=", "), ".")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if (treeprior %in% allowed_models)
	
	
	# If CalibratedYule, a nodes_df object is required
	if (treeprior == "starBeast2_CalibratedYule")
		{
		if (is.null(nodes_df) == TRUE)
			{
			txt = "STOP ERROR in make_BD_model_for_starBeast2(): You have specified a CalibratedYule tree prior, with 'starBeast2_CalibratedYule' in worksheet 'treemodel', column 'treeModel', but this requires that you give make_BD_model_for_starBeast2() a 'nodes_df' input. Instead, currently you have the default, nodes_df=NULL"
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			}
		}
	

	

	
	# Initialize the StarBeast2 tree
	geneTreeNames = unique(seqs_df$geneTreeName)
	geneTree_XML_list = list()
	for (i in 1:length(geneTreeNames))
		{
		tmpXML = xmlNode(name="geneTree", attrs=list(idref=geneTreeNames[i]))
		geneTree_XML_list = c(geneTree_XML_list, list(tmpXML))
		}
	
	# Add the population size model
	if (treemodel_df$popModel == "old2.4.2")
		{
		popModel_id = "constantPopModel.Species"
		}
	if (treemodel_df$popModel == "old2.4.2")
		{
		popmodel_XML = xmlNode(name="populationModel", attrs=list(id=popModel_id, spec="starbeast2.ConstantPopulation", populationSizes="@popSizes.Species"))
		}

	
	if (popModel == "constant")
		{
		childModel_XML = xmlNode(name="childModel", attrs=list(id=popModel_id, spec="starbeast2.ConstantPopulations", populationSizes=popSizes_idref, speciesTree=speciesTree_default_idref))
		
		popmodel_XML = xmlNode(name="populationModel", attrs=list(id=popModelBridge_id, spec="starbeast2.PassthroughModel"), .children=list(childModel_XML))
		}

	if (popModel == "analytic")
		{
		childModel_XML = xmlNode(name="childModel", attrs=list(id=popModel_id, spec="starbeast2.DummyModel"))
		
		popmodel_XML = xmlNode(name="populationModel", attrs=list(id=popModelBridge_id, spec="starbeast2.PassthroughModel"), .children=list(childModel_XML))
		}

	if (popModel == "linear")
		{
		# lwcr = LinearWithConstantRoot
		childModel_XML = xmlNode(name="childModel", attrs=list(id=popModel_id, spec="starbeast2.LinearWithConstantRoot", speciesTree=speciesTree_default_idref, tipPopulationSizes=tip_popSizes_idref, topPopulationSizes=top_popSizes_idref))
		
		popmodel_XML = xmlNode(name="populationModel", attrs=list(id=popModelBridge_id, spec="starbeast2.PassthroughModel"), .children=list(childModel_XML))
		}
	
	
	geneTree_XML_list = c(geneTree_XML_list, list(popmodel_XML))	
	
	speciesTree_idref = paste0("@", tree_name)
	birthRate_idref = paste0("@", birthRate_id)
	deathRate_idref = paste0("@", deathRate_id)
	
	# StarBeast2 tree initializer
	if ((popModel == "constant") || (popModel == "linear")  || (popModel == "analytic"))
		{
		if (random == TRUE)
			{
			starBeast2_initializer_XML = xmlNode(name="init", attrs=list(id="SBI", spec="starbeast2.StarBeastInitializer", birthRate=birthRate_idref, estimate="true", speciesTree=speciesTree_idref), .children=geneTree_XML_list)
			}
		if (random == FALSE)
			{
			# Use the Newick string
			starBeast2_initializer_XML = xmlNode(name="init", attrs=list(id="SBI", spec="starbeast2.StarBeastInitializer",  birthRate=birthRate_idref, estimate="true", speciesTree=tree_name_idref, initial=trstr), .children=geneTree_XML_list)
			}
		} else {
		starBeast2_initializer_XML = xmlNode(name="init", attrs=list(id="SBI", spec="starbeast2.StarBeastInitializer",  birthRate=birthRate_idref, estimate="false", speciesTree=tree_name_idref), .children=geneTree_XML_list)
		} # END if ((popModel == "constant") || (popModel == "linear")  || (popModel == "analytic"))
	
	
	# Min and max on birthRate - defaults
	if (isblank_TF(treemodel_df$birthRate_min[1]))
		{
		birthRate_min = 0
		} else {
		birthRate_min = treemodel_df$birthRate_min[1]
		}
	if (isblank_TF(treemodel_df$birthRate_max[1]))
		{
		birthRate_max = 10
		} else {
		birthRate_max = treemodel_df$birthRate_max[1]
		}

	# stateNodes for the netDiversificationRate and ExtinctionFraction
	birthRate_XML = xmlNode(name="parameter", as.numeric(treemodel_df$birthRate_starting_vals), attrs=list(id=birthRate_id, lower=birthRate_min, upper=birthRate_max, name="stateNode"))

	# Prior on the birthRate and deathRate
	birthRate_XML_priorLog = make_generic_XML_prior(dfline=treemodel_df, colname_prefix="birthRate", param_name=birthRate_id, xml=NULL)

	# Operators on the birthRate
	tmp_XML_comment = xmlCommentNode(" Operators on birthRate and deathRate of species tree ")
	
	operator_id = paste0("ScaleOperator_on_", birthRate_id)
	birthRate_op1_XML = xmlNode(name="operator", attrs=list(id=operator_id, spec="ScaleOperator", parameter=birthRate_idref, scaleFactor="0.5", weight="3"))

	# Min and max on deathRate - defaults
	if (isblank_TF(treemodel_df$deathRate_min[1]))
		{
		deathRate_min = 0
		} else {
		deathRate_min = treemodel_df$deathRate_min[1]
		}
	if (isblank_TF(treemodel_df$deathRate_max[1]))
		{
		deathRate_max = 10
		} else {
		deathRate_max = treemodel_df$deathRate_max[1]
		}

	if (treeprior != "starBeast2_CalibratedYule")
		{
		# Max deathRate is 1 (relative extinction rate)
		# (for StarBEAST2 analyses, this is relative extinction rate, between 0 and 1)
		deathRate_max = as.numeric(deathRate_max)
		if (deathRate_max > 1)
			{
			deathRate_max = 1
			}

		deathRate_XML = xmlNode(name="parameter", as.numeric(treemodel_df$deathRate_starting_vals), attrs=list(id=deathRate_id, lower=deathRate_min, upper=deathRate_max, name="stateNode"))
		deathRate_XML_priorLog = make_generic_XML_prior(dfline=treemodel_df, colname_prefix="deathRate", param_name=deathRate_id, xml=NULL)
		# Operators on the deathRate
		operator_id = paste0("ScaleOperator_on_", deathRate_id)
		deathRate_op1_XML = xmlNode(name="operator", attrs=list(id=operator_id, spec="ScaleOperator", parameter=deathRate_idref, scaleFactor="0.5", weight="1"))
		operator_id = paste0("UniformOperator_on_", deathRate_id)
		deathRate_op2_XML = xmlNode(name="operator", attrs=list(id=operator_id, spec="UniformOperator", parameter=deathRate_idref, weight="1"))

		} else  {
		deathRate_XML = NULL
		deathRate_XML_priorLog  = NULL
		deathRate_op1_XML = NULL
		deathRate_op2_XML = NULL		
		} # END if (treeprior != "starBeast2_CalibratedYule")
	
	
	operators_XMLs = list(bl(), tmp_XML_comment, birthRate_op1_XML, deathRate_op1_XML, deathRate_op2_XML)
	
	
	
	if (treeprior == "starBeast2_BD")
		{	
		# Tree Model (tree prior) for species tree
		tree_idref = paste0("@", tree_name)
		treemodel_id = paste0("BirthDeathModel_", tree_name)
	
		tree_prior_XML = xmlNode(name="distribution", attrs=list(id=treemodel_id, spec="beast.evolution.speciation.BirthDeathGernhard08Model", birthDiffRate=birthRate_idref, relativeDeathRate=deathRate_idref, tree=tree_idref))
		} # END if (treeprior == "starBeast2_BD")
	
	if (treeprior == "starBeast2_CalibratedYule")
		{
		# Tree Model (tree prior) for species tree
		tree_idref = paste0("@", tree_name)
		treemodel_id = paste0("CalibratedYuleModel_", tree_name)
		
		# Set up the calibrations, inside the Calibrated Yule setup
		XML_children = make_cladePrior_XMLs(nodes_df, xml=NULL, list_of_empty_taxa=xml$list_of_empty_taxa, CalibratedYule_txt="CalibratedYule")
		
		# Write out to XML
		tree_prior_XML = xmlNode(name="distribution", attrs=list(id=treemodel_id, spec="beast.evolution.speciation.CalibratedYuleModel", birthRate=birthRate_idref, tree=tree_idref, type="full"), .children=XML_children)
		} # END if (treeprior == "starBeast2_CalibratedYule")	
	
	
	#######################################################
	# Logs on the species tree
	#######################################################
	#######################################################
	# treeLog
	#######################################################
	TreeWithMetaDataLogger_id = paste0("TreeWithMetaDataLogger_", tree_name)
	clockModel_name_idref = paste0("@", clockModel_name)
	TreeWithMetaDataLogger_XML = xmlNode(name="log", attrs=list(id=TreeWithMetaDataLogger_id, tree=tree_name_idref, branchratemodel=clockModel_name_idref, substitutions="false", spec="beast.evolution.tree.TreeWithMetaDataLogger") )
	TreeWithMetaDataLogger_XMLs = list(
	xmlCommentNode(" Log the tree and branch rates "),
	TreeWithMetaDataLogger_XML
	)
	
	#######################################################
	# subsLog -- tree of 
	#######################################################
	Subs_TreeWithMetaDataLogger_id = paste0("Subs_TreeWithMetaDataLogger_", tree_name)
	Subs_TreeWithMetaDataLogger_XML = xmlNode(name="log", attrs=list(id=Subs_TreeWithMetaDataLogger_id, tree=tree_name_idref, substitutions="true", spec="beast.evolution.tree.TreeWithMetaDataLogger") )
	Subs_TreeWithMetaDataLogger_XMLs = list(
	xmlCommentNode(" Log the tree and branch rates "),
	Subs_TreeWithMetaDataLogger_XML
	)

	
	
	#######################################################
	# speciesTreeLogger
	#######################################################
	speciesTreeLoggerX_id = paste0("speciesTreeLoggerX_for_", tree_name)
	fn = paste0("species_trees_for_", tree_name, ".trees.txt")
	speciesTree_default_idref = paste0("@", speciesTree_default)
	
	popModel_idref = paste0("@", popModel_id)
	child = xmlNode(name="log", attrs=list(id=speciesTreeLoggerX_id, spec="starbeast2.SpeciesTreeLogger", populationmodel=popModel_idref, speciesTree=speciesTree_default_idref))
	
	speciesTreeLogger_XML = xmlNode(name="logger", attrs=list(id="speciesTreeLogger", fileName=fn, logEvery=logEvery, mode="tree"), .children=list(child))
	speciesTreeLogger_XMLs = list(bl(), xmlCommentNode(" starBeast2 logger for the speciesTree object "), speciesTreeLogger_XML)
	
	

	#######################################################
	# Logs on the gene trees
	#######################################################
	geneTree_treeHeights_XMLs_list = list()
	geneTree_treeLogs_XMLs = list()
	for (i in 1:length(geneTreeNames))
		{
		id = paste0("TreeHeight_of_", geneTreeNames[i])
		geneTreeName_idref = paste0("@", geneTreeNames[i])
		tmpXML = xmlNode(name="log", attrs=list(id=id, spec="beast.evolution.tree.TreeHeightLogger", tree=geneTreeName_idref))
		geneTree_treeHeights_XMLs_list = c(geneTree_treeHeights_XMLs_list, list(tmpXML))
		
		# Tree log
		id = paste0("TreeWithMetaDataLogger_for_", geneTreeNames[i])
		child = xmlNode(name="log", attrs=list(id=id, spec="beast.evolution.tree.TreeWithMetaDataLogger", tree=geneTreeName_idref))
		id = paste0("treeLog_for_", geneTreeNames[i])
		fn = paste0("treeLog_for_", geneTreeNames[i], ".trees.txt")
		tmpXML = xmlNode(name="logger", attrs=list(id=id, fileName=fn, logEvery=logEvery, mode="tree"), .children=list(child))
		geneTree_treeLogs_XMLs = c(geneTree_treeLogs_XMLs, list(tmpXML))
		}
	
	
	
	geneTree_treeHeights_XMLs_list = c(list(bl(), xmlCommentNode(" Loggers for genetree heights ")), geneTree_treeHeights_XMLs_list)
	geneTree_treeLogs_XMLs = c(list(bl(), xmlCommentNode(" Loggers for genetrees ")), geneTree_treeLogs_XMLs)
	
	
	if (is.null(xml))
		{
		res = NULL
		res$init = starBeast2_initializer_XML
		res$state = cl(res$state, speciesTree_statenode_XML, birthRate_XML, deathRate_XML)
		res$operators = cl(res$operators, operators_XMLs)
		res$priors = cl(res$priors, tree_prior_XML, birthRate_XML_priorLog$priors, deathRate_XML_priorLog$priors)
		res$tracelog = cl(res$tracelog, birthRate_XML_priorLog$tracelog, deathRate_XML_priorLog$tracelog, geneTree_treeHeights_XMLs_list)
		res$screenlog = cl(res$screenlog, birthRate_XML_priorLog$tracelog, deathRate_XML_priorLog$tracelog)
		res$treelog_XMLs = TreeWithMetaDataLogger_XMLs
		res$subslog_XMLs = Subs_TreeWithMetaDataLogger_XMLs
		res$speciesTreeLogger_XMLs = speciesTreeLogger_XMLs
		res$geneTree_treeLogs_XMLs = geneTree_treeLogs_XMLs
		tracelog
		return(res)
		} else {
		xml$init = cl(xml$init, starBeast2_initializer_XML)
		xml$state = cl(xml$state, speciesTree_statenode_XML, birthRate_XML, deathRate_XML)
		xml$operators = cl(xml$operators, operators_XMLs)
		xml$priors = cl(xml$priors, tree_prior_XML, birthRate_XML_priorLog$priors, deathRate_XML_priorLog$priors)
		xml$tracelog = cl(xml$tracelog, birthRate_XML_priorLog$tracelog, deathRate_XML_priorLog$tracelog, geneTree_treeHeights_XMLs_list)
		xml$screenlog = cl(xml$screenlog, birthRate_XML_priorLog$tracelog, deathRate_XML_priorLog$tracelog)
		xml$treelog = cl(xml$treelog, TreeWithMetaDataLogger_XMLs)
		xml$subslog = cl(xml$subslog, Subs_TreeWithMetaDataLogger_XMLs)
		xml$speciesTreeLogger_XMLs = speciesTreeLogger_XMLs
		xml$geneTree_treeLogs_XMLs = geneTree_treeLogs_XMLs
		return(xml)
		}
	
	return(stop("make_BD_model_for_starBeast2() ERROR: Shouldn't get here!"))
	} # END make_BD_model_for_starBeast2 <- function(treemodel_df, tree_name="shared_tree", xml=NULL)



#######################################################
# ONLY FOR NON-STARBEAST2 ANALYSES
#######################################################

make_BDSKY_model <- function(treemodel_df, tree_name="shared_tree", clockModel_name="shared_clock", clock_type="ucld", alignment_name_w_taxa, taxonset_w_taxa="list_of_OTUs", tipDates_id="tipDates", OTUs=NULL, xml=NULL, printall="short", relscale_fraction=0.05, XML_mod_for_cont_chars=FALSE, trstr=NULL, xlsfn=NULL)
	{
	defaults='
	treemodel_df = readWorksheetFromFile(file=xlsfn, sheet="treemodel", startRow=15, startCol=1, header=TRUE)
	tree_name = "shared_tree"
	alignment_name_w_taxa="seqs_morph_morph2_unordered"
	tipDates_id="tipDates"
	OTUs=NULL
	xml = NULL
	printall="short"
	relscale_fraction=0.05
	trstr=NULL
	'
	###############################################
	# starting Tree: Can be user-specified, or randomly determined
	###############################################
	# Setup / error fixes
	if (is.null(treemodel_df$use))
		{
		treemodel_df$use = rep("yes", times=nrow(treemodel_df))
		}
	treemodel_df$use[isblank_TF(treemodel_df$use) == TRUE] = "yes"
	treemodel_df$use = tolower(treemodel_df$use)
	treemodel_df = treemodel_df[treemodel_df$use == "yes"]
	
	# file_type: Newick or NEXUS
	if (is.null(treemodel_df$file_type))
		{
		treemodel_df$file_type = rep("newick", times=nrow(treemodel_df))
		}
	treemodel_df$file_type[isblank_TF(treemodel_df$file_type) == TRUE] = "newick"
	treemodel_df$file_type = tolower(treemodel_df$file_type)
	file_type = treemodel_df$file_type[1]
	
	if ((file_type %in% c("newick", "nexus")) == FALSE)
		{
		txt = paste0("STOP ERROR in make_starting_tree_for_starBeast2(): The starting tree type must be 'newick' or 'nexus', in treemodel_df$file_type.  Blank means 'newick'. However, you have: '", file_type, "'. Fix this in the Excel settings file, and re-run.")
		
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		}
	
	startingTree_option = treemodel_df$starting_tree[1]
	treeModel_option = treemodel_df$treeModel[1]
	
	if (length(startingTree_option) != 1)
		{
		errortxt = "\n\nSTOP ERROR in make_BDSKY_model(): there must be one and only entry in 'starting_tree' in worksheet 'treemodel'\n\n"
		cat(errortxt)
		stop(errortxt)
		} # END if (length(startingTree_option) != 1)

	if (length(treeModel_option) != 1)
		{
		errortxt = "\n\nSTOP ERROR in make_BDSKY_model(): there must be one and only entry in 'treeModel' in worksheet 'treemodel'\n\n"
		cat(errortxt)
		stop(errortxt)
		} # END if (length(treeModel_option) != 1)
	
	
	# Error check
	treeModel_options = c("BD", "BDSKY", "SABDSKY")
	if ( (treeModel_option %in% treeModel_options) == FALSE)
		{
		error_msg = "STOP ERROR in make_BDSKY_model(): treeModel_option, from treemodel_df$treeModel, must be one of the options allowed for a non-starBeast2 analysis."
		
		cat("\n\n")
		cat("error_msg")
		
		cat("\n\nAllowed options for treemodel_df$treeModel:\n\n")
		print(treeModel_options)
		cat("\n\n")
		
		stop(error_msg)
		}

	# Note: If you have a Newick file, you could use the same XML scheme as 
	# used for when there is a continuous trait -- however, this only seems
	# to work for BDSKY. It looks like you cannot combine
	# SABD and continuous traits.
	if ((XML_mod_for_cont_chars == TRUE) && (treeModel_option == "SABDSKY"))
		{
		txt = "STOP ERROR in make_BDSKY_model(): You are attempting to include continuous \ncharacters in an SABDSKY analysis. I cannot figure out how to make this work, it may\nrequire manual edits of the XML, or perhaps a change to the BEAST2 source code. -- Nick Matzke, 2015-04-02"
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if ((XML_mod_for_cont_chars == TRUE) && (treeModel_option == "SABDSKY"))
	


	# Tree XML method for continuous characters, or not?
	# Putting tree into XML in Remco's (RRB) way (necessary for continuous traits)
	if (XML_mod_for_cont_chars == TRUE)
		{
		# Putting tree into XML in Remco's (RRB) way (necessary for continuous traits)
		treetag = "init"
		# Tree XML method for continuous characters
		tree_name_init = paste0(tree_name, "_init")
		tree_name_idref = paste0("@", tree_name)
		tree_name_idref_for_init = paste0("@", tree_name)
		} else {
		# Putting tree into XML in normal way
		treetag = "stateNode"
		# Tree XML method, normal
		tree_name_init = paste0(tree_name, "")
		tree_name_idref = paste0("@", tree_name)
		tree_name_idref_for_init = NULL
		} # END if (XML_mod_for_cont_chars == TRUE)
	
	
	# idref for taxa list (and taxonset list for SABDSKY)
	alignment_name_w_taxaref = paste0("@", alignment_name_w_taxa)
	taxonset_w_taxaref = paste0("@", taxonset_w_taxa)
	
	
	

	if ((startingTree_option == "construct_starting_tree") && (is.null(xlsfn)) )
		{
		txt = "ERROR in make_BDSKY_model(): the 'startingTree_option' was 'construct_starting_tree', but xlsfn=NULL. Supply xlsfn (Excel filename (fn)) in order for make_BDSKY_model to use the 'construct_starting_tree' option. (Note: The 'startingTree_option' setting is derived from the column 'starting_tree' in treemodel_df, from the Excel worksheet 'treemodel' in the BEASTmasteR settings spreadsheet. The options are 'random', 'upgma', a newick filename, or 'construct_starting_tree'."
		
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)  
		}

	# We have an xlsfn, so...Construct the starting tree
	if (startingTree_option == "construct_starting_tree")
		{
		# Set the output filename
		outfn = "constructed_starting_tree.newick"
		
		# Change the option, so that this filename is now the starting tree
		startingTree_option = outfn
		
		# Construct the starting tree, write to outfn
		construct_starting_tree(xlsfn=xlsfn, min_brlen=0.0001, outfn=outfn)
		} # END if (startingTree_option == "construct_starting_tree")
	
	
	if ((startingTree_option == "random") || (startingTree_option == "upgma"))
		{
		# Generate a starting tree
		if (startingTree_option == "random")
			{
			# UPGMA starting tree (commented out)
			tipDates_idref = paste0("@", tipDates_id)
			UPGMA_startingTree_XML = xmlNode(name=treetag, attrs=list(id=tree_name_init, initial=tree_name_idref_for_init, spec="beast.util.ClusterZBSATree", nodetype="beast.evolution.tree.ZeroBranchSANode", clusterType="neighborjoining2", taxa=alignment_name_w_taxaref, trait=tipDates_idref) )		
			XMLstring = saveXML(doc=UPGMA_startingTree_XML, file=NULL, prefix=NULL)

			txt0_XML = xmlCommentNode(" RRB (Remco) suggested: add initialiser ")
			UPGMA_startingTree_XMLs = list(bl(), txt0_XML, xmlCommentNode(" Starting tree is constructed using UPGMA method (actually neighborjoining2) "), xmlCommentNode(XMLstring))


			# Parameter for the random starting tree population is 1.0 (?)
			popsize_id = paste0("randomPopSize_", tree_name)
			popsize_XML = xmlNode(name="parameter", 1.0, attrs=list(id=popsize_id, name="popSize") )
			
			# Population Model
			popModel_id = paste0("ConstantPopulation_", tree_name)
			popModel_XML = xmlNode(name="populationModel", attrs=list(id=popModel_id, spec="ConstantPopulation"), .children=list(popsize_XML) )
			
			# Starting Tree
			tipDates_idref = paste0("@", tipDates_id)
			
			if (treeModel_option == "BDSKY")
				{
				startingTree_XML = xmlNode(name=treetag, attrs=list(id=tree_name_init, initial=tree_name_idref_for_init, estimate="true", spec="beast.evolution.tree.RandomTree", taxa=alignment_name_w_taxaref, trait=tipDates_idref), .children=list(popModel_XML) )
				} # END if (treeModel_options == "BDSKY")
			if (treeModel_option == "SABDSKY")
				{
				startingTree_XML = xmlNode(name=treetag, attrs=list(id=tree_name_init, initial=tree_name_idref_for_init, estimate="true", spec="beast.evolution.tree.ZeroBranchSARandomTree", taxonset=taxonset_w_taxaref, trait=tipDates_idref), .children=list(popModel_XML) )
				} # END if (treeModel_options == "SABDSKY")			
				
			txt0_XML = xmlCommentNode(" RRB (Remco) suggested: add initialiser ")
			startingTree_XMLs = c(UPGMA_startingTree_XMLs, list(bl(), txt0_XML, xmlCommentNode(" Starting tree is constructed using random coalescence "), startingTree_XML))
			} # END if (startingTree_option == "random")

		if (startingTree_option == "upgma")
			{
			# Random Starting Tree (commented out)
			tipDates_idref = paste0("@", tipDates_id)
			
			if (treeModel_option == "BDSKY")
				{
				random_startingTree_XML = xmlNode(name=treetag, attrs=list(id=tree_name_init, initial=tree_name_idref_for_init, estimate="true", spec="beast.evolution.tree.RandomTree", taxa=alignment_name_w_taxaref, trait=tipDates_idref), .children=list(popModel_XML) )
				} # END if (treeModel_options == "BDSKY")
			if (treeModel_option == "SABDSKY")
				{
				random_startingTree_XML = xmlNode(name=treetag, attrs=list(id=tree_name_init, initial=tree_name_idref_for_init, estimate="true", spec="beast.evolution.tree.ZeroBranchSARandomTree", taxonset=taxonset_w_taxaref, trait=tipDates_idref), .children=list(popModel_XML) )
				} # END if (treeModel_options == "SABDSKY")			

			XMLstring = saveXML(doc=random_startingTree_XML, file=NULL, prefix=NULL)

			txt0_XML = xmlCommentNode(" RRB (Remco) suggested: add initialiser ")
			random_startingTree_XMLs = list(bl(), txt0_XML, xmlCommentNode(" Starting tree is constructed using random coalescence "), xmlCommentNode(XMLstring))


			# BEAST v2.2.0 Documentation: beast.util.ClusterZBSATree
			# Create initial beast.tree by hierarchical clustering, either through one of 
			# the classic link methods or by neighbor joining. The following link methods 
			# are supported: 
			# o single link, 
			# o complete link, 
			# o UPGMA=average link, 
			# o mean link, 
			# o centroid, 
			# o Ward and 
			# o adjusted complete link 
			# o neighborjoining 
			# o neighborjoining2 - corrects tree for tip data, unlike plain neighborjoining

			# Starting Tree
			tipDates_idref = paste0("@", tipDates_id)
			startingTree_XML = xmlNode(name=treetag, attrs=list(id=tree_name_init, initial=tree_name_idref_for_init, spec="beast.util.ClusterZBSATree", nodetype="beast.evolution.tree.ZeroBranchSANode", clusterType="neighborjoining2", taxa=alignment_name_w_taxaref, trait=tipDates_idref) )	

			txt0_XML = xmlCommentNode(" RRB (Remco) suggested: add initialiser ")
			startingTree_XMLs = c(random_startingTree_XMLs, list(bl(), txt0_XML, xmlCommentNode(" Starting tree is constructed using UPGMA method (actually neighborjoining2) "), startingTree_XML))
			} # END if (startingTree_option == "upgma")
		
		
		} else {
		# Use a user-specified starting tree
		if (isblank_TF(treemodel_df$tree_dir[1]) == TRUE)
			{
			trfn = startingTree_option
			} else {
			# Directory not blank, paste to filename
			trfn = slashslash(paste0(addslash(treemodel_df$tree_dir[1]), startingTree_option))
			} # END if (isblank_TF(treemodel_df$tree_dir[1]))
		
		if (printall != "none")
			{
			cat("\n\nmake_BDSKY_model() is running 'ape::read.tree()' on:\n'", getwd(), "\n", trfn, "'...\n\n", sep="")
			} # END if (printall != "none")
		
		
		#######################################################
		# Take the tree from a prior string, if specified
		# If not, take from a filename
		#######################################################
		if (is.null(trstr) == FALSE)
			{
			# Read the Newick string
			tr = read.tree(file="", text=trstr)
			} else {
			# Read the Newick or NEXUS file
			# (trfn comes from startingTree_option)
			# (which comes from startingTree_option = treemodel_df$starting_tree[1])
			#tr = read.tree(file=trfn)
			if (file_type == "newick")
				{
				tr = read.tree(file=trfn)
				} 
			if (file_type == "nexus")
				{
				tr = read.nexus(file=trfn)
				}
			} # END if (is.null(trstr) == FALSE)
		
		# Error check: Make sure there is just one (1) tree!
		if (class(tr) != "phylo")
			{
			txt = paste0("STOP ERROR in make_BDSKY_model(): there should be only 1 tree in '", trfn, "', but BEASTmasteR detects a class(tr)='", class(tr), "' object of length=", length(tr), ".")
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			cat("Printing 'tr':\n")
			print(tr)
			cat("Printing 'tr' string:\n")
			print(write_tree_noSci(tr, file=""))
			cat("\n\n")
			stop(txt)
			} # END if (length(tr) != 1)
		
		# Subset to remove missing OTUs
		if (is.null(OTUs) != TRUE)
			{
			# Remove tree tips not in OTUs 
			# (avoids manual subsetting of tree)
			tips_in_OTUs_TF = tr$tip.label %in% OTUs
			
			if (sum(tips_in_OTUs_TF == FALSE) > 0)
				{
				tips_to_drop = tr$tip.label[tips_in_OTUs_TF == FALSE]
				txt = paste0("WARNING in make_BDSKY_model(): Dropping ", length(tips_to_drop), " tips, that were not in 'OTUs', from the starting tree file\n", getwd(), "\n", trfn, "")
				cat("\n\n")
				cat(txt)
				cat("\n\n")
				cat("List of the dropped tips:\n")
				cat(tips_to_drop, sep="\n")
				cat("\n\n")
				warning(txt)
			
				# Perform the dropping
				tr2 = drop.tip(phy=tr, tip=tips_to_drop)
				tr = tr2
				} # END if (sum(tips_in_OTUs_TF == FALSE) > 0)
			
			# Check that no OTUs are missing from the tree
			OTUs_in_tr_TF = OTUs %in% tr$tip.label
			if (sum(OTUs_in_tr_TF == FALSE) > 0)
				{
				txt = paste0("STOP ERROR in make_BDSKY_model(): some specified 'OTUs' are nowhere in the Newick tree '", trfn, "'.")
				cat("\n\n")
				cat(txt)
				cat("\n\n")
				cat("List of the OTUs missing from Newick tree:\n")
				cat(OTUs[OTUs_in_tr_TF == FALSE], sep="\n")
				cat("\n\n")
				stop(txt)
				} # END if (sum(OTUs_in_tr_TF) > 0)
			} # END if (is.null(OTUs) != TRUE)
		
		# Write the tree string, for inclusion in the Newick file
		trstr = write_tree_noSci(tr, file="")
		
		if (treeModel_option == "BDSKY")
			{
			startingTree_XML = xmlNode(name=treetag, attrs=list(id=tree_name_init, initial=tree_name_idref_for_init, IsLabelledNewick="true", adjustTipHeights="false", estimate="true", spec="beast.util.TreeParser", taxa=alignment_name_w_taxaref, threshold="0.001", newick=trstr) )	
			} # END if (treeModel_option == "BDSKY")
			
		if (treeModel_option == "SABDSKY")
			{
			startingTree_XML = xmlNode(name=treetag, attrs=list(id=tree_name_init, initial=tree_name_idref_for_init, IsLabelledNewick="true", singlechild="true", estimate="true", threshold="0.001", spec="beast.util.ZeroBranchSATreeParser", taxa=alignment_name_w_taxaref, newick=trstr) )	
			} # END if (treeModel_option == "SABDSKY")


		#######################################################
		# Write code for random or UPGMA trees (for ease of 
		# use by user, it will be commented out)
		#######################################################
		# XML for random starting tree (commented out)
		# Parameter for the random starting tree population is 1.0 (?)
		popsize_id = paste0("randomPopSize_", tree_name)
		popsize_XML = xmlNode(name="parameter", 1.0, attrs=list(id=popsize_id, name="popSize") )

		# Population Model
		popModel_id = paste0("ConstantPopulation_", tree_name)
		popModel_XML = xmlNode(name="populationModel", attrs=list(id=popModel_id, spec="ConstantPopulation"), .children=list(popsize_XML) )

		# Starting Tree
		tipDates_idref = paste0("@", tipDates_id)
		random_startingTree_XML = xmlNode(name=treetag, attrs=list(id=tree_name_init, initial=tree_name_idref_for_init, estimate="true", spec="beast.evolution.tree.RandomTree", taxa=alignment_name_w_taxaref, trait=tipDates_idref), .children=list(popModel_XML) )

		XMLstring = saveXML(doc=random_startingTree_XML, file=NULL, prefix=NULL)

		txt0_XML = xmlCommentNode(" RRB (Remco) suggested: add initialiser ")
		random_startingTree_XMLs = list(bl(), txt0_XML, xmlCommentNode(" Starting tree is constructed using random coalescence "), xmlCommentNode(XMLstring))

		# XML for UPGMA starting tree (commented out)
		# Starting Tree
		# BEAST v2.2.0 Documentation: beast.util.ClusterZBSATree
		# Create initial beast.tree by hierarchical clustering, either through one of 
		# the classic link methods or by neighbor joining. The following link methods 
		# are supported: 
		# o single link, 
		# o complete link, 
		# o UPGMA=average link, 
		# o mean link, 
		# o centroid, 
		# o Ward and 
		# o adjusted complete link 
		# o neighborjoining 
		# o neighborjoining2 - corrects tree for tip data, unlike plain neighborjoining
		
		tipDates_idref = paste0("@", tipDates_id)
		UPGMA_startingTree_XML = xmlNode(name=treetag, attrs=list(id=tree_name_init, initial=tree_name_idref_for_init, spec="beast.util.ClusterZBSATree", nodetype="beast.evolution.tree.ZeroBranchSANode", clusterType="neighborjoining2", taxa=alignment_name_w_taxaref, taxonset="list_of_OTUs", trait=tipDates_idref) )	
		XMLstring = saveXML(doc=UPGMA_startingTree_XML, file=NULL, prefix=NULL)

		txt0_XML = xmlCommentNode(" RRB (Remco) suggested: add initialiser ")
		UPGMA_startingTree_XMLs = list(bl(), txt0_XML, xmlCommentNode(" Starting tree is constructed using UPGMA method (actually neighborjoining2) "), xmlCommentNode(XMLstring))
		#######################################################
		# END writing code for random or UPGMA trees
		#######################################################

		txt0_XML = xmlCommentNode(" RRB (Remco) suggested: add initialiser ")
		txt1_XML = xmlCommentNode(" Starting tree is user-specified Newick file from: ")
		txt2_XML = xmlCommentNode(paste0(" ", trfn, " "))
		startingTree_XMLs = c(random_startingTree_XMLs, UPGMA_startingTree_XMLs, list(bl(), txt0_XML, txt1_XML, txt2_XML, startingTree_XML))
		startingTree_XMLs
		} # END if ((startingTree_option == "random") || (startingTree_option == "upgma"))
	
	#######################################################
	# Copy the starting tree into the states (helps make 
	# continuous characters not crash; see Remco's (RRB's) advice)
	#######################################################
	
	if (XML_mod_for_cont_chars == TRUE)
		{
		# Remco says: RRB: create a separate tree
		# <tree name='stateNode' id="shared_tree" taxonset="@list_of_OTUs" />
		tree_RRM_comment_XML = xmlCommentNode(" Remco says: RRB: create a separate tree ")
		tree_RRM_statenode_XML = xmlNode(name="tree", attrs=list(id=tree_name, name="stateNode", taxonset="@list_of_OTUs"))
		tree_RRM_statenode_XMLs = list(bl(), tree_RRM_comment_XML, tree_RRM_statenode_XML)
		} else {
		# Normal, for situations with non-continuous traits
		# Put the starting tree straight in as a stateNode in <state>:
		tree_RRM_statenode_XMLs = startingTree_XMLs
		startingTree_XMLs = NULL
		}# END if (XML_mod_for_cont_chars == TRUE)
	
	###############################################
	# originTime: the time_bp at which the sampling process starts
	# (a nuisance parameter)
	###############################################
	originTime = treemodel_df$originTime[isblank_TF(treemodel_df$originTime)==FALSE]
	originTime_id = "originTime"
	originTime_idref = "@originTime"
	originTime_XML = xmlNode(name="parameter", originTime, attrs=list(id=originTime_id, lower="0.0", upper="4600", name="origin") )
		
	###############################################
	# birthRates, from oldest to youngest
	###############################################
	birthRates = treemodel_df$birthRate_starting_vals[isblank_TF(treemodel_df$birthRate_times_tops)==FALSE]
	birthRateChangeTimes_tops = treemodel_df$birthRate_times_tops[isblank_TF(treemodel_df$birthRate_times_tops)==FALSE]
	
	# Check lengths
	if (length(birthRates) != length(birthRateChangeTimes_tops))
		{
		txt = paste0("\n\nERROR in make_BDSKY_model():\nLength of treemodel_df$birthRate_starting_vals = ", length(birthRates), "\nLength of treemodel_df$birthRate_times_tops = ", length(birthRateChangeTimes_tops), "\nThese must be equal in the Excel spreadsheet.\n\n")
		cat(txt)
		stop(txt)
		}
	
	# Ages in youngest to oldest order
	order_indices = order(birthRateChangeTimes_tops)
	
	# We must print to XML in oldest to youngest order
	ages_txt = paste(birthRateChangeTimes_tops[rev(order_indices)], sep=" ", collapse=" ")
	rates_txt = paste(birthRates[rev(order_indices)], sep=" ", collapse=" ")
	rates_txt
	
	# Make the birthRate parameters
		
	# Min and max on birth rate - defaults
	if (isblank_TF(treemodel_df$birthRate_min[1]))
		{
		birthRate_min = NULL
		} else {
		birthRate_min = treemodel_df$birthRate_min[1]
		}
	if (isblank_TF(treemodel_df$birthRate_max[1]))
		{
		birthRate_max = NULL
		} else {
		birthRate_max = treemodel_df$birthRate_max[1]
		}
	
	birthRate_params_id = "birthRate"
	birthRate_params_idref = "@birthRate"
	birthRate_params_XML = xmlNode(name="parameter", attrs=list(id=birthRate_params_id, name="birthRate", lower=birthRate_min, upper=birthRate_max, value=rates_txt))
	birthRateChangeTimes_tops_id = "birthRateChangeTimes_tops"
	birthRateChangeTimes_tops_idref = paste0("@", birthRateChangeTimes_tops_id)
	birthRateChangeTimes_tops_XML = xmlNode(name="parameter", attrs=list(id=birthRateChangeTimes_tops_id, name="birthRateChangeTimes", value=ages_txt))
	
	#######################################################
	# Extract the priors on one or more birthRates
	#######################################################
	#treemodel_df = readWorksheetFromFile(file=xlsfn, sheet="treemodel", startRow=15, startCol=1, header=TRUE)
	
	birthRate_defn_XMLs = NULL
	birthRate_prior_XMLs = NULL
	birthRate_log_XMLs = NULL
	
	dflines = treemodel_df[isblank_TF(treemodel_df$birthRate_times_tops)==FALSE, ]
	num_birthRates = nrow(dflines)
	num_birthRates
	for (i in 1:num_birthRates)
		{
		dfline = dflines[i,]

		# Define the parameter names
		if (num_birthRates == 1)
			{
			param_name = "birthRate"
			} else {
			param_name = paste0("birthRate_bin", i-1)
			tmp_expression = paste0("birthRate[", i-1, "] * 1")
			
			child_XML = xmlNode(name="x", attrs=list(idref="birthRate") )
			
			birthRate_XML = xmlNode(name="parameter", attrs=list(id=param_name, spec="beast.util.Script", expression=tmp_expression), .children=list(child_XML))
			xmls = list(bl(), birthRate_XML)
			birthRate_defn_XMLs = c(birthRate_defn_XMLs, xmls)
			}
	
		# If distribution is fixed, make a very tight uniform prior (0.1%)
		if (dfline$birthRate_function == "fixed")
			{
			pass = 1
			# 		if ( (dfline$birthRate_function == "fixed") || (dfline$birthRate_starting_vals > 0))

# 			start_val = dfline$birthRate_starting_vals[i]
# 			param1_val = start_val - 0.001*start_val
# 			param2_val = start_val + 0.001*start_val
# 			tmpXML = make_generic_XML_prior(dfline, colname_prefix="birthRate", param_name=param_name, distrib="uniform", param1=param1_val, param2=param2_val, tmp_offset=0, meanInRealSpace="yes")
# 			tmpXML
# 			birthRate_prior_XMLs = c(birthRate_prior_XMLs, tmpXML$priors)
			} # END if ((dfline$ == "estimated") || (dfline$ == "estimated_base"))
	
	
	
		# Reasons to skip to next iteration of the loop
		if ( (dfline$birthRate_starting_vals==0 ) || 
			 (dfline$birthRate_function=="fixed") )
			{
			next()
			} # END if ((dfline$birthRate_starting_vals==0)||(dfline$birthRate_function=="fixed"))
		
		
		# Write priors, and log the prior probs, for estimated birthRates
		if ((dfline$birthRate_function == "estimated") || (dfline$birthRate_function == "estimated_base"))
			{
			print(dfline)
			tmpXML = make_generic_XML_prior(dfline, colname_prefix="birthRate", param_name=param_name)
			tmpXML
			birthRate_prior_XMLs = c(birthRate_prior_XMLs, tmpXML$priors)
			birthRate_log_XMLs = c(birthRate_log_XMLs, tmpXML$tracelog)
			} # END if ((dfline$ == "estimated") || (dfline$ == "estimated_base"))

		# Write priors, and log the prior probs, for scaling factors on birthRates
		if ( grepl(pattern="scaled", x=dfline$birthRate_function) == TRUE )
			{
			# Find which one is the base rate
			TF = dflines$birthRate_function[1:num_birthRates] == "estimated_base"
			if (sum(TF, na.rm=TRUE) != 1)
				{
				errortxt = "\n\nSTOP ERROR in make_BDSKY_model(): when scaling, you need to have have one and only one 'estimated_base' rate in the 'birthRate_function' column\n\n"
				cat(errortxt)
				stop(errortxt)
				} # END if (sum(TF, na.rm=TRUE) != 1)
			
			base_rownum = (1:num_birthRates)[TF]
			scaled_birthRate_equation = paste0("birthRate_bin", i-1, "")
			base_birthRate_equation = paste0("birthRate_bin", base_rownum-1, "")
			base_param_name = paste0("birthRate_bin", base_rownum-1)
			
			calculation_result_name = paste0("scale_difference_in_", scaled_birthRate_equation, "_vs_", base_birthRate_equation)
			scaled_prior_name = paste0("prior_on_", calculation_result_name)
			
			tmptxt = paste0(" Set up a prior on the SCALE DIFFERENCE in ", scaled_birthRate_equation, " vs. ", scaled_birthRate_equation, " ")
			txt_XML = xmlCommentNode(tmptxt)

			# scaled_min or scaled_max
			if ((grepl(pattern="min", x=dfline$birthRate_function) == TRUE) || (grepl(pattern="max", x=dfline$birthRate_function) == TRUE))
				{
				# Do the equation
				equation_txt = paste0(scaled_birthRate_equation, " / ", base_birthRate_equation)
				x_XML1 = xmlNode(name="x", attrs=list(idref=scaled_birthRate_equation) )
				x_XML2 = xmlNode(name="x", attrs=list(idref=base_birthRate_equation) )
				eqn_XML = xmlNode(name="x", attrs=list(id=calculation_result_name, spec="beast.util.Script", expression=equation_txt), .children=list(x_XML1, x_XML2) )
			
				# Do the probability density
				# scaled_min
				if (grepl(pattern="min", x=dfline$birthRate_function) == TRUE)
					{
					min_of_scaled = dfline$birthRate_scaling_relative_to_estimated_base
					min_of_uniform_id = paste0("min_of_UniformDistrib_on_", calculation_result_name)
					
					max_of_scaled = 10 * min_of_scaled
					max_of_uniform_id = paste0("max_of_UniformDistrib_on_", calculation_result_name)
					} # END if (grepl(pattern="min", x=dfline$birthRate_function) == TRUE)

				# scaled_max
				if (grepl(pattern="max", x=dfline$birthRate_function) == TRUE)
					{
					min_of_scaled = 0
					min_of_uniform_id = paste0("min_of_UniformDistrib_on_", calculation_result_name)
					
					max_of_scaled = dfline$birthRate_scaling_relative_to_estimated_base
					max_of_uniform_id = paste0("max_of_UniformDistrib_on_", calculation_result_name)
					} # END if (grepl(pattern="min", x=dfline$birthRate_function) == TRUE)			
				# Do the probability density
				uniform_id = paste0("UniformDistrib_on_", calculation_result_name)
				uniform_on_scaled_XML = xmlNode(name="Uniform", attrs=list(id=uniform_id, lower=min_of_scaled, upper=max_of_scaled, offset=0, spec="beast.math.distributions.Uniform", name="distr") )
			
				# Do the prior
				scaled_prior_XML = xmlNode(name="prior", attrs=list(id=scaled_prior_name, spec="beast.math.distributions.Prior", name="distribution"), .children=list(eqn_XML, uniform_on_scaled_XML))
				scaled_prior_XMLs = list(bl(), txt_XML, scaled_prior_XML)
			
				# Log the equation result, and the prior
				idref_of_scaled_prior_name = paste0(scaled_prior_name)
				tmplog_XMLtxt = xmlCommentNode(paste0(" Log of the prior probability density of '", scaled_prior_name, "' "))
				tmplog_XML = xmlNode(name="log", attrs=list(idref=idref_of_scaled_prior_name))
				tmplog_XMLs1 = list(bl(), tmplog_XMLtxt, tmplog_XML)

				idref_of_calculation_result_name = paste0(calculation_result_name)
				tmplog_XMLtxt = xmlCommentNode(paste0(" Log of the value of '", calculation_result_name, "' "))
				tmplog_XML = xmlNode(name="log", attrs=list(idref=idref_of_calculation_result_name))
				tmplog_XMLs2 = list(bl(), tmplog_XMLtxt, tmplog_XML)

				tmplog_XMLs = c(tmplog_XMLs1, tmplog_XMLs2)
			
				# Store
				birthRate_prior_XMLs = c(birthRate_prior_XMLs, scaled_prior_XMLs)
				birthRate_log_XMLs = c(birthRate_log_XMLs, tmpXML$tracelog, tmplog_XMLs)
				} else {
				# Normal distribution on scaled multiplier
			
				# Do the equation
				scaled_birthRate_equation = paste0("birthRate_bin", i-1, "")
				base_birthRate_equation = paste0("birthRate_bin", base_rownum-1, "")

				equation_txt = paste0(scaled_birthRate_equation, " / ", base_birthRate_equation)
				x_XML1 = xmlNode(name="x", attrs=list(idref=scaled_birthRate_equation) )
				x_XML2 = xmlNode(name="x", attrs=list(idref=base_birthRate_equation) )
				eqn_XML = xmlNode(name="x", attrs=list(id=calculation_result_name, spec="beast.util.Script", expression=equation_txt), .children=list(x_XML1, x_XML2) )
			
				# Do the probability density
				# Mean
				mean_of_scaled = dfline$birthRate_scaling_relative_to_estimated_base
				mean_of_normal_id = paste0("mean_of_NormalDistrib_on_", calculation_result_name)
				param1_XML = xmlNode(name="parameter", value=mean_of_scaled, attrs=list(id=mean_of_normal_id, name="mean") )
			
				# Standard deviation
				sd_of_scaled = relscale_fraction * mean_of_scaled
				sd_of_normal_id = paste0("sd_of_NormalDistrib_on_", calculation_result_name)
				param2_XML = xmlNode(name="parameter", value=sd_of_scaled, attrs=list(id=sd_of_normal_id, name="sigma") )
			
				# Do the probability density
				normal_id = paste0("NormalDistrib_on_", calculation_result_name)
				normal_on_scaled_XML = xmlNode(name="Normal", attrs=list(id=normal_id, spec="beast.math.distributions.Normal", name="distr"), .children=list(param1_XML, param2_XML) )
			
				# Do the prior
				scaled_prior_XML = xmlNode(name="prior", attrs=list(id=scaled_prior_name, spec="beast.math.distributions.Prior", name="distribution"), .children=list(eqn_XML, normal_on_scaled_XML))
				scaled_prior_XMLs = list(bl(), txt_XML, scaled_prior_XML)
			
				# Log the equation result, and the prior
				idref_of_scaled_prior_name = paste0(scaled_prior_name)
				tmplog_XMLtxt = xmlCommentNode(paste0(" Log of the prior probability density of '", scaled_prior_name, "' "))
				tmplog_XML = xmlNode(name="log", attrs=list(idref=idref_of_scaled_prior_name))
				tmplog_XMLs1 = list(bl(), tmplog_XMLtxt, tmplog_XML)

				idref_of_calculation_result_name = paste0(calculation_result_name)
				tmplog_XMLtxt = xmlCommentNode(paste0(" Log of the value of '", calculation_result_name, "' "))
				tmplog_XML = xmlNode(name="log", attrs=list(idref=idref_of_calculation_result_name))
				tmplog_XMLs2 = list(bl(), tmplog_XMLtxt, tmplog_XML)

				tmplog_XMLs = c(tmplog_XMLs1, tmplog_XMLs2)
			
				# Store
				birthRate_prior_XMLs = c(birthRate_prior_XMLs, scaled_prior_XMLs)
				birthRate_log_XMLs = c(birthRate_log_XMLs, tmpXML$tracelog, tmplog_XMLs)
				} # END if ( grepl(pattern="scaled", x=dfline$birthRate_function) == TRUE )
			} # END if ((grepl(pattern="min", x=dfline$birthRate_function) == TRUE) || (grepl(pattern="max", x=dfline$birthRate_function) == TRUE))		
		} # END for (i in 1:num_birthRates)
	
	birthRate_defn_XMLs
	birthRate_prior_XMLs
	birthRate_log_XMLs
	#######################################################
	# End absolute and relative priors on birthRates
	#######################################################
	
	
	
	###############################################
	# deathRates, from oldest to youngest
	###############################################
	deathRates = treemodel_df$deathRate_starting_vals[isblank_TF(treemodel_df$deathRate_times_tops)==FALSE]
	deathRateChangeTimes_tops = treemodel_df$deathRate_times_tops[isblank_TF(treemodel_df$deathRate_times_tops)==FALSE]
	
	# Check lengths
	if (length(deathRates) != length(deathRateChangeTimes_tops))
		{
		txt = paste0("\n\nERROR in make_BDSKY_model():\nLength of treemodel_df$deathRate_starting_vals = ", length(deathRates), "\nLength of treemodel_df$deathRate_times_tops = ", length(deathRateChangeTimes_tops), "\nThese must be equal in the Excel spreadsheet.\n\n")
		cat(txt)
		stop(txt)
		}
	
	# Ages in youngest to oldest order
	order_indices = order(deathRateChangeTimes_tops)
	
	# We must print to XML in oldest to youngest order
	ages_txt = paste(deathRateChangeTimes_tops[rev(order_indices)], sep=" ", collapse=" ")
	rates_txt = paste(deathRates[rev(order_indices)], sep=" ", collapse=" ")
	rates_txt
	
	# Make the deathRate parameters
	# Min and max on death rate - defaults
	if (isblank_TF(treemodel_df$deathRate_min[1]))
		{
		deathRate_min = NULL
		} else {
		deathRate_min = treemodel_df$deathRate_min[1]
		}
	if (isblank_TF(treemodel_df$deathRate_max[1]))
		{
		deathRate_max = NULL
		} else {
		deathRate_max = treemodel_df$deathRate_max[1]
		}

	deathRate_params_id = "deathRate"
	deathRate_params_idref = "@deathRate"
	deathRate_params_XML = xmlNode(name="parameter", attrs=list(id=deathRate_params_id, name="deathRate", lower=deathRate_min, upper=deathRate_max, value=rates_txt))
	deathRateChangeTimes_tops_id = "deathRateChangeTimes_tops"
	deathRateChangeTimes_tops_idref = paste0("@", deathRateChangeTimes_tops_id)
	deathRateChangeTimes_tops_XML = xmlNode(name="parameter", attrs=list(id=deathRateChangeTimes_tops_id, name="deathRateChangeTimes", value=ages_txt))



	#######################################################
	# Extract the priors on one or more deathRates
	#######################################################
	#treemodel_df = readWorksheetFromFile(file=xlsfn, sheet="treemodel", startRow=15, startCol=1, header=TRUE)
	
	deathRate_defn_XMLs = NULL
	deathRate_prior_XMLs = NULL
	deathRate_log_XMLs = NULL
	
	dflines = treemodel_df[isblank_TF(treemodel_df$deathRate_times_tops)==FALSE, ]
	num_deathRates = nrow(dflines)
	num_deathRates
	for (i in 1:num_deathRates)
		{
		dfline = dflines[i,]

		# Define the parameter names
		if (num_deathRates == 1)
			{
			param_name = "deathRate"
			} else {
			param_name = paste0("deathRate_bin", i-1)
			tmp_expression = paste0("deathRate[", i-1, "] * 1")

			child_XML = xmlNode(name="x", attrs=list(idref="deathRate") )
			
			deathRate_XML = xmlNode(name="parameter", attrs=list(id=param_name, spec="beast.util.Script", expression=tmp_expression), .children=list(child_XML))
			
			xmls = list(bl(), deathRate_XML)
			deathRate_defn_XMLs = c(deathRate_defn_XMLs, xmls)
			}
	
		# If distribution is fixed, make a very tight uniform prior (0.1%)
		if (dfline$deathRate_function == "fixed")
			{
			pass = 1

# 			start_val = dfline$deathRate_starting_vals[i]
# 			param1_val = start_val - 0.001*start_val
# 			param2_val = start_val + 0.001*start_val
# 			tmpXML = make_generic_XML_prior(dfline, colname_prefix="deathRate", param_name=param_name, distrib="uniform", param1=param1_val, param2=param2_val, tmp_offset=0, meanInRealSpace="yes")
# 			tmpXML
# 			deathRate_prior_XMLs = c(deathRate_prior_XMLs, tmpXML$priors)
			} # END if ((dfline$ == "estimated") || (dfline$ == "estimated_base"))
	
	
	
		# Reasons to skip to next iteration of the loop
		if ( (dfline$deathRate_starting_vals==0 ) || 
			 (dfline$deathRate_function=="fixed") )
			{
			next()
			} # END if ((dfline$deathRate_starting_vals==0)||(dfline$deathRate_function=="fixed"))
		
		
		# Write priors, and log the prior probs, for estimated deathRates
		if ((dfline$deathRate_function == "estimated") || (dfline$deathRate_function == "estimated_base"))
			{
			print(dfline)
			tmpXML = make_generic_XML_prior(dfline, colname_prefix="deathRate", param_name=param_name)
			tmpXML
			deathRate_prior_XMLs = c(deathRate_prior_XMLs, tmpXML$priors)
			deathRate_log_XMLs = c(deathRate_log_XMLs, tmpXML$tracelog)
			} # END if ((dfline$ == "estimated") || (dfline$ == "estimated_base"))

		# Write priors, and log the prior probs, for scaling factors on deathRates
		if ( grepl(pattern="scaled", x=dfline$deathRate_function) == TRUE )
			{
			# Find which one is the base rate
			TF = dflines$deathRate_function[1:num_deathRates] == "estimated_base"

			if (sum(TF, na.rm=TRUE) != 1)
				{
				errortxt = "\n\nSTOP ERROR in make_BDSKY_model(): when scaling, you need to have have one and only one 'estimated_base' rate in the 'deathRate_function' column\n\n"
				cat(errortxt)
				stop(errortxt)
				} # END if (sum(TF, na.rm=TRUE) != 1)

			base_rownum = (1:num_deathRates)[TF]
			scaled_deathRate_equation = paste0("deathRate_bin", i-1, "")
			base_deathRate_equation = paste0("deathRate_bin", base_rownum-1, "")
			base_param_name = paste0("deathRate_bin", base_rownum-1)
			
			calculation_result_name = paste0("scale_difference_in_", scaled_deathRate_equation, "_vs_", base_deathRate_equation)
			scaled_prior_name = paste0("prior_on_", calculation_result_name)
			
			tmptxt = paste0(" Set up a prior on the SCALE DIFFERENCE in ", scaled_deathRate_equation, " vs. ", base_deathRate_equation, " ")
			txt_XML = xmlCommentNode(tmptxt)
			

			# scaled_min or scaled_max
			if ((grepl(pattern="min", x=dfline$deathRate_function) == TRUE) || (grepl(pattern="max", x=dfline$deathRate_function) == TRUE))
				{
				# Do the equation
				equation_txt = paste0(scaled_deathRate_equation, " / ", base_deathRate_equation)
				x_XML1 = xmlNode(name="x", attrs=list(idref=scaled_deathRate_equation) )
				x_XML2 = xmlNode(name="x", attrs=list(idref=base_deathRate_equation) )
				eqn_XML = xmlNode(name="x", attrs=list(id=calculation_result_name, spec="beast.util.Script", expression=equation_txt), .children=list(x_XML1, x_XML2) )
			
				# Do the probability density
				# scaled_min
				if (grepl(pattern="min", x=dfline$deathRate_function) == TRUE)
					{
					min_of_scaled = dfline$deathRate_scaling_relative_to_estimated_base
					min_of_uniform_id = paste0("min_of_UniformDistrib_on_", calculation_result_name)
					
					max_of_scaled = 10 * min_of_scaled
					max_of_uniform_id = paste0("max_of_UniformDistrib_on_", calculation_result_name)
					} # END if (grepl(pattern="min", x=dfline$deathRate_function) == TRUE)

				# scaled_max
				if (grepl(pattern="max", x=dfline$deathRate_function) == TRUE)
					{
					min_of_scaled = 0
					min_of_uniform_id = paste0("min_of_UniformDistrib_on_", calculation_result_name)
					
					max_of_scaled = dfline$deathRate_scaling_relative_to_estimated_base
					max_of_uniform_id = paste0("max_of_UniformDistrib_on_", calculation_result_name)
					} # END if (grepl(pattern="min", x=dfline$deathRate_function) == TRUE)			
				# Do the probability density
				uniform_id = paste0("UniformDistrib_on_", calculation_result_name)
				uniform_on_scaled_XML = xmlNode(name="Uniform", attrs=list(id=uniform_id, lower=min_of_scaled, upper=max_of_scaled, offset=0, spec="beast.math.distributions.Uniform", name="distr") )
			
				# Do the prior
				scaled_prior_XML = xmlNode(name="prior", attrs=list(id=scaled_prior_name, spec="beast.math.distributions.Prior", name="distribution"), .children=list(eqn_XML, uniform_on_scaled_XML))
				scaled_prior_XMLs = list(bl(), txt_XML, scaled_prior_XML)
			
				# Log the equation result, and the prior
				idref_of_scaled_prior_name = paste0(scaled_prior_name)
				tmplog_XMLtxt = xmlCommentNode(paste0(" Log of the prior probability density of '", scaled_prior_name, "' "))
				tmplog_XML = xmlNode(name="log", attrs=list(idref=idref_of_scaled_prior_name))
				tmplog_XMLs1 = list(bl(), tmplog_XMLtxt, tmplog_XML)

				idref_of_calculation_result_name = paste0(calculation_result_name)
				tmplog_XMLtxt = xmlCommentNode(paste0(" Log of the value of '", calculation_result_name, "' "))
				tmplog_XML = xmlNode(name="log", attrs=list(idref=idref_of_calculation_result_name))
				tmplog_XMLs2 = list(bl(), tmplog_XMLtxt, tmplog_XML)

				tmplog_XMLs = c(tmplog_XMLs1, tmplog_XMLs2)
			
				# Store
				deathRate_prior_XMLs = c(deathRate_prior_XMLs, scaled_prior_XMLs)
				deathRate_log_XMLs = c(deathRate_log_XMLs, tmpXML$tracelog, tmplog_XMLs)
				} else {
				# Normal distribution on scaled multiplier

				# Do the equation
				scaled_deathRate_equation = paste0("deathRate_bin", i-1, "")
				base_deathRate_equation = paste0("deathRate_bin", base_rownum-1, "")
				equation_txt = paste0(scaled_deathRate_equation, " / ", base_deathRate_equation)
				x_XML1 = xmlNode(name="x", attrs=list(idref=scaled_deathRate_equation) )
				x_XML2 = xmlNode(name="x", attrs=list(idref=base_deathRate_equation) )
				eqn_XML = xmlNode(name="x", attrs=list(id=calculation_result_name, spec="beast.util.Script", expression=equation_txt), .children=list(x_XML1, x_XML2) )
			
				# Do the probability density
				# Mean
				mean_of_scaled = dfline$deathRate_scaling_relative_to_estimated_base
				mean_of_normal_id = paste0("mean_of_NormalDistrib_on_", calculation_result_name)
				param1_XML = xmlNode(name="parameter", value=mean_of_scaled, attrs=list(id=mean_of_normal_id, name="mean") )
			
				# Standard deviation
				sd_of_scaled = relscale_fraction * mean_of_scaled
				sd_of_normal_id = paste0("sd_of_NormalDistrib_on_", calculation_result_name)
				param2_XML = xmlNode(name="parameter", value=sd_of_scaled, attrs=list(id=sd_of_normal_id, name="sigma") )
			
				# Do the probability density
				normal_id = paste0("NormalDistrib_on_", calculation_result_name)
				normal_on_scaled_XML = xmlNode(name="Normal", attrs=list(id=normal_id, spec="beast.math.distributions.Normal", name="distr"), .children=list(param1_XML, param2_XML) )
			
				# Do the prior
				scaled_prior_XML = xmlNode(name="prior", attrs=list(id=scaled_prior_name, spec="beast.math.distributions.Prior", name="distribution"), .children=list(eqn_XML, normal_on_scaled_XML))
				scaled_prior_XMLs = list(bl(), txt_XML, scaled_prior_XML)
			
				# Log the equation result, and the prior
				idref_of_scaled_prior_name = paste0(scaled_prior_name)
				tmplog_XMLtxt = xmlCommentNode(paste0(" Log of the prior probability density of '", scaled_prior_name, "' "))
				tmplog_XML = xmlNode(name="log", attrs=list(idref=idref_of_scaled_prior_name))
				tmplog_XMLs1 = list(bl(), tmplog_XMLtxt, tmplog_XML)

				idref_of_calculation_result_name = paste0(calculation_result_name)
				tmplog_XMLtxt = xmlCommentNode(paste0(" Log of the value of '", calculation_result_name, "' "))
				tmplog_XML = xmlNode(name="log", attrs=list(idref=idref_of_calculation_result_name))
				tmplog_XMLs2 = list(bl(), tmplog_XMLtxt, tmplog_XML)

				tmplog_XMLs = c(tmplog_XMLs1, tmplog_XMLs2)
			
				# Store
				deathRate_prior_XMLs = c(deathRate_prior_XMLs, scaled_prior_XMLs)
				deathRate_log_XMLs = c(deathRate_log_XMLs, tmpXML$tracelog, tmplog_XMLs)
				} # END if ( grepl(pattern="scaled", x=dfline$deathRate_function) == TRUE )
			} # END if ((grepl(pattern="min", x=dfline$deathRate_function) == TRUE) || (grepl(pattern="max", x=dfline$deathRate_function) == TRUE))		
		} # END for (i in 1:num_deathRates)
	
	deathRate_defn_XMLs
	deathRate_prior_XMLs
	deathRate_log_XMLs
	#######################################################
	# End absolute and relative priors on deathRates
	#######################################################







	

	#######################################################
	# Check: the deathRate cannot exceed the birthRate
	#######################################################
	example='
	originTime = 20
	birthRates = c(0.1, 0.2, 0.3)
	deathRates = c(0.2)
	birthRateChangeTimes_tops = c(0, 5, 10)
	deadRateChangeTimes_tops = c(0)
	deathRateChangeTimes_tops = c(0)
	'
	
	times = seq(0, originTime, (originTime / 50))
	birthTime_tops_tmp = c(birthRateChangeTimes_tops, originTime)
	deathTime_tops_tmp = c(deathRateChangeTimes_tops, originTime)
	birth_death_rates_table = NULL
	for (ii in 1:(length(times)-1))
		{
		tmptime = times[ii]
		bTF1 = tmptime >= birthTime_tops_tmp[1:(length(birthTime_tops_tmp)-1)]
		bTF2 = tmptime < birthTime_tops_tmp[2:(length(birthTime_tops_tmp)-0)]
		bTF = (bTF1 + bTF2) == 2
		dTF1 = tmptime >= deathTime_tops_tmp[1:(length(deathTime_tops_tmp)-1)]
		dTF2 = tmptime < deathTime_tops_tmp[2:(length(deathTime_tops_tmp)-0)]
		dTF = (dTF1 + dTF2) == 2
		brate = birthRates[bTF]
		drate = deathRates[dTF]
		OK = (brate > drate)
		
		tmprow = c(tmptime, brate, drate, OK)
		
		birth_death_rates_table = rbind(birth_death_rates_table, tmprow)
		}
	birth_death_rates_table = as.data.frame(birth_death_rates_table, stringsAsFactors=FALSE, row.names=NULL)
	row.names(birth_death_rates_table) = NULL
	names(birth_death_rates_table) = c("time", "brate", "drate", "OK")
	TF = birth_death_rates_table$OK == 1
	birth_death_rates_table$OK = TF
	birth_death_rates_table
	
	if (any(birth_death_rates_table$OK == FALSE))
		{
		error_txt = paste0("\n\nSTOP ERROR in make_BDSKY_model(): In Beast2 tree models, it appears that the birthRate can never exceed the deathRate. This is probably due to some simplification in equations used by Beast2 in the code calculating the likelihood of a tree given the birthRates, deathRates, etc. If deathRate exceeds birthRate, Beast2 will crash with a line like: 'P(BDSKY_prior_shared_tree) = NaN'.\n\n")
		
		txt2 = paste0("Yes, this is somewhat silly, since it is perfectly possible in real life that the deathRate could exceed the birthRate for a period of time, and still produce an observed tree with some probability. In a purely molecular tree of purely living taxa, it may be that the ML estimate of birthRate must always be bigger than the deathRate, but in a tree with fossils this is not the case. So, we should lobby the Beast2 developers to change the equation (or perhaps I am misunderstanding something).\n\n")
		txt3 = paste0("Until then, make sure your starting birthRates always exceed your starting deathRates.\n\nThe table printed below shows approximately where you have this problem (where the 'OK' column equals 'FALSE').\n\n")
		
		cat(error_txt)
		cat(txt2)
		cat(txt3)
		
		print(birth_death_rates_table)
		
		stop(error_txt)
		} # END if (any(birth_death_rates_table$OK) == FALSE)
	
	
	
	
	###############################################
	# samplingRate dates of the tops of bins, sorted from youngest to oldest, but with zero at 
	# the end, as specified by Denise on the beast-users group
	# https://groups.google.com/d/msg/beast-users/UyAS0sYDkKk/PB0ktWnssn0J
	# (reverseTimeArrays is, I guess, about whether time itself is backwards or forwards)
	###############################################
	samplingRates = treemodel_df$samplingRate_starting_vals[isblank_TF(treemodel_df$samplingRate_times_tops)==FALSE]
	samplingRateTimes_tops = treemodel_df$samplingRate_times_tops[isblank_TF(treemodel_df$samplingRate_times_tops)==FALSE]
	
	# Check lengths
	if (length(samplingRates) != length(samplingRateTimes_tops))
		{
		txt = paste0("\n\nERROR in make_BDSKY_model():\nLength of treemodel_df$samplingRate_starting_vals = ", length(samplingRates), "\nLength of treemodel_df$samplingRate_times_tops = ", length(samplingRateTimes_tops), "\nThese must be equal in the Excel spreadsheet.\n\n")
		cat(txt)
		stop(txt)
		}
	
	# Add a zero, if there isn't one
	# and, make its rate 0
	if ((0 %in% samplingRateTimes_tops) == FALSE)
		{
		samplingRateTimes_tops = c(samplingRateTimes_tops, 0)
		samplingRates = c(samplingRates, 0)
		} # END if ((0 %in% samplingRateTimes_tops) == FALSE)
	
	# Ages in youngest to oldest order
	order_indices = order(samplingRateTimes_tops)
	
	# samplingRateTimes_tops in Denise's order:
	Denises_samplingRateTimes_tops_order = TRUE
	if (Denises_samplingRateTimes_tops_order == TRUE)
		{
		# Skip the first zero...
		everything_but_zero = samplingRateTimes_tops[order_indices][-1]
		# ...and put it at the end...
		samplingRateTimes_tops_to_use = c(everything_but_zero, samplingRateTimes_tops[order_indices][1])
		} else {
		# Original default; seemed to not work with most recent samplingRate being 0
		samplingRateTimes_tops_to_use = samplingRateTimes_tops[rev(order_indices)]		
		} # END if (Denises_samplingRateTimes_tops_order == TRUE)
		
	
	# We must print to XML in oldest to youngest order
	ages_txt = paste(samplingRateTimes_tops_to_use, sep=" ", collapse=" ")
	rates_txt = paste(samplingRates[rev(order_indices)], sep=" ", collapse=" ")
	rates_txt
	
	# Make the samplingRate parameters
	# Min and max on sampling rate - defaults
	if (isblank_TF(treemodel_df$samplingRate_min[1]))
		{
		samplingRate_min = NULL
		} else {
		samplingRate_min = treemodel_df$samplingRate_min[1]
		}
	if (isblank_TF(treemodel_df$samplingRate_max[1]))
		{
		samplingRate_max = NULL
		} else {
		samplingRate_max = treemodel_df$samplingRate_max[1]
		}
	
	samplingRate_params_id = "samplingRate"
	samplingRate_params_idref = "@samplingRate"
	samplingRate_params_XML = xmlNode(name="parameter", attrs=list(id=samplingRate_params_id, name="samplingRate", lower=samplingRate_min, upper=samplingRate_max, value=rates_txt))
	samplingRateTimes_tops_id = "samplingRateChangeTimes_tops"
	samplingRateTimes_tops_idref = paste0("@", samplingRateTimes_tops_id)
	samplingRateTimes_tops_XML = xmlNode(name="parameter", attrs=list(id=samplingRateTimes_tops_id, name="samplingRateChangeTimes", value=ages_txt))
	


	#######################################################
	# Extract the priors on one or more samplingRates
	#######################################################
	#treemodel_df = readWorksheetFromFile(file=xlsfn, sheet="treemodel", startRow=15, startCol=1, header=TRUE)
	
	samplingRate_defn_XMLs = NULL
	samplingRate_prior_XMLs = NULL
	samplingRate_log_XMLs = NULL
	
	dflines = treemodel_df[isblank_TF(treemodel_df$samplingRate_times_tops)==FALSE, ]
	num_samplingRates = nrow(dflines)
	num_samplingRates
	for (i in 1:num_samplingRates)
		{
		dfline = dflines[i,]

		# Define the parameter names
		if (num_samplingRates == 1)
			{
			param_name = "samplingRate"
			} else {
			# Weird Beast2 numbering scheme -- last entry, e.g. 0-based 3, 1-based 4,
			# is actually the rate at time 0
			if (i == num_samplingRates)
				{
				param_name = paste0("samplingRate_bin", 0)
				tmp_expression = paste0("samplingRate[", i-1, "] * 1")
				} else {
				param_name = paste0("samplingRate_bin", i)
				tmp_expression = paste0("samplingRate[", i-1, "] * 1")
				}
			
			child_XML = xmlNode(name="x", attrs=list(idref="samplingRate") )
						
			samplingRate_XML = xmlNode(name="parameter", attrs=list(id=param_name, spec="beast.util.Script", expression=tmp_expression), .children=list(child_XML))
			xmls = list(bl(), samplingRate_XML)
			
			
			
			samplingRate_defn_XMLs = c(samplingRate_defn_XMLs, xmls)
			}
	
		# If distribution is fixed, make a very tight uniform prior (0.1%)
		if (dfline$samplingRate_function == "fixed")
			{
			pass = 1

# 			start_val = dfline$samplingRate_starting_vals[i]
# 			param1_val = start_val - 0.001*start_val
# 			param2_val = start_val + 0.001*start_val
# 			tmpXML = make_generic_XML_prior(dfline, colname_prefix="samplingRate", param_name=param_name, distrib="uniform", param1=param1_val, param2=param2_val, tmp_offset=0, meanInRealSpace="yes")
# 			tmpXML
# 			samplingRate_prior_XMLs = c(samplingRate_prior_XMLs, tmpXML$priors)
			} # END if ((dfline$ == "estimated") || (dfline$ == "estimated_base"))
	
	
	
		# Reasons to skip to next iteration of the loop
		if ( (dfline$samplingRate_starting_vals==0 ) || 
			 (dfline$samplingRate_function=="fixed") )
			{
			next()
			} # END if ((dfline$samplingRate_starting_vals==0)||(dfline$samplingRate_function=="fixed"))
		
		
		# Write priors, and log the prior probs, for estimated samplingRates
		if ((dfline$samplingRate_function == "estimated") || (dfline$samplingRate_function == "estimated_base"))
			{
			print(dfline)
			tmpXML = make_generic_XML_prior(dfline, colname_prefix="samplingRate", param_name=param_name)
			tmpXML
			samplingRate_prior_XMLs = c(samplingRate_prior_XMLs, tmpXML$priors)
			samplingRate_log_XMLs = c(samplingRate_log_XMLs, tmpXML$tracelog)
			} # END if ((dfline$ == "estimated") || (dfline$ == "estimated_base"))

		# Write priors, and log the prior probs, for scaling factors on samplingRates
		if ( grepl(pattern="scaled", x=dfline$samplingRate_function) == TRUE )
			{
			# Find which one is the base rate
			TF = dflines$samplingRate_function[1:num_samplingRates] == "estimated_base"

			if (sum(TF, na.rm=TRUE) != 1)
				{
				errortxt = "\n\nSTOP ERROR in make_BDSKY_model(): when scaling, you need to have have one and only one 'estimated_base' rate in the 'samplingRate_function' column\n\n"
				cat(errortxt)
				stop(errortxt)
				} # END if (sum(TF, na.rm=TRUE) != 1)

			base_rownum = (1:num_samplingRates)[TF]
			scaled_samplingRate_equation = paste0("samplingRate_bin", i-1, "")
			base_samplingRate_equation = paste0("samplingRate_bin", base_rownum-1, "")
			base_param_name = paste0("samplingRate_bin", base_rownum-1)
			
			calculation_result_name = paste0("scale_difference_in_", scaled_samplingRate_equation, "_vs_", base_samplingRate_equation)
			scaled_prior_name = paste0("prior_on_", calculation_result_name)
			
			tmptxt = paste0(" Set up a prior on the SCALE DIFFERENCE in ", scaled_samplingRate_equation, " vs. ", base_samplingRate_equation, " ")
			txt_XML = xmlCommentNode(tmptxt)
			
			# scaled_min or scaled_max
			if ((grepl(pattern="min", x=dfline$samplingRate_function) == TRUE) || (grepl(pattern="max", x=dfline$samplingRate_function) == TRUE))
				{
				# Do the equation
				equation_txt = paste0(scaled_samplingRate_equation, " / ", base_samplingRate_equation)
				x_XML1 = xmlNode(name="x", attrs=list(idref=scaled_samplingRate_equation) )
				x_XML2 = xmlNode(name="x", attrs=list(idref=base_samplingRate_equation) )
				eqn_XML = xmlNode(name="x", attrs=list(id=calculation_result_name, spec="beast.util.Script", expression=equation_txt), .children=list(x_XML1, x_XML2) )
			
				# Do the probability density
				# scaled_min
				if (grepl(pattern="min", x=dfline$samplingRate_function) == TRUE)
					{
					min_of_scaled = dfline$samplingRate_scaling_relative_to_estimated_base
					min_of_uniform_id = paste0("min_of_UniformDistrib_on_", calculation_result_name)
					
					max_of_scaled = 10 * min_of_scaled
					max_of_uniform_id = paste0("max_of_UniformDistrib_on_", calculation_result_name)
					} # END if (grepl(pattern="min", x=dfline$samplingRate_function) == TRUE)

				# scaled_max
				if (grepl(pattern="max", x=dfline$samplingRate_function) == TRUE)
					{
					min_of_scaled = 0
					min_of_uniform_id = paste0("min_of_UniformDistrib_on_", calculation_result_name)
					
					max_of_scaled = dfline$samplingRate_scaling_relative_to_estimated_base
					max_of_uniform_id = paste0("max_of_UniformDistrib_on_", calculation_result_name)
					} # END if (grepl(pattern="min", x=dfline$samplingRate_function) == TRUE)			
				# Do the probability density
				uniform_id = paste0("UniformDistrib_on_", calculation_result_name)
				uniform_on_scaled_XML = xmlNode(name="Uniform", attrs=list(id=uniform_id, lower=min_of_scaled, upper=max_of_scaled, offset=0, spec="beast.math.distributions.Uniform", name="distr") )
			
				# Do the prior
				scaled_prior_XML = xmlNode(name="prior", attrs=list(id=scaled_prior_name, spec="beast.math.distributions.Prior", name="distribution"), .children=list(eqn_XML, uniform_on_scaled_XML))
				scaled_prior_XMLs = list(bl(), txt_XML, scaled_prior_XML)
			
				# Log the equation result, and the prior
				idref_of_scaled_prior_name = paste0(scaled_prior_name)
				tmplog_XMLtxt = xmlCommentNode(paste0(" Log of the prior probability density of '", scaled_prior_name, "' "))
				tmplog_XML = xmlNode(name="log", attrs=list(idref=idref_of_scaled_prior_name))
				tmplog_XMLs1 = list(bl(), tmplog_XMLtxt, tmplog_XML)

				idref_of_calculation_result_name = paste0(calculation_result_name)
				tmplog_XMLtxt = xmlCommentNode(paste0(" Log of the value of '", calculation_result_name, "' "))
				tmplog_XML = xmlNode(name="log", attrs=list(idref=idref_of_calculation_result_name))
				tmplog_XMLs2 = list(bl(), tmplog_XMLtxt, tmplog_XML)

				tmplog_XMLs = c(tmplog_XMLs1, tmplog_XMLs2)
			
				# Store
				samplingRate_prior_XMLs = c(samplingRate_prior_XMLs, scaled_prior_XMLs)
				samplingRate_log_XMLs = c(samplingRate_log_XMLs, tmpXML$tracelog, tmplog_XMLs)
				} else {
				# Normal distribution on scaled multiplier
				
				# Do the equation
				scaled_samplingRate_equation = paste0("samplingRate[", i-1, "]")
				base_samplingRate_equation = paste0("samplingRate[", base_rownum-1, "]")
				equation_txt = paste0(scaled_samplingRate_equation, " / ", base_samplingRate_equation)
				x_XML = xmlNode(name="x", attrs=list(idref="samplingRate") )
				eqn_XML = xmlNode(name="x", attrs=list(id=calculation_result_name, spec="beast.util.Script", expression=equation_txt), .children=list(x_XML) )
			
				# Do the probability density
				# Mean
				mean_of_scaled = dfline$samplingRate_scaling_relative_to_estimated_base
				mean_of_normal_id = paste0("mean_of_NormalDistrib_on_", calculation_result_name)
				param1_XML = xmlNode(name="parameter", value=mean_of_scaled, attrs=list(id=mean_of_normal_id, name="mean") )
			
				# Standard deviation
				sd_of_scaled = relscale_fraction * mean_of_scaled
				sd_of_normal_id = paste0("sd_of_NormalDistrib_on_", calculation_result_name)
				param2_XML = xmlNode(name="parameter", value=sd_of_scaled, attrs=list(id=sd_of_normal_id, name="sigma") )
			
				# Do the probability density
				normal_id = paste0("NormalDistrib_on_", calculation_result_name)
				normal_on_scaled_XML = xmlNode(name="Normal", attrs=list(id=normal_id, spec="beast.math.distributions.Normal", name="distr"), .children=list(param1_XML, param2_XML) )
			
				# Do the prior
				scaled_prior_XML = xmlNode(name="prior", attrs=list(id=scaled_prior_name, spec="beast.math.distributions.Prior", name="distribution"), .children=list(eqn_XML, normal_on_scaled_XML))
				scaled_prior_XMLs = list(bl(), txt_XML, scaled_prior_XML)
			
				# Log the equation result, and the prior
				idref_of_scaled_prior_name = paste0(scaled_prior_name)
				tmplog_XMLtxt = xmlCommentNode(paste0(" Log of the prior probability density of '", scaled_prior_name, "' "))
				tmplog_XML = xmlNode(name="log", attrs=list(idref=idref_of_scaled_prior_name))
				tmplog_XMLs1 = list(bl(), tmplog_XMLtxt, tmplog_XML)

				idref_of_calculation_result_name = paste0(calculation_result_name)
				tmplog_XMLtxt = xmlCommentNode(paste0(" Log of the value of '", calculation_result_name, "' "))
				tmplog_XML = xmlNode(name="log", attrs=list(idref=idref_of_calculation_result_name))
				tmplog_XMLs2 = list(bl(), tmplog_XMLtxt, tmplog_XML)

				tmplog_XMLs = c(tmplog_XMLs1, tmplog_XMLs2)
			
				# Store
				samplingRate_prior_XMLs = c(samplingRate_prior_XMLs, scaled_prior_XMLs)
				samplingRate_log_XMLs = c(samplingRate_log_XMLs, tmpXML$tracelog, tmplog_XMLs)
				} # END if ( grepl(pattern="scaled", x=dfline$samplingRate_function) == TRUE )
			} # END if ((grepl(pattern="min", x=dfline$samplingRate_function) == TRUE) || (grepl(pattern="max", x=dfline$samplingRate_function) == TRUE))
		
		} # END for (i in 1:num_samplingRates)
	

	birthRate_defn_XMLs
	birthRate_prior_XMLs
	birthRate_log_XMLs
	deathRate_defn_XMLs
	deathRate_prior_XMLs
	deathRate_log_XMLs
	samplingRate_defn_XMLs
	samplingRate_prior_XMLs
	samplingRate_log_XMLs
	#######################################################
	# End absolute and relative priors on samplingRates
	#######################################################




	###############################################
	# rhoSamplingProbs, from oldest to youngest
	###############################################
	rhoSamplingProbs = treemodel_df$rho[isblank_TF(treemodel_df$rho)==FALSE]
	rhoSamplingTimes = treemodel_df$rhoSamplingTimes[isblank_TF(treemodel_df$rhoSamplingTimes)==FALSE]
	
	# Check lengths
	if (length(rhoSamplingProbs) != length(rhoSamplingTimes))
		{
		txt = paste0("\n\nERROR in make_BDSKY_model():\nLength of treemodel_df$rho = ", length(rhoSamplingProbs), "\nLength of treemodel_df$rhoSamplingTimes = ", length(rhoSamplingTimes), "\nThese must be equal in the Excel spreadsheet.\n\n")
		cat(txt)
		stop(txt)
		}
	
	# Ages in youngest to oldest order
	order_indices = order(rhoSamplingTimes)
	
	# We must print to XML in oldest to youngest order
	ages_txt = paste(rhoSamplingTimes[rev(order_indices)], sep=" ", collapse=" ")
	rates_txt = paste(rhoSamplingProbs[rev(order_indices)], sep=" ", collapse=" ")
	rates_txt
	
	# Make the rhoSamplingRate parameters
	rhoSamplingRate_params_id = "rho"
	rhoSamplingRate_params_idref = "@rho"
	rhoSamplingRate_params_XML = xmlNode(name="parameter", attrs=list(id=rhoSamplingRate_params_id, name="rho", lower="0.0", value=rates_txt))
	rhoSamplingTimes_id = "rhoSamplingTimes"
	rhoSamplingTimes_idref = "@rhoSamplingTimes"
	rhoSamplingTimes_XML = xmlNode(name="parameter", attrs=list(id=rhoSamplingTimes_id, name="rhoSamplingTimes", value=ages_txt))
	
	
	#######################################################
	# reverseTimeArrays: which sampling times are reversed
	# (seems to only matter in that samplingRateChangeTimes is true)
	# allegedly, order is times for birth, death, psi, rho, R
	#######################################################
	reverseTimeArrays_id = "reverseTimeArrays"
	reverseTimeArrays_idref = "@reverseTimeArrays"
	reverseTimeArrays_XML = xmlNode(name="reverseTimeArrays", attrs=list(id=reverseTimeArrays_id, spec="beast.core.parameter.BooleanParameter", value="true true true true true") )

	
	
	####################################################
	# Assemble into the BirthDeathSkylineModel node
	####################################################
	kids_XML = list(
	bl(),
	xmlCommentNode(" Starting parameters of the BirthDeathSkylineModel for the observed tree "), 
	
	bl(),
	xmlCommentNode(" originTime: the start of the sampling+tree growth process "), 
	xmlCommentNode(" (This is a nuisance parameter, but its starting value must be older than the start tree root or Beast2 will crash) "), 
	originTime_XML,
	
	bl(),
	xmlCommentNode(" birthRates (lambda: speciation or lineage splitting rates, oldest to youngest) "), 
	birthRate_params_XML,
	birthRateChangeTimes_tops_XML,
	
	bl(),
	xmlCommentNode(" deathRates (mu: lineage extinction rates, oldest to youngest) "), 
	deathRate_params_XML,
	deathRateChangeTimes_tops_XML,
	
	bl(),
	xmlCommentNode(" samplingRates (psi: rate of sampling per lineage per my)                  "), 
	xmlCommentNode(" sorted from youngest to oldest, but with zero at the end,                 "), 
	xmlCommentNode(" as specified by Denise on the beast-users group:                          "), 
	xmlCommentNode(" https://groups.google.com/d/msg/beast-users/UyAS0sYDkKk/PB0ktWnssn0J          "), 
	xmlCommentNode(" The actual psi samplingRates, though, are in order from oldest to youngest "), 
	samplingRate_params_XML,
	samplingRateTimes_tops_XML,
	
	bl(),
	xmlCommentNode(" rho (probability of sampling each lineage extant at a particular timepoint, typically the present, 0 Ma) "), 
	rhoSamplingRate_params_XML,
	rhoSamplingTimes_XML,
	
	bl(),
	xmlCommentNode(" Are time arrays reversed, for birth, death, psi, rho, and R. Regardless, the parameter values are listed from oldest to youngest time bin. "),
	reverseTimeArrays_XML
	) # END kids_XML list
	kids_XML


	# If SABDSKY, add Removal Probability
	if ( treeModel_option == "SABDSKY" )
		{

# 	<!-- removalProbability (probability of a lineage becoming non-infectious) -->
# 	<!-- (this is relevant only in diseases) -->
# 	
# 	<!-- http://blog.beast2.org/2014/09/02/sampled-ancestor-trees-in-beast/ -->
# 	<!-- "For the fossilisation process, r=0 because any sample will never
# 	      be removed, so any of its descendants can be sampled." 
# 	-->
# 	
# 	<parameter id="removalProbability" name="removalProbability" lower="0" upper="1" value="0.0"/>
# 	
		removalProbability_XML = xmlNode(name="parameter", attrs=list(id="removalProbability", name="removalProbability", lower="0", upper="1", value="0.0", estimate="false") )
		
		# Add to kids_XML
		kids_XML = c(kids_XML, 
		list(bl()), 
		list(xmlCommentNode(" removalProbability (probability of a lineage becoming non-infectious) ")), 
		list(xmlCommentNode(" (this is estimated only in diseases) ")), 
		list(bl()), 
		list(xmlCommentNode(" http://blog.beast2.org/2014/09/02/sampled-ancestor-trees-in-beast/ ")), 
		list(xmlCommentNode(" For the fossilisation process, r=0 because any sample will never ")), 
		list(xmlCommentNode(" be removed, so any of its descendants can be sampled. ")), 
		list(removalProbability_XML))
		} # END if ( treeModel_option == "SABDSKY" )


	
	###########################################
	# Set up the tree model
	###########################################
	treeModel_id = tree_name
	treeModel_idref = paste0("@", treeModel_id)
	
	# BDSKY
	if ( treeModel_option == "BDSKY" )
		{
		treeModel_XML = xmlNode(name="BirthDeathSkylineModel", attrs=list(id=treeModel_idref, tree=treeModel_idref, contemp="false", spec="beast.evolution.speciation.BirthDeathSkylineModel"), .children=kids_XML)
	
		treeModel_XMLs = list(bl(), xmlCommentNode(" Set up the BDSKY tree model: BirthDeathSkylineModel "), xmlCommentNode(" (Other models, e.g. Yule (pure-birth), BD (birth-death), BDSS (BD, constant serial sampling) are special cases of BDSKY) "), xmlCommentNode(" (Except the Birth-Death-Skyline-Sampled-Ancestor model) "), treeModel_XML)
		treeModel_XMLs
		} # END if ( treeModel_option == "BDSKY" )



	
	
	# SABDSKY
	if ( treeModel_option == "SABDSKY" )
		{
		treeModel_XML = xmlNode(name="SABDSkylineModel", attrs=list(id=treeModel_idref, tree=treeModel_idref, contemp="false", spec="beast.evolution.speciation.SABDSkylineModel"), .children=kids_XML)
	
		treeModel_XMLs = list(bl(), xmlCommentNode(" Set up the SABDSKY tree model: SABDSkylineModel "), xmlCommentNode(" (Other models, e.g. Yule (pure-birth), BD (birth-death), BDSS (BD, constant serial sampling) are special cases of BDSKY) "), xmlCommentNode(" (Except the Birth-Death-Skyline-Sampled-Ancestor model) "), treeModel_XML)
		treeModel_XMLs
		} # END if ( treeModel_option == "SABDSKY" )
	

	
	


	
	
	#######################################################
	# operators on the tree model parameters
	#######################################################

	# birthRate operators
	birthRate_operator_id = paste0("rateScaler_", birthRate_params_id, "_", tree_name)
	birthRate_operator_XML = xmlNode(name="operator", attrs=list(id=birthRate_operator_id, parameter=birthRate_params_idref, scaleFactor="0.75", weight="10.0", spec="ScaleOperator") )

	# If any of these are "estimated", put in the operator
	birthRate_function = treemodel_df$birthRate_function[isblank_TF(treemodel_df$birthRate_function)==FALSE]
	TF = grepl(pattern="estimated", x=birthRate_function)
	if (any(TF) == FALSE)
		{
		xmltxt = saveXML(doc=birthRate_operator_XML, file=NULL, prefix=NULL)
		birthRate_operator_XML = xmlCommentNode(paste0(" ", xmltxt, " "))
		} # END if (any(TF))
	
	
	# deathRate operators
	deathRate_operator_id = paste0("rateScaler_", deathRate_params_id, "_", tree_name)
	deathRate_operator_XML = xmlNode(name="operator", attrs=list(id=deathRate_operator_id, parameter=deathRate_params_idref, scaleFactor="0.75", weight="10.0", spec="ScaleOperator") )

	# If any of these are "estimated", put in the operator
	deathRate_function = treemodel_df$deathRate_function[isblank_TF(treemodel_df$deathRate_function)==FALSE]
	TF = grepl(pattern="estimated", x=deathRate_function)
	if (any(TF) == FALSE)
		{
		xmltxt = saveXML(doc=deathRate_operator_XML, file=NULL, prefix=NULL)
		deathRate_operator_XML = xmlCommentNode(paste0(" ", xmltxt, " "))
		} # END if (any(TF))
	
	
	# samplingRate operators
	samplingRate_operator_id = paste0("rateScaler_", samplingRate_params_id, "_", tree_name)
	samplingRate_operator_XML = xmlNode(name="operator", attrs=list(id=samplingRate_operator_id, parameter=samplingRate_params_idref, scaleFactor="0.75", weight="10.0", spec="ScaleOperator") )

	# If any of these are "estimated", put in the operator
	samplingRate_function = treemodel_df$samplingRate_function[isblank_TF(treemodel_df$samplingRate_function)==FALSE]
	TF = grepl(pattern="estimated", x=samplingRate_function)
	if (any(TF) == FALSE)
		{
		xmltxt = saveXML(doc=samplingRate_operator_XML, file=NULL, prefix=NULL)
		samplingRate_operator_XML = xmlCommentNode(paste0(" ", xmltxt, " "))
		} # END if (any(TF))
	
	
	# rhoSamplingRate operators
	rhoSamplingRate_operator_id = paste0("rateScaler_", rhoSamplingRate_params_id, "_", tree_name)
	rhoSamplingRate_operator_XML = xmlNode(name="operator", attrs=list(id=rhoSamplingRate_operator_id, parameter=rhoSamplingRate_params_idref, scaleFactor="0.75", weight="10.0", spec="ScaleOperator") )

	# If any of these are "estimated", put in the operator
	rhoSamplingRate_function = treemodel_df$rho_function[isblank_TF(treemodel_df$rho_function)==FALSE]
	TF = ((rhoSamplingRate_function == "estimated") || (rhoSamplingRate_function == "estimated_base"))
	if (any(TF) == FALSE)
		{
		xmltxt = saveXML(doc=rhoSamplingRate_operator_XML, file=NULL, prefix=NULL)
		rhoSamplingRate_operator_XML = xmlCommentNode(paste0(" ", xmltxt, " "))
		} # END if (any(TF))
	
	
	# originTime operators
	originTime_operator_id = paste0("rateScaler_", originTime_id, "_", tree_name)
	originTime_operator_XML = xmlNode(name="operator", attrs=list(id=originTime_operator_id, parameter=originTime_idref, scaleFactor="0.75", weight="1.0", spec="ScaleOperator") )
	
	
	
	
	#######################################################
	# operators on the tree
	#######################################################
	tree_name_idref = paste0("@", tree_name)
	
	#######################################################
	# operators on BirthDeathSkylineModel (BDSKY) tree
	#######################################################
	if ( treeModel_option == "BDSKY" )
		{
		treeScaler_id = paste0("treeScaler_", tree_name)
		treeScaler_XML = xmlNode(name="operator", attrs=list(id=treeScaler_id, tree=tree_name_idref, scaleFactor="0.5", weight="3.0", spec="ScaleOperator") )
	
		treeRootScaler_id = paste0("treeRootScaler_", tree_name)
		treeRootScaler_XML = xmlNode(name="operator", attrs=list(id=treeRootScaler_id, tree=tree_name_idref, scaleFactor="0.5", weight="3.0", spec="ScaleOperator") )
	
		UniformOperator_id = paste0("UniformOperator_", tree_name)
		UniformOperator_XML = xmlNode(name="operator", attrs=list(id=UniformOperator_id, tree=tree_name_idref, weight="30.0", spec="Uniform") )

		SubtreeSlide_id = paste0("SubtreeSlide_", tree_name)
		SubtreeSlide_XML = xmlNode(name="operator", attrs=list(id=SubtreeSlide_id, tree=tree_name_idref, weight="15.0", spec="SubtreeSlide") )
	
		narrowExchange_id = paste0("narrowExchange_", tree_name)
		narrowExchange_XML = xmlNode(name="operator", attrs=list(id=narrowExchange_id, tree=tree_name_idref, weight="15.0", spec="Exchange") )

		wideExchange_id = paste0("wideExchange_", tree_name)
		wideExchange_XML = xmlNode(name="operator", attrs=list(id=wideExchange_id, tree=tree_name_idref, isNarrow="false", weight="3.0", spec="Exchange") )

		WilsonBalding_id = paste0("WilsonBalding_", tree_name)
		WilsonBalding_XML = xmlNode(name="operator", attrs=list(id=WilsonBalding_id, tree=tree_name_idref, weight="3.0", spec="WilsonBalding") )

		# Make the list of tree operators
		tree_operators_XMLs = list(
			bl(),
			xmlCommentNode(" Operators on the tree model parameters "),
			xmlCommentNode(" (commented out unless specified as estimated in Excel settings file) "),
			birthRate_operator_XML,
			deathRate_operator_XML,
			samplingRate_operator_XML,
			rhoSamplingRate_operator_XML,
			originTime_operator_XML,
			bl(),
			xmlCommentNode(" Operators on the phylogeny (topology and node dates; tip dates could be added) "),
			treeScaler_XML,
			treeRootScaler_XML,
			UniformOperator_XML,
			SubtreeSlide_XML,
			narrowExchange_XML,
			wideExchange_XML,
			bl(),
			xmlCommentNode(" Wilson-Balding Operator: for sampled ancestor trees, can change "),
			xmlCommentNode(" the number of branches in the tree: "),
			xmlCommentNode(" https://figshare.com/articles/_The_Wilson_Balding_operator_/1260483 "),
			WilsonBalding_XML,
			bl()
			) # END tree_operators_XMLs
		} # END if ( treeModel_option == "BDSKY" )
	

	if ( treeModel_option == "SABDSKY" )
		{
		#######################################################
		# operators on SampledAncestorBirthDeathSkylineModel (SABDSKY) tree
		#######################################################
		treeDimension_id = paste0("TreeDimensionJump_", tree_name)
		rateCategories_idref = paste0("@rateCategories_of_", clockModel_name)
		removalProbability_idref = paste0("@removalProbability")
		
		# WITH rateCategories in operator -- seems to cause this error under
		# SABD trees, after a few hundred generations:
		# 
		# java.lang.ArrayIndexOutOfBoundsException: -1
		# at beast.evolution.branchratemodel.UCRelaxedClockModel.getRateForBranch(Unknown Source)
		# at beast.evolution.likelihood.TreeLikelihood.traverse(Unknown Source)
		# at beast.evolution.likelihood.TreeLikelihood.traverse(Unknown Source)
		# 
		#treeDimension_XML = xmlNode(name="operator", attrs=list(id=treeDimension_id, tree=tree_name_idref, removalProbability=removalProbability_idref, rateCategories=rateCategories_idref, weight="10.0", spec="TreeDimensionJump") )
		treeDimension_XML = xmlNode(name="operator", attrs=list(id=treeDimension_id, tree=tree_name_idref, removalProbability=removalProbability_idref, weight="10.0", spec="TreeDimensionJump") )

		treeScaler_id = paste0("treeScaler_", tree_name)
		treeScaler_XML = xmlNode(name="operator", attrs=list(id=treeScaler_id, tree=tree_name_idref, scaleFactor="0.95", weight="3.0", spec="ScaleOperatorForZeroBranchSATrees") )
	
		treeRootScaler_id = paste0("treeRootScaler_", tree_name)
		treeRootScaler_XML = xmlNode(name="operator", attrs=list(id=treeRootScaler_id, tree=tree_name_idref, scaleFactor="0.95", weight="1.0", rootOnly="true", spec="ScaleOperatorForZeroBranchSATrees") )
	
		UniformOperator_id = paste0("UniformOperator_", tree_name)
		UniformOperator_XML = xmlNode(name="operator", attrs=list(id=UniformOperator_id, tree=tree_name_idref, weight="20.0", spec="UniformForZeroBranchSATrees") )

	# 	SubtreeSlide_id = paste0("SubtreeSlide_", tree_name)
	# 	SubtreeSlide_XML = xmlNode(name="operator", attrs=list(id=SubtreeSlide_id, tree=tree_name_idref, weight="15.0", spec="SubtreeSlide") )
		
		clock_rate_idref = paste0("@", clock_type, "_Mean_of_", clockModel_name)
		updownTree_id = paste0("updownTree_", tree_name)
		updownTree_XML = xmlNode(name="operator", attrs=list(id=updownTree_id, up=clock_rate_idref, down=tree_name_idref, scaleFactor="0.95", weight="20.0", spec="UpDownOperator") )

	
		narrowExchange_id = paste0("narrowExchange_", tree_name)
		narrowExchange_XML = xmlNode(name="operator", attrs=list(id=narrowExchange_id, tree=tree_name_idref, weight="15.0", spec="Exchange") )

		wideExchange_id = paste0("wideExchange_", tree_name)
		wideExchange_XML = xmlNode(name="operator", attrs=list(id=wideExchange_id, tree=tree_name_idref, isNarrow="false", weight="3.0", spec="Exchange") )

		WilsonBalding_id = paste0("WilsonBalding_", tree_name)
		WilsonBalding_XML = xmlNode(name="operator", attrs=list(id=WilsonBalding_id, tree=tree_name_idref, weight="3.0", spec="WilsonBalding") )

		# Make the list of tree operators
		tree_operators_XMLs = list(
			bl(),
			xmlCommentNode(" Operators on the tree model parameters "),
			xmlCommentNode(" (commented out unless specified as estimated in Excel settings file) "),
			birthRate_operator_XML,
			deathRate_operator_XML,
			samplingRate_operator_XML,
			rhoSamplingRate_operator_XML,
			originTime_operator_XML,
			bl(),
			xmlCommentNode(" Operators on the phylogeny (topology and node dates; tip dates could be added) "),
			xmlCommentNode(" Changing tree dimension suggests moving specimens between tips and direct ancestors "),
			treeDimension_XML,
			treeScaler_XML,
			treeRootScaler_XML,
			UniformOperator_XML,
			updownTree_XML,
			narrowExchange_XML,
			wideExchange_XML,
			bl(),
			xmlCommentNode(" The sampled-ancestor operator, IIRC described in arXiv 2014, is WilsonBaldingForZeroBranchSampledAncestorTrees "),
			WilsonBalding_XML
			) # END tree_operators_XMLs
		} # END if ( treeModel_option == "SABDSKY" )
	
	tree_operators_XMLs
	
	#######################################################
	# XML for the states section of MCMC run
	#######################################################
	
	epidemiology_comment1_XML = xmlCommentNode(" Epidemiology version of a tree model ")
	
    txt2 = '<parameter dimension="10" id="samplingProportionS.t:shared_tree" lower="0.0" name="stateNode" upper="1.0">0.01</parameter>\n        <parameter dimension="10" id="becomeUninfectiousRateS.t:shared_tree" lower="0.0" name="stateNode" upper="10.0">1.0</parameter>\n        <parameter dimension="10" id="R0S.t:shared_tree" lower="0.0" name="stateNode" upper="10.0">2.0</parameter>'
    epidemiology_comment2_XML = xmlCommentNode(txt2)
        
	
	state_of_treeModel_XMLs = list(
	bl(),
	xmlCommentNode(" State of the tree & tree model parameters during MCMC search "), 
	xmlNode(name="stateNode", attrs=list(idref=tree_name) ),
	xmlNode(name="stateNode", attrs=list(idref=originTime_id) ),
	xmlNode(name="stateNode", attrs=list(idref=birthRate_params_id) ),
	xmlNode(name="stateNode", attrs=list(idref=deathRate_params_id) ),
	xmlNode(name="stateNode", attrs=list(idref=samplingRate_params_id) ),
	xmlNode(name="stateNode", attrs=list(idref=birthRateChangeTimes_tops_id) ),
	xmlNode(name="stateNode", attrs=list(idref=deathRateChangeTimes_tops_id) ),
	xmlNode(name="stateNode", attrs=list(idref=samplingRateTimes_tops_id) ),
	xmlNode(name="stateNode", attrs=list(idref=rhoSamplingRate_params_id) ),
	xmlNode(name="stateNode", attrs=list(idref=rhoSamplingTimes_id) ),
	xmlNode(name="stateNode", attrs=list(idref=reverseTimeArrays_id) ),
	bl(),
	epidemiology_comment1_XML,
	epidemiology_comment2_XML
	)
	# Add Remco's (RRB's) tree-copying procedure
	state_of_treeModel_XMLs = c(state_of_treeModel_XMLs, tree_RRM_statenode_XMLs)
	state_of_treeModel_XMLs
	
	#######################################################
	# XML for the parameters of the tree prior
	#######################################################
	# Prior on originTime
	originTime_id
	originTime_prior_id = paste0("prior_", originTime_id)
	originTime_prior_Uniform_id = paste0("prior_", originTime_id, "_Uniform_distrib")
	originTime_prior_Uniform_XML = xmlNode(name="Uniform", attrs=list(id=originTime_prior_Uniform_id, lower="0.0", upper="4600.0", name="distr") )
	originTime_prior_XML = xmlNode(name="prior", attrs=list(id=originTime_prior_id, name="distribution", x=originTime_idref), .children=list(originTime_prior_Uniform_XML) )
	
	OLD_PRIORS='	
	# Prior on birthRate
	birthRate_prior_id = paste0("prior_", birthRate_params_id, "_", tree_name)
	birthRate_prior_logNormal_mean_id = paste0("prior_", birthRate_params_id, "_", tree_name, "_logNormal_mean")
	birthRate_prior_logNormal_mean_XML = xmlNode(name="parameter", 0.1, attrs=list(id=birthRate_prior_logNormal_mean_id, name="M", estimate="false") )
	birthRate_prior_logNormal_stdev_id = paste0("prior_", birthRate_params_id, "_", tree_name, "_logNormal_stdev")
	birthRate_prior_logNormal_stdev_XML = xmlNode(name="parameter", 0.33, attrs=list(id=birthRate_prior_logNormal_stdev_id, name="S", estimate="false") )
	
	birthRate_prior_logNormal_dist_id = paste0("prior_", birthRate_params_id, "_", tree_name, "_logNormal_distrib")
	birthRate_prior_logNormal_dist_XML = xmlNode(name="LogNormal", attrs=list(id=birthRate_prior_logNormal_dist_id, name="distr", meanInRealSpace="true"), .children=list(birthRate_prior_logNormal_mean_XML, birthRate_prior_logNormal_stdev_XML))
	
	birthRate_prior_XML = xmlNode(name="prior", attrs=list(id=birthRate_prior_id, name="distribution", x=birthRate_params_idref), .children=list(birthRate_prior_logNormal_dist_XML) )

	
	# Prior on deathRate
	deathRate_prior_id = paste0("prior_", deathRate_params_id, "_", tree_name)
	deathRate_prior_logNormal_mean_id = paste0("prior_", deathRate_params_id, "_", tree_name, "_logNormal_mean")
	deathRate_prior_logNormal_mean_XML = xmlNode(name="parameter", 0.1, attrs=list(id=deathRate_prior_logNormal_mean_id, name="M", estimate="false") )
	deathRate_prior_logNormal_stdev_id = paste0("prior_", deathRate_params_id, "_", tree_name, "_logNormal_stdev")
	deathRate_prior_logNormal_stdev_XML = xmlNode(name="parameter", 0.33, attrs=list(id=deathRate_prior_logNormal_stdev_id, name="S", estimate="false") )
	
	deathRate_prior_logNormal_dist_id = paste0("prior_", deathRate_params_id, "_", tree_name, "_logNormal_distrib")
	deathRate_prior_logNormal_dist_XML = xmlNode(name="LogNormal", attrs=list(id=deathRate_prior_logNormal_dist_id, name="distr", meanInRealSpace="true"), .children=list(deathRate_prior_logNormal_mean_XML, deathRate_prior_logNormal_stdev_XML))
	
	deathRate_prior_XML = xmlNode(name="prior", attrs=list(id=deathRate_prior_id, name="distribution", x=deathRate_params_idref), .children=list(deathRate_prior_logNormal_dist_XML) )

	
	# Prior on samplingRate (exponential prior)
	samplingRate_prior_id = paste0("prior_", samplingRate_params_id, "_", tree_name)
	samplingRate_prior_Exponential_mean_id = paste0("prior_", samplingRate_params_id, "_", tree_name, "_Exponential_mean")
	samplingRate_prior_Exponential_mean_XML = xmlNode(name="parameter", 0.01, attrs=list(id=samplingRate_prior_Exponential_mean_id, name="mean", estimate="false") )
	
	samplingRate_prior_Exponential_dist_id = paste0("prior_", samplingRate_params_id, "_", tree_name, "_Exponential_distrib")
	samplingRate_prior_Exponential_dist_XML = xmlNode(name="Exponential", attrs=list(id=samplingRate_prior_Exponential_dist_id, name="distr"), .children=list(samplingRate_prior_Exponential_mean_XML))
	
	samplingRate_prior_XML = xmlNode(name="prior", attrs=list(id=samplingRate_prior_id, name="distribution", x=samplingRate_params_idref), .children=list(samplingRate_prior_Exponential_dist_XML) )
	' # END OLD_PRIORS=
	# rho sampling probability is fixed for now
	
	
	
	#######################################################
	# XML for the tree prior
	#######################################################
	tree_prior_id = paste0("BDSKY_prior_", tree_name)
	tree_prior_idref = paste0("@", tree_prior_id)
	tree_name_idref = paste0("@", tree_name)
	
	tree_prior_XML = xmlNode(name="distribution", attrs=list(
	id=tree_prior_id,
	spec="beast.evolution.speciation.BirthDeathSkylineModel",
	tree=tree_name_idref,
	origin=originTime_idref,
	birthRate=birthRate_params_idref,
	birthRateChangeTimes=birthRateChangeTimes_tops_idref,
	deathRate=deathRate_params_idref,
	deathRateChangeTimes=deathRateChangeTimes_tops_idref,
	samplingRate=samplingRate_params_idref,
	samplingRateChangeTimes=samplingRateTimes_tops_idref,
	rho=rhoSamplingRate_params_idref,
	rhoSamplingTimes=rhoSamplingTimes_idref,
	reverseTimeArrays=reverseTimeArrays_idref
	) )
	
	tree_prior_XMLs = c(
	list(
	bl(),
	xmlCommentNode(" Prior probability distributions of tree model parameters "),
	bl(),
	xmlCommentNode(" Uniform prior on the time of origin of the process (assumes my) "),
	xmlCommentNode(" x means whatever parameter you are assessing via the distribution "),
	originTime_prior_XML,
	bl(),
	xmlCommentNode(" Priors on birthRate(s) ")),
	birthRate_prior_XMLs,
	list( bl(),
	xmlCommentNode(" Priors on deathRate(s) ")),
	deathRate_prior_XMLs,
	list( bl(),
	xmlCommentNode(" Prior on samplingRate(s) ")),
	samplingRate_prior_XMLs,
	list(bl(),
	xmlCommentNode(" Also: we are assuming for now that rho at time 0 is fixed "),
	bl(),
	xmlCommentNode(" Probability of MCMC sampled tree given treeModel parameters "),
	tree_prior_XML)
	)
	tree_prior_XMLs




	#######################################################
	# Log the tree prior, parameter values, and parameter priors
	#######################################################
	tracelog_XMLs = list(
	bl(), 
	xmlCommentNode(" Log the probability of the tree given the treeModel parameters"), 
	xmlNode(name="log", attrs=list(idref=tree_prior_id) ),
	bl(), 
	xmlCommentNode(" Log the prior probability of the treeModel parameters"), 
	xmlNode(name="log", attrs=list(idref=originTime_id) ),
	xmlNode(name="log", attrs=list(idref=birthRate_params_id) ),
	xmlNode(name="log", attrs=list(idref=deathRate_params_id) ),
	xmlNode(name="log", attrs=list(idref=samplingRate_params_id) ),
	xmlNode(name="log", attrs=list(idref=rhoSamplingRate_params_id) )		
	)
	tracelog_XMLs
	# Add the new rate parameters
	tracelog_XMLs = c(tracelog_XMLs, birthRate_log_XMLs, deathRate_log_XMLs, samplingRate_log_XMLs)

	#######################################################
	# screenLog the tree prior, parameter values, and parameter priors
	#######################################################
	screenlog_XMLs = list(
	bl(), 
	xmlCommentNode(" Log the probability of the tree given the treeModel parameters"), 
	xmlNode(name="log", attrs=list(idref=tree_prior_id) ),
	bl(), 
	xmlCommentNode(" Log the prior probability of the treeModel parameters"), 
	xmlNode(name="log", attrs=list(idref=originTime_id) ),
	xmlNode(name="log", attrs=list(idref=birthRate_params_id) ),
	xmlNode(name="log", attrs=list(idref=deathRate_params_id) ),
	xmlNode(name="log", attrs=list(idref=samplingRate_params_id) ),
	xmlNode(name="log", attrs=list(idref=rhoSamplingRate_params_id) )
	)
	# Nah don't need this...
	#screenlog_XMLs = c(screenlog_XMLs, birthRate_log_XMLs, deathRate_log_XMLs, samplingRate_log_XMLs)

	
	#######################################################
	# treeLog
	#######################################################
	TreeWithMetaDataLogger_id = paste0("TreeWithMetaDataLogger_", tree_name)
	clockModel_name_idref = paste0("@", clockModel_name)
	TreeWithMetaDataLogger_XML = xmlNode(name="log", attrs=list(id=TreeWithMetaDataLogger_id, tree=tree_name_idref, branchratemodel=clockModel_name_idref, substitutions="false", spec="beast.evolution.tree.TreeWithMetaDataLogger") )
	TreeWithMetaDataLogger_XMLs = list(
	xmlCommentNode(" Log the tree and branch rates "),
	TreeWithMetaDataLogger_XML
	)
	
	#######################################################
	# subsLog -- tree of 
	#######################################################
	Subs_TreeWithMetaDataLogger_id = paste0("Subs_TreeWithMetaDataLogger_", tree_name)
	Subs_TreeWithMetaDataLogger_XML = xmlNode(name="log", attrs=list(id=Subs_TreeWithMetaDataLogger_id, tree=tree_name_idref, substitutions="true", spec="beast.evolution.tree.TreeWithMetaDataLogger") )
	Subs_TreeWithMetaDataLogger_XMLs = list(
	xmlCommentNode(" Log the tree and branch rates "),
	Subs_TreeWithMetaDataLogger_XML
	)
	



	
	if (is.null(xml))
		{
		res = NULL
		res$misc = c(birthRate_defn_XMLs, deathRate_defn_XMLs, samplingRate_defn_XMLs)

		res$startingTree_XMLs = startingTree_XMLs
		res$treeModel_XMLs = treeModel_XMLs
		res$tree_operators_XMLs = tree_operators_XMLs
		res$state_of_treeModel_XMLs = state_of_treeModel_XMLs
		res$tree_prior_XMLs = tree_prior_XMLs
		res$tracelog_XMLs = tracelog_XMLs
		res$screenlog_XMLs = screenlog_XMLs
		res$treelog_XMLs = TreeWithMetaDataLogger_XMLs
		res$subslog_XMLs = Subs_TreeWithMetaDataLogger_XMLs
		
		extract='
		misc_XMLs = res$misc
		startingTree_XMLs = res$startingTree_XMLs
		treeModel_XMLs = res$treeModel_XMLs
		tree_operators_XMLs = res$tree_operators_XMLs
		state_of_treeModel_XMLs = res$state_of_treeModel_XMLs
		tree_prior_XMLs = res$tree_prior_XMLs
		tracelog_XMLs = res$tracelog_XMLs
		screenlog_XMLs = res$screenlog_XMLs
		TreeWithMetaDataLogger_XMLs = res$treelog_XMLs
		Subs_TreeWithMetaDataLogger_XMLs = res$subslog_XMLs
		'
		
		return(res)
		} else {
		xml$misc = c(xml$misc, birthRate_defn_XMLs, deathRate_defn_XMLs, samplingRate_defn_XMLs)
		xml$starting_tree = c(xml$starting_tree, startingTree_XMLs)
		xml$tree = c(xml$tree, treeModel_XMLs)
		xml$state = c(xml$state, state_of_treeModel_XMLs)
		xml$operators = c(xml$operators, tree_operators_XMLs)
		xml$priors = c(xml$priors, tree_prior_XMLs)
		xml$tracelog = c(xml$tracelog, tracelog_XMLs)
		xml$screenlog = c(xml$screenlog, screenlog_XMLs)
		xml$treelog = c(xml$treelog, TreeWithMetaDataLogger_XMLs)
		xml$subslog = c(xml$subslog, Subs_TreeWithMetaDataLogger_XMLs)
		return(xml)
		} # END if (is.null(xml))
	
	cat("\n\n")
	stop("ERROR in make_BDSKY_model(): Shouldn't get here.")
	} # END make_BDSKY_model
	
	
	
	
	
make_taxon_superset <- function(taxonsets_df, speciesTree_taxonset_listname="taxonsuperset")
	{
	# Filter taxonsets_df to the "use" column
	taxonsets_df = readWorksheetFromFile(xlsfn, sheet="taxonsets", startRow=15)
	taxonsets_df$use[isblank_TF(taxonsets_df$use)] = "yes"
	keepTF = (taxonsets_df$use != "no")
	taxonsets_df = taxonsets_df[keepTF, ]
	
	list_of_taxonset_children = list(
	bl(), 
	xmlCommentNode(" Taxon list for the Species Tree "))
	
	# There may be many specimen names for each speciesName!!
	list_of_taxonset_children = list()
	uniq_speciesNames = unique(taxonsets_df$speciesName)
	for (i in 1:length(uniq_speciesNames))
		{
		rows_TF = taxonsets_df$speciesName == uniq_speciesNames[i]
		specimenStrings = taxonsets_df$specimenString[rows_TF]
		specimen_children = list()
		for (j in 1:length(specimenStrings))
			{
			specimenString = specimenStrings[j]
			xml_child = xmlNode(name="taxon", attrs=list(id=specimenString, spec="Taxon") )
			specimen_children = c(specimen_children, list(xml_child))
			}
		xml_specimens_with_a_species = xmlNode(name="taxon", attrs=list(id=uniq_speciesNames[i], spec="TaxonSet"), .children=specimen_children)
		list_of_taxonset_children = c(list_of_taxonset_children, list(xml_specimens_with_a_species))
		} # END for (i in 1:length(uniq_speciesNames))
	taxonset_XML = xmlNode(name="taxonset", attrs=list(id=speciesTree_taxonset_listname, spec="TaxonSet"), .children=list_of_taxonset_children)
	return(taxonset_XML)
	}
	

make_speciesTree_statenode <- function(tree_name="shared_tree", taxonset_XML, starting_tree_string="")
	{
	if (isblank_TF(starting_tree_string) == TRUE)
		{
		# Default: tree will be random
		speciesTree_statenode_XML = xmlNode(name="tree", attrs=list(id=tree_name, name="stateNode"), .children=list(taxonset_XML))
		} else {
		# Use a user-specified starting tree Newick string
		speciesTree_statenode_XML = xmlNode(name="tree", attrs=list(id=tree_name, name="stateNode", initial="NULL", IsLabelledNewick="true", adjustTipHeights="false", estimate="true", spec="beast.util.TreeParser", threshold="0.001", newick=starting_tree_string), .children=list(taxonset_XML))
		}
	return(speciesTree_statenode_XML)
	}


make_speciesTree_statenode2 <- function(tree_name="shared_tree", taxonset_XML)
	{
	speciesTree_statenode_XML = xmlNode(name="stateNode", attrs=list(id=tree_name, spec="starbeast2.SpeciesTree"), .children=list(taxonset_XML))

	return(speciesTree_statenode_XML)
	}
	
	
	
make_speciesTree_starting <- function(tree_name="shared_tree", taxonset_XML, starting_tree_string="")
	{
	if (isblank_TF(starting_tree_string) == TRUE)
		{
		# Default: tree will be random
		speciesTree_statenode_XML = xmlNode(name="tree", attrs=list(id=tree_name, name="tree"), .children=list(taxonset_XML))
		} else {
		# Use a user-specified starting tree Newick string
		speciesTree_statenode_XML = xmlNode(name="tree", attrs=list(id=tree_name, name="tree", initial="NULL", IsLabelledNewick="true", adjustTipHeights="false", estimate="true", spec="beast.util.TreeParser", threshold="0.001", newick=starting_tree_string), .children=list(taxonset_XML))
		}
	return(speciesTree_statenode_XML)
	}
	
	




treeClock_operators <- function(tree_name, clock_name, clock_type)
	{
	# Starter
	operator_id_suffix = paste0("operator_for_TREE_", tree_name, "_CLOCK_", clock_name)
	clockName_id = clock_name
	genetree_id = tree_name
	clockName_idref = paste0("@", clockName_id)
	genetree_idref = paste0("@", genetree_id)
	
	clockRate_id = paste0(clock_type, "_Mean_of_", clockName_id)
	clockRate_idref = paste0("@", clockRate_id)
	
	# Clock operators
	clockRateScaler_id = paste0("clockRateScaler_", operator_id_suffix)
	clockRateScaler_XML = xmlNode(name="operator", attrs=list(id=clockRateScaler_id, spec="ScaleOperator", parameter=clockRate_idref, scaleFactor="0.5", weight="3.0"))
	
	clockUpDownOperator_id = paste0("clockUpDownOperator_", operator_id_suffix)
	up_XML = xmlNode(name="up", attrs=list(idref=clockRate_id))
	down_XML = xmlNode(name="down", attrs=list(idref=genetree_id))
	updown_XML_list = list(up_XML, down_XML)
	clockUpDownOperator_XML = xmlNode(name="operator", attrs=list(id=clockUpDownOperator_id, spec="UpDownOperator", scaleFactor="0.95", weight="3.0"), .children=updown_XML_list)	
	
	# Tree operators
	treeScaler_id = paste0("treeScaler_", operator_id_suffix)
	treeScaler_XML = xmlNode(name="operator", attrs=list(id=treeScaler_id, tree=genetree_idref, scaleFactor="0.95", weight="3.0", spec="ScaleOperator") )

	treeRootScaler_id = paste0("treeRootScaler_", operator_id_suffix)
	treeRootScaler_XML = xmlNode(name="operator", attrs=list(id=treeRootScaler_id, tree=genetree_idref, scaleFactor="0.7", weight="3.0", spec="ScaleOperator", rootOnly="true") )

	UniformOperator_id = paste0("UniformOperator_", operator_id_suffix)
	UniformOperator_XML = xmlNode(name="operator", attrs=list(id=UniformOperator_id, tree=genetree_idref, weight="15.0", spec="Uniform") )

	SubtreeSlide_id = paste0("SubtreeSlide_", operator_id_suffix)
	SubtreeSlide_XML = xmlNode(name="operator", attrs=list(id=SubtreeSlide_id, tree=genetree_idref, weight="15.0", spec="SubtreeSlide", size="0.002") )

	narrowExchange_id = paste0("narrowExchange_", operator_id_suffix)
	narrowExchange_XML = xmlNode(name="operator", attrs=list(id=narrowExchange_id, tree=genetree_idref, weight="15.0", spec="Exchange", isNarrow="true") )

	wideExchange_id = paste0("wideExchange_", operator_id_suffix)
	wideExchange_XML = xmlNode(name="operator", attrs=list(id=wideExchange_id, tree=genetree_idref, isNarrow="false", weight="15.0", spec="Exchange") )

	WilsonBalding_id = paste0("WilsonBalding_", operator_id_suffix)
	WilsonBalding_XML = xmlNode(name="operator", attrs=list(id=WilsonBalding_id, tree=genetree_idref, weight="15.0", spec="WilsonBalding") )

	# Make the list of tree operators
	treeClock_operators_XMLs = list(bl(),
		xmlCommentNode(paste0(" Operators on the geneTree clock/tree model parameters for: ", operator_id_suffix, " ")),
		xmlCommentNode(" (commented out unless specified as estimated in Excel settings file) "),
		bl(),
		xmlCommentNode(" Operators on the clock model "),
		clockRateScaler_XML, 
		clockUpDownOperator_XML,
		xmlCommentNode(" Operators on the tree/phylogeny (topology and node dates; tip dates could be added) "),
		treeScaler_XML,
		treeRootScaler_XML,
		UniformOperator_XML,
		SubtreeSlide_XML,
		narrowExchange_XML,
		wideExchange_XML,
		bl(),
		xmlCommentNode(" Wilson-Balding Operator: for sampled ancestor trees, can change "),
		xmlCommentNode(" the number of branches in the tree: "),
		xmlCommentNode(" https://figshare.com/articles/_The_Wilson_Balding_operator_/1260483 "),
		WilsonBalding_XML,
		bl()
		) # END treeClock_operators_XMLs
		
	return(treeClock_operators_XMLs)
	} # END treeclock_operators <- function(tree_name, clock_name)
	

	

genetrees_operators <- function(seqs_df, taxonsets_df, treemodel_df, tree_name="shared_tree", speciesTree_taxonset_listname="taxonsuperset", speciesTree_default="speciesTree", xml=NULL, birthRate_id="netDiversificationRate")
	{
	defaults='
	speciesTree_taxonset_listname="taxonsuperset"
	taxonsets_df = readWorksheetFromFile(xlsfn, sheet="taxonsets", startRow=15)
	seqs_df = readWorksheetFromFile(xlsfn, sheet="data", startRow=15)
	birthRate_id="netDiversificationRate"
	'
	
	tree_name_idref = paste0("@", tree_name)
	
	# Check for BEASTmasteR-programmed population models
	if (is.null(treemodel_df$popModel))
		{
		treemodel_df$popModel = "constant"
		} # END if (is.null(treemodel_df$popModel))
	if (isblank_TF(treemodel_df$popModel))
		{
		treemodel_df$popModel = "constant"
		} # END if (is.null(treemodel_df$popModel))
	popModel = treemodel_df$popModel

	if (popModel == "constant")
		{
		popModel_prefix = "constantPopModel"
		}
	if (popModel == "analytic")
		{
		popModel_prefix = "constantPopIOModel"
		}
	if ((treemodel_df$popModel == "linear"))
		{
		popModel_prefix = "lwcrPopModel"
		}
	
	popModel_id = paste0(popModel_prefix, "_for_", speciesTree_default)
	popModel_idref = paste0("@", popModel_id)
	popSizes_id = paste0(popModel_prefix, "_popSizes.Species")
	popSizes_idref = paste0("@", popSizes_id)
	popSizesMean_id = paste0(popModel_prefix, "_popSizesMean.Species")
	popSizesShape_id = paste0(popModel_prefix, "_popSizesShape.Species")
	popSizes_idref = paste0("@", popModel_prefix, "_popSizes.Species")
	popSizesMean_idref = paste0("@", popModel_prefix, "_popSizesMean.Species")
	popSizesShape_idref = paste0("@", popModel_prefix, "_popSizesShape.Species")
	popModelBridge_id = "popModelBridge.Species"
			
	tip_popSizes_id = paste0("tip_", popModel_prefix, "_popSizes.Species")
	tip_popSizes_idref = paste0("@", popSizes_id)
	top_popSizes_id = paste0("top_", popModel_prefix, "_popSizes.Species")
	top_popSizes_idref = paste0("@", popSizes_id)
	lwcrPopScale_id = paste0(popModel_prefix, "_popScale.Species")
	lwcrPopScale_idref = paste0("@", lwcrPopScale_id)
	
	
	# Make the stateNodes and operators for the genetree / species tree
	
	# Get the gene tree names
	gene_tree_names = unique(seqs_df$geneTreeName)
	rownums = match(x=gene_tree_names, table=seqs_df$geneTreeName)
	
	tree_XMLs = list()
	genetree_idref_XMLs = list()
	geneTree_idref_XMLs = list()
	for (i in 1:length(gene_tree_names))
		{
		gene_tree_name = gene_tree_names[i]
		alignment_source_id = seqs_df$datasetName[rownums[i]]
		
		child = xmlNode(name="alignment", attrs=list(idref=alignment_source_id))
		taxonset_id = paste0("taxonset_", gene_tree_name)
		taxonset_XML = xmlNode(name="taxonset", attrs=list(id=taxonset_id, spec="TaxonSet"), .children=list(child))
		
		tree_XML = xmlNode(name="tree", attrs=list(id=gene_tree_name, name="stateNode"), .children=list(taxonset_XML))
		tree_XMLs = c(tree_XMLs, list(tree_XML))
		
		# Operators
		genetree_idref_XML = xmlNode(name="genetree", attrs=list(idref=gene_tree_name))
		genetree_idref_XMLs = c(genetree_idref_XMLs, list(genetree_idref_XML))

		# Some like "genetree", some like "geneTree"
		geneTree_idref_XML = xmlNode(name="geneTree", attrs=list(idref=gene_tree_name))
		geneTree_idref_XMLs = c(geneTree_idref_XMLs, list(geneTree_idref_XML))
		}
	tree_XMLs
	
	
	
	#######################################################
	# Operators for the gene trees
	#######################################################
	operator_name = paste0("Reheight_operator_for_", tree_name)
	taxonset_idref = paste0("@", speciesTree_taxonset_listname)
	tree_idref = paste0("@", tree_name)
	weight = 8 * length(gene_tree_names)	# Rough weight, Nick's guess
	NodeReheight_operator_XML = xmlNode(name="operator", attrs=list(id=operator_name, spec="NodeReheight", taxonset=taxonset_idref, tree=tree_idref, weight=weight), .children=genetree_idref_XMLs)

	operator_name = paste0("coordinatedUniform_operator_for_", tree_name)
	species_tree_name_idref = paste0("@", tree_name)
	speciesTree_idref = paste0("@speciesTree")
	weight = 2 * length(gene_tree_names)	# Rough weight, Nick's guess
	
	if ((popModel == "constant") || (popModel == "linear")  || (popModel == "analytic"))
		{
		coordinatedUniform_operator_XML = xmlNode(name="operator", attrs=list(id=operator_name, spec="starbeast2.CoordinatedUniform", speciesTree=tree_name_idref, weight=weight), .children=geneTree_idref_XMLs)
		} else {
		coordinatedUniform_operator_XML = xmlNode(name="operator", attrs=list(id=operator_name, spec="starbeast2.CoordinatedUniform", speciesTree=speciesTree_idref, tree=species_tree_name_idref, weight=weight), .children=geneTree_idref_XMLs)
		}


	operator_name = paste0("coordinatedExponential_operator_for_", tree_name)
	species_tree_name_idref = paste0("@", tree_name)
	speciesTree_idref = paste0("@speciesTree")
	weight = 2 * length(gene_tree_names)	# Rough weight, Nick's guess
	if ((popModel == "constant") || (popModel == "linear")  || (popModel == "analytic"))
		{
		coordinatedExponential_operator_XML = xmlNode(name="operator", attrs=list(id=operator_name, spec="starbeast2.CoordinatedExponential", speciesTree=tree_name_idref, weight=weight), .children=geneTree_idref_XMLs)
		} else {
		coordinatedExponential_operator_XML = xmlNode(name="operator", attrs=list(id=operator_name, spec="starbeast2.CoordinatedExponential", speciesTree=speciesTree_idref, tree=species_tree_name_idref, weight=weight), .children=geneTree_idref_XMLs)
		}


	if ((popModel == "constant") || (popModel == "linear")  || (popModel == "analytic"))
		{	
		operator_name = paste0("TreeScaler_operator_for_", tree_name)
		species_tree_name_idref = paste0("@", tree_name)
		weight = 0.5 * length(gene_tree_names)	# Rough weight, Nick's guess
		TreeScaler_operator_XML = xmlNode(name="operator", attrs=list(id=operator_name, spec="ScaleOperator", scaleFactor="0.95", tree=tree_name_idref, weight=weight))

		operator_name = paste0("TreeRootScaler_operator_for_", tree_name)
		species_tree_name_idref = paste0("@", tree_name)
		weight = 0.5 * length(gene_tree_names)	# Rough weight, Nick's guess
		TreeRootScaler_operator_XML = xmlNode(name="operator", attrs=list(id=operator_name, spec="ScaleOperator", rootOnly="true", scaleFactor="0.7", tree=tree_name_idref, weight=weight))
	
		operator_name = paste0("UniformOperator_operator_for_", tree_name)
		species_tree_name_idref = paste0("@", tree_name)
		weight = 2 * length(gene_tree_names)	# Rough weight, Nick's guess
		UniformOperator_operator_XML = xmlNode(name="operator", attrs=list(id=operator_name, spec="Uniform", tree=tree_name_idref, weight=weight))
	
		operator_name = paste0("SubtreeSlide_operator_for_", tree_name)
		species_tree_name_idref = paste0("@", tree_name)
		weight = 2 * length(gene_tree_names)	# Rough weight, Nick's guess
		SubtreeSlide_operator_XML = xmlNode(name="operator", attrs=list(id=operator_name, spec="SubtreeSlide", size="0.002", tree=tree_name_idref, weight=weight))

		operator_name = paste0("Narrow_operator_for_", tree_name)
		species_tree_name_idref = paste0("@", tree_name)
		weight = 2 * length(gene_tree_names)	# Rough weight, Nick's guess
		Narrow_operator_XML = xmlNode(name="operator", attrs=list(id=operator_name, spec="Exchange", isNarrow="true", tree=tree_name_idref, weight=weight))
		
		operator_name = paste0("Wide_operator_for_", tree_name)
		species_tree_name_idref = paste0("@", tree_name)
		weight = 2 * length(gene_tree_names)	# Rough weight, Nick's guess
		Wide_operator_XML = xmlNode(name="operator", attrs=list(id=operator_name, spec="Exchange", isNarrow="false", tree=tree_name_idref, weight=weight))

		operator_name = paste0("WilsonBalding_operator_for_", tree_name)
		species_tree_name_idref = paste0("@", tree_name)
		weight = 2 * length(gene_tree_names)	# Rough weight, Nick's guess
		WilsonBalding_operator_XML = xmlNode(name="operator", attrs=list(id=operator_name, spec="WilsonBalding", tree=tree_name_idref, weight=weight))
		} else {
		# Old method
		operator_name = paste0("TreeScaler_operator_for_", tree_name)
		species_tree_name_idref = paste0("@", tree_name)
		weight = 0.5 * length(gene_tree_names)	# Rough weight, Nick's guess
		TreeScaler_operator_XML = xmlNode(name="operator", attrs=list(id=operator_name, spec="ScaleOperator", scaleFactor="0.95", tree=species_tree_name_idref, weight=weight))

		operator_name = paste0("TreeRootScaler_operator_for_", tree_name)
		species_tree_name_idref = paste0("@", tree_name)
		weight = 0.5 * length(gene_tree_names)	# Rough weight, Nick's guess
		TreeRootScaler_operator_XML = xmlNode(name="operator", attrs=list(id=operator_name, spec="ScaleOperator", rootOnly="true", scaleFactor="0.7", tree=species_tree_name_idref, weight=weight))
	
		operator_name = paste0("UniformOperator_operator_for_", tree_name)
		species_tree_name_idref = paste0("@", tree_name)
		weight = 2 * length(gene_tree_names)	# Rough weight, Nick's guess
		UniformOperator_operator_XML = xmlNode(name="operator", attrs=list(id=operator_name, spec="Uniform", tree=species_tree_name_idref, weight=weight))
	
		operator_name = paste0("SubtreeSlide_operator_for_", tree_name)
		species_tree_name_idref = paste0("@", tree_name)
		weight = 2 * length(gene_tree_names)	# Rough weight, Nick's guess
		SubtreeSlide_operator_XML = xmlNode(name="operator", attrs=list(id=operator_name, spec="SubtreeSlide", size="0.002", tree=species_tree_name_idref, weight=weight))

		operator_name = paste0("Narrow_operator_for_", tree_name)
		species_tree_name_idref = paste0("@", tree_name)
		weight = 2 * length(gene_tree_names)	# Rough weight, Nick's guess
		Narrow_operator_XML = xmlNode(name="operator", attrs=list(id=operator_name, spec="Exchange", isNarrow="true", tree=species_tree_name_idref, weight=weight))
		
		operator_name = paste0("Wide_operator_for_", tree_name)
		species_tree_name_idref = paste0("@", tree_name)
		weight = 2 * length(gene_tree_names)	# Rough weight, Nick's guess
		Wide_operator_XML = xmlNode(name="operator", attrs=list(id=operator_name, spec="Exchange", isNarrow="false", tree=species_tree_name_idref, weight=weight))

		operator_name = paste0("WilsonBalding_operator_for_", tree_name)
		species_tree_name_idref = paste0("@", tree_name)
		weight = 2 * length(gene_tree_names)	# Rough weight, Nick's guess
		WilsonBalding_operator_XML = xmlNode(name="operator", attrs=list(id=operator_name, spec="WilsonBalding", tree=species_tree_name_idref, weight=weight))
		}


	# UpDown for all clockrates, and gene trees
	updown_XMLs_list = list()
	clock_rownums = match(x=clockModel_names, table=seqs_df$clockModel_name)
	for (i in 1:length(clockModel_names))
		{
		clock_type = seqs_df$clockmodel_type[clock_rownums[i]]
		
		if (grepl(pattern="SpeciesTree", x=clock_type) == TRUE)
			{
			# Collect the geneTrees within this SpeciesTree clock model

			# SPECIESTREE RELAXED CLOCK STARBEAST2 ANALYSIS
			
			TFs = seqs_df$clockModel_name == clockModel_names[i]
			
			for (j in 1:sum(TFs))
				{
				tree_name = seqs_df$geneTreeName[TFs][j]
				geneTree_clockModel_rate_id = paste0(seqs_df$clockmodel_type[clock_rownums[i]], "_Mean_of_clockRate_for_geneTree_", tree_name)
				geneTree_clockModel_rate_idref = paste0("@", geneTree_clockModel_rate_id)
			
				clock_rate_id = paste0(seqs_df$clockmodel_type[clock_rownums[i]], "_", "Mean_of_clockRate_", clockModel_names[i])

				up_XML = xmlNode(name="up", attrs=list(idref=clock_rate_id))
				updown_XMLs_list = c(updown_XMLs_list, list(up_XML))
				} # END for (j in 1:sum(TFs))
			} else {
			
			# REGULAR STARBEAST2 ANALYSIS
			clock_rate_id = paste0(seqs_df$clockmodel_type[clock_rownums[i]], "_Mean_of_", clockModel_names[i])

			up_XML = xmlNode(name="up", attrs=list(idref=clock_rate_id))
			updown_XMLs_list = c(updown_XMLs_list, list(up_XML))
			} # END if (grepl(pattern="SpeciesTree", x=clock_type) == TRUE)
		} # END for (i in 1:length(clockModel_names))
	
	# Add netDiversificationRate
	up_XML = xmlNode(name="up", attrs=list(idref=birthRate_id))
	updown_XMLs_list = c(updown_XMLs_list, list(up_XML))
	
	# Down moves
	down_XML = xmlNode(name="down", attrs=list(idref=tree_name))
	updown_XMLs_list = c(updown_XMLs_list, list(down_XML))

	for (i in 1:length(gene_tree_names))
		{
		down_XML = xmlNode(name="down", attrs=list(idref=gene_tree_names[i]))
		updown_XMLs_list = c(updown_XMLs_list, list(down_XML))
		}
	
	if ((treemodel_df$popModel == "constant"))
		{
		down_XML = xmlNode(name="down", attrs=list(idref=popSizes_id))
		updown_XMLs_list = c(updown_XMLs_list, list(down_XML))

		down_XML = xmlNode(name="down", attrs=list(idref=popSizesMean_id))
		updown_XMLs_list = c(updown_XMLs_list, list(down_XML))
		}

	if ((treemodel_df$popModel == "linear"))
		{
		down_XML = xmlNode(name="down", attrs=list(idref=tip_popSizes_id))
		updown_XMLs_list = c(updown_XMLs_list, list(down_XML))

		down_XML = xmlNode(name="down", attrs=list(idref=top_popSizes_id))
		updown_XMLs_list = c(updown_XMLs_list, list(down_XML))

		down_XML = xmlNode(name="down", attrs=list(idref=lwcrPopScale_id))
		updown_XMLs_list = c(updown_XMLs_list, list(down_XML))
		}

	
	weight = 1 * length(gene_tree_names)
	updown_all_species_XML = xmlNode(name="operator", attrs=list(id="updown.all.Species", spec="UpDownOperator", scaleFactor="0.75", weight=weight), .children=updown_XMLs_list)
		
	

	
	# Scalers for clock rates and gene trees
	# Do it for every unique combination of gene tree and clock rate
	if (grepl(pattern="SpeciesTree", x=clock_type) == TRUE)
		{
		clockModel_names = paste0("clockRate_for_geneTree_", seqs_df$geneTreeName)
		genetree_clockrates_fused = paste0(seqs_df$geneTreeName, "_", clockModel_names)
		} else {
		genetree_clockrates_fused = paste0(seqs_df$geneTreeName, "_", seqs_df$clockModel_name)
		clockModel_names = seqs_df$clockModel_name
		} # END if (grepl(pattern="SpeciesTree", x=clock_type) == TRUE)
	
	genetree_clockrates_fused_unique = unique(genetree_clockrates_fused)
	rownums = match(x=genetree_clockrates_fused_unique, table=genetree_clockrates_fused)
	treeClock_operators_XMLs = list()
	previous_clock_type = NA
	for (r in rownums)
		{
		if (grepl(pattern="SpeciesTree", x=clock_type) == TRUE)
			{
			tree_name = seqs_df$geneTreeName[r]
			clock_name = paste0("clockRate_for_geneTree_", tree_name)
			} else {
			clock_name = seqs_df$clockModel_name[r]
			tree_name = seqs_df$geneTreeName[r]
			} # END if (grepl(pattern="SpeciesTree", x=clock_type) == TRUE)
		clock_type = seqs_df$clockmodel_type[r]
		
		# Store the previous clock_type, use if this one is NA
		if (isblank_TF(clock_type) == FALSE)
			{
			previous_clock_type = clock_type
			} else {
			clock_type = previous_clock_type
			}
		# Error check
		if (isblank_TF(clock_type) == TRUE)
			{
			txt = paste0("STOP ERROR in genetrees_operators(). For clock_name '", clock_name, "', tree_name '", tree_name, "', the clock_type variable came out as blank/NA. Please check the 'clockmodel_type' column in the 'data' worksheet of the Excel settings file.")
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			}
		
		treeClock_operators_XML = treeClock_operators(tree_name=tree_name, clock_name=clock_name, clock_type=clock_type)
		treeClock_operators_XMLs = c(treeClock_operators_XMLs, treeClock_operators_XML)
		}
	
	treeClock_operators_XMLs[1:17]
	treeClock_operators_XMLs[(length(treeClock_operators_XMLs)-17+1):length(treeClock_operators_XMLs)]
	length(treeClock_operators_XMLs)
	
	length(treeClock_operators_XMLs)/17



	# Put the tree stuff in the stateNodes

	genetree_operators_XMLs = c(list(NodeReheight_operator_XML),
	list(coordinatedUniform_operator_XML),
	list(coordinatedExponential_operator_XML),
	list(TreeScaler_operator_XML),
	list(TreeRootScaler_operator_XML),
	list(UniformOperator_operator_XML),
	list(SubtreeSlide_operator_XML),
	list(Narrow_operator_XML),
	list(Wide_operator_XML),
	list(WilsonBalding_operator_XML),
	list(updown_all_species_XML),
	treeClock_operators_XMLs)

	if (is.null(xml))
		{
		res = NULL
		res$tree_XMLs = tree_XMLs
		res$genetree_operators_XMLs = genetree_operators_XMLs

		extract='
		tree_XMLs = res$tree_XMLs
		genetree_operators_XMLs = res$genetree_operators_XMLs
		'
		return(res)
		} else {
		xml$state = c(xml$state, tree_XMLs)
		xml$operators = c(xml$operators, genetree_operators_XMLs)
		return(xml)
		} # END if (is.null(xml))
	} # END genetrees_operators <- function(seqs_df, taxonsets_df, tree_name="shared_tree", speciesTree_taxonset_listname="taxonsuperset", speciesTree_default="speciesTree", xml=NULL)


	
	