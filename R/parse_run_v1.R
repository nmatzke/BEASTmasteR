

#######################################################
# Is the XML item an XMLCommentNode?
#######################################################
is_XMLcomment_TF <- function(item)
	{
	classes = class(item)
	tmpTF = classes == "XMLCommentNode"
	if (sum(tmpTF) > 0)
		{
		return(TRUE)
		} else {
		return(FALSE)
		}
	} # END is_XMLcomment_TF <- function(item)




#######################################################
# Remove duplicate XMLS
#######################################################
remove_duplicate_xmls <- function(xml)
	{
	XML_names = names(xml)
	for (i in 1:length(XML_names))
		{
		#print("Here")
		#print(i)
		#print(class(xml[[i]]))
		# Only proceed if it's a list
		# Extract this XML name
		cmdtxt = paste0("tmpXML = xml$", XML_names[i])
		eval(parse(text=cmdtxt))
		
		if (class(tmpXML) == "list")
			{
			x=1 # pass
			} else {
			# If not
			next()
			}
		
	
		# Proceed
		
		# Remove comments and NULLs
		comment_TF = sapply(X=tmpXML, FUN=is_XMLcomment_TF)
		null_TF = sapply(X=tmpXML, FUN=is.null)
		comment_TF2 = comment_TF + null_TF
		tmpXML2 = tmpXML[comment_TF2 == FALSE]
		
		# Of the remaining, which ones are duplicates?
		duplicated_TF = duplicated(tmpXML2)
		
		if (sum(duplicated_TF) == 0)
			{
			# No duplicates
			newXML = tmpXML
			} else {
			# Find the unique versions of the duplicates
			unique_xmls = unique(tmpXML2[duplicated_TF])

			for (j in 1:length(unique_xmls))
				{
				# Check for identical matches
				matches_TF = sapply(X=tmpXML, FUN=identical, y=unique_xmls[[j]])
				matchnums = (1:length(matches_TF))[matches_TF]
				matchnums_to_remove = matchnums[-1] # keep the first instance
				matchnums_to_remove = matchnums_to_remove * -1 # minus removes them
				tmpXML = tmpXML[matchnums_to_remove]
				} # END for (j in 1:length(unique_xmls))

			newXML = tmpXML
			} # END if (sum(duplicated_TF) == 0)
	
		# Back into original XML
		cmdtxt = paste0("xml$", XML_names[i], " = NULL")
		eval(parse(text=cmdtxt))

		cmdtxt = paste0("xml$", XML_names[i], " = newXML")
		eval(parse(text=cmdtxt))
		} # END for (i in 1:length(XML_names))
	
	return(xml)
	} # END remove_duplicate_xmls <- function(xml)




#######################################################
# Save the stats to default filenames
#######################################################
# Specifically, saves the completeness statistics 
# for morphological data
#
# And the overall morphological matrix summary stats
# 
# And the "morphList"
#
save_data_stats <- function(data_XMLs, datastats_dir=NULL)
	{
	# Output filenames
	outfns = NULL
	
	# Print stats out
	names(data_XMLs)
	
	if (length(data_XMLs$matrices_stats) == 0)
		{
		txt = "WARNING from save_data_stats(): Nothing was found in data_XMLs$matrices_stats. Probably this means your dataset has no morphological characters. No statistics will be calculated on 0 characters. save_data_stats() will now return NULL." 
		cat("\n\n")
		warning(txt)
		cat("\n\n")
		return(NULL)
		}
	
	
	# Modify datastats_dir
	if (is.null(datastats_dir) == FALSE)
		{
		datastats_dir = paste0(slashslash(datastats_dir), "/")
		}
	
	# Show the stats
	data_XMLs$matrices_stats

	# Show the number of data matrices
	length(data_XMLs$charsdf_list)

	for (i in 1:length(data_XMLs$matrices_stats))
		{
		completeness_df = data_XMLs$matrices_stats[[i]]$completeness_df
		outfn = paste0(datastats_dir, "completeness_stats_matrix", i, ".txt")
		write.table(completeness_df, file=outfn, quote=FALSE, sep="\t")
		outfns = c(outfns, outfn)

		matrix_stats_df = data_XMLs$matrices_stats[[i]]$matrix_stats_df
		outfn = paste0(datastats_dir, "stats_matrix", i, ".txt")
		write.table(matrix_stats_df, file=outfn, quote=FALSE, sep="\t")
		outfns = c(outfns, outfn)
		
		if (is.null(data_XMLs$matrices_stats[[i]]$autapomorphies_desc_df) == FALSE)
			{
			autapomorphies_desc_df = data_XMLs$matrices_stats[[i]]$autapomorphies_desc_df
			
			# Write autapomorphies description matrix to file
			outfn = paste0(datastats_dir, "autapomorphies_desc_morph", i, ".txt")
			write.table(autapomorphies_desc_df, file=outfn, quote=FALSE, sep="\t")
			outfns = c(outfns, outfn)
			} # END if (is.null(data_XMLs$matrices_stats[[i]]$char_is_autapomorphic) == FALSE)
		
		} # END for (i in 1:length(data_XMLs$matrices_stats))

	morphList = data_XMLs$morphList
	outfn = paste0(datastats_dir, "morphList.txt")
	outfns = c(outfns, outfn)
	write.table(morphList, file=outfn, quote=FALSE, sep="\t")
	
	txt = paste(outfns, sep="", collapse="\n")
	cat("\n")
	cat("Saved morphology stats to: \n\n", txt)
	cat("\n")
	
	return(outfns)
	} # END save_data_stats <- function(data_XMLs)





#######################################################
# Parse the "run" worksheet, and convert "xml" to the full XML file
#######################################################

parse_run <- function(run_df, xml, outfn="run_df", burnin_fraction=0.5, dataset_source="", printall="short", BEAST2_version="2.4.3")
	{
	defaults='
	outfn=NULL
	outfn="run_df"
	printall="short"
	'
	
	# Error check that only 1 set of setting is used
	TF = run_df$use == "yes"
	if (sum(TF) != 1)
		{
		errortxt = paste0("\n\nERROR in parse_run(): worksheet 'run' must have only 1 line where 'use'=='yes'. Instead, you have ", sum(TF), " such lines.\n")
		cat(errortxt)
		cat("\n")
		stop(errortxt)
		}

	#######################################################
	# Remove any duplicates that were accidentally created
	# (excluding XMLCommentNodes, blanks, etc.)
	#######################################################
	xml = remove_duplicate_xmls(xml)


	
	# Get the line you are using
	settings = run_df[TF,]
	
	# Fill in blanks
	if (isblank_TF(settings$chainLength))
		{
		settings$chainLength = 10000
		}
	if (isblank_TF(settings$numsamps))
		{
		settings$numsamps = 2000
		}
	if (isblank_TF(settings$state_store))
		{
		settings$state_store = floor(settings$chainLength / settings$numsamps)
		}
	if (isblank_TF(settings$tracelog_store))
		{
		settings$tracelog_store = floor(settings$chainLength / settings$numsamps)
		}
	if (isblank_TF(settings$screenlog_store))
		{
		settings$screenlog_store = floor(settings$chainLength / settings$numsamps)
		}
	if (isblank_TF(settings$treelog_store))
		{
		settings$treelog_store = floor(settings$chainLength / settings$numsamps)
		}
	
	
	# Error correction (/0, divide-by-zero error, if logEvery is 0)
	if (settings$state_store < 1) 
		{
		settings$state_store = 100
		}
	if (settings$tracelog_store < 1)
		{
		settings$tracelog_store = 100
		}
	if (settings$screenlog_store < 1)
		{
		settings$screenlog_store = 100
		}
	if (settings$treelog_store < 1)
		{
		settings$treelog_store = 100
		}
	
	# Extract the filename from run_df
	if (outfn == "run_df")
		{
		outfn = settings$outfn
		
		if (isblank_TF(outfn))
			{
			txt = paste("\n\nERROR in parse_run(): You specified outfn='run_df' (the default), which means\n  that outfn (=output filename) will be read from the 'run' worksheet of the Excel\n  setting file. However, the column 'outfn' is blank for the row you are using. To fix:\n(1) Fill it in, (2) set outfn= to an actual filename, or (3) set outfn=NULL to just return\n  the XML to R.\n\n", sep="")
			cat(txt)
			stop(txt)
			}
		}

	
	# Make the header XML_namespaces
	# (XML_namespaces will be included in the <beast> tag)
	XML_namespaces = "beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood"
	
	# Make the map XML
	map_XML_list = list(
	bl(),
	xmlCommentNode(" A list of class mappings used in this Beast2 XML "),
	xmlNode(name="map", "beast.evolution.sitemodel.SiteModel", attrs=list(name="SiteModel") ),
	xmlNode(name="map", "beast.math.distributions.Beta", attrs=list(name="Beta") ),
	xmlNode(name="map", "beast.math.distributions.Exponential", attrs=list(name="Exponential") ),
	xmlNode(name="map", "beast.math.distributions.InverseGamma", attrs=list(name="InverseGamma") ),
	xmlNode(name="map", "beast.math.distributions.LogNormalDistributionModel", attrs=list(name="LogNormal") ),
	xmlNode(name="map", "beast.math.distributions.Gamma", attrs=list(name="Gamma") ),
	xmlNode(name="map", "beast.math.distributions.Uniform", attrs=list(name="Uniform") ),
	xmlNode(name="map", "beast.math.distributions.Prior", attrs=list(name="prior") ),
	xmlNode(name="map", "beast.math.distributions.LaplaceDistribution", attrs=list(name="LaplaceDistribution") ),
	xmlNode(name="map", "beast.math.distributions.OneOnX", attrs=list(name="OneOnX") ),
	xmlNode(name="map", "beast.math.distributions.Normal", attrs=list(name="Normal") )
	) # END map_XML_list
	xml$header = c(xml$header, map_XML_list)
	
	
	#######################################################
	# Make a note about the input dataset
	# (e.g., reference, source, etc.)
	#######################################################
	
	code_source1 = paste0(" This XML file was written using BEASTmasteR at ", Sys.time())
	code_source2 = " BEASTmasteR took a lot of work.  Please cite if you use! Citation info: "
	code_source3 = ' Matzke, Nicholas J. (2014). "BEASTmasteR: Automated conversion of NEXUS data to BEAST2 XML format, for fossil tip-dating and other uses." Online at PhyloWiki, http://phylo.wikidot.com/beastmaster. Accessed (access date). ' 
	
	txt1 = " Information & citation on the data source(s) that went into this XML: "
	txt2 = " (it is crucial to track, and cite, your data sources!) "
	
	citation_XMLs = list(
	bl(),
	bl(),
	xmlCommentNode(code_source1),
	bl(),
	xmlCommentNode(code_source2),
	xmlCommentNode(code_source3),
	bl(),
	xmlCommentNode(txt1),
	xmlCommentNode(txt2),
	bl(),
	xmlCommentNode(dataset_source),
	bl()
	)
	
	xml$header = c(xml$header, citation_XMLs)
	
	
	# Assemble the prior, likelihood, posterior
	
	# Prior
	prior_id = "prior"
	prior_XML = xmlNode(name="distribution", attrs=list(id=prior_id, spec="util.CompoundDistribution"), .children=xml$priors )

	# Likelihood
	like_id = "likelihood"
	like_XML = xmlNode(name="distribution", attrs=list(id=like_id, spec="util.CompoundDistribution"), .children=xml$likes )
	
	# Posterior
	posterior_id = "posterior"
	posterior_idref = paste0("@", posterior_id)
	posterior_XML = xmlNode(name="distribution", attrs=list(id=posterior_id, spec="util.CompoundDistribution"), .children=list(prior_XML, like_XML) )
	
	# Logs of prior, likelihood, posterior
	logs_XML = list(
	xmlNode(name="log", attrs=list(idref=posterior_id) ),
	xmlNode(name="log", attrs=list(idref=prior_id) ),
	xmlNode(name="log", attrs=list(idref=like_id) )
	)
	
	# Add to xml
	xml$tracelog = c(logs_XML, xml$tracelog)
	xml$screenlog = c(logs_XML, xml$screenlog)
	
	# Add these to the screen log and trace log:
	#    <log id="ESS.0" spec="util.ESS" arg="@posterior"/>
    #    <log id="ESS.0" spec="util.ESS" arg="@likelihood"/>
    #    <log id="ESS.0" spec="util.ESS" arg="@prior"/>
	
	# Adding these to tracelog kills Tracer because the first row has
	# "N" as the output
	# xml$tracelog = cl(ESS_log1, ESS_log2, ESS_log3, xml$tracelog)

	ESS_log1 = xmlNode(name="log", attrs=list(id="ESS.po", spec="util.ESS", arg="@prior"))
	ESS_log2 = xmlNode(name="log", attrs=list(id="ESS.pr", spec="util.ESS", arg="@likelihood"))
	ESS_log3 = xmlNode(name="log", attrs=list(id="ESS.l", spec="util.ESS", arg="@posterior"))
	xml$screenlog = cl(ESS_log1, ESS_log2, ESS_log3, xml$screenlog)
	
	
	# Assemble the state
	storeEvery_statelog = sprintf("%0.0f", as.numeric(settings$state_store))
	state_id = "state"
	state_XML = xmlNode(name="state", attrs=list(id=state_id, storeEvery=storeEvery_statelog), .children=xml$state )
	
	# Assemble the trace log
	storeEvery_tracelog = sprintf("%0.0f", as.numeric(settings$tracelog_store))
	fileName = settings$tracelog_fn
	traceLog_id = "tracelog"
	traceLog_XML = xmlNode(name="logger", attrs=list(id=traceLog_id, fileName=fileName, logEvery=storeEvery_tracelog, model=posterior_idref, sanitiseHeaders="true", sort="smart"), .children=xml$tracelog )

	# Assemble the screen log
	storeEvery_screenlog = sprintf("%0.0f", as.numeric(settings$screenlog_store))
	screenLog_id = "screenlog"
	screenLog_XML = xmlNode(name="logger", attrs=list(id=screenLog_id, fileName="", logEvery=storeEvery_screenlog, sanitiseHeaders="true", sort="smart"), .children=xml$screenlog )

	# Assemble the tree log
	storeEvery_treelog = sprintf("%0.0f", as.numeric(settings$treelog_store))
	fileName = settings$treelog_fn
	treeLog_id = "treelog"
	# *NO* , substitutions="false" on plain treelog
	treeLog_XML = xmlNode(name="logger", attrs=list(id=treeLog_id, fileName=fileName, logEvery=storeEvery_treelog, mode="tree"), .children=xml$treelog )

	# Assemble the "subs" log (tree of substitutions rather than time)
	storeEvery_treelog = sprintf("%0.0f", as.numeric(settings$treelog_store))
	fileName = settings$subslog_fn
	if (isblank_TF(fileName) == FALSE)
		{
		# *NO* , substitutions="true" on plain treelog
		subsLog_id = "subslog"
		subsLog_XML = xmlNode(name="logger", attrs=list(id=subsLog_id, fileName=fileName, logEvery=storeEvery_treelog, mode="tree"), .children=xml$subslog )
		} else {
		subsLog_XML = xmlCommentNode(" Note: substitutions treelog (subslog.txt) not saved, as 'subslog_fn' in worksheet 'run' was left blank. ")
		} # END if (isblank_TF(fileName) == FALSE)
	
	
	# SpeciesTree/GeneTree loggers
	if (is.null(xml$speciesTreeLogger_XMLs) == FALSE)
		{
		speciesTreeLogger_XMLs = xml$speciesTreeLogger_XMLs
		} else {
		speciesTreeLogger_XMLs = NULL
		}
	if (is.null(xml$geneTree_treeLogs_XMLs) == FALSE)
		{
		geneTree_treeLogs_XMLs = xml$geneTree_treeLogs_XMLs
		} else {
		geneTree_treeLogs_XMLs = NULL
		}


	#######################################################
	# Remove any duplicates that were accidentally created
	# (excluding XMLCommentNodes, blanks, etc.)
	#######################################################
	xml = remove_duplicate_xmls(xml)


	# Assemble the run tag
	run_id = "mcmc"
	chainLength = settings$chainLength
	tmp_children = c( 
		list(bl()), 
		list(state_XML), 
		xml$starting_tree,	# Put the starting tree in here
		list(bl()), 
		xml$init, 			# Initializers
		list(bl()), 
		list(posterior_XML), 
		list(bl()), 
		list(bl()), 
		xml$operators, 
		list(bl()), 
		list(bl()), 
		list(traceLog_XML), 
		list(bl()), 
		list(bl()), 
		list(screenLog_XML), 
		list(bl()), 
		list(bl()), 
		list(treeLog_XML), 
		list(bl()), 
		list(bl()), 
		list(subsLog_XML),
		list(bl()), 
		list(bl()), 
		speciesTreeLogger_XMLs,
		list(bl()), 
		list(bl()), 
		geneTree_treeLogs_XMLs
		)
	
	# Assemble the big XML
	chainLength_txt = sprintf("%0.0f", as.numeric(chainLength))
	run_XML = xmlNode(name="run", attrs=list(id=run_id, chainLength=chainLength_txt, spec="MCMC"), .children=tmp_children )
	run_XMLs = list(bl(), bl(), xmlCommentNode(" Defining and logging the MCMC search chain "), run_XML)
	
	# Number of trees sampled
	numtrees_sampled = round(chainLength / settings$treelog_store)
	burnin = round(numtrees_sampled * burnin_fraction)
	
	#######################################################
	# LIST ALL XML TAGS
	#######################################################
	ALL_XML = c(
	xml$header,
	xml$taxa,
	xml$sequences,
	xml$partitions,
	xml$misc,
	xml$sitemodels,
	xml$clock,
	xml$tree,	
	run_XMLs
	)	
	
	# Make the <beast> node!
	beast_XML = xmlNode(name="beast", attrs=list(beautitemplate="Standard", beautistatus="", namespace=XML_namespaces, version="2.0"), .children=ALL_XML)
	
	
	# Write to output file?
	if (is.null(outfn))
		{
		return(beast_XML)
		} else {
		# Save the XML
		saveXML(doc=beast_XML, file=outfn, prefix='<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
		
		#cmdstr = paste('perl -e "s/\\(A\\ C\\)/M/g;" -pi $(find ', file, ' -type f)', sep="")
		str_to_change = "<!-- blank line -->"
		cmdstr1 = paste('perl -e "s/', str_to_change, '//g;" -pi.bak $(find ', outfn, ' -type f)', sep="")
		system(cmdstr1)
		
		str_to_change = "&quot;"
		cmdstr2 = paste('perl -e "s/', str_to_change, '/\\"/g;" -pi.bak $(find ', outfn, ' -type f)', sep="")
		system(cmdstr2)

		str_to_change = "&gt;"
		cmdstr3 = paste('perl -e "s/', str_to_change, '/>/g;" -pi.bak $(find ', outfn, ' -type f)', sep="")
		system(cmdstr3)

		str_to_change = "&lt;"
		cmdstr4 = paste('perl -e "s/', str_to_change, '/</g;" -pi.bak $(find ', outfn, ' -type f)', sep="")
		system(cmdstr4)

		str_to_change = "&amp;"
		cmdstr5 = paste('perl -e "s/', str_to_change, '/&/g;" -pi.bak $(find ', outfn, ' -type f)', sep="")
		system(cmdstr5)

		# change apostrophes
		str_to_change = "&apos;"
		cmdstr6 = paste('perl -e "s/', str_to_change, '/\\`/g;" -pi.bak $(find ', outfn, ' -type f)', sep="")
		system(cmdstr6)
		
		if (printall != "none")
			{
			cat("\n\nModifying '", outfn, "' by replacing 'blank line' with blank lines for formatting. Terminal command:\n ", cmdstr1, " \n\n", sep="")
			cat("\n\nModifying '", outfn, "' by replacing '&quot;' with double quotes. Terminal command:\n ", cmdstr2, " \n\n", sep="")
			cat("\n\nModifying '", outfn, "' by replacing '&gt;' with '>'. Terminal command:\n ", cmdstr3, " \n\n", sep="")
			cat("\n\nModifying '", outfn, "' by replacing '&lt;' with '<'. Terminal command:\n ", cmdstr4, " \n\n", sep="")
			cat("\n\nModifying '", outfn, "' by replacing '&amp;' with '&'. Terminal command:\n ", cmdstr5, " \n\n", sep="")
			cat("\n\nModifying '", outfn, "' by replacing '&apos;' with `. Terminal command:\n ", cmdstr6, " \n\n", sep="")
			}
		

		
		if (printall != "none")
			{
		cat("\n\n==================================================")
		cat("\nNote: BEAST analysis gets easier when you learn the command-line commands.\n\nHere are some likely commands, assuming you have Mac OS X 10.9, and are using Terminal (usually found in: /Applications/Utilities/Terminal.app). These may not be perfect -- it all depends on the locations of each file/program. And YMMV (your mileage may vary) on Windows etc.\n\n")

		cat("(1) cd (change directory) to wherever you saved outfn. Probably your working directory:\n")
		
		run_cmd1 = paste0(">>> cd ", addslash(getwd()))
		cat(run_cmd1)
		cat("\n\n")
		cat("(2) Run Beast2.\n")
		run_cmd2 = paste0(">>> java -Xms512m -Xmx512m -jar /Applications/BEAST_", BEAST2_version, "/lib/beast.jar -java -seed 754321 -overwrite ", outfn)
		cat(run_cmd2)
		cat("\n\n")
		cat("-java means you don't need Beagle, but makes it slower.\n-Xms512m and -Xmx512 control how much RAM is used; double etc., if you need more.\n-seed is the starting seed so you can replicate the analysis.\n-overwrite means you will overwrite the output file. Rename it if you don't want to do this!\n\n")

		cat("(3) This runs the help for Beast. See also google and the Beast2 website and beast-users google group!\n")

		run_cmd3 = paste0(">>> java -Xms512m -Xmx512m -jar /Applications/BEAST_", BEAST2_version, "/lib/beast.jar -help")
		cat(run_cmd3)
		cat("\n\n")

		cat("(4) This runs TreeAnnotator, which goes through your post-burnin tree samples and calculates an MCC (Maximum Clade Credibility) tree. Note that this is not the same thing as e.g. a majority-rule 'consensus' tree. See: http://en.wikipedia.org/wiki/Maximum_clade_credibility_tree\n")
		
		mcc_fn = paste0(get_fn_prefix(settings$treelog_fn), ".mcc")
		mcc_fn
		run_cmd4 = paste0(">>> java -Xms512m -Xmx512m -jar /Applications/BEAST_", BEAST2_version, "/treeannotator.jar -heights median -burnin ", burnin, " -limit 0 ", settings$treelog_fn, " ", mcc_fn)
		cat(run_cmd4)
		cat("\n\n")
		
		cat("(5) This runs the help for TreeAnnotator. See also google and the Beast2 website and beast-users google group!\n")
		run_cmd5 = paste0(">>> java -Xms512m -Xmx512m -jar /Applications/BEAST_", BEAST2_version, "/lib/treeannotator.jar -help")
		cat(run_cmd5)
		cat("\n\n")
		
		cat("(6) This opens the MCC tree in FigTree. You should also inspect the trace log files with Tracer.app!!\n")
		mcc_fn_wDir = slashslash(paste0(addslash(getwd()), mcc_fn))
		run_cmd6 = paste0(">>> open -n /Applications/beast/FigTree_v1.4.0.app --args ", mcc_fn)
		cat(run_cmd6)
		cat("\n==================================================\n")

		cat("\n==================================================\n")
		cat("All commands, without comments:\n")
		txt = paste(c(run_cmd1, run_cmd2, run_cmd3, run_cmd4, run_cmd5, run_cmd6), sep="\n", collapse="\n")
		txt = gsub(pattern=">>> ", replacement="", x=txt)
		cat(txt)
		cat("\n==================================================\n")

	
		txt = paste("\n==================================================\nBEASTmasteR is writing Beast2 XML to:\n\n", outfn, "\n\nYour current working directory (getwd()) is:\n\n", getwd(), "\n==================================================\n", "\n====================================================================================================\n====================================================================================================\nIF YOU USE BEASTmasteR IN YOUR PUBLICATIONS, PLEASE CITE IT: http://phylo.wikidot.com/beastmaster#citation\n====================================================================================================\n====================================================================================================\n\n")
		cat(txt)
		
		} # END if (printall != "none")
		

		
		return(outfn)
		} # END if (is.null(outfn))
	
	txt = "ERROR in parse_run(): Shouldn't get here."
	stop(txt)
	} # END parse_run