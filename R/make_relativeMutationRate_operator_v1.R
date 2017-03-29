
pick_a_partitionName_for_taxa_list <- function(seqs_df, morphList=NULL)
	{
	partitionNames = unique(seqs_df$filteredAlignmentName)
	rownums = match(x=partitionNames, table=seqs_df$filteredAlignmentName)
	dataTypes_for_uniq_partitionNames = seqs_df$type[rownums]
	geneTreeNames = seqs_df$geneTreeName[rownums]
	
	# Sequence partitions do NOT have dataType "morph"
	TF1 = dataTypes_for_uniq_partitionNames != "morph"
	TF2 = dataTypes_for_uniq_partitionNames != "continuous"
	# Also exclude partitions that are for gene trees, not species trees
	TF3 = isblank_TF(geneTreeNames)
	
	seqs_TF = (TF1 + TF2 + TF3) == 3
	partitionNames_seq = partitionNames[seqs_TF]

	if (is.null(morphList))
		{
		# Put sequences first, then morphology
		alignment_names_w_taxa = partitionNames_seq
		} else {
		# Put sequences first, then morphology
		partitionNames_morph = gsub(pattern="datatype_", replacement="", x=morphList$name_of_userDataType)

		alignment_names_w_taxa = c(partitionNames_seq, partitionNames_morph)
		} # END if (is.null(morphLengths))
	
	if (length(alignment_names_w_taxa) < 1)
		{
		txt = "\n\nERROR in pick_a_partitionName_for_taxa_list(): length(alignment_names_w_taxa) < 1, meaning neither seqs_df nor morphList provided a partition name. If this happened because all of your sequences are for gene trees ('data' worksheet, 'geneTreeName' column), i.e. you are running starBeast2, then use pick_a_partitionName_for_taxa_list_starBEAST().\n\n"
		cat(txt)
		warning(txt)
		return(NULL)
		} # END if (length(alignment_names_w_taxa) < 1)
	
	# Otherwise, pick the first usable alignment
	alignment_name_w_taxa = alignment_names_w_taxa[1]
	return(alignment_name_w_taxa)
	} # END pick_a_partitionName_for_taxa_list




# Call this, by looping through each clockModel_name
make_relativeMutationRate_operator <- function(seqs_df, partitionLengths_df=NULL, morphLengths=NULL, clockModel_name=NULL, relrate_suffix = "_relRate", force_relRates=FALSE, xml=NULL)
	{	
	# Remove use="no" rows...
	seqs_df = seqs_df[seqs_df$use != "no", ]
	
	# Error check
# 	if (is.null(partitionLengths_df) == FALSE)
# 		{
# 		if (nrow(partitionLengths_df) != nrow(seqs_df))
# 			{
# 			txt = paste0("ERROR in make_relativeMutationRate_operator(): nrow(partitionLengths_df)='", nrow(partitionLengths_df), " != nrow(seqs_df)='", nrow(seqs_df), "'. They must be equal (after seqs_df$use are removed).")
# 			cat("\n\n")
# 			cat(txt)
# 			cat("\n\n")
# 			} # END if (nrow(partitionLengths_df) != nrow(seqs_df))
# 		} # END if (is.null(partitionLengths_df) == FALSE)
	
	# If multiple clocks are specified
	if (is.null(clockModel_name) == FALSE)
		{
		# Remove those with the wrong clock
		TF = seqs_df$clockModel_name == clockModel_name
		tmp_seqs_df = seqs_df[TF, ]
		FixMeanMutationRatesOperator_id = paste0("FixMeanMutationRatesOperator_for_partitions_within_", clockModel_name)
		} else {
		FixMeanMutationRatesOperator_id = "FixMeanMutationRatesOperator"
		tmp_seqs_df = seqs_df
		}# END if (is.null(clockModel_name) == FALSE)
	
	# Get the list of relRates
	relRates_list = unique(tmp_seqs_df$clockModel_relRates)
	relRates_list = relRates_list[isblank_TF(relRates_list) == FALSE]
	

	
	# If there are no relative rates
	if ( (length(relRates_list) == 0) || (length(relRates_list) == 1) )
		{
		if (force_relRates == TRUE)
			{
			relRates_list = paste0(tmp_seqs_df$clockModel_name, "_oneRate")
			} else {
			if (is.null(xml) == TRUE)
				{
				res = NULL
				return(res)
				} else {
				return(xml)
				}
			}
		} # END if (length(relRates_list) == 0)
	
	rownums = match(x=relRates_list, table=tmp_seqs_df$clockModel_relRates)
	dataTypes = tmp_seqs_df$type[rownums]
	TF1 = dataTypes != "morph"
	TF2 = dataTypes != "continuous"
	seqs_TF = (TF1 + TF2) == 2
	relRates_list_seqs = relRates_list[seqs_TF]


	# For each relRates_list item, get the partitionName
	partitionNames = NULL
	for (i in 1:length(relRates_list_seqs))
		{
		TF = tmp_seqs_df$clockModel_relRates == relRates_list_seqs[i]
		partitionNames = c(partitionNames, tmp_seqs_df$filteredAlignmentName[TF][1])
		}
		
	
	# Get sequence lengths within each relRate
	partitionNames_seq = relRates_list_seqs
	partitionLengths_seq = NULL
	if (length(relRates_list_seqs) < 1)
		{
		nada=NULL
		} else {
		for (r in 1:length(relRates_list_seqs))
			{
			TF = tmp_seqs_df$clockModel_relRates == relRates_list_seqs[r]
			tmp_seqs_df2 = tmp_seqs_df[TF,]
		
			# Get the data type for this relRate
			tmp_datatype = tmp_seqs_df2$type[1]
		
			# Subset the clock to this datatype
			TF2 = tmp_seqs_df$type == tmp_datatype
			tmp_seqs_df3 = tmp_seqs_df[TF2,]
		
			if (tmp_datatype == "morph")
				{
				TF3 = tmp_seqs_df3$clockModel_relRates == relRates_list_seqs[r]
				morph_rownum = (1:length(TF3))[TF3]
				tmp_length = morphLengths[morph_rownum]
				} else {
				TF3 = tmp_seqs_df3$clockModel_relRates == relRates_list_seqs[r]
				DNA_partitionName = tmp_seqs_df3$filteredAlignmentName[TF3][1]
				tmp_length = partitionLengths_df$partitionLengths[partitionLengths_df$partitionNames==DNA_partitionName]
				} # END if (tmp_datatype == "morph")

	# 		BROKEN: tmp_partitionLengths = partitionLengths_df$partitionLengths[TF]
	# 		startvals = tmp_seqs_df$startchar[TF]
	# 		stopvals = tmp_seqs_df$stopchar[TF]
	# 		byvals = tmp_seqs_df$by[TF]
	# 		listvals = tmp_seqs_df$list[TF]
	# 		
	# 		tmp_length = 0
	# 		for (q in 1:sum(TF))
	# 			{
	# 			if (isblank_TF(listvals[q]) == FALSE)
	# 				{
	# 				char_positions = trim(strsplit(x=listvals[q], split=",")[[1]])
	# 				} else {
	# 				char_positions = seq(startvals[q], stopvals[q], by=byvals[q])
	# 				}
	# 			}
	# 		
	# 		tmp_partitionLengths = ceiling(tmp_partitionLengths / byvals)
			# Add up all of the partitionLengths in this particular relRate category
			partitionLengths_seq = c(partitionLengths_seq, tmp_length)
			}
		} # END if (length(relRates_list_seqs) < 1)
	
	#partitionLengths_seq = partitionLengths_df$partitionLengths
	
	# Get the partitionNames for morphology
	morph_TF = dataTypes == "morph"
	partitionNames_morph = relRates_list[morph_TF]
	if (length(partitionNames_morph) == 0)
		{
		morphLengths = partitionNames_morph
		}
	
	
	# Check length of partitionLengths, morphLengths
	if (length(partitionLengths) != length(partitionNames_seq))
		{
		errortxt = paste0("\n\nSTOP ERROR in make_relativeMutationRate_operator():\nLength of 'partitionLengths' is ", length(partitionLengths), " but the number of non-morphology\npartitions in 'tmp_seqs_df', partitionNames_seq, is ", length(partitionNames_seq), ".\n\n")
		cat(errortxt)
		cat("PartitionNames extracted for clockModel_name '", clockModel_name, "', from tmp_seqs_df$filteredAlignmentName:\n")
		print(partitionNames_seq)
		cat("\n\nInput 'partitionLengths_df':\n")
		print(partitionLengths_df)
		stop(errortxt)
		}

	if ( (length(partitionNames_morph) > 0) && (length(morphLengths) != length(partitionNames_morph)) )
		{
		#######################################################
		# Special case: ordered/unordered data from 
		# the same file
		# Try to fix by re-doing morphLengths
		#######################################################
		morphLengths = NULL
		for (ii in 1:nrow(tmp_seqs_df))
			{
			if (isblank_TF(tmp_seqs_df$list[ii]) == FALSE)
				{
				tmp_length = length(strsplit(tmp_seqs_df$list[ii], split=",")[[1]])
				} else {
				if (isblank_TF(tmp_seqs_df$endchar[ii]) == TRUE)
					{
					tmp_items = seq(tmp_seqs_df$startchar[ii], tmp_seqs_df$endchar[ii], tmp_seqs_df$by[ii])
					tmp_length = length(tmp_items)
					} else {
					warning(paste0("WARNING in make_relativeMutationRate_operator(): you should specify column 'endchar' or 'list' for a morphology character. As a random guess, BEASTmasteR assumes that the number of characters in morphology row #", ii, " is 100."))
					tmp_length = 100
					}
				}
			morphLengths = c(morphLengths, tmp_length)
			} # END for (ii in 1:nrow(tmp_seqs_df))
		} # END if


	if ( (length(partitionNames_morph) > 0) && (length(morphLengths) != length(partitionNames_morph)) )
		{
		errortxt = paste0("\n\nSTOP ERROR in make_relativeMutationRate_operator():\nLength of 'morphLengths' is ", length(morphLengths), " but the number of morphology\npartitions in 'tmp_seqs_df', partitionNames_morph, is ", length(partitionNames_morph), ".\n\n")
		cat(errortxt)
		cat("Unique partitionNames_morph extracted from tmp_seqs_df$filteredAlignmentName:\n")
		print(partitionNames_morph)
		cat("\n\nInput 'morphLengths':\n")
		print(morphLengths)
		stop(errortxt)
		}
	
	# Put sequences first, then morphology
	relative_mutationRate_ids = c(partitionNames_seq, partitionNames_morph)
	if (force_relRates == TRUE)
		{
		relative_mutationRate_ids = paste0(relative_mutationRate_ids, relrate_suffix)
		}
	all_partition_lengths = c(partitionLengths_seq, morphLengths)
	
	# Make the operators and logs for relative_mutationRate_ids
	screenLog_mutationRates_XML = list(bl(), xmlCommentNode(" Log the relative mutation rates (really, substitution rates!) ") )
	operator_mutationRates_XML = list(bl(), xmlCommentNode(" Operator on the relative mutation rates (really, substitution rates!) ") )
	
	
	if (length(relative_mutationRate_ids) < 1)
		{
		txt = "\n\nERROR in make_relativeMutationRate_operator(): (length(partitionNames_seq) + length(partitionNames_morph)) is LESS THAN ONE -- i.e., you have NO PARTITIONS to put rates on!\n\n"
		cat(txt)
		stop(txt)
		}

	if (length(relative_mutationRate_ids) == 1)
		{
		mutationRate_id = relative_mutationRate_ids[1]

		# Just one partition rate
		operator_mutationRate_id = paste0(mutationRate_id, "_Scaler")
		mutationRate_idref = paste0("@", mutationRate_id)
		operator_mutationRate_XML = xmlNode(name="operator", attrs=list(id=operator_mutationRate_id, parameter=mutationRate_idref, scaleFactor="0.5", weight="0.0", spec="ScaleOperator") )
		
		operator_mutationRate_XMLs = list(bl(), xmlCommentNode(" Only one partition input, so the DeltaExchangeOperator and FixMeanMutationRatesOperator are not needed, just a simple Scaler "), xmlCommentNode(" Weight is 0.0 since otherwise it's not identifiable with the clock rate "), operator_mutationRate_XML)
		
		relativeMutationRates_operators_XML = operator_mutationRate_XMLs

		# Log the mutation rate
		screenLog_mutationRate_XML = xmlNode(name="parameter", attrs=list(idref=mutationRate_id, name="log"))
		screenLog_mutationRates_XML = c(screenLog_mutationRates_XML, list(screenLog_mutationRate_XML))

		} else {
		# Multiple partitions with relative rates
		# Multiple partition rates, need to be scaled relative to each other
	
		# Loop through the partitions 
		for (i in 1:length(relative_mutationRate_ids))
			{
			mutationRate_id = relative_mutationRate_ids[i]

			screenLog_mutationRate_XML = xmlNode(name="parameter", attrs=list(idref=mutationRate_id, name="log"))
			screenLog_mutationRates_XML = c(screenLog_mutationRates_XML, list(screenLog_mutationRate_XML))

			operator_mutationRate_XML = xmlNode(name="parameter", attrs=list(idref=mutationRate_id))
			operator_mutationRates_XML = c(operator_mutationRates_XML, list(operator_mutationRate_XML))
			} # END for (i in 1:length(relative_mutationRate_ids))
		
		# Build the operator
		weightvector_list = all_partition_lengths
		weightvector_txt = paste(weightvector_list, sep="", collapse=" ")
		weightvector_XML = xmlNode(name="weightvector", weightvector_txt, attrs=list(id="weightparameter", dimension=length(relative_mutationRate_ids), estimate="false", lower="0", upper="0", spec="parameter.IntegerParameter") )
		weightvector_XML = c(list(bl()), list(xmlCommentNode(" The weights should be the number of characters in each partition ")), list(weightvector_XML))
	
		children_XMLs = c(operator_mutationRates_XML, weightvector_XML)
		relativeMutationRates_operators_XML = xmlNode(name="operator", attrs=list(id=FixMeanMutationRatesOperator_id, delta="0.75", weight="2.0", spec="DeltaExchangeOperator"), .children=children_XMLs)
		relativeMutationRates_operators_XML = c(list(bl()), list(xmlCommentNode(" Operators to shift the relative rates of the partitions about the overall (clock) mean ")), list(relativeMutationRates_operators_XML))
		} # END if (length(relative_mutationRate_ids) == 1)

	
	# Output the results
	if (is.null(xml) == TRUE)
		{
		res = NULL
		# Output to res
		res$relativeMutationRates_operators_XML = relativeMutationRates_operators_XML
		res$traceLog_mutationRates_XML = screenLog_mutationRates_XML
		res$screenLog_mutationRates_XML = screenLog_mutationRates_XML

		extract='
		relativeMutationRates_operators_XML = res$relativeMutationRates_operators_XML
		traceLog_mutationRates_XML = res$traceLog_mutationRates_XML 
		screenLog_mutationRates_XML = res$screenLog_mutationRates_XML
		'
		return(res)
		} else {
		# Update the 'xml' list
		xml$operators = c(xml$operators, relativeMutationRates_operators_XML)
		xml$tracelog = c(xml$tracelog, screenLog_mutationRates_XML)
		xml$screenlog = c(xml$screenlog, screenLog_mutationRates_XML)
		return(xml)
		}

	stop("ERROR in make_relativeMutationRate_operator(): Shouldn't get here.")
	} # END make_relativeMutationRate_operator <- function(seqs_df, partitionLengths=NULL, morphLengths=NULL, xml=NULL)