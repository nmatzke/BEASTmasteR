
pick_a_partitionName_for_taxa_list <- function(seqs_df, morphList=NULL)
	{
	partitionNames = unique(seqs_df$filteredAlignmentName)
	rownums = match(x=partitionNames, table=seqs_df$filteredAlignmentName)
	dataTypes_for_uniq_partitionNames = seqs_df$type[rownums]

	# Sequence partitions do NOT have dataType "morph"
	TF1 = dataTypes_for_uniq_partitionNames != "morph"
	TF2 = dataTypes_for_uniq_partitionNames != "continuous"
	seqs_TF = (TF1 + TF2) == 2
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
		txt = "\n\nERROR in pick_a_partitionName_for_taxa_list(): length(alignment_names_w_taxa) < 1, meaning neither seqs_df nor morphList provided a partition name\n\n"
		cat(txt)
		stop(txt)
		} # END if (length(alignment_names_w_taxa) < 1)
	
	# Otherwise, pick the first usable alignment
	alignment_name_w_taxa = alignment_names_w_taxa[1]
	return(alignment_name_w_taxa)
	} # END pick_a_partitionName_for_taxa_list

make_relativeMutationRate_operator <- function(seqs_df, partitionLengths=NULL, morphLengths=NULL, clockModel_name=NULL, xml=NULL)
	{
	# Remove use="no" rows...
	seqs_df = seqs_df[seqs_df$use == "yes", ]
	
	# If multiple clocks are specified
	if (is.null(clockModel_name) == FALSE)
		{
		# Remove those with the wrong clock
		TF = seqs_df$clockmodel_name == clockModel_name
		seqs_df = seqs_df[TF, ]
		FixMeanMutationRatesOperator_id = paste0("FixMeanMutationRatesOperator_for_partitions_within_", clockModel_name)
		} else {
		FixMeanMutationRatesOperator_id = "FixMeanMutationRatesOperator"
		}# END if (is.null(clockModel_name) == FALSE)
	
	partitionNames = unique(seqs_df$filteredAlignmentName)
	rownums = match(x=partitionNames, table=seqs_df$filteredAlignmentName)
	dataTypes = seqs_df$type[rownums]
	
	# Sequence partitions do NOT have dataType "morph" *or* "continuous
	TF1 = dataTypes != "morph"
	TF2 = dataTypes != "continuous"
	seqs_TF = (TF1 + TF2) == 2
	partitionNames_seq = partitionNames[seqs_TF]
	
	morph_TF = dataTypes == "morph"
	partitionNames_morph = partitionNames[morph_TF]
	
	# Check length of partitionLengths, morphLengths
	if (length(partitionLengths) != length(partitionNames_seq))
		{
		errortxt = paste0("\n\nSTOP ERROR in make_relativeMutationRate_operator():\nLength of 'partitionLengths' is ", length(partitionLengths), " but the number of non-morphology\npartitions in 'seqs_df' is ", length(partitionNames_seq), ".\n\n")
		cat(errortxt)
		cat("Unique partitionNames extracted from seqs_df$filteredAlignmentName:\n")
		print(partitionNames)
		cat("\n\nInput 'partitionLengths':\n")
		print(partitionNames)
		stop(errortxt)
		}
	if (length(morphLengths) != length(partitionNames_morph))
		{
		errortxt = paste0("\n\nSTOP ERROR in make_relativeMutationRate_operator():\nLength of 'morphLengths' is ", length(morphLengths), " but the number of morphology\npartitions in 'seqs_df' is ", length(partitionNames_morph), ".\n\n")
		cat(errortxt)
		cat("Unique partitionNames extracted from seqs_df$filteredAlignmentName:\n")
		print(partitionNames)
		cat("\n\nInput 'morphLengths':\n")
		print(partitionNames)
		stop(errortxt)
		}
	
	# Put sequences first, then morphology
	relative_mutationRate_ids = c(partitionNames_seq, partitionNames_morph)
	relative_mutationRate_ids = paste0(relative_mutationRate_ids, "_relRate")
	all_partition_lengths = c(partitionLengths, morphLengths)
	
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