
get_DNA_from_sequence_tag <- function(xmltag)
	{
	dnaval = xmlAttrs(xmltag)$value
	dnalist = strsplit(dnaval, split="")[[1]]
	return(dnalist)
	} # END get_DNA_from_sequence_tag <- function(xmltag)

get_name_from_sequence_tag <- function(xmltag, remove_taxon=TRUE, string_to_delete=string_to_delete)
	{
	nameval = xmlAttrs(xmltag)$id
	taxonval = xmlAttrs(xmltag)$taxon
	
	if (remove_taxon == TRUE)
		{
		nameval = gsub(pattern=taxonval, replacement="", x=nameval)
		}
	if (isblank_TF(string_to_delete) == FALSE)
		{
		nameval = gsub(pattern=string_to_delete, replacement="", x=nameval)
		}
	return(nameval)
	} # END get_DNA_from_sequence_tag <- function(xmltag)

get_taxon_from_sequence_tag <- function(xmltag)
	{
	taxonval = xmlAttrs(xmltag)$taxon
	return(taxonval)
	} # END get_DNA_from_sequence_tag <- function(xmltag)




convert_data_XML_to_dnalists <- function(data_XML, remove_taxon=TRUE, string_to_delete="_seqfrom_", name_seqs_by_taxonname=TRUE, use_single_unique_seqName_for_alignment_name=TRUE)
 	{

	defaults = '
	data_XML = data_XMLs$data_XML
	remove_taxon=TRUE
	string_to_delete = "_seqfrom_"
	name_seqs_by_taxonname=TRUE
	use_single_unique_seqName_for_alignment_name=TRUE
	'
	
	##########################################################################
	######### START EXAMPLE 
	##########################################################################
	example='
	runslow = TRUE
	data_XML_fn = "data_XML.Rdata"
	if (runslow)
		{
		data_XMLs = parse_datasets(seqs_df, xml=NULL, add_morphLength=TRUE, OTUs=OTUs, add_morphList=TRUE, printall="short", xlsfn=xlsfn, return_charsdf=TRUE, convert_ambiguous_to_IUPAC=FALSE)
		save(data_XMLs, file=data_XML_fn)
		} else {
		# Loads to "data_XML"
		load(file=data_XML_fn)
		} # END if (runslow)

	# For StarBeast2, prune out sequences from specimenNames not listed in taxonsets_df
	# (no action taken when not a StarBeast2 analysis)
	data_XML = prune_seqs_based_on_taxonsets(data_XML=data_XMLs$data_XML, taxonsets_df=taxonsets_df, StarBeast2_TF=StarBeast2_TF)
	data_XMLs$data_XML = data_XML
	
	# Extract DNA, use seqNames for alignment name, use taxon (specimen) names for sequence names
	#source("/drives/GDrive/__github/BEASTmasteR/R/extract_from_XML_v1.R")
	string_to_delete = "_seqfrom_"
	output = convert_data_XML_to_dnalists(data_XML=data_XMLs$data_XML, remove_taxon=TRUE, string_to_delete=string_to_delete, name_seqs_by_taxonname=TRUE, use_single_unique_seqName_for_alignment_name=TRUE)
	names(output)
	output$dna_seqNames_list
	output$dna_taxa_list
	output$dna_dfs_list

	# Extract DNA, use numbers for alignment names, use original sequence names
	string_to_delete = ""
	output = convert_data_XML_to_dnalists(data_XML=data_XMLs$data_XML, remove_taxon=FALSE, string_to_delete=string_to_delete, name_seqs_by_taxonname=FALSE, use_single_unique_seqName_for_alignment_name=FALSE)
	names(output)
	output$dna_seqNames_list
	output$dna_taxa_list
	output$dna_dfs_list
	'
	##########################################################################
	######### END EXAMPLE
	##########################################################################
	
	
	require(ape)
	require(BioGeoBEARS)
	require(gdata)	# for trim
	require(seqinr)
	require(stringr)
	require(pegas)


	numXMLs = length(data_XML)
	dna_dfs_list = NULL
	dna_seqNames_list = NULL
	dna_taxa_list = NULL
	alignment_names = NULL
	datatypes = NULL
	dnum = 0
	dna_alignment_lengths = NULL
	
	for (i in 1:numXMLs)
		{
		# Escape, if not data
		if (xmlName(data_XML[[i]]) != "data")
			{
			next()
			} # END if (xmlName(data_XML[[i]]) != "data")
		
		# Increment dnum
		dnum = dnum + 1

		# Datatypes
		datatype = xmlAttrs(data_XMLs$data_XML[[i]])$dataType
		datatypes = c(datatypes, datatype)
		
		# If it is data, let's copy it out to a DNAbin
		children = xmlChildren(data_XML[[i]])
	
		# Find the sequences
		seqs_TF = names(children) == "sequence"
		children2 = children[seqs_TF]
	
		namevals = sapply(X=children2, FUN=get_name_from_sequence_tag, remove_taxon=remove_taxon, string_to_delete=string_to_delete)
		names(namevals) = NULL
		dna_seqNames_list[[dnum]] = namevals
		
		# Extract the DNA	
		dnavals = lapply(X=children2, FUN=get_DNA_from_sequence_tag)

		# Get the taxon values
		if (name_seqs_by_taxonname == TRUE)
			{
			taxonvals = lapply(X=children2, FUN=get_taxon_from_sequence_tag)
			names(taxonvals) = namevals
			
			dna_taxa_list[[dnum]] = taxonvals
			} # END if (name_seqs_by_taxonname == TRUE)
	
		# If you want to have one unique name (e.g., the locus name) to name the alignment
		if (use_single_unique_seqName_for_alignment_name == TRUE)
			{
			unique_nameval = unique(namevals)
	
			if (length(unique_nameval) > 1)
				{
				strs = paste(unique_nameval, collapse=", ", sep="")
				strs = paste0("'", strs, "'")
			
				txt = paste0("STOP ERROR in convert_data_XML_to_dnadf():\n\nIn XML #", i, "/", numXMLs, ", there was more than one unique sequence name, given options remove_taxon=", remove_taxon, " and string_to_delete='", string_to_delete, "'. They were:\n\n", strs, ". Fix the inputs, or change use_single_unique_seqName_for_alignment_name to FALSE.")
			
				cat("\n\n")
				cat(txt)
				cat("\n\n")
				stop(txt)
				} # END if (length(unique_nameval) == FALSE)
		
			# Use the taxon names, not the sequence names
			if (name_seqs_by_taxonname == TRUE)
				{
				names(dnavals) = taxonvals
				} else {
				names(dnavals) = namevals
				} # END if (name_seqs_by_taxonname == TRUE)
		
			alignment_names = c(alignment_names, unique_nameval)
			cmdstr = paste0("dna_dfs_list$", unique_nameval, " = dnavals")
			eval(parse(text=cmdstr))
			
			tmp_alignment_length = length(dnavals[[1]])
			dna_alignment_lengths = c(dna_alignment_lengths, tmp_alignment_length)
			} else {
			# If you want to use the individual sequence names, not just the specimen names
			names(dnavals) = namevals
			dna_dfs_list[[dnum]] = dnavals
			} # END if (unique_nameval == TRUE)
	
		} # END for (i in 1:numXMLs)

	length(dna_dfs_list)

	output = NULL
	output$dna_seqNames_list = dna_seqNames_list
	if (name_seqs_by_taxonname == TRUE)
		{
		output$dna_taxa_list = dna_taxa_list
		}
	output$dna_dfs_list = dna_dfs_list
	
	# Unlist the results
	output$dna_seqNames_list = lapply(X=output$dna_seqNames_list, FUN=unlist)
	output$dna_taxa_list = lapply(X=output$dna_taxa_list, FUN=unlist)
	output$ntaxa = unlist(lapply(X=dna_dfs_list, FUN=length))
	output$dna_alignment_lengths = unlist(dna_alignment_lengths)
	output$datatypes = datatypes

	
	extract='
	dna_seqNames_list = output$dna_seqNames_list
	dna_taxa_list = output$dna_taxa_list
	dna_dfs_list = output$dna_dfs_list
	ntaxa = output$ntaxa
	dna_alignment_lengths = output$dna_alignment_lengths
	datatypes = output$datatypes
	'
	
	return(output)
	} # END convert_data_XML_to_dnalists <- function(data_XML, remove_taxon=TRUE, string_to_delete="_seqfrom_", name_seqs_by_taxonname=TRUE, use_single_unique_seqName_for_alignment_name=TRUE)



















charseqs_list_to_partition_txt <- function(dna_dfs_list, returnval="text")
	{
	pos = 0
	partition_nums_txt = ""
	num_alignments = length(dna_dfs_list)
	alignment_names = names(dna_dfs_list)
	alignment_lengths = NULL
	startpositions = NULL
	endpositions = NULL
	super_seqs_list = NULL
	
	for (i in 1:num_alignments)
		{
		tmp_alignment_charseqs = dna_dfs_list[[i]]
		seqlengths = sapply(X=tmp_alignment_charseqs, FUN=length)
		alignment_length = unique(seqlengths)
		if (length(alignment_length) > 1)
			{
			alignment_lengths_txt = paste(alignment_lengths, collapse=", ", sep="")
			txt = paste0("STOP ERROR in charseqs_list_to_partition(): not all sequences in alignment #", i, "/", num_alignments, " (alignment name: '", alignment_names[i], "') are equal. The unique alignment_lengths are printed below: ", alignment_lengths_txt)
			
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			paste(alignment_lengths)
			cat("\n\n")
			stop(txt)
			} # END if (length(alignment_length) > 1)
		
		# OK, if we have one alignment 
		startpos = pos+1
		endpos = pos+alignment_length
		pos = endpos
		
		startpositions = c(startpositions, startpos)
		endpositions = c(endpositions, endpos)
		alignment_lengths = c(alignment_lengths, alignment_length)
		newtxt = paste0(startpos, "-", endpos)
		if (partition_nums_txt != "")
			{
			partition_nums_txt = paste(partition_nums_txt, newtxt, sep=", ")
			} else {
			partition_nums_txt = newtxt
			}

		# Add the sequences together
		if (i == 1)
			{
			super_seqs_list = tmp_alignment_charseqs
			} else {
			# Match names
			match_newseqs_to_superseqs = match(x=names(tmp_alignment_charseqs), table=names(super_seqs_list))
			newTF = is.na(match_newseqs_to_superseqs)
			nums_of_seqnames_that_are_new = (1:length(tmp_alignment_charseqs))[newTF==TRUE]
			nums_of_seqnames_that_are_old = (1:length(tmp_alignment_charseqs))[newTF==FALSE]
			place_to_put_seqnames_that_are_old = match_newseqs_to_superseqs[newTF==FALSE]
			nums_of_superseqnames_not_matched_TF = ( ((1:length(super_seqs_list)) %in% match_newseqs_to_superseqs) == FALSE )
			nums_of_superseqnames_not_matched = (1:length(super_seqs_list))[nums_of_superseqnames_not_matched_TF]
			
			if (length(nums_of_seqnames_that_are_old) > 0)
				{
				for (j in 1:length(nums_of_seqnames_that_are_old))
					{
					seqnum_to_take_from = nums_of_seqnames_that_are_old[j]
					seqnum_to_splice_at = place_to_put_seqnames_that_are_old[j]
					super_seqs_list[[seqnum_to_splice_at]] = c(super_seqs_list[[seqnum_to_splice_at]], tmp_alignment_charseqs[[seqnum_to_take_from]])
					} # END for (j in 1:length(nums_of_seqnames_that_are_old))
				} # END if (length(nums_of_seqnames_that_are_old) > 0)

			# For the new sequence names
			listnum = length(super_seqs_list)
			if (length(nums_of_seqnames_that_are_new) > 0)
				{
				for (k in 1:length(nums_of_seqnames_that_are_new))
					{
					listnum = listnum+1
					seqnum_to_take_from = nums_of_seqnames_that_are_new[k]
					seqname_to_add = names(tmp_alignment_charseqs)[seqnum_to_take_from]
					# Add the new sequence
					cmdstr = paste0("super_seqs_list$", seqname_to_add, " = tmp_alignment_charseqs[[seqnum_to_take_from]]")
					eval(parse(text=cmdstr))
					} # END for (k in 1:length(nums_of_seqnames_that_are_new))
				} # END if (length(nums_of_seqnames_that_are_new) > 0)
			
			if (length(nums_of_superseqnames_not_matched) > 0)
				{
				for (l in 1:length(nums_of_superseqnames_not_matched))
					{
					print("Adding Ns")
					newseq = rep("N", times=alignment_length)
					super_seqs_list[[nums_of_superseqnames_not_matched[l]]] = c(super_seqs_list[[nums_of_superseqnames_not_matched[l]]], newseq)
					}
				} # END if (length(nums_of_superseqnames_not_matched) > 0)
				
			} # END if (i == 1)


		} # END for (i in 1:num_alignments)
	
	partitions_df = as.data.frame(cbind(startpositions, endpositions, alignment_lengths), stringsAsFactors=FALSE)
	names(partitions_df) = c("start", "end", "length")
	
	if (returnval == "super_seqs_list")
		{
		return(super_seqs_list)
		} 
	
	if (returnval == "text")
		{
		return(partition_nums_txt)
		} 
	if (returnval == "df")
		{
		return(partitions_df)
		}
	} # END charseqs_list_to_partition


# *DO* count "N"/"n" (IUPAC ambiguous)
count_ambig_DNA <- function(charslist)
	{
	TF1 = charslist == "?"
	TF2 = charslist == "-"
	TF3 = charslist == "N"
	TF4 = charslist == "n"
	TF = TF1 + TF2 + TF3 + TF4
	num_ambig = sum(TF)
	return(num_ambig)
	}

# Don't count "N"/"n" (Asparagine)
count_ambig_nonDNA <- function(charslist)
	{
	TF1 = charslist == "?"
	TF2 = charslist == "-"
	TF = TF1 + TF2
	num_ambig = sum(TF)
	return(num_ambig)
	}

# *DO* count "N"/"n" (IUPAC ambiguous)
count_ambig_DNA_in_string <- function(seq_string)
	{
	count1 = str_count(string=seq_string, pattern="\\?")
	count2 = str_count(string=seq_string, pattern="-")
	count3 = str_count(string=seq_string, pattern="N")
	count4 = str_count(string=seq_string, pattern="n")
	num_ambig = count1 + count2 + count3 + count4
	return(num_ambig)
	}

# Don't count "N"/"n" (Asparagine)
count_ambig_nonDNA_in_string <- function(seq_string)
	{
	count1 = str_count(string=seq_string, pattern="\\?")
	count2 = str_count(string=seq_string, pattern="-")
	count3 = str_count(string=seq_string, pattern="Z")
	count4 = str_count(string=seq_string, pattern="z")
	num_ambig = count1 + count2 + count3 + count4
	return(num_ambig)
	}










write_DNAdf_to_files <- function(dna_taxa_list, dna_dfs_list, datatypes, prefix="", suffix="", write_unmerged_also_TF=TRUE, outdir="z_saved")
	{
	defaults='
	suffix="SavedOut"
	write_unmerged_also_TF=TRUE
	outdir="z_saved"
	'
	# The list of files you've created
	outfns = NULL
	
	# Create the directory if it doesn't exist (and suppress the warning, if it does)
	outdir1 = paste0(outdir, "1")
	outdir2 = paste0(outdir, "2")
	dir.create(path=outdir1, showWarnings=FALSE)
	
	# End the outdir with a slash
	if (endsWith(x=outdir1, suffix="/") == FALSE)
		{
		outdir1 = paste0(outdir1, "/")
		} # END if (endsWith(x=outdir, suffix="/") == FALSE)	
	if (endsWith(x=outdir2, suffix="/") == FALSE)
		{
		outdir2 = paste0(outdir2, "/")
		} # END if (endsWith(x=outdir, suffix="/") == FALSE)	
	
	# Only works for a single datatype; with multiple datatypes, run as a subset
	uniq_datatypes = unique(datatypes)
	if (length(uniq_datatypes) > 1)
		{
		outfns = NULL
		# Split into multiple groups
		for (i in 1:length(uniq_datatypes))
			{
			datatype = uniq_datatypes[i]
			TF = datatypes == datatype
			nums = (1:length(TF))[TF]
			tmp_datatype = gsub(pattern=" ", replacement="_", x=datatype)
			tmp_suffix = paste0(suffix, "_", tmp_datatype)
			tmp_outfns = write_DNAdf_to_files(dna_taxa_list=dna_taxa_list[nums], dna_dfs_list=dna_dfs_list[nums], datatypes=datatypes[nums], prefix=prefix, suffix=tmp_suffix, write_unmerged_also_TF=write_unmerged_also_TF, outdir=outdir)
			outfns = c(outfns, tmp_outfns)
			} # END for (i in 1:length(uniq_datatypes))	
		return(outfns)
		} # END if (length(uniq_datatypes) == 1)
	
	# Datatype is unique, so...
	if (unique(datatypes) == "nucleotide")
		{
		datatype = "dna"
		} else if (unique(datatypes) == "protein") {
		datatype = "protein"
		} else if (unique(datatypes) == "aminoacid") {
		datatype = "protein"
		} else {
		datatype = "standard"
		}
	
		
	min_numseqs = 0
	
	if (prefix != "")
		{
		if (endsWith(prefix, suffix="_"))
			{
			prefix = prefix
			} else {
			prefix = paste0(prefix, "_")
			} # END if (endsWith(prefix, suffix="_"))
		} # END if (prefix != "")

	# Get list of OTUs (specimen names)
	locus_names = NULL
	list_of_OTUs = NULL
	for (i in 1:length(dna_taxa_list))
		{
		list_item = dna_taxa_list[[i]]
		locus_names = c(locus_names, names(list_item))
		list_of_OTUs = c(list_of_OTUs, unname(list_item))
		} # END for (i in 1:length(output$dna_taxa_list))
	loci_all = sort(unique(locus_names))
	OTUs_all = sort(unique(list_of_OTUs))
	length(OTUs_all)


	#######################################################
	# Get info on each locus
	#######################################################
	numchars_list = NULL
	for (i in 1:length(dna_dfs_list))
		{
		dna_df = dna_dfs_list[[i]]
		numchars_list = c(numchars_list, length(dna_df[[1]]))
		}

	# Merge merge merge
	out_seqs_charslist = list()
	for (i in 1:length(OTUs_all))
		{
		cmdstr = paste0("out_seqs_charslist$", OTUs_all[i], " = list()")
		eval(parse(text=cmdstr))
		#out_seqs_charslist[[i]] = list()
		}
	out_seqs_charslist

	#######################################################
	# Merge loci, save in outdir1
	#######################################################
	# Remove "__" (double underscore)
	tmptxt = paste0(prefix, "concatenated_", suffix, ".nex")
	tmptxt = gsub(pattern="__", replacement="_", x=tmptxt)
	out_nexfn = suppressWarnings(slashslash(paste0(outdir1, tmptxt)))

	tmptxt = paste0(prefix, "concatenated_", suffix, "_partitions.nex")
	tmptxt = gsub(pattern="__", replacement="_", x=tmptxt)
	partitions_nexfn = suppressWarnings(slashslash(paste0(outdir1, tmptxt)))
	
	tmptxt = paste0(prefix, "concatenated_", suffix, "_partitions_table.txt")
	tmptxt = gsub(pattern="__", replacement="_", x=tmptxt)
	partitions_table_fn = suppressWarnings(slashslash(paste0(outdir1, tmptxt)))

	tmptxt = paste0(prefix, "concatenated_", suffix, "+partitions.nex")
	tmptxt = gsub(pattern="__", replacement="_", x=tmptxt)
	both_nexfn = suppressWarnings(slashslash(paste0(outdir1, tmptxt)))

	tmptxt = paste0(prefix, "taxonsets_specimenNames_to_speciesNames_", suffix, ".txt")
	tmptxt = gsub(pattern="__", replacement="_", x=tmptxt)
	taxonsets_fn = suppressWarnings(slashslash(paste0(outdir1, tmptxt)))
	
	
	cat("\nConcatenating loci for output into file '", out_nexfn, "'...", sep="")
	cat("\nOf ", length(dna_dfs_list), " alignments, adding in #:\n", sep="")
	for (i in 1:length(dna_dfs_list))
		{
		cat(i, " ", sep="")
		names_to_add = names(dna_dfs_list[[i]])
		for (j in 1:length(out_seqs_charslist))
			{
			name_to_add_to = names(out_seqs_charslist)[j]
			pos_to_take_new_data = match(x=name_to_add_to, table=names_to_add)
	
			if (is.na(pos_to_take_new_data) == FALSE)
				{
				newdata = dna_dfs_list[[i]][[pos_to_take_new_data]]
				} else {
				newdata = rep(x="N", times=numchars_list[[i]])
				}
	
			# Stick on the new data
			out_seqs_charslist[[j]] = unlist(c(unlist(out_seqs_charslist[[j]]), newdata))
			}
		}
	cat("\n...done.\n")

	out_seqs_charslist
	lapply(X=out_seqs_charslist, FUN=length)


	dna_all_df = as.data.frame(out_seqs_charslist, stringsAsFactors=FALSE)

	#outfn = "dna_all.nex"
	write_nexus_data2(x=out_seqs_charslist, file=out_nexfn, format=datatype)


	# Write the partitioning scheme
	charset_strings = list()
	inum = 1
	charset_strings[[(inum=inum+1)]] = "begin mrbayes;"
	partitionNames = names(dna_dfs_list)
	partitionNames

	startnum = 1
	startnums = NULL
	endnums = NULL
	byval = 1
	for (i in 1:length(partitionNames))
		{
		endnum = startnum + numchars_list[i] - 1
		charset_strings[[(inum=inum+1)]] = paste0("charset ", partitionNames[i], " = ", startnum, "-", endnum, ";")
		startnums = c(startnums, startnum)
		endnums = c(endnums, endnum)
		startnum = endnum + 1
		} # END for (i in 1:length(partitionNames))

	partition_comma_txt = paste(partitionNames, collapse=", ", sep="")
	charset_strings[[(inum=inum+1)]] = paste0("partition userPartitions = ", length(partitionNames), ": ", partition_comma_txt, ";")
	charset_strings[[(inum=inum+1)]] = "END;"
	charset_strings
	
	# Write partitions as a table
	tmptable = cbind(partitionNames, startnums, endnums, byval)
	tmptable_df = as.data.frame(tmptable, stringsAsFactors=FALSE)
	tmptable_df = dfnums_to_numeric(tmptable_df)
	write.table(x=tmptable_df, file=partitions_table_fn, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
	
	# Write to a partitions-only file
	write.table(x=unlist(charset_strings), file=partitions_nexfn, append=FALSE, quote=FALSE, sep="\n", row.names=FALSE, col.names=FALSE)
	
	# Copy the concatentated NEXUS file
	file.copy(from=out_nexfn, to=both_nexfn)
	# Append the two files
	file.append(file1=both_nexfn, file2=partitions_nexfn)




	#######################################################
	# Write out the un-merged dataset to outdir2
	#######################################################
	numseqs_list = NULL
	specimenNames_kept = NULL
	cat("\nSaving individual loci filename(s): \n", sep="")

	if (write_unmerged_also_TF == TRUE)
		{
		#cat("\nOf ", length(dna_dfs_list), " alignments, processing/saving #:\n", sep="")
		dir.create(path=outdir2, showWarnings=FALSE)

		for (i in 1:length(dna_dfs_list))
			{
			tmpseqs = dna_dfs_list[[i]]
			names_to_add = names(tmpseqs)
			names_to_add
		
			numQs = sapply(X=tmpseqs, FUN=count_ambig_DNA)
			lengths = sapply(X=tmpseqs, FUN=length)
			fracts = 1-(numQs / lengths)

			#keepTF1 = fracts >= min_fraction_to_keep_sequence
			#keepTF2 = (names_to_add %in% specimenNames_to_keep) == TRUE
			keepTF1 = fracts >= 0
			keepTF2 = rep(TRUE, times=length(keepTF1))
		
			keepTF = (keepTF1 + keepTF2) == 2
			keepNums = (1:length(tmpseqs))[keepTF]
			tmpseqs2 = tmpseqs[keepNums]
		
			numseqs = sum(keepTF)
			if (numseqs >= min_numseqs)
				{
				names_to_add2 = names_to_add[keepNums]
				specimenNames_kept = c(specimenNames_kept, names_to_add2)
			
				# Write out to NEXUS
				tmptxt = paste0(prefix, get_fn_prefix(fn=partitionNames[i]), "_", suffix, "_locus.nexus")
				tmptxt = gsub(pattern="__", replacement="_", x=tmptxt)
				outfn = suppressWarnings(slashslash(paste0(outdir2, tmptxt)))
				outfns = c(outfns, outfn)

				cat(i, ":\t", outfn, "\n", sep="")
				write_nexus_data2(x=tmpseqs2, file=outfn, format=datatype)
				} # END if (numseqs >= min_numseqs)
			} # END for (i in 1:length(dna_dfs_list))
	
		specimenNames_kept_unique = unique(specimenNames_kept)	
		} # END if (merge_loci_TF == FALSE)
	specimenNames_kept_unique
	length(specimenNames_kept_unique)


	# Make the taxonsets file
	corresponding_species_list = NULL
	for (i in 1:length(specimenNames_kept_unique))
		{
		tmpname = specimenNames_kept_unique[i]
		matchnum = match(x=tmpname, table=taxonsets_df$specimenString)
		corresponding_species = taxonsets_df$speciesName[matchnum]
		corresponding_species_list = c(corresponding_species_list, corresponding_species)
		}
	taxonsets_output_df = cbind(specimenNames_kept_unique, corresponding_species_list)
	taxonsets_output_df = as.data.frame(taxonsets_output_df, stringsAsFactors=FALSE)
	if (length(names(taxonsets_output_df) ) == 2)
		{
		names(taxonsets_output_df) = c("traits", "species")
		} else if (length(names(taxonsets_output_df) ) == 1)
		{
		names(taxonsets_output_df) = c("species")
		} else {
		abcd = NULL
		}
		
	write.table(x=taxonsets_output_df, file=taxonsets_fn, quote=FALSE, sep="\t", append=FALSE, row.names=FALSE, col.names=TRUE)
	#moref(taxonsets_fn)
	outfns = c(outfns, taxonsets_fn)
	return(outfns)
	} # END write_DNAdf_to_files <- function(dna_taxa_list, dna_dfs_list, datatype="dna", prefix="", suffix="", write_unmerged_also_TF=TRUE)



get_DNA_numQs_from_data_XML <- function(xml_items_list)
	{
	defaults='
	xml_items_list = data_XMLs$data_XML
	'
	tdf = NULL
	
	# Get the number of question marks per taxon
	# in just the DNA alignment
	#xml_items_list = data_XML
	for (j in 1:length(xml_items_list))
		{
		xml_item = xml_items_list[[j]]
	
		tmp_datatype = xmlAttrs(xml_item)$dataType
		if ( (isblank_TF(tmp_datatype) == FALSE) && (tmp_datatype=="nucleotide") )
			{
			tmprows = NULL
			for (i in 1:length(xml_item))
				{
				xmltag = xml_item[[i]]
				taxon = xmlAttrs(xmltag)$taxon
				id = xmlAttrs(xmltag)$id
				dna = get_DNA_from_sequence_tag(xmltag)
				numQs = count_ambig_DNA(dna)
	
				tmprow = c(taxon, id, numQs)
				tmprows = rbind(tmprows, tmprow)
				} # END for (i in 1:length(xml_item))
			tdf = as.data.frame(tmprows, stringsAsFactors=FALSE)
			names(tdf) = c("species", "locus", "numQs")
			row.names(tdf) = NULL
			#head(tdf)
			} # END if ( (isblank_TF(tmp_datatype) == FALSE) && (tmp_datatype=="nucleotide") )
		} # END for (j in 1:length(xml_items_list))	

	return(tdf)
	} # END get_DNA_numQs_from_data_XML <- function(xml_items_list)







# After parse_datasets() has read in the data matrices,
# this function:

# (a) saves them out to individual matrices
# (b) concatenates different DNA loci into a supermatrix
#     (in case one wants to run a separate concatenated 
#      analysis in a different run or program)
# (c) saves tables with completeness information
# (d) generates new datasets based on some completeness
#     filters
# (e) generates .Rdata files for easy loading if
#     runslow=FALSE in the main script\
#
# The directories created are based on z_saved, z_saved_preCompleteness
# The 1, 2, 3, etc., refer to the different data matrices

save_out_matrices_and_completeness <- function(data_XMLs_preCompleteness, datastats_dir="z_data_stats", min_numseqs=0, min_fraction_to_keep_sequence=0.5, StarBeast2_TF=FALSE, taxonsets_df=NULL,  data_XMLs_fn="data_XMLs.Rdata", data_XMLs_preCompleteness_fn="data_XMLs_preCompleteness.Rdata", dnalist_conversion_output_fn="dnalist_conversion_dnalist_conversion_output_preCompleteness.Rdata", dnalist_conversion_output_preCompleteness_fn="dnalist_conversion_dnalist_conversion_output.Rdata")
	{
	defaults = '
	datastats_dir="z_data_stats"
	min_numseqs=0
	min_fraction_to_keep_sequence=0.5
	StarBeast2_TF=FALSE
	taxonsets_df=NULL
	data_XMLs = parse_datasets(seqs_df, add_morphLength=TRUE, add_morphList=TRUE, OTUs=OTUs,  printall="short", convert_ambiguous_to_IUPAC=FALSE, xml=NULL, xlsfn=xlsfn, return_charsdf=TRUE)
	
	data_XMLs_preCompleteness = data_XMLs
	data_XMLs_fn = "data_XMLs.Rdata"
data_XMLs_preCompleteness_fn = "data_XMLs_preCompleteness.Rdata"
dnalist_conversion_output_fn = "dnalist_conversion_dnalist_conversion_output_preCompleteness.Rdata"
dnalist_conversion_output_preCompleteness_fn = "dnalist_conversion_dnalist_conversion_output.Rdata"

	save_out_matrices_and_completeness_files = save_out_matrices_and_completeness(data_XMLs_preCompleteness, datastats_dir=datastats_dir, min_numseqs=0, min_fraction_to_keep_sequence=0.5, StarBeast2_TF=StarBeast2_TF, taxonsets_df=taxonsets_df,  data_XMLs_fn=data_XMLs_fn, data_XMLs_preCompleteness_fn=data_XMLs_preCompleteness_fn, dnalist_conversion_output_fn=dnalist_conversion_output_fn, dnalist_conversion_output_preCompleteness_fn=dnalist_conversion_output_preCompleteness_fn)
	'
	
	# Save to an .Rdata file for easy loading
	save(data_XMLs_preCompleteness, file=data_XMLs_preCompleteness_fn)


	#######################################################
	# Read the data_XML to a dnalist, calculate completeness statistics PRE-completeness filtering
	#######################################################
	# For StarBeast2, prune out sequences from specimenNames not listed in taxonsets_df
	# (no action taken when not a StarBeast2 analysis)
	data_XML_preCompleteness = prune_seqs_based_on_taxonsets(data_XML=data_XMLs$data_XML, taxonsets_df=taxonsets_df, StarBeast2_TF=StarBeast2_TF)

	# remove_taxon = remove the taxon name from the name for a starBeast2 locus renaming
	dnalist_conversion_output_preCompleteness = convert_data_XML_to_dnalists(data_XML=data_XML_preCompleteness, remove_taxon=TRUE, string_to_delete="_seqfrom_", name_seqs_by_taxonname=TRUE, use_single_unique_seqName_for_alignment_name=TRUE)
	save(dnalist_conversion_output_preCompleteness, file=dnalist_conversion_output_preCompleteness_fn)

	# Calculate completeness, and write stats and datasets, before completeness filtering
	completeness_df_preCompleteness = completeness_stats_from_data_XML(data_XML=data_XML_preCompleteness)
	data_XMLs_preCompleteness$completeness_df = completeness_df_preCompleteness
	preCompleteness_df_outfn = paste0(datastats_dir, "/", "completeness_df_preCompleteness.txt")
	write.table(x=completeness_df_preCompleteness, file=preCompleteness_df_outfn, append=FALSE, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
	
	
	# Write out the data files BEFORE completeness filtering
	outfns_preCompleteness = write_DNAdf_to_files(dna_taxa_list=dnalist_conversion_output_preCompleteness$dna_taxa_list, dna_dfs_list=dnalist_conversion_output_preCompleteness$dna_dfs_list, datatypes=dnalist_conversion_output_preCompleteness$datatypes, prefix="", suffix="SavedOut", write_unmerged_also_TF=TRUE, outdir="z_saved_preCompleteness")
	
	cat("\n\nwrite_DNAdf_to_files() has written...\n")
	cat(outfns_preCompleteness, sep="\n")

	#######################################################
	# Prune based on completeness, and get final completeness statistics
	#######################################################
	data_XML = prune_seqs_based_on_completeness(data_XML=data_XML_preCompleteness, min_numseqs=0, min_fraction_to_keep_sequence=0.5)
	data_XMLs$data_XML = data_XML

	completeness_df = completeness_stats_from_data_XML(data_XMLs$data_XML)
	data_XMLs$completeness_df = completeness_df
	completeness_df_outfn = paste0(datastats_dir, "/", "completeness_df.txt")
	write.table(x=completeness_df, file=completeness_df_outfn, append=FALSE, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

	# Read the data_XML to a dnalist, for filtering based on completeness etc.
	dnalist_conversion_output = convert_data_XML_to_dnalists(data_XML=data_XMLs$data_XML, remove_taxon=TRUE, string_to_delete="_seqfrom_", name_seqs_by_taxonname=TRUE, use_single_unique_seqName_for_alignment_name=TRUE)

	# Write out the data to a bunch of files
	outfns = write_DNAdf_to_files(dna_taxa_list=dnalist_conversion_output$dna_taxa_list, dna_dfs_list=dnalist_conversion_output$dna_dfs_list, datatypes=dnalist_conversion_output$datatypes, suffix="SavedOut", write_unmerged_also_TF=TRUE, outdir="z_saved")
	cat("\n\nwrite_DNAdf_to_files() has written...\n")
	cat(outfns, sep="\n")

	save(data_XMLs, file=data_XMLs_fn)
	save(dnalist_conversion_output, file=dnalist_conversion_output_fn)
	
	save_out_matrices_and_completeness_files = list()
	save_out_matrices_and_completeness_files$datastats_dir = datastats_dir

	save_out_matrices_and_completeness_files$outfns_preCompleteness = outfns_preCompleteness
	save_out_matrices_and_completeness_files$outfns = outfns
	
	save_out_matrices_and_completeness_files$data_XMLs_fn = data_XMLs_fn
	save_out_matrices_and_completeness_files$data_XMLs_preCompleteness_fn = data_XMLs_preCompleteness_fn
	save_out_matrices_and_completeness_files$dnalist_conversion_output_fn = dnalist_conversion_output_fn
	save_out_matrices_and_completeness_files$dnalist_conversion_output_preCompleteness_fn = dnalist_conversion_output_preCompleteness_fn
	return(save_out_matrices_and_completeness_files)
	} # END save_out_matrices_and_completeness


#######################################################
# Which taxa have data?
#######################################################
# Produces text files of yes/no under different conditions,
# based on what is in the XML, described in data_XMLs$matrices_stats
# 
# These are useful for pasting back into the "OTUs" column
# of the BEASTmasteR settings Excel file, for different Beast2 analyses
# (morphology-only, DNA-only, etc.)
#######################################################

which_taxa_have_data <- function(matrices_stats, count_constraints_df, OTUs_raw=NULL)
	{
	setup = '
	OTUs_df = readWorksheetFromFile(xlsfn, sheet="OTUs", startRow=15)
	head(OTUs_df)
	OTUs_raw = OTUs_df$OTUs	# raw list, for safekeeping / pasting back to Excel

	matrices_stats = data_XMLs$matrices_stats
	count_constraints_df = xml$count_constraints_df
	matrices_stats
	count_constraints_df
	
	z_Excel_OTUs_yesNo_fns = which_taxa_have_data(matrices_stats=matrices_stats, count_constraints_df=count_constraints_df, OTUs_raw=OTUs_raw)
z_Excel_OTUs_yesNo_fns
	' # END setup
	
	
	if (is.null(OTUs_raw))
		{
		OTUs_raw = matrices_stats[[1]]$completeness_df$OTU
		} # END if (is.null(OTUs_raw))
	
	if (is.null(OTUs_raw))
		{
		stoptxt = paste0("ERROR in which_taxa_have_data(): input OTUs_raw has not been specified, and matrices_stats[[1]]$completeness_df$OTU was NULL. One of these has to provide the list of OTUs.")
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		stop(stoptxt)
		} # END if (is.null(OTUs_raw))
	
	indices_to_OTUs_raw = match(x=count_constraints_df$OTUs, table=OTUs_raw)
	constraints_count = count_constraints_df$count_constraints[indices_to_OTUs_raw]

	taxa_wData_df = cbind(OTUs_raw, constraints_count)
	taxa_wData_df = as.data.frame(taxa_wData_df, stringsAsFactors=FALSE)

	dir.create(path="z_Excel_OTUs_yesNo", showWarnings=FALSE)

	for (i in 1:length(matrices_stats))
		{
		dataset_name = data_XMLs$dataset_lengths_df$dataset_names[[i]]
		indices_to_OTUs_raw = match(x=matrices_stats[[i]]$completeness_df$OTU, table=OTUs_raw)
	
		# Add a column of 0s
		cmdstr = paste0("taxa_wData_df$", dataset_name, " = rep(0, nrow(taxa_wData_df))")
		eval(parse(text=cmdstr))

		# Fill it in
		notNAs = is.na(indices_to_OTUs_raw) == FALSE
		indices_to_OTUs_raw = indices_to_OTUs_raw[notNAs]
		vals = matrices_stats[[i]]$completeness_df$PctData[notNAs]
	
		cmdstr = paste0("taxa_wData_df$", dataset_name, "[indices_to_OTUs_raw] = vals")
		eval(parse(text=cmdstr))
		}
	taxa_wData_df = dfnums_to_numeric(taxa_wData_df)
	save(taxa_wData_df, file="z_Excel_OTUs_yesNo/taxa_wData_df.Rdata")
	write.table(taxa_wData_df, file="z_Excel_OTUs_yesNo/taxa_wData_df.txt", quote=FALSE, sep="\t", row.names=FALSE)

	cols_named_morph_TF = grepl(pattern="morph", x=tolower(names(taxa_wData_df)))
	if (sum(cols_named_morph_TF) == 0)
		{
		dataYes_TFmorph = rep(FALSE, nrow(taxa_wData_df))
		} else if (sum(cols_named_morph_TF) > 1) {
		dataYes_TFmorph = rowSums(taxa_wData_df[,cols_named_morph_TF]) > 0
		} else {
		dataYes_TFmorph = taxa_wData_df[,cols_named_morph_TF] > 0
		}

	cols_named_DNA_TF = grepl(pattern="dna", x=tolower(names(taxa_wData_df)))
	if (sum(cols_named_DNA_TF) == 0)
		{
		dataYes_TFdna = rep(FALSE, nrow(taxa_wData_df))
		} else if (sum(cols_named_DNA_TF) > 1) {
		dataYes_TFdna = rowSums(taxa_wData_df[,cols_named_DNA_TF]) > 0
		} else {
		dataYes_TFdna = taxa_wData_df[,cols_named_DNA_TF] > 0
		}

	dataYes_TF1 = taxa_wData_df$constraints_count > 1

	cols_wChars_TF = rep(FALSE, ncol(taxa_wData_df))
	cols_wChars_TF[3:ncol(taxa_wData_df)] = TRUE
	if (sum(cols_wChars_TF) == 0)
		{
		dataYes_TFmore = rep(FALSE, nrow(taxa_wData_df))
		} else if (sum(cols_wChars_TF) > 1) {
		dataYes_TFmore = rowSums(taxa_wData_df[,cols_wChars_TF]) > 0
		} else {
		dataYes_TFmore = taxa_wData_df[,cols_wChars_TF] > 0
		}

	dataYes_TF = cbind(dataYes_TF1, dataYes_TFmore)

	cols_not_named_morph = names(taxa_wData_df)[cols_named_morph_TF==FALSE]
	cols_not_named_morph = cols_not_named_morph[cols_not_named_morph != "constraints_count"]
	cols_not_named_morph = cols_not_named_morph[cols_not_named_morph != "OTUs_raw"]

	if (length(cols_not_named_morph) == 0)
		{
		dataYes_TFmore_noMorph = rep(FALSE, nrow(taxa_wData_df))
		} else if (length(cols_not_named_morph) > 1) {
		dataYes_TFmore_noMorph = rowSums(taxa_wData_df[,cols_not_named_morph]) > 0
		} else {
		dataYes_TFmore_noMorph = taxa_wData_df[,cols_not_named_morph] > 0
		}




	# OTUs with any data
	use_Yes_No = rep("yes", nrow(taxa_wData_df))
	OTUs_with_any_data_TF = rowSums(dataYes_TF) > 0
	use_Yes_No[OTUs_with_any_data_TF==FALSE] = "no"
	use_Yes_No_df = as.data.frame(use_Yes_No, stringsAsFactors=FALSE)
	names(use_Yes_No_df) = "Excel_OTUs_yesNo_any_data"
	write.table(use_Yes_No_df, file="z_Excel_OTUs_yesNo/Excel_OTUs_yesNo_any_data.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)


	# OTUs with any morphological data
	use_Yes_No = rep("yes", nrow(taxa_wData_df))
	use_Yes_No[dataYes_TFmorph==FALSE] = "no"
	use_Yes_No_df = as.data.frame(use_Yes_No, stringsAsFactors=FALSE)
	names(use_Yes_No_df) = "Excel_OTUs_yesNo_any_data_named_morph"
	write.table(use_Yes_No_df, file="z_Excel_OTUs_yesNo/Excel_OTUs_yesNo_any_data_named_morph.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

	# OTUs with any character data, not named morph
	use_Yes_No = rep("yes", nrow(taxa_wData_df))
	use_Yes_No[dataYes_TFmore_noMorph==FALSE] = "no"
	use_Yes_No_df = as.data.frame(use_Yes_No, stringsAsFactors=FALSE)
	names(use_Yes_No_df) = "Excel_OTUs_yesNo_any_non-morph_data"
	write.table(use_Yes_No_df, file="z_Excel_OTUs_yesNo/Excel_OTUs_yesNo_any_non-morph_data.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

	# OTUs with any DNA data
	use_Yes_No = rep("yes", nrow(taxa_wData_df))
	use_Yes_No[dataYes_TFdna==FALSE] = "no"
	use_Yes_No_df = as.data.frame(use_Yes_No, stringsAsFactors=FALSE)
	names(use_Yes_No_df) = "Excel_OTUs_yesNo_any_data_named_DNA"
	write.table(use_Yes_No_df, file="z_Excel_OTUs_yesNo/Excel_OTUs_yesNo_any_data_named_DNA.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

	# OTUs with any character data whatsoever
	use_Yes_No = rep("yes", nrow(taxa_wData_df))
	use_Yes_No[dataYes_TFmore==FALSE] = "no"
	use_Yes_No_df = as.data.frame(use_Yes_No, stringsAsFactors=FALSE)
	names(use_Yes_No_df) = "Excel_OTUs_yesNo_any_char_data"
	write.table(use_Yes_No_df, file="z_Excel_OTUs_yesNo/Excel_OTUs_yesNo_any_char_data.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
	
	z_Excel_OTUs_yesNo_fns = list.files(path="z_Excel_OTUs_yesNo")
	
	cat("\n\nwhich_taxa_have_data() wrote these files to subdirectory 'z_Excel_OTUs_yesNo':\n\n")
	cat(z_Excel_OTUs_yesNo_fns, sep="\n")
	cat("\n")
	return(z_Excel_OTUs_yesNo_fns)
	} # END which_taxa_have_data <- function(matrices_stats, count_constraints_df, OTUs_raw=NULL)

