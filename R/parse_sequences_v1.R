#######################################################
# Functions for parsing sequence data and adding to XML
#######################################################

convert_NEXUS_to_FASTA <- function(fasta_fn, out_nexus_aln_fn=NULL)
	{
	require(ape)
	require(seqinr)
	
	if (is.null(out_nexus_aln_fn))
		{
		out_nexus_aln_fn = paste0(fasta_fn, ".nex")
		}
		
	# Write out to NEXUS
	seqs = seqinr::read.fasta(fasta_fn)
	ape::write.nexus.data(x=seqs, file=out_nexus_aln_fn)
	txt = paste0("\nconvert_NEXUS_to_FASTA() converted a FASTA file named '", fasta_fn, "' to a NEXUS data file named '", out_nexus_aln_fn, "'.\n")
	cat(txt)
	
	return(out_nexus_aln_fn)
	} # END convert_NEXUS_to_FASTA <- function(fasta_fn, out_nexus_aln_fn=NULL)


collapse_seqs <- function(dataset)
	{
	txtseqs = unlist(lapply(X=dataset, FUN=paste0, collapse=""))
	return(txtseqs)
	}


# Make XML for a single sequence
# For use in sapply
make_seq_XML <- function(txtseq, OTU_name, dataset_name, totalcount)
	{
	id = paste0(dataset_name, "_seqfrom_", OTU_name)
	seq_XML = xmlNode(name="sequence", attrs=list(id=id, taxon=OTU_name, totalcount=totalcount, value=txtseq))
	
	return(seq_XML)
	}
	
# 	tipdates_txt = make_txt_tipdates_list(OTUs, tipdates, leading_space="\t\t\t", trailing=",\n")
# 	tipdates_txt = paste0(tipdates_txt, "\n", sep="")
# 	txt3 = " Specify the taxonSet which the tipdate names refer to. "
# 	XML_comment3 = xmlCommentNode(txt3)
# 	taxa_names_source_XML = xmlNode(name="taxa", attrs=list(idref=taxon_name))
# 	
# 	XML_tipdates = xmlNode(name="trait", tipdates_txt, attrs=list(id=name, spec="beast.evolution.tree.TraitSet", traitname=traitname))
# 	XML_tipdates
# 



print_unique_numstates <- function(numstates_morph_list, printall="short")
	{
	# Now, go through each number of characters, and get the ambiguities
	# present in each.
	unique_numstates = unique(numstates_morph_list)
	unique_numstates = unique_numstates[order(unique_numstates)]
	
	for (numstates in unique_numstates)
		{
		match_TF = numstates_morph_list == numstates
		match_indices = c(1:length(numstates_morph_list))[match_TF]

		if (printall != "none")
			{
			if (numstates <= 2)
				{
				txt = " (this includes invariant and autapomorphic characters)"
				} else {
				txt = ""
				}
			
			cat("\nThese characters have ", numstates, " unique states (skipping '?' or '-')", txt, ":\n", sep="")
			cat(match_indices, "\n", sep=" ")
			} # END if (printall != "none")
		}
	return(unique_numstates)
	}




make_BEAST2_userDataTypes <- function(numstates_morph_list, nexd6, dataset_name="", ordering="unordered", morph_transition_rates="equal", baseFreqs="equal", numGammaCat=1, clockModel_name="shared_clock", clockModel_relRate=NA, gammaShape_suffix=NA, add_morphList=TRUE, printall="short")
	{
	defaults='
	numstates_morph_list = numstates_morph_list_subset
	nexd6 = nexd6
	printall="short"
	'
	
	# Warning for when shared_clock in clockModel_name is lost somehow...
	if (is.null(clockModel_name) == TRUE)
		{
		txt = "WARNING from make_BEAST2_userDataTypes(): input 'clockModel_name' was NULL, so setting it to 'shared_clock'."
		cat("\n")
		cat(txt)
		cat("\n")
		warning(txt)
		
		clockModel_name = "shared_clock"
		}
	
	
	# Error check on ordering
	if ( (ordering != "unordered") && (ordering != "ordered") )
		{
		txt = paste0("\n\nERROR in make_BEAST2_userDataTypes(): 'ordering' must be either\n'unordered' (default) or 'ordered'. Instead, you have: ", ordering, "\n\n2016-09-24 update: The 'ordering' column is deprecated, instead the 'model' column is checked for 'Mk_unord' or 'Mk_ord'.\n\n")
		cat(txt)
		stop(txt)
		}
	
	# Keep track of the list of morphology sections
	# (for later use in making the likelihood models etc.)
	morphList = NULL


	
	# Error check on characters with only 0 or 1 observed states.
	#print(numstates_morph_list)
	#print("here1")
	
	chars_wLT_2_states_df = NULL
	TF = numstates_morph_list < 2
	if (sum(TF) > 0)
		{
		txt = paste0("\n\nWARNING in make_BEAST2_userDataTypes(): input 'numstates_morph_list' has\n", sum(TF), "  characters with invariant states.\n\n")
		txt1 = paste0("Printing numstates_morph_list[numstates_morph_list < 2] (saving as chars_wLT_2_states_df):\n")
		cat(txt)
		cat(txt1)
		
		names_as_nums = 1:length(numstates_morph_list)
		chars_wLT_2_states_df = as.data.frame(matrix(data=numstates_morph_list, nrow=1), row.names=NULL)
		names(chars_wLT_2_states_df) = names_as_nums
		chars_wLT_2_states_df = chars_wLT_2_states_df[TF]
		print(chars_wLT_2_states_df)
		cat("\n\n")
		
		#stop(txt)
		} # END if (sum(TF) > 0)
	
	# Ordered characters must have 3 or more states
	if (ordering == "ordered")
		{
		TF = numstates_morph_list == 2
		if (sum(TF) > 0)
			{
			errortxt = paste0("\n\nERROR in make_BEAST2_userDataTypes(): You have included ", sum(TF), " characters with only 2 states\nin your setup of the ordered morphology characters. Ordered characters must have\nthree or more states, if you think about it.\n")
			cat(errortxt)

			txt = paste0("These are the offending characters, of the total of ", length(TF), " characters you specified to be ordered:\n")
			cat(txt)

			txt = paste0("(may or may not correspond to original character numbers; these are the index numbers of just this subset)\n")
			cat(txt)

			charnums = (1:length(TF))[TF]
			charnums_txt = paste(charnums, sep="", collapse=", ")
			cat(charnums_txt)
			cat("\n\n")
			stop(errortxt)
			} # END if (sum(TF) > 0)
		} # END if (ordering == "ordered")
	
	# Correct the 
	numstates_morph_list[numstates_morph_list < 2] = 2
	unique_numstates = print_unique_numstates(numstates_morph_list, printall=printall)
	
	# subset the nexus dataframe to just characters (morphological, not DNA) under consideration
	nexd4 = nexd6[1:length(numstates_morph_list), ]
	
	# output nexus df
	nexd5 = nexd4
	
	# Set up a null list
	userDataType_nodes_list = NULL
	gnum = 0

	# Get the NEXUS/TNT data codes
	tmp_charstate_codes = charstate_codes()
	
	# set up the multistate BEAST codes (Y backwards; skip Z)
	multistate_BEAST_codes = rev(LETTERS)[-1]
	
	for (i in 1:length(unique_numstates))
		{
		numstates = unique_numstates[i]
		
		# Get all the characters for all of the sites with
		# unique_numstates[i] states:
		all_cells = as.character(unlist(nexd4[numstates_morph_list==numstates, ]))
		
		unique_all_cells = unique(all_cells)
		
		# Morph2
		statenodes_list = NULL
		snum = 0

		# put in the single-state characters
		# put in 0, 1, etc.
		for (j in 1:numstates)
			{
			codeMap_txt = paste0(tmp_charstate_codes[j], "=", tmp_charstate_codes[j])
			if (j == 1)
				{
				codeMap_txts = codeMap_txt
				} else {
				codeMap_txts = paste0(codeMap_txts, ", ", codeMap_txt)
				} # END if (j == 1)
			}
		
		# get the multi-state characters
		# (which have nchar > 1)
		multistate_characters_TF = sapply(X=unique_all_cells, FUN=nchar) > 1
		multistate_characters = unique_all_cells[multistate_characters_TF]
		multistate_characters = multistate_characters[order(multistate_characters)]

		# Set up the "Z" state (totally ambiguous; also covers "?" and "-")
		Zstate = paste(tmp_charstate_codes[1:numstates], collapse="", sep="")
		
		# Remove the Zstate from multistate_characters
		multistate_characters = multistate_characters[multistate_characters != Zstate]



		#######################################################
		# Remove the states that will be coded as "Z" and convert to Z
		# (do this BEFORE you do the smaller ambiguous characters, or life will SUCK!)
		#######################################################

		# Edit the datamatrix to have the BEAST Zstate coding
		tmppattern = Zstate
		tmpreplacement = "Z"
		nexd5[numstates_morph_list==numstates, ] = mapply(gsub, tmppattern, tmpreplacement, nexd5[numstates_morph_list==numstates, ])
		tmptxt = paste0(strsplit(Zstate, split="")[[1]], collapse=" ", sep="")
		Ztxt = paste0("Z=", tmptxt) 
		
		# Edit the datamatrix to have the BEAST Zstate coding
		tmppattern = "\\?"
		tmpreplacement = "Z"
		nexd5[numstates_morph_list==numstates, ] = mapply(gsub, tmppattern, tmpreplacement, nexd5[numstates_morph_list==numstates, ])
		Qtxt = paste0("?=", tmptxt)
		
		# Edit the datamatrix to have the BEAST Zstate coding
		tmppattern = "-"
		tmpreplacement = "Z"
		nexd5[numstates_morph_list==numstates, ] = mapply(gsub, tmppattern, tmpreplacement, nexd5[numstates_morph_list==numstates, ])
		dashtxt = paste0("-=", tmptxt)
		
		ZQdtxt = paste(dashtxt, Qtxt, Ztxt, collapse="", sep=", ")
		ZQdtxt

		#######################################################
		# Order by REVERSE NCHAR LENGTH
		#######################################################
		multistate_characters = sort(multistate_characters)		# alphabetical sorting, but will impose size sorting on this
		
		# sorting by nchar
		nchars_multistate_characters = nchar(multistate_characters)
		reverse_order_by_size = rev(order(nchars_multistate_characters))
		multistate_characters = multistate_characters[reverse_order_by_size]
		multistate_characters


		#######################################################
		# Now recode non-Z ambiguous characters, by reverse nchar() size
		#######################################################		
		# If there are still multistate characters, code those up from "Y" backwards
		if (length(multistate_characters) > 0)
			{
			for (j in 1:length(multistate_characters))
				{
				
				multistate_character = multistate_characters[j]
				multistate_character_split = strsplit(multistate_character, split="")[[1]]
				multistate_character_spaced = paste(multistate_character_split, sep="", collapse=" ")
				multistate_BEAST_code = multistate_BEAST_codes[j]

				# Put the multistate into XML
				#statenodes_list[[(snum=snum+1)]] = xmlNode("ambiguity", attrs=c(code=multistate_BEAST_code, states=multistate_character))
				codeMap_txt = paste0(multistate_BEAST_code, "=", multistate_character_spaced, sep="") 
				codeMap_txt
				codeMap_txts = paste0(codeMap_txts, ", ", codeMap_txt)
				codeMap_txts
				
				# Check that matrix cells are being converted to length 1
				tc = unlist(nexd5)
				tc = tc[nchar(tc) > 1]
				if (printall != "none")
					{
					cat("\nAll characters of length >1:\n", sep="")
					cat(paste0(unique(tc), sep=", ", collapse=""))
					} # END if (printall != "none")

				tc = unlist(nexd5[numstates_morph_list==numstates, ])
				tc = tc[nchar(tc) > 1]
				if (printall != "none")
					{
					cat("\nRemaining characters of length >1:\n", sep="")
					cat(paste0(unique(tc), sep=", ", collapse=""))
					} # END if (printall != "none")


				# Edit the datamatrix to have the BEAST multistate coding
				tmppattern = multistate_character
				tmpreplacement = multistate_BEAST_code
				
				# Extract just the row with the number of states specified by "i"
				#tmprows = nexd5[numstates_morph_list==numstates, ]
				# Now find the EXACT cell matches (NOT just the text-replace substitutions, which leads to e.g. 012 --> W2 (?)
				#match_TF = tmprows == tmppattern
				#tmprow[match_TF] = tmpreplacement
				
				if (printall != "none")
					{
					cat("i=", i, "	numstates=", numstates, "	find=", tmppattern, "	replace=", tmpreplacement, "\n", sep="")
					} # END if (printall != "none")

				
				# This was probably BAD
				nexd5[numstates_morph_list==numstates, ] = mapply(gsub, tmppattern, tmpreplacement, nexd5[numstates_morph_list==numstates, ])
				
				# Better:
				#nexd5[numstates_morph_list==numstates, ] = tmprow
				}
				
			# Error check:
			txt_states = unlist(nexd5[numstates_morph_list==numstates, ])
			cells_nchar_GT0_TF = nchar(txt_states) > 1
			if (sum(cells_nchar_GT0_TF))
				{
				errortxt = paste("ERROR -- you size have cells with nchar() > 1, for unique_numstates[", i, "] = ", unique_numstates[i], sep="")
				print(txt_states[cells_nchar_GT0_TF])
				stop(errortxt)
				}
			}
		

		# Put the Z state into XML
		#statenodes_list[[(snum=snum+1)]] = xmlNode("ambiguity", attrs=c(code="Z", states=Zstate))
		codeMap_txts = paste(codeMap_txts, ", ", ZQdtxt, sep="")
		codeMap_txts
		
		# Set up the datatype name
		name_of_userDataType = paste("datatype_", dataset_name, "_morph", numstates, "_", ordering, sep="")
		txt = paste0(" Datatype for a ", numstates, "-state ", ordering, " character ")
		
		# Save in morphList
		numchars = sum(numstates_morph_list==numstates)
		
		# 2016-09-24: added morph_transition_rates, baseFreqs
		# 2016-12-19: added clockModel_relRate
		morphList_row = c(dataset_name, name_of_userDataType, numstates, ordering, morph_transition_rates, baseFreqs, numchars, numGammaCat, clockModel_name, clockModel_relRate, gammaShape_suffix)
		morphList = rbind(morphList, morphList_row)
		
		# Add this 
		userDataType_node = xmlNode("userDataType", attrs=list(id=name_of_userDataType, spec="beast.evolution.datatype.UserDataType", states=numstates, codeMap=codeMap_txts, codelength="1"))
		
		if (numstates == 2)
			{
			txt2 = "(all 2-state characters are unordered, duh!)"
			userDataType_nodes_to_add = list(bl(), xmlCommentNode(txt), xmlCommentNode(txt2), userDataType_node)
			} else {
			userDataType_nodes_to_add = list(bl(), xmlCommentNode(txt), userDataType_node)
			}
		
		userDataType_nodes_list = c(userDataType_nodes_list, userDataType_nodes_to_add)

		} # END for (i in 1:length(unique_numstates))
	
	result = NULL
	result$nexd5 = nexd5
	result$userDataType_nodes_list = userDataType_nodes_list
# 	print("Checkpoint_1 in make_BEAST2_userDataTypes()...")
	
	if (add_morphList == TRUE)
		{
		row.names(morphList) = NULL
		morphList = as.data.frame(morphList, row.names=NULL, stringsAsFactors=FALSE)
# 		print(names(morphList))
# 		print(morphList)
		
		
		names(morphList) = c("dataset_name", "name_of_userDataType", "numstates", "ordering", "morph_transition_rates", "baseFreqs", "numchars", "numGammaCat", "clockModel_name", "clockModel_relRate", "gammaShape_suffix")
		class(morphList$numstates) = "numeric"
		class(morphList$numchars) = "numeric"
		result$morphList = morphList
		result$chars_wLT_2_states_df = chars_wLT_2_states_df
		}
	
	extract='
	nexd5 = result$nexd5
	userDataType_nodes_list = result$userDataType_nodes_list
	morphList = result$morphList
	chars_wLT_2_states_df = result$chars_wLT_2_states_df

	nexd7 = result2$nexd5
	userDataType_nodes_list = result2$userDataType_nodes_list
	morphList = result$morphList
	chars_wLT_2_states_df = result2$chars_wLT_2_states_df
	'
	
	return(result)
	}



# Options for ascertainment model:
# ascertainment="Mk": Markov-k model for k states, no correction for ascertainment bias
# ascertainment="Mkv": Markov-k model, conditioning on unobservability of invariant sites
# Now implemented:
# ascertainment="MkParsInf": Markov-k model, conditioning on all sites being 
#                         parsimony-informative (all invariant and autapomorphic sites
#                         excluded)
# ascertainment="noabsencesites": Markov-k model, sites of all-0 are disallowed
#                                 (e.g. for restriction-site data, where all-absence
#                                 -- all-0 -- site data is disallowed). See e.g. MrBayes.
# ascertainment="nopresencesites": Markov-k model, sites of all-1 are disallowed
#                                 (e.g. for restriction-site data, where all-absence
#                                 -- all-0 -- site data is disallowed). See e.g. MrBayes.
# 
# MkA: Pyron (2016)'s F81-like model for binary characters. 
# pi0 = basefrequency of state 0, and rate of 1->0
# pi1 = basefrequency of state 1, and rate of 0->1
# 
# 
write_BEAST2_morphology_characters <- function(nexd7, numstates_morph_list, dataset_name="", ordering="unordered", morph_transition_rates="equal", baseFreqs="equal", ascertainment="Mk", max_num_patterns=2500, assumed_nstates="")
	{
	junk='
	nexd7=nexdf
	numstates_morph_list=numstates_morph_list
	ordering="unordered"
	'
	
	# Correct numstates_morph_list to set characters with 0 or 1 states to 2
	numstates_morph_list[numstates_morph_list<2] = 2
	
	allowed_ascertainment_models = c("Mk", "Mkv", "noabsencesites", "nopresencesites", "MkInf", "MkParsInf")
	
	# Error check on ordering
	if ( (ordering != "unordered") && (ordering != "ordered") )
		{
		txt = paste0("\n\nSTOP ERROR in write_BEAST2_morphology_characters(): 'ordering' must be either\n'unordered' (default) or 'ordered'. Instead, you have: ", ordering, "\n\n")
		cat(txt)
		stop(txt)
		}
	
	
	# Error check on characters with only 0 or 1 observed states.
	TF = numstates_morph_list < 2
	if (sum(TF) > 0)
		{
		txt = paste0("\n\nWARNING in write_BEAST2_morphology_characters(): input 'numstates_morph_list' has\n", sum(TF), "  characters with only 0 or 1 states.\n\n")
		txt1 = paste0("Printing numstates_morph_list:\n")
		cat(txt)
		cat(txt1)
		
		names_as_nums = 1:length(numstates_morph_list)
		numstates_morph_list_df = as.data.frame(matrix(data=numstates_morph_list, nrow=1), row.names=NULL)
		names(numstates_morph_list_df) = names_as_nums
		print(numstates_morph_list_df)
		cat("\n\n")
		
		#stop(txt)
		}
		
	
	unique_numstates = unique(numstates_morph_list)
	unique_numstates = unique_numstates[order(unique_numstates)]
	
	# Initialize a list to contain all of the alignment nodes, one for each
	# type of morphology character
	morph_alignment_nodes_list = NULL
	anum = 0
	
	
	for (i in 1:length(unique_numstates))
		{
		# Get the list of character numbers (in these ordered/unordered characters)
		matchTF = numstates_morph_list %in% unique_numstates[i]
		list_of_charnums = (1:length(numstates_morph_list))[matchTF]
		list_of_charnums_txt = paste(list_of_charnums, sep=", ", collapse=", ")
		list_of_charnums_txt = paste0(" ", list_of_charnums_txt, " ", sep="")
		list_of_charnums_caption = paste0("Indices of the characters in the ", length(numstates_morph_list), " ", ordering, " characters subset from dataset: ", dataset_name, ": ")
		
		
		
		partition_name = paste0(dataset_name, "_morph", unique_numstates[i], "_", ordering)
		name_of_userDataType = paste("datatype_", dataset_name, "_morph", unique_numstates[i], "_", ordering, sep="")
		
		
		numstates = unique_numstates[i]
		tmpnexd = nexd7[numstates_morph_list==numstates, ]

		# Create the string of sequence data for each taxon,
		# for the characters with unique_numstates[i] states
		tmpseqs = apply(X=tmpnexd, MARGIN=2, FUN=paste, collapse="")
		
		if (isblank_TF(ascertainment))
			{
			txt = "WARNING in write_BEAST2_morphology_characters(): 'ascertainment' was blank, which probably means you have 'ascertainment' blank in the 'data' worksheet of the Excel settings file, for the morphology row(s). Defaulting to 'Mk' model (no ascertainment bias correction), which you may not want."
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			#stop(txt)
			
			ascertainment = "Mk"
			} # END if (isblank_TF(ascertainment))
		
		
		# Add the sites to condition on
		if (ascertainment == "Mk")
			{
			ascertainment_txt = " The ascertainment-bias model for these characters has been set to 'Mk' (Markov-k). This means *NO* correction for unobservable character patterns will be performed.\n"
			} # END if (ascertainment == "Mk")
		

		if (ascertainment == "Mkv")
			{
			states_to_exclude = seq(0, (numstates-1), 1)
			states_to_exclude_txt = paste(states_to_exclude, sep="", collapse="")
			states_to_exclude_txt
			ascertainment_txt = paste0(" The ascertainment-bias model for these characters has been set to 'Mkv' (Markov-k-variable). This corrects for the unobservability of invariant characters.\n\nOther ways of describing the ascertainment-bias correction:\n\n- The likelihood of these characters is conditioned on them being variable characters.\n- The characters are ascertained in a biased way, in that only variable characters are ascertained.\n\nIn Beast2, this 'ascertainment bias' is corrected for by setting up ", numstates, " dummy characters and adding them to the front of the alignment. In this case, for a ", numstates, "-state character, that means adding '", states_to_exclude_txt, "' to the beginning of the character sequence for each taxon. During calculation of the likelihood, the likelihood of those ", numstates, " invariant 'character patterns' is calculated, and the likelihood of each observed site is normalized by dividing by (1-likelihood of unobservable patterns).\n\nThe basics of ascertainment bias correction are described in p. 163, equations 8-9, of:\n\nFelsenstein, Joseph (1992). Phylogenies from Restriction Sites: A Maximum-Likelihood Approach. Evolution, 46(1), 159-173. February 1992. http://www.jstor.org/stable/2409811\n\n See also: https://groups.google.com/forum/#!topic/raxml/QOgSN8jRA2o\n ")
			
			# Make the dummy characters, and add them
			dummy_characters = rep(states_to_exclude_txt, length(tmpseqs))			
			tmpseqs = paste(dummy_characters, tmpseqs, sep="")
			} # END if (ascertainment == "Mkv")

		# MkInf = old version of ParsInf name
		if ((ascertainment == "MkParsInf") || (ascertainment == "MkInf"))
			{
			allowed_ascertainment_models_txt = paste(allowed_ascertainment_models, sep=", ", collapse=", ")
			error_txt = paste0("\n\nSTOP ERROR IN write_BEAST2_morphology_characters(): The ascertainment-bias model for these characters has been set to 'MkParsInf', but this model is not yet implemented in BEASTmasteR. Allowed models: ", allowed_ascertainment_models_txt, "\n\n")
			ascertainment_txt = " The ascertainment-bias model for these characters has been set to 'MkParsInf' (Markov-k, with correction for observing only parsimony-informative characters). This correction will be attempted, but will fail if there are too many possible unobservable site patterns.\n"
			#ascertainment_txt = paste0(" The ascertainment-bias model for these characters has been set to 'MkParsInf' (Markov-k, only parsimony-informative characters allowed, i.e. no autapomorphies). This model is not yet implemented in BEASTmasteR because ", "(1) It's a pain (you have to list in the XML all of the patterns you are *not* able to observe, i.e., an autapomorphy on each possible taxon; this soon gets ridiculous for many taxa/many character states). ", "(2) You should probably describe/code autapomorphies in morphological characters anyway: these are (slightly) informative in a likelihood/Bayesian framework, they can distinguish OTUs that would otherwise be identical, ", "they give a little more information about the amount of morphological divergence and thus time elapsed (I credit Randy Irmis for pointing this out to me). Finally, it's always possible that new taxa/new descriptions will make autapomorphies into synapomorphies. (I credit Randy Irmis for convincing me of this.) ", "(3) As far as I can tell from the MrBayes 3.1.2 code, MrBayes is only set up to do parsimony-informative ascertainment bias correction for 2-state (binary) characters anyway: i.e., if you thought you were doing it for 3-state, 4-state characters etc., you probably weren't actually. This is also hinted at in the MrBayes manual: https://groups.google.com/forum/#!topic/raxml/QOgSN8jRA2o . ", "(4) The MrBayes manual suggests that parsimony-informative ascertainment bias doesn't matter much anyway for large phylogenies. (5) Even if it does matter for the likelihoods, it's debatable whether or not the dating and topology would be effected much, ", " except in extreme cases, because the rates should adjust to match date calibrations regardless of whether the rate is absolutely correct. All of this constitutes the casual opinion of Nick Matzke and I Might Be Wrong.\n\n")
			#cat(error_txt)
			cat("\n\n")
			cat(ascertainment_txt)
			cat("\n\n")
			#stop(error_txt)
			
			if (ordering == "ordered")
				{
				error_txt = "STOP ERROR in write_BEAST2_morphology_characters(): The ascertainment-bias model for these characters has been set to 'MkParsInf', but the ordering is 'ordered'. MkParsInf for ordered characters will be different than for unordered characters, and I have not implemented it yet.\n"
				cat("\n\n")
				cat(error_txt)
				cat("\n\n")
				stop(error_txt)
				} # END if (ordering == "ordered")
			
			# Otherwise, implement the characters and fuse them to the beginning of the data matrix:
			ntaxa = length(tmpseqs)
			nstates = numstates
			# Returns a matrix
			newcols = list_unobservable_patterns_ParsInf(ntaxa=ntaxa, nstates=nstates, printflag=TRUE, max_num_patterns=max_num_patterns, ordering=ordering, assumed_nstates=assumed_nstates)
			
			num_parsInf_dummy_characters = ncol(newcols)
			
			# Fuse the matrix columns
			dummy_characters = apply(X=newcols, MARGIN=1, FUN=paste, collapse="")
			tmpseqs = paste(dummy_characters, tmpseqs, sep="")
			} # END if (ascertainment == "MkParsInf")


		if (ascertainment == "noabsencesites")
			{
			states_to_exclude = 0
			states_to_exclude_txt = paste(states_to_exclude, sep="", collapse="")
			states_to_exclude_txt
			ascertainment_txt = paste0(" The ascertainment-bias model for these characters has been set to 'noabsencesites'. This corrects for the unobservability of an all-0 column (all absent). This is typically used for e.g. restriction site data, where the character will not be ascertained (observed) at all if the restriction enzyme fails to recognize and cut anywhere.\n\nOther ways of describing the ascertainment-bias correction:\n\n- The likelihood of these characters is conditioned on them having at least one non-zero state.\n- The characters are ascertained in a biased way.\n\nIn Beast2, this 'ascertainment bias' is corrected for by setting up ", 1, " dummy character and adding it to the front of the alignment. In this case, for a ", numstates, "-state character, that means adding '", states_to_exclude_txt, "' to the beginning of the character sequence for each taxon. During calculation of the likelihood, the likelihood of that all-0 'character patterns' is calculated, and the likelihood of each observed site is normalized by dividing by (1-likelihood of unobservable pattern).\n\nThe basics of ascertainment bias correction are described in p. 163, equations 8-9, of:\n\nFelsenstein, Joseph (1992). Phylogenies from Restriction Sites: A Maximum-Likelihood Approach. Evolution, 46(1), 159-173. February 1992. http://www.jstor.org/stable/2409811\n\n See also: https://groups.google.com/forum/#!topic/raxml/QOgSN8jRA2o\n ")
			
			# Make the dummy characters, and add them
			dummy_characters = rep(states_to_exclude_txt, length(tmpseqs))			
			tmpseqs = paste(dummy_characters, tmpseqs, sep="")
			} # END if (ascertainment == "noabsencesites")
		if (ascertainment == "nopresencesites")
			{
			states_to_exclude = 1
			states_to_exclude_txt = paste(states_to_exclude, sep="", collapse="")
			states_to_exclude_txt
			ascertainment_txt = paste0(" The ascertainment-bias model for these characters has been set to 'nopresencesites'. This corrects for the unobservability of an all-1 column (all present). Correction for ascertainment bias against all-absence", " is typically used for e.g. restriction site data, where the character will not be ascertained (observed) at all if the restriction enzyme fails to recognize and cut anywhere. I'm not sure where all-presence ascertainment bias occurs, but I'm just replicating MrBayes's options here.\n\nOther ways of describing the ascertainment-bias correction:\n\n- The likelihood of these characters is conditioned on them having at least one non-1 state.\n- The characters are ascertained in a biased way.\n\nIn Beast2, this 'ascertainment bias' is corrected for by setting up ", 1, " dummy character and adding it to the front of the alignment. In this case, for a ", numstates, "-state character, that means adding '", states_to_exclude_txt, " ' to the beginning of the character sequence for each taxon. During calculation of the likelihood, the likelihood of that all-1 'character patterns' is calculated, and the likelihood of each observed site is normalized by dividing by (1-likelihood of unobservable pattern).\n\nThe basics of ascertainment bias correction are described in p. 163, equations 8-9, of:\n\nFelsenstein, Joseph (1992). Phylogenies from Restriction Sites: A Maximum-Likelihood Approach. Evolution, 46(1), 159-173. February 1992. http://www.jstor.org/stable/2409811\n\n See also: https://groups.google.com/forum/#!topic/raxml/QOgSN8jRA2o\n ")
			
			# Make the dummy characters, and add them
			dummy_characters = rep(states_to_exclude_txt, length(tmpseqs))			
			tmpseqs = paste(dummy_characters, tmpseqs, sep="")
			} # END if (ascertainment == "nopresencesites")		
		
		# Check for incorrect ascertainment pattern:
		if ((ascertainment %in% allowed_ascertainment_models) == FALSE)
			{
			allowed_ascertainment_models_txt = paste(allowed_ascertainment_models, sep=", ", collapse=", ")
			error_txt = paste0("\n\nSTOP ERROR IN write_BEAST2_morphology_characters(): The ascertainment-bias model for these characters has been set to ascertainment='", ascertainment, "', but this model option does not exist. Allowed models: ", allowed_ascertainment_models_txt, "\n\nAlso, check that, for the 'ascertainment' column in the 'data' worksheet, the morphology rows in use each have the 'ascertainment' cell filled in.\n\n")
			stop(error_txt)
			} # END if ((ascertainment %in% allowed_ascertainment_models) == FALSE)
		
		# For each sequence in the list (also, dataType node)
		sequence_nodes_list = NULL
		snum = 0
		
		
		# Create a dataType node
		#  <userDataType idref="datatype_Morph6" />
		sequence_nodes_list[[(snum=snum+1)]] = xmlNode(name="userDataType", attrs=list(idref=name_of_userDataType) )
			
		# This creates a txt node
		# s3 = xmlNode("sequence", .children=list(xmlNode("taxon", attrs=c(idref="Osteolepiformes")), "ZZZZZZZZ"))
		
		for (j in 1:length(tmpseqs))
			{
			taxon_name = colnames(nexd7)[j]
			idval = paste0("seq", unique_numstates[i], "_", ordering, "_from_", dataset_name, "_for_", taxon_name, sep="")
			sequence_nodes_list[[(snum=snum+1)]] = xmlNode(name="sequence", attrs=list(id=idval, taxon=taxon_name, totalcount=unique_numstates[i], value=tmpseqs[j]) )
			}
		#sequence_nodes_list
		
		# Put the dataType node and the sequence nodes as children of an alignment node
		
		morph_alignment_nodes_list[[(anum=anum+1)]] = bl()		
		morph_alignment_nodes_list[[(anum=anum+1)]] = xmlCommentNode(list_of_charnums_caption)
		morph_alignment_nodes_list[[(anum=anum+1)]] = xmlCommentNode(list_of_charnums_txt)
		morph_alignment_nodes_list[[(anum=anum+1)]] = bl()
		morph_alignment_nodes_list[[(anum=anum+1)]] = xmlCommentNode(ascertainment_txt)
		morph_alignment_nodes_list[[(anum=anum+1)]] = bl()
		
		# Write the "data" / ascertained alignment XML
		if (ascertainment == "Mk")
			{
			# Original
			#morph_alignment_nodes_list[[(anum=anum+1)]] = xmlNode(name="data", attrs=list(id=partition_name, dataType="user defined"), .children=sequence_nodes_list)
			
			# DO NOT STRIP invariant sites from the likelihood (strip="true") when ascertainment is Mk
			# (so: strip=false)
			morph_alignment_nodes_list[[(anum=anum+1)]] = xmlNode(name="data", attrs=list(id=partition_name, dataType="user defined", ascertained="false", statecount=numstates, strip="false", spec="beast.evolution.alignment.AscertainedAlignment"), .children=sequence_nodes_list)
			} # END if (ascertainment == "Mk")

		if (ascertainment == "Mkv")
			{
			# DO STRIP invariant sites from the likelihood (strip="true") when ascertainment is Mkv
			# (so: strip=true)
			# But: the R-script/user should do this

			# Exclude from 0 to the first state that is NOT excluded
			# Using 0-based counting
			excludefrom = 0
			excludeto = nchar(states_to_exclude_txt) - 0 
			excludeevery = 1
			morph_alignment_nodes_list[[(anum=anum+1)]] = xmlNode(name="data", attrs=list(id=partition_name, dataType="user defined", ascertained="true", statecount=numstates, excludefrom=excludefrom, excludeto=excludeto, excludeevery=excludeevery, strip="false", spec="beast.evolution.alignment.AscertainedAlignment"), .children=sequence_nodes_list)
			} # END if (ascertainment == "Mkv")


		if ((ascertainment == "MkParsInf") || (ascertainment == "MkInf"))
			{
			# DO STRIP invariant sites from the likelihood (strip="true") when ascertainment is MkParsInf
			# (so: strip=true)
			# But: the R-script/user should do this

			# Exclude from 0 to the first state that is NOT excluded
			# Using 0-based counting
			excludefrom = 0
			excludeto = num_parsInf_dummy_characters - 0 
			excludeevery = 1
			morph_alignment_nodes_list[[(anum=anum+1)]] = xmlNode(name="data", attrs=list(id=partition_name, dataType="user defined", ascertained="true", statecount=numstates, excludefrom=excludefrom, excludeto=excludeto, excludeevery=excludeevery, strip="false", spec="beast.evolution.alignment.AscertainedAlignment"), .children=sequence_nodes_list)
			} # END if (ascertainment == "Mkv")


		if (ascertainment == "noabsencesites")
			{
			# Exclude from 0 to the first state that is NOT excluded
			# Using 0-based counting
			excludefrom = 0
			excludeto = nchar(states_to_exclude_txt) - 0 
			excludeevery = 1
			morph_alignment_nodes_list[[(anum=anum+1)]] = xmlNode(name="data", attrs=list(id=partition_name, dataType="user defined", ascertained="true", statecount=numstates, excludefrom=excludefrom, excludeto=excludeto, excludeevery=excludeevery, strip="false", spec="beast.evolution.alignment.AscertainedAlignment"), .children=sequence_nodes_list)
			} # END if (ascertainment == "noabsencesites")

		if (ascertainment == "nopresencesites")
			{
			# Exclude from 0 to the first state that is NOT excluded
			# Using 0-based counting
			excludefrom = 0
			excludeto = nchar(states_to_exclude_txt) - 0 
			excludeevery = 1
			morph_alignment_nodes_list[[(anum=anum+1)]] = xmlNode(name="data", attrs=list(id=partition_name, dataType="user defined", ascertained="true", statecount=numstates, excludefrom=excludefrom, excludeto=excludeto, excludeevery=excludeevery, strip="false", spec="beast.evolution.alignment.AscertainedAlignment"), .children=sequence_nodes_list)
			} # END if (ascertainment == "nopresencesites")


		}
	return(morph_alignment_nodes_list)
	} # END write_BEAST2_morphology_characters





# Read in sequences, correct them to meet Beast2 requirements, 
# format for XML, output to XML tag 

# If add_morphLength=TRUE, add the length of the morphology data to the xml output
# CUT: check_numstates If TRUE (default), then conflicts between the observed number of states, and the theoretical number
parse_datasets <- function(seqs_df, add_morphLength=TRUE, add_morphList=TRUE, OTUs=NULL, printall="short", convert_ambiguous_to_IUPAC=FALSE, xml=NULL, xlsfn=NULL, return_charsdf=TRUE)#, check_numstates=TRUE)
	{
	defaults='
	# Example script defaults
	add_morphLength=TRUE
	add_morphList=TRUE
	OTUs=OTUs
	printall="short"
	convert_ambiguous_to_IUPAC=FALSE
	xml=NULL
	xlsfn=xlsfn
	return_charsdf=TRUE	
	'
	
	# Get the list of OTUs you actually want to use
	if (is.null(OTUs) == TRUE)
		{
		tmpOTUs_df = readWorksheetFromFile(file=xlsfn, sheet="OTUs", startRow=15, startCol=1, header=TRUE)
		tmpOTUs_df$use[isblank_TF(tmpOTUs_df$use)] = ""
		OTUs_to_keep_TF = tmpOTUs_df$use != "no"
		OTUs = tmpOTUs_df$OTUs[OTUs_to_keep_TF]
		tmpOTUs_df = tmpOTUs_df[OTUs_to_keep_TF,]
		}
	
	# Default is use="yes"
	seqs_df$use[isblank_TF(seqs_df$use)] = "yes"
	
	
	# Initialize
	chars_wLT_2_states_df_list = list()
	count_per_charstate_per_char = NULL
	
	# Return the characters data frame for morphology, if desired
	if (return_charsdf == TRUE)
		{
		cnum = 0
		charsdf_list = list()
		} # END if (return_charsdf == TRUE)
	
	# Subset to just the sequence datasets you intend to use
	seqs_df = seqs_df[seqs_df$use == "yes", ]
	
	
	# Get unique NEXUS filenames
	
	# Standardize the filenames, to include the directory etc.
	# Convert blank directories to ""
	# If there is only one slash, fix that
	seqs_df$dir[isblank_TF(seqs_df$dir)] = ""
	dirs = seqs_df$dir
	fns = seqs_df$filename
	for (fnum in 1:length(fns))
		{
		if (isblank_TF(dirs[fnum]) == TRUE)
			{
			next()
			}
		
		# If the prefix in the filename is redundant with 
		# the one in dir, remove it from fn
		prefix = addslash(get_fn_prefix(fns[fnum]))
		tmpdir = addslash(dirs[fnum])
		if (grepl(pattern=prefix, x=tmpdir) == TRUE)
			{
			fns[fnum] = gsub(pattern=prefix, replacement="", x=fns[fnum])
			}

		# Paste on the directory, unless the filename is "traits"
		if (fns[fnum] != "traits")
			{
			fns[fnum] = slashslash(paste(seqs_df$dir[fnum], fns[fnum], sep="/"))
			} # END if (fns[fnum] != "traits")
		
		# If only 1 slash, at the beginning, fix that
		hits = gregexpr(pattern="/", text=fns[fnum])[[1]]
		if ( (length(hits) == 1) && (hits[[1]] == 1) )
			{
			fns[fnum] = gsub(pattern="/", replacement="", x=fns[fnum])
			} # END if (length(hits) == 1)
		} # END for (fnum in 1:fns)
	
	uniq_fns = unique(fns)
	uniq_fn_rownums = match(uniq_fns, fns)
	datatypes = seqs_df$type[uniq_fn_rownums]
	uniq_fn_rownums
	
	# Keep track of the number of (used) morphology characters in 
	# each morphology DATASET (not partition)
	morphLengths = NULL
	
	# Keep track of the morphology sections
	# (used for making the morphology models and likelihood tags)
	morphList = NULL
	
	matrices_stats = list()
	
	# First, see if there are two versions of morph, ordered and unordered
	# Revise the rows to include accordingly...
	unique_rows_revised = NULL
	for (i in 1:length(uniq_fn_rownums))
		{
		current_rownum = uniq_fn_rownums[i]
		if (toupper(datatypes[i]) != toupper("morph"))
			{
			unique_rows_revised = c(unique_rows_revised, current_rownum)
			next()	# Skip to next item in the loop
			}
	
		if (toupper(datatypes[i]) == toupper("morph"))
			{
			dataset_name = seqs_df$datasetName[current_rownum]
			TF = seqs_df$datasetName == dataset_name
			rownums_matching = (1:nrow(seqs_df))[TF]
			for (j in 1:length(rownums_matching))
				{
				current_rownum = rownums_matching[j]
				unique_rows_revised = c(unique_rows_revised, current_rownum)
				}
			}
		} # END for (i in 1:length(uniq_fn_rownums))
	uniq_fn_rownums = unique_rows_revised
	uniq_fn_rownums
	uniq_fns = fns[uniq_fn_rownums]
	
	
	# There should be one dataset_name for each file, 
	# or for e.g. unordered vs. ordered morph from a single file
	dataset_names = seqs_df$datasetName[uniq_fn_rownums]
	
	# For DNA or AAs, choose which function to use to read DNA/AAs
	# (read_nexus_data2 is SLOOOOOW for large files)
	seqs_read_functions = seqs_df$seqs_read_function[uniq_fn_rownums]
	
	if (is.null(seqs_df$file_types) == FALSE)
		{
		file_types = seqs_df$file_types[uniq_fn_rownums]
		} else {
		file_types = rep("", times=length(uniq_fn_rownums))
		} # END if (is.null(seqs_df$file_types) == FALSE)
	
	
	# Number of characters in each raw dataset (each file)
	dataset_lengths = rep(0, length(dataset_names))
	dataset_ntaxa = rep(0, length(dataset_names))
	
	# Completeness statistics
# 	meanQs
# 	maxQs
# 	minQs
# 	pctData
# 	nMajData

	meanQs = rep(0, length(dataset_names))
	maxQs = rep(0, length(dataset_names))
	minQs = rep(0, length(dataset_names))
	pctData = rep(0, length(dataset_names))
	nMajData = rep(0, length(dataset_names))

	
	# filteredAlignment_names can be whole files, or 
	# subsets of those files (e.g. partitions)
	filteredAlignment_names = seqs_df$filteredAlignmentName[uniq_fn_rownums]
	datatypes = seqs_df$type[uniq_fn_rownums]
	geneTreeNames = seqs_df$geneTreeName[uniq_fn_rownums]
	seqNames_to_cut = seqs_df$seqNames_to_cut[uniq_fn_rownums]
	
	# Old versions had an "order_type" column; this is deprecated, but 
	# being filled in from the "model" column
	models = seqs_df$model[uniq_fn_rownums]
	# default is unordered
	order_types = rep("unordered", times=length(models))
	order_types[models == "Mk_ord"] = "ordered"
	
	# Added 2016-09-24:
	# Different kinds of morphological models (MkA for e.g. asymmetric)
	morph_transition_rates = seqs_df$morph_transition_rates[uniq_fn_rownums]
	baseFreqs = seqs_df$baseFreqs[uniq_fn_rownums]
	
	# Ascertainment bias models
	ascertainment = seqs_df$ascertainment[uniq_fn_rownums]
	
	# Defaults
	# (irrelevant for DNA)
	morph_transition_rates[isblank_TF(morph_transition_rates)] = "equal"
	baseFreqs[isblank_TF(baseFreqs)] = "equal"
	ascertainment[isblank_TF(ascertainment)] = "Mkv"
	
	
	# Maximum number of unobservable patterns to list for ascertainment bias correction: 
	# mostly to keep your R script and/or BEAST from crashing, when they try to 
	# write out, or read in, the potentially BILLIONS of unobservable site patterns
	# that there are for, say, a 5-state unordered character with 100 taxa
	max_num_patterns = seqs_df$max_num_patterns[uniq_fn_rownums]
	max_num_patterns[isblank_TF(max_num_patterns)] = 2500
	max_num_patterns = as.numeric(max_num_patterns)
	
	
	# Assumed nstates -- if e.g. you want to run a 2-state correction on a 3-state
	# character (mostly because this might be what MrBayes has actually been doing).
	assumed_nstates = seqs_df$assumed_nstates[uniq_fn_rownums]
	assumed_nstates[isblank_TF(assumed_nstates)] = ""
	
	#print("assumed_nstates2:")
	#print(assumed_nstates)
	
	
	charlists = seqs_df$list[uniq_fn_rownums]
	startchars = seqs_df$startchar[uniq_fn_rownums]
	endchars = seqs_df$endchar[uniq_fn_rownums]
	bynums = seqs_df$by[uniq_fn_rownums]
	numGammaCats = seqs_df$gammaNum[uniq_fn_rownums]
	clockModel_names = seqs_df$clockModel_name[uniq_fn_rownums]
	clockModel_relRates = seqs_df$clockModel_relRates[uniq_fn_rownums]
	gammaShape_suffixes = seqs_df$gammaShape_suffix[uniq_fn_rownums]
	
	if (length(dataset_names) != length(uniq_fns) )
		{
		txt = paste("\n\nERROR in parse_datasets(): In worksheet 'data', the number of unique filenames (", length(uniq_fns),")\ndoes not match the number of unique dataset names (", length(dataset_names), ").\n\n", sep="")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if (length(dataset_names) != length(fns) )

	
	# Are there continuous traits in the dataset?
	traits_TF = "continuous" %in% datatypes
	traits_TF
	
	# Initialize XMLs for holding the datasets
	data_XML = list()
	prior_XML = list()
	misc_XML = list()
	sitemodels_XML = list()
	state_XML = list()
	
	if (traits_TF == TRUE)
		{
		likes_XMLs = list(bl(), xmlCommentNode(" idrefs for additional likelihoods, e.g. for continuous traits "))
		operator_XMLs = list(bl(), xmlCommentNode(" Additional operators, e.g. for continuous traits "))
		tracelog_XMLs = list(bl(), xmlCommentNode(" Additional logs, e.g. for continuous traits "))
		screenlog_XMLs = list(bl(), xmlCommentNode(" Additional logs, e.g. for continuous traits "))
		} else {
		likes_XMLs = list(bl(), xmlCommentNode(" idrefs for additional likelihoods "))
		operator_XMLs = list(bl(), xmlCommentNode(" Additional operators "))
		tracelog_XMLs = list(bl(), xmlCommentNode(" Additional logs "))
		screenlog_XMLs = list(bl(), xmlCommentNode(" Additional logs "))
		} # END if (traits_TF == TRUE)
	
	old_dataset_name = "not_used_yet"	# to avoid reloading morphology data
	
	####################################################
	####################################################
	# LOOP THROUGH DATA FILES
	####################################################
	####################################################
	for (i in 1:length(uniq_fns))
		{
		cat("\nparse_datasets() is starting read of file ", i, "/", length(uniq_fns), ": '", uniq_fns[i], "'...\n", sep="")
		
		dataset_name = dataset_names[i]
		fn = uniq_fns[i]
		file_type = tolower(file_types[i])
		datatype = toupper(datatypes[i])
		numGammaCat = numGammaCats[i]
		clockModel_name = clockModel_names[i]
		clockModel_relRate = clockModel_relRates[i]
		gammaShape_suffix = gammaShape_suffixes[i]
		
		# Error check: Check if the file exists
		if (file.exists(fn) == FALSE)
			{
			stoptxt = paste0("STOP ERROR in parse_datasets(): File #", i, "/", length(uniq_fns), " could not be found. This is the file that was sought:\n\n", fn, "\n\n")
			
			cat("\n\n")
			cat(stoptxt)
			cat("\n\n")
			
			stop(stoptxt)
			} # END Error check: Check if the file exists

		
		
		# If file_type is blank, try to guess the file type from the filename
		if (isblank_TF(file_type) == TRUE)
			{
			if (endsWith(x=fn, suffix="fas") == TRUE)
				{
				file_type = "fasta"
				}
			if (endsWith(x=fn, suffix="fasta") == TRUE)
				{
				file_type = "fasta"
				}
			if (endsWith(x=fn, suffix="nexus") == TRUE)
				{
				file_type = "nexus"
				}
			if (endsWith(x=fn, suffix="nex") == TRUE)
				{
				file_type = "nexus"
				}
			if (endsWith(x=fn, suffix="nxs") == TRUE)
				{
				file_type = "nexus"
				}
			# If it's STILL blank, say "nexus"
			if (isblank_TF(file_type) == TRUE)
				{
				file_type = "nexus"
				} # END if (isblank_TF(file_type) == TRUE)

			txt = paste0("\nNOTE: for data file '", fn, "' BEASTmasteR is guessing that the file_type is '", file_type, "'.\n")
			cat(txt)
			} # END if (isblank_TF(file_type) == TRUE)
		
		
		if (file_type == "fasta")
			{
			newfn = paste0(fn, ".nex")
			txt = paste0("\nNOTE: for data file '", fn, "' the file_type is '", file_type, "', so BEASTmasteR will copy and convert to a NEXUS file with name '", newfn, "'.\n")
			cat(txt)
			
			# Save out to NEXUS, make that the filename
			out_nexus_aln_fn = convert_NEXUS_to_FASTA(fasta_fn=fn, out_nexus_aln_fn=newfn)
			file_type = "nexus"
			fn = out_nexus_aln_fn
			} # END if (file_type == "fasta")
		
		
		if (any(toupper(c("DNA", "AA", "morph", "continuous")) == datatype) == FALSE)
			{
			txt = paste0("\n\nERROR in parse_dataset(): datatype '", datatype, "' not recognized.\n\n")
			stop(txt)
			}
	
		if (datatype == toupper("DNA"))
			{
			if (seqs_read_functions[i] == "read_nex_phyloch_to_list")
				{
				cmdstr = paste(dataset_name, " = read_nex_phyloch_to_list(fn=fn)")
				} else {
				cmdstr = paste(dataset_name, " = read_nexus_data2(file=fn, convert_ambiguous_to_IUPAC=convert_ambiguous_to_IUPAC, printall=printall, check_ambig_chars=FALSE)")
				} # END if (seqs_read_functions[i] == "read_nex_phyloch_to_list")

			# Loads to whatever variable given by "dataset_name"
			eval(parse(text=cmdstr))
			
			# Store the length (number of characters) of the raw dataset
			cmdstr = paste0("dataset_lengths[i] = length(", dataset_name, "[[1]])")
			eval(parse(text=cmdstr))

			# Store the number of taxa of the raw dataset
			cmdstr = paste0("dataset_ntaxa[i] = length(", dataset_name, ")")
			eval(parse(text=cmdstr))
			
			
			# Statistics on the DNA matrix
			# typically seqs_DNA
			cmdstr = paste0("DNAstats = DNA_matrix_stats(charslist=", dataset_name, ", charsdf=NULL)")
			eval(parse(text=cmdstr))
			matrices_stats_tmp = DNAstats
			matrices_stats = c(matrices_stats, list(matrices_stats_tmp))
						
			
			# Get the completeness stats by OTU
			# numQs (number of nodata)
			
			# For this dataset, find numQs for the dataset
			cmdstr = paste0("numQs_for_each_OTU = sapply(X=", dataset_name, ", FUN=count_ambig_DNA)")
			eval(parse(text=cmdstr))
			
			meanQs[i] = mean(numQs_for_each_OTU, na.rm=TRUE)
			maxQs[i] = max(numQs_for_each_OTU, na.rm=TRUE)
			minQs[i] = min(numQs_for_each_OTU, na.rm=TRUE)
			pctData_for_each_OTU = (1-(numQs_for_each_OTU / dataset_lengths[i])) * 100
			pctData[i] = mean(pctData_for_each_OTU, na.rm=TRUE)
			nMajData[i] = sum( pctData_for_each_OTU >= 50)
			} # END if (datatype == "DNA")
		
		
		if (datatype == toupper("AA"))
			{
			if (seqs_read_functions[i] == "read_nex_phyloch_to_list")
				{
				cmdstr = paste(dataset_name, " = read_nex_phyloch_to_list(fn=fn)")
				} else {
				cmdstr = paste(dataset_name, " = read_nexus_data2(file=fn, convert_ambiguous_to_IUPAC=convert_ambiguous_to_IUPAC, printall=printall, check_ambig_chars=FALSE)")
				} # END if (seqs_read_functions[i] == "read_nex_phyloch_to_list")

			# Loads to whatever variable given by "dataset_name"
			eval(parse(text=cmdstr))

			# Store the length (number of characters) of the raw dataset
			cmdstr = paste0("dataset_lengths[i] = length(", dataset_name, "[[1]])")
			eval(parse(text=cmdstr))

			# Store the number of taxa of the raw dataset
			cmdstr = paste0("dataset_ntaxa[i] = length(", dataset_name, ")")
			eval(parse(text=cmdstr))

			# Statistics on the AA matrix
			cmdstr = paste0("AAstats = AA_matrix_stats(charslist=", dataset_name, ", charsdf=NULL)")
			eval(parse(text=cmdstr))

			matrices_stats_tmp = AAstats
			matrices_stats = c(matrices_stats, list(matrices_stats_tmp))
			


			# Get the completeness stats
			# numQs (number of nodata)
			
			# For this dataset, find numQs for the dataset
			cmdstr = paste0("numQs_for_each_OTU = sapply(X=", dataset_name, ", FUN=count_ambig_nonDNA)")
			eval(parse(text=cmdstr))
			
			meanQs[i] = mean(numQs_for_each_OTU, na.rm=TRUE)
			maxQs[i] = max(numQs_for_each_OTU, na.rm=TRUE)
			minQs[i] = min(numQs_for_each_OTU, na.rm=TRUE)
			pctData_for_each_OTU = (1-(numQs_for_each_OTU / dataset_lengths[i])) * 100
			pctData[i] = mean(pctData_for_each_OTU, na.rm=TRUE)
			nMajData[i] = sum( pctData_for_each_OTU >= 50)
			} # END if (datatype == "AA")
		
		# Continuous morphology data
		# Assumes tab-delimited file
		if (datatype == toupper("continuous"))
			{
			defaults='
			fn="/drives/Dropbox/_njm/__packages/BEASTmasteR_permahelp/examples/ex_basic_venerid_morphDNA_v4/SABD_tipsVary_noOutg_v1/traits_table.tab"
			dataset_name="traits_table"
			'
			
			# If the filename is "traits", extract from the Excel file,
			if (fn == "traits")
				{
				# Read from Excel file, "traits" tab
				#print(xlsfn)
				#print(fn)
				
				if (is.null(xlsfn))
					{
					errortxt = "\n\nSTOP ERROR in parse_datasets(): seqs_df, derived from the Excel 'data' worksheet, indicates that you have a continuous dataset in the 'traits' worksheet. However, this requires that the Excel filename be given in xlsfn= in parse_datasets(). You have xlsfn=NULL or not specified, that is the cause of the error.\n\n"
					cat(errortxt)
					stop(errortxt)
					} # END if (xlsfn == NULL)
				
				tmpdf = readWorksheetFromFile(file=xlsfn, sheet=fn, startRow=15, startCol=1, header=TRUE)
				} else {
				# Read from tab-delimited text file
				# Row names in column 1
				#print(xlsfn)
				#print(fn)
				cmdstr=paste0("tmpdf = read.table(file=fn, header=TRUE, sep='\t', strip.white=TRUE, stringsAsFactors=FALSE)")
				# Loads to 'tmpdf'
				eval(parse(text=cmdstr))
				}
			
			# Cut out unused OTUs
			# Subset this table based on the master OTUs
			TF = tmpdf$OTUs %in% OTUs
			tmpdf = tmpdf[TF,]
			
			# Apply the OTU names as row.names
			row.names(tmpdf) = tmpdf$OTUs

			
			# Store the continuous trait data in dataset and in e.g. "traits"
			dataset = tmpdf[,2:ncol(tmpdf)]
			cmdstr = paste0(dataset_name, " = dataset")
			# Loads to whatever variable given by "dataset_name"
			eval(parse(text=cmdstr))
			
			# Store the length (number of characters) of the raw dataset
			cmdstr = paste0("dataset_lengths[i] = ncol(", dataset_name, ")")
			eval(parse(text=cmdstr))

			# Store the number of taxa of the raw dataset
			cmdstr = paste0("dataset_ntaxa[i] = nrow(", dataset_name, ")")
			eval(parse(text=cmdstr))
			
			cmdstr = paste0("TF1 = as.matrix(", dataset_name, ") == '?'")
			eval(parse(text=cmdstr))
			
			cmdstr = paste0("TF2 = as.matrix(", dataset_name, ") == '-'")
			eval(parse(text=cmdstr))

			cmdstr = paste0("TF3 = isblank_TF(as.matrix(", dataset_name, "))")
			eval(parse(text=cmdstr))
			Q_TF = (TF1 + TF2 + TF3) > 0
			meanQs = mean(Q_TF, na.rm=TRUE)
			maxQs = max(Q_TF, na.rm=TRUE)
			minQs = min(Q_TF, na.rm=TRUE)
			
			numQs_by_taxon = rowSums(Q_TF)
			pctData_by_taxon = (1 - (numQs_by_taxon / ncol(Q_TF))) * 100
			pctData = mean(pctData_by_taxon, na.rm=TRUE)
			nMajData = sum(pctData_by_taxon >= 50)
			
			# Each different partition of continuous data 
			# should come from a different file!
			seqs_df_continuous_TF = seqs_df$type == "continuous"
			
			if ( (length(unique(seqs_df$datasetName[seqs_df_continuous_TF])) != 1) || (length(unique(seqs_df$filename[seqs_df_continuous_TF])) != 1) )
				{
				error_txt = "STOP ERROR in parse_datasets(). Each continuous dataset ('alignment') in the 'data' worksheet needs its own filename & dataset. Instead, you have repeats..."
				cat("\n\n")
				cat(error_txt)
				cat("\n\n")
				
				cat("Printing repeats in seqs_df:\n\n")
				print(seqs_df[seqs_df_continuous_TF,])
				cat("\n\n")
				stop(error_txt)
				} # END if (more than one use of continuous filename or datasetName
			} # END if (datatype == "continuous")


		if (datatype == toupper("morph"))
			{
			# Avoid re-loading twice, if two morph rows are adjacent
			if (dataset_name == old_dataset_name)
				{
				# Dataset is already in eval(parse(text=dataset_name))
				old_dataset_name = dataset_name
				dataset_previously_loaded = TRUE
				# Here, we are just using the old dataset name

				# Get statistics for the morphology matrix
				cmdstr = paste0("charslist = ", dataset_name)
				eval(parse(text=cmdstr))
				matrices_stats_tmp = morphology_matrix_stats(charslist, charsdf=NULL)
				matrices_stats = c(matrices_stats, list(matrices_stats_tmp))

				} else {
				cmdstr = paste0(dataset_name, " = read_nexus_data2(file=fn, convert_ambiguous_to_IUPAC=FALSE, check_ambig_chars=TRUE, printall=printall, convert_ambiguous_to=NULL)")
				old_dataset_name = dataset_name
				dataset_previously_loaded = FALSE
				# Loads to whatever variable given by "dataset_name"
				eval(parse(text=cmdstr))

				# Get statistics for the morphology matrix
				cmdstr = paste0("charslist = ", dataset_name)
				eval(parse(text=cmdstr))
				matrices_stats_tmp = morphology_matrix_stats(charslist, charsdf=NULL)
				matrices_stats = c(matrices_stats, list(matrices_stats_tmp))
				
				if (return_charsdf == TRUE)
					{
					charsdf_list[[(cnum=cnum+1)]] = as.data.frame(charslist, stringsAsFactors=FALSE)
					}
				} # END if (dataset_name == old_dataset_name)

			# Store the length (number of characters) of the raw dataset
			cmdstr = paste0("dataset_lengths[i] = length(", dataset_name, "[[1]])")
			eval(parse(text=cmdstr))

			# Store the number of taxa of the raw dataset
			cmdstr = paste0("dataset_ntaxa[i] = length(", dataset_name, ")")
			eval(parse(text=cmdstr))

			# Get the completeness stats
			# numQs (number of nodata)
			
			# For this dataset, find numQs for the dataset
			cmdstr = paste0("numQs_for_each_OTU = sapply(X=", dataset_name, ", FUN=count_ambig_nonDNA)")
			eval(parse(text=cmdstr))
			
			meanQs[i] = mean(numQs_for_each_OTU, na.rm=TRUE)
			maxQs[i] = max(numQs_for_each_OTU, na.rm=TRUE)
			minQs[i] = min(numQs_for_each_OTU, na.rm=TRUE)
			pctData_for_each_OTU = (1-(numQs_for_each_OTU / dataset_lengths[i])) * 100
			pctData[i] = mean(pctData_for_each_OTU, na.rm=TRUE)
			nMajData[i] = sum( pctData_for_each_OTU >= 50)
			} # END if (datatype == "morph")


		# Which characters (split "list" cell on comma)
		charlist = charlists[i]
		startchar = startchars[i]
		endchar = endchars[i]
		bynum = bynums[i]
		
		# Fill in the charnums/charlist 
		# (the list of which
		if (isblank_TF(charlist) == TRUE)
			{
			# Corrections
			if (isblank_TF(startchar) == TRUE)
				{
				startchar = 1
				}
			if (isblank_TF(bynum) == TRUE)
				{
				bynum = 1
				}
			if (isblank_TF(endchar) == TRUE)
				{
				endchar = dataset_lengths[i]
				} # END if (isblank_TF(endchar) == TRUE)
		
			# If the charnums are missing, use all characters
# 			if ( isblank_TF(startchar) || isblank_TF(endchar) || isblank_TF(bynum) )
# 				{
# 				charnums = 1:nrow(morph_df2_corrected)
# 				} else {
				# Otherwise, use the specified characters
			charnums = seq(from=startchar, to=endchar, by=bynum)
# 				} # END check for start/end/by
			} else {
			# If the list is user-specified
			charnums = as.numeric(strsplit(x=charlist, split=",")[[1]])
			} # END if (isblank_TF(charnums))


		

		#########################################
		# Write seqs to XML
		#########################################
		# Datasets now exist in:
		dataset_names
		# seqs_DNA
		# seqs_morph
		# Store the dataset in a common variable, "dataset"
		cmdstr = paste0("dataset = ", dataset_names[i])
		eval(parse(text=cmdstr))
		#print(dataset)
		# now stored in "dataset"


		#########################################
		# Write DNA to XML
		#########################################
		if (datatypes[i] == "DNA")
			{
			totalcount = 4
			dataType_txt = "nucleotide"
			txtseqs = toupper(collapse_seqs(dataset))
			OTU_names = names(txtseqs)
			
			cat("\nLength of DNA seqs list before cutting unused OTUs: ", length(txtseqs))
			cat("\n")
			
			# Subset DNA sequences, if OTUs list specified
			# But, add back in any loci that are for gene trees, rather
			# than species tree
			if ((is.null(OTUs) == FALSE) && (isblank_TF(geneTreeNames[i]) == TRUE) )
				{
				keepTF = OTU_names %in% OTUs
				OTU_names = OTU_names[keepTF]
				
				txtseqs = as.list(txtseqs[keepTF])
				} # END if (is.null(OTUs) == FALSE)

			cat("Length of DNA seqs list after cutting unused OTUs: ", length(txtseqs))
			cat("\n")

			# For a geneTree, you might want to cut some specific sequences from an
			# alignment. These are listed in "seqNames_to_cut"
			if (isblank_TF(geneTreeNames[i]) == FALSE)
				{
				if (isblank_TF(seqNames_to_cut[i]) == FALSE)
					{
					# Cut the listed sequences
					words = trim(strsplit(seqNames_to_cut[i], split=",")[[1]])
					keepTF = rep(TRUE, length(OTU_names))
					for (w in 1:length(words))
						{
						cutTF = grepl(pattern=words[w], x=OTU_names)
						keepTF[cutTF==TRUE] = FALSE
						}
					txtseqs = as.list(txtseqs[keepTF])
					OTU_names = OTU_names[keepTF]
					} # END if (isblank_TF(seqNames_to_cut[i]) == FALSE)
				cat("Length of DNA seqs list after cutting unused seqNames_to_cut: ", length(txtseqs))
				cat("\n(applies only if worksheet 'data', column 'geneTreeName' is not blank)\n\n")
				} # END if (isblank_TF(geneTreeNames[i]) == FALSE)

			
			seq_XML = make_seq_XML(txtseqs[1], names(txtseqs[1]), dataset_name, totalcount)
			seqs_XML = mapply(FUN=make_seq_XML, txtseq=txtseqs, OTU_name=OTU_names, MoreArgs=list(dataset_name=dataset_name, totalcount=totalcount))
		
			# Comment
			txt = paste0(" ", datatypes[i], " sequence. Name: ", dataset_name, ". From: ", uniq_fns[i], " ")
			XML_comment = xmlCommentNode(txt)
	
			# Data node
			seqs_node = xmlNode(name="data", attrs=list(id=dataset_name, dataType=dataType_txt), .children=seqs_XML)
			seqs_node = list(bl(), XML_comment, seqs_node)
		
			# Store in big data XML list
			data_XML = c(data_XML, seqs_node)
			} # END if (datatypes[i] == "DNA")


		#########################################
		# Write AAs to XML
		#########################################
		if (datatypes[i] == "AA")
			{
			totalcount = 20
			dataType_txt = "aminoacid"
			txtseqs = toupper(collapse_seqs(dataset))
			OTU_names = names(txtseqs)
			
			# Subset AA sequences, if OTUs list specified
			if (is.null(OTUs) == FALSE)
				{
				TF = OTU_names %in% OTUs
				OTU_names = OTU_names[TF]
				txtseqs = as.list(txtseqs[TF])
				} # END if (is.null(OTUs) == FALSE)
			
			seq_XML = make_seq_XML(txtseqs[1], names(txtseqs[1]), dataset_name, totalcount)
			seqs_XML = mapply(FUN=make_seq_XML, txtseq=txtseqs, OTU_name=OTU_names, MoreArgs=list(dataset_name=dataset_name, totalcount=totalcount))
		
			# Comment
			txt = paste0(" ", datatypes[i], " sequence. Name: ", dataset_name, ". From: ", uniq_fns[i], " ")
			XML_comment = xmlCommentNode(txt)
	
			# Data node
			seqs_node = xmlNode(name="data", attrs=list(id=dataset_name, dataType=dataType_txt), .children=seqs_XML)
			seqs_node = list(bl(), XML_comment, seqs_node)

			# Store in big data XML list
			data_XML = c(data_XML, seqs_node)
			} # END if (datatypes[i] == "AA")

	
		#########################################
		# Write continuous data to XML (and lots of other stuff)
		#########################################
		if (datatypes[i] == "continuous")
			{
			OTU_names = row.names(dataset)
			
			# Subset continuous sequences, if OTUs list specified
			if (is.null(OTUs) == FALSE)
				{
				TF = OTU_names %in% OTUs
				OTU_names = OTU_names[TF]
				dataset = dataset[TF,]
				} # END if (is.null(OTUs) == FALSE)
			
			# Each different partition of continuous data 
			# should come from a different file!
			
			# Subset columns by charnums
			#print(dataset)
			#print(charnums)
			dataset = dataset[, charnums]
			
			# Header on data XML
			txt = paste0(" Continuous trait data from '", dataset_name, "'-'", filteredAlignment_names[i], "', encoded 1 trait at a time, which assumes independence. ")
			txt_XML = xmlCommentNode(txt)
			txt_XMLs = list(bl(), bl(), txt_XML)
			data_XML = c(data_XML, txt_XMLs)
			
			# Set up elements/models used for ALL continuous characters
			# in this partition
			dataset_name = dataset_names[i]
			userDataType_id = paste0("userDataType_for_", dataset_name)
			userDataType_XML = xmlNode(name="userDataType", attrs=list(id=userDataType_id, spec="beast.evolution.datatype.ContinuousDataType"))
			userDataType_XML
			txt = " The continuous datatype; same element used by all continuous characters in the partition (assumes independence) "
			userDataType_XMLs = list(bl(), xmlCommentNode(txt), userDataType_XML)
			misc_XML = c(misc_XML, userDataType_XMLs)
			
			
	
			# Set up the siteModel used for ALL continuous characters
			# in this partition
			numGammaCat = numGammaCats[i]
			if (isblank_TF(numGammaCat))
				{
				numGammaCat = 1
				}
			if (numGammaCat=="no")
				{
				numGammaCat = 1
				}
			if (numGammaCat!="no")
				{
				numGammaCat = 1
				}
			
			
			# Make a gamma shape parameter for among-site rate variation for
			# continuous traits
			
			# stateNode for gammaShape of trait dataset
			trait_gammaShape_id = paste0("gammaShape_for_", dataset_name)
			trait_gammaShape_idref = paste0("@", trait_gammaShape_id)
			trait_gammaShape_XML = xmlNode(name="parameter", value=1, attrs=list(id=trait_gammaShape_id, name="stateNode"))
			state_XML = c(state_XML, list(trait_gammaShape_XML))
			
			# Same prior distribution for gammaShape parameter 
			# for traits in a trait dataset
			# Add to misc
			trait_gammaShape_Exponential_id = paste0("exponential_on_gammaShape_for_", dataset_name)
			trait_gammaShape_Exponential_idref = paste0("@", trait_gammaShape_Exponential_id)
			gammaExponential_param_id = paste0("exponential_mean_of_1_for_", dataset_name)
			gammaExponential_param_XML = xmlNode(name="parameter", value=1, attrs=list(id=gammaExponential_param_id, lower="0.0", upper="0.0", name="mean"))
			trait_gammaShape_Exponential_XML = xmlNode(name="Exponential", attrs=list(id=trait_gammaShape_Exponential_id, name="distr"), .children=list(gammaExponential_param_XML))
			trait_gammaShape_Exponential_XMLs = list(bl(), xmlCommentNode(" This is an exponential distribution used as the prior on each gamma shape parameter "), trait_gammaShape_Exponential_XML)
			misc_XML = c(misc_XML, trait_gammaShape_Exponential_XMLs)
			misc_XML
			
			# Actual prior on that gammaShape parameter for traits in a trait dataset
			gammaShape_prior_id = paste0("gammaShape_for_", dataset_name, "_prior")
			gammaShape_prior_XML = xmlNode(name="prior", attrs=list(id=gammaShape_prior_id, x=trait_gammaShape_idref, distr=trait_gammaShape_Exponential_idref, name="distribution"))
			txt = paste0(" Prior on gamma shape parameter for traits in ", dataset_name, " ")
			gammaShape_prior_XMLs = list(bl(), xmlCommentNode(txt), gammaShape_prior_XML)
			prior_XML = c(prior_XML, gammaShape_prior_XMLs)
			
			trait_gammaShape_Exponential_XMLlog = xmlNode(name="log", attrs=list(idref=trait_gammaShape_id))
			trait_gammaShape_prior_XMLlog = xmlNode(name="log", attrs=list(idref=gammaShape_prior_id))
			tracelog_XMLs = c(tracelog_XMLs, list(trait_gammaShape_Exponential_XMLlog, trait_gammaShape_prior_XMLlog))
			
			# Operator on the gammaShape parameter for traits in a trait dataset
			gammaShape_operator_id = paste0("gammaShapeScaler_for_", dataset_name)
			
			# Set operator to 0 weight, if there is only one gamma category
			if (numGammaCat > 1)
				{
				gamma_operator_weight = 0.1
				} else {
				gamma_operator_weight = 0.0
				}
			gammaShape_operator_XML = xmlNode(name="operator", attrs=list(id=gammaShape_operator_id, parameter=trait_gammaShape_idref, scaleFactor="0.5", spec="ScaleOperator", weight=gamma_operator_weight))
			operator_XMLs = c(operator_XMLs, list(gammaShape_operator_XML))
			

			
			# Write out each column separately
			# (independent characters)
			# (yes, this is kind of ridiculous, but better than
			#  dealing with a huge covariance matrix when we 
			#  are assuming the characters are independent)
			for (colnum in 1:ncol(dataset))
				{
				traitName = names(dataset)[colnum]
				traitmap_id = paste0(dataset_name, "_", traitName, "_traitmap")
				traitmap_idref = paste0("@", traitmap_id)
				parameter_id = paste0(dataset_name, "_", traitName, "_parameter")

				trait_data_vector = apply(X=cbind(row.names(dataset), dataset[,colnum]), MARGIN=1, FUN=paste, sep="", collapse="=")
				trait_data_txt = paste(trait_data_vector, sep="", collapse=",\n")
				#cat(trait_data_txt)
				


				###############################################
				# Separate substitution model, sitemodel, and precisionMatrix (rate)
				# for each continuous trait
				###############################################

				# Add the precision matrix (NOT same for every character in the partition,
				# as the different traits are getting different rates)
				# it's just 1x1)
				# Add it to stateNode
				precisionMatrix_id = paste0("precisionMatrix_for_", traitName, "_in_", dataset_name)
				precisionMatrix_idref = paste0("@", precisionMatrix_id)
				startvalue = 1
				precisionMatrix_XML = xmlNode(name="parameter", startvalue, attrs=list(id=precisionMatrix_id, dimension="1", minordimension="1", name="stateNode"))
				txt = paste0(" Using the same precisionMatrix (1x1) for all continuous characters from dataset '", dataset_name, "'. ")
				precisionMatrix_XMLs = list(bl(), xmlCommentNode(txt), precisionMatrix_XML)
				state_XML = c(state_XML, precisionMatrix_XMLs)
				precisionMatrix_XMLlog = xmlNode(name="log", attrs=list(idref=precisionMatrix_id))
				tracelog_XMLs = c(tracelog_XMLs, list(precisionMatrix_XMLlog))

				# Add the trait node values
				# Add it to stateNode
				# number of nodes = ntips + (ntips-1)
				numnodes = nrow(dataset) + nrow(dataset) - 1
				meanval = mean(as.numeric(dataset[,colnum]), na.rm=TRUE)
				traitNodeVals_id = paste0("traitNodeVals_for_", traitName, "_in_", dataset_name)
				traitNodeVals_idref = paste0("@", traitNodeVals_id)
				startvalue = 1
				traitNodeVals_XML = xmlNode(name="parameter", meanval, attrs=list(id=traitNodeVals_id, dimension=numnodes, name="stateNode"))
				txt = paste0(" Trait values at each node for '", traitName, "' in '", dataset_name, "'.\n(Tree has ", nrow(dataset), " tips, ", nrow(dataset)-1, " internal nodes, ", numnodes, " total nodes.) ")
				traitNodeVals_XMLs = list(bl(), xmlCommentNode(txt), traitNodeVals_XML)
				state_XML = c(state_XML, traitNodeVals_XMLs)
				#tracelog_XMLs = c(tracelog_XMLs, precisionMatrix_XMLs)





				
				# Settings for interpretation of the traitmap values
				# NJM: These are important!  
				# This error is easy to get and very frustrating:
				##########################################################################
				# java.lang.NullPointerException
				# 	at beast.evolution.tree.TreeTraitMap.initAndValidate(Unknown Source)
				# 	at beast.util.XMLParser.initPlugins(Unknown Source)
				# 	at beast.util.XMLParser.parse(Unknown Source)
				# 	at beast.util.XMLParser.parseFile(Unknown Source)
				# 	at beast.app.BeastMCMC.parseArgs(Unknown Source)
				# 	at beast.app.beastapp.BeastMain.main(Unknown Source)
				# 
				# Error 110 parsing the xml input file
				# 
				# validate and intialize error: null
				# 
				# Error detected about here:
				#   <beast>
				#       <run id='mcmc' spec='MCMC'>
				#           <distribution id='posterior' spec='util.CompoundDistribution'>
				#               <distribution id='prior' spec='util.CompoundDistribution'>
				#                   <prior id='Prior_on_traits_table_trait1_rootTrait' name='distribution'>
				#                       <x id='traits_table_trait1_rootTrait' spec='beast.evolution.tree.RootTrait'>
				#                           <traitmap id='traitmap1' spec='beast.evolution.tree.TreeTraitMap'>
				##########################################################################
				# 
				# The following things cause this error:
				# 
				# 1. Not having 2 trees initialized outside and inside state
				# (See Remco (RRB)'s help:
				# https://groups.google.com/forum/#!topic/beast-users/7KFzWd4PVeQ
				#
				# 2. Also, <traitmap> has to be *inside* <x> *inside* the root prior
				#    inside <prior> inside <posterior> inside <run>.  No references!
				#
				# 3. Similarly, <siteModel> for the continuous trait has to be
				#    inside the likelihood, and the likelhood has to be inside
				#    <likelihood> inside <posterior> inside <run>
				# 
				# 4. Finally, initByMean="true" inside <traitmap>
				# 
				# 5. Also, no extra taxa inside "list_of_OTUs"
				#
				# This message was brought to you by NJM's entire weekend, 3/14/15-3/15/15
				# (happy Pi Day!)
				##########################################################################
				#
				# Settings for traitmap
				
				# Make randomize lower 1/10000 of lowest value, or 0
				randomizelower_val = 1/10000 * min(as.numeric(dataset[,colnum]), na.rm=TRUE)

				# Make randomize upper 10 times highest value, or 0
				randomizeupper_val = 10 * max(as.numeric(dataset[,colnum]), na.rm=TRUE)

				# Make the jitter value 1/100 of max of data values
				jitterval = (0.01 * max(as.numeric(dataset[,colnum]), na.rm=TRUE))
				
				#######################################
				# Make traitmap XML
				#######################################
				traitmap_XML = xmlNode(name="traitmap", trait_data_txt, attrs=list(id=traitmap_id, traitName=traitName, initByMean="true", randomizelower=randomizelower_val, randomizeupper=randomizeupper_val, jitter=jitterval, parameter=traitNodeVals_idref, tree="@shared_tree", spec="beast.evolution.tree.TreeTraitMap") )
				comment_txt = paste0(" Trait data for dataset: ", dataset_name, ", trait #", colnum, " (#", charnums[colnum], " of original matrix); trait name: ", traitName, ". ")
				comment_XML = xmlCommentNode(comment_txt)
				traitmap_XMLs = list(bl(), comment_XML, traitmap_XML)
				
				# Store in big data XML list - NO
				# Store inside prior inside priors inside run
				#data_XML = c(data_XML, traitmap_XMLs)
				
				###########################################
				# Might as well do the prior here, as well
				# Reference the traitmap(= trait data) which is 
				# stored earlier in the XML
				###########################################
				# Make "x" XML
				# Store traitmap here
				rootTrait_id = paste0(dataset_name, "_", traitName, "_rootTrait")
				
				# Old, doesn't work
				# x_XML = xmlNode(name="x", attrs=list(id=rootTrait_id, traitmap=traitmap_idref, spec="beast.evolution.tree.RootTrait") )
				# Works, put traitmap inside prior
				x_XML = xmlNode(name="x", attrs=list(id=rootTrait_id, spec="beast.evolution.tree.RootTrait"), .children=traitmap_XMLs )
				
				comment_txt = paste0(" rootTrait for dataset: ", dataset_name, ", trait #", colnum, " (#", charnums[colnum], " of original matrix); trait name: ", traitName, ". ")
				comment_XML = xmlCommentNode(comment_txt)
				x_XMLs = list(comment_XML, x_XML)
				
				########################################
				# Prior on root trait value
				# (just have a very broad prior, i.e. normal with 
				#  mean = uncorrected mean and
				#  sd = 10 * mean )
				meanval = mean(as.numeric(dataset[,colnum]), na.rm=TRUE)
				sdval = 10 * meanval
				Normal_id = paste0("Normal_on_", rootTrait_id)
				Normal_XML = xmlNode(name="Normal", attrs=list(id=Normal_id, name="distr", mean=meanval, sigma=sdval))
				
				rootTrait_prior_id = paste0("Prior_on_", rootTrait_id)
				rootTrait_prior_XML = xmlNode(name="prior", attrs=list(id=rootTrait_prior_id, name="distribution"), .children=c(list(Normal_XML), x_XMLs))
				txt1 = paste0(" Prior distribution on", comment_txt)
				txt2 = paste0(" By default, we are putting a very flat prior on the root value of the trait,\nbut if you had a reason (e.g. known ancestral value) you could change this manually. ")
				rootTrait_prior_XMLs = list(bl(), xmlCommentNode(txt1), xmlCommentNode(txt2), rootTrait_prior_XML)
				prior_XML = c(prior_XML, rootTrait_prior_XMLs)


						

				# Set up the substModel used for ALL continuous characters
				# in this partition
				substModel_id = paste0("substModel_for_", traitName, "_in_", dataset_name)
				substModel_idref = paste0("@", substModel_id)
				substModel_XML = xmlNode(name="substModel", attrs=list(id=substModel_id, precisionMatrix=precisionMatrix_idref, spec="beast.continuous.MultivariateDiffusionModel"))
		
		
				continuous_siteModel_id = paste0("siteModel_for_", traitName, "_in_", dataset_name)
				continuous_siteModel_idref = paste0("@", continuous_siteModel_id)
				continuous_siteModel_XML = xmlNode(name="siteModel", attrs=list(id=continuous_siteModel_id, shape=trait_gammaShape_idref, gammaCategoryCount=numGammaCat, spec="SiteModel"), .children=list(substModel_XML))
				continuous_siteModel_txt = paste0(" A different siteModel element is used for each continuous character in the trait '", traitName, "' in the dataset '", dataset_name, "' (assumes independence) ")
				continuous_siteModel_XMLs = list(bl(), xmlCommentNode(continuous_siteModel_txt), continuous_siteModel_XML)

				#sitemodels_XML = c(sitemodels_XML, continuous_siteModel_XMLs)




				
				
				#######################################################
				# ADD LIKELIHOODS FOR EACH TRAIT
				# AN OPERATOR FOR PRECISION MATRIX
				# etc...
				#######################################################
				# scaleByTime *MIGHT* mean this:
				# https://code.google.com/p/beast-classic/source/browse/trunk/src/beast/continuous/AbstractMultivariateTraitLikelihood.java?spec=svn81&r=81
				# if (scaleByTime) {
				# 	length /= treeLength;
				# }
				# Going with  the settings here, basically:
				# ITS_continuous_phylogeography_1trait.xml :
				# 
				# <distribution id="locationtreeLikelihood.location"
				#  reciprocalRates="true" scaleByTime="false"
				#  spec="beast.continuous.SampledMultivariateTraitLikelihood"
				#  reportAsMultivariate="false" traitParameter="@location.location"
				#  tree="@Tree.t:tree" useTreeLength="false">
				# 
				#userDataType_id = paste0(dataset_name, "_", traitName, "_userDataType")
				userDataType_XML = xmlNode(name="userDataType", attrs=list(idref=userDataType_id))
				dataref_id = paste0(dataset_name, "_", traitName, "_dataref")
				dataref_for_trait_XML = xmlNode(name="data", attrs=list(id=dataref_id, spec="AlignmentFromTraitMap", traitMap=traitmap_idref), .children=list(userDataType_XML))
				
				# Don't use siteModel reference in likelihood, use actual siteModel
				#siteModelref_XML = xmlNode(name="siteModel", attrs=list(idref=continuous_siteModel_id))
				
				#branchRateModel_id = 
				# Specifying a strict clock here; but we want to use the
				# overall clock model!!
				# <branchRateModel clock.rate="@clockRate.c:location"
				# id="StrictClock.c:location"
				# spec="beast.evolution.branchratemodel.StrictClockModel"/>
				branchRateModel_XML = xmlNode(name="branchRateModel", attrs=list(idref=clockModel_names[i]))
				
				tmp_children = c(list(dataref_for_trait_XML), continuous_siteModel_XMLs, list(branchRateModel_XML))
				trait_likelihood_id = paste0(dataset_name, "_", traitName, "_likelihood")
				trait_likelihood_idref = paste0("@", trait_likelihood_id)
				

				######################################################
				# BEGIN guessing options of SampledMultivariateTraitLikelihood
				######################################################				
				# Reciprocal Rates, documentation says:
				# https://code.google.com/p/beast-classic/source/browse/trunk/src/beast/continuous/AbstractMultivariateTraitLikelihood.java?spec=svn81&r=81
				# if (reciprocalRates) {
                # true means:
                # length /= rateModel.getRateForBranch(node); // branch rate scales as precision (inv-time)
				# NJM Translation: false means branch lengths divided  
				# by Brownian rate(s) (variance)
				# false means:
				# length *= rateModel.getRateForBranch(node); // branch rate scales as variance (time)
				# NJM Translation: false means branch lengths multiplied 
				# by Brownian rate(s) (variance)
				#
				# if (scaleByTime) { length /= treeLength; }
				# NJM translation: divide branches so that tree is length 1
				#
				# useTreeLength
				# NJM best guess: something about recalculating tree length;
				# going with default, which is false
				######################################################
				# END guessing options of SampledMultivariateTraitLikelihood
				######################################################				
				
				trait_likelihood_XML = xmlNode(name="distribution", attrs=list(id=trait_likelihood_id, traitParameter=traitNodeVals_idref, tree="@shared_tree", reciprocalRates="false", scaleByTime="false", reportAsMultivariate="false", useTreeLength="false", spec="beast.continuous.SampledMultivariateTraitLikelihood"), .children=tmp_children)
				
				# Save the actual likelihood calculator in misc
				txt = paste0(" Calculation of likelihood on continuous ", traitName, " in dataset ", dataset_name, " ")
				trait_likelihood_XMLs = list(bl(), xmlCommentNode(txt), trait_likelihood_XML)
				misc_XML = c(misc_XML, trait_likelihood_XMLs)
				
				# Refer to the likelihood calculators in the Likelihoods section
				like_XML = xmlNode(name="distribution", attrs=list(idref=trait_likelihood_id))
				likes_XMLs = c(likes_XMLs, list(like_XML))
				
				# Refer to the likelihood calculators in the log
				log_XML = xmlNode(name="log", attrs=list(idref=trait_likelihood_id))
				tracelog_XMLs = c(tracelog_XMLs, list(log_XML))



				#######################################################
				# AN OPERATOR FOR PRECISION MATRIX FOR THIS TRAIT 
				#######################################################
				# Gibbs operator
				scaleMatrix_parameter_id = paste0("scaleMatrix_parameter_for_precisionPrior_for_", precisionMatrix_id)
				scaleMatrix_parameter_XML = xmlNode(name="parameter", attrs=list(id=scaleMatrix_parameter_id, value=1.0, dimension=1, minordimension=1, name="scaleMatrix"))
				precisionPrior_for_GibbsOperator_id = paste0("precisionPrior_for_GibbsOperator_for_", precisionMatrix_id)
				precisionPrior_for_GibbsOperator_XML = xmlNode(name="prior", attrs=list(id=precisionPrior_for_GibbsOperator_id, arg=precisionMatrix_idref, df="1.0", spec="beast.math.distributions.WishartDistribution"), .children=list(scaleMatrix_parameter_XML))


				precisionMatrix_GibbsOperator_id = paste0("precisionGibbsOperator_for_", precisionMatrix_id)
				precisionMatrix_GibbsOperator_XML = xmlNode(name="operator", attrs=list(id=precisionMatrix_GibbsOperator_id, parameter=precisionMatrix_idref, likelihood=trait_likelihood_idref, traitmap=traitmap_idref, tree="@shared_tree", weight="15.0", spec="PrecisionMatrixGibbsOperator"), .children=list(precisionPrior_for_GibbsOperator_XML))


				traitGibbsOperator_id = paste0("traitGibbsOperator_for_", precisionMatrix_id)
				traitGibbsOperator_XML = xmlNode(name="operator", attrs=list(id=traitGibbsOperator_id, precisionMatrix=precisionMatrix_idref, likelihood=trait_likelihood_idref, traitmap=traitmap_idref, tree="@shared_tree", weight="50.0", spec="TraitGibbsOperator"))

				RootTraitRandowWalkOperator_id = paste0("RootTraitRandowWalkOperator_for_", precisionMatrix_id)
				RootTraitRandowWalkOperator_XML = xmlNode(name="operator", attrs=list(id=RootTraitRandowWalkOperator_id, parameter=traitNodeVals_idref, traitmap=traitmap_idref, weight="5.0", windowSize="10.0", spec="RootTraitRandowWalkOperator"))
				
				# Add to operators
				operator_XMLs = c(operator_XMLs, list(precisionMatrix_GibbsOperator_XML, traitGibbsOperator_XML, RootTraitRandowWalkOperator_XML))
				
				} # END for (colnum in 1:ncol(dataset))
			} # END if (datatypes[i] == "continuous")

	
		#########################################
		# Write discrete morphology to XML (and lots of other stuff)
		#########################################
		if (toupper(datatypes[i]) == toupper("morph"))
			{
		
			# Comment
			txt = paste0(" ", datatypes[i], " ", order_types[i], " characters. Dataset: ", dataset_name, ". From: ", uniq_fns[i], " ")
			XML_comment = xmlCommentNode(txt)
	
			# Data node
			#seqs_node = xmlNode(name="data", attrs(id=dataset_name), .children=seqs_XML)
			seqs_node = NULL
	
	
				
			#print("Printing dataset:")
			#print(dataset)
			#print(dim(dataset))
			#print(length(dataset))
			#print(class(dataset))
			#print("Printing dataset done.")
			#dataset2 = dataset
			
			# Convert dataset to data.frame
			morph_df = as.data.frame(x=dataset, col.names=names(dataset), stringsAsFactors=FALSE)
			morph_df
			
			#print(dim(morph_df))
			#print(OTUs)
			
			# Subset morphological OTUs, if OTUs list specified
			if (is.null(OTUs) == FALSE)
				{
				# dataset (named list of sequences)
				OTU_names = names(dataset)
				TF = OTU_names %in% OTUs
				OTU_names = OTU_names[TF]
				dataset = as.list(dataset[TF])
				#print(TF)
				
				# morph_df (data.frame with OTUs as columns)
				OTU_names = names(morph_df)
				TF = OTU_names %in% OTUs
				OTU_names = OTU_names[TF]
				morph_df = morph_df[,TF]
				#print(TF)
				} # END if (is.null(OTUs) == FALSE)
			
		
			# Check for characters that have missing states, and 
			# correct
			#print(morph_df)
			#print(morph_df[353,])
			# NO: print(morph_df[,353])
		
			res = get_numstates_per_char(morph_df, ambig_to_remove=c("\\(", "\\)", " ", ","), return_missing_chars="correct", printall=printall, count_autapomorphies=TRUE, sort_ambigs=TRUE)
			missing_charstates_list = res$missing_charstates_list
			allQs_TF = res$allQs_TF
			morph_df2_corrected = res$nexdf
			count_per_charstate_per_char = res$count_per_charstate_per_char
			char_is_autapomorphic = res$char_is_autapomorphic
			autapomorphies_desc_df = res$autapomorphies_desc_df
			matrices_stats[[i]]$autapomorphies_desc_df = autapomorphies_desc_df
			
			# Keep track of the number of (used) morphology characters in 
			# each morphology DATASET (not partition)
			number_of_morphological_characters = nrow(morph_df2_corrected)
			if (dataset_previously_loaded == FALSE)
				{
				morphLengths = c(morphLengths, number_of_morphological_characters)
				}
			
			# Convert any "1 state" characters to 2 states
			# (probably people should correct/revise their data)
			numstates_morph_list = res$numstates_morph_list
			numstates_morph_list[numstates_morph_list == 1] = 2
			
			if (dataset_previously_loaded == FALSE)
				{
				# Text results for the XML comments
				txt1 = paste0(" Results of get_numstates_per_char() on dataset: ", dataset_name, " ")
				txt2 = paste0(" Filename: ", fn)
				txt3 = paste0(" There were ", length(missing_charstates_list), " characters missing an implied state. ")
				txt4 = paste0(" They were characters numbered: ", paste(missing_charstates_list, sep="", collapse=", "), " ")
				txt5 = paste0(" E.g., states 013 were found, but not state 2. This character would get re-coded as 012. ")
				txt6 = paste0(" (Beast complains if states are missing; this might be OK with ordered characters,  ")
				txt7 = paste0("  but I have not explored it...NJM 2014-10-22)                                        ")
				txt8 = paste0(" Listing ORIGINAL data in PHYLIP format: ")
				orig_matrix = dtf_to_phylip(char_dtf=t(morph_df), txtTF=TRUE)
				orig_matrix = paste("\n", orig_matrix, "\n", sep="")
			
				# Have to convert "-" to "?" for XML display
				orig_matrix = gsub(pattern="-", replacement="?", x=orig_matrix)
			
			
				txt9 = paste0(" Listing CORRECTED data in PHYLIP format: ")
				morph_df2_corrected_wParens = t(add_parens_to_multistate_chars(char_dtf=t(morph_df2_corrected)))
				corrected_matrix = dtf_to_phylip(char_dtf=t(morph_df2_corrected_wParens), txtTF=TRUE)
				corrected_matrix = paste("\n", corrected_matrix, "\n", sep="")
				# Have to convert "-" to "?" for XML display
				corrected_matrix = gsub(pattern="-", replacement="?", x=corrected_matrix)
		
				txt10 = paste0(" Number of characters that were all-unknown (?, dashes, etc.): ", sum(allQs_TF), " ")
				txt11 = paste0(" They were: ", paste((1:nrow(morph_df))[allQs_TF], sep=", "), " ")
				txt12 = paste0(" Any all-unknown characters will be cut from the Beast XML ")
		
				data_XML = c(data_XML, list(bl()))
				data_XML = c(data_XML, list(xmlCommentNode(txt1)))
				data_XML = c(data_XML, list(xmlCommentNode(txt2)))
				data_XML = c(data_XML, list(bl()))
				data_XML = c(data_XML, list(xmlCommentNode(txt3)))
				data_XML = c(data_XML, list(xmlCommentNode(txt4)))
				data_XML = c(data_XML, list(xmlCommentNode(txt5)))
				data_XML = c(data_XML, list(xmlCommentNode(txt6)))
				data_XML = c(data_XML, list(xmlCommentNode(txt7)))
				data_XML = c(data_XML, list(bl()))
				data_XML = c(data_XML, list(xmlCommentNode(txt10)))
				data_XML = c(data_XML, list(xmlCommentNode(txt11)))
				data_XML = c(data_XML, list(xmlCommentNode(txt12)))
				data_XML = c(data_XML, list(bl()))
				data_XML = c(data_XML, list(xmlCommentNode(txt8)))
				data_XML = c(data_XML, list(xmlCommentNode(orig_matrix)))
				data_XML = c(data_XML, list(bl()))
				data_XML = c(data_XML, list(xmlCommentNode(txt9)))
				data_XML = c(data_XML, list(xmlCommentNode(corrected_matrix)))
				data_XML = c(data_XML, list(bl()))
				data_XML
				} # END if (dataset_previously_loaded == FALSE)
		

			# Now that you have the charnums, subset both the 
			# characters, and numstates_morph_list
			rownums = 1:nrow(morph_df2_corrected)
			TFs = rep(FALSE, times=nrow(morph_df2_corrected))
			TFs[charnums] = TRUE
			TFs[allQs_TF] = FALSE # Remove all "?", "-"
			charnums_subset = rownums[TFs]
			numstates_morph_list_subset = numstates_morph_list[charnums_subset]
			morph_df2_corrected_subset = morph_df2_corrected[charnums_subset,]
		
			charnums_txt = paste(charnums_subset, sep="", collapse=",")
			charnums_txt
			txt1 = paste0(" Character numbers in original matrix: ")
			txt2 = paste0(" ", charnums_txt, " ")
			data_XML = c(data_XML, list(bl()))	
			data_XML = c(data_XML, list(xmlCommentNode(txt1)))
			data_XML = c(data_XML, list(xmlCommentNode(txt2)))
	
			
			# Store morphology characters in the XML
			# Unordered *or* ordered characters
			nexd6 = morph_df2_corrected_subset
			
			if (dataset_previously_loaded == FALSE)
				{
				result2 = make_BEAST2_userDataTypes(numstates_morph_list_subset, nexd6=nexd6, dataset_name=dataset_name, ordering=order_types[i], morph_transition_rates=morph_transition_rates[i], baseFreqs=baseFreqs[i], numGammaCat=numGammaCat, clockModel_name=clockModel_name, clockModel_relRate=clockModel_relRate, gammaShape_suffix=gammaShape_suffix, add_morphList=add_morphList, printall=printall)
				
				# NJM edit for 2 different morphology matrices
				if ( !is.null(morphList) && (dim(morphList) == dim(result2$morphList)) && all(morphList == result2$morphList) )
					{
					morphList = result2$morphList
					} else {
					morphList = rbind(morphList, result2$morphList) 
					}
				chars_wLT_2_states_df_list = c(chars_wLT_2_states_df_list, list(result2$chars_wLT_2_states_df))
				nexd7 = result2$nexd5
				} else {
				result2 = make_BEAST2_userDataTypes(numstates_morph_list_subset, nexd6=nexd6, dataset_name=dataset_name, ordering=order_types[i], morph_transition_rates=morph_transition_rates[i], baseFreqs=baseFreqs[i], numGammaCat=numGammaCat, clockModel_name=clockModel_name, clockModel_relRate=clockModel_relRate, gammaShape_suffix=gammaShape_suffix, add_morphList=add_morphList, printall=printall)
				morphList = rbind(morphList, result2$morphList)
				chars_wLT_2_states_df_list = c(chars_wLT_2_states_df_list, list(result2$chars_wLT_2_states_df))
				nexd7 = result2$nexd5
				}
		
			# Add the datatypes to the XML
			userDataType_nodes_list = result2$userDataType_nodes_list
			data_XML = c(data_XML, userDataType_nodes_list)
		
			# Add the morphology data to the XML
			morph_alignment_nodes_list = write_BEAST2_morphology_characters(nexd7, numstates_morph_list=numstates_morph_list_subset, dataset_name=dataset_name, ordering=order_types[i], morph_transition_rates=morph_transition_rates[i], baseFreqs=baseFreqs[i], ascertainment=ascertainment[i], max_num_patterns=max_num_patterns[i], assumed_nstates=assumed_nstates[i])
		
			seqs_node = c(list(bl()), list(XML_comment), morph_alignment_nodes_list)

			# Store in big data XML list
			data_XML = c(data_XML, seqs_node)
			} # END if (datatypes[i] == "morph")

		cat("\nparse_datasets() has finished reading of file ", i, "/", length(uniq_fns), ".\n", sep="")

		} # END for (i in 1:length(uniq_fns))
		#########################################
		# END writing out characters
		#########################################
	
	dataset_lengths_df = cbind(dataset_names, dataset_lengths, dataset_ntaxa, nMajData, round(pctData,1), round(minQs,2), round(meanQs,2), round(maxQs,2))
	dataset_lengths_df = as.data.frame(dataset_lengths_df, stringsAsFactors=FALSE)
	names(dataset_lengths_df) = c("dataset_names", "dataset_lengths", "dataset_ntaxa", "nMajData", "pctData", "minQs", "meanQs", "maxQs")
	
	# Output nodes, or xml
	if (is.null(xml))
		{
		data_XMLs = NULL
		# If desired, add the number of (used!) characters
		# of each morphology DATASET (not partition)
		data_XMLs$data_XML = data_XML
		data_XMLs$prior_XML = prior_XML
		data_XMLs$misc_XML = misc_XML
		data_XMLs$state_XML = state_XML
		data_XMLs$sitemodels_XML = sitemodels_XML
		data_XMLs$likes_XMLs = likes_XMLs
		data_XMLs$operator_XMLs = operator_XMLs
		data_XMLs$tracelog_XMLs = tracelog_XMLs
		data_XMLs$screenlog_XMLs = screenlog_XMLs
		data_XMLs$chars_wLT_2_states_df_list = chars_wLT_2_states_df_list
		data_XMLs$dataset_lengths_df = dataset_lengths_df
		if (add_morphLength == TRUE)
			{
			data_XMLs$morphLengths = morphLengths
			}
		if (add_morphLength == TRUE)
			{
			data_XMLs$morphList = morphList
			}
		if (!is.null(matrices_stats))
			{
			if (is.null(count_per_charstate_per_char) == FALSE)
				{
				#matrices_stats[[length(matrices_stats)]]$autapomorphies_desc_df = autapomorphies_desc_df
				} # END if (is.null(count_per_charstate_per_char) == FALSE)
			
			data_XMLs$matrices_stats = matrices_stats
			} # END if (!is.null(matrices_stats))
		if (return_charsdf == TRUE)
			{
			data_XMLs$charsdf_list = charsdf_list
			} # END if (return_charsdf == TRUE)
		return(data_XMLs)
		} else {
		# Store the sequences
		xml$sequences = c(xml$sequences, data_XML)
		xml$prior = c(xml$prior, prior_XML)
		xml$misc = c(xml$misc, misc_XML)
		xml$state = c(xml$state, state_XML)
		xml$likes = c(xml$likes, likes_XMLs)
		xml$operators = c(xml$sitemodels, operator_XMLs)
		xml$tracelog = c(xml$sitemodels, tracelog_XMLs)
		xml$screenlog = c(xml$sitemodels, screenlog_XMLs)
		# If desired, add the number of (used!) characters
		# of each morphology DATASET (not partition)
		xml$chars_wLT_2_states_df_list = chars_wLT_2_states_df_list
		if (add_morphLength == TRUE)
			{
			xml$morphLengths = morphLengths
			xml$morphList = morphList
			}
		xml$dataset_lengths_df = dataset_lengths_df
			
			
		morphList
		if (!is.null(matrices_stats))
			{
			if (is.null(count_per_charstate_per_char) == FALSE)
				{
				matrices_stats[[length(matrices_stats)]]$autapomorphies_desc_df = autapomorphies_desc_df
				} # END if (is.null(count_per_charstate_per_char) == FALSE)

			xml$matrices_stats = c(xml$matrices_stats, list(matrices_stats))
			}		
		if (return_charsdf == TRUE)
			{
			xml$charsdf_list = c(xml$charsdf_list, charsdf_list)
			} # END if (return_charsdf == TRUE)

		return(xml)
		}
	
	stop("parse_datasets(): shouldn't get here")
	} # END parse datasets



# XML SECTION 4: Partitions 
# (filters applied to the sequence alignments produce partitions)
make_partitions_XML <- function(seqs_df, xml=NULL, add_partitionLength=TRUE, dataset_lengths_df=NULL)
	{
	nonmorph_TF = toupper(seqs_df$type) != toupper("morph")
	noncontinuous_TF = toupper(seqs_df$type) != toupper("continuous")
	keepTF = (nonmorph_TF + noncontinuous_TF) == 2
	
	# If morphology-only, return null
	if (sum(keepTF) == 0)
		{
		if (is.null(xml))
			{
			return(NULL)
			} else {
			return(xml)
			} # END if (is.null(xml))
		} # END if (length(seqs_df$type != "morph") == 0)

	# If no dataset_lengths, get them from the xml
	if (is.null(dataset_lengths_df))
		{
		dataset_lengths_df = xml$dataset_lengths_df
		}

	# For non-morphology datasets, extract partitions from the alignments
	# Remove morphology from seqs
	seqs_df_original = seqs_df
	seqs_df = seqs_df[keepTF, ]
	
	# Gather the partition lengths, if desired
	partitionLengths = NULL
	partitionNames = NULL
	
	# Make filters
	uniq_partitions = unique(seqs_df$filteredAlignmentName)
	uniq_partitions
	
	# Get the corresponding dataset names
	rows_TF = seqs_df$filteredAlignmentName %in% uniq_partitions
	uniq_partitions_datasetNames = seqs_df$datasetName[rows_TF]
	
	data_filtered_to_partitions_XML = NULL

	for (i in 1:length(uniq_partitions))
		{
		partitionName = uniq_partitions[i]
		TF = seqs_df$filteredAlignmentName == partitionName
		rownums = (1:nrow(seqs_df))[TF]
		partition_df = seqs_df[rownums, ]
		
		# Zero out the length counter for this partition
		partitionLengths_tmp = 0
		
		# Loop through the partition subsections
		for (j in 1:length(rownums))
			{
			rownum = rownums[j]
		
			# Get the source dataset
			datasetName = seqs_df$datasetName[rownum]
			
			# Get the partition information
			# Which characters (split "list" cell on comma)
			charlist = seqs_df$list[rownum]
			startchar = seqs_df$startchar[rownum]
			endchar = seqs_df$endchar[rownum]
			bynum = seqs_df$by[rownum]
	
			if (isblank_TF(bynum) == TRUE)
				{
				bynum = 1
				}

			# Fill in the charnums/charlist 
			# (the list of which
			if (isblank_TF(charlist) == TRUE)
				{
				# Corrections
				if (isblank_TF(startchar) == TRUE)
					{
					startchar = 1
					}
				if (isblank_TF(bynum) == TRUE)
					{
					bynum = 1
					}
				if (isblank_TF(endchar) == TRUE)
					{
					#print(dataset_lengths_df$dataset_names)
					#print(datasetName)
					
					TF = dataset_lengths_df$dataset_names == datasetName
					#print(TF)
					endchar = as.numeric(dataset_lengths_df$dataset_lengths[TF])
					
					# Error check if this didn't work
					if (isblank_TF(endchar))
						{
						txt = paste0("STOP ERROR in make_partitions_XML(). For dataset '", datasetName, "' (partition '", partitionName, "'), rownum=", rownum, ", column 'endchar' was not specified, nor was make_partitions_XML() given the data.frame 'dataset_lengths_df' with columns 'dataset_names' and 'dataset_lengths'. Fix one of these so that make_partitions_XML() can know how long the dataset is.")
						cat("\n\n")
						cat(txt)
						cat("\n\n")
						stop(txt)
						}
					} # END if (isblank_TF(endchar) == TRUE)
		
				# If the charnums are missing, use all characters
	# 			if ( isblank_TF(startchar) || isblank_TF(endchar) || isblank_TF(bynum) )
	# 				{
	# 				charnums = 1:nrow(morph_df2_corrected)
	# 				} else {
					# Otherwise, use the specified characters
				charnums = seq(from=startchar, to=endchar, by=bynum)
	# 				} # END check for start/end/by
				} else {
				# If the list is user-specified
				charnums = as.numeric(strsplit(x=charlist, split=",")[[1]])
				} # END if (isblank_TF(charnums))
			

	
			if (isblank_TF(charlist) == TRUE)
				{
				# If the charnums are missing, use all characters
				if ( isblank_TF(startchar) || isblank_TF(endchar) || isblank_TF(bynum) )
					{
					# Use all characters
					filtertxt = paste0("1", "-", nrow(morph_df2_corrected))
					charlist_txt_for_length = seq(1, nrow(morph_df2_corrected), by=1)
					} else {
					# Otherwise, use the specified characters
					filtertxt = paste(startchar, "-", endchar, "\\", bynum, sep="")
					charlist_txt_for_length = seq(startchar, endchar, by=bynum)
					#filtertxt = gsub(pattern="\\", replacement="\\", x=filtertxt)
					} # END check for start/end/by
				} else {
				# Use a comma-delimited list of the characters
				charnums_txt = strsplit(x=charlist, split=",")[[1]]
				charlist_txt_for_length = charnums_txt
				filtertxt = paste(charnums_txt, sep=",", collapse="")
				} # END if (isblank_TF(charnums))
			
			# Gather the lengths of this subset of this partition
			partitionLengths_tmp = partitionLengths_tmp + length(charlist_txt_for_length)
			
			# Store, or add to, the filtertxt
			if (j == 1)
				{
				stored_filtertxt = filtertxt
				} else {
				stored_filtertxt = paste(stored_filtertxt, ",", filtertxt, sep="")
				}
	
			} # END for (j in 1:length(rownums))

		# Gather the partition lengths, if desired
		partitionLengths = c(partitionLengths, partitionLengths_tmp)
		partitionNames = c(partitionNames, partitionName)
		
		txt = paste(" Partition ", partitionName, " filter text: ", stored_filtertxt, sep="")
		cat(" \n", txt)
	
		partitionNode_comment = xmlCommentNode(txt)
	
		# Make the partition nodes
		# data reference
		dataref_XML = xmlNode(name="data", attrs=list(idref=datasetName) )
		partitionNode = xmlNode(name="alignment", attrs=list(id=partitionName, filter=stored_filtertxt, spec="FilteredAlignment"), .children=list(dataref_XML) )

	
		partitionNode_XMLlist = list(bl(), partitionNode_comment, partitionNode)
		data_filtered_to_partitions_XML = c(data_filtered_to_partitions_XML, partitionNode_XMLlist)
		} # END for (i in 1:length(uniq_partitions))
	data_filtered_to_partitions_XML
	
	partitionLengths_df = NULL
	partitionLengths_df$partitionNames = partitionNames
	partitionLengths_df$partitionLengths = partitionLengths
	partitionLengths_df = as.data.frame(partitionLengths_df, stringsAsFactors=FALSE)
		
	if (is.null(xml))
		{
		# Just return the data_filtered_to_partitions_XML 
		
		# Gather the partition lengths, if desired
		if (add_partitionLength == TRUE)
			{
			data_filtered_to_partitions_XML$partitionLengths_df = partitionLengths_df
			}
		return(data_filtered_to_partitions_XML)
		} else {
		xml$partitions = c(xml$partitions, data_filtered_to_partitions_XML)
		
		# Gather the partition lengths, if desired
		if (add_partitionLength == TRUE)
			{
			xml$partitionLengths_df = partitionLengths_df
			}
		
		return(xml)
		} # END if (is.null(xml))
	
	stop("ERROR in make_partitions_XML(): shouldn't get here")
	} # END make_partitions_XML <- function(seqs_df, xml=NULL)


# For StarBeast2 analyses:
# Prune the sequences for just those that match the 
# "taxonsets" worksheet, if this is a StarBeast2 run	
prune_seqs_based_on_taxonsets <- function(data_XML, taxonsets_df, StarBeast2_TF)
	{
	# Exit, if not a StarBeast2 run
	if (StarBeast2_TF == FALSE)
		{
		return(data_XML)
		}
	
	# Filter taxonsets_df to the "use" column
	taxonsets_df = readWorksheetFromFile(xlsfn, sheet="taxonsets", startRow=15)
	taxonsets_df$use[isblank_TF(taxonsets_df$use)] = "yes"
	keepTF = (taxonsets_df$use != "no")
	taxonsets_df = taxonsets_df[keepTF, ]
	
	for (i in 1:length(data_XML))
		{
		namevals = names(data_XML[[i]])
		if ( (length(namevals) > 0) && (namevals[1] == "sequence") )
			{
			# Go through the sequences
			seqs = data_XML[[i]]
			
			for (j in 1:length(seqs))
				{
				taxon_val = xmlAttrs(seqs[[j]])$taxon
				TF = taxon_val %in% taxonsets_df$specimenString
				if (TF == FALSE)
					{
					txt = paste0(" Note: Sequence from specimen '", taxon_val, "' was removed as specimen was not found in taxonsets_df$specimenString. ")
					seqs[[j]] = xmlCommentNode(txt)
					} # END if (TF == FALSE)
				} # END for (j in 1:length(seqs))
			
			# Provide edited sequences
			data_XML[[i]] = seqs
			} # END if ( (length(namevals) > 0) && (namevals[1] == "sequence") )
		} # END for (i in 1:length(data_XML)
	return(data_XML)
	}



get_numQs_from_XMLtag <- function(xmltag, isDNA=TRUE)
	{
	# Return NA for non-sequence xmltags
	TF = is_XML_comment(xmltag)
	if (TF == TRUE)
		{
		return(NA)
		}
	# E.g., exclude sequences named "text"
# 	if (names(xmltag) != "sequence")
# 		{
# 		return(NA)
# 		}
	# OK, you have sequence, so extract it
	# Get the DNA (or other sequence)
	chars = xmlAttrs(xmltag)$value
	if (isDNA == TRUE)
		{
		numQs = count_ambig_DNA_in_string(chars)
		}
	if (isDNA == FALSE)
		{
		numQs = count_ambig_nonDNA_in_string(chars)
		}
	return(numQs)
	} # END get_numQs_from_XMLtag <- function(xmltag, isDNA=TRUE)

get_lengths_from_XMLtag <- function(xmltag)
	{
	# Return NA for non-sequence xmltags
	TF = is_XML_comment(xmltag)
	if (TF == TRUE)
		{
		return(NA)
		}
	# E.g., exclude sequences named "text"
# 	if (names(xmltag) != "sequence")
# 		{
# 		return(NA)
# 		}
	# OK, you have sequence, so extract it
	# Get the DNA (or other sequence)
	chars = xmlAttrs(xmltag)$value
	lengths = nchar(chars)

	return(lengths)
	} # END get_lengths_from_XMLtag <- function(xmltag, isDNA=TRUE)


is_XML_comment <- function(xmltag)
	{
	classes = class(xmltag)
	if (("XMLCommentNode" %in% classes) == TRUE)
		{
		return(TRUE)
		} else {
		return(FALSE)
		}
	} # END is_XML_comment <- function(xmltag)


prune_seqs_based_on_completeness <- function(data_XML, min_numseqs=0, min_fraction_to_keep_sequence=0.5)
	{
	XML_comments_TF = sapply(X=data_XML, FUN=is_XML_comment)
	alignment_nums = (1:length(data_XML))[XML_comments_TF==FALSE]
	
#	for (i in 1:4)
	for (i in 1:length(alignment_nums))
		{
		tmp_alignment = data_XML[[alignment_nums[i]]]
		namevals = names(tmp_alignment)
		
		# Get the dataType
		tmp_datatype = xmlAttrs(tmp_alignment)$dataType
		if (isblank_TF(tmp_datatype) == TRUE)
			{
			# Then, it is morphology or something, don't cut
			cat("\nNote: in dataXML, xmlAttrs(tmp_alignment)$dataType returns blank. This is probably morphology; prune_seqs_based_on_completeness() is tested only on DNA at the moment. Skipping...\n")
			next()
			#return(data_XML)
			}
		if (tmp_datatype == "nucleotide")
			{
			isDNA = TRUE
			} else {
			isDNA = FALSE
			} # END if (tmp_datatype == "nucleotide")
		
		
		# If this alignment (XML tag) has more than 0 names, and 
		# if some of them are sequence, then this is an alignment
		seqs_TF = "sequence" %in% namevals
		xmltag_is_seqs_TF = namevals == "sequence"
		dataset_ntaxa = sum(xmltag_is_seqs_TF)

		if ( (length(namevals) > 0) && (sum(seqs_TF) > 0) )
			{
			# Get the lengths of each DNA sequence, and the number
			# of nondata in each
			seqnums = (1:length(namevals))[xmltag_is_seqs_TF]
			children = xmlChildren(tmp_alignment)
			seqs = children
			#sapply(X=seqs, FUN=is_XML_comment)
			lengths = unlist(sapply(X=seqs, FUN=get_lengths_from_XMLtag))
			if (isDNA == TRUE)
				{
				numQ_vals = unlist(sapply(X=seqs, FUN=get_numQs_from_XMLtag, isDNA=TRUE))
				}
			if (isDNA == FALSE)
				{
				numQ_vals = unlist(sapply(X=seqs, FUN=get_numQs_from_XMLtag, isDNA=FALSE))
				}
			names(numQ_vals) = NULL
			names(lengths) = NULL
			
			fract_data_vals = 1 - (numQ_vals / lengths)
			fract_data_vals
			
			# Go through the sequences
			seqnums_eval_TF = is.na(fract_data_vals) == FALSE
			seqnums2 = (1:length(fract_data_vals))[seqnums_eval_TF]
			
			# Cut this locus entirely, if there are not enough sequences with enough data
			num_seqs_that_are_complete_enough = sum( fract_data_vals[seqnums_eval_TF] >= min_fraction_to_keep_sequence )
			if (num_seqs_that_are_complete_enough < min_numseqs)
				{
				cut_all_seqs = TRUE
				} else {
				cut_all_seqs = FALSE
				}
				
			# For sequences that are being cut, replace the XML
			for (j in 1:length(seqnums2))
				{
				# Get the DNA (or other sequence)
				seqid = xmlAttrs(seqs[[seqnums2[j]]])$id
				taxon_val = xmlAttrs(seqs[[seqnums2[j]]])$taxon
				if (is.null(taxon_val) == TRUE)
					{
					taxon_val = "NULL"
					}
				chars = xmlAttrs(seqs[[seqnums2[j]]])$value
				fract_data = fract_data_vals[seqnums2[j]]
				
				# If too much is missing, replace sequence with a comment
				if (fract_data < min_fraction_to_keep_sequence)
					{
					txt = paste0(" Note: Sequence with id '", seqid, "' from specimen '", taxon_val, "' was removed as its fraction of data was only ", fract_data, ", and min_fraction_to_keep_sequence=", min_fraction_to_keep_sequence, ". ")
					seqs[[seqnums2[j]]] = xmlCommentNode(txt)
					next()
					} # END if (fract_data < min_fraction_to_keep_sequence)
				
				if (cut_all_seqs == TRUE)
					{
					txt = paste0(" Note: Sequence with id '", seqid, "' from specimen '", taxon_val, "' was removed as this whole locus is being cut, because num_seqs_that_are_complete_enough=", num_seqs_that_are_complete_enough, ", and min_numseqs=", min_numseqs, ". ")
					seqs[[seqnums2[j]]] = xmlCommentNode(txt)
					next()
					}
				} # END for (j in 1:length(seqnums2))
			
			# Replace the XML children
			xmlChildren(tmp_alignment) = seqs
			} # END if ( (length(namevals) > 0) && (sum(seqs_TF) > 0) )
		
		# Save it; not sure what to do if ALL sequences were eliminated -- probably cut locus manually!
		data_XML[[alignment_nums[i]]] = tmp_alignment
		} # END for (i in 1:length(data_XML)
	return(data_XML)
	} # END prune_seqs_based_on_completeness <- function(data_XML, min_numseqs=0, min_fraction_to_keep_sequence=0.5, isDNA=TRUE)





completeness_stats_from_data_XML <- function(data_XML)
	{
	XML_comments_TF = sapply(X=data_XML, FUN=is_XML_comment)
	alignment_nums = (1:length(data_XML))[XML_comments_TF==FALSE]
	dataset_names = rep(NA, length(alignment_nums))
	
	# Exclude XML tags with no children
	keepTF = rep(TRUE, length(alignment_nums))
	for (i in 1:length(alignment_nums))
		{
		tmp_alignment = data_XML[[alignment_nums[i]]]
		dataset_names[i] = xmlAttrs(tmp_alignment)$id
		namevals = names(tmp_alignment)
		if (length(namevals) == 0)
			{
			keepTF[i] = FALSE
			} # END if (length(namevals) == 0)
		} # END for (i in 1:length(alignment_nums))
	alignment_nums = alignment_nums[keepTF]
	
	dataset_names = rep(0, length(alignment_nums))
	dataset_ntaxa = rep(0, length(alignment_nums))
	dataset_lengths = rep(0, length(alignment_nums))
	allSameLength = rep("", length(alignment_nums)) 
	meanQs = rep(0, length(alignment_nums))
	maxQs = rep(0, length(alignment_nums))
	minQs = rep(0, length(alignment_nums))
	pctData = rep(0, length(alignment_nums))
	nMajData = rep(0, length(alignment_nums))
	
	for (i in 1:length(alignment_nums))
		{
		tmp_alignment = data_XML[[alignment_nums[i]]]
		dataset_names[i] = xmlAttrs(tmp_alignment)$id

		namevals = names(tmp_alignment)
		
		# Get the dataType
		tmp_datatype = xmlAttrs(tmp_alignment)$dataType
		if (isblank_TF(tmp_datatype) == TRUE)
			{
			# Then, it is morphology or something, don't cut
			cat("\nNote: in dataXML, xmlAttrs(tmp_alignment)$dataType returns blank. This is probably morphology; completeness_stats_from_data_XML() is tested only on DNA at the moment. Skipping...\n")
			#completeness_df = NULL
			#return(completeness_df)
			next()
			}
		if (tmp_datatype == "nucleotide")
			{
			isDNA = TRUE
			} else {
			isDNA = FALSE
			} # END if (tmp_datatype == "nucleotide")
		
		# If this alignment (XML tag) has more than 0 names, and 
		# if some of them are sequence, then this is an alignment
		seqs_TF = "sequence" %in% namevals
		xmltag_is_seqs_TF = namevals == "sequence"
		dataset_ntaxa[i] = sum(xmltag_is_seqs_TF)
		if ( (length(namevals) > 0) && (sum(seqs_TF) > 0) )
			{
			# Get the lengths of each DNA sequence, and the number
			# of nondata in each
			seqnums = (1:length(namevals))[xmltag_is_seqs_TF]
			children = xmlChildren(tmp_alignment)
			seqs = children[seqnums]
			#sapply(X=seqs, FUN=is_XML_comment)
			lengths = unlist(sapply(X=seqs, FUN=get_lengths_from_XMLtag))
			
			if (length(unique(lengths)) == 1)
				{
				dataset_lengths[i] = unique(lengths)
				allSameLength[i] = "yes"
				} else {
				dataset_lengths[i] = mean(lengths, na.rm=TRUE)
				allSameLength[i] = "no"
				}
			
			if (isDNA == TRUE)
				{
				numQ_vals = unlist(sapply(X=seqs, FUN=get_numQs_from_XMLtag, isDNA=TRUE))
				}
			if (isDNA == FALSE)
				{
				numQ_vals = unlist(sapply(X=seqs, FUN=get_numQs_from_XMLtag, isDNA=FALSE))
				}
			names(numQ_vals) = NULL
			names(lengths) = NULL
			
			fract_data_vals = 1 - (numQ_vals / lengths)
			fract_data_vals

			meanQs[i] = round(mean(numQ_vals, na.rm=TRUE), 2)
			maxQs[i] = round(max(numQ_vals, na.rm=TRUE), 2)
			minQs[i] = round(min(numQ_vals, na.rm=TRUE), 2)
			pctData_for_each_OTU = round(fract_data_vals * 100, 1)
			pctData[i] = round(mean(pctData_for_each_OTU, na.rm=TRUE),1)
			nMajData[i] = sum( pctData_for_each_OTU >= 50)				
			} # END if ( (length(namevals) > 0) && (sum(seqs_TF) > 0) )
		} # END for (i in 1:length(alignment_nums))

	completeness_df = cbind(dataset_names, dataset_lengths, allSameLength, dataset_ntaxa, nMajData, pctData, minQs, meanQs, maxQs)
	completeness_df = as.data.frame(completeness_df, stringsAsFactors=FALSE)
	names(completeness_df) = c("dataset_names", "dataset_lengths", "allSameLength", "dataset_ntaxa", "nMajData", "pctData", "minQs", "meanQs", "maxQs")
	completeness_df = dfnums_to_numeric(completeness_df)	
	return(completeness_df)
	} # END completeness_stats_from_data_XML_seqs_based_on_completeness(data_XML)
