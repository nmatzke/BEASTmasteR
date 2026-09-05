# Load the R libraries for dealing with phylogenies & XML
library(XML)
library(ape)   # for read.tree
library(gdata) # for read.xls
library(BioGeoBEARS) # for sourceall
library(XLConnect)	# for readWorksheetFromFile
library(stringr)	# for str_count
library(adegenet)	# for seqinr
library(pegas)		# for seqinr
library(seqinr)		# for seqinr

#####################################################
# Source BEASTmasteR code via the internet
# (You could also save the R files locally and source() their locations on your hard disk.
# This will be especially handy if your internet sucks. I have archived the R code in a
# dated zipfile, see "Files" at the bottom of the page):
#####################################################

# On Nick's development computer:
#sourceall("/GDrive/__github/BEASTmasteR/R/")
sourceall("~/GitHub/BEASTmasteR/R/")


#######################################################
# CHANGE SCRIPT HERE: USE YOUR OWN WORKING DIRECTORY
# NOTE: Windows uses "\\" instead of "/"
#######################################################
#wd = "/drives/Dropbox/_njm/__packages/BEASTmasteR_permahelp/examples/ex_basic_venerid_morphDNA_v4/SABD_tipsVary_wOutg_v1/"

wd = "~/GitHub/BEASTmasteR/ex/rails/"
setwd(wd)


# Create a directory to save statistics about your data
datastats_dir = "z_data_stats"
dir.create(datastats_dir)


# The name of the Excel settings file
xlsfn = "settings_v9b_Rails.xlsx"



#######################################################
# Double-check your working directory with getwd()
#######################################################
getwd()

#######################################################
# Double-check that you have the right files with list.files()
#######################################################
list.files()

#######################################################
# YOU SHOULD BE ABLE TO CUT-N-PASTE CODE FROM HERE DOWN
# (But, do it line-by-line, so you can see what is going on!)
#######################################################

#######################################################
# Required change to R settings (run every time)
#######################################################
# Turn off R's default stringsAsFactors silliness
options(stringsAsFactors = FALSE)

# When you have a starting tree named tree.newick, you could view it by un-commenting the code below:
# trfn = "tree.newick"
# tr = read.tree(trfn)
# tr_table = prt(tr)
# tr_table[,c("label","time_bp")]
# names(tr_table)

#######################################################
# Setup the overall xml structure
#######################################################
xml = setup_master_XML()
print_master_XML_section_headers(xml)


#######################################################
# Get the name of the species tree -- the main
# target of the analysis
#######################################################
#tree_name = "shared_tree"
treemodel_df = readWorksheetFromFile(xlsfn, sheet="treemodel", startRow=15)
data_df = readWorksheetFromFile(xlsfn, sheet="data", startRow=15)
data_df$use[isblank_TF(data_df$use)] = "yes"
data_df = data_df[data_df$use == "yes",]
tree_names = unique(c(treemodel_df$speciesTreeName, data_df$speciesTreeName))
if (length(tree_names) > 1)
	{
	txt = "Error in BEASTmasteR: BEASTmasteR is not currently programed to accept more than one speciesTreeName. You have these:\n\n"
	cat("\n\n")
	cat(txt)
	cat("treemodel_df$speciesTreeName:")
	print(treemodel_df$speciesTreeName)
	cat("\n\n")
	cat("data_df$speciesTreeName:")
	print(data_df$speciesTreeName)
	cat("\n\n")
	stop(txt)
	} else {
	tree_name = tree_names[1]
	}
if (is.null(tree_name))
	{
	tree_name = "shared_tree"
	}
tree_name

#clockModel_name = "shared_clock"


# Decide if this is a StarBeast2 analysis
if (grepl(pattern="starBeast2", x=treemodel_df$treeModel[1]) == TRUE)
	{
	StarBeast2_TF = TRUE

	# This is required, if StarBeast2 is being used
	taxonsets_df = readWorksheetFromFile(xlsfn, sheet="taxonsets", startRow=15)

	# Error check
	if (nrow(treemodel_df) > 1)
		{
		txt = paste0("STOP ERROR in check_for_starBeast2(): in worksheet 'treemodel', you have specified that treemodel_df$treeModel='starBeast2', meaning you want to do a starBeast2 analysis. However, starBeast2 assumes a simple birth-death model (constant rate, not skyline models), so you should only have one row in treemodel. However, you have ", nrow(treemodel_df), " rows.  Printing treemodel_df:")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		}
	} else {
	StarBeast2_TF = FALSE
	taxonsets_df = NULL
	}



# XML SECTION 1: Header and Beast2 Java class references / mappings
# Done below

# XML SECTION 2: Taxa and clade definitions (and tip-dates if desired)
#######################################################
# Get the OTUs from the Excel settings file
#######################################################
OTUs_df = readWorksheetFromFile(xlsfn, sheet="OTUs", startRow=15)
head(OTUs_df)
OTUs_raw = OTUs_df$OTUs	# raw list, for safekeeping / pasting back to Excel

# Remove OTUs with "no" in "use" column of "OTUs" worksheet
OTUs_df$use[isblank_TF(OTUs_df$use)] = "yes"
keepTF = OTUs_df$use != "no"
OTUs_df = OTUs_df[keepTF,]
OTUs = OTUs_df$OTUs
OTUs

# Counting samples by time-bin (irrelevant if everything is living in the present)
# (this is just a rough assessment of sampling-through-time -- delete if
#  not interested)
empirical_sampling_info = analyze_tipdates_for_gaps_in_sampling_EXAMPLES(OTUs_df)


#######################################################
# Write OTUs to ="list_of_OTUs"
#######################################################

# Not needed really
tmpXML = make_XML_tipdate_priors(OTUs_df, min_precision=0.01, tree_name="shared_tree", xml=NULL)
tmpXML$priors
tmpXML$operators
tmpXML$tracelog
tmpXML$tipdatelog
# Not needed really

XML_list_of_OTUs = make_XMLs_for_OTUs(OTUs, OTU_idref=FALSE, StarBeast2_TF=StarBeast2_TF)
list_of_OTUs_XML = make_XML_taxon_block(taxon_name="list_of_OTUs", XML_list_of_OTUs=XML_list_of_OTUs)
list_of_OTUs_XML

# Add to the list of taxa/clades (manually, as an example)
xml$taxa = c(xml$taxa, list_of_OTUs_XML)



#######################################################
# Add taxonSets additional clades/taxa specified in "taxa"
#######################################################
taxa_df = readWorksheetFromFile(xlsfn, sheet="taxa", startRow=15)
head(taxa_df)

# Add the taxon groups (monophyletic clades and perhaps non-monophyletic groups of interest)
xml = make_taxa_groups(taxa_df, OTUs=OTUs, xml=xml)

# Add the genera if desired (for logging monophyly, if possible; if not, log ages)
#genera = get_genera_from_OTUs(OTUs, mintaxa=2, split="_")
#genera
# Don't add genera to XML (you would only do it if you were SURE they were monophyletic)
#xml = make_XML_for_genera_from_OTUs(xml, OTUs, mintaxa=2, split="_")
#xml

#######################################################
# Add priors for the node constraints
#######################################################
nodes_df = readWorksheetFromFile(xlsfn, sheet="nodes", startRow=15)
head(nodes_df)

# Make the clade priors (for the general XML, not inside the Calibrated Yule prior;
# that is set up elsewhere if desired)
xml = make_cladePrior_XMLs(nodes_df, xml=xml, list_of_empty_taxa=xml$list_of_empty_taxa, CalibratedYule_txt="")

# Add to the logs
xml = make_cladePrior_logs(nodes_df, xml=xml, tree_name=tree_name, trace_or_screen="both", list_of_empty_taxa=xml$list_of_empty_taxa)

# XML SECTION 3: Sequence alignments (e.g. DNA, morphology); filtered in later section to produce partitions
seqs_df = readWorksheetFromFile(xlsfn, sheet="data", startRow=15)
seqs_df = seqs_df[seqs_df$use == "yes", ]
head(seqs_df)

# See starting length of the XML
orig_length_xml = length(xml)


#######################################################
# Add priors on the stem ages of terminal branches, if any
# (specified in OTUs_df)
#######################################################
# # Default is no, so only pull out "yes"
# stem_ages_df = OTUs_df[OTUs_df$stem_make_age_prior=="yes",]
# for (i in 1:nrow(stem_ages_df))
# 	{
# 	dfline = stem_ages_df[i,]
# 	param_name = paste0("stem_age_of_", dfline$OTUs)
# 	stemXML = make_generic_XML_prior(dfline, colname_prefix="stem", param_name=param_name, header_scheme=1, stem=TRUE, distrib=NULL, param1=NULL, param2=NULL, tmp_offset=NULL, meanInRealSpace=NULL)
# 	stemXML
# 	}
# stemXML



#######################################################
# Add relative priors, if any
# (specified in OTUs_df and/or nodes_df)
#######################################################
xml = make_generic_difference_statistic(nodes_df=nodes_df, OTUs_df=OTUs_df, xml=xml)
xml



################################################################
# PARSING NEXUS FILES WITH parse_datasets()
################################################################
#
# NOTE: THIS FUNCTION WILL PRODUCE A LOT OF SCREEN
# OUTPUT, UNLESS YOU SET printall="none". THE SCREEN
# OUTPUT IS GOOD, IT SHOWS YOU ALL THE REFORMATTING
# OF THE CHARACTERS TO MAKE THEM ACCEPTABLE TO
# BEAST.
#
# E.g., a particular character must have its
# lowest state be state 0, it must not skip states
# (no characters with only states 0, 3, and 4), and
# the observed number of distinct character states
# must match the number of claimed character states.
# These are typically true for original datasets,
# but after researchers cut out taxa, re-code
# certain states, etc., sometimes it is no longer
# true. The script counts the observed character
# states, re-numbers them if needed, classifies them
# by number of character states, and then reformats
# for Beast2 (e.g., for a 2-state character, "?",
# "-", "(0 1)", "(1 0)", "{0 1}", and {1 0} would
# all be converted to "Z").
#
# For your own data, you may discover yet more weird
# features of NEXUS data file formatting that I have
# not coded for. Email them to the
# beastmaster_package Google Group and I will try to
# fix them.
#
# NEXUS files must be: (1) Simplified NEXUS format
# (as exported from Mesquite), and (2) the taxon
# names must have no spaces, quotes, or other weird
# characters. Use underscores ("_") instead.
#
# Tracking characters in the NEXUS file: as a
# double-check, just to make sure the data is not
# being fundamentally changed, I have the script
# write into XML: (1) the original morphology data
# matrix, (2) the modified data matrix, and (3) the
# character numbers (original numbers or indices
# with each data section, depending) that have been
# placed into each morphology data section. This
# occurs because Beast2 requires that the 2-state,
# 3-state, 4-state, etc. characters each be coded in
# the XML in separate sections.
################################################################
# Run the parsing of NEXUS file(s)

# Save old xml
pre_parsing_xml = xml
# xml = pre_parsing_xml

# Once you've run parse_datasets successfully, you can
# change runslow=FALSE for future runs to save time
# (as long as you have data_XML.Rdata in the working
# directory).
runslow = TRUE
data_XMLs_fn = "data_XMLs.Rdata"
data_XMLs_preCompleteness_fn = "data_XMLs_preCompleteness.Rdata"
dnalist_conversion_output_fn = "dnalist_conversion_output.Rdata"
dnalist_conversion_output_preCompleteness_fn = "dnalist_conversion_output_preCompleteness.Rdata"

if (runslow)
	{
	data_XMLs = parse_datasets(seqs_df, add_morphLength=TRUE, add_morphList=TRUE, OTUs=OTUs,  printall="short", convert_ambiguous_to_IUPAC=FALSE, xml=NULL, xlsfn=xlsfn, return_charsdf=TRUE, add_allQ_DNA_seqs_in_TF=TRUE)
	data_XMLs_preCompleteness = data_XMLs
	# Save loaded data to various .Rdata files, and NEXUS separate and concatenated
	# Also produces a version of the data after completeness filters are applied
	# (more for multilocus StarBeast2 analyses)d
	save_out_matrices_and_completeness_files = save_out_matrices_and_completeness(data_XMLs_preCompleteness, datastats_dir=datastats_dir, min_numseqs=0, min_fraction_to_keep_sequence=0.5, StarBeast2_TF=StarBeast2_TF, taxonsets_df=taxonsets_df,  data_XMLs_fn=data_XMLs_fn, data_XMLs_preCompleteness_fn=data_XMLs_preCompleteness_fn, dnalist_conversion_output_fn=dnalist_conversion_output_fn, dnalist_conversion_output_preCompleteness_fn=dnalist_conversion_output_preCompleteness_fn)
	save(save_out_matrices_and_completeness_files, file="save_out_matrices_and_completeness_files.Rdata")
	} # END if (runslow)

load(data_XMLs_fn)									# Loads to "data_XML"
load(data_XMLs_preCompleteness_fn)					# Loads to "data_XMLs_preCompleteness"
load(dnalist_conversion_output_fn)					# Loads to "dnalist_conversion_output"
load(dnalist_conversion_output_preCompleteness_fn)	# Loads to "dnalist_conversion_output_preCompleteness"
load("save_out_matrices_and_completeness_files.Rdata") 	# Loads to "save_out_matrices_and_completeness_files"

# Do you want to use the dataset as loaded (pre-completeness filtering?) (default TRUE)
# (mostly use post-completeness for complex multilocus datasets)
use_preCompleteness = TRUE
if (use_preCompleteness == TRUE)
	{
	data_XMLs = data_XMLs_preCompleteness
	data_XML = data_XMLs_preCompleteness$data_XML
	dnalist_conversion_output = dnalist_conversion_output_preCompleteness
	} else {
	nothing = "blah"
	}

# Print lengths without locus names
unname(dnalist_conversion_output$dna_alignment_lengths)


#data_XMLs[[1]][3]
# If you have problems, read the help paragraphs above,
# carefully look at the messages printed to screen, the E R R O R
# messages, and carefully look at your input NEXUS file.
#
# It can also help to run read_nexus_data2 by itself:
# E.g. get first used NEXUS file in 'data' worksheet:
#
# Un-comment to run:
# nexus_filename_to_read = seqs_df$filename[seqs_df$use==TRUE][1]
# read_nexus_data2(file=nexus_filename_to_read, check_ambig_chars=TRUE, convert_ambiguous_to=NULL, printall="short", convert_ambiguous_to_IUPAC=FALSE)

# Show the length of each data matrix (and completeness statistics)
data_XMLs$dataset_lengths_df

# Completeness statistics for the first matrix
data_XMLs$matrices_stats[[1]]$completeness_df

# Show morphological characters with <2 morphological states
data_XMLs$chars_wLT_2_states_df_list

# Show the stats (both DNA and morphology)
data_XMLs$matrices_stats

# Show the number of morphological data matrices
length(data_XMLs$charsdf_list)
# Size of matrix #1
if (length(data_XMLs$charsdf_list) > 0)
	{
	data_XMLs$charsdf_list[[1]]
	data_XMLs$charsdf_list[[1]][1:5,1:5]
	}


# Get the number of question marks per taxon
# in just the DNA alignment
tdf = get_DNA_numQs_from_data_XML(xml_items_list=data_XMLs$data_XML)
outfn = paste0(datastats_dir, "/", "DNA_numQs_v1.txt")
write.table(tdf, file=outfn, sep="\t", quote=FALSE)



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
matrices_stats = data_XMLs$matrices_stats
count_constraints_df = xml$count_constraints_df
#OTUs_raw
z_Excel_OTUs_yesNo_fns = which_taxa_have_data(matrices_stats=matrices_stats, count_constraints_df=count_constraints_df, OTUs_raw=OTUs_raw)
z_Excel_OTUs_yesNo_fns

#######################################################
# Continue with adding morphology, DNA, etc. to the XML
#######################################################
xml$sequences = c(xml$sequences, data_XMLs$data_XML)
xml$misc = c(xml$misc, data_XMLs$misc_XML)
xml$sitemodels = c(xml$sitemodels, data_XMLs$sitemodels_XML)
xml$state = c(xml$state, data_XMLs$state_XML)
xml$priors = c(xml$priors, data_XMLs$prior_XML)
xml$likes = c(xml$likes, data_XMLs$likes_XMLs)
xml$operators = c(xml$operators, data_XMLs$operator_XMLs)
xml$tracelog = c(xml$tracelog, data_XMLs$tracelog_XMLs)
xml$screenlog = c(xml$screenlog, data_XMLs$screenlog_XMLs)

xml$morphLengths = data_XMLs$morphLengths
xml$morphList = data_XMLs$morphList
xml$dataset_lengths_df = data_XMLs$dataset_lengths_df
xml$morphList


# The number of (used!) characters
# of each morphology DATASET (not partition)
morphLengths = xml$morphLengths
morphLengths

# The number of (used!) characters in each morphology "partition" (section)
morphList = xml$morphList
sum(xml$morphList$numchars)
morphList

# Check that xml changed:
length(orig_length_xml)
length(xml)


####################################################################
# XML SECTION 4: Partitions
# (filters applied to the sequence alignments produce partitions)
####################################################################

# For non-morphology datasets, extract partitions from the alignments
xml = make_partitions_XML(seqs_df, xml=xml, add_partitionLength=TRUE, dataset_lengths_df=NULL)
xml$partitions

partitionLengths = xml$partitionLengths_df$partitionLengths
partitionLengths

DNA_columns_total_number = sum(partitionLengths)
DNA_columns_txt_fn = paste0(datastats_dir, "/", "DNA_columns_total_number.txt")
write.table(x=DNA_columns_total_number, file=DNA_columns_txt_fn, append=FALSE, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

#######################################################
# Add the tip-dates (assume age = 0 if none)
# All blanks/NAs are converted to 0
#######################################################
OTUs = OTUs_df$OTUs[OTUs_df$use=="yes"]


if (StarBeast2_TF == FALSE)
	{
	tipdates = OTUs_df$tipdate

	# Get the alignment that will serve as a source of OTUs to reference with "@"
	# in the tipdates 'taxa' tag
	alignment_name_w_taxa = pick_a_partitionName_for_taxa_list(seqs_df=seqs_df, morphList=morphList)

	# Make the tipdates
	OTUs_df$use[isblank_TF(OTUs_df$use)] = "yes"
	XML_tipdates = make_XML_tipdates(name="tipDates", OTUs_df, alignment_name_w_taxa=alignment_name_w_taxa, backward=TRUE, xml=NULL)

	# Add to the list of taxa/clades (here, manually, just to show the process)
	xml$taxa = c(xml$taxa, XML_tipdates)
	xml$taxa

	#######################################################
	# Add a taxonSet for the fossils, if there are any
	#######################################################
	xml = add_fossils_taxon_to_xml(OTUs, tipdates, xml)

	#######################################################
	# Add tipdate uncertainty, sampling, and logging, if desired
	#######################################################
	xml = make_XML_tipdate_priors(OTUs_df, min_precision=0.01, tree_name=tree_name, xml=xml)
	} else {

	if ((exists("alignment_name_w_taxa") == FALSE) || (is.null(alignment_name_w_taxa) == TRUE) )
		{
		taxa_list = taxa_list_for_starBEAST(xlsfn=xlsfn, sheet="taxonsets")
		}
	} # END if (StarBeast2_TF == FALSE)


# XML SECTION 5: Miscellaneous
# (none yet)



#####################################################
# SECTIONS 6-7: Partitions and clock models
#####################################################
# Load the sequence partitions and their clocks;
# Remember to subset to just the ones where use=yes
seqs_df = readWorksheetFromFile(xlsfn, sheet="data", startRow=15)
seqs_df = seqs_df[seqs_df$use != "no", ]
seqs_df = check_seqs_df(seqs_df)

# Get the number of taxa
ntaxa = length(OTUs_df$OTUs)

# XML SECTION 6: Site Models (sequence/morphology evolution models for each partition)
xml = parse_DNA_AA_siteModels(seqs_df, xml=xml, tree_name=tree_name, StarBeast2_TF=StarBeast2_TF)
xml$dataset_lengths_df


######################################################
# XML SECTION 7: Shared clock model (relaxed ucld currently; also, the strict clock would be where stdev=~0)
######################################################

# Get the names of each clock (usually there is just 1)
clockModel_names = get_clock_names(seqs_df)
for (i in 1:length(clockModel_names))
	{
	# Get the name of the clock model
	clockModel_name = clockModel_names[i]

	if (StarBeast2_TF==TRUE)
		{
		# Get corresponding name of the gene tree (for StarBeast2 analysis)
		# (and the dataset)
		TF = seqs_df$clockModel_name==clockModel_name
		gene_tree_name = seqs_df$geneTreeName[TF][1]
		dataset_name = seqs_df$datasetName[TF][1]
		TF2 = xml$dataset_lengths_df$dataset_names == dataset_name
		ntaxa = as.numeric(xml$dataset_lengths_df$dataset_ntaxa[TF2][1])
		clock_type = seqs_df$clockmodel_type[TF][1]

		# StarBeast2 analysis
		# SpeciesTreeRelaxedClock
		if (grepl(pattern=tolower("SpeciesTree"), x=tolower(clock_type)) == TRUE)
			{
			# For a SpeciesTreeRelaxedClock, you have to go through
			# all the gene trees
# 			for (j in 1:sum(TF))
# 				{
# 				# Is it the first?
# 				if (j == 1)
# 					{
# 					first = TRUE
# 					} else {
# 					first = FALSE
# 					}
# 				# Each gene tree inherits relative rate from SpeciesTree;
# 				# But here, we are just defining the SpeciesTree Relaxed Clock
# 				gene_tree_name = seqs_df$geneTreeName[TF][j]
# 				dataset_name = seqs_df$datasetName[TF][j]
# 				TF2 = xml$dataset_lengths_df$dataset_names == dataset_name
# 				ntaxa = as.numeric(xml$dataset_lengths_df$dataset_ntaxa[TF2][1])
# 				#geneTree_clockRate_id = paste0("clockRate_for_geneTree_", gene_tree_name)
				xml = define_a_shared_clock(seqs_df, ntaxa=ntaxa, clockModel_name=clockModel_name, species_tree_name=tree_name, gene_tree_name=gene_tree_name, xml=xml, first=NULL)
# 				} # END for (j in 1:sum(TF))
			} else {
			# Individual genetree relaxed clock
			# (i.e., just do it once, for this clockModel)
			xml = define_a_shared_clock(seqs_df, ntaxa=ntaxa, clockModel_name=clockModel_name, species_tree_name=tree_name, gene_tree_name=gene_tree_name, xml=xml, first=NULL)
			} # END if (grepl(pattern=tolower("SpeciesTree"), x=tolower(clock_type)) == TRUE)

		} else {
		# Regular Beast2 analysis
		tree_name = tree_name
		xml = define_a_shared_clock(seqs_df, ntaxa, clockModel_name=clockModel_name, species_tree_name=tree_name, xml=xml)
		} # END if (StarBeast2_TF==TRUE)
	} # END for (i in 1:length(clockModel_names))
xml$clock



# Do the XML for relativeMutationRates
# (the relative rate of each partition)

# For multiple clocks: Put this in a loop with
# different partitionLengths (DNA) and morphLengths (morphology)
# for each clock
partitionLengths_df = xml$partitionLengths_df
morphLengths = xml$morphLengths

clockModel_names = unique(seqs_df$clockModel_name)
clockModel_name=clockModel_names[1]
# Test making of mutationRate operators
make_relativeMutationRate_operator(seqs_df, partitionLengths_df=partitionLengths_df, morphLengths=morphLengths, clockModel_name=clockModel_name, relrate_suffix = "", force_relRates=FALSE, xml=NULL)

for (i in 1:length(clockModel_names))
	{
	clockModel_name = clockModel_names[i]
	xml = make_relativeMutationRate_operator(seqs_df, partitionLengths_df=partitionLengths_df, morphLengths=morphLengths, clockModel_name=clockModel_name, relrate_suffix = "", xml=xml)
	}
#xml$operators


#######################################################
# Make relative mutation rates, ONLY for the situation
# where there are multiple partitions
#######################################################

# XML for morphology models, priors, likelihoods
if (!is.null(morphList))
	{
	# Find unique clocks or relRates
	txts = paste(morphList$clockModel_name, morphList$clockModel_relRate, collapse=NULL, sep=".")
	morph_types = unique(txts)
	for (i in 1:length(morph_types))
		{
		rowsTF = txts == morph_types[i]
		clockModel_name = morphList$clockModel_name[rowsTF][1]
		morphRate_name = morphList$clockModel_relRate[rowsTF][1]
		morphGamma_name = morphList$gammaShape_suffix[rowsTF][1]
		xml = make_Beast2_morph_models(morphList[rowsTF,], morphRate_name=morphRate_name, morphGamma_name=morphGamma_name, clockModel_name=clockModel_name, tree_name=tree_name, xml=xml)
		} # END for (i in 1:length(morph_types))
	#xml = make_Beast2_morph_models(morphList[1:4,], morphRate_name="morph_relRate", morphGamma_name="morph_gammaShape", clockModel_name=clockModel_name, tree_name=tree_name, xml=xml)
	# NJM: for 2nd morph matrix
	#xml = make_Beast2_morph_models(morphList[5:7,], morphRate_name="morphSl_relRate", morphGamma_name="morphSl_gammaShape", clockModel_name=clockModel_name, tree_name=tree_name, xml=xml)
	} # END if (!is.null(morphList))


# Get the logEvery for the speciesTree logger
run_df = readWorksheetFromFile(xlsfn, sheet="run", startRow=15)
run_df


# XML SECTION 8: Shared tree model (Yule, Birth-Death (BD), Birth-Death Skyline (BDSKY), etc.)
treemodel_df = readWorksheetFromFile(xlsfn, sheet="treemodel", startRow=15)

# Try this
#constructed_starting_tree_fn = "constructed_starting_tree.newick"
#min_brlen = 0.001
#starting_tree = construct_starting_tree(xlsfn=xlsfn, min_brlen=min_brlen, outfn=constructed_starting_tree_fn)
#write.tree(starting_tree, file="")

# Does the dataset include continuous characters?  The XML for the tree needs to
# be different depending on whether or not this is so!  See:
# https://groups.google.com/forum/#!topic/beast-users/7KFzWd4PVeQ
TF = seqs_df$use != "no"
if (("continuous" %in% seqs_df$type[TF]) == TRUE)
	{
	XML_mod_for_cont_chars = TRUE
	} else {
	XML_mod_for_cont_chars = FALSE
	} # END if (("continuous" %in% seqs_df$type[TF]) == TRUE)


if (StarBeast2_TF == FALSE)
	{
	# Traditional concatenated analysis

	# Get a partition for the taxa when building the starting tree
	alignment_name_w_taxa = pick_a_partitionName_for_taxa_list(seqs_df=seqs_df, morphList=morphList)

	#xmltmp = make_BDSKY_model(treemodel_df, tree_name=tree_name, clockModel_name=clockModel_name, alignment_name_w_taxa=alignment_name_w_taxa, tipDates_id="tipDates", xml=NULL, XML_mod_for_cont_chars=FALSE)

	xml = make_BDSKY_model(treemodel_df, tree_name=tree_name, clockModel_name=clockModel_name, alignment_name_w_taxa=alignment_name_w_taxa, taxonset_w_taxa="list_of_OTUs", tipDates_id="tipDates", xml=xml, XML_mod_for_cont_chars=XML_mod_for_cont_chars, xlsfn=xlsfn)

	} else {
	#######################################################
	# StarBeast2 setup
	#######################################################
	# Make a taxonsuperset for the species tree
	speciesTree_taxonset_listname="taxonsuperset"
	taxonsets_df = readWorksheetFromFile(xlsfn, sheet="taxonsets", startRow=15)
	taxonset_XML = make_taxon_superset(taxonsets_df=taxonsets_df, speciesTree_taxonset_listname=speciesTree_taxonset_listname)


	xml = make_BD_model_for_starBeast2(treemodel_df=treemodel_df, seqs_df=seqs_df, taxonset_XML=taxonset_XML, tree_name=tree_name, speciesTree_default=tree_name, birthRate_id="netDiversificationRate", logEvery=run_df$treelog_store, xml=xml, nodes_df=nodes_df, min_brlen=0.001)

	xml$likes
	xml = make_speciescoalescent_for_starBeast2(OTUs_df=OTUs_df, seqs_df=seqs_df, treemodel_df=treemodel_df, tree_name=tree_name, speciesTree_default=tree_name, xml=xml)
	xml$likes


	###################################################
	# Make the tree stateNode for the species tree
	###################################################
	#xml = make_starting_tree_for_starBeast2(treemodel_df, taxonset_XML, xlsfn=xlsfn, speciesTree_default="speciesTree", min_brlen=0.001, xml=xml)

	#######################################################
	# Operators on the gene trees
	#######################################################
	xml = genetrees_operators(seqs_df=seqs_df, taxonsets_df=taxonsets_df, treemodel_df=treemodel_df, tree_name=tree_name, speciesTree_taxonset_listname=speciesTree_taxonset_listname, speciesTree_default=tree_name, xml=xml)
	} # END if (StarBeast2_TF == FALSE)



# These are all done above, with further additions by parse_run()
# XML SECTION 9: Starting tree (random, user-fixed, etc.)
# XML SECTION 10: MCMC states: the state of the MCMC chain is saved at various timepoints,
# XML SECTION 11: Priors
# XML SECTION 12: Likelihoods
# XML SECTION 13: Operators. These specify how the parameter values can change at MCMC time-steps.
# XML SECTION 14: Trace log. This logs the parameters in a .log file, which may be viewed in Tracer.
# XML SECTION 15: Screen log. Which parameters print to screen, and how often?
# XML SECTION 16: Tree log. How often should trees from the MCMC chain be saved, and what data should they have?

# STICK IT ALL TOGETHER
#run_df = readWorksheetFromFile(xlsfn, sheet="run", startRow=15)
run_df

# Save to outfn
# This is where it will get saved:
getwd()
# outfn="run_df" means the output filename will is specified in the run_df worksheet

# Input the dataset source
dataset_source = ' Dataset: Juans Rallidae+outgroups dataset, DNA+morph; BEASTmasteR by Nick Matzke '

outfn = parse_run(run_df, xml, outfn="run_df", dataset_source=dataset_source,  printall="short")

# Here is the filename
outfn

# Glance at the XML file in R (will take up many screens)
#moref(outfn)

# Save the morphology matrix stats to text file(s):
outfns = save_data_stats(data_XMLs, datastats_dir=datastats_dir)


