
#######################################################
# Based on: _which_OTUs_have_data_v2.R
#######################################################

tally_data <- function(dtf_list=NULL, fns=NULL, dtf_list_fn=NULL, xlsfn=NULL, OTUs_df=NULL, taxa_df=NULL, nodes_df=NULL, outfn=NULL)
{
defaults='
fns = c(
"canidae_DNA_wQs_v3.nex",
"canidae_simp2_wQs_v2.nex",
"CLA_298_sm_Appendix_S3_discrete_v5_morph_simp2_wQs.nex"
)

dtf_list=dtf_list
fns=NULL
dtf_list_fn=NULL
xlsfn=xlsfn
OTUs_df=NULL
taxa_df=NULL
nodes_df=NULL
outfn=NULL

'

require(XLConnect)

gather_OTUs = FALSE
if (is.null(xlsfn) && is.null(OTUs_df))
	{
	gather_OTUs = TRUE
	} # END if (is.null(xlsfn) && is.null(OTUs_df))

if (is.null(xlsfn)==FALSE)
	{
	OTUs_df = readWorksheetFromFile(xlsfn, sheet="OTUs", startRow=15)
	head(OTUs_df)

	# Not needed: instead, record "use" column
	# Remove OTUs with "no" in "use" column of "OTUs" worksheet
	remove_OTUs = FALSE
	if (remove_OTUs)
		{
		OTUs_df$use[isblank_TF(OTUs_df$use)] = "yes"
		keepTF = OTUs_df$use != "no"
		OTUs_df = OTUs_df[keepTF,]
		}
	# INSTEAD, record this in the table
	OTUs = OTUs_df$OTUs
	OTUs

	# 
	OTUs_df = readWorksheetFromFile(xlsfn, sheet="OTUs", startRow=15)
	taxa_df = readWorksheetFromFile(xlsfn, sheet="taxa", startRow=15)
	nodes_df = readWorksheetFromFile(xlsfn, sheet="nodes", startRow=15)
	}


if (is.null(taxa_df))
	{
	txt = "STOP ERROR in tally_data(): at least one of taxa_df, xlsfn must not be NULL."
	cat("\n\n")
	cat(txt)
	cat("\n\n")
	stop(txt)
	}


if (is.null(nodes_df))
	{
	txt = "STOP ERROR in tally_data(): at least one of nodes_df, xlsfn must not be NULL."
	cat("\n\n")
	cat(txt)
	cat("\n\n")
	stop(txt)
	}


if (is.null(dtf_list) && is.null(dtf_list_fn) && is.null(fns))
	{
	txt = "STOP ERROR in tally_data(): at least one of dtf_list, fns, dtf_list_fn must not be NULL."
	
	cat("\n\n")
	cat(txt)
	cat("\n\n")

	stop(txt)
	} # END if (is.null(dtf_list) && is.null(dtf_list_fn) && is.null(fns))


# Load the data
data_loaded = FALSE

if (is.null(dtf_list) == FALSE)
	{
	# OK, you have data, you are good
	data_loaded = TRUE
	}

if (data_loaded==FALSE && is.null(dtf_list_fn)==FALSE)
	{
	# Loads to "dtf_list"
	load(file=dtf_list_fn)
	data_loaded = TRUE
	}


if (data_loaded==FALSE && is.null(fns)==FALSE)
	{
	####################################################
	# Load the NEXUS files... this can be slow
	####################################################

	runslow = TRUE
	dtf_list_fn = "dtf_list.Rdata"
	dtf_list = NULL
	if (runslow)
		{
# 		# DNA - don't check for ambig characters for this
# 		dtf_list[[1]] = read_nexus_data2(file=fns[1], convert_ambiguous_to_IUPAC=FALSE, check_ambig_chars=FALSE, printall=TRUE, convert_ambiguous_to=NULL)
# 	
# 		# DNA - don't check for ambig characters for this
# 		dtf_list[[2]] = read_nexus_data2(file=fns[2], convert_ambiguous_to_IUPAC=FALSE, check_ambig_chars=TRUE, printall=TRUE, convert_ambiguous_to=NULL)
# 		dtf_list[[3]] = read_nexus_data2(file=fns[3], convert_ambiguous_to_IUPAC=FALSE, check_ambig_chars=TRUE, printall=TRUE, convert_ambiguous_to=NULL)

		for (i in 1:length(fns))
			{
			dtf_list[[i]] = read_nexus_data2(file=fns[i], convert_ambiguous_to_IUPAC=FALSE, check_ambig_chars=FALSE, printall=TRUE, convert_ambiguous_to=NULL)
			}
		
		save(dtf_list, file=dtf_list_fn)
		} else {
		# Loads to "dtf_list"
		load(file=dtf_list_fn)
		} # END if (runslow)
	data_loaded = TRUE
	}



if (gather_OTUs == TRUE)
	{
	OTUs_to_gather = NULL

	for (i in 1:length(fns))
		{
		OTUs_to_gather = c(OTUs_to_gather, names(dtf_list[[i]]))
		}
	
	OTUs = unique(OTUs_to_gather)
	OTUs = OTUs[isblank_TF(OTUs)==FALSE]

	OTUs_df = NULL
	OTUs_df$use = rep("yes", length(OTUs))
	}

OTUs = OTUs[isblank_TF(OTUs)==FALSE]
	
if (is.null(fns))
	{
	fns_cols = paste("file", 1:length(dtf_list), sep="")
	fns_cols
	} else {
	fns_cols = fns
	fns_cols
	}


#######################################################
# Now, look for sequence/morphology data
#######################################################

####################################################
# Query for all "?"
####################################################
data_tally = matrix(0, nrow=length(OTUs), ncol=3+1+length(fns_cols))

# Check if OTUs are in meaningful clade constraints
clade_names = nodes_df$Taxon[nodes_df$use == "yes"]
data_tally[,1] = OTUs
data_tally[,2] = OTUs_df$use

# Cut the meaningless ones
clade_names_orig = clade_names
if (length(clade_names_orig) >= 3)
	{
	clade_names = clade_names[3:length(clade_names)]

	# For each clade, get the list of OTUs, add to data_tally
	for (c in 1:length(clade_names))
		{
		clade_name = clade_names[c]
		cmdstr = paste0("OTUs_w_constraint = taxa_df$", clade_name)
		eval(parse(text=cmdstr))
	
		TF = OTUs %in% OTUs_w_constraint
		data_tally[TF,3] = 1
		} # END for (c in 1:length(clade_names))
	} else {
	# Otherwise, mark everything as 0
	data_tally[,3] = 0
	}
data_tally


# For each sequence dataset, get the list of OTUs that are 
# *not* all "?"
for (f in 1:length(fns_cols))
	{
	dtf = dtf_list[[f]]
	OTUs_in_NEXUS = names(dtf)
	
	# Check for allQs
	TF = rep(FALSE, length(OTUs_in_NEXUS))
	for (i in 1:length(dtf))
		{
		if ( all(dtf[[i]] == "?") == FALSE)
			{
			TF[i] = TRUE
			} # END if ( all(dtf[[i]] == "?") == FALSE)
		} # END for (i in length(dtf):1)
		
	data_tally[TF,3+f] =  1
	} # END for (f in 1:length(fns_cols))

data_tally_orig = data_tally

data_tally = as.data.frame(x=data_tally, row.names=NULL, stringsAsFactors=FALSE)
names(data_tally) = c("OTUs", "use", "clade_constr", fns_cols, "any_data")
data_tally = dfnums_to_numeric(data_tally)

# Fill the last column with whether or not any of the preceding are true
data_tally[, ncol(data_tally)] = rowSums(data_tally[,3:(ncol(data_tally)-1)]) > 0

data_tally = dfnums_to_numeric(data_tally)
data_tally$any_data = as.numeric(data_tally$any_data)
cls.df(data_tally)

# Count up each column
colSums(data_tally[,3:ncol(data_tally)])
nrow(data_tally)

#######################################################
# Taxa for which we have no data/constraint at the moment:
#######################################################
TF = data_tally$any_data==FALSE
OTUs_with_NO_data = data_tally$OTUs[TF]
cat(OTUs_with_NO_data, sep="\n")


#######################################################
# USED taxa for which we have no data/constraint at the moment:
#######################################################
TF = ((data_tally$any_data==FALSE) + (data_tally$use=="yes")) == 2
OTUs_with_NO_data = data_tally$OTUs[TF]
cat(OTUs_with_NO_data, sep="\n")




#######################################################
# Output the table of data presence/absence
#######################################################
if (is.null(outfn))
	{
	outfn = "data_tally_v1.txt"
	}

write.table(data_tally, file=outfn, quote=FALSE, sep="\t")
cat("\nData completeness data for each dataset written to: '", outfn, "'.\n", sep="")

# data_tally = read.table(outfn)
# data_tally
# colSums(data_tally[,3:ncol(data_tally)])
# nrow(data_tally)

return(data_tally)
} # END tally_data <- function ()
