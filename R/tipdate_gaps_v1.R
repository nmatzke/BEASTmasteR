

#######################################################
# Run some basic checks on tipdates
#######################################################

check_tipdates <- function(OTUs_df)
	{
	# Are any of the distributions NA/missing?
	blankTF = isblank_TF(OTUs_df$distribution)
	if ( any(blankTF) == TRUE )
		{
		sumEmpty = sum(blankTF)
		names_empty = OTUs_df$OTUs[blankTF]
		nums_empty = (1:length(OTUs_df$OTUs))[blankTF]
		empty_df = as.data.frame(cbind(nums_empty, names_empty))
		names(empty_df) = c("num", "OTU")
		
		cat("\n\n")
		error_txt = paste0("ERROR in check_tipdates(OTUs_df): You have ", sumEmpty, " tips with no 'distribution' specified in your Excel 'OTUs' worksheet.\nSet these entries to 'fixed', at least.")
		cat(error_txt)
		cat("\n\n")
		
		cat("OTUs with blank cells in 'distribution' column:\n\n")
		print(empty_df)
		cat("\n\n")
		
		stop(error_txt)
		} # END if ( any(blankTF) == TRUE )


	# Check for allowed distributions
	allowed_distributions = c("fixed", "uniform", "normal", "lognormal")
	
	nomatchTF = (OTUs_df$distribution %in% allowed_distributions) == FALSE
	if ( any(nomatchTF) == TRUE )
		{
		sum_nomatch = sum(nomatchTF)
		names_nomatch = OTUs_df$OTUs[nomatchTF]
		nums_nomatch = (1:length(OTUs_df$OTUs))[nomatchTF]
		distributions_nomatch = OTUs_df$distribution[nomatchTF]
		nomatch_df = as.data.frame(cbind(nums_nomatch, names_nomatch, distributions_nomatch))
		names(nomatch_df) = c("num", "OTU", "distribution")
		
		cat("\n\n")
		error_txt = paste0("ERROR in check_tipdates(OTUs_df): You have ", sum_nomatch, " tips with no allowed 'distribution' specified in your Excel 'OTUs' worksheet.\nSet these entries to 'fixed', at least.")
		cat(error_txt)
		cat("\n\n")
		
		cat("Allowed distributions:\n\n")
		print(allowed_distributions)
		
		cat("\n\nOTUs with cells in 'distribution' column that don't match the allowed distributions:\n")
		print(nomatch_df)
		cat("\n\n")
		
		stop(error_txt)
		} # END if ( any(nomatchTF) == TRUE )
		
	# Check for missing dates
	TF = (OTUs_df$distribution == "fixed") == FALSE
	OTUs_df = OTUs_df[TF,]

	nomatchTF = isblank_TF(OTUs_df$tipdate)
	if ( any(nomatchTF) == TRUE )
		{
		sum_nomatch = sum(nomatchTF)
		names_nomatch = OTUs_df$OTUs[nomatchTF]
		nums_nomatch = (1:length(OTUs_df$OTUs))[nomatchTF]
		tipdates_nomatch = OTUs_df$tipdate[nomatchTF]
		nomatch_df = as.data.frame(cbind(nums_nomatch, names_nomatch, tipdates_nomatch))
		names(nomatch_df) = c("num", "OTU", "tipdate")
		
		cat("\n\n")
		error_txt = paste0("ERROR in check_tipdates(OTUs_df): You have ", sum_nomatch, " tips with no 'tipdate' specified in your Excel 'OTUs' worksheet.\nSet these dates to '0', at least.")
		cat(error_txt)
		
		cat("\n\n")
		cat("OTUs with blank cells in 'tipdate' column:\n\n")
		print(nomatch_df)
		cat("\n\n")
		
		stop(error_txt)
		} # END if ( any(nomatchTF) == TRUE )
	
	
	# All of the non-fixed tipdates need to have 2 parameters
	fossil_TF = OTUs_df$tipdate > 0
	if (sum(fossil_TF) > 0)
		{
		fossils_df = OTUs_df[fossil_TF, ]
		
		param1_blankTF = isblank_TF(fossils_df$param1)
		param2_blankTF = isblank_TF(fossils_df$param2)
		
		params_blankTF = (param1_blankTF + param2_blankTF) > 0
		
		if (sum(params_blankTF) > 0)
			{
			sum_noparam = sum(params_blankTF)
			names_noparam = fossils_df$OTUs[params_blankTF]
			nums_noparam = (1:length(fossils_df$OTUs))[params_blankTF]
			tipdates_noparam = fossils_df$tipdate[params_blankTF]
			param1_noparam = fossils_df$param1[params_blankTF]
			param2_noparam = fossils_df$param2[params_blankTF]
			noparam_df = as.data.frame(cbind(nums_noparam, names_noparam, tipdates_noparam, param1_noparam, param2_noparam))
			names(noparam_df) = c("num", "OTU", "tipdate", "param1", "param2")
			
			cat("\n\n")
			error_txt = paste0("ERROR in check_tipdates(OTUs_df): You have ", sum_noparam, " tips with a missing 'param1' and/or 'param2' in your Excel 'OTUs' worksheet.")
			cat(error_txt)

			cat("\n\n")
			cat("OTUs missing param1 or param2:\n\n")
			print(noparam_df)
			cat("\n\n")

			stop(error_txt)
			} # END if (sum(params_blankTF) > 0)
		} # END if (sum(fossil_TF) > 0)
	
	# All of the non-fixed tipdates need to have 2 parameters
	uniform_TF = OTUs_df$distribution == "uniform"
	if (sum(uniform_TF) > 0)
		{
		uniform_df = OTUs_df[uniform_TF, ]
		uniform_param2_too_small_TF = uniform_df$param2 < uniform_df$param1
		
		if (sum(uniform_param2_too_small_TF) > 0)
			{
		
			sum_param2_too_small = sum(uniform_param2_too_small_TF)
			names_param2_too_small = uniform_df$OTUs[uniform_param2_too_small_TF]
			nums_param2_too_small = (1:length(uniform_df$OTUs))[uniform_param2_too_small_TF]
			tipdates_param2_too_small = uniform_df$tipdate[uniform_param2_too_small_TF]
			param1_param2_too_small = uniform_df$param1[uniform_param2_too_small_TF]
			param2_param2_too_small = uniform_df$param2[uniform_param2_too_small_TF]
			param2_too_small_df = as.data.frame(cbind(nums_param2_too_small, names_param2_too_small, tipdates_param2_too_small, param1_param2_too_small, param2_param2_too_small))
			names(param2_too_small_df) = c("num", "OTU", "tipdate", "param1", "param2")
			
			cat("\n\n")
			error_txt = paste0("ERROR in check_tipdates(OTUs_df): You have ", sum_param2_too_small, " tips with a 'uniform' distribution where param2 (older age) is less than or equal to param1 (which should be younger, i.e. smaller age in Ma). 'param2' should be the older (i.e. larger age in Ma). Please edit your Excel 'OTUs' worksheet.")
			cat(error_txt)

			cat("\n\n")
			cat("OTUs with param1 (which should be younger, i.e. smaller age in Ma) > param2 (which should be older, i.e. larger age in Ma) column:\n\n")
			print(param2_too_small_df)
			cat("\n\n")

			stop(error_txt)
			} # END if (sum(uniform_param2_too_small_TF) > 0)
		} # END if (sum(uniform_TF) > 0)
		
	
	# If no problems
	return("passes")
	}



# Fine the gaps amongst the fossil tipdates
find_sampling_gaps_from_tipdates <- function(OTUs_df, bintops=seq(0,200,10), bot_of_bins=1000, collapse=FALSE, psiSamplingStopsAt=0.01)
	{
	defaults='
	bintops=seq(0,200,by=10)
	bot_of_bins=1000
	psiSamplingStopsAt=0.1
	'
	
	# Get tipdates (midpoints, means, or starting values)
	tipdates_orig = OTUs_df$tipdate
	
	# Remove 0s (actually, anything close enough to zero to be
	# caught in the "tolerance" value)
	TF = tipdates_orig >= psiSamplingStopsAt
	# Add psiSamplingStopsAt as the top bin (before 0, added at the end)
	tipdates = sort(tipdates_orig[TF])
	
	# Gaps between them
	bintops = c(bintops, bot_of_bins)
	tops = bintops[1:(length(bintops)-1)]
	tops[1] = psiSamplingStopsAt
	bots = bintops[2:(length(bintops))]
	gapsize = bots-tops
	count = rep(0, length(tops))
	counts_bybin = as.data.frame(cbind(tops, bots, gapsize, count))
	names(counts_bybin) = c("tops", "bots", "gapsize", "count")
	
	# Number of fossils in each bin
	for (i in 1:length(tipdates))
		{
		date_below_top_TF = tipdates[i] > tops
		date_above_bot_TF = tipdates[i] <= bots
		in_bin_TF = (date_below_top_TF + date_above_bot_TF) == 2
		counts_bybin$count[in_bin_TF] = counts_bybin$count[in_bin_TF] + 1
		}
	

	if (collapse == TRUE)
		{
		counts_bybin2 = NULL
		TF = counts_bybin$count > 0
		
		for (i in 1:nrow(counts_bybin))
			{
			if (i == 1)
				{
				state1 = TF[i]
				counts_bybin2 = counts_bybin[i,]
				j = 1
				next()
				}
			state2 = TF[i]
			
			# Add the rows if the same,
			# Add new row if different
			if (state2 == state1)
				{
				counts_bybin2$count[j] = counts_bybin2$count[j] + counts_bybin$count[i]
				} else {
				counts_bybin2 = rbind(counts_bybin2, counts_bybin[i,])
				j = j+1
				}
			state1 = state2
			} # END for (i in 1:nrow(counts_bybin))
		
		# Recalculate gapsize
		counts_bybin2$bots[1:(nrow(counts_bybin2)-1)] = counts_bybin2$tops[2:(nrow(counts_bybin2))]
		counts_bybin2$bots[nrow(counts_bybin2)] = bot_of_bins
		
		counts_bybin2$gapsize = counts_bybin2$bots - counts_bybin2$tops
		return(counts_bybin2)
		} # END if (collapse == TRUE)

	return(counts_bybin)
	}


convert_nonUniform_dates_to_uniform <- function(OTUs_df, CI=0.999)
	{
	defaults='
	CI=0.999
	'
	
#	OTUs_df is just a data.frame with these columns:
# 	tipdate
# 	distribution
# 	param1
# 	param2
# 	offset
# 	meanInRealSpace
	
	# Error check
	if ((CI <= 0) || (CI >=1))
		{
		errortxt = paste0("STOP ERROR in convert_nonUniform_dates_to_uniform(): CI must be greater than 0, and less than 1. Your CI=", CI)
		
		cat("\n\n")
		cat(errortxt)
		cat("\n\n")
		stop(errortxt)
		} # END if ((CI <= 0) || (CI >=1))

	# Get the lower and upper p intervals
	lowerCI = (1-CI) / 2
	upperCI = 1-lowerCI
	
	
	# Find the tipdates that are not fixed, that have non-uniform distributsion
	nonfixed_TF = OTUs_df$distribution != "fixed"
	nonuniform_TF = OTUs_df$distribution != "uniform"
	distribs_to_convert_TF = (nonfixed_TF + nonuniform_TF) == 2
	distribs_to_convert_nums = (1:length(OTUs_df$distribution))[distribs_to_convert_TF]
	
	if (sum(distribs_to_convert_TF) > 0)
		{
		#distribs_to_convert_df = OTUs_df[distribs_to_convert_TF, ]
		
		for (i in 1:sum(distribs_to_convert_TF))
			{
			rownum = distribs_to_convert_nums[i]
			rowdf = OTUs_df[rownum, ]

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
				} else {
				meanInRealSpace = TRUE
				}
			
			if (rowdf$distribution == "normal")
				{
				meanval = rowdf$param1
				sdval = rowdf$param2
				
				# Calculate the approximate bounds of the uniform 
				# distribution
				minval = qnorm(p=lowerCI, mean=meanval, sd=sdval) + offset
				maxval = qnorm(p=upperCI, mean=meanval, sd=sdval) + offset
				} # END if (rowdf$distribution == "normal")

			if (rowdf$distribution == "lognormal")
				{
				meanval = rowdf$param1
				sdval = rowdf$param2

				# Convert to real space, if needed
				if (meanInRealSpace == FALSE)
					{
					meanval = exp(meanval)
					sdval = exp(meanval)
					}
				 
				# Calculate the approximate bounds of the uniform 
				# distribution
				# Inputs means SHOULD BE MEANS IN REAL SPACE, possibly with offset
				minval = qlnorm(p=lowerCI, meanlog=log(meanval), sdlog=log(sdval)) + offset
				maxval = qlnorm(p=upperCI, meanlog=log(meanval), sdlog=log(sdval)) + offset
				} # END if (rowdf$distribution == "normal")

			if (rowdf$distribution == "exponential")
				{
				meanval = rowdf$param1

				# Convert to real space, if needed
				if (meanInRealSpace == FALSE)
					{
					meanval = 1/(meanval)
					}
				 
				# Calculate the approximate bounds of the uniform 
				# distribution
				minval = qexp(p=lowerCI, rate=(1/meanval)) + offset
				maxval = qexp(p=upperCI, rate=(1/meanval)) + offset
				} # END if (rowdf$distribution == "normal")


			rowdf$param1 = minval
			rowdf$param2 = maxval
			rowdf$distribution = "uniform"
			
			# Put back into OTUs_df
			OTUs_df[rownum, ] = rowdf
			} # END for (i in 1:nrow(distribs_to_convert_df))
		} # END if (sum(uniform_TF) > 0)
		
	return(OTUs_df)		
	}



# Convert all tipdates to normal
convert_dates_to_normal <- function(OTUs_df, sd_fraction=0.01, convert_fixed=TRUE)
	{
	defaults='
	sd_fraction=0.01
	convert_fixed=TRUE
	'
	
	# Find the tipdates that are not fixed, that have non-uniform distributsion
	if (convert_fixed == FALSE)
		{
		nonfixed_TF = OTUs_df$distribution != "fixed"
		distribs_to_convert_TF = nonfixed_TF
		} else {
		nonfixed_TF = rep(TRUE, times=length(OTUs_df$distribution))
		}
	
	
	nonnormal_TF = OTUs_df$distribution != "normal"
	distribs_to_convert_TF = (nonfixed_TF + nonnormal_TF) == 2
	distribs_to_convert_nums = (1:length(OTUs_df$distribution))[distribs_to_convert_TF]
	
	if (sum(distribs_to_convert_TF) > 0)
		{
		#distribs_to_convert_df = OTUs_df[distribs_to_convert_TF, ]
		
		for (i in 1:sum(distribs_to_convert_TF))
			{
			rownum = distribs_to_convert_nums[i]
			rowdf = OTUs_df[rownum, ]

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
				} else {
				meanInRealSpace = TRUE
				}
			
			# Change fixed dates, EXCEPT if they are 0
			if (rowdf$distribution == "fixed")
				{
				if (rowdf$tipdate == 0)
					{
					# Don't modify anything
					next()
					} # END if (rowdf$tipdate == 0)
				
				meanval = rowdf$tipdate  + offset
				sdval = sd_fraction * meanval
				} # END if (rowdf$distribution == "normal")


			if (rowdf$distribution == "uniform")
				{
				meanval = ((rowdf$param1 + rowdf$param2) / 2)  + offset
				sdval = sd_fraction * meanval
				} # END if (rowdf$distribution == "normal")

			if (rowdf$distribution == "lognormal")
				{
				meanval = rowdf$param1
				sdval = rowdf$param2
				
				# Convert to real space, if needed
				if (meanInRealSpace == FALSE)
					{
					meanval = exp(meanval) + offset
					sdval = sd_fraction * meanval
					}
				} # END if (rowdf$distribution == "normal")

			if (rowdf$distribution == "exponential")
				{
				meanval = rowdf$param1
				 
				# Convert to real space, if needed
				if (meanInRealSpace == FALSE)
					{
					meanval = (1/meanval) + offset
					sdval = sd_fraction * meanval
					}
				} # END if (rowdf$distribution == "normal")


			rowdf$param1 = meanval
			rowdf$param2 = sdval
			rowdf$distribution = "normal"
			
			# Put back into OTUs_df
			OTUs_df[rownum, ] = rowdf
			} # END for (i in 1:nrow(distribs_to_convert_df))
		} # END if (sum(uniform_TF) > 0)
		
	return(OTUs_df)		
	}



# Find the gaps amongst the fossil tipdates
# (Dates of 0 are removed)
# (non-uniform distributions are converted to uniform ones, given "CI")
find_sampling_gaps_from_tipdate_ranges <- function(OTUs_df, CI=0.999, bot_of_bins=1000, psiSamplingStopsAt=0.01, min_precision=0.1 )
	{
	defaults='
	bintops=seq(0,200,by=10)
	bot_of_bins=1000
	psiSamplingStopsAt=0.1
	min_precision = 0.1
	'

	OTUs_df_orig = OTUs_df
	# OTUs_df = OTUs_df_orig

	# Get tipdates (midpoints, means, or starting values)
	tipdates_orig = OTUs_df$tipdate

	# Remove 0s (actually, anything close enough to zero to be
	# caught in the "tolerance" value)
	fossil_TF = tipdates_orig >= psiSamplingStopsAt
	# Add psiSamplingStopsAt as the top bin (before 0, added at the end)
	tipdates = sort(tipdates_orig[fossil_TF])

	# Tipdate ranges
	# Convert non-uniform tipdates to uniform
	unifOTUs_df = convert_nonUniform_dates_to_uniform(OTUs_df[fossil_TF, ], CI=0.999)
	unifOTUs_df

	# Fixed fossil tipdates
	# (Not needed, already in unifOTUs_df)
# 	fixed_fossils_TF = OTUs_df$distribution[fossil_TF] == "fixed"
# 	fixed_OTUs_df = OTUs_df[fossil_TF, ][fixed_fossils_TF, ]
# 	fixed_OTUs_df$param1 = fixed_OTUs_df$tipdate
# 	fixed_OTUs_df$param2 = fixed_OTUs_df$tipdate
# 	unifOTUs_df = cbind(unifOTUs_df, fixed_OTUs_df)
	fixed_fossils_TF = unifOTUs_df$distribution == "fixed"
	unifOTUs_df$param1[fixed_fossils_TF] = unifOTUs_df$tipdate[fixed_fossils_TF] - min_precision/2
	unifOTUs_df$param2[fixed_fossils_TF] = unifOTUs_df$tipdate[fixed_fossils_TF] + min_precision/2

	# Tipdate range tops
	tipdate_rangetops = unifOTUs_df$param1
	tipdate_rangebots = unifOTUs_df$param2
	cbind(tipdate_rangetops, tipdate_rangebots)

	# Gaps with/without sampling
	timesteps = seq(0, bot_of_bins, min_precision)

	timestep_below_unif_rangetops_TF = rep(NA, length(tipdate_rangetops))
	timestep_above_unif_rangebots_TF = rep(NA, length(tipdate_rangebots))
	num_timesteps_inside_unif = rep(0, length(timesteps))
	for (i in 1:length(tipdate_rangetops))
		{
		timesteps_below_unif_rangetops_TF = timesteps > tipdate_rangetops[i]
		timesteps_above_unif_rangebots_TF = timesteps <= tipdate_rangebots[i]
		inside_unifs_TF = (timesteps_below_unif_rangetops_TF + timesteps_above_unif_rangebots_TF) == 2
		num_timesteps_inside_unif = num_timesteps_inside_unif + inside_unifs_TF
		} # END for (i in 1:length(timesteps))
	timesteps_inside_a_unif_TF = num_timesteps_inside_unif > 0
	cbind(num_timesteps_inside_unif, timesteps_inside_a_unif_TF)
	plot(timesteps[1:2000], num_timesteps_inside_unif[1:2000], pch=".")
	lines(timesteps[1:2000], num_timesteps_inside_unif[1:2000])
	
	counts_bybin = cbind(timesteps[1:(length(timesteps)-1)], timesteps[2:length(timesteps)], rep(0, (length(timesteps)-1)), rep(1, (length(timesteps)-1)), timesteps_inside_a_unif_TF[1:(length(timesteps)-1)])
	counts_bybin = as.data.frame(counts_bybin)
	names(counts_bybin) = c("tops", "bots", "gapsize", "count", "insideTF")
	head(counts_bybin)

	# Collapse the gaps
	counts_bybin2 = NULL
	for (i in 1:nrow(counts_bybin))
		{
		if (i == 1)
			{
			state1 = counts_bybin$insideTF[i]
			counts_bybin2 = counts_bybin[i,]
			j = 1
			next()
			}
		state2 = counts_bybin$insideTF[i]
		
		# Add the rows if the same,
		# Add new row if different
		if (state2 == state1)
			{
			counts_bybin2$count[j] = counts_bybin2$count[j] + counts_bybin$count[i]
			} else {
			counts_bybin2 = rbind(counts_bybin2, counts_bybin[i,])
			j = j+1
			}
		state1 = state2
		} # END for (i in 1:nrow(counts_bybin))
	
	# Recalculate gapsize
	counts_bybin2$bots[1:(nrow(counts_bybin2)-1)] = counts_bybin2$tops[2:(nrow(counts_bybin2))]
	counts_bybin2$bots[nrow(counts_bybin2)] = bot_of_bins
	counts_bybin2$gapsize = counts_bybin2$bots - counts_bybin2$tops
	counts_bybin2
	
	# Re-do counts, based on tipdates
	for (i in 1:nrow(counts_bybin2))
		{
		older_than_top_TF = tipdates > counts_bybin2$tops[i]
		younger_than_bot_TF = tipdates <= counts_bybin2$bots[i]
		inside_TF = (older_than_top_TF + younger_than_bot_TF) == 2
		counts_bybin2$count[i] = sum(inside_TF)
		} # END for (i in 1:nrow(counts_bybin2))
	counts_bybin2
	
	return(counts_bybin2)
	} # END find_sampling_gaps_from_tipdate_ranges




# Fine the gaps amongst the fossil tipdates
find_bins_sampled_from_tipdate_ranges <- function(OTUs_df, bintops=seq(0,200,10), bot_of_bins=1000, collapse=FALSE, psiSamplingStopsAt=0.01, min_precision=0.1)
	{
	defaults='
	bintops=seq(0,200,by=10)
	bot_of_bins=1000
	psiSamplingStopsAt=0.1
	collapse = FALSE
	collapse = "count_tipdates"
	collapse = "count_overlap"
	collapse = FALSE
	'
	
	# Get tipdates (midpoints, means, or starting values)
	tipdates_orig = OTUs_df$tipdate
	
	# Remove 0s (actually, anything close enough to zero to be
	# caught in the "tolerance" value)
	fossil_TF = tipdates_orig >= psiSamplingStopsAt
	# Add psiSamplingStopsAt as the top bin (before 0, added at the end)
	tipdates = sort(tipdates_orig[fossil_TF])

	# Tipdate ranges
	# Convert non-uniform tipdates to uniform
	unifOTUs_df = convert_nonUniform_dates_to_uniform(OTUs_df[fossil_TF, ], CI=0.999)
	unifOTUs_df

	# Fixed fossil tipdates
	fixed_fossils_TF = unifOTUs_df$distribution == "fixed"
	unifOTUs_df$param1[fixed_fossils_TF] = unifOTUs_df$tipdate[fixed_fossils_TF] - min_precision/2
	unifOTUs_df$param2[fixed_fossils_TF] = unifOTUs_df$tipdate[fixed_fossils_TF] + min_precision/2

	# Tipdate range tops
	tipdate_rangetops = unifOTUs_df$param1
	tipdate_rangebots = unifOTUs_df$param2
	cbind(tipdate_rangetops, tipdate_rangebots)

	# Gaps with/without sampling
	timesteps = seq(0, bot_of_bins, min_precision)

	timestep_below_unif_rangetops_TF = rep(NA, length(tipdate_rangetops))
	timestep_above_unif_rangebots_TF = rep(NA, length(tipdate_rangebots))
	num_timesteps_inside_unif = rep(0, length(timesteps))
	for (i in 1:length(tipdate_rangetops))
		{
		timesteps_below_unif_rangetops_TF = timesteps > tipdate_rangetops[i]
		timesteps_above_unif_rangebots_TF = timesteps <= tipdate_rangebots[i]
		inside_unifs_TF = (timesteps_below_unif_rangetops_TF + timesteps_above_unif_rangebots_TF) == 2
		num_timesteps_inside_unif = num_timesteps_inside_unif + inside_unifs_TF
		} # END for (i in 1:length(timesteps))
	timesteps_inside_a_unif_TF = num_timesteps_inside_unif > 0

	
	# Sampling bins
	bintops = c(bintops, bot_of_bins)
	tops = bintops[1:(length(bintops)-1)]
	tops[1] = psiSamplingStopsAt
	bots = bintops[2:(length(bintops))]
	gapsize = bots-tops
	count_tipdates = rep(0, length(tops))
	count_overlap = rep(0, length(tops))
	counts_bybin = as.data.frame(cbind(tops, bots, gapsize, count_tipdates, count_overlap))
	names(counts_bybin) = c("tops", "bots", "gapsize", "count_tipdates", "count_overlap")
	
	
	
	
	# Number of fossil ranges overlapping with each bin
	for (i in 1:nrow(counts_bybin))
		{
		top = counts_bybin$tops[i]
		bot = counts_bybin$bots[i]

		timesteps_below_sampling_bin_top_TF = timesteps > top
		timesteps_above_sampling_bin_bot_TF = timesteps <= bot
		timesteps_inside_sampling_bin_TF = (timesteps_below_sampling_bin_top_TF + timesteps_above_sampling_bin_bot_TF) == 2
		
		timesteps_inside_both_sampling_bin_and_unif_TF = (timesteps_inside_a_unif_TF + timesteps_inside_sampling_bin_TF) == 2
		if (sum(timesteps_inside_both_sampling_bin_and_unif_TF) > 0)
			{
			sampling_bin_contains_unif_TF = TRUE
			} else {
			sampling_bin_contains_unif_TF = FALSE
			}
		
		# Count the number of tipdates
		tipdates_below_sampling_bin_top_TF = tipdates > top
		tipdates_above_sampling_bin_bot_TF = tipdates <= bot
		tipdates_in_sampling_bin_TF = (tipdates_below_sampling_bin_top_TF + tipdates_above_sampling_bin_bot_TF) == 2
		counts_bybin$count_tipdates[i] = sum(tipdates_in_sampling_bin_TF)
		
		# Count the number of overlaps
		if (sampling_bin_contains_unif_TF == TRUE)
			{
			counts_bybin$count_overlap[i] = max(num_timesteps_inside_unif[timesteps_inside_both_sampling_bin_and_unif_TF])
			} else {
			counts_bybin$count_overlap[i] = 0
			} # END if (sampling_bin_contains_unif_TF == TRUE)
		} # END for (i in 1:nrow(counts_bybin))
	counts_bybin

	if (collapse == "count_overlap")
		{
		TF = counts_bybin$count_overlap > 0
		counts_bybin$count = counts_bybin$count_overlap		
		}
	if (collapse == "count_tipdates")
		{
		TF = counts_bybin$count_tipdates > 0
		counts_bybin$count = counts_bybin$count_tipdates
		}

	if (collapse != FALSE)
		{
		counts_bybin2 = NULL
		
		for (i in 1:nrow(counts_bybin))
			{
			if (i == 1)
				{
				state1 = TF[i]
				counts_bybin2 = counts_bybin[i,]
				j = 1
				next()
				}
			state2 = TF[i]
			
			# Add the rows if the same,
			# Add new row if different
			if (state2 == state1)
				{
				counts_bybin2$count[j] = counts_bybin2$count[j] + counts_bybin$count[i]
				} else {
				counts_bybin2 = rbind(counts_bybin2, counts_bybin[i,])
				j = j+1
				}
			state1 = state2
			} # END for (i in 1:nrow(counts_bybin))
		
		# Recalculate gapsize
		counts_bybin2$bots[1:(nrow(counts_bybin2)-1)] = counts_bybin2$tops[2:(nrow(counts_bybin2))]
		counts_bybin2$bots[nrow(counts_bybin2)] = bot_of_bins
		
		counts_bybin2$gapsize = counts_bybin2$bots - counts_bybin2$tops
		
		
		# If the collapsed table ends up shorter, re-do counts, 
		# this time DON'T collapse
		if (nrow(counts_bybin2) < nrow(counts_bybin))
			{
			counts_bybin2 = find_bins_sampled_from_tipdate_ranges(OTUs_df, bintops=counts_bybin2$tops, bot_of_bins=bot_of_bins, collapse=FALSE, psiSamplingStopsAt=psiSamplingStopsAt, min_precision=min_precision)
			} # END if (nrow(counts_bybin2) < nrow(counts_bybin))
		
		# The "count_overlap" column is really a sum( max within each sub-bin )
		
		return(counts_bybin2)
		} # END if (collapse == TRUE)
	
	return(counts_bybin)
	}


