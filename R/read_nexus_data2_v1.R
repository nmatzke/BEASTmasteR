#######################################################
# read_nexus_data2: Reading in NEXUS data files that
# have e.g. morphology data
#
# Requirements:
# 1. NEXUS data files should be "simplified NEXUS" format, 
# e.g. by exporting data matrices from Mesquite as "simplified NEXUS"
#
# 2. Taxon names should have no spaces. Use "_" instead.
#
# 3. Taxon names should have no ' or " marks.
#
# 4. Sometimes, if you manually remove ' or ", there is no space
#    between the taxon name and the start of the sequence data.
#    Make sure there is whitespace there.
# 
#######################################################



# For Debugging NEXUS input -- see TNT R stuff 
# for better functions.
defaults='
file=nexus_fn
check_ambig_chars=TRUE

file=infn
check_ambig_chars=TRUE
convert_ambiguous_to=NULL
printall="short"
convert_ambiguous_to_IUPAC=FALSE
'

# read_nexus_data2 is based on APE's "read.data.nexus", but has been 
# substantially modified to take into account morphology data, common 
# ambiguities, etc.

read_nexus_data2 <- function(file, check_ambig_chars=TRUE, convert_ambiguous_to=NULL, printall="short", convert_ambiguous_to_IUPAC=FALSE) 
	{
	
	
	if (printall != "none")
		{
		cat("\n\nRunning read_nexus_data2(). This function modifies APE's read.nexus.data() in order to\nsuccessfully parse ambiguous morphological characters and other issues.  It might not \nwork perfectly, as NEXUS is a very complex format. Also, just because it reads the file\ndoes not mean that functions in other packages will magically be able to handle\nyour weird morphology data.\n")

		cat("\nPlease read TIPS (below). If you can't figure it out, or find a bug, you can ask the author, Nick Matzke, at the BEASTmasteR Google Group\n")
		
		cat("\nTIPS:")
		cat("\n- You should use Mesquite to export your data to 'simplified\n  NEXUS' format before using this function.")
		cat("\n- Taxon/OTU names should include *NO* spaces or ' characters!")		
		cat("\n- Filenames should also have *NO* spaces or ' characters!")		
		cat("\n- If convert_ambiguous_to==TRUE (default), the script will check for\n  spaces in the character rows in the data matrix, e.g. due to '(0 1)'.", sep="")
		cat("\n- But your names still need to have no spaces, no apostrophes (')...\n\n", sep="")
		cat("\n\n====================================================================\n", sep="")
		cat("\n- DNA data: convert_ambiguous_to_IUPAC=TRUE runs some Perl scripts to\n  attempt to convert non-IUPAC stuff like e.g. {A C}\n  to IUPAC like M\n  See DNA IUPAC codes at: http://www.dna.affrc.go.jp/misc/MPsrch/InfoIUPAC.html\n")
		cat("\n- Those Perl command may fail on e.g. Windows. If so change\n  convert_ambiguous_to_IUPAC=FALSE and make sure your data is IUPAC!\n")	
		cat("\n- Amino acids: *NO* '(' or '{' ambiguities!  You should convert these\n  manually to '?', or figure how how to let Beast2 read these.\n")	
		cat("\n\n====================================================================\n", sep="")
		cat("Processing ", file, "\n", sep="")
		cat("(notes may follow; set printall='none' to turn off)\n", sep="")
		cat("====================================================================\n\n", sep="")
		} # END if (printall != "none")


	#######################################################
	# If desired, use GREP to convert all ambiguous characters to their IUPAC codes, 
	# according to:
	# http://www.boekhoff.info/?pid=data&dat=fasta-codes
	#######################################################
	if (convert_ambiguous_to_IUPAC == TRUE)
		{
		strs = NULL
		sn = 0
		
		prefix = get_fn_prefix(file) 
		suffix = get_fn_suffix(file) 
		new_fn = paste0(prefix, "_noAmbig", ".", suffix)
		file = new_fn
		
		# After running this, you (SHOULDN'T!) need to check ambiguous
		# characters in DNA
		check_ambig_chars = FALSE
		
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ C\\)/M/g;" -pi.bak $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ G\\)/R/g;" -pi.bak $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ T\\)/W/g;" -pi.bak $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(C\\ G\\)/S/g;" -pi.bak $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(C\\ T\\)/Y/g;" -pi.bak $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(G\\ T\\)/K/g;" -pi.bak $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ C\\ G\\)/V/g;" -pi.bak $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ C\\ T\\)/H/g;" -pi.bak $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ G\\ T\\)/D/g;" -pi.bak $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(C\\ G\\ T\\)/B/g;" -pi.bak $(find ', file, ' -type f)', sep="")
		strs[[(sn=sn+1)]] = paste('perl -e "s/\\(A\\ C\\ G\\ T\\)/N/g;" -pi.bak $(find ', file, ' -type f)', sep="")
		
		if (printall != "none")
			{
			cat("\n\n")
			cat("read_nexus_data2() is running perl to convert DNA from non-IUPAC to IUPAC...\ne.g. '(A C)' to 'M'...\n\n")
			} # END if (printall != "none")
		
		for (i in 1:sn)
			{
			if (printall != "none")
				{
				cat(paste("Running: ", strs[[i]], "...\n", sep=""))
				}
			system(strs[[i]])
			}

		} # END if (convert_ambiguous_to_IUPAC == TRUE)
	
	if (printall != 'none')
		{
		cat("\n\nReading in data...\n\n")
		}
	
	# Find the number of taxa
    "find.ntax" <- function(x)
    	{
        for (i in 1:NROW(x))
        	{
            if (any(f <- grep("\\bntax", x[i], ignore.case = TRUE)))
            	{
                ntax <- as.numeric(sub("(.+?)(ntax\\s*\\=\\s*)(\\d+)(.+)", 
                  "\\3", x[i], perl = TRUE, ignore.case = TRUE))
                break
            	} # END if (any(f <- grep("\\bntax", x[i], ignore.case = TRUE)))
        	} # END for (i in 1:NROW(x))
        return(ntax)
    	} # END "find.ntax" <- function(x)

	# Find the number of characters
    "find.nchar" <- function(x)
    	{
        for (i in 1:NROW(x))
        	{
            if (any(f <- grep("\\bnchar", x[i], ignore.case = TRUE)))
            	{
                nchar <- as.numeric(sub("(.+?)(nchar\\s*\\=\\s*)(\\d+)(.+)", 
                  "\\3", x[i], perl = TRUE, ignore.case = TRUE))
                break
				} # END if (any(f <- grep("\\bnchar", x[i], ignore.case = TRUE)))
			} # END  for (i in 1:NROW(x))
		return(nchar)
		}

	# Find the data matrix startline
    "find.matrix.line" <- function(x)
    	{
        for (i in 1:NROW(x))
        	{
            if (any(f <- grep("\\bmatrix\\b", x[i], ignore.case = TRUE)))
            	{
                matrix.line <- as.numeric(i)
                break
				} # END if (any(f <- grep("\\bmatrix\\b", x[i], ignore.case = TRUE)))
			} # END for (i in 1:NROW(x))
        matrix.line
	    } # END "find.matrix.line" <- function(x)
    
    
    "trim.whitespace" <- function(x)
    	{
        gsub("\\s+", "", x)
    	}
    
    "trim.semicolon" <- function(x)
    	{
        gsub(";", "", x)
    	}
    
    # Check if you can find the file
    if (file.access(file, mode = 4))
    	{
        stop("file could not be found")
    	}
    
    # Scan in the datafile
    X <- scan(file=file, what=character(), sep = "\n", quiet = TRUE, 
        comment.char = "[", strip.white = TRUE)
    ntax <- find.ntax(X)
    nchar <- find.nchar(X)
    matrix.line <- find.matrix.line(X)
    start.reading <- matrix.line + 1
    Obj <- list()
    length(Obj) <- ntax
    i <- 1
    pos <- 0
    tot.nchar <- 0
    tot.ntax <- 0
    





	# Some lines that have caused problems in the past.
	bugcheck = '
	Xj = "y2008Wola 2(0 1)002???91011101001002100511?200?1??01?11?0?30030100000100?1??1???10201?10?1?10?1030101201000?10?1???10?1?1?1?????0???1??00&00201?"
	# What if you have ambiguous characters with more than one state?
	Xj = "Gafrarium_tumidum_FMNH307858+FMNH312294+IM200741731+IM200741732+IM200741733 0000(1 2)0000(0 1)001000022110(1 3 4)0(0 2 3)0(0 2)13000-01(0 2)11(0 2)1120-10010101011210000100??????(0 1 2)"
	Xj = "H._ergaster         (0 1)(0 1)(0 1)(0 1)???????00110111111011110?101101??????????????????????????1????1???1??????11??1???1?1?????????????0??????????0????0?(0 1)????11?????????????????????????????????????????????????0??????????????????????????????????????1?0??10(0 1)10?1?00??0??01?????00100010???0001(0 1)110000001011(0 1)(0 1)01001?0(0 1)110?11100????????????0100(0 1)(0 1)0101100001000001001011011000011(0 1)00(0 1)(0 1)01(0 1)000(0 1)(0 1)(0 1)10(0 1)1(0 1)10(0 1)(0 1)(0 1)0112?(1 2)???0211(0 1)0100(0 1 2)(0 1)0(1 2)(1 2)(0 2)2?00013200?22101?????2???????????????????????????????1?21002112???2122????1?1221102?021222(1 2)2121(0 1)1(0 1)21122?(0 1 2)(0 1)(1 2)0(0 1)0(0 1 2)02(0 1)0(0 1 2)(1 2)(1 2)12331021(0 1)0(0 2)13(1 2)11"
	
	# Squiggely brackets:
	Xj = "Byronosaurus_jaffei                 ?????101???101?1100110101011?00??????20120????????????1?100??????0000001??11????00021??01?0?0???????010110????????0??02??????????????????????????????????????????????????????????????1????0?0??0????????1??????????0????0????100???0?{23}0?100?101100??0100???10?00110?????????????????????????????????????????????????????????????11????????????????????0??11?????0??????????00?00??01?1"

	ts <- strsplit(Xj, "[ \t]+")[[1]]
	Name = ts[1]
	
	tmpSeq = list2str_fast_nosplit(ts[-1], spacer=" ")
	'
	
	
	# Go through the characters, check for ambiguities if desired
    if ((check_ambig_chars == TRUE) && (printall != "none"))
    	{
    	cat("\n\nread_nexus_data2() is checking for ambiguous characters in line...\n", sep="")
    	} # END  if ((check_ambig_chars == TRUE) && (printall != "none"))
    
    
    #######################################################
	# LOOP THROUGH ROWS j OF THE NEXUS FILE
	#######################################################
    
    #for (j in start.reading:22)
	for (j in start.reading:NROW(X))
    	{
    	if ( (check_ambig_chars == TRUE) && (printall != "none") )
    		{
		    cat(j, ", ", sep="")
		    }
        Xj <- trim.semicolon(X[j])
        if (Xj == "")
        	{
            break
        	}
        # Stop once you hit "END"
        if (any(jtmp <- grep("\\bend\\b", X[j], perl = TRUE, ignore.case = TRUE)))
            {
            break
	        }
        
        # Convert character string line with e.g. (0 1) to Z (or whatever)
        # count the number of spaces
        # number_of_spaces = count_chars(char_to_count=" ", Xj)
        
	
		
		need_to_correct_ambiguous = FALSE
        if (check_ambig_chars == TRUE)
        	{
        	# split on whitespace
        	ts <- strsplit(Xj, "[ \t]+")[[1]]
        	Name = ts[1]
        	
        	# Take everything that wasn't the name and merge back into old string
        	tmpSeq = list2str_fast_nosplit(ts[-1], spacer=" ")
        	tmpSeq_orig = tmpSeq
        	
        	
        	# Convert any "{" or "}" to "(" and ")"
        	if (grepl(pattern="\\{", x=tmpSeq) == TRUE)
        		{
				tmpSeq = gsub(pattern="\\{", replacement="(", x=tmpSeq)
				tmpSeq = gsub(pattern="\\}", replacement=")", x=tmpSeq)
				}
        	
        	# convert the ambiguous characters to "~" for later processing
        	# locate e.g. (0 1), (0 1 2), etc...
        	ambig_chars_paren_locations = gregexpr('(\\(.*?\\))', tmpSeq, perl=TRUE)[[1]]
			
			# If there ARE ambiguous characters, convert them to ~
			# otherwise, don't
			ambig_chars_TF = FALSE
			if ( length(ambig_chars_paren_locations) > 1)
				{
				ambig_chars_TF = TRUE
				} else {
				if (ambig_chars_paren_locations != -1)
					{
					ambig_chars_TF = TRUE
					} # END if (ambig_chars_paren_locations != -1)
				} # END if ( length(ambig_chars_paren_locations) > 1)
			
			if (ambig_chars_TF == TRUE)
				{
				need_to_correct_ambiguous = TRUE
				ambig_chars_paren_list = NULL
				
				# Count down backwards, to avoid the changing of the sequence length
				for (k in length(ambig_chars_paren_locations):1)
					{
					
					startnum = ambig_chars_paren_locations[k]
					ambig_length = attr(ambig_chars_paren_locations, "match.length")[k]
					endnum = startnum + ambig_length - 1
					tmpSeq_split = strsplit(tmpSeq, split="")[[1]]

					# Do you want to put the full ambiguous coding into the NEXUS
					# object, or replace with e.g. a question mark ("?") ?
					if (is.null(convert_ambiguous_to))
						{
						# Put in original ambiguous coding into NEXUS object
						ambig_charstring = list2str_fast(tmpSeq_split[startnum:endnum])
						ambig_charstring
						
						# Check for ambig_charstring of length 4,
						# e.g. "(23)"
						# Things SHOULD have spaces, but maybe
						# {23} doesn't need them
						if (nchar(ambig_charstring) == 4)
							{
							ambig_charstring_chars = strsplit(ambig_charstring, split="")[[1]]
							ambig_charstring_chars
							ambig_charstring_chars2 = c(ambig_charstring_chars[1:2], " ", ambig_charstring_chars[3:4])
							ambig_charstring = list2str_fast(ambig_charstring_chars2, spacer="")
							}
						
						
						tmprow = c(startnum, endnum, ambig_length, ambig_charstring)
						} else {
						# Put ? into NEXUS object
						#print(Name)
						tmprow = c(startnum, endnum, ambig_length, "?")
						#print(tmprow)
						} # END if (is.null(convert_ambiguous_to))
			
					# ambiguous character = ~
# 					if ((j == 22) && (k==73))
# 						{
# 						print("Sequence before recoding ambiguous")
# 						print(tmpSeq)
# 						#print(ambig_chars_paren_locations)
# 						cat("\n", startnum, endnum, startnum-1, "\n")
# 						}
					
					
					# If you are NOT at the end of the string -- 
					# catches where end of string is e.g. (0 1)
					if (endnum < length(tmpSeq_split))
						{
						# This catches the case where the first character is e.g. (0 1)
						if (startnum-1 > 0)
							{
							newSeq = list2str_fast(c(tmpSeq_split[1:(startnum-1)], "~", tmpSeq_split[(endnum+1):length(tmpSeq_split)] ))
							} else {
							newSeq = list2str_fast(c("~", tmpSeq_split[(endnum+1):length(tmpSeq_split)] ))
							}
						} else {
						# If you ARE at the end of the string: -- catches where end of string is e.g. (0 1)
						# This catches the case where the first character is e.g. (0 1)
						if (startnum-1 > 0)
							{
							newSeq = list2str_fast(c(tmpSeq_split[1:(startnum-1)], "~") )
							} else {
							newSeq = list2str_fast(c("", "~") )
							}
						} # END if (endnum < length(tmpSeq_split))
	
					# add to the list of ambiguous characters
					ambig_chars_paren_list = rbind(ambig_chars_paren_list, tmprow)
					
					tmpSeq = newSeq
					##############################
					# END for-loop for characters within this taxon
					##############################
					} # END for (k in length(ambig_chars_paren_locations):1)
				
				# We will have to store this in a permanent place somewhere
				ambig_chars_paren_list = as.data.frame(as.matrix(ambig_chars_paren_list), row.names=NULL, stringsAsFactors=FALSE)
				row.names(ambig_chars_paren_list) = NULL
				names(ambig_chars_paren_list) = c("startnum", "endnum", "ambig_length", "ambig_charstring")
        		} # END if (ambig_chars_TF == TRUE)
        	
        	        	
        	# locate e.g. [0 1], [0 1 2], etc...
        	#ambig_locations2 = gregexpr('(\\[(.*?)\\])', tmpSeq, perl=TRUE)[[1]]

        	# locate e.g. 0&1, etc...
        	#ambig_locations3 = gregexpr('((.&.))', tmpSeq, perl=TRUE)[[1]]
        	
        	# store the revised sequence in Seq
        	Seq = tmpSeq
        	# endif if (check_ambig_chars == TRUE)
        	
        	} else {
        	# So: if (check_ambig_chars == FALSE)
			ts <- unlist(strsplit(Xj, "(?<=\\S)(\\s+)(?=\\S)", perl = TRUE))
			if (length(ts) > 2)
				{
				cat("\n\nERROR: length of ts (data string, name-space-characters) is >2.  printing ts:\n", sep="")
				print(ts)
				stop("\nOriginal error message: nexus parser does not handle spaces in sequences or taxon names (ts>2)")
				} # END if (length(ts) > 2)
			if (length(ts) != 2)
				{
				cat("\n\nERROR: length of ts (data string, name-space-characters) is !=2.  printing ts:\n", sep="")
				print(ts)
				stop("\nOriginal error message: nexus parser failed to read the sequences (ts!=2)")
				} # END if (length(ts) != 2)
			Seq <- trim.whitespace(ts[2])
			Name <- trim.whitespace(ts[1])
			} # END if (check_ambig_chars == TRUE)
		
		# Get the name so that you can 
		# extract it from the line of the NEXUS file
        nAME <- paste(c("\\b", Name, "\\b"), collapse = "")
		
		# Here, this if() didnt run on a standard simplified MrBayes NEXUS file
		l <- grep(nAME, names(Obj))
		#print("Printing l")
		#print(l)
		
		if (any(l))
        	{
            tsp <- strsplit(Seq, NULL)[[1]]
            tsp_orig = tsp
            
			# Correct for ambiguous sequences
            if (need_to_correct_ambiguous == TRUE)
            	{
            	# Add rev() to as ambig chars were recognized in reverse order!
				tsp[tsp=="~"] = rev(ambig_chars_paren_list$ambig_charstring)
				
				if ((j == 22) && (k==73))
					{
					ambig_chars_paren_list
					}
				
				
				# Tell the user what happened if ambiguous characters are detected
				if (is.null(convert_ambiguous_to))
					{
					if (printall != "none")
						{
						if (printall == "short")
							{
							cat("\n", j, ": read_nexus_data2() detected & stored ambiguous characters in: ", Name, "\n", sep="")
							} else {
							cat("\n", j, ": read_nexus_data2() detected & stored ambiguous characters in: ", Name, "\n", tmpSeq_orig, "\n", sep="")
							} # END if (printall == "short")
						} # END if (printall != "none")
					} else {
					if (printall != "none")
						{
						if (printall == "short")
							{
							cat("\n", j, ": read_nexus_data2() converted ambiguous characters in: '", Name, "' to '", convert_ambiguous_to, "'\n", sep="")
							} else {
							cat("\n", j, ": read_nexus_data2() converted ambiguous characters in: '", Name, "' to '", convert_ambiguous_to, "':", "\n", tmpSeq_orig, " --> \n", sep="")
							cat(list2str_fast_nosplit(tsp), "\n", sep="")				
							} # END if (printall == "short")
						} # END if (printall != "none")
					} # END if (is.null(convert_ambiguous_to))
				} # END if (need_to_correct_ambiguous == TRUE)
            
            # input the characters into the object
            for (k2 in 1:length(tsp))
            	{
                p <- k2 + pos
                try_result = try(Obj[[l]][p] <- tsp[k2])
                if (class(try_result) == "try-error")
                	{
                 	cat("\n\nA parsing error happened at read_nexus_data2, position 1.\n\n")
                	txt = "ERROR in parsing detected by read_nexus_data2(): printing j, Name, nAME, l, p, k2, pos"
                	cat("\n\n")
                	cat(txt)
                	cat("\n")
                	print(j)
                	cat("\n\n")
                	print("Name:")
                	print(Name)
                	cat("\n\n")
                	print("nAME:")
                	print(nAME)
                	cat("\n\n")
                	txt2 = "STOP ERROR NOTE: If you see characters (01011 etc..) where 'Name' or 'nAME' are, this probably means you have INAPPROPRIATE LINE BREAKS in your NEXUS file. R's scan() function reads files line-by-line, so each row of the data matrix HAS TO BE ON A SINGLE LINE.\n\n"
                	cat(txt2)

                	#print(Obj)
                	print(l)
                	print(p)
                	print(k2)
                	print(pos)
                	#print(tsp[[k2]])
                	cat("\n\n")
                	txt3 = paste0(txt, "\n\n", txt2)
                	stop(txt3)
                	}
                chars.done <- k2
	            }
	        # endif if (any(l))
    	    } else {
			# Here, this else() DID run on a standard simplified MrBayes NEXUS file
			tsp <- strsplit(Seq, NULL)[[1]]
            names(Obj)[i] <- Name
			tsp_orig = tsp
			l <- grep(nAME, names(Obj))	# Needed 2015-10-08


			# Correct for ambiguous sequences
            if (need_to_correct_ambiguous == TRUE)
            	{
				# Add rev() to as ambig chars were recognized in reverse order!
				tsp[tsp=="~"] = rev(ambig_chars_paren_list$ambig_charstring)

				# Tell the user what happened if ambiguous characters are detected
				if (is.null(convert_ambiguous_to))
					{
				    if (printall != "none")
    					{
						if (printall == "short")
							{
							cat("\n", j, ": read_nexus_data2() detected & stored ambiguous characters in: ", Name, "\n", sep="")
							} else {
							cat("\n", j, ": read_nexus_data2() detected & stored ambiguous characters in: ", Name, "\n", tmpSeq_orig, "\n", sep="")
							}
						} # END if (printall != "none")
					} else {
					if (printall == "short")
						{
						if (printall == "short")
							{
							cat("\n", j, ": read_nexus_data2() converted ambiguous characters in: '", Name, "' to '", convert_ambiguous_to, "'\n", sep="")
							} else {
							cat("\n", j, ": read_nexus_data2() converted ambiguous characters in: '", Name, "' to '", convert_ambiguous_to, "':", "\n", tmpSeq_orig, " --> \n", sep="")
							cat(list2str_fast_nosplit(tsp), "\n", sep="")				
							}
						} # END if (printall != "none")
					} # END if (is.null(convert_ambiguous_to))
				} # END if (need_to_correct_ambiguous == TRUE)
            
            for (k2 in 1:length(tsp))
            	{
                p <- k2 + pos
                Obj[[i]][p] <- tsp[k2]
                try_result = try(Obj[[l]][p] <- tsp[k2])
                if (class(try_result) == "try-error")
                	{
                 	cat("\n\nA parsing error happened at read_nexus_data2, position 2.\n\n")
                	txt = "ERROR in parsing detected by read_nexus_data2(): printing j, Name, nAME, l, p, k2, pos"
                	cat("\n\n")
                	cat(txt)
                	cat("\n")
                	print(j)
                	cat("\n\n")
                	print("Name:")
                	print(Name)
                	cat("\n\n")
                	print("nAME:")
                	print(nAME)
                	cat("\n\n")
                	txt2 = "STOP ERROR NOTE: If you see characters (01011 etc..) where 'Name' or 'nAME' are, this probably means you have INAPPROPRIATE LINE BREAKS in your NEXUS file. R's scan() function reads files line-by-line, so each row of the data matrix HAS TO BE ON A SINGLE LINE.\n\n"
                	cat(txt2)

                	#print(Obj)
                	print(l)
                	print(p)
                	print(k2)
                	print(pos)
                	#print(tsp[[k2]])
                	cat("\n\n")
                	txt3 = paste0(txt, "\n\n", txt2)
                	stop(txt3)
                	}
                chars.done <- k2
            	}
	        }
        
        tot.ntax <- tot.ntax + 1
        if (tot.ntax == ntax)
        	{
            i <- 1
            tot.ntax <- 0
            tot.nchar <- tot.nchar + chars.done
            
            # Make sure the total number of characters = 
            # = total # of taxa * # of characters per taxon
            if (tot.nchar == nchar * ntax)
            	{
                print("Warning: ntot was more than nchar*ntax")
                break
				} # END if (tot.nchar == nchar * ntax)
				pos <- tot.nchar
			} else {
            i <- i + 1
			} # END if (tot.ntax == ntax)
		##################################
		# ENDING THE BIG FOR-LOOP
		##################################
		} # END for (j in start.reading:NROW(X))

	# Obj
	# ambig_chars_paren_list
	# ambig_chars_paren_list$ambig_charstring
	# tsp_orig[tsp_orig == "~"]
	# cbind(ambig_chars_paren_list, tsp_orig[tsp_orig == "~"])
	# cbind(ambig_chars_paren_list, tsp[tsp_orig == "~"])

    if (tot.ntax != 0)
    	{
        cat("ntax:", ntax, "differ from actual number of taxa in file?\n")
        stop("read_nexus_data2 parser did not read names correctly (tot.ntax!=0)")
    	} # END if (tot.ntax != 0)
    
    for (i in 1:length(Obj))
    	{
        if (length(Obj[[i]]) != nchar)
        	{
            cat(names(Obj[i]), "has", length(Obj[[i]]), "characters\n")
            stop("nchar differ from sequence length (length(Obj[[i]])!=nchar)")
        	} # END if (length(Obj[[i]]) != nchar)
	    } # END for (i in 1:length(Obj))
    Obj <- lapply(Obj, tolower)
    
    if (printall != "none")
    	{
		cat("\n...read_nexus_data2() done.\n\n")
		} # END if (printall != "none")
    return(Obj)
	} # END read_nexus_data2





# Slight modification of the APE function, to allow the output of characters
# listed as DATATYPE=STANDARD
# (defaults allow only "DNA" and "PROTEIN")
write_nexus_data2 <- function(x, file, format = "dna", datablock = TRUE, interleaved = TRUE, charsperline = NULL, gap = NULL, missing = NULL) 
	{
	# Convert the "format" string to lowercase
	format = tolower(format)
	
    indent <- "  "
    maxtax <- 5
    defcharsperline <- 80
    defgap <- "-"
    defmissing <- "?"
    ntax <- length(x)
    nchars <- length(x[[1]])
    zz <- file(file, "w")
    if (is.null(names(x))) {
        names(x) <- as.character(1:ntax)
    }
    "fcat" <- function(..., file = zz) {
        cat(..., file = file, sep = "", append = TRUE)
    }
    "find.max.length" <- function(x) {
        max <- 0
        for (i in 1:length(x)) {
            val <- length((strsplit(x[i], split = NULL))[[1]])
            if (val > max) {
                max <- val
            }
        }
        max
    }
    "print.matrix" <- function(x, dindent = "    ") {
        Names <- names(x)
        printlength <- find.max.length(Names) + 2
        if (interleaved == FALSE) {
            for (i in 1:length(x)) {
                sequence <- paste(x[[i]], collapse = "")
                taxon <- Names[i]
                thestring <- sprintf("%-*s%s%s", printlength, 
                  taxon, dindent, sequence)
                fcat(indent, indent, thestring, "\n")
            }
        }
        else {
            ntimes <- ceiling(nchars/charsperline)
            start <- 1
            end <- charsperline
            for (j in 1:ntimes) {
                for (i in 1:length(x)) {
                  sequence <- paste(x[[i]][start:end], collapse = "")
                  taxon <- Names[i]
                  thestring <- sprintf("%-*s%s%s", printlength, 
                    taxon, dindent, sequence)
                  fcat(indent, indent, thestring, "\n")
                }
                if (j < ntimes) {
                  fcat("\n")
                }
                start <- start + charsperline
                end <- end + charsperline
                if (end > nchars) {
                  end <- nchars
                }
            }
        }
    }
    fcat("#NEXUS\n[Data written by write_nexus_data2.R,", " ", 
        date(), "]\n")
    NCHAR <- paste("NCHAR=", nchars, sep = "")
    NTAX <- paste("NTAX=", ntax, sep = "")
    if (format == "dna")
    	{
        DATATYPE <- "DATATYPE=DNA"
	    }
    if (format == "protein")
    	{
        DATATYPE <- "DATATYPE=PROTEIN"
	    }
    if (format == "standard")
    	{
        DATATYPE <- "DATATYPE=STANDARD"
	    }
    if (format == "restriction")
    	{
        DATATYPE <- "DATATYPE=RESTRICTION"
	    }
    if (is.null(charsperline)) {
        if (nchars < defcharsperline) {
            charsperline <- nchars
            interleaved <- FALSE
        }
        else {
            if (nchars > defcharsperline) {
                charsperline <- defcharsperline
            }
        }
    }
    if (is.null(missing)) {
        MISSING <- paste("MISSING=", defmissing, sep = "")
    }
    else {
        MISSING <- paste("MISSING=", missing, sep = "")
    }
    if (is.null(gap)) {
        GAP <- paste("GAP=", defgap, sep = "")
    }
    else {
        GAP <- paste("GAP=", gap, sep = "")
    }
    if (interleaved == TRUE) {
        INTERLEAVE <- "INTERLEAVE=YES"
    }
    if (interleaved == FALSE) {
        INTERLEAVE <- "INTERLEAVE=NO"
    }
    if (datablock == TRUE) {
        fcat("BEGIN DATA;\n")
        fcat(indent, "DIMENSIONS", " ", NTAX, " ", NCHAR, ";\n")
        if (format %in% c("dna", "protein", "standard", "restriction"))
        	{
            fcat(indent, "FORMAT", " ", DATATYPE, " ", MISSING, 
                " ", GAP, " ", INTERLEAVE, ";\n")
	        }
        fcat(indent, "MATRIX\n")
        print.matrix(x)
        fcat(indent, ";\n")
        fcat("END;\n\n")
    }
    else {
        fcat("BEGIN TAXA;\n")
        fcat(indent, "DIMENSIONS", " ", NTAX, ";\n")
        fcat(indent, "TAXLABELS\n")
        fcat(indent, indent)
        j <- 0
        for (i in 1:ntax) {
            fcat(names(x[i]), " ")
            j <- j + 1
            if (i == ntax) {
                fcat("\n", indent, ";\n")
            }
            else {
                if (j == maxtax) {
                  fcat("\n", indent, indent)
                  j <- 0
                }
            }
        }
        fcat("END;\n\n")
        fcat("BEGIN CHARACTERS;\n")
        fcat(indent, "DIMENSIONS", " ", NCHAR, ";\n")
        if (format %in% c("dna", "protein", "standard", "restriction"))
        	{
            fcat(indent, "FORMAT", " ", MISSING, " ", GAP, " ", 
                DATATYPE, " ", INTERLEAVE, ";\n")
	        }
        fcat(indent, "MATRIX\n")
        print.matrix(x)
        fcat(indent, ";")
        fcat("\nEND;\n\n")
    }
    close(zz)
	}







# Pull out and code observed ambiguous states for each collection of
# observed character states
# for morphology only
# Hannah's data: 69 morphological characters
# return_missing_chars = "list", "correct", or "none"
# (for listing the sites with missing character states, correcting them, or neither)
get_numstates_per_char <- function(nexd, ambig_to_remove=c("\\(", "\\)", " ", ","), return_missing_chars="list", printall="short")
	{
	defaults='
	ambig_to_remove=c("\\(", "\\)", " ", ",")
	return_missing_chars="list"
	printall="short"
	'
	
	if (is.null(nexd))
		{
		txt = "STOP ERROR in get_numstates_per_char(): the input 'nexd' is NULL"
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if (is.null(nexd))
	
	
	# Make a list of all ? (allQs)
	ntaxa = ncol(nexd)
	allQs_TF = rep(FALSE, times=nrow(nexd))
	invariant_TF = rep(FALSE, times=nrow(nexd))
	
	# Convert data to data.frame
	#cat("\n\nPrinting nexd:\n\n")
	#print(nexd)
	
	nexdf = as.data.frame(x=nexd, row.names=NULL, stringsAsFactors=FALSE)
	#cat("\n\nPrinting nexdf:\n\n")
	#print(nexdf)

	#column_strings_morph = NULL
	numstates_morph_list = NULL

	if (printall != "none")
		{
		cat("\n\nget_numstates_per_char() is checking for missing charstates:\n", sep="")
		} # END if (printall != "none")
	
	# Initialize storage variables
	allQs_characters_list = NULL
	invariant_characters_list = NULL
	missing_charstates_list = NULL
	tmp_charstate_codes = charstate_codes()
	missing_charstates = FALSE
	
	
	# Go through each column, count the number of states
	# (but remove "(", ")" )
	cat("\n\nCounting the number of states in each column:\n")
	for (rownum in 1:nrow(nexdf))
		{
		cat(rownum, " ", sep="")
		
		#print(rownum)
		#print(nexdf[rownum,])
		
		# Check for all "?" characters
		if (countQs(nexdf[rownum,]) == ntaxa)
			{
			if (printall != "none")
				{
				cat("\n\nWARNING: Data column #", rownum, " is ALL ? CHARACTERS, i.e. all missing data. PLEASE CORRECT.\n", sep="")
				cat("  (...adding to list 'allQs_colnums'...)   \n", sep="")
				} # END if (printall != "none")
			# Add to list of allQs rows
			allQs_TF[rownum] = TRUE
			allQs_characters_list = c(allQs_characters_list, rownum)
			}
		
		
		# Is the column a "standard" character?
		# (contains some digits, i.e. "0s" (or 1s!) == "\\d" covers all digits)
		standard_TF = sum(grepl("\\d", nexdf[rownum,])) > 0
		#print(standard_TF)
		
		#cat("\n")
		if (standard_TF != TRUE)
			{
			warning_txt = paste0("\nWARNING from get_numstates_per_char(): rownum '", rownum, "' seems to contain no digits, i.e. it might not be entirely of type 'standard': ", rownum, "\n")
			cat(warning_txt)
			cat("\nPrinting the row:\n\n", sep="")
			print(nexdf[rownum,])
			numstates = 0
			numstates_morph_list = c(numstates_morph_list, numstates)
			
			} else {
			# are there any ambiguous characters?
			ambiguous_TF = grepl("\\(", nexdf[rownum,])
			
			# If there are ambiguous characters, convert them 
			# from e.g. (0 1) to 01
			if (sum(ambiguous_TF) > 0)
				{
				ambig_chars = nexdf[rownum,][ambiguous_TF]
				
				for (i in 1:length(ambig_to_remove))
					{
					# replace "(", ")", " "
					#gsub("\\(", "", c("(0 1)", "(0 1)", "(0 1)"))
					ambig_chars = gsub(ambig_to_remove[i], "", ambig_chars)
					ambig_chars = gsub(ambig_to_remove[i], "", ambig_chars)
					ambig_chars = gsub(ambig_to_remove[i], "", ambig_chars)
					}
				
				nexdf[rownum,][ambiguous_TF] = ambig_chars
				}
			# add these to the list
			unique_characters = unique(strsplit(list2str_fast(nexdf[rownum,]), split="")[[1]])
			
			# drop "-" and "?"
			unique_characters = unique_characters[unique_characters != "-"]
			unique_characters = unique_characters[unique_characters != "?"]
			numstates = length(unique_characters)
			
			# Check for characters missing states (the highest character
			# should equal the one you get from counting up states)
			unique_characters = unique_characters[order(unique_characters)]
			
			# if (rownum == 353)
			# 	{
			# 	print("row 353")
			# 	print(nexd[rownum,])
			# 	print(nexdf[rownum,])
			# 	print(unique_characters)
			# 	print(numstates)
			# 	
			# 	print(strsplit(list2str_fast(nexdf[rownum,]), split="")[[1]])
			# 	stop()
			# 	}
			
			
			# If more than one character state, check for
			# missing character states in a multistate character
			if ( (numstates > 0) && (unique_characters[length(unique_characters)] != tmp_charstate_codes[length(unique_characters)]) )
				{

				# Note characters with only 1 state; and make sure they convert!
				if (numstates == 1)
					{
					# Is the character absolutely invariant, or not?
					TF = nexdf[rownum,] == unique_characters
					if (sum(TF) == length(TF))
						{
						# Absolutely invariant character
						txt1 = ""
					
						} else {
						# Absolutely invariant character, except for 
						# question marks / dashes
						txt1 = ", that is not '?' or '-'"
					
						} # END if (sum(TF) == length(TF))
				
				
					if (printall != "none")
						{
						cat("\n\nWarning from get_numstates_per_char(): in rownum #", rownum, " there is only one character state, '", unique_characters, "'", txt1, ". This can sometimes happen e.g. if you have cut some taxa.  You should probably revise your dataset or exclude this character.", sep="")
						}

					# Add to invariant characters list
					invariant_TF[rownum] = TRUE
					invariant_characters_list = c(invariant_characters_list, rownum)
				
					# Add to list of characters missing implied character states
					missing_charstates_list = c(missing_charstates_list, rownum)
					} # END if (numstates == 1)


				if (printall != "none")
					{
					cat("\n\nWarning from get_numstates_per_char(): in rownum #", rownum, " the maximum charstate is\nnot the same as the maximum charstate derived by counting up charstates...\nmax charstate: ", unique_characters[length(unique_characters)], "\nmax code: ", tmp_charstate_codes[length(unique_characters)], "\n", sep="")
					print(unique_characters)
					#print(tmp_charstate_codes)
					cat("Old row #", rownum, "\n", sep="")
					print(as.character(nexdf[rownum,]))
					} # END if (printall != "none")

				
				missing_charstates_list = c(missing_charstates_list, rownum)
				missing_charstates = TRUE
				if (printall != "none")
					{
					cat(rownum, ", ", sep="")
					} # END if (printall != "none")

				
				# If you want to correct these, do it here
				if (return_missing_chars == "correct")
					{
					#cat("correcting, ", sep="")
					TFs_to_change_list = NULL
					TFs_to_change_num = 0
					unique_chars_to_change_list = NULL
					new_character_values_list = NULL
					
					for (j in 1:numstates)
						{
						unique_char = unique_characters[j]
						if (unique_char != tmp_charstate_codes[j])
							{
							unique_chars_to_change_list = c(unique_chars_to_change_list, unique_char)
							new_character_values_list = c(new_character_values_list, tmp_charstate_codes[j])
							TFs_to_change_list[[(TFs_to_change_num=TFs_to_change_num+1)]] = grepl(unique_char, nexdf[rownum,])
							} # END if (unique_char != tmp_charstate_codes[j])
						} # END for (j in 1:numstates)
					# Now fix these, ONLY in the appropriate cells
					
					if (printall != "none")
						{
						cat("\n...correcting: ", rownum, "...\n", sep="")
						} # END if (printall != "none")
					for (j in 1:length(unique_chars_to_change_list))
						{
						unique_char = unique_chars_to_change_list[j]
						TF_to_change = TFs_to_change_list[[j]] 
						nexdf[rownum, ][TF_to_change] = gsub(unique_char, new_character_values_list[j], nexdf[rownum,][TF_to_change])
						} # END for (j in 1:length(unique_chars_to_change_list))
					if (printall != "none")
						{
						cat("New row #", rownum, "\n", sep="")
						} # END if (printall != "none")

					if (printall != "none")
						{
						print(as.character(nexdf[rownum,]))
						} # END if (printall != "none")
					} # END if (return_missing_chars == "correct")
				} # END if ( (numstates > 1) && (unique_characters...
			numstates_morph_list = c(numstates_morph_list, numstates)			
			} # END if (standard_TF != TRUE)
		# end forloop
		} #END for (rownum in 1:nrow(nexdf))
	cat("\n...done counting the number of states in each column.\n")
		
	# If missing character states were detected, 
	# print "...done. No missing character states detected."
	if (missing_charstates == FALSE)
		{
		if (printall != "none")
			{
			cat("\n\nThe function get_numstates_per_char() says: No missing character states detected.\n", sep="")
			cat("\n")
			} # END if (printall != "none")
		} # END if (missing_charstates == FALSE)

	if (sum(invariant_TF) > 0)
		{
		if (printall != "none")
			{
			cat("\n\nWarning from get_numstates_per_char(): Detected ", sum(invariant_TF), " characters that are invariant (except for '?' or '-'). They are: \n", sep="")
			cat(invariant_characters_list, sep=", ")
			cat("\n")
			} # END if (printall != "none")
		} # END if (sum(invariant_TF) > 0)

	if (sum(allQs_TF) > 0)
		{
		if (printall != "none")
			{
			cat("\n\nWarning from get_numstates_per_char(): Detected ", sum(allQs_TF), " characters that are all '?' or '-'. They are: \n", sep="")
			cat(allQs_characters_list, sep=", ")
			cat("\n")
			} # END if (printall != "none")
		} # END if (sum(allQs_TF) > 0)
		
	if (printall != "none")
		{
		cat("\n\nThe function get_numstates_per_char() is finished; returning 'results' object.\n\n")
		}

	# Also return the columns that have missing states, if desired.
	if (return_missing_chars == "list")
		{
		results = NULL
		results$numstates_morph_list = numstates_morph_list
		results$missing_charstates_list = missing_charstates_list
		results$allQs_TF = allQs_TF
		results$invariant_TF = invariant_TF
		results$invariant_characters_list = invariant_characters_list
		return(results)
		}

	if (return_missing_chars == "correct")
		{
		results = NULL
		results$nexdf = nexdf
		results$numstates_morph_list = numstates_morph_list
		results$missing_charstates_list = missing_charstates_list
		results$allQs_TF = allQs_TF
		results$invariant_TF = invariant_TF
		results$invariant_characters_list = invariant_characters_list
		return(results)
		}
	
	retrieve_cmds='
	numstates_morph_list = res$numstates_morph_list
	missing_charstates_list = res$missing_charstates_list
	allQs_TF = res$allQs_TF
	morph_df2_corrected = res$nexdf
	'
			
	return(numstates_morph_list)
	}


# Count question marks
countQs <- function(tmprow)
	{
	numQs = sum(tmprow == "?", na.rm=TRUE)
	return(numQs)
	}



# Ambiguous single-letter codes for morphology
charstate_codes <- function(max_numstates=32)
	{
	# TNT can do 32 character states
	statenums = 1:max_numstates
	
	alphabet = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z")
	
	alphabet = LETTERS
	
	statecodes = c(0:9, LETTERS[1:22])
	return(statecodes)
	}







dtf_to_phylip <- function(char_dtf, outfn="tmpdata.phylip", txtTF=FALSE, max_namelength="max")
	{
	default='
	char_dtf = t(morph_df)
	outfn="tmpdata.phylip"
	txtTF=TRUE
	max_namelength="max"
	multistate_chars_as="paren"
	'
	
	# write a data.frame to phylip
	# you may want to put columns with different numbers of states into different columns!!
	

	
	# convert columns to character
	d = apply(char_dtf, 2, paste)
	d2 = as.data.frame(d)
	tmplevels = as.character(sort(unique(unlist(char_dtf))))
	numcols = ncol(d2)
	
	
	
	str1 = paste(nrow(d2), numcols, sep=" ")
	list_of_strings = NULL
	list_of_strings = c(list_of_strings, str1)
	
	# Maximum numer of characters
	if (max_namelength == "max")
		{
		max_namelength = max(nchar(row.names(char_dtf))) + 1
		}
	
	for (i in 1:nrow(d2))
		{
		tmpname = row.names(char_dtf)[i]
		numchars = nchar(tmpname)
		if(numchars > max_namelength)
			{
			tmpname = strtrim(tmpname, max_namelength)
			numspaces_to_add = 0
			} else {
			numspaces_to_add = max_namelength-numchars
			}
		
		spaces_to_add = paste(rep(" ", numspaces_to_add), collapse="")
		
		charstring = paste(d2[i, ], collapse="")
		outrowstr = paste(tmpname, spaces_to_add, " ", charstring, sep="")
		
		list_of_strings = c(list_of_strings, outrowstr)
		}
	
	
	if (txtTF == FALSE)
		{
		# Start writing to file
		write(list_of_strings[[1]], file=outfn, ncolumns=1, append=FALSE, sep="")
		
		# Continue writing
		for (i in 2:length(list_of_strings))
			{
			outrowstr = list_of_strings[[i]]
			write(outrowstr, file=outfn, ncolumns=1, append=TRUE, sep="")
			} # END for (i in 1:nrow(d2))
		return(outfn)
		} # END if (txtTF == FALSE)
	
	# Otherwise, paste everything together into a string
	seqs_txt = paste(unlist(list_of_strings), collapse="\n", sep="")
	
	return(seqs_txt)
	}


# Get statistics for a morphology matrix
# morphstats = morphology_matrix_stats(charslist, charsdf=NULL)
morphology_matrix_stats <- function(charslist=NULL, charsdf=NULL)
	{
	if (is.null(charsdf))
		{
		charsdf = as.data.frame(charslist, stringsAsFactors=FALSE)
		charsdf
		}
	numchars = nrow(charsdf)
	numtaxa = ncol(charsdf)

	# invariant_characters
	uniq_chars = apply(X=charsdf, MAR=1, FUN=unique)
	num_uniq = NULL
	for (u in 1:length(uniq_chars))
		{
		tmpchars = uniq_chars[[u]]
		TF = grepl(pattern="\\?", x=tmpchars)
		tmpchars = tmpchars[TF == FALSE]
		TF = grepl(pattern="\\(", x=tmpchars)
		tmpchars = tmpchars[TF == FALSE]
		num_uniq = c(num_uniq, length(tmpchars))
		}
	num_uniq
	LT_2states_TF = num_uniq < 2

	numQs = NULL
	numParens = NULL
	nQ_woInvar = NULL
	for (i in 1:numtaxa)
		{
		tmpcol = charsdf[,i]
		paren_TF = grepl(pattern="\\(", x=tmpcol)
		Qs_TF = grepl(pattern="\\?", x=tmpcol)
		Qs_TF
	
		numParens = c(numParens, sum(paren_TF))
		numQs = c(numQs, sum(Qs_TF))
	
		# Convert invariant characters to all Qs, and recount
		#print(sum(Qs_TF))
		tmpcol[LT_2states_TF] = "?"
		Qs_TF_wo_Invar = grepl(pattern="\\?", x=tmpcol)
		#print(sum(Qs_TF_wo_Invar))
	
		nQ_woInvar = c(nQ_woInvar, sum(Qs_TF_wo_Invar))
		}
	nQ_woInvar
	completeness_df = as.data.frame(cbind(names(charsdf), numParens, numQs, nQ_woInvar), stringsAsFactors=FALSE)
	names(completeness_df) = c("OTU", "nPar", "nQ", "nQ_woInvar")
	completeness_df$nPar = as.numeric(completeness_df$nPar)
	completeness_df$nQ = as.numeric(completeness_df$nQ)
	completeness_df$nQ_woInvar = as.numeric(completeness_df$nQ_woInvar)
	completeness_df$numchars = numchars - (completeness_df$nPar + completeness_df$nQ)
	completeness_df$completePercent = round(100*(completeness_df$numchars/numchars), digits=2)

	completeness_df$numchars2 = numchars - (completeness_df$nPar + completeness_df$nQ_woInvar)
	completeness_df$completePercent2 = round(100*(completeness_df$numchars2/numchars), digits=2)
	completeness_df

	ttl_nPar = sum(completeness_df$nPar)
	ttl_nQ = sum(completeness_df$nQ)
	ttl_nQ_woInvar = sum(completeness_df$nQ_woInvar)

	ttl_nchar = numchars * numtaxa
	complete_chars = sum(completeness_df$numchars)
	complete_percent = round(complete_chars / ttl_nchar * 100, digits=2)
	complete_chars2 = sum(completeness_df$numchars2)
	complete_percent2 = round(complete_chars2 / ttl_nchar * 100, digits=2)

	matrix_stats = c(numtaxa, numchars, ttl_nchar, ttl_nPar, ttl_nQ, ttl_nQ_woInvar, complete_chars, complete_percent, complete_chars2, complete_percent2)
	matrix_stats_df = as.data.frame(as.matrix(matrix_stats, ncol=1), stringsAsFactors=FALSE)
	row.names(matrix_stats_df) = c("numtaxa", "numchars", "ttl_nchar", "ttl_nPar", "ttl_nQ", "ttl_nQ_woInvar", "complete_chars", "complete_percent", "complete_chars2", "complete_percent2")
	names(matrix_stats_df) = "matrix_stats"
	matrix_stats_df
	
	cat("\n\nStatistics on the completeness of the morphology data matrix, by taxon:\n\n")
	print(completeness_df)

	cat("\n\nStatistics on the completeness of the morphology data matrix, overall:\n\n")
	print(matrix_stats_df)
	
	morphstats = NULL
	morphstats$completeness_df = completeness_df
	morphstats$matrix_stats_df = matrix_stats_df
	
	return(morphstats)
	} # END morphology_matrix_stats <- function(charslist=NULL, charsdf=NULL)
