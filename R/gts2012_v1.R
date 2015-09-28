
get_GTS2012 <- function(level="all", quiet=FALSE)
	{
	defaults='
	level="Stage.Age"
	quiet=TRUE
	'
	cols = c("Precambrian", "Eon", "Era", "Period", "Series.Epoch")
	
	# Gradstein, F.M., Ogg, J.G., Schmitz, M.D. and Ogg, G.M. (2012). 
	# "The Geologic Time Scale 2012", Elsevier, 2 Volume set.
	# https://engineering.purdue.edu/Stratigraphy/charts/chart.html
	
	gts2012 = structure(list(Precambrian = structure(c(2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 
3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 1L), .Label = c("age_of_Earth", 
"Phanerozoic", "Precambrian"), class = "factor"), Eon = structure(c(4L, 
4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
4L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 2L, 2L, 2L, 
2L, 3L, 1L), .Label = c("age_of_Earth", "Archean", "Hadean", 
"Phanerozoic", "Proterozoic"), class = "factor"), Era = structure(c(2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 
7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 
7L, 7L, 7L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 
12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 
12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 
12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 8L, 
8L, 8L, 5L, 5L, 5L, 10L, 10L, 10L, 10L, 9L, 6L, 11L, 3L, 4L, 
1L), .Label = c("age_of_Earth", "Cenozoic", "Eoarchean", "Hadean", 
"Meso-proterozoic", "Mesoarchean", "Mesozoic", "Neo-proterozoic", 
"Neoarchean", "Paleo-proterozoic", "Paleoarchean", "Paleozoic"
), class = "factor"), Period = structure(c(18L, 18L, 18L, 18L, 
18L, 13L, 13L, 13L, 13L, 13L, 13L, 13L, 13L, 16L, 16L, 16L, 16L, 
16L, 16L, 16L, 16L, 16L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 
6L, 6L, 6L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 
12L, 25L, 25L, 25L, 25L, 25L, 25L, 25L, 17L, 17L, 17L, 17L, 17L, 
17L, 17L, 17L, 17L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 8L, 8L, 8L, 8L, 
8L, 8L, 8L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 14L, 14L, 
14L, 14L, 14L, 14L, 14L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
4L, 10L, 7L, 24L, 23L, 9L, 3L, 22L, 15L, 19L, 20L, 2L, 2L, 2L, 
2L, 11L, 1L), .Label = c("age_of_Earth", "blank", "Calymmian", 
"Cambrian", "Carboniferous", "Cretaceous", "Cryogenian", "Devonian", 
"Ectasian", "Ediacaran", "Hadean", "Jurassic", "Neogene", "Ordovician", 
"Orosirian", "Paleogene", "Permian", "Quaternary", "Rhyacian", 
"Siderian", "Silurian", "Statherian", "Stenian", "Tonian", "Triassic"
), class = "factor"), Series.Epoch = structure(c(7L, 27L, 27L, 
27L, 27L, 28L, 28L, 24L, 24L, 24L, 24L, 24L, 24L, 25L, 25L, 4L, 
4L, 4L, 4L, 26L, 26L, 26L, 33L, 33L, 33L, 33L, 33L, 33L, 10L, 
10L, 10L, 10L, 10L, 10L, 35L, 35L, 35L, 19L, 19L, 19L, 19L, 12L, 
12L, 12L, 12L, 39L, 39L, 39L, 23L, 23L, 16L, 16L, 9L, 9L, 6L, 
6L, 6L, 3L, 3L, 3L, 3L, 38L, 38L, 22L, 15L, 36L, 20L, 13L, 34L, 
34L, 18L, 18L, 11L, 11L, 11L, 29L, 17L, 17L, 40L, 40L, 8L, 8L, 
8L, 37L, 37L, 37L, 21L, 21L, 14L, 14L, 5L, 5L, 5L, 31L, 31L, 
31L, 30L, 30L, 32L, 32L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 2L), .Label = c("", "age_of_Earth", "Cisuralian", 
"Eocene", "Furongian", "Guadalupian", "Holocene", "Llandovery", 
"Lopingian", "Lower Cretaceous", "Lower Devonian", "Lower Jurassic", 
"Lower Mississippian", "Lower Ordovician", "Lower Pennsylvanian", 
"Lower Triassic", "Ludlow", "Middle Devonian", "Middle Jurassic", 
"Middle Mississippian", "Middle Ordovician", "Middle Pennsylvanian", 
"Middle Triassic", "Miocene", "Oligocene", "Paleocene", "Pleistocene", 
"Pliocene", "Pridoli", "Series 2", "Series 3", "Terreneuvian", 
"Upper Cretaceous", "Upper Devonian", "Upper Jurassic", "Upper Mississippian", 
"Upper Ordovician", "Upper Pennsylvanian", "Upper Triassic", 
"Wenlock"), class = "factor"), Stage.Age = structure(c(45L, 96L, 
48L, 18L, 37L, 66L, 102L, 60L, 92L, 80L, 55L, 17L, 8L, 25L, 74L, 
69L, 13L, 58L, 101L, 89L, 78L, 27L, 59L, 20L, 77L, 26L, 95L, 
23L, 5L, 7L, 12L, 42L, 97L, 16L, 90L, 52L, 64L, 19L, 15L, 11L, 
2L, 91L, 67L, 82L, 43L, 71L, 62L, 22L, 54L, 6L, 63L, 47L, 24L, 
100L, 21L, 99L, 73L, 53L, 9L, 75L, 10L, 41L, 50L, 61L, 14L, 79L, 
98L, 93L, 33L, 36L, 38L, 31L, 32L, 68L, 56L, 70L, 57L, 39L, 46L, 
81L, 88L, 3L, 72L, 44L, 51L, 76L, 29L, 28L, 34L, 94L, 83L, 49L, 
65L, 40L, 30L, 87L, 86L, 85L, 84L, 35L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 4L), .Label = c("", "Aalenian", 
"Aeronian", "age_of_Earth", "Albian", "Anisian", "Aptian", "Aquitanian", 
"Artinskian", "Asselian", "Bajocian", "Barremian", "Bartonian", 
"Bashkirian", "Bathonian", "Berriasian", "Burdigalian", "Calabrian", 
"Callovian", "Campanian", "Capitanian", "Carnian", "Cenomanian", 
"Changhsingian", "Chattian", "Coniacian", "Danian", "Dapingian", 
"Darriwillan", "Drumian", "Eifelian", "Emsian", "Famennian", 
"Floian", "Fortunian", "Frasnian", "Gelasian", "Givetian", "Gorstian", 
"Guzhangian", "Gzhelian", "Hauterivian", "Hettangian", "Hirnantian", 
"Holocene", "Homerian", "Induan", "Ionian", "Jiangshanian", "Kasimovian", 
"Katian", "Kimmeridgian", "Kungurian", "Ladinian", "Langhian", 
"Lochkovian", "Ludfordian", "Lutetian", "Maastrictian", "Messinian", 
"Moscovian", "Norian", "Olenekian", "Oxfordian", "Paibian", "Piacenzian", 
"Pliensbachian", "Pragian", "Priabonian", "Pridoli", "Rhaetian", 
"Rhuddanian", "Roadian", "Rupelian", "Sakmarian", "Sandbian", 
"Santonian", "Selandian", "Serpukhovian", "Serravallian", "Sheinwoodian", 
"Sinemurian", "Stage 10", "Stage 2", "Stage 3", "Stage 4", "Stage 5", 
"Telychian", "Thanetian", "Tithonian", "Toarcian", "Tortonian", 
"Toumaisian", "Tremadocian", "Turonian", "Upper Pleistocene", 
"Valanginian", "Visean", "Wordian", "Wuchiapingian", "Ypresian", 
"Zanclean"), class = "factor"), tops = c(0, 0.0118, 0.126, 0.781, 
1.806, 2.588, 3.6, 5.333, 7.246, 11.63, 13.82, 15.97, 20.44, 
23.03, 28.1, 33.9, 37.8, 41.2, 47.8, 56, 59.2, 61.6, 66, 72.1, 
83.6, 86.3, 89.8, 93.9, 100.5, 113, 123.6, 130.8, 133.9, 139.4, 
145, 152.1, 157.3, 163.5, 166.1, 168.3, 170.3, 174.1, 182.7, 
190.8, 199.3, 201.3, 209.5, 228.4, 237, 241.5, 247.1, 250, 252.2, 
254.2, 259.8, 265.1, 268.8, 272.3, 279.3, 290.1, 295.5, 298.9, 
303.7, 307, 315.2, 323.2, 330.9, 346.7, 358.9, 372.2, 382.7, 
387.7, 393.3, 407.6, 410.8, 419.2, 423, 425.6, 427.4, 430.5, 
433.4, 438.5, 440.8, 443.8, 445.2, 453, 458.4, 467.3, 470, 477.7, 
485.4, 489.5, 494, 497, 500.5, 504.5, 509, 514, 521, 529, 541, 
635, 850, 1000, 1200, 1400, 1600, 1800, 2050, 2300, 2500, 2800, 
3200, 3600, 4000, 4560), CI = c(NA, NA, NA, NA, NA, NA, NA, NA, 
NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5, 0.2, NA, NA, NA, 0.05, 
0.2, 0.2, 0.5, 0.3, 0.2, 0.4, 0.4, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
1, 1.1, 1.2, 1.3, 1.4, 1, 0.7, 1, 0.3, 0.2, NA, NA, 1, 1, 0.2, 
0.5, 0.5, 0.3, 0.4, 0.4, 0.5, 0.5, 0.6, 0.2, 0.4, 0.2, 0.1, 0.2, 
0.2, 0.4, 0.3, 0.4, 0.4, 1.6, 1.6, 0.8, 1.2, 2.6, 2.8, 3.2, 2.3, 
0.9, 0.5, 0.7, 0.8, 1.1, 1.2, 1.5, 1.4, 0.7, 0.9, 1.1, 1.4, 1.4, 
1.9, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1, NA, NA, NA, NA, NA, 
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), Tilde = structure(c(1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 2L), .Label = c("", "~"), class = "factor"), Notes = structure(c(2L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L), .Label = c("", "stage blank"), class = "factor")), .Names = c("Precambrian", 
"Eon", "Era", "Period", "Series.Epoch", "Stage.Age", "tops", 
"CI", "Tilde", "Notes"), class = "data.frame", row.names = c(NA, 
-116L))
	
	if (quiet != TRUE)
		{
		txt = paste0("\n\nThis gts2012 dataset was typed into Excel by Nick Matzke on 2014-12-07. \n\nOriginal source:\n\nGradstein, F.M., Ogg, J.G., Schmitz, M.D. and Ogg, G.M. (2012). 'The Geologic Time Scale 2012', Elsevier, 2 Volume set.\nhttps://engineering.purdue.edu/Stratigraphy/charts/chart.html\n\n")
		cat(txt)
		}

# 	if (level != "all")
# 		{
# 		cmdstr = paste0("return(gts2012$", level, ")")
# 		eval(parse(text=cmdstr))
# 		}
	
	if (level != "all")
		{
		# Evaluates to "factors"
		cmdstr = paste0("factors = unique(gts2012$", level, ")")
		eval(parse(text=cmdstr))
		
		gts2012_subset = NULL
		for (i in 1:length(factors))
			{
			f = factors[i]
			
			# Evaluates to TF
			cmdstr = paste0("TF = gts2012$", level, " == f")
			eval(parse(text=cmdstr))
			
			tmp = gts2012[TF,]
			toprowTF = tmp$tops == min(tmp$tops)
			#botrowTF = tmp$tops == max(tmp$tops)
			
			gts2012_subset = rbind(gts2012_subset, tmp[toprowTF, ])
			#gts2012_subset = rbind(gts2012_subset, tmp[botrowTF, ])
			} # END for (i in 1:length(factors))
		
		return(gts2012_subset)
		} # END if (level != "all")
	
	# Otherwise, return the full table
	return(gts2012)
	} # END get_GTS2012()

junk='	
get_GTS2012(quiet=TRUE)
get_GTS2012(level="Series.Epoch", quiet=TRUE)
get_GTS2012(level="Period", quiet=TRUE)
get_GTS2012(level="Era", quiet=TRUE)
'

