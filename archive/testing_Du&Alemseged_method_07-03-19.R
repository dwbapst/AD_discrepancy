library(paleotree)

setwd("~/AD_discrepancy")

# let's do a bunch and concatanate
nRun <- 10000

# invert avg species duration to get ext rate
extRate <- 1/0.97
# set ext rate equal to orig rate
origRate <- extRate
# set sampling rate maybe equal to ext rate (??)
sampRate <- extRate
# shouldn't matter as sampling rate falls out of A&D's equations...

# set a seed
set.seed(44)
isExtant_all <- AD_overlap_all <- AD_diff_all <- list()
# iterate over nRun
for(run in 1:nRun){
	# get a simulated fossil record
	record <- simFossilRecord(
		p = origRate,
		q = extRate, 
		r = sampRate, 
		nruns = 1,
		# controls
		totalTime = c(0, 1000), 
		nTotalTaxa = c(1, 1000),
		# basically only conditioned on at least 1 extant taxon
		nExtant = c(1, 1000),
		# and 10 to 100 sampled fossil taxa
		nSamp = c(10,100),
		plot = TRUE
		)
	# convert to taxa data
	taxa <- fossilRecord2fossilTaxa(record)
	# get ranges
	ranges <- fossilRecord2fossilRanges(record)
	# drop NAs
	ranges <- ranges[!is.na(ranges[,1]),]
	# get sampled taxa
	sampledTaxaNames <- rownames(ranges)
	# get last sampled ancestor for each taxon
	sampAnc <- character()
	for(i in 1:length(sampledTaxaNames)){
		taxon <- sampledTaxaNames[i]
		isSampled <- FALSE
		while(!isSampled){
			ancID <- record[[taxon]]$taxa.data["ancestor.id"]
			anc <- names(record)[ancID]
			if(is.na(anc)){
				isSampled <- TRUE
			}else{
				isSampled <- any(anc == sampledTaxaNames)
				taxon <- anc
				}
			}
		sampAnc[i] <- anc
		}
	# filter to those with sampled ancestors
	names(sampAnc) <- sampledTaxaNames
	sampAnc <- sampAnc[!is.na(sampAnc)]
	# are there any with sampled ancestors
	if(length(sampAnc)<1){
		isExtant <- AD_overlap <- AD_diff <- NA
	}else{
		# need to filter so each sampled ancestor
			# is only used once for independence concerns
		resample <- function(x, ...) x[sample.int(length(x), ...)]
		uniqDesc <- sapply(unique(sampAnc),function(x) 
			names(sampAnc)[
				resample(which(sampAnc == x),1)
				]
			)
		sampAnc <- unique(sampAnc)
		names(sampAnc) <- uniqDesc
		#
		# whats the anc-desc FAD difference
		ancFADs <- ranges[sampAnc,1]
		descFADs <- ranges[names(sampAnc),1]
		AD_diff <- ancFADs - descFADs	
		#
		# get the true overlap between taxon ranges
			# start - desc orig time
			# end - and or desc ext, whichever comes first
		desc_orig <- sapply(names(sampAnc),function(taxon)
			record[[taxon]]$taxa.data["orig.time"]
			)
		desc_ext <- sapply(names(sampAnc),function(taxon)
			record[[taxon]]$taxa.data["ext.time"]
			)
		anc_ext <- sapply(sampAnc,function(taxon)
			record[[taxon]]$taxa.data["ext.time"]
			)
		end_date <- max(anc_ext,desc_ext)
		AD_overlap <- desc_orig - end_date
		# 
		# are either the anc or desc extant
		desc_extant <- sapply(names(sampAnc),function(taxon)
			record[[taxon]]$taxa.data["still.alive"] == 1
			)
		anc_extant <- sapply(sampAnc,function(taxon)
			record[[taxon]]$taxa.data["still.alive"] == 1
			)
		isExtant <- desc_extant | anc_extant
		}
	AD_diff_all[[run]] <- AD_diff
	AD_overlap_all[[run]] <- AD_overlap
	isExtant_all[[run]] <- isExtant
	}

save.image("saved_AD_discrepancy_workspace.Rda")

}
	
############################################

# discrepancy
AD_diff_all <- unlist(AD_diff_all)
AD_diff_all <- AD_diff_all[!is.na(AD_diff_all)]
#hist(AD_diff_all)
#summary(AD_diff_all)

# overlap
AD_overlap_all <- unlist(AD_overlap_all)
AD_overlap_all <- AD_overlap_all[!is.na(AD_diff_all)]
#hist(AD_overlap_all)
#summary(AD_overlap_all)

# isExtant
isExtant_all <- unlist(isExtant_all)
isExtant_all <- isExtant_all[!is.na(AD_diff_all)]
#sum(isExtant_all)
#sum(isExtant_all)/length(isExtant_all)

###############################################

# but how to compare to D&A results? 
	# These simulations have variable T_R, T_o
	# and T_R and T_o differ between anc and desc
# basically noncomparable

# relationship of true overlap to observed discrepancy

plot(AD_diff_all, AD_overlap_all)
# negative values are due to taxa having no observed overlap

# discrepancy only
plot(AD_diff_all[AD_diff_all <= -1], 
	AD_overlap_all[AD_diff_all <= -1])
# very small negative overlaps are probably due to having extant taxa

# does relationship of T_o and P(T_d) hold up even with variable T_R?
	# would need to pick a particular discrepancy, then look at
	# cumulative freq of T_o in that set (I think?)

T_d <- 0

########################

# condition overlap on having a discrepancy greater than threshold
cond_overlap <- AD_overlap_all[AD_diff_all <= -(T_d+0.001)]
# sort 
cond_overlap <- sort(cond_overlap)
# calculate how many overlaps in total dataset are as big
	# as each of these
nOverlap <- sapply(cond_overlap, 
	function(x) sum(AD_overlap_all < x)
	)
# how many cond_overlap are less than each cond_overlap
nCondOverlap <- sapply(cond_overlap, 
	function(x) sum(cond_overlap <= x)
	)
# 
# cum prob that 
	# an anc-desc pair with this overlap would have a discrepancy
	# greater than T_d, relative to total number of
	# overlaps observed greater than that value
Prob_overlap <- nCondOverlap/nOverlap
#plot(nCondOverlap,nOverlap)
#
plot(-cond_overlap, Prob_overlap,
	xlab = "True Overlap between Anc and Desc Greater Than...",
	ylab = paste0(
		"Proportion Having a Discrepancy >= ",
		T_d, " Ma")
	)

# number of discrepancies greater than set T_d
length(cond_overlap)

#######################################################
# does being extant have any affect on these

## FULLY EXTINCT PAIR ONLY

AD_overlap_dead <- AD_overlap_all[!isExtant]
AD_diff_dead <- AD_diff_all[!isExtant]

# condition overlap on having a discrepancy greater than threshold
cond_overlap <- AD_overlap_dead[AD_diff_dead <= -(T_d+0.001)]
# sort 
cond_overlap <- sort(cond_overlap)
# calculate how many overlaps in total dataset are as big
	# as each of these
nOverlap <- sapply(cond_overlap, 
	function(x) sum(AD_overlap_dead < x)
	)
# how many cond_overlap are less than each cond_overlap
nCondOverlap <- sapply(cond_overlap, 
	function(x) sum(cond_overlap <= x)
	)
# 
# cum prob that 
	# an anc-desc pair with this overlap would have a discrepancy
	# greater than T_d, relative to total number of
	# overlaps observed greater than that value
Prob_overlap <- nCondOverlap/nOverlap
#plot(nCondOverlap,nOverlap)
#
plot(-cond_overlap, Prob_overlap,
	xlab = "True Overlap between Anc and Desc Greater Than...",
	ylab = paste0(
		"Proportion Having a Discrepancy >= ",
		T_d, " Ma")
	)

# number of discrepancies greater than set T_d
length(cond_overlap)

#####################################

## PAIRS WHERE AT LEAST ONE IS EXTANT

AD_overlap_live <- AD_overlap_all[isExtant]
AD_diff_live <- AD_diff_all[isExtant]

# condition overlap on having a discrepancy greater than threshold
cond_overlap <- AD_overlap_live[AD_diff_live <= -(T_d+0.001)]
# sort 
cond_overlap <- sort(cond_overlap)
# calculate how many overlaps in total dataset are as big
	# as each of these
nOverlap <- sapply(cond_overlap, 
	function(x) sum(AD_overlap_live < x)
	)
# how many cond_overlap are less than each cond_overlap
nCondOverlap <- sapply(cond_overlap, 
	function(x) sum(cond_overlap <= x)
	)
# 
# cum prob that 
	# an anc-desc pair with this overlap would have a discrepancy
	# greater than T_d, relative to total number of
	# overlaps observed greater than that value
Prob_overlap <- nCondOverlap/nOverlap
#plot(nCondOverlap, nOverlap)
#
plot(-cond_overlap, Prob_overlap,
	xlab = "True Overlap between Anc and Desc Greater Than...",
	ylab = paste0(
		"Proportion Having a Discrepancy >= ",
		T_d, " Ma")
	)

# number of discrepancies greater than set T_d
length(cond_overlap)





###########################################################
















T_d_obs <- AD_diff_all[AD_diff_all < 0]

AD_diff_sort <- sort(AD_diff_all)
cumProb <- cumsum(rep(1,length(AD_diff_sort))) / length(AD_diff_sort)

# reverse cumulative prob distribution
plot(AD_diff_sort, rev(cumProb))

# cumulative probability distribution
plot(AD_diff_sort, cumProb,
	xlab = "Ancestral First Appearance Time - Descendant First Appearance Date",
	ylab = "Empirical Cumulative Probability for Difference in FADs",
	main = "")
	
# the proportion of A-D pairs with Desc FAD before Anc Fad 
	# should be equal to the P(T_D < 0)
	# so let's look at that
# and look at discrepancies greater than 1 Ma

P_D_0 <- sum(AD_diff_all < 0) / length(AD_diff_all)
P_D_0

P_D_1 <- sum(AD_diff_all <= -1) / length(AD_diff_all)
P_D_1

P_D_2 <- sum(AD_diff_all <= -2) / length(AD_diff_all)
P_D_2

