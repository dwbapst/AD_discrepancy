
# variables

# duration of the entire range
	# D&A treat this as a fixed variable based on avg species duration
	# for now let's say 1
T_R <- 1

# MINIMUM amount of time that an
	# ancestor's FAD (H_A) post-dates a descendant's FAD (H_D)
	# this is a fixed variable of interest in D&A
		# a lower bound of observed discrepancy
	# e.g. 0.1 
T_d <- 0.1

# amount of range overlap between anc and desc
	# this is a MINIMUM and represents range extention
		# in either direction (ancestor sampled later or
		# descendant originated earlier)
	#
	# the time between the true origin of the ancestor
		# and the first observed find of the descendant (H_D)
	#
	# D&A treat as identical to the time between the
		# true origin of the descendant and the 
		# first observed find of the ancestor (H_A)
		# see (Fig 2)
	#
	# this is a random variable in D&A
	# for now let's say 0.5
T_o <- 0.5

# for desc FAD (H_D) to post-date anc FAD (H_A)
# T_R > T_o > T_d

# Thus, T_R places a maximum limit on T_o values
#  and, T_d places a minimum limit on T_o values

# The difference between T_d and T_o is X_A and X_D

#####
# X_A and X_D

# for T_o to be greater than T_d, the ancestor 
	# must be sampled some time (H_A) *after*
	# the first sampled occurrence of descendant (H_D) + T_d
# this is X_A

# OR 

# for T_o to be greater than T_d, the descendant 
	# must be sampled some time (H_D) *before*
	# the first sampled occurrence of ancestor (H_A) - T_d
# this is X_A	

# X_A and X_D are thus random variables, with minimums at 0

# for overlap to occur at all, 
	# X_A > X_D
# And
	# X_A + X_D = T_o - T_d
# (sort of - remember that T_o is a random variable with a minimum)




#####################################################
# equation 1
# probability that you'd sample the descendant before the ancestor
	# given the total ranges and the amount of overlap

P_endD <- (T_o - T_d) / T_R

# and because we assume that the
	# ancestor and descendant have identical durations
	# we assume that the potential amount of time between
	# T_o and T_d for H_A and H_D is distributed similarly
P_endA <- P_endD

############################################
# equ 2
# joint probabilities for P_endD and P_endA

P_endA_endD <- ((T_o - T_d)^2) / (T_R^2)

####
# skipping ahead

#########

# equation 3g
# prob (X_A > XD) = lambda / 2*lambda = 1/2

# equation 5
# prob (H_A - H_D > T_D, given P_endA_endD) is 0.5

# equation 5b

# prob (H_A - H_D > T_D) = prob (H_A - H_D > T_D, given P_endA_endD) * P_endA_endD
# prob (H_A - H_D > T_D) = 1/2 * ((T_o - T_d)^2) / (T_R^2)

# let's use notation:
# P_H_diff_greater = prob (H_A - H_D > T_D) 

P_H_diff_greater <- ((T_o - T_d)^2) / (2*(T_R^2)) # 5b

# prob (H_A - H_D > T_D) must be zero when T_d > T_o
# thus 5b only holds when T_d > T_o

# D&A then used
T_R <- 0.97
T_d <- 0.8
# and tried many values of T_o
	# presumably from T_o -> T_R

# fig 3 actually shows lots of values with T_d > T_o, and Prob at 0
	# why even show that
# obviously the presumed overlap must be greater than the
	# age discrepancy between FADs of Anc & Desc
	
##########################################

# let's try to recreate their figure 3

# note they put both an absolute and relative axes for range
	# overlap on this plot
#Relative overlap = total overlap (T_o) / total taxon range (T_R)
	# but since T_R is so close to 1, they aren't visually distinct

# equation 5b
# P_H_diff_greater <- ((T_o - T_d)^2) / (2*(T_R^2)) 

T_R <- 0.97
T_d <- 0.8

T_o <- seq(T_d, T_R, 
	length.out = 1000)

P_H_diff_greater <- sapply(T_o,
	function(T_o_i) ((T_o_i - T_d)^2) / (2*(T_R^2)) 
	)
	
plot(T_o, P_H_diff_greater,
	xlab = "Range Overlap (Possible value of T_o)",
	ylab = "Prob of ((Desc's FAD - Anc's FAD) > T_d)",
	main = ""
	)

# okay, that was almost too easy to recreate...

###################################################

# what if we constrain T_o to 0.9 (Homo originated at most 100kyr before sediba)
# and asked how Prob of Discrepancy varied with T_R

T_d <- 0.8
T_o <- 0.97

# let's say that Homo and Sediba might be at most 3 million years old
max_T_R <- 2

T_R <- seq(T_o, max_T_R, 
	length.out = 1000)

P_H_diff_greater <- ((T_o - T_d)^2) / (2*(T_R^2))
	
plot(T_R, P_H_diff_greater,
	xlab = "Total Taxon Range (T_R)",
	ylab = "Prob of ((Desc's FAD - Anc's FAD) > T_d)",
	main = ""
	)

# very low probabilities??
	# worse if T_o is larger
	
##########################################################

library(paleotree)

# invert avg species duration to get ext rate
extRate <- 1/0.97
# set ext rate equal to orig rate
origRate <- extRate
# set sampling rate maybe equal to ext rate (??)
sampRate <- extRate

# shouldn't matter as sampling rate falls out of A&D's equations...

# set a seed
set.seed(44)
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
	nSamp = c(100,1000),
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
	AD_diff <- NA
}else{
	# whats the anc-desc FAD difference
	ancFADs <- ranges[sampAnc,1]
	descFADs <- ranges[names(sampAnc),1]
	AD_diff <- ancFADs - descFADs	
	}
	
hist(AD_diff)

##################################

library(paleotree)

# let's do a bunch and concatanate
nRun <- 100

# invert avg species duration to get ext rate
extRate <- 1/0.97
# set ext rate equal to orig rate
origRate <- extRate
# set sampling rate maybe equal to ext rate (??)
sampRate <- extRate
# shouldn't matter as sampling rate falls out of A&D's equations...

# set a seed
set.seed(44)
AD_diff_all <- list()
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
		AD_diff <- NA
	}else{
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
		}
	AD_diff_all[[run]] <- AD_diff
	}
	
AD_diff_all <- unlist(AD_diff_all)
AD_diff_all <- AD_diff_all[!is.na(AD_diff_all)]
hist(AD_diff_all)

T_d_obs <- AD_diff_all[AD_diff_all < 0]

AD_diff_sort <- sort(AD_diff_all)
cumProb <- cumsum(rep(1,length(AD_diff_sort)))/length(AD_diff_sort)

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

# but how to compare to D&A results? 
	# These simulations have variable T_R, T_o
	# and T_R and T_o differ between anc and desc
# basically noncomparable

# does relationship of T_o and P(T_d) hold up even with variable T_R?




###########################################################

# does being extant have any affect on these



