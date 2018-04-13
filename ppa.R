#!/usr/bin/env Rscript


# load required packages
suppressMessages(suppressWarnings(library(seqinr, quietly=T)))
suppressMessages(suppressWarnings(library(ape, quietly=T)))
suppressMessages(suppressWarnings(library(optparse, quietly=T)))
suppressMessages(suppressWarnings(library(dplyr, quietly=T)))

###-----------------------------------------
###               Options 
###-----------------------------------------

option_list = list(
  make_option(c("-s", "--simulations"), type="character", default=NULL,
              help="posterior predictive simulations from P4 script", metavar="path"),
  make_option(c("-e", "--empirical_alignment"), type="character", default=NULL,
              help="empirical alignment - should be same dims as simulations", metavar="path"),
  make_option(c("-f", "--empirical_alignment_format"), type="character", default=NULL,
              help="format of empirical alignment: can be one of fasta or phylip", metavar="character"),
  make_option(c("-d", "--distribution"), type="character", default=FALSE,
              help="print all values of the posterior statitic (i.e. for each alignment)", metavar="logical"),
  make_option(c("-l", "--line_number"), type="numeric", default=1,
              help="number of lines in file between simulated alignments, including things like phylip headers", metavar="number"),
  make_option(c("-ry", "--RY_coding"), type="character", default=FALSE,
              help="is the third codon position coded as puRines and pYrimidines? assumes the empirical alignment has every 3rd site coded 
              as R and Y, but the simulations are coded a,g,c,t (as they are from P4).", metavar="logical"),
  make_option(c("-D", "--PPADIV"), type="character", default=TRUE,
              help="PPADIV: mean diversity per site across alignment", metavar="logical"),
  make_option(c("-C", "--PPAX2"), type="character", default= TRUE,
              help="PPAX2: Chi-squared test of compositional homogeneity between taxa", metavar="logical"),
  make_option(c("-M", "--PPAMULTI"), type="character", default= TRUE,
              help="PPAMULTI: Multinomial test of alignment site patterns sensu Bollback (2002) eqn. 7", metavar="logical")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$simulations) ){
  print_help(opt_parser)
  stop("Can't do PPA without simulations!", call.=FALSE)
}

if (is.null(opt$empirical_alignment) ){
  print_help(opt_parser)
  stop("Can't do PPA without an alignment against which to compare simulations!", call.=FALSE)
}

if (is.null(opt$empirical_alignment_format) ){
  print_help(opt_parser)
  stop("Please provide format of empirical alignment", call.=FALSE)
}


###-------------------------------
###      Read empirical data
###-------------------------------

empirical_alignment = read.alignment(opt$empirical_alignment,format=opt$empirical_alignment_format) # read in alignment
empirical_matrix = as.matrix.alignment(empirical_alignment) # convert to matrix
empirical_matrix[empirical_matrix=="x"]<-"-"   # just in case there happen to be a run of X in an AA alignment

if(opt$RY_coding==TRUE){                                 # See below re sims for RYcoding P4 sims. Same deal here, R and Y would be
                                                        # counted as 5th and 6th 'states' without changing them. 
  RYcode3rd = function(x){
    temp = x[,seq(3, ncol(x), 3)]
    temp2 = recode(temp,r="a", y="t", n="-")
    x[,seq(3, ncol(x), 3)] = temp2
    print(x)
  }
  empirical_matrix = RYcode3rd(empirical_matrix)
}


if (opt$PPADIV==TRUE){div_empirical = empirical_matrix} 
if (opt$PPAX2==TRUE){chix_empirical = empirical_matrix}
if (opt$PPAMULTI==TRUE){multi_empirical = empirical_matrix}

gap_idx = which(empirical_matrix=="-", arr.ind=T) # find out where the gaps are
gaps_cols = as.numeric(levels(as.factor(gap_idx[,2]))) # gap containing columns


###-----------------------------------------------
###         Read and reformat simulations
###-----------------------------------------------

# Here we're reading and then converting the flat text file (containing phylip formatted simulated alignments, each on a new line) 
# into a matrix the same as empirical_matrix. Couple of seconds for ~1000 simulations

P4_sim_sorter = function(x){
sims = scan(x, what="", sep="\n", quiet=T)
sims2 = split(sims, ceiling(seq_along(sims)/(nrow(empirical_matrix)+opt$line_number))) # split into separate alignments 
sims3 = lapply(sims2, function(x) {x[-opt$line_number]}) # drop the phylip header
sims4 = lapply(sims3, function(x) gsub(" ", "", x, fixed = TRUE)) #remove whitespace between name and sequence
sims5 = lapply(sims4, function(x) unlist(strsplit(x, "(?<=[[:upper:]])(?=[[:lower:]])", perl=TRUE))) # split name and sequence and unlist
sims6 = lapply(sims5, function(x) matrix(x,ncol=2, byrow=T)) # make into matrices
sims7 = lapply(sims6, function(x) t(sapply(strsplit(x[,2],""), tolower))) # split sites
sims8 = lapply(sims7, function(x) {row.names(x) <- as.character(sims6[[1]][,1]);x})} # name rows

sims = P4_sim_sorter(opt$simulations)

if(opt$RY_coding==TRUE){                                   #### P4 simulations have the correct proportion of purines and pyrimidines
                                                          #### at 3rd sites but coded as a,g,c,t. This will swap one of each for the other
  RYcode3rd = function(x){                                #### thus turning the 3rd sites binary again. 
    temp = x[,seq(3, ncol(x), 3)]
    temp2 = recode(temp,g="a", c="t", n="-")
    x[,seq(3, ncol(x), 3)] = temp2
    print(x)
  }
  sims = lapply(sims, function(x){RYcode3rd(x)}) 
} 
                           
               
if (opt$PPADIV== TRUE){div_sims = sims} 
if (opt$PPAX2==TRUE){chix_sims = sims}
if (opt$PPAMULTI==TRUE){multi_sims = sims}

###------------------------
###        PPADIV
###------------------------

if (opt$PPADIV==TRUE){
# 1. First we need to 'fill in' gaps in the empirical and simulated matrices. 
# This works because PPADIV is the diversity per column, so changing the gaps to another base already
# occuring in that column doesn't have any effect on the calculation, and its bettter than deleting all
# the columns... I think!

# 2. Data sorting
uncomplete_rows =  as.numeric(as.character(levels(as.factor(gap_idx[,1])))) # gap containing rows  
fill_in_gaps_row =  c(1:nrow(div_empirical))[-uncomplete_rows] # find all rows/seqs in empirical that have no gaps 

div_empirical[gap_idx[,1:2]] <- div_empirical[fill_in_gaps_row[1], c(gap_idx[,2])] # change the gaps to the same residue as in the first ungapped sequence

adjusted_div_empirical = div_empirical
adjusted_div_sims = lapply(div_sims, function(x) { x[gap_idx[,1:2]] <- x[fill_in_gaps_row[1], c(gap_idx[,2])] ; x}   ) # adjust the same sites in the sims

# 3. Define function 
# straight up no gaps version
ppaDiv <- function(x){
  mean(
    apply(x, 2, function(i)length(unique(i)))
  )
}

# 4. Run it on empirical and sims 
empiricalDiv = ppaDiv(adjusted_div_empirical)
simsDiv = as.numeric(unlist(lapply(adjusted_div_sims,ppaDiv)))

# 5. z-score
ZDiv = abs(empiricalDiv-median(simsDiv))/sd(simsDiv)

# 6. write to file
PPAs <- file("ppaScores.txt","w")
writeLines(".........PPADIV...........", PPAs)
writeLines(paste("The emprical value is: ", empiricalDiv, sep=""), PPAs)
writeLines(paste("The Z score is: ", ZDiv, sep=""), PPAs)
writeLines(paste("The range of the sims is: ", round(min(simsDiv),2), " to ", round(max(simsDiv),2), sep=""), PPAs)
close(PPAs)

if (opt$distribution==TRUE){
  PPA_dist <- file("PPADIV-distribution.txt","w")
  writeLines(".....PPADIV-distribution.....")
  writeLines(paste(simsDiv,collapse=","), PPA_dist)
  close(PPA_dist)
}
} # close off PPADIV

###-------------------------------
###         PPA Chi-Squared
###-------------------------------

if (opt$PPAX2==TRUE){
# Here gap containing sites contribute nothing to the final stat so the same gap cotaining index in the sims as in the empirical need to be made NA  
# In some papers the sites will be removed from the empirical dataset before simulation, however, I'm not sure how this would work with 
# a partitioned data set. The point of the PPAs here is to compare partitioned models with CAT models so I'm doing the sims, then 
# removing the data. Otherwise, we would have to redefine the partitions, which might (at worst) exclude whole partitions that are 
# actually contributing to the estimate of phylogeny, or at least change the partitions so much that a more appropraiate substitution model 
# may be required.Hmmmmmm not sure about this but at least the empirical and simulated datasets are treated the same, and the sims are 
# done for a partitioned model. Moving on:


chix_empirical[chix_empirical=="-"] <- NA # make empirical NA
adjusted_chix_empirical<- chix_empirical
adjusted_chix_empirical_table = table(row(adjusted_chix_empirical), as.matrix(adjusted_chix_empirical)) # make contingency table of bases


adjusted_chix_sims = lapply(chix_sims, function(x) { x[gap_idx[,1:2]] <- NA ; x}   ) # make sims NA in same place as empirical NA
adjusted_chix_sims_table = lapply(adjusted_chix_sims, function (x) table(row(x), as.matrix(x))) # calc contingency table for each sim- takes ~15s for 1000! 


empiricalX2 = chisq.test(adjusted_chix_empirical_table, simulate.p.value = F)$statistic
simsX2 = lapply(adjusted_chix_sims_table, function(x) chisq.test(x, simulate.p.value = F)$statistic)
simsX2 = as.numeric(unlist(simsX2))
ZChi = abs(empiricalX2-median(simsX2))/sd(simsX2)


PPAs <- file("ppaScores.txt","a")
writeLines(".........PPAX2...........",PPAs)
writeLines(paste("The emprical value is: ", empiricalX2, sep=""), PPAs)
writeLines(paste("The Z score is: ", ZChi, sep=""), PPAs)
writeLines(paste("The range of the sims is: ", round(min(simsX2),2), " to ", round(max(simsX2),2), sep=""), PPAs)
close(PPAs)

if (opt$distribution==TRUE){
  PPAX2_dist <- file("PPAX2-distributions.txt","w")
  writeLines(".....PPAX2-distribution.....")
  writeLines(paste(simsX2,collapse=","), PPAX2_dist)
  close(PPAX2_dist)
}
} # close off PPAX2

###--------------------------------
###        PPA-Multinomial
###--------------------------------

if (opt$PPAMULTI==TRUE){
# Multinomial proposed by Goldman 1993, follows eqn.7 in Bollback 2002. The idea is to measure site patterns across the 
# alignment. From Bolback, if we have an alignment 10 bases long, where three different site patterns each occur twice, and the 
# remaining four sites each occur once, then:
#T(x) = (3*2*(log(2)) + 4*(log(1))) - 10*(log(10))
#T(x) = -18.87
# To do this calculation we need to remove all the gap containing columns in the empirical, and corresponding sites in the sims.

adjusted_multi_empirical = multi_empirical[,-gaps_cols] # remove the gap containing columns

site_patterns = apply(adjusted_multi_empirical, 2, paste, collapse = "") #collapse columns to strings
unique_site_patterns = as.vector(table(site_patterns)) # get the number of occurences of each patterns
site_pattern_occurences = sort(unique(unique_site_patterns)) # and 
site_pattern_numbers = as.vector(table(unique_site_patterns))
empirical_multi_DF = as.data.frame(cbind(site_pattern_occurences, site_pattern_numbers))

sims2m = lapply(multi_sims, function(x) x[,-gaps_cols])
sims3m = lapply(sims2m, function(x) apply(x, 2, paste, collapse = ""))
sims4m = lapply(sims3m, function(x) as.vector(table(x)))
sims5m = lapply(sims4m, function(x) sort(unique(x)))
sims6m = lapply(sims5m, function(x) as.vector(table(x)))
sims_multi_DF = list()
for (i in 1:length(multi_sims)){
 sims_multi_DF[[i]]= as.data.frame(cbind(sims5m[[i]], sims6m[[i]]))
}

ppaMultinomial = function(x) {sum(x[,2]*x[,1]*(log(x[,1]))) - ncol(empirical_matrix)*(log(ncol(empirical_matrix)))} # define function

empiricalMulti = ppaMultinomial(empirical_multi_DF)
simsMulti = as.numeric(unlist(lapply(sims_multi_DF,ppaMultinomial)))


ZMulti = abs(empiricalMulti-median(simsMulti))/sd(simsMulti)


PPAs <- file("ppaScores.txt","a")
writeLines(".........PPAMULTI...........", PPAs)
writeLines(paste("The emprical value is: ", empiricalMulti, sep=""), PPAs)
writeLines(paste("The Z score is: ", ZMulti, sep=""), PPAs)
writeLines(paste("The range of the sims is: ", round(min(simsMulti),2), " to ", round(max(simsMulti),2), sep=""), PPAs)
close(PPAs)

if (opt$distribution==TRUE){
  PPAMULTI_dist <- file("PPAMUTLI-distributions.txt","w")
  writeLines(".....PPAMUTLI-distribution.....")
  writeLines(paste(simsMulti,collapse=","), PPAMUTLI_dist)
  close(PPAMUTLI_dist)
}                      
} # close off PPAMULTI

message("PPAs done, have a nice day")
