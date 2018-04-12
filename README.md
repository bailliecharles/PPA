# PPA
Tools for doing posterior predictive analyses of Bayesian phylogenetics in R

I'll add some graphics scripts to this repo eventually but for now this single script will do three different types of PPA. There is an example P4 script to do the simulations, and if you want to write them to file without too much clutter I've also provided a tweaked version of the alignment_readwrite.py script that comes with P4 in which I've added append. Just swap the version installed for this one. 

As input it needs:
1. the empirical alignment (all ambiguities, like N's, or uncertain characters, ?, need to coded as "-"). It also needs to know whether 3rd codon positions are RY-coded with the -ry option. Assumes empirical is actually 'R' and 'Y', but the sims are {a,g,c,t} albeit in the correct proportions of purines to pyrimidines. 
2. the format of that alignment (only fasta or phylip allowed)
3. a file containing alignments generated from the posterior parameter estimates of a Bayesian phylogenetic analysis. This file is a little awkward: all the alignments need to be in one file with the same number of lines between each alignment (every line that isn't a sequence containing line, i.e. phylip headers). This number of lines needs to be specified and with the -l option.   

I've been using P4 to do these simluations from MrBayes runs, mainly because it can handle simulation on partitioned models. Maybe there are other better things out there? Eventually I'll also write a script to parse MrBayes NEXUS file into a P4 script. The kind of neat thing is that the script will also read in the simulations produced by PhyloBayes. Just run readpb_mpi -ppred on one of your chains to simulate a bunch of post burnin alignments, then cat \*.ali > pb_sims.ali. The -l option for this file will be 2 (the phylip header, and one blank space). PhyloBayes does also perform these tests but I have no idea how they are calculated and I prefer to use the same metrics across different methods.  

The PPAs performed are all 'data' checks rather than 'inference' based:
1. PPADIV - the mean nucleotide/aminoacid diversity per site
2. PPAX2 - a Chi-squared test for homogeneity of base frequencies between taxa
3. PPAMULTI - the multinomial likelihood test of site patterns from Bollback (2002)

With default settings it will calculate all three stats and output for each, in one file:
1. The empirical statistic
2. The range of the simulated statistics
3. The z-score (calculated as the median of the simulated values subtracted from the empirical value, divided by the standard deviation of the simulated values)
There is also an option to print out all values of each test for each simulated alignment. 

This is not a very fast script - there are a lot of \*apply functions and data manipulations going on! For 1,000 simulated alignments (which is probably more than representative) of around 10K sites, it will run in approximately 3 minutes. 

