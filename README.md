# PPA
Tools for doing posterior predictive analyses of Bayesian phylogenetics in R

I'll add some graphics scripts to this repo eventually but for now this single script will do three different types of PPA. There is an example P4 script for do the simulations, and if you want to write them to file without too much clutter I've also provided a tweaked version of the alignment_readwrite.py script that comes with P4 in which I've added append. Just swap the version installed for this one. 

As input it needs the empirical alignment, the format of that alignment (only fasta or phylip allowed), and a file containing alignments generated from the posterior parameter estimates of Bayesian phylogenetics analysis. This last file is a little awkward: all the alignment need to be one file with the same number of lines between each alignment (every line that isn't a sequence containing line, i.e. phylip headers). This number of lines needs to be specified and with the -l option.   I've been using P4 to do these simluations from MrBayes runs, mainly because it can handle simulation on partitioned models. Maybe there are other better things out there? Eventually I'll also write a script to parse MrBayes NEXUS file into a P4 script. The kind of neat thing is that with a tweak the script will also read in the simulations produced by PhyloBayes (readpb_mpi -ppred). 

The PPAs performed are all 'data' checks rather than 'inference' based:
1. PPADIV - the mean nucleotide/aminoacid diversity per site
2. PPAX2 - a Chi-squared test for homogeneity of base frequencies between taxa
3. PPAMULTI - the multinomial likelihood test of site patterns from Bollback (2002)

With default settings it will calculate all three stats and output for each, in one file:
1. The empirical statistic
2. The range of the simulated statistics
3. The z-score (calculated as the median of the simulated values subtracted from the empirical value, divided by the standard deviation of the simulated values)
There is also an option to print out all values of each test for each simulated alignment. 

This is not a very fast script - there are a lot of \*apply functions and data manipulations going on! For 1,000 simulated alignments (which is probably more than representative) the script will run in approximately 3 minutes. 


Please let me know if there are any issues and it would be great if someone could double check the calculations... 
