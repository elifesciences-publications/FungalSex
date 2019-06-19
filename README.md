# FungalSex
Code and data for the article entitled: 


"A direct ecological benefit of outcrossing in an obligate pathogen" 

### This code is associated with the paper from Laine et al., "Variable opportunities for outcrossing result in hotspots of novel genetic variation in a pathogen metapopulation". eLife, 2019. http://dx.doi.org/10.7554/eLife.47091



In particular: 

AllDataFungalSex.csv   contains the data used for the statistical models

run_final_models.R     has the code to run the models in the paper (also other models presented)

fungalsexfunctions.R   contains some helper functions to run the models or to generate the figures. 




Variables in the data: 


patch 				the id of the local host population
year				the year of sampling
numberofstrains			the number of strains observed in the local population at the year of sampling
numberofstrainsnextyear		the number of strains observed in total in the same location the next year
newstrainsnextyear		the number of strains observed the next year that are new to the location 	
pathogenconnectivity		the computed pathogen connectivity (i.e. the expected infection load from outside the local population)
hostconnectivity		the computed host connectivity based on the host coverages in other host locations
nsamples			number of samples collected from the local population at the year of sampling
proportionofcoinfected		the proportion of the collected samples that were found out to be coinfected 
numberofcoinfections     	the absolute number of coinfections in the population
absoluteabundance		The categorical variable describing the abundance of the powdery mildew in the local population
x				The longitude of the patch from where the samples were collected from
y				The latitude of the patch from where the samples were collected from
pl				The log-host covarage in the patch (m²)
area				The geographical area of the patch, within which the host plants are distributed
coinfpresence			A binary variable describing whether there was a coinfection observed among the samples
abundance1			The categorical abundance variable transformed into a three binary variables. 
abundance2
abundance3
mildewpresencenextyear		A binary variable descrbing mildew presence next year. 
Intercept			A vector of ones (used in the modelling)
