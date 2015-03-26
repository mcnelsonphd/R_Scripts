########################### SCRIPT NOTES #################################
##                                                                      ##
## R-script to take a QIIME v1.8 generated OTU biom file that has       ##
## meta-data attached to it and perform NMDS ordination and then pretty ##
## plot creation using ggplot. Ellipses, sample colorings and plotting  ##
## characters are determined by one or more meta-data fields in the     ##
## biom file.                                                           ##
##                                                                      ##
## NB: Currently, this is not an as-is executable script. The user will ##
## need to change some of the hardcoded variables for proper execution. ##
##                                                                      ##
## Items to change prior to running are:                                ##
## PATH_TO_BIOM_FILE to the absolute path of the input biom file.       ##
## The meta-data fields and names to use for ellipse generation and     ##
## plotting aesthetics, currently set as METANAME1, METANAME2, META1,   ##
## META2.
#########################################################################

require(vegan)
require(ggplot2)
require(biom)
require(grid) # Required if putting on env. fit vectors to get arrowheads a-la https://oliviarata.wordpress.com/2014/04/17/ordinations-in-ggplot2/
# Env fitting not currently done.
# Read in the biom file, which should have meta-data attached to it.
mybiom <- read_biom("PATH_TO_BIOM_FILE")
# Get some basic info on the file
mybiom
# Put the species/sites data into a matrix and then transpose it to put into the correct "orientation" and make it a data frame
tmp <- as(biom_data(mybiom), "matrix")
mybiom.data <- as.data.frame(t(tmp))
rm(tmp)
# Create a sample metadata matrix for easy referencing
mybiom.meta <- sample_metadata(mybiom)
# Print out the names of the meta-data fields
colnames(mybiom.meta)
# Now lets ask the user which meta-data field they want to use for grouping
#while(n<1){  This assumes script is being run interactively, probably want to do this programmatically instead.
#  meta <- readline("Please type the name of the meta-data you wish to use for sample grouping: ")
#  meta <- ifelse(grepl("\\D",n),-1,as.integer(n))
#  if(is.na(n)){break}  # breaks when hit enter
#}

# Now for the fun, run NMDS on the observation matrix
mybiom.nmds <- metaMDS(mybiom.data, k=2, trymax=500)
mybiom.nmds
stressplot(mybiom.nmds)

# Create a quick and dirty plot for initial analysis, also prevents ordiellipse from plotting onto the stressplot.
ordiplot(mybiom.nmds)

# Create new data frames with the values for the two NMDS vectors and ellipse values
# Note that METANAME1=mybiom.meta$META1 referes to a meta-data column of the metadata data frame, and that METANAME will be used in the legend to identify the meta-data field.
# Same for METANAME2
mybiom.nmds.sites <- data.frame(NMDS1=mybiom.nmds$points[,1],NMDS2=mybiom.nmds$points[,2],METANAME1=mybiom.meta$META1,METANAME2=mybiom.meta$META2)
mybiom.nmds.species <- data.frame(NMDS1=mybiom.nmds$species[,1],NMDS2=mybiom.nmds$species[,2])
mybiom.ellipse <- ordiellipse(mybiom.nmds,mybiom.meta$META1,display="sites",kind="sd",conf=0.95,label=T)

# Create a function to get the data used for creating the ellipses
# Taken from http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
veganCovEllipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

ellipse.data <- data.frame()
for(g in levels(mybiom.nmds.sites$META1)){
  ellipse.data <- rbind(ellipse.data, cbind(as.data.frame(with(mybiom.nmds.sites[mybiom.nmds.sites$META1==g,],
                  veganCovEllipse(mybiom.ellipse[[g]]$cov,mybiom.ellipse[[g]]$center,mybiom.ellipse[[g]]$scale)))
                  ,METANAME1=g))
}
rm(g)
# Now lets construct the plot by layer
mybiom.plot <- ggplot(data=mybiom.nmds.sites, aes(NMDS1, NMDS2)) +
               theme_minimal() +
               geom_point(aes(color=METANAME2, shape=METANAME1), size=5) + # Plots the sample points and colors according to temp, point according to feedstock.
               geom_path(data=ellipse.data,aes(x=NMDS1, y=NMDS2,linetype=META1), size=0.75) + # Plots the ellipses, using different line types for the different feedstocks.
               geom_point(data=mybiom.nmds.species, color='grey10', size=0.75) # Plots the species as small grey dots to give more context.
# And now display it
mybiom.plot 

