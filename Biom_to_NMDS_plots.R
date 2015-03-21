require(vegan)
require(ggplot2)
require(biom)
require(grid) # Required if putting on env. vectors to get arrowheads a-la https://oliviarata.wordpress.com/2014/04/17/ordinations-in-ggplot2/
# Read in the biom file, which should have meta-data attached to it.
mybiom <- read_biom("~/Desktop/Sandbox/Qiime_Dev/Jacqui_table.biom")
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
# Note that the Feedstock=mybiom.meta$Feedstock referes to the Feedstock column of the metadata data frame
mybiom.nmds.sites <- data.frame(NMDS1=mybiom.nmds$points[,1],NMDS2=mybiom.nmds$points[,2],Feedstock=mybiom.meta$Feedstock,Temp=mybiom.meta$Temperature)
mybiom.nmds.species <- data.frame(NMDS1=mybiom.nmds$species[,1],NMDS2=mybiom.nmds$species[,2])
mybiom.ellipse <- ordiellipse(mybiom.nmds,mybiom.meta$Feedstock,display="sites",kind="sd",conf=0.95,label=T)

# Create a function to get the data used for creating the ellipses
# Taken from http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
veganCovEllipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

ellipse.data <- data.frame()
for(g in levels(mybiom.nmds.sites$Feedstock)){
  ellipse.data <- rbind(ellipse.data, cbind(as.data.frame(with(mybiom.nmds.sites[mybiom.nmds.sites$Feedstock==g,],
                  veganCovEllipse(mybiom.ellipse[[g]]$cov,mybiom.ellipse[[g]]$center,mybiom.ellipse[[g]]$scale)))
                  ,Feedstock=g))
}
rm(g)
# Now lets construct the plot by layer
mybiom.plot <- ggplot(data=mybiom.nmds.sites, aes(NMDS1, NMDS2)) +
               theme_minimal() +
               geom_point(aes(color=Temp, shape=Feedstock), size=5) + # Plots the sample points and colors according to temp, point according to feedstock.
               geom_path(data=ellipse.data,aes(x=NMDS1, y=NMDS2,linetype=Feedstock), size=0.75) + # Plots the ellipses, using different line types for the different feedstocks.
               geom_point(data=mybiom.nmds.species, color='grey10', size=0.75) # Plots the species as small grey dots to give more context.
# And now display it
mybiom.plot 

