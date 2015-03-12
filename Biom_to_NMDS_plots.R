require(biom)
# Read in the biom file
mybiom <- read_biom("~/Desktop/Sandbox/Qiime_Dev/even_table_meta.biom")
# Get some basic info on the file
mybiom
# Put the species/sites data into a matrix and then transpose it to put into the correct "orientation"
mybiom.matrix <- as(biom_data(mybiom), "matrix")
mybiom.matrixt <- t(mybiom.matrix)
# Create a sample metadata matrix for easy referencing
mybiom.meta <- as.data.frame(sample_metadata(mybiom))


# Now for the fun, run NMDS on the observation matrix
mybiom.nmds <- metaMDS(mybiom.matrixt, k=2, trymax=500)
mybiom.nmds
stressplot(mybiom.nmds)

# Create a quick and dirty plot for initial analysis
ordiplot(mybiom.nmds)

mybiom.nmds.data <- data.frame(MDS1=mybiom.nmds$points[,1],MDS2=mybiom.nmds$points[,2],group=mybiom.meta$Feed)
ord <- ordiellipse(NB.nmds,NB.env$FS,display="sites",kind="sd",conf=0.95,label=T)

