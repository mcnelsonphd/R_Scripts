# R Scripts

These are a collection of R scripts to perform various ecological analyses of 16S datasets. They are still a work in progress and are in no way guaranteed to work correctly.

It is recommended that users use R Studio 

###Biom_to_NMDS_plots.R
An R "script" that takes a users biom formatted OTU table or taxonomy table and runs NMDS, plotting the results with ggplot. This script assumes that sample meta-data has been added to the biom file, which is default action of the Qiime_Process script since version 13 (2015-03-12). Meta-data is not added by default to the taxonomy biom files.

In order to take advantage of sample meta-data for plotting and statistical grouping analyses, the user will need to modify the script to note the meta-data fields that they wish to use. At present, only two fields are utilized for determining plotting characters and color.

__Required R packages:__ biom, ggplot2, vegan