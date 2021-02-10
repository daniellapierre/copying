##Plotting Script##
##Uses environment from main script##
#source("plots.R")

### Density Histograms - Geographic Distance ###

#all sites#
layout(matrix(1,1,1,byrow = T))
plotDensityHistogram(geodist_global,prob = T,adjust = 1.5,main = "Pairwise Intersite Distances - All Sites",xlab="Distance (km)")

#all sites by region#
layout(matrix(c(1,2,3,4,5),5,1,byrow = T))
plotDensityHistogram(jeff_geodist,prob = T,adjust = 1,main = "Pairwise Intersite Distances - Jefferson County",xlab="Distance (km)")
plotDensityHistogram(erie_geodist,prob = T,adjust = 1,main = "Pairwise Intersite Distances - Erie / Niagara",xlab="Distance (km)")
plotDensityHistogram(ny_geodist,prob = T,adjust = 1,main = "Pairwise Intersite Distances - New York",xlab="Distance (km)")
plotDensityHistogram(on_geodist,prob = T,adjust = 1,main = "Pairwise Intersite Distances - Ontario",xlab="Distance (km)")
plotDensityHistogram(stl_geodist,prob = T,adjust = 1,main = "Pairwise Intersite Distances - St. Lawrence",xlab="Distance (km)")

#all regions by timeblock#
layout(matrix(c(1,2,3,4,5),5,1,byrow = T))
plotDensityHistogram(timeblock1.geodist,prob = T,adjust = 1,main = "Pairwise Intersite Distances - Timeblock 1 (1350-1450)",xlab="Distance (km)")
plotDensityHistogram(timeblock2.geodist,prob = T,adjust = 1,main = "Pairwise Intersite Distances - Timeblock 2 (1400-1500)",xlab="Distance (km)")
plotDensityHistogram(timeblock3.geodist,prob = T,adjust = 1,main = "Pairwise Intersite Distances - Timeblock 3 (1450-1550)",xlab="Distance (km)")
plotDensityHistogram(timeblock4.geodist,prob = T,adjust = 1,main = "Pairwise Intersite Distances - Timeblock 4 (1500-1600)",xlab="Distance (km)")
plotDensityHistogram(timeblock5.geodist,prob = T,adjust = 1,main = "Pairwise Intersite Distances - Timeblock 5 (1550-1650)",xlab="Distance (km)")

#regions by timeblock#
layout(matrix(c(1,2),2,1,byrow = T))
plotDensityHistogram(jeff_geodist.tb2,prob = T,adjust = 1,main = "Pairwise Intersite Distances - Jefferson County Timeblock 2 (1400-1500)",xlab="Distance (km)")
plotDensityHistogram(jeff_geodist.tb3,prob = T,adjust = 1,main = "Pairwise Intersite Distances - Jefferson County Timeblock 3 (1450-1550)",xlab="Distance (km)")

layout(matrix(c(1,2),2,1,byrow = T))
plotDensityHistogram(erie_geodist.tb4,prob = T,adjust = 1,main = "Pairwise Intersite Distances - Erie / Niagara Timeblock 4 (1500-1600)",xlab="Distance (km)")
plotDensityHistogram(erie_geodist.tb5,prob = T,adjust = 1,main = "Pairwise Intersite Distances - Erie / Niagara Timeblock 5 (1550-1650)",xlab="Distance (km)")

layout(matrix(c(1,2,3,4,5),5,1,byrow = T))
plotDensityHistogram(ny_geodist.tb1,prob = T,adjust = 1,main = "Pairwise Intersite Distances - New York Timeblock 1 (1350-1450)",xlab="Distance (km)")
plotDensityHistogram(ny_geodist.tb2,prob = T,adjust = 1,main = "Pairwise Intersite Distances - New York Timeblock 2 (1400-1500)",xlab="Distance (km)")
plotDensityHistogram(ny_geodist.tb3,prob = T,adjust = 1,main = "Pairwise Intersite Distances - New York Timeblock 3 (1450-1550)",xlab="Distance (km)")
plotDensityHistogram(ny_geodist.tb4,prob = T,adjust = 1,main = "Pairwise Intersite Distances - New York Timeblock 4 (1500-1600)",xlab="Distance (km)")
plotDensityHistogram(ny_geodist.tb5,prob = T,adjust = 1,main = "Pairwise Intersite Distances - New York Timeblock 5 (1550-1650)",xlab="Distance (km)")

layout(matrix(c(1,2,3,4,5),5,1,byrow = T))
plotDensityHistogram(on_geodist.tb1,prob = T,adjust = 1,main = "Pairwise Intersite Distances - Ontario Timeblock 1 (1350-1450)",xlab="Distance (km)")
plotDensityHistogram(on_geodist.tb2,prob = T,adjust = 1,main = "Pairwise Intersite Distances - Ontario Timeblock 2 (1400-1500)",xlab="Distance (km)")
plotDensityHistogram(on_geodist.tb3,prob = T,adjust = 1,main = "Pairwise Intersite Distances - Ontario Timeblock 3 (1450-1550)",xlab="Distance (km)")
plotDensityHistogram(on_geodist.tb4,prob = T,adjust = 1,main = "Pairwise Intersite Distances - Ontario Timeblock 4 (1500-1600)",xlab="Distance (km)")
plotDensityHistogram(on_geodist.tb5,prob = T,adjust = 1,main = "Pairwise Intersite Distances - Ontario Timeblock 5 (1550-1650)",xlab="Distance (km)")

layout(matrix(c(1,2,3),3,1,byrow = T))
plotDensityHistogram(stl_geodist.tb1,prob = T,adjust = 1,main = "Pairwise Intersite Distances - St. Lawrence Timeblock 1 (1350-1450)",xlab="Distance (km)")
plotDensityHistogram(stl_geodist.tb2,prob = T,adjust = 1,main = "Pairwise Intersite Distances - St. Lawrence Timeblock 2 (1400-1500)",xlab="Distance (km)")
plotDensityHistogram(stl_geodist.tb3,prob = T,adjust = 1,main = "Pairwise Intersite Distances - St. Lawrence Timeblock 3 (1450-1500)",xlab="Distance (km)")

### Density Histograms - BR Distance ###

#all sites
layout(matrix(1,1,1,byrow = T))
plotDensityHistogram(brdist_global,prob = T,adjust = 1.5,main = "Pairwise Interassemblage Similarity - All Sites",xlab=expression("BR"["Dist"]))

#all timeblocks
layout(matrix(c(1,2,3,4,5),5,1,byrow = T))
plotDensityHistogram(timeblock1.brdist,prob = T,adjust = 1,main = "Pairwise Interassemblage Similarity - Timeblock 1 (1350-1450)",xlab=expression("BR"["Dist"]))
plotDensityHistogram(timeblock2.brdist,prob = T,adjust = 1,main = "Pairwise Interassemblage Similarity  - Timeblock 2 (1400-1500)",xlab=expression("BR"["Dist"]))
plotDensityHistogram(timeblock3.brdist,prob = T,adjust = 1,main = "Pairwise Interassemblage Similarity - Timeblock 3 (1450-1550)",xlab=expression("BR"["Dist"]))
plotDensityHistogram(timeblock4.brdist,prob = T,adjust = 1,main = "Pairwise Interassemblage Similarity - Timeblock 4 (1500-1600)",xlab=expression("BR"["Dist"]))
plotDensityHistogram(timeblock5.brdist,prob = T,adjust = 1,main = "Pairwise Interassemblage Similarity - Timeblock 5 (1550-1650)",xlab=expression("BR"["Dist"]))

#all regions
layout(matrix(c(1,2,3,4,5),5,1,byrow = T))
plotDensityHistogram(jeff_brdist,prob = T,adjust = 1,main = "Pairwise Interassemblage Similarity - Jefferson County",xlab=expression("BR"["Dist"]))
plotDensityHistogram(erie_brdist,prob = T,adjust = 1,main = "Pairwise Interassemblage Similarity  - Erie / Niagara",xlab=expression("BR"["Dist"]))
plotDensityHistogram(ny_brdist,prob = T,adjust = 1,main = "Pairwise Interassemblage Similarity - New York",xlab=expression("BR"["Dist"]))
plotDensityHistogram(on_brdist,prob = T,adjust = 1,main = "Pairwise Interassemblage Similarity - Ontario",xlab=expression("BR"["Dist"]))
plotDensityHistogram(stl_brdist,prob = T,adjust = 1,main = "Pairwise Interassemblage Similarity - St. Lawrence",xlab=expression("BR"["Dist"]))

### Mantel Correlograms ###

#all sites#
layout(matrix(1,1,1,byrow = T))
plot(global_mgram)
title(main="Mantel Correlogram - All Sites")  

#all sites by timeblock#
layout(matrix(c(1,2,3,4,5),5,1,byrow = T))
plot(tb1.mgram)
title(main="Mantel Correlogram - Timeblock 1 (1350-1450)")
plot(tb2.mgram)
title(main="Mantel Correlogram - Timeblock 2 (1400-1500)")
plot(tb3.mgram)
title(main="Mantel Correlogram - Timeblock 3 (1450-1550)")
plot(tb4.mgram)
title(main="Mantel Correlogram - Timeblock 4 (1500-1600)")
plot(tb5.mgram)
title(main="Mantel Correlogram - Timeblock 5 (1550-1650)")

#regions all timeblocks#
layout(matrix(c(1,2,3,4,5),5,1,byrow = T))
plot(jeff.mgram)
title(main="Mantel Correlogram - Jefferson County")
plot(erie.mgram)
title(main="Mantel Correlogram - Erie / Niagara")
plot(ny.mgram)
title(main="Mantel Correlogram - New York")
plot(on.mgram)
title(main="Mantel Correlogram - Ontario")
plot(stl.mgram)
title(main="Mantel Correlogram - St. Lawrence")

#region by timeblock#
layout(matrix(c(1,2),2,1,byrow = T))
plot(jeff.tb2.mgram)
title(main="Mantel Correlogram - Jefferson County Timeblock 2 (1400-1500)")
plot(jeff.tb3.mgram)
title(main="Mantel Correlogram - Jefferson County Timeblock 3 (1450-1550)")

layout(matrix(c(1,2),2,1,byrow = T))
plot(erie.tb4.mgram)
title(main="Mantel Correlogram - Erie / Niagara Timeblock 4 (1500-1600)")
plot(erie.tb5.mgram)
title(main="Mantel Correlogram - Erie / Niagara Timeblock 5 (1550-1650)")

layout(matrix(c(1,2,3,4,5),5,1,byrow = T))
plot(ny.tb1.mgram)
title(main="Mantel Correlogram - New York Timeblock 1 (1350-1450)")
plot(ny.tb2.mgram)
title(main="Mantel Correlogram - New York Timeblock 2 (1400-1550)")
plot(ny.tb3.mgram)
title(main="Mantel Correlogram - New York Timeblock 3 (1450-1550)")
plot(ny.tb4.mgram)
title(main="Mantel Correlogram - New York Timeblock 4 (1500-1600)")
plot(ny.tb5.mgram)
title(main="Mantel Correlogram - New York Timeblock 5 (1550-1650)")

layout(matrix(c(1,2,3,4,5),5,1,byrow = T))
plot(on.tb1.mgram)
title(main="Mantel Correlogram - Ontario Timeblock 1 (1350-1450)")
plot(on.tb2.mgram)
title(main="Mantel Correlogram - Ontario Timeblock 2 (1400-1550)")
plot(on.tb3.mgram)
title(main="Mantel Correlogram - Ontario Timeblock 3 (1450-1550)")
plot(on.tb4.mgram)
title(main="Mantel Correlogram - Ontario Timeblock 4 (1500-1600)")
plot(on.tb5.mgram)
title(main="Mantel Correlogram - Ontario Timeblock 5 (1550-1650)")

layout(matrix(c(1,2,3),3,1,byrow = T))
plot(stl.tb1.mgram)
title(main="Mantel Correlogram - St. Lawrence Timeblock 1 (1350-1450)")
plot(stl.tb2.mgram)
title(main="Mantel Correlogram - St. Lawrence Timeblock 2 (1400-1500)")
plot(stl.tb3.mgram)
title(main="Mantel Correlogram - St. Lawrence Timeblock 3 (1450-1550)")

##BR Boxplots JC vs EN vs NY vs ON vs STL by Timeblock##
BR <- read_csv("brbycgroup.csv",col_names = T) #Table of BR corrected similarity coefficients labelled by combined group and timeblock

#Create plot
BR %>% 
  ggplot(data = BR,mapping = aes(x=TB,y=,log(BR),fill=CG)) + #X axis = Timeblock, Y axis = log(BR)
  labs(y="log(BR)", x = "Timeblock") +
  geom_boxplot2(width.errorbar = 0.1,position = "dodge2") + #install_github("kongdd/Ipaper") better version of geom_boxplot
  scale_x_discrete(labels=c("1350-1450","1400-1500","1450-1550","1500-1600","1550-1650")) +
  scale_fill_discrete(name = "Regional Group", labels = c("EN","JC","NY","ON","STL")) +
  theme_classic()

###Corrected BR Coefficients
#BR coefficient correction code by GM Alberti
# for (s1 in 1:rd) {
#   for (s2 in 1:rd) {
#     zero.categ.a <- length(which(x[s1, ] == 0)) #number of unrepresented categories in site 1
#     zero.categ.b <- length(which(x[s2, ] == 0)) #number of unrepresented categories in site 2
#     joint.absence <- sum(colSums(rbind(x[s1, ], x[s2, ])) == 0) #number of shared absennces (which do not effect BR)
#     if (zero.categ.a == zero.categ.b) { 
#     divisor.final <- 1 #if shared absences are equal BR is not affected
#     }
#     else {
#       divisor.final <- max(zero.categ.a, zero.categ.b) - joint.absence + 0.5 # Highest number of unrepresented categories (NOT count size of category) in either site 1 or 2 - number of shared absences
#     }
#     results[s1, s2] <- round((1 - (sum(abs(x[s1, 
#     ]/sum(x[s1, ]) - x[s2, ]/sum(x[s2, ]))))/2)/divisor.final, #normal BR coefficient formula / number of UNshared absences. GM Alberti has demonstration of correction on his webpage cainarchaeology.weebly.com
#     digits = 3)
#   }
# }
# #correction removes false positives from shared pcmcs with large sample sizes between sites (especially simples PCMC #3)
# library(GmAMisc)
# library(reshape2)
# library(tidyverse)
# 
# #all_pcmcs <- dt %>% select(m2:m30) From R Script
# testsample <- all_pcmcs %>% 
#   sample_frac(1/3,replace = T) #take sample of all counts
# test.corr <- BRsim(testsample,correction = T,clust = F,rescale = T) #corrected sample
# test.uncorr <- BRsim(testsample,correction = F,clust = F,rescale = T) #uncorrected sample
# 
# testdata <- data.frame("corrected" = as.vector(test.corr$BR_similarity_matrix),"uncorrected" = as.vector(test.uncorr$BR_similarity_matrix))
# 
# #test if the two samples are significantly different
# perm.t.test(testdata,"short","corrected","uncorrected",10000) 
# #2-sided p-value 0.0002 shows the samples are very different
# #on average uncorrected BR coefficients are 0.32 larger than corrected ones
# #mean corrected BR = 0.1 mean uncorrected BR = 0.42 +- 0.01
# #plot the distributions
# testdata %>%
#   melt() %>%
#   ggplot(aes(x=value,color = variable)) +
#   geom_freqpoly(size=2,bins=25) +
#   labs(title = "Corrected vs. Uncorrected BR Coefficients",y="Frequency",x="BR Similarity Coefficient") +
#   scale_color_discrete(name = "Test Sample",labels = c("Corrected","Uncorrected")) +
#   theme(plot.title=element_text(hjust = 0.5))
# #The corrected coefficient is much more conservative in assigning similarity to two sites, which is important because the dataset is sparse with many unshared categories.
