############################################################################################################################
############################################################################################################################
##Community Coalescence and Regional Geospatial Trends of Ceramic Decorative Variation in Late Woodland Northern Iroquoia ##
############################################################################################################################
###################################################Daniel LaPierre##########################################################
############################################################################################################################
rm(list=ls())
set.seed(133769420)
source("utilityfunctions_Crema.R") #binaryCulture and distMat functions and pairwise PhiST for-loops written by Enrico Crema, Shennan et al. 2015
library(GmAMisc) #BRsim function and BR coefficient correction equation written by GM Alberti (2020)
library(usedist)
library(tidyverse)
library(lemon)
library(rcompanion)
library(vegan)
library(pegas)
library(Ipaper) #boxplot2
library(otuSummary) #matrix collapse tool

#Load data
dt <- read_tsv("SiteData.txt",col_types = cols( #define variable types so sorting works properly
  siteid = "i",
  sitename = "c",
  cgroup_code = col_factor(c("A","B","C","D","E")),
  groupcode = col_factor(c("1","2","3","4","5","6","7","8","9","10","11","12","13")),
  timeslicecode = col_factor(c("I","II","III","IV","V","VI")),
  edate = "i",
  ldate = "i",
  m1 = "i",
  m2 = "i",
  m3 = "i",
  m4 = "i",
  m5 = "i",
  m6 = "i",
  m7 = "i",
  m8 = "i",
  m9 = "i",
  m10 = "i",
  m11 = "i",
  m12 = "i",
  m13 = "i",
  m14 = "i",
  m15 = "i",
  m16 = "i",
  m17 = "i",
  m18 = "i",
  m19 = "i",
  m20 = "i",
  m21 = "i",
  m22 = "i",
  m23 = "i",
  m24 = "i",
  m25 = "i",
  m26 = "i",
  m27 = "i",
  m28 = "i",
  m29 = "i",
  m30 = "i",
  totalnoplain = "i",
  ncombo = "i"),trim_ws = TRUE)

#Set Rownames/subgroups
dt$sitename <- abbreviate(dt$sitename,minlength = 4,use.classes = F,dot = F) #abbreviate/trim site names
dt <- dt %>%  column_to_rownames("sitename") #add rownames to df
ny.dt <- dt[dt$cgroup_code=="C",] #subset NY region
on.dt <- dt[dt$cgroup_code=="D",] #subset ON region
stl.dt <- dt[dt$cgroup_code=="E",] #subset STL region
all_pcmcs <- dt %>% select(m2:m30) #select pcmc counts only

#Load All Sites Distance Matrix
geodist_global <- readRDS("geodistglobal.rds")
geodist_global <- dist_setNames(geodist_global,rownames(dt)) #set rownames of global matrix based on dt for subsetting exactly as above

#Subset For Case Scenarios
### Regions ###
jeff_geodist <- dist_subset(geodist_global,rownames(dt[dt$cgroup_code=="A",]))
erie_geodist <- dist_subset(geodist_global,rownames(dt[dt$cgroup_code=="B",]))
ny_geodist <- dist_subset(geodist_global,rownames(dt[dt$cgroup_code=="C",]))
on_geodist <- dist_subset(geodist_global,rownames(dt[dt$cgroup_code=="D",]))
stl_geodist <- dist_subset(geodist_global,rownames(dt[dt$cgroup_code=="E",]))

### TIME PERIODS ###
timeblock1.geodist <- dist_subset(geodist_global,rownames(dt[dt$timeslicecode == "I" | dt$timeslicecode == "II",]))
timeblock2.geodist <- dist_subset(geodist_global,rownames(dt[dt$timeslicecode == "II" | dt$timeslicecode == "III",]))
timeblock3.geodist <- dist_subset(geodist_global,rownames(dt[dt$timeslicecode == "III" | dt$timeslicecode == "IV",]))
timeblock4.geodist <- dist_subset(geodist_global,rownames(dt[dt$timeslicecode == "IV" | dt$timeslicecode == "V",]))
timeblock5.geodist <- dist_subset(geodist_global,rownames(dt[dt$timeslicecode == "V" | dt$timeslicecode == "VI",]))

### Regions by Time Period ###
##Jefferson
jeff_geodist.tb2 <- dist_subset(jeff_geodist,rownames(dt[dt$cgroup_code=="A" & dt$timeslicecode %in% c("II", "III"),]))
jeff_geodist.tb3 <- dist_subset(jeff_geodist,rownames(dt[dt$cgroup_code=="A" & dt$timeslicecode %in% c("III", "IV"),]))
##Erie
erie_geodist.tb4 <- dist_subset(erie_geodist,rownames(dt[dt$cgroup_code=="B" & dt$timeslicecode %in% c("IV", "V"),]))
erie_geodist.tb5 <- dist_subset(erie_geodist,rownames(dt[dt$cgroup_code=="B" & dt$timeslicecode %in% c("V", "VI"),]))
#New York
ny_geodist.tb1 <- dist_subset(ny_geodist,rownames(ny.dt[ny.dt$timeslicecode == "I" | ny.dt$timeslicecode == "II",]))
ny_geodist.tb2 <- dist_subset(ny_geodist,rownames(ny.dt[ny.dt$timeslicecode == "II" | ny.dt$timeslicecode == "III",]))
ny_geodist.tb3 <- dist_subset(ny_geodist,rownames(ny.dt[ny.dt$timeslicecode == "III" | ny.dt$timeslicecode == "IV",]))
ny_geodist.tb4 <- dist_subset(ny_geodist,rownames(ny.dt[ny.dt$timeslicecode == "IV" | ny.dt$timeslicecode == "V",]))
ny_geodist.tb5 <- dist_subset(ny_geodist,rownames(ny.dt[ny.dt$timeslicecode == "V" | ny.dt$timeslicecode == "VI",]))
#Ontario
on_geodist.tb1 <- dist_subset(on_geodist,rownames(on.dt[on.dt$timeslicecode == "I" | on.dt$timeslicecode == "II",]))
on_geodist.tb2 <- dist_subset(on_geodist,rownames(on.dt[on.dt$timeslicecode == "II" | on.dt$timeslicecode == "III",]))
on_geodist.tb3 <- dist_subset(on_geodist,rownames(on.dt[on.dt$timeslicecode == "III" | on.dt$timeslicecode == "IV",]))
on_geodist.tb4 <- dist_subset(on_geodist,rownames(on.dt[on.dt$timeslicecode == "IV" | on.dt$timeslicecode == "V",]))
on_geodist.tb5 <- dist_subset(on_geodist,rownames(on.dt[on.dt$timeslicecode == "V" | on.dt$timeslicecode == "VI",]))
#St. Lawrence
stl_geodist.tb1 <- dist_subset(stl_geodist,rownames(dt[dt$cgroup_code=="E" & dt$timeslicecode %in% c("I", "II"),]))
stl_geodist.tb2 <- dist_subset(stl_geodist,rownames(dt[dt$cgroup_code=="E" & dt$timeslicecode %in% c("II", "III"),]))
stl_geodist.tb3 <- dist_subset(stl_geodist,rownames(dt[dt$cgroup_code=="E" & dt$timeslicecode %in% c("III", "IV"),]))

#Generate BR Similarity Coefficients and Distance Matrices
### ALL SITES ###
br_global <- BRsim(all_pcmcs,correction = T,rescale = T,clust = T,oneplot = F) #correction removes unwanted false positives from rich single category matches
brdist_global <- as.dist(m = br_global$BR_distance_matrix,upper = F) #extract distance matrix
brdist_global <- dist_setNames(brdist_global,rownames(dt)) #define the rownames of the matrix, used to subset complete sample into subgroups

### Regions ###
jeff_brdist <- dist_subset(brdist_global,rownames(dt[dt$cgroup_code=="A",])) #subset global matrix by rownames of rows in dt grouped in CG A (Jefferson County)
erie_brdist <- dist_subset(brdist_global,rownames(dt[dt$cgroup_code=="B",])) #subset global matrix by rownames of rows in dt grouped in CG B (Lake Erie Niagara)
ny_brdist <- dist_subset(brdist_global,rownames(dt[dt$cgroup_code=="C",])) #subset global matrix by rownames of rows in dt grouped in CG C (New York)
on_brdist <- dist_subset(brdist_global,rownames(dt[dt$cgroup_code=="D",])) #subset global matrix by rownames of rows in dt grouped in CG D (Ontario)
stl_brdist <- dist_subset(brdist_global,rownames(dt[dt$cgroup_code=="E",])) #subset global matrix by rownames of rows in dt grouped in CG A (St. Lawrence)

### TIME PERIODS ###
timeblock1.brdist <- dist_subset(brdist_global,rownames(dt[dt$timeslicecode == "I" | dt$timeslicecode == "II",])) #subset global matrix by rownames of rows in dt grouped in timeslices I or II
timeblock2.brdist <- dist_subset(brdist_global,rownames(dt[dt$timeslicecode == "II" | dt$timeslicecode == "III",])) #subset global matrix by rownames of rows in dt grouped in timeslices II or III
timeblock3.brdist <- dist_subset(brdist_global,rownames(dt[dt$timeslicecode == "III" | dt$timeslicecode == "IV",])) #subset global matrix by rownames of rows in dt grouped in timeslices III or IV
timeblock4.brdist <- dist_subset(brdist_global,rownames(dt[dt$timeslicecode == "IV" | dt$timeslicecode == "V",])) #subset global matrix by rownames of rows in dt grouped in timeslices IV or V
timeblock5.brdist <- dist_subset(brdist_global,rownames(dt[dt$timeslicecode == "V" | dt$timeslicecode == "VI",])) #subset global matrix by rownames of rows in dt grouped in timeslices V or VI

### Regions by Time Period ###
##Jefferson
jeff_brdist.tb2 <- dist_subset(jeff_brdist,rownames(dt[dt$cgroup_code=="A" & dt$timeslicecode %in% c("II", "III"),])) #subset global matrix by rownames of rows in dt grouped in combined group A AND timeslices II or III
jeff_brdist.tb3 <- dist_subset(jeff_brdist,rownames(dt[dt$cgroup_code=="A" & dt$timeslicecode %in% c("III", "IV"),])) #subset global matrix by rownames of rows in dt grouped in combined group A AND timeslices III or IV
##Erie
erie_brdist.tb4 <- dist_subset(erie_brdist,rownames(dt[dt$cgroup_code=="B" & dt$timeslicecode %in% c("IV", "V"),])) #subset global matrix by rownames of rows in dt grouped in combined group B AND timeslices IV or V
erie_brdist.tb5 <- dist_subset(erie_brdist,rownames(dt[dt$cgroup_code=="B" & dt$timeslicecode %in% c("V", "VI"),])) #subset global matrix by rownames of rows in dt grouped in combined group B AND timeslices V or VI
#New York
ny_brdist.tb1 <- dist_subset(ny_brdist,rownames(ny.dt[ny.dt$timeslicecode == "I" | ny.dt$timeslicecode == "II",])) #subset global matrix by rownames of rows in ny.dt grouped in timeslices I or II
ny_brdist.tb2 <- dist_subset(ny_brdist,rownames(ny.dt[ny.dt$timeslicecode == "II" | ny.dt$timeslicecode == "III",])) #subset global matrix by rownames of rows in ny.dt grouped in timeslices II or III
ny_brdist.tb3 <- dist_subset(ny_brdist,rownames(ny.dt[ny.dt$timeslicecode == "III" | ny.dt$timeslicecode == "IV",])) #subset global matrix by rownames of rows in ny.dt grouped in timeslices III or IV
ny_brdist.tb4 <- dist_subset(ny_brdist,rownames(ny.dt[ny.dt$timeslicecode == "IV" | ny.dt$timeslicecode == "V",])) #subset global matrix by rownames of rows in ny.dt grouped in timeslices IV or V
ny_brdist.tb5 <- dist_subset(ny_brdist,rownames(ny.dt[ny.dt$timeslicecode == "V" | ny.dt$timeslicecode == "VI",])) #subset global matrix by rownames of rows in ny.dt grouped in timeslices V or VI
#Ontario
on_brdist.tb1 <- dist_subset(on_brdist,rownames(on.dt[on.dt$timeslicecode == "I" | on.dt$timeslicecode == "II",])) #subset global matrix by rownames of rows in on.dt grouped in timeslices I or II
on_brdist.tb2 <- dist_subset(on_brdist,rownames(on.dt[on.dt$timeslicecode == "II" | on.dt$timeslicecode == "III",])) #subset global matrix by rownames of rows in on.dt grouped in timeslices II or III
on_brdist.tb3 <- dist_subset(on_brdist,rownames(on.dt[on.dt$timeslicecode == "III" | on.dt$timeslicecode == "IV",])) #subset global matrix by rownames of rows in on.dt grouped in timeslices III or IV
on_brdist.tb4 <- dist_subset(on_brdist,rownames(on.dt[on.dt$timeslicecode == "IV" | on.dt$timeslicecode == "V",])) #subset global matrix by rownames of rows in on.dt grouped in timeslices IV or V
on_brdist.tb5 <- dist_subset(on_brdist,rownames(on.dt[on.dt$timeslicecode == "V" | on.dt$timeslicecode == "VI",])) #subset global matrix by rownames of rows in on.dt grouped in timeslices V or VI
#St. Lawrence
stl_brdist.tb1 <- dist_subset(stl_brdist,rownames(dt[dt$cgroup_code=="E" & dt$timeslicecode %in% c("I", "II"),])) #subset global matrix by rownames of rows in combined group E AND grouped in timeslices I or II
stl_brdist.tb2 <- dist_subset(stl_brdist,rownames(dt[dt$cgroup_code=="E" & dt$timeslicecode %in% c("II", "III"),])) #subset global matrix by rownames of rows in combined group E AND grouped in timeslices II or III
stl_brdist.tb3 <- dist_subset(stl_brdist,rownames(dt[dt$cgroup_code=="E" & dt$timeslicecode %in% c("III", "IV"),])) #subset global matrix by rownames of rows in combined group E AND grouped in timeslices III or IV

###################################
###################################
#####Step 1: Mantel Tests##########
###################################
###################################

### ALL SITES ###
mantel(brdist_global,geodist_global,method = "spearman",permutations = 10000,parallel = 4) #spearman rho to handle skewed nonnormal archaeological data

### Regions All Time Periods ###
mantel(jeff_brdist,jeff_geodist,method = "spearman",permutations = 10000,parallel = 4)
mantel(erie_brdist,erie_geodist,method = "spearman",permutations = 10000,parallel = 4)
mantel(ny_brdist,ny_geodist,method = "spearman",permutations = 10000,parallel = 4)
mantel(on_brdist,on_geodist,method = "spearman",permutations = 10000,parallel = 4)
mantel(stl_brdist,stl_geodist,method = "spearman",permutations = 10000,parallel = 4)

### Time Periods All Regions ###
mantel(timeblock1.brdist,timeblock1.geodist,method = "spearman",permutations = 10000,parallel = 4)
mantel(timeblock2.brdist,timeblock2.geodist,method = "spearman",permutations = 10000,parallel = 4)
mantel(timeblock3.brdist,timeblock3.geodist,method = "spearman",permutations = 10000,parallel = 4)
mantel(timeblock4.brdist,timeblock4.geodist,method = "spearman",permutations = 10000,parallel = 4)
mantel(timeblock5.brdist,timeblock5.geodist,method = "spearman",permutations = 10000,parallel = 4)

### Regions by Time Period ###
mantel(jeff_brdist.tb2,jeff_geodist.tb2,method = "spearman",permutations = 10000,parallel = 4)
mantel(jeff_brdist.tb3,jeff_geodist.tb3,method = "spearman",permutations = 10000,parallel = 4)

mantel(erie_brdist.tb4,erie_geodist.tb4,method = "spearman",permutations = 10000,parallel = 4)
mantel(erie_brdist.tb5,erie_geodist.tb5,method = "spearman",permutations = 10000,parallel = 4)

mantel(ny_brdist.tb1,ny_geodist.tb1,method = "spearman",permutations = 10000,parallel = 4)
mantel(ny_brdist.tb2,ny_geodist.tb2,method = "spearman",permutations = 10000,parallel = 4)
mantel(ny_brdist.tb3,ny_geodist.tb3,method = "spearman",permutations = 10000,parallel = 4)
mantel(ny_brdist.tb4,ny_geodist.tb4,method = "spearman",permutations = 10000,parallel = 4)
mantel(ny_brdist.tb5,ny_geodist.tb5,method = "spearman",permutations = 10000,parallel = 4)

mantel(on_brdist.tb1,on_geodist.tb1,method = "spearman",permutations = 10000,parallel = 4)
mantel(on_brdist.tb2,on_geodist.tb2,method = "spearman",permutations = 10000,parallel = 4)
mantel(on_brdist.tb3,on_geodist.tb3,method = "spearman",permutations = 10000,parallel = 4)
mantel(on_brdist.tb4,on_geodist.tb4,method = "spearman",permutations = 10000,parallel = 4)
mantel(on_brdist.tb5,on_geodist.tb5,method = "spearman",permutations = 10000,parallel = 4)

mantel(stl_brdist.tb1,stl_geodist.tb1,method = "spearman",permutations = 10000,parallel = 4)
mantel(stl_brdist.tb2,stl_geodist.tb2,method = "spearman",permutations = 10000,parallel = 4)
mantel(stl_brdist.tb3,stl_geodist.tb3,method = "spearman",permutations = 10000,parallel = 4)

### Mantel Correlograms ###
brkpts = seq(from=0,to=1100,by=10) #define sequence of 10 km distance classes to cover all possible distances in global sample (slow but thorough)

#All Sites#
global_mgram <- mantel.correlog(brdist_global,geodist_global,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000) #cutoff = F ensures all possible pairs of sites are compared 

#all sites by time period#
tb1.mgram <- mantel.correlog(timeblock1.brdist,timeblock1.geodist,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm=10000)
tb2.mgram <- mantel.correlog(timeblock2.brdist,timeblock2.geodist,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm=10000)
tb3.mgram <- mantel.correlog(timeblock3.brdist,timeblock3.geodist,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm=10000)
tb4.mgram <- mantel.correlog(timeblock4.brdist,timeblock4.geodist,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm=10000)
tb5.mgram <- mantel.correlog(timeblock5.brdist,timeblock5.geodist,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm=10000)

#regions all time periods#
jeff.mgram <- mantel.correlog(jeff_brdist,jeff_geodist,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)
erie.mgram <- mantel.correlog(erie_brdist,erie_geodist,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)
ny.mgram <- mantel.correlog(ny_brdist,ny_geodist,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)
on.mgram <- mantel.correlog(on_brdist,on_geodist,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)
stl.mgram <- mantel.correlog(stl_brdist,stl_geodist,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)

#region by time period#
jeff.tb2.mgram <- mantel.correlog(jeff_brdist.tb2,jeff_geodist.tb2,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)
jeff.tb3.mgram <- mantel.correlog(jeff_brdist.tb3,jeff_geodist.tb3,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)

erie.tb4.mgram <- mantel.correlog(erie_brdist.tb4,erie_geodist.tb4,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)
erie.tb5.mgram <- mantel.correlog(erie_brdist.tb5,erie_geodist.tb5,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)

ny.tb1.mgram <- mantel.correlog(ny_brdist.tb1,ny_geodist.tb1,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)
ny.tb2.mgram <- mantel.correlog(ny_brdist.tb2,ny_geodist.tb2,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)
ny.tb3.mgram <- mantel.correlog(ny_brdist.tb3,ny_geodist.tb3,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)
ny.tb4.mgram <- mantel.correlog(ny_brdist.tb4,ny_geodist.tb4,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)
ny.tb5.mgram <- mantel.correlog(ny_brdist.tb5,ny_geodist.tb5,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)

on.tb1.mgram <- mantel.correlog(on_brdist.tb1,on_geodist.tb1,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)
on.tb2.mgram <- mantel.correlog(on_brdist.tb2,on_geodist.tb2,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)
on.tb3.mgram <- mantel.correlog(on_brdist.tb3,on_geodist.tb3,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)
on.tb4.mgram <- mantel.correlog(on_brdist.tb4,on_geodist.tb4,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)
on.tb5.mgram <- mantel.correlog(on_brdist.tb5,on_geodist.tb5,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)

stl.tb1.mgram <- mantel.correlog(stl_brdist.tb1,stl_geodist.tb1,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)
stl.tb2.mgram <- mantel.correlog(stl_brdist.tb2,stl_geodist.tb2,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)
stl.tb3.mgram <- mantel.correlog(stl_brdist.tb3,stl_geodist.tb3,break.pts = brkpts,cutoff = F,r.type = "spearman",nperm = 10000)

###################################
###################################
########Step 2: AMOVAs#############
###################################
###################################

### AMOVA squares distances so sqrt() maintains BR coefficients (sqrt(BR)^2 = BR) ###
sqbrdist_global <- sqrt(brdist_global) #all sites

timeblock1.sqbrdist <- sqrt(timeblock1.brdist) #timeblock 1
timeblock2.sqbrdist <- sqrt(timeblock2.brdist) #timeblock 2
timeblock3.sqbrdist <- sqrt(timeblock3.brdist) #timeblock 3
timeblock4.sqbrdist <- sqrt(timeblock4.brdist) #timeblock 4
timeblock5.sqbrdist <- sqrt(timeblock5.brdist) #timeblock 5

ny_sqbrdist <- sqrt(ny_brdist) #new york all sites
on_sqbrdist <- sqrt(on_brdist) #ontario all sites
stl_sqbrdist <- sqrt(stl_brdist) #st lawrence all sites

ny_sqbrdist.tb1 <- sqrt(ny_brdist.tb1) #new york tb1
ny_sqbrdist.tb2 <- sqrt(ny_brdist.tb2) #new york tb2
ny_sqbrdist.tb3 <- sqrt(ny_brdist.tb3) #new york tb3
ny_sqbrdist.tb4 <- sqrt(ny_brdist.tb4) #new york tb4
ny_sqbrdist.tb5 <- sqrt(ny_brdist.tb5) #new york tb5

on_sqbrdist.tb1 <- sqrt(on_brdist.tb1) #ontario tb1
on_sqbrdist.tb2 <- sqrt(on_brdist.tb2) #ontario tb2
on_sqbrdist.tb3 <- sqrt(on_brdist.tb3) #ontario tb3
on_sqbrdist.tb4 <- sqrt(on_brdist.tb4) #ontario tb4
on_sqbrdist.tb5 <- sqrt(on_brdist.tb5) #ontario tb5

stl_sqbrdist.tb1 <- sqrt(stl_brdist.tb1) #stl tb1
stl_sqbrdist.tb2 <- sqrt(stl_brdist.tb2) #stl tb2
stl_sqbrdist.tb3 <- sqrt(stl_brdist.tb3) #stl tb3

### Create Grouping Factors ###
groups_global <-  as.factor(dt$groupcode) #by spatial groups
cgroups_global <-  as.factor(dt$cgroup_code) #by regional groups

groups_ny <- as.factor(dt$groupcode[dt$cgroup_code=="C"]) #spatial groups in NY all timeblocks
groups_on <- as.factor(dt$groupcode[dt$cgroup_code=="D"]) #spatial groups in ON all timeblocks
groups_stl <- as.factor(dt$groupcode[dt$cgroup_code=="E"]) #spatial groups in STL all timeblocks

groups_tb1 <- as.factor(dt$groupcode[dt$timeslicecode == "I" | dt$timeslicecode == "II"]) #all spatial groups timeblock 1
groups_tb2 <- as.factor(dt$groupcode[dt$timeslicecode == "II" | dt$timeslicecode == "III"]) #all spatial groups timeblock 2
groups_tb3 <- as.factor(dt$groupcode[dt$timeslicecode == "III" | dt$timeslicecode == "IV"]) #all spatial groups timeblock 3
groups_tb4 <- as.factor(dt$groupcode[dt$timeslicecode == "IV" | dt$timeslicecode == "V"]) #all spatial groups timeblock 4
groups_tb5 <- as.factor(dt$groupcode[dt$timeslicecode == "V" | dt$timeslicecode == "VI"]) #all spatial groups timeblock 5
cgroups_tb1 <- as.factor(dt$cgroup_code[dt$timeslicecode == "I" | dt$timeslicecode == "II"]) #all regional groups in timeblock 1
cgroups_tb2 <- as.factor(dt$cgroup_code[dt$timeslicecode == "II" | dt$timeslicecode == "III"]) #all regional groups in timeblock 2
cgroups_tb3 <- as.factor(dt$cgroup_code[dt$timeslicecode == "III" | dt$timeslicecode == "IV"]) #all regional groups in timeblock 3
cgroups_tb4 <- as.factor(dt$cgroup_code[dt$timeslicecode == "IV" | dt$timeslicecode == "V"]) #all regional groups in timeblock 4
cgroups_tb5 <- as.factor(dt$cgroup_code[dt$timeslicecode == "V" | dt$timeslicecode == "VI"]) #all regional groups in timeblock 5

groups_ny.tb1 <- as.factor(ny.dt$groupcode[ny.dt$timeslicecode == "I" | ny.dt$timeslicecode == "II"]) #all spatial groups in NY timeblock 1
groups_ny.tb2 <- as.factor(ny.dt$groupcode[ny.dt$timeslicecode == "II" | ny.dt$timeslicecode == "III"]) #all spatial groups in NY timeblock 2
groups_ny.tb3 <- as.factor(ny.dt$groupcode[ny.dt$timeslicecode == "III" | ny.dt$timeslicecode == "IV"]) #all spatial groups in NY timeblock 3
groups_ny.tb4 <- as.factor(ny.dt$groupcode[ny.dt$timeslicecode == "IV" | ny.dt$timeslicecode == "V"]) #all spatial groups in NY timeblock 4
groups_ny.tb5 <- as.factor(ny.dt$groupcode[ny.dt$timeslicecode == "V" | ny.dt$timeslicecode == "VI"]) #all spatial groups in NY timeblock 5

groups_on.tb1 <- as.factor(on.dt$groupcode[on.dt$timeslicecode == "I" | on.dt$timeslicecode == "II"]) #all spatial groups in ON timeblock 1
groups_on.tb2 <- as.factor(on.dt$groupcode[on.dt$timeslicecode == "II" | on.dt$timeslicecode == "III"]) #all spatial groups in ON timeblock 2
groups_on.tb3 <- as.factor(on.dt$groupcode[on.dt$timeslicecode == "III" | on.dt$timeslicecode == "IV"]) #all spatial groups in ON timeblock 3
groups_on.tb4 <- as.factor(on.dt$groupcode[on.dt$timeslicecode == "IV" | on.dt$timeslicecode == "V"]) #all spatial groups in ON timeblock 4
groups_on.tb5 <- as.factor(on.dt$groupcode[on.dt$timeslicecode == "V" | on.dt$timeslicecode == "VI"]) #all spatial groups in ON timeblock 5

groups_stl.tb1 <- as.factor(stl.dt$groupcode[stl.dt$timeslicecode == "I" | stl.dt$timeslicecode == "II"]) #all spatial groups in STL timeblock 1
groups_stl.tb2 <- as.factor(stl.dt$groupcode[stl.dt$timeslicecode == "II" | stl.dt$timeslicecode == "III"]) #all spatial groups in STL timeblock 2
groups_stl.tb3 <- as.factor(stl.dt$groupcode[stl.dt$timeslicecode == "III" | stl.dt$timeslicecode == "IV"]) #all spatial groups in STL timeblock 3


### Single AMOVA ###
AMOVAglobal <- amova(sqbrdist_global~groups_global,nperm = 10000) #all sites by spatial
AMOVAglobal.cgroups <- amova(sqbrdist_global~cgroups_global,nperm = 10000) #all sites by combined group

AMOVAtb1 <- amova(timeblock1.sqbrdist~groups_tb1,nperm = 10000) #timeblock 1 sites by spatial group
AMOVAtb2 <- amova(timeblock2.sqbrdist~groups_tb2,nperm = 10000) #timeblock 2 sites by spatial group
AMOVAtb3 <- amova(timeblock3.sqbrdist~groups_tb3,nperm = 10000) #timeblock 3 sites by spatial group
AMOVAtb4 <- amova(timeblock4.sqbrdist~groups_tb4,nperm = 10000) #timeblock 4 sites by spatial group
AMOVAtb5 <- amova(timeblock5.sqbrdist~groups_tb5,nperm = 10000) #timeblock 5 sites by spatial group
AMOVAtb1.cgroups <- amova(timeblock1.sqbrdist~cgroups_tb1,nperm = 10000) #timeblock 1 sites by combined group
AMOVAtb2.cgroups <- amova(timeblock2.sqbrdist~cgroups_tb2,nperm = 10000) #timeblock 2 sites by combined group
AMOVAtb3.cgroups <- amova(timeblock3.sqbrdist~cgroups_tb3,nperm = 10000) #timeblock 3 sites by combined group
AMOVAtb4.cgroups <- amova(timeblock4.sqbrdist~cgroups_tb4,nperm = 10000) #timeblock 4 sites by combined group
AMOVAtb5.cgroups <- amova(timeblock5.sqbrdist~cgroups_tb5,nperm = 10000) #timeblock 5 sites by combined group

AMOVAny <- amova(ny_sqbrdist~groups_ny,nperm = 10000) #all sites in NY combined group by spatial
AMOVAon <- amova(on_sqbrdist~groups_on,nperm = 10000) #all sites in ON combined group by spatial
AMOVAstl <- amova(stl_sqbrdist~groups_stl,nperm = 10000) #all sites in STL combined group by spatial

AMOVAny.tb1 <- amova(ny_sqbrdist.tb1~groups_ny.tb1,nperm = 10000) #timeblock 1 sites in NY combined group by spatial
AMOVAny.tb2 <- amova(ny_sqbrdist.tb2~groups_ny.tb2,nperm = 10000) #timeblock 2 sites in NY combined group by spatial
AMOVAny.tb3 <- amova(ny_sqbrdist.tb3~groups_ny.tb3,nperm = 10000) #timeblock 3 sites in NY combined group by spatial
AMOVAny.tb4 <- amova(ny_sqbrdist.tb4~groups_ny.tb4,nperm = 10000) #timeblock 4 sites in NY combined group by spatial
AMOVAny.tb5 <- amova(ny_sqbrdist.tb5~groups_ny.tb5,nperm = 10000) #timeblock 5 sites in NY combined group by spatial

AMOVAon.tb1 <- amova(on_sqbrdist.tb1~groups_on.tb1,nperm = 10000)  #timeblock 1 sites in ON combined group by spatial
AMOVAon.tb2 <- amova(on_sqbrdist.tb2~groups_on.tb2,nperm = 10000)  #timeblock 2 sites in ON combined group by spatial
AMOVAon.tb3 <- amova(on_sqbrdist.tb3~groups_on.tb3,nperm = 10000)  #timeblock 3 sites in ON combined group by spatial
AMOVAon.tb4 <- amova(on_sqbrdist.tb4~groups_on.tb4,nperm = 10000)  #timeblock 4 sites in ON combined group by spatial
AMOVAon.tb5 <- amova(on_sqbrdist.tb5~groups_on.tb5,nperm = 10000)  #timeblock 5 sites in ON combined group by spatial

AMOVAstl.tb1 <- amova(stl_sqbrdist.tb1~groups_stl.tb1,nperm = 10000) #timeblock 1 sites in STL combined group by spatial
AMOVAstl.tb2 <- amova(stl_sqbrdist.tb2~groups_stl.tb2,nperm = 10000) #timeblock 2 sites in STL combined group by spatial
AMOVAstl.tb3 <- amova(stl_sqbrdist.tb3~groups_stl.tb3,nperm = 10000) #timeblock 3 sites in STL combined group by spatial

### Single AMOVA Results ###
AMOVAglobal$varcomp$sigma2[1]/sum(AMOVAglobal$varcomp$sigma2) #Phi all sites by spatial group. within groups / across groups
AMOVAglobal$varcomp$P.value[1] #p value all sites by spatial group

AMOVAglobal.cgroups$varcomp$sigma2[1]/sum(AMOVAglobal.cgroups$varcomp$sigma2) #Phi all sites by regional group
AMOVAglobal.cgroups$varcomp$P.value[1] #p value all sites by regional group

AMOVAtb1$varcomp$sigma2[1]/sum(AMOVAtb1$varcomp$sigma2) #Phi all sites by spatial group timeblock 1
AMOVAtb1$varcomp$P.value[1] #p value all sites by spatial group timeblock 1
AMOVAtb2$varcomp$sigma2[1]/sum(AMOVAtb2$varcomp$sigma2) #Phi all sites by spatial group timeblock 2
AMOVAtb2$varcomp$P.value[1] #p value all sites by spatial group timeblock 2
AMOVAtb3$varcomp$sigma2[1]/sum(AMOVAtb3$varcomp$sigma2) #Phi all sites by spatial group timeblock 3
AMOVAtb3$varcomp$P.value[1] #p value all sites timeblock 3
AMOVAtb4$varcomp$sigma2[1]/sum(AMOVAtb4$varcomp$sigma2) #Phi all sites by spatial group timeblock 4
AMOVAtb4$varcomp$P.value[1] #p value all sites by spatial group timeblock 4
AMOVAtb5$varcomp$sigma2[1]/sum(AMOVAtb5$varcomp$sigma2) #Phi all sites by spatial group timeblock 5
AMOVAtb5$varcomp$P.value[1] #p value all sites by spatial group timeblock 5

AMOVAtb1.cgroups$varcomp$sigma2[1]/sum(AMOVAtb1.cgroups$varcomp$sigma2) #Phi all sites by regional group timeblock 1 
AMOVAtb1.cgroups$varcomp$P.value[1] #p value all sites by regional group timeblock 1
AMOVAtb2.cgroups$varcomp$sigma2[1]/sum(AMOVAtb2.cgroups$varcomp$sigma2) #Phi all sites by regional group timeblock 1
AMOVAtb2.cgroups$varcomp$P.value[1] #p value all sites by regional group timeblock 1
AMOVAtb3.cgroups$varcomp$sigma2[1]/sum(AMOVAtb3.cgroups$varcomp$sigma2) #Phi all sites by regional group timeblock 1
AMOVAtb3.cgroups$varcomp$P.value[1] #p value all sites by regional group timeblock 1
AMOVAtb4.cgroups$varcomp$sigma2[1]/sum(AMOVAtb4.cgroups$varcomp$sigma2) #Phi all sites by regional group timeblock 1
AMOVAtb4.cgroups$varcomp$P.value[1] #p value all sites by regional group timeblock 1
AMOVAtb5.cgroups$varcomp$sigma2[1]/sum(AMOVAtb5.cgroups$varcomp$sigma2) #Phi all sites by regional group timeblock 1
AMOVAtb5.cgroups$varcomp$P.value[1] #p value all sites by regional group timeblock 1

AMOVAny$varcomp$sigma2[1]/sum(AMOVAny$varcomp$sigma2) #Phi NY combined group all timeblocks
AMOVAny$varcomp$P.value[1] #p value NY combined group all timeblocks
AMOVAon$varcomp$sigma2[1]/sum(AMOVAon$varcomp$sigma2) #Phi ON combined group all timeblocks
AMOVAon$varcomp$P.value[1] #p value ON combined group all timeblocks
AMOVAstl$varcomp$sigma2[1]/sum(AMOVAstl$varcomp$sigma2) #Phi STL combined group all timeblocks
AMOVAstl$varcomp$P.value[1] #p value STL combined group all timeblocks

AMOVAon.tb1$varcomp$sigma2[1]/sum(AMOVAon.tb1$varcomp$sigma2) #Phi ON combined group timeblock 1
AMOVAon.tb1$varcomp$P.value[1] #p value ON combined group timeblock 1
AMOVAon.tb2$varcomp$sigma2[1]/sum(AMOVAon.tb2$varcomp$sigma2) #Phi ON combined group timeblock 2
AMOVAon.tb2$varcomp$P.value[1] #p value ON combined group timeblock 2
AMOVAon.tb3$varcomp$sigma2[1]/sum(AMOVAon.tb3$varcomp$sigma2) #Phi ON combined group timeblock 3
AMOVAon.tb3$varcomp$P.value[1] #p value ON combined group timeblock 3
AMOVAon.tb4$varcomp$sigma2[1]/sum(AMOVAon.tb4$varcomp$sigma2) #Phi ON combined group timeblock 4
AMOVAon.tb4$varcomp$P.value[1] #p value ON combined group timeblock 4
AMOVAon.tb5$varcomp$sigma2[1]/sum(AMOVAon.tb5$varcomp$sigma2) #Phi ON combined group timeblock 5
AMOVAon.tb5$varcomp$P.value[1] #p value ON combined group timeblock 5

AMOVAny.tb1$varcomp$sigma2[1]/sum(AMOVAny.tb1$varcomp$sigma2) #Phi NY combined group timeblock 1
AMOVAny.tb1$varcomp$P.value[1] #p value NY combined group timeblock 1
AMOVAny.tb2$varcomp$sigma2[1]/sum(AMOVAny.tb2$varcomp$sigma2) #Phi NY combined group timeblock 2
AMOVAny.tb2$varcomp$P.value[1] #p value NY combined group timeblock 2
AMOVAny.tb3$varcomp$sigma2[1]/sum(AMOVAny.tb3$varcomp$sigma2) #Phi NY combined group timeblock 3
AMOVAny.tb3$varcomp$P.value[1] #p value NY combined group timeblock 3
AMOVAny.tb4$varcomp$sigma2[1]/sum(AMOVAny.tb4$varcomp$sigma2) #Phi NY combined group timeblock 4
AMOVAny.tb4$varcomp$P.value[1] #p value NY combined group timeblock 4
AMOVAny.tb5$varcomp$sigma2[1]/sum(AMOVAny.tb5$varcomp$sigma2) #Phi NY combined group timeblock 5
AMOVAny.tb5$varcomp$P.value[1] #p value NY combined group timeblock 5

AMOVAstl.tb1$varcomp$sigma2[1]/sum(AMOVAstl.tb1$varcomp$sigma2) #Phi STL combined group timeblock 1
AMOVAstl.tb1$varcomp$P.value[1] #p value STL combined group timeblock 1
AMOVAstl.tb2$varcomp$sigma2[1]/sum(AMOVAstl.tb2$varcomp$sigma2) #Phi STL combined group timeblock 2
AMOVAstl.tb2$varcomp$P.value[1] #p value STL combined group timeblock 2
AMOVAstl.tb3$varcomp$sigma2[1]/sum(AMOVAstl.tb3$varcomp$sigma2) #Phi STL combined group timeblock 3
AMOVAstl.tb3$varcomp$P.value[1] #p value STL combined group timeblock 3

### Pairwise AMOVA Tests ###

### Define Group Names ###
spgrnames_global<-unique(as.character(groups_global))
spgrnames_global.cgroups<-unique(as.character(cgroups_global)) #regional group names
spgrnames_tb1<-unique(as.character(groups_tb1))
spgrnames_tb2<-unique(as.character(groups_tb2))
spgrnames_tb3<-unique(as.character(groups_tb3))
spgrnames_tb4<-unique(as.character(groups_tb4))
spgrnames_tb5<-unique(as.character(groups_tb5))
spgrnames_tb1.cgroups<-unique(as.character(cgroups_tb1))
spgrnames_tb2.cgroups<-unique(as.character(cgroups_tb2))
spgrnames_tb3.cgroups<-unique(as.character(cgroups_tb3))
spgrnames_tb4.cgroups<-unique(as.character(cgroups_tb4))
spgrnames_tb5.cgroups<-unique(as.character(cgroups_tb5))
spgrnames_stl<-unique(as.character(groups_stl))
spgrnames_ny<-unique(as.character(groups_ny))
spgrnames_on<-unique(as.character(groups_on))
spgrnames_ny.tb1<-unique(as.character(groups_ny.tb1))
spgrnames_ny.tb2<-unique(as.character(groups_ny.tb2))
spgrnames_ny.tb3<-unique(as.character(groups_ny.tb3))
spgrnames_ny.tb4<-unique(as.character(groups_ny.tb4))
spgrnames_ny.tb5<-unique(as.character(groups_ny.tb5))
spgrnames_on.tb1<-unique(as.character(groups_on.tb1))
spgrnames_on.tb2<-unique(as.character(groups_on.tb2))
spgrnames_on.tb3<-unique(as.character(groups_on.tb3))
spgrnames_on.tb4<-unique(as.character(groups_on.tb4))
spgrnames_on.tb5<-unique(as.character(groups_on.tb5))
spgrnames_stl.tb1<-unique(as.character(groups_stl.tb1))
spgrnames_stl.tb2<-unique(as.character(groups_stl.tb2))
spgrnames_stl.tb3<-unique(as.character(groups_stl.tb3))

### Define Number of Groups ###
numspgr_global<-length(spgrnames_global)
numspgr_global.cgroups<-length(spgrnames_global.cgroups) #regional groups
numspgr_tb1<-length(spgrnames_tb1)
numspgr_tb2<-length(spgrnames_tb2)
numspgr_tb3<-length(spgrnames_tb3)
numspgr_tb4<-length(spgrnames_tb4)
numspgr_tb5<-length(spgrnames_tb5)  
numspgr_tb1.cgroups<-length(spgrnames_tb1.cgroups)
numspgr_tb2.cgroups<-length(spgrnames_tb2.cgroups)
numspgr_tb3.cgroups<-length(spgrnames_tb3.cgroups)
numspgr_tb4.cgroups<-length(spgrnames_tb4.cgroups)
numspgr_tb5.cgroups<-length(spgrnames_tb5.cgroups)
numspgr_stl<-length(spgrnames_stl)
numspgr_ny<-length(spgrnames_ny)
numspgr_on<-length(spgrnames_on)
numspgr_ny.tb1<-length(spgrnames_ny.tb1)
numspgr_ny.tb2<-length(spgrnames_ny.tb2)
numspgr_ny.tb3<-length(spgrnames_ny.tb3)
numspgr_ny.tb4<-length(spgrnames_ny.tb4)
numspgr_ny.tb5<-length(spgrnames_ny.tb5)
numspgr_on.tb1<-length(spgrnames_on.tb1)
numspgr_on.tb2<-length(spgrnames_on.tb2)
numspgr_on.tb3<-length(spgrnames_on.tb3)
numspgr_on.tb4<-length(spgrnames_on.tb4)
numspgr_on.tb5<-length(spgrnames_on.tb5)
numspgr_stl.tb1<-length(spgrnames_stl.tb1)
numspgr_stl.tb2<-length(spgrnames_stl.tb2)
numspgr_stl.tb3<-length(spgrnames_stl.tb3)

###Create Symmetrical Square Matrix ###
sqbrdist_global.matrix <- as.matrix(sqbrdist_global)
sqbrdist_tb1.matrix <- as.matrix(timeblock1.sqbrdist)
sqbrdist_tb2.matrix <- as.matrix(timeblock2.sqbrdist)
sqbrdist_tb3.matrix <- as.matrix(timeblock3.sqbrdist)
sqbrdist_tb4.matrix <- as.matrix(timeblock4.sqbrdist)
sqbrdist_tb5.matrix <- as.matrix(timeblock5.sqbrdist)

stl_sqbrdist.matrix <- as.matrix(stl_sqbrdist)
ny_sqbrdist.matrix <- as.matrix(ny_sqbrdist)
on_sqbrdist.matrix <- as.matrix(on_sqbrdist)

ny_sqbrdist.tb1.matrix <- as.matrix(ny_sqbrdist.tb1)
ny_sqbrdist.tb2.matrix <- as.matrix(ny_sqbrdist.tb2)
ny_sqbrdist.tb3.matrix <- as.matrix(ny_sqbrdist.tb3)
ny_sqbrdist.tb4.matrix <- as.matrix(ny_sqbrdist.tb4)
ny_sqbrdist.tb5.matrix <- as.matrix(ny_sqbrdist.tb5)

on_sqbrdist.tb1.matrix <- as.matrix(on_sqbrdist.tb1)
on_sqbrdist.tb2.matrix <- as.matrix(on_sqbrdist.tb2)
on_sqbrdist.tb3.matrix <- as.matrix(on_sqbrdist.tb3)
on_sqbrdist.tb4.matrix <- as.matrix(on_sqbrdist.tb4)
on_sqbrdist.tb5.matrix <- as.matrix(on_sqbrdist.tb5)

stl_sqbrdist.tb1.matrix <- as.matrix(stl_sqbrdist.tb1)
stl_sqbrdist.tb2.matrix <- as.matrix(stl_sqbrdist.tb2)
stl_sqbrdist.tb3.matrix <- as.matrix(stl_sqbrdist.tb3)

### Define list of Groups for Looping ###
spatiallist_global <- groups_global
spatiallist_global.cgroups <- cgroups_global

spatiallist_tb1 <- groups_tb1
spatiallist_tb2 <- groups_tb2
spatiallist_tb3 <- groups_tb3
spatiallist_tb4 <- groups_tb4
spatiallist_tb5 <- groups_tb5 

spatiallist_tb1.cgroups <- cgroups_tb1
spatiallist_tb2.cgroups <- cgroups_tb2
spatiallist_tb3.cgroups <- cgroups_tb3
spatiallist_tb4.cgroups <- cgroups_tb4
spatiallist_tb5.cgroups <- cgroups_tb5

spatiallist_stl <- groups_stl
spatiallist_ny <- groups_ny
spatiallist_on <- groups_on

spatiallist_ny.tb1 <- groups_ny.tb1
spatiallist_ny.tb2 <- groups_ny.tb2
spatiallist_ny.tb3 <- groups_ny.tb3
spatiallist_ny.tb4 <- groups_ny.tb4
spatiallist_ny.tb5 <- groups_ny.tb5

spatiallist_on.tb1 <- groups_on.tb1
spatiallist_on.tb2 <- groups_on.tb2
spatiallist_on.tb3 <- groups_on.tb3
spatiallist_on.tb4 <- groups_on.tb4
spatiallist_on.tb5 <- groups_on.tb5

spatiallist_stl.tb1 <- groups_stl.tb1
spatiallist_stl.tb2 <- groups_stl.tb2
spatiallist_stl.tb3 <- groups_stl.tb3

###Define PhiST and pval output matrices ###
PhiMatrix_global<-matrix(NA,nrow=numspgr_global,ncol=numspgr_global)
PhiMatrix_global.cgroups<-matrix(NA,nrow=numspgr_global.cgroups,ncol=numspgr_global.cgroups)

PhiMatrix_tb1<-matrix(NA,nrow=numspgr_tb1,ncol=numspgr_tb1)
PhiMatrix_tb2<-matrix(NA,nrow=numspgr_tb2,ncol=numspgr_tb2)
PhiMatrix_tb3<-matrix(NA,nrow=numspgr_tb3,ncol=numspgr_tb3)
PhiMatrix_tb4<-matrix(NA,nrow=numspgr_tb4,ncol=numspgr_tb4)
PhiMatrix_tb5<-matrix(NA,nrow=numspgr_tb5,ncol=numspgr_tb5)

PhiMatrix_tb1.cgroups<-matrix(NA,nrow=numspgr_tb1.cgroups,ncol=numspgr_tb1.cgroups)
PhiMatrix_tb2.cgroups<-matrix(NA,nrow=numspgr_tb2.cgroups,ncol=numspgr_tb2.cgroups)
PhiMatrix_tb3.cgroups<-matrix(NA,nrow=numspgr_tb3.cgroups,ncol=numspgr_tb3.cgroups)
PhiMatrix_tb4.cgroups<-matrix(NA,nrow=numspgr_tb4.cgroups,ncol=numspgr_tb4.cgroups)
PhiMatrix_tb5.cgroups<-matrix(NA,nrow=numspgr_tb5.cgroups,ncol=numspgr_tb5.cgroups)

PhiMatrix_stl<-matrix(NA,nrow=numspgr_stl,ncol=numspgr_stl)
PhiMatrix_ny<-matrix(NA,nrow=numspgr_ny,ncol=numspgr_ny)
PhiMatrix_on<-matrix(NA,nrow=numspgr_on,ncol=numspgr_on)

PhiMatrix_ny.tb1<-matrix(NA,nrow=numspgr_ny.tb1,ncol=numspgr_ny.tb1)
PhiMatrix_ny.tb2<-matrix(NA,nrow=numspgr_ny.tb2,ncol=numspgr_ny.tb2)
PhiMatrix_ny.tb3<-matrix(NA,nrow=numspgr_ny.tb3,ncol=numspgr_ny.tb3)
PhiMatrix_ny.tb4<-matrix(NA,nrow=numspgr_ny.tb4,ncol=numspgr_ny.tb4)
PhiMatrix_ny.tb5<-matrix(NA,nrow=numspgr_ny.tb5,ncol=numspgr_ny.tb5)

PhiMatrix_on.tb1<-matrix(NA,nrow=numspgr_on.tb1,ncol=numspgr_on.tb1)
PhiMatrix_on.tb2<-matrix(NA,nrow=numspgr_on.tb2,ncol=numspgr_on.tb2)
PhiMatrix_on.tb3<-matrix(NA,nrow=numspgr_on.tb3,ncol=numspgr_on.tb3)
PhiMatrix_on.tb4<-matrix(NA,nrow=numspgr_on.tb4,ncol=numspgr_on.tb4)
PhiMatrix_on.tb5<-matrix(NA,nrow=numspgr_on.tb5,ncol=numspgr_on.tb5)

PhiMatrix_stl.tb1<-matrix(NA,nrow=numspgr_stl.tb1,ncol=numspgr_stl.tb1)
PhiMatrix_stl.tb2<-matrix(NA,nrow=numspgr_stl.tb2,ncol=numspgr_stl.tb2)
PhiMatrix_stl.tb3<-matrix(NA,nrow=numspgr_stl.tb3,ncol=numspgr_stl.tb3)

PvalMatrix_global<-matrix(NA,nrow=numspgr_global,ncol=numspgr_global)
PvalMatrix_global.cgroups<-matrix(NA,nrow=numspgr_global.cgroups,ncol=numspgr_global.cgroups)

PvalMatrix_tb1.cgroups<-matrix(NA,nrow=numspgr_tb1.cgroups,ncol=numspgr_tb1.cgroups)
PvalMatrix_tb2.cgroups<-matrix(NA,nrow=numspgr_tb2.cgroups,ncol=numspgr_tb2.cgroups)
PvalMatrix_tb3.cgroups<-matrix(NA,nrow=numspgr_tb3.cgroups,ncol=numspgr_tb3.cgroups)
PvalMatrix_tb4.cgroups<-matrix(NA,nrow=numspgr_tb4.cgroups,ncol=numspgr_tb4.cgroups)
PvalMatrix_tb5.cgroups<-matrix(NA,nrow=numspgr_tb5.cgroups,ncol=numspgr_tb5.cgroups)  

PvalMatrix_tb1<-matrix(NA,nrow=numspgr_tb1,ncol=numspgr_tb1)
PvalMatrix_tb2<-matrix(NA,nrow=numspgr_tb2,ncol=numspgr_tb2)
PvalMatrix_tb3<-matrix(NA,nrow=numspgr_tb3,ncol=numspgr_tb3)
PvalMatrix_tb4<-matrix(NA,nrow=numspgr_tb4,ncol=numspgr_tb4)
PvalMatrix_tb5<-matrix(NA,nrow=numspgr_tb5,ncol=numspgr_tb5)

PvalMatrix_stl<-matrix(NA,nrow=numspgr_stl,ncol=numspgr_stl)
PvalMatrix_ny<-matrix(NA,nrow=numspgr_ny,ncol=numspgr_ny)
PvalMatrix_on<-matrix(NA,nrow=numspgr_on,ncol=numspgr_on)

PvalMatrix_ny.tb1<-matrix(NA,nrow=numspgr_ny.tb1,ncol=numspgr_ny.tb1)
PvalMatrix_ny.tb2<-matrix(NA,nrow=numspgr_ny.tb2,ncol=numspgr_ny.tb2)
PvalMatrix_ny.tb3<-matrix(NA,nrow=numspgr_ny.tb3,ncol=numspgr_ny.tb3)
PvalMatrix_ny.tb4<-matrix(NA,nrow=numspgr_ny.tb4,ncol=numspgr_ny.tb4)
PvalMatrix_ny.tb5<-matrix(NA,nrow=numspgr_ny.tb5,ncol=numspgr_ny.tb5)

PvalMatrix_on.tb1<-matrix(NA,nrow=numspgr_on.tb1,ncol=numspgr_on.tb1)
PvalMatrix_on.tb2<-matrix(NA,nrow=numspgr_on.tb2,ncol=numspgr_on.tb2)
PvalMatrix_on.tb3<-matrix(NA,nrow=numspgr_on.tb3,ncol=numspgr_on.tb3)
PvalMatrix_on.tb4<-matrix(NA,nrow=numspgr_on.tb4,ncol=numspgr_on.tb4)
PvalMatrix_on.tb5<-matrix(NA,nrow=numspgr_on.tb5,ncol=numspgr_on.tb5)

PvalMatrix_stl.tb1<-matrix(NA,nrow=numspgr_stl.tb1,ncol=numspgr_stl.tb1)
PvalMatrix_stl.tb2<-matrix(NA,nrow=numspgr_stl.tb2,ncol=numspgr_stl.tb2)
PvalMatrix_stl.tb3<-matrix(NA,nrow=numspgr_stl.tb3,ncol=numspgr_stl.tb3)

### Begin For Looping###
### All For Loops Identical except for groups and data ###

##All Sites
XYT_global<-data.frame(x=numeric(length=numspgr_global),y=numeric(length=numspgr_global),t=numeric(length=numspgr_global)) #define spatial group centroids input dataframe based on length of number of groups
for (i in 1:numspgr_global) #loop over each group in the case scenario
{
  print(paste(i," out of ",numspgr_global))
  for (j in 1:numspgr_global)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_global)==as.character(spgrnames_global[i])) #select first group to pair
      indexJ=which(as.character(spatiallist_global)==as.character(spgrnames_global[j])) #select second group to pair
      #Define subset of all sites in each group
      subDistMat<-sqbrdist_global.matrix[c(indexI,indexJ),c(indexI,indexJ)] #fill subset distance matrix with all BR distances from all pairs of sites in the two groups
      subDistMat<-as.dist(subDistMat) #change class back to distance matrix
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ)))) #character matrix defining whether each site is grouped into group 1 or 2
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000) #compute differentiation between group 1 and 2, permute to assess significance
      N0=tmp$varcoef #variance coefficient of paired groups
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef) #BR distance variation constrained between groups -- mean square deviation between groups - mean square deviation across groups / variance coefficient
      WP=(tmp$tab$MSD[2]) #BR distance variation within groups
      PhiMatrix_global[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)) #sum of variation constrained between group 1 and 2 -- greater differentiation of either pair makes this value larger
      if(PhiMatrix_global[i,j]<0){PhiMatrix_global[i,j]=0} #end loop once all pairs are assessed
      PvalMatrix_global[i,j]=tmp$varcomp$P.value[1] #add p value to matrix
    }
  }
  XYT_global$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_global[i]))],na.rm=TRUE) #calculate mean longitude of group centroid
  XYT_global$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_global[i]))],na.rm=TRUE) #calculate mean latitude of group centroid
}
row.names(PhiMatrix_global)<-spgrnames_global #set row and column names of matrix equal to names of paired spatial groups assessed
colnames(PhiMatrix_global)<-spgrnames_global
row.names(PvalMatrix_global)<-spgrnames_global
colnames(PvalMatrix_global)<-spgrnames_global

#Timeblock I
XYT_tb1<-data.frame(x=numeric(length=numspgr_tb1),y=numeric(length=numspgr_tb1),t=numeric(length=numspgr_tb1)) #different length
for (i in 1:numspgr_tb1)
{
  print(paste(i," out of ",numspgr_tb1))
  for (j in 1:numspgr_tb1)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_tb1)==as.character(spgrnames_tb1[i])) #different list of groups
      indexJ=which(as.character(spatiallist_tb1)==as.character(spgrnames_tb1[j]))
      #subsetDistance Matrix
      subDistMat<-sqbrdist_tb1.matrix[c(indexI,indexJ),c(indexI,indexJ)] #new square matrix
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_tb1[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_tb1[i,j]<0){PhiMatrix_tb1[i,j]=0} #different matrix output
      PvalMatrix_tb1[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_tb1$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_tb1[i]))],na.rm=TRUE)
  XYT_tb1$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_tb1[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_tb1)<-spgrnames_tb1
colnames(PhiMatrix_tb1)<-spgrnames_tb1
row.names(PvalMatrix_tb1)<-spgrnames_tb1
colnames(PvalMatrix_tb1)<-spgrnames_tb1

#Timeblock 2
XYT_tb2<-data.frame(x=numeric(length=numspgr_tb2),y=numeric(length=numspgr_tb2),t=numeric(length=numspgr_tb2))
for (i in 1:numspgr_tb2)
{
  print(paste(i," out of ",numspgr_tb2))
  for (j in 1:numspgr_tb2)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_tb2)==as.character(spgrnames_tb2[i]))
      indexJ=which(as.character(spatiallist_tb2)==as.character(spgrnames_tb2[j]))
      #subsetDistance Matrix
      subDistMat<-sqbrdist_tb2.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_tb2[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_tb2[i,j]<0){PhiMatrix_tb2[i,j]=0}
      PvalMatrix_tb2[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_tb2$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_tb2[i]))],na.rm=TRUE)
  XYT_tb2$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_tb2[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_tb2)<-spgrnames_tb2
colnames(PhiMatrix_tb2)<-spgrnames_tb2
row.names(PvalMatrix_tb2)<-spgrnames_tb2
colnames(PvalMatrix_tb2)<-spgrnames_tb2

#Timeblock 3
XYT_tb3<-data.frame(x=numeric(length=numspgr_tb3),y=numeric(length=numspgr_tb3),t=numeric(length=numspgr_tb3))
for (i in 1:numspgr_tb3)
{
  print(paste(i," out of ",numspgr_tb3))
  for (j in 1:numspgr_tb3)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_tb3)==as.character(spgrnames_tb3[i]))
      indexJ=which(as.character(spatiallist_tb3)==as.character(spgrnames_tb3[j]))
      #subsetDistance Matrix
      subDistMat<-sqbrdist_tb3.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_tb3[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_tb3[i,j]<0){PhiMatrix_tb3[i,j]=0}
      PvalMatrix_tb3[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_tb3$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_tb3[i]))],na.rm=TRUE)
  XYT_tb3$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_tb3[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_tb3)<-spgrnames_tb3
colnames(PhiMatrix_tb3)<-spgrnames_tb3
row.names(PvalMatrix_tb3)<-spgrnames_tb3
colnames(PvalMatrix_tb3)<-spgrnames_tb3

#Timeblock 4
XYT_tb4<-data.frame(x=numeric(length=numspgr_tb4),y=numeric(length=numspgr_tb4),t=numeric(length=numspgr_tb4))
for (i in 1:numspgr_tb4)
{
  print(paste(i," out of ",numspgr_tb4))
  for (j in 1:numspgr_tb4)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_tb4)==as.character(spgrnames_tb4[i]))
      indexJ=which(as.character(spatiallist_tb4)==as.character(spgrnames_tb4[j]))
      #subsetDistance Matrix
      subDistMat<-sqbrdist_tb4.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_tb4[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_tb4[i,j]<0){PhiMatrix_tb4[i,j]=0}
      PvalMatrix_tb4[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_tb4$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_tb4[i]))],na.rm=TRUE)
  XYT_tb4$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_tb4[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_tb4)<-spgrnames_tb4
colnames(PhiMatrix_tb4)<-spgrnames_tb4
row.names(PvalMatrix_tb4)<-spgrnames_tb4
colnames(PvalMatrix_tb4)<-spgrnames_tb4

#Timeblock 5
XYT_tb5<-data.frame(x=numeric(length=numspgr_tb5),y=numeric(length=numspgr_tb5),t=numeric(length=numspgr_tb5))
for (i in 1:numspgr_tb5)
{
  print(paste(i," out of ",numspgr_tb5))
  for (j in 1:numspgr_tb5)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_tb5)==as.character(spgrnames_tb5[i]))
      indexJ=which(as.character(spatiallist_tb5)==as.character(spgrnames_tb5[j]))
      #subsetDistance Matrix
      subDistMat<-sqbrdist_tb5.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_tb5[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_tb5[i,j]<0){PhiMatrix_tb5[i,j]=0}
      PvalMatrix_tb5[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_tb5$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_tb5[i]))],na.rm=TRUE)
  XYT_tb5$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_tb5[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_tb5)<-spgrnames_tb5
colnames(PhiMatrix_tb5)<-spgrnames_tb5
row.names(PvalMatrix_tb5)<-spgrnames_tb5
colnames(PvalMatrix_tb5)<-spgrnames_tb5

#All sites grouped by region
XYT_global.cgroups<-data.frame(x=numeric(length=numspgr_global.cgroups),y=numeric(length=numspgr_global.cgroups),t=numeric(length=numspgr_global.cgroups))
for (i in 1:numspgr_global.cgroups)
{
  print(paste(i," out of ",numspgr_global.cgroups))
  for (j in 1:numspgr_global.cgroups)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_global.cgroups)==as.character(spgrnames_global.cgroups[i]))
      indexJ=which(as.character(spatiallist_global.cgroups)==as.character(spgrnames_global.cgroups[j]))
      #subsetDistance Matrix
      subDistMat<-sqbrdist_global.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_global.cgroups[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_global.cgroups[i,j]<0){PhiMatrix_global.cgroups[i,j]=0}
      PvalMatrix_global.cgroups[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_global.cgroups$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_global.cgroups[i]))],na.rm=TRUE)
  XYT_global.cgroups$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_global.cgroups[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_global.cgroups)<-spgrnames_global.cgroups
colnames(PhiMatrix_global.cgroups)<-spgrnames_global.cgroups
row.names(PvalMatrix_global.cgroups)<-spgrnames_global.cgroups
colnames(PvalMatrix_global.cgroups)<-spgrnames_global.cgroups

#timeblock 1 region

XYT_tb1.cgroups<-data.frame(x=numeric(length=numspgr_tb1.cgroups),y=numeric(length=numspgr_tb1.cgroups),t=numeric(length=numspgr_tb1.cgroups))
for (i in 1:numspgr_tb1.cgroups)
{
  print(paste(i," out of ",numspgr_tb1.cgroups))
  for (j in 1:numspgr_tb1.cgroups)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_tb1.cgroups)==as.character(spgrnames_tb1.cgroups[i]))
      indexJ=which(as.character(spatiallist_tb1.cgroups)==as.character(spgrnames_tb1.cgroups[j]))
      #subsetDistance Matrix
      subDistMat<-sqbrdist_tb1.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_tb1.cgroups[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_tb1.cgroups[i,j]<0){PhiMatrix_tb1.cgroups[i,j]=0}
      PvalMatrix_tb1.cgroups[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_tb1.cgroups$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_tb1.cgroups[i]))],na.rm=TRUE)
  XYT_tb1.cgroups$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_tb1.cgroups[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_tb1.cgroups)<-spgrnames_tb1.cgroups
colnames(PhiMatrix_tb1.cgroups)<-spgrnames_tb1.cgroups
row.names(PvalMatrix_tb1.cgroups)<-spgrnames_tb1.cgroups
colnames(PvalMatrix_tb1.cgroups)<-spgrnames_tb1.cgroups

#Timeblock 2 region
XYT_tb2.cgroups<-data.frame(x=numeric(length=numspgr_tb2.cgroups),y=numeric(length=numspgr_tb2.cgroups),t=numeric(length=numspgr_tb2.cgroups))
for (i in 1:numspgr_tb2.cgroups)
{
  print(paste(i," out of ",numspgr_tb2.cgroups))
  for (j in 1:numspgr_tb2.cgroups)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_tb2.cgroups)==as.character(spgrnames_tb2.cgroups[i]))
      indexJ=which(as.character(spatiallist_tb2.cgroups)==as.character(spgrnames_tb2.cgroups[j]))
      #subsetDistance Matrix
      subDistMat<-sqbrdist_tb2.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_tb2.cgroups[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_tb2.cgroups[i,j]<0){PhiMatrix_tb2.cgroups[i,j]=0}
      PvalMatrix_tb2.cgroups[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_tb2.cgroups$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_tb2.cgroups[i]))],na.rm=TRUE)
  XYT_tb2.cgroups$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_tb2.cgroups[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_tb2.cgroups)<-spgrnames_tb2.cgroups
colnames(PhiMatrix_tb2.cgroups)<-spgrnames_tb2.cgroups
row.names(PvalMatrix_tb2.cgroups)<-spgrnames_tb2.cgroups
colnames(PvalMatrix_tb2.cgroups)<-spgrnames_tb2.cgroups

#Timeblock 3 region
XYT_tb3.cgroups<-data.frame(x=numeric(length=numspgr_tb3.cgroups),y=numeric(length=numspgr_tb3.cgroups),t=numeric(length=numspgr_tb3.cgroups))
for (i in 1:numspgr_tb3.cgroups)
{
  print(paste(i," out of ",numspgr_tb3.cgroups))
  for (j in 1:numspgr_tb3.cgroups)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_tb3.cgroups)==as.character(spgrnames_tb3.cgroups[i]))
      indexJ=which(as.character(spatiallist_tb3.cgroups)==as.character(spgrnames_tb3.cgroups[j]))
      #subsetDistance Matrix
      subDistMat<-sqbrdist_tb3.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_tb3.cgroups[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_tb3.cgroups[i,j]<0){PhiMatrix_tb3.cgroups[i,j]=0}
      PvalMatrix_tb3.cgroups[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_tb3.cgroups$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_tb3.cgroups[i]))],na.rm=TRUE)
  XYT_tb3.cgroups$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_tb3.cgroups[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_tb3.cgroups)<-spgrnames_tb3.cgroups
colnames(PhiMatrix_tb3.cgroups)<-spgrnames_tb3.cgroups
row.names(PvalMatrix_tb3.cgroups)<-spgrnames_tb3.cgroups
colnames(PvalMatrix_tb3.cgroups)<-spgrnames_tb3.cgroups

#Timeblock 4 region
XYT_tb4.cgroups<-data.frame(x=numeric(length=numspgr_tb4.cgroups),y=numeric(length=numspgr_tb4.cgroups),t=numeric(length=numspgr_tb4.cgroups))
for (i in 1:numspgr_tb4.cgroups)
{
  print(paste(i," out of ",numspgr_tb4.cgroups))
  for (j in 1:numspgr_tb4.cgroups)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_tb4.cgroups)==as.character(spgrnames_tb4.cgroups[i]))
      indexJ=which(as.character(spatiallist_tb4.cgroups)==as.character(spgrnames_tb4.cgroups[j]))
      #subsetDistance Matrix
      subDistMat<-sqbrdist_tb4.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_tb4.cgroups[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_tb4.cgroups[i,j]<0){PhiMatrix_tb4.cgroups[i,j]=0}
      PvalMatrix_tb4.cgroups[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_tb4.cgroups$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_tb4.cgroups[i]))],na.rm=TRUE)
  XYT_tb4.cgroups$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_tb4.cgroups[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_tb4.cgroups)<-spgrnames_tb4.cgroups
colnames(PhiMatrix_tb4.cgroups)<-spgrnames_tb4.cgroups
row.names(PvalMatrix_tb4.cgroups)<-spgrnames_tb4.cgroups
colnames(PvalMatrix_tb4.cgroups)<-spgrnames_tb4.cgroups

#Timeblock 5 region
XYT_tb5.cgroups<-data.frame(x=numeric(length=numspgr_tb5.cgroups),y=numeric(length=numspgr_tb5.cgroups),t=numeric(length=numspgr_tb5.cgroups))
for (i in 1:numspgr_tb5.cgroups)
{
  print(paste(i," out of ",numspgr_tb5.cgroups))
  for (j in 1:numspgr_tb5.cgroups)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_tb5.cgroups)==as.character(spgrnames_tb5.cgroups[i]))
      indexJ=which(as.character(spatiallist_tb5.cgroups)==as.character(spgrnames_tb5.cgroups[j]))
      #subsetDistance Matrix
      subDistMat<-sqbrdist_tb5.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_tb5.cgroups[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_tb5.cgroups[i,j]<0){PhiMatrix_tb5.cgroups[i,j]=0}
      PvalMatrix_tb5.cgroups[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_tb5.cgroups$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_tb5.cgroups[i]))],na.rm=TRUE)
  XYT_tb5.cgroups$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_tb5.cgroups[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_tb5.cgroups)<-spgrnames_tb5.cgroups
colnames(PhiMatrix_tb5.cgroups)<-spgrnames_tb5.cgroups
row.names(PvalMatrix_tb5.cgroups)<-spgrnames_tb5.cgroups
colnames(PvalMatrix_tb5.cgroups)<-spgrnames_tb5.cgroups


#New York All Sites
XYT_ny<-data.frame(x=numeric(length=numspgr_ny),y=numeric(length=numspgr_ny),t=numeric(length=numspgr_ny))
for (i in 1:numspgr_ny)
{
  print(paste(i," out of ",numspgr_ny))
  for (j in 1:numspgr_ny)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_ny)==as.character(spgrnames_ny[i]))
      indexJ=which(as.character(spatiallist_ny)==as.character(spgrnames_ny[j]))
      #subsetDistance Matrix
      subDistMat<-ny_sqbrdist.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_ny[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_ny[i,j]<0){PhiMatrix_ny[i,j]=0}
      PvalMatrix_ny[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_ny$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_ny[i]))],na.rm=TRUE)
  XYT_ny$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_ny[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_ny)<-spgrnames_ny
colnames(PhiMatrix_ny)<-spgrnames_ny
row.names(PvalMatrix_ny)<-spgrnames_ny
colnames(PvalMatrix_ny)<-spgrnames_ny

#Ontario all sites
XYT_on<-data.frame(x=numeric(length=numspgr_on),y=numeric(length=numspgr_on),t=numeric(length=numspgr_on))
for (i in 1:numspgr_on)
{
  print(paste(i," out of ",numspgr_on))
  for (j in 1:numspgr_on)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_on)==as.character(spgrnames_on[i]))
      indexJ=which(as.character(spatiallist_on)==as.character(spgrnames_on[j]))
      #subsetDistance Matrix
      subDistMat<-on_sqbrdist.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_on[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_on[i,j]<0){PhiMatrix_on[i,j]=0}
      PvalMatrix_on[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_on$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_on[i]))],na.rm=TRUE)
  XYT_on$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_on[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_on)<-spgrnames_on
colnames(PhiMatrix_on)<-spgrnames_on
row.names(PvalMatrix_on)<-spgrnames_on
colnames(PvalMatrix_on)<-spgrnames_on

#New York Timeblock 1
XYT_ny.tb1<-data.frame(x=numeric(length=numspgr_ny.tb1),y=numeric(length=numspgr_ny.tb1),t=numeric(length=numspgr_ny.tb1))
for (i in 1:numspgr_ny.tb1)
{
  print(paste(i," out of ",numspgr_ny.tb1))
  for (j in 1:numspgr_ny.tb1)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_ny.tb1)==as.character(spgrnames_ny.tb1[i]))
      indexJ=which(as.character(spatiallist_ny.tb1)==as.character(spgrnames_ny.tb1[j]))
      #subsetDistance Matrix
      subDistMat<-ny_sqbrdist.tb1.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_ny.tb1[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_ny.tb1[i,j]<0){PhiMatrix_ny.tb1[i,j]=0}
      PvalMatrix_ny.tb1[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_ny.tb1$x[i]<-mean(ny.dt$lon_dd[which(dt$groupcode==as.character(spgrnames_ny.tb1[i]))],na.rm=TRUE)
  XYT_ny.tb1$y[i]<-mean(ny.dt$lat_dd[which(dt$groupcode==as.character(spgrnames_ny.tb1[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_ny.tb1)<-spgrnames_ny.tb1
colnames(PhiMatrix_ny.tb1)<-spgrnames_ny.tb1
row.names(PvalMatrix_ny.tb1)<-spgrnames_ny.tb1
colnames(PvalMatrix_ny.tb1)<-spgrnames_ny.tb1

#NY Timeblock 2
XYT_ny.tb2<-data.frame(x=numeric(length=numspgr_ny.tb2),y=numeric(length=numspgr_ny.tb2),t=numeric(length=numspgr_ny.tb2))
for (i in 1:numspgr_ny.tb2)
{
  print(paste(i," out of ",numspgr_ny.tb2))
  for (j in 1:numspgr_ny.tb2)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_ny.tb2)==as.character(spgrnames_ny.tb2[i]))
      indexJ=which(as.character(spatiallist_ny.tb2)==as.character(spgrnames_ny.tb2[j]))
      #subsetDistance Matrix
      subDistMat<-ny_sqbrdist.tb2.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_ny.tb2[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_ny.tb2[i,j]<0){PhiMatrix_ny.tb2[i,j]=0}
      PvalMatrix_ny.tb2[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_ny.tb2$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_ny.tb2[i]))],na.rm=TRUE)
  XYT_ny.tb2$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_ny.tb2[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_ny.tb2)<-spgrnames_ny.tb2
colnames(PhiMatrix_ny.tb2)<-spgrnames_ny.tb2
row.names(PvalMatrix_ny.tb2)<-spgrnames_ny.tb2
colnames(PvalMatrix_ny.tb2)<-spgrnames_ny.tb2

#NY Timeblock 3
XYT_ny.tb3<-data.frame(x=numeric(length=numspgr_ny.tb3),y=numeric(length=numspgr_ny.tb3),t=numeric(length=numspgr_ny.tb3))
for (i in 1:numspgr_ny.tb3)
{
  print(paste(i," out of ",numspgr_ny.tb3))
  for (j in 1:numspgr_ny.tb3)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_ny.tb3)==as.character(spgrnames_ny.tb3[i]))
      indexJ=which(as.character(spatiallist_ny.tb3)==as.character(spgrnames_ny.tb3[j]))
      #subsetDistance Matrix
      subDistMat<-ny_sqbrdist.tb3.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_ny.tb3[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_ny.tb3[i,j]<0){PhiMatrix_ny.tb3[i,j]=0}
      PvalMatrix_ny.tb3[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_ny.tb3$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_ny.tb3[i]))],na.rm=TRUE)
  XYT_ny.tb3$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_ny.tb3[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_ny.tb3)<-spgrnames_ny.tb3
colnames(PhiMatrix_ny.tb3)<-spgrnames_ny.tb3
row.names(PvalMatrix_ny.tb3)<-spgrnames_ny.tb3
colnames(PvalMatrix_ny.tb3)<-spgrnames_ny.tb3
#NY Timeblock 4
XYT_ny.tb4<-data.frame(x=numeric(length=numspgr_ny.tb4),y=numeric(length=numspgr_ny.tb4),t=numeric(length=numspgr_ny.tb4))
for (i in 1:numspgr_ny.tb4)
{
  print(paste(i," out of ",numspgr_ny.tb4))
  for (j in 1:numspgr_ny.tb4)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_ny.tb4)==as.character(spgrnames_ny.tb4[i]))
      indexJ=which(as.character(spatiallist_ny.tb4)==as.character(spgrnames_ny.tb4[j]))
      #subsetDistance Matrix
      subDistMat<-ny_sqbrdist.tb4.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_ny.tb4[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_ny.tb4[i,j]<0){PhiMatrix_ny.tb4[i,j]=0}
      PvalMatrix_ny.tb4[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_ny.tb4$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_ny.tb4[i]))],na.rm=TRUE)
  XYT_ny.tb4$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_ny.tb4[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_ny.tb4)<-spgrnames_ny.tb4
colnames(PhiMatrix_ny.tb4)<-spgrnames_ny.tb4
row.names(PvalMatrix_ny.tb4)<-spgrnames_ny.tb4
colnames(PvalMatrix_ny.tb4)<-spgrnames_ny.tb4
#NY Timeblock 5
XYT_ny.tb5<-data.frame(x=numeric(length=numspgr_ny.tb5),y=numeric(length=numspgr_ny.tb5),t=numeric(length=numspgr_ny.tb5))
for (i in 1:numspgr_ny.tb5)
{
  print(paste(i," out of ",numspgr_ny.tb5))
  for (j in 1:numspgr_ny.tb5)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_ny.tb5)==as.character(spgrnames_ny.tb5[i]))
      indexJ=which(as.character(spatiallist_ny.tb5)==as.character(spgrnames_ny.tb5[j]))
      #subsetDistance Matrix
      subDistMat<-ny_sqbrdist.tb5.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_ny.tb5[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_ny.tb5[i,j]<0){PhiMatrix_ny.tb5[i,j]=0}
      PvalMatrix_ny.tb5[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_ny.tb5$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_ny.tb5[i]))],na.rm=TRUE)
  XYT_ny.tb5$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_ny.tb5[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_ny.tb5)<-spgrnames_ny.tb5
colnames(PhiMatrix_ny.tb5)<-spgrnames_ny.tb5
row.names(PvalMatrix_ny.tb5)<-spgrnames_ny.tb5
colnames(PvalMatrix_ny.tb5)<-spgrnames_ny.tb5

#Ontario Timeblock 1
XYT_on.tb1<-data.frame(x=numeric(length=numspgr_on.tb1),y=numeric(length=numspgr_on.tb1),t=numeric(length=numspgr_on.tb1))
for (i in 1:numspgr_on.tb1)
{
  print(paste(i," out of ",numspgr_on.tb1))
  for (j in 1:numspgr_on.tb1)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_on.tb1)==as.character(spgrnames_on.tb1[i]))
      indexJ=which(as.character(spatiallist_on.tb1)==as.character(spgrnames_on.tb1[j]))
      #subsetDistance Matrix
      subDistMat<-on_sqbrdist.tb1.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_on.tb1[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_on.tb1[i,j]<0){PhiMatrix_on.tb1[i,j]=0}
      PvalMatrix_on.tb1[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_on.tb1$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_on.tb1[i]))],na.rm=TRUE)
  XYT_on.tb1$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_on.tb1[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_on.tb1)<-spgrnames_on.tb1
colnames(PhiMatrix_on.tb1)<-spgrnames_on.tb1
row.names(PvalMatrix_on.tb1)<-spgrnames_on.tb1
colnames(PvalMatrix_on.tb1)<-spgrnames_on.tb1

#Ontario Timeblock 2
XYT_on.tb2<-data.frame(x=numeric(length=numspgr_on.tb2),y=numeric(length=numspgr_on.tb2),t=numeric(length=numspgr_on.tb2))
for (i in 1:numspgr_on.tb2)
{
  print(paste(i," out of ",numspgr_on.tb2))
  for (j in 1:numspgr_on.tb2)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_on.tb2)==as.character(spgrnames_on.tb2[i]))
      indexJ=which(as.character(spatiallist_on.tb2)==as.character(spgrnames_on.tb2[j]))
      #subsetDistance Matrix
      subDistMat<-on_sqbrdist.tb2.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_on.tb2[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_on.tb2[i,j]<0){PhiMatrix_on.tb2[i,j]=0}
      PvalMatrix_on.tb2[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_on.tb2$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_on.tb2[i]))],na.rm=TRUE)
  XYT_on.tb2$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_on.tb2[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_on.tb2)<-spgrnames_on.tb2
colnames(PhiMatrix_on.tb2)<-spgrnames_on.tb2
row.names(PvalMatrix_on.tb2)<-spgrnames_on.tb2
colnames(PvalMatrix_on.tb2)<-spgrnames_on.tb2

#Ontario Timeblock 3
XYT_on.tb3<-data.frame(x=numeric(length=numspgr_on.tb3),y=numeric(length=numspgr_on.tb3),t=numeric(length=numspgr_on.tb3))
for (i in 1:numspgr_on.tb3)
{
  print(paste(i," out of ",numspgr_on.tb3))
  for (j in 1:numspgr_on.tb3)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_on.tb3)==as.character(spgrnames_on.tb3[i]))
      indexJ=which(as.character(spatiallist_on.tb3)==as.character(spgrnames_on.tb3[j]))
      #subsetDistance Matrix
      subDistMat<-on_sqbrdist.tb3.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_on.tb3[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_on.tb3[i,j]<0){PhiMatrix_on.tb3[i,j]=0}
      PvalMatrix_on.tb3[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_on.tb3$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_on.tb3[i]))],na.rm=TRUE)
  XYT_on.tb3$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_on.tb3[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_on.tb3)<-spgrnames_on.tb3
colnames(PhiMatrix_on.tb3)<-spgrnames_on.tb3
row.names(PvalMatrix_on.tb3)<-spgrnames_on.tb3
colnames(PvalMatrix_on.tb3)<-spgrnames_on.tb3

#Ontario Timeblock 4
XYT_on.tb4<-data.frame(x=numeric(length=numspgr_on.tb4),y=numeric(length=numspgr_on.tb4),t=numeric(length=numspgr_on.tb4))
for (i in 1:numspgr_on.tb4)
{
  print(paste(i," out of ",numspgr_on.tb4))
  for (j in 1:numspgr_on.tb4)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_on.tb4)==as.character(spgrnames_on.tb4[i]))
      indexJ=which(as.character(spatiallist_on.tb4)==as.character(spgrnames_on.tb4[j]))
      #subsetDistance Matrix
      subDistMat<-on_sqbrdist.tb4.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_on.tb4[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_on.tb4[i,j]<0){PhiMatrix_on.tb4[i,j]=0}
      PvalMatrix_on.tb4[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_on.tb4$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_on.tb4[i]))],na.rm=TRUE)
  XYT_on.tb4$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_on.tb4[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_on.tb4)<-spgrnames_on.tb4
colnames(PhiMatrix_on.tb4)<-spgrnames_on.tb4
row.names(PvalMatrix_on.tb4)<-spgrnames_on.tb4
colnames(PvalMatrix_on.tb4)<-spgrnames_on.tb4

#Ontario Timeblock 5
XYT_on.tb5<-data.frame(x=numeric(length=numspgr_on.tb5),y=numeric(length=numspgr_on.tb5),t=numeric(length=numspgr_on.tb5))
for (i in 1:numspgr_on.tb5)
{
  print(paste(i," out of ",numspgr_on.tb5))
  for (j in 1:numspgr_on.tb5)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_on.tb5)==as.character(spgrnames_on.tb5[i]))
      indexJ=which(as.character(spatiallist_on.tb5)==as.character(spgrnames_on.tb5[j]))
      #subsetDistance Matrix
      subDistMat<-on_sqbrdist.tb5.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_on.tb5[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_on.tb5[i,j]<0){PhiMatrix_on.tb5[i,j]=0}
      PvalMatrix_on.tb5[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_on.tb5$x[i]<-mean(dt$lon_dd[which(dt$groupcode==as.character(spgrnames_on.tb5[i]))],na.rm=TRUE)
  XYT_on.tb5$y[i]<-mean(dt$lat_dd[which(dt$groupcode==as.character(spgrnames_on.tb5[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_on.tb5)<-spgrnames_on.tb5
colnames(PhiMatrix_on.tb5)<-spgrnames_on.tb5
row.names(PvalMatrix_on.tb5)<-spgrnames_on.tb5
colnames(PvalMatrix_on.tb5)<-spgrnames_on.tb5


#St. Lawrence All Sites
XYT_stl<-data.frame(x=numeric(length=numspgr_stl),y=numeric(length=numspgr_stl),t=numeric(length=numspgr_stl))
for (i in 1:numspgr_stl)
{
  print(paste(i," out of ",numspgr_stl))
  for (j in 1:numspgr_stl)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_stl)==as.character(spgrnames_stl[i]))
      indexJ=which(as.character(spatiallist_stl)==as.character(spgrnames_stl[j]))
      #subsetDistance Matrix
      subDistMat<-stl_sqbrdist.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_stl[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_stl[i,j]<0){PhiMatrix_stl[i,j]=0}
      PvalMatrix_stl[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_stl$x[i]<-mean(stl.dt$lon_dd[which(stl.dt$groupcode==as.character(spgrnames_stl[i]))],na.rm=TRUE)
  XYT_stl$y[i]<-mean(stl.dt$lat_dd[which(stl.dt$groupcode==as.character(spgrnames_stl[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_stl)<-spgrnames_stl
colnames(PhiMatrix_stl)<-spgrnames_stl
row.names(PvalMatrix_stl)<-spgrnames_stl
colnames(PvalMatrix_stl)<-spgrnames_stl

#St. Lawrence Timeblock 1
XYT_stl.tb1<-data.frame(x=numeric(length=numspgr_stl.tb1),y=numeric(length=numspgr_stl.tb1),t=numeric(length=numspgr_stl.tb1))
for (i in 1:numspgr_stl.tb1)
{
  print(paste(i," out of ",numspgr_stl.tb1))
  for (j in 1:numspgr_stl.tb1)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_stl.tb1)==as.character(spgrnames_stl.tb1[i]))
      indexJ=which(as.character(spatiallist_stl.tb1)==as.character(spgrnames_stl.tb1[j]))
      #subsetDistance Matrix
      subDistMat<-stl_sqbrdist.tb1.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_stl.tb1[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_stl.tb1[i,j]<0){PhiMatrix_stl.tb1[i,j]=0}
      PvalMatrix_stl.tb1[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_stl.tb1$x[i]<-mean(stl.dt$lon_dd[which(stl.dt$groupcode==as.character(spgrnames_stl.tb1[i]))],na.rm=TRUE)
  XYT_stl.tb1$y[i]<-mean(stl.dt$lat_dd[which(stl.dt$groupcode==as.character(spgrnames_stl.tb1[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_stl.tb1)<-spgrnames_stl.tb1
colnames(PhiMatrix_stl.tb1)<-spgrnames_stl.tb1
row.names(PvalMatrix_stl.tb1)<-spgrnames_stl.tb1
colnames(PvalMatrix_stl.tb1)<-spgrnames_stl.tb1

#St. Lawrence Timeblock 2
XYT_stl.tb2<-data.frame(x=numeric(length=numspgr_stl.tb2),y=numeric(length=numspgr_stl.tb2),t=numeric(length=numspgr_stl.tb2))
for (i in 1:numspgr_stl.tb2)
{
  print(paste(i," out of ",numspgr_stl.tb2))
  for (j in 1:numspgr_stl.tb2)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_stl.tb2)==as.character(spgrnames_stl.tb2[i]))
      indexJ=which(as.character(spatiallist_stl.tb2)==as.character(spgrnames_stl.tb2[j]))
      #subsetDistance Matrix
      subDistMat<-stl_sqbrdist.tb2.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_stl.tb2[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_stl.tb2[i,j]<0){PhiMatrix_stl.tb2[i,j]=0}
      PvalMatrix_stl.tb2[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_stl.tb2$x[i]<-mean(stl.dt$lon_dd[which(stl.dt$groupcode==as.character(spgrnames_stl.tb2[i]))],na.rm=TRUE)
  XYT_stl.tb2$y[i]<-mean(stl.dt$lat_dd[which(stl.dt$groupcode==as.character(spgrnames_stl.tb2[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_stl.tb2)<-spgrnames_stl.tb2
colnames(PhiMatrix_stl.tb2)<-spgrnames_stl.tb2
row.names(PvalMatrix_stl.tb2)<-spgrnames_stl.tb2
colnames(PvalMatrix_stl.tb2)<-spgrnames_stl.tb2

#St Lawrence Timeblock 3
XYT_stl.tb3<-data.frame(x=numeric(length=numspgr_stl.tb3),y=numeric(length=numspgr_stl.tb3),t=numeric(length=numspgr_stl.tb3))
for (i in 1:numspgr_stl.tb3)
{
  print(paste(i," out of ",numspgr_stl.tb3))
  for (j in 1:numspgr_stl.tb3)
  {
    if (i!=j&i>j)
    { 
      indexI=which(as.character(spatiallist_stl.tb3)==as.character(spgrnames_stl.tb3[i]))
      indexJ=which(as.character(spatiallist_stl.tb3)==as.character(spgrnames_stl.tb3[j]))
      #subsetDistance Matrix
      subDistMat<-stl_sqbrdist.tb3.matrix[c(indexI,indexJ),c(indexI,indexJ)]
      subDistMat<-as.dist(subDistMat)
      classes<-as.factor(c(rep("a",length(indexI)),rep("b",length(indexJ))))
      #compute Phi:
      tmp<-amova(subDistMat~classes,nperm = 10000)
      N0=tmp$varcoef
      AP=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)
      WP=(tmp$tab$MSD[2])
      PhiMatrix_stl.tb3[i,j]=((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef)/((tmp$tab$MSD[2])+((tmp$tab$MSD[1]-tmp$tab$MSD[2])/tmp$varcoef))
      if(PhiMatrix_stl.tb3[i,j]<0){PhiMatrix_stl.tb3[i,j]=0}
      PvalMatrix_stl.tb3[i,j]=tmp$varcomp$P.value[1]
    }
  }
  XYT_stl.tb3$x[i]<-mean(stl.dt$lon_dd[which(stl.dt$groupcode==as.character(spgrnames_stl.tb3[i]))],na.rm=TRUE)
  XYT_stl.tb3$y[i]<-mean(stl.dt$lat_dd[which(stl.dt$groupcode==as.character(spgrnames_stl.tb3[i]))],na.rm=TRUE)
}
row.names(PhiMatrix_stl.tb3)<-spgrnames_stl.tb3
colnames(PhiMatrix_stl.tb3)<-spgrnames_stl.tb3
row.names(PvalMatrix_stl.tb3)<-spgrnames_stl.tb3
colnames(PvalMatrix_stl.tb3)<-spgrnames_stl.tb3

############################################
############################################
########Step 3: Evaluate Biases#############
############################################
############################################

### Stratified Mantel tests ###
#stratified tests only work for populations with more than one group
#purpose of this test is to remove the potential bias of geographic clustering on apparent Mantel correlation
#A persistent result after stratification means there is a robust overall correlation with distance

mantel(brdist_global,geodist_global,method = "spearman",strata = groups_global,permutations = 10000) #permute rows by the same grouping factor as AMOVA

mantel(timeblock1.brdist,timeblock1.geodist,method = "spearman",strata = groups_tb1,permutations = 10000)
mantel(timeblock2.brdist,timeblock2.geodist,method = "spearman",strata = groups_tb2,permutations = 10000)
mantel(timeblock3.brdist,timeblock3.geodist,method = "spearman",strata = groups_tb3,permutations = 10000)
mantel(timeblock4.brdist,timeblock4.geodist,method = "spearman",strata = groups_tb4,permutations = 10000)
mantel(timeblock5.brdist,timeblock5.geodist,method = "spearman",strata = groups_tb5,permutations = 10000)

mantel(ny_brdist,ny_geodist,method = "spearman",strata = groups_ny,permutations = 10000)
mantel(on_brdist,on_geodist,method = "spearman",strata = groups_on,permutations = 10000)
mantel(stl_brdist,stl_brdist,method = "spearman",strata = groups_stl,permutations = 10000)

mantel(ny_brdist.tb1,ny_geodist.tb1,method = "spearman",strata = groups_ny.tb1,permutations = 10000)
mantel(ny_brdist.tb2,ny_geodist.tb2,method = "spearman",strata = groups_ny.tb2,permutations = 10000)
mantel(ny_brdist.tb3,ny_geodist.tb3,method = "spearman",strata = groups_ny.tb3,permutations = 10000)
mantel(ny_brdist.tb4,ny_geodist.tb4,method = "spearman",strata = groups_ny.tb4,permutations = 10000)
mantel(ny_brdist.tb5,ny_geodist.tb5,method = "spearman",strata = groups_ny.tb5,permutations = 10000)

mantel(on_brdist.tb1,on_geodist.tb1,method = "spearman",strata = groups_on.tb1,permutations = 10000)
mantel(on_brdist.tb2,on_geodist.tb2,method = "spearman",strata = groups_on.tb2,permutations = 10000)
mantel(on_brdist.tb3,on_geodist.tb3,method = "spearman",strata = groups_on.tb3,permutations = 10000)
mantel(on_brdist.tb4,on_geodist.tb4,method = "spearman",strata = groups_on.tb4,permutations = 10000)
mantel(on_brdist.tb5,on_geodist.tb5,method = "spearman",strata = groups_on.tb5,permutations = 10000)

mantel(stl_brdist.tb1,stl_geodist.tb1,method = "spearman",strata = groups_stl.tb1,permutations = 10000)
mantel(stl_brdist.tb2,stl_geodist.tb2,method = "spearman",strata = groups_stl.tb2,permutations = 10000)
mantel(stl_brdist.tb3,stl_geodist.tb3,method = "spearman",strata = groups_stl.tb3,permutations = 10000)


### Pairwise PhiST and Geospatial Distance ###
#Purpose of this test is to quantify how strongly apparent differentiation between groups is influenced by their geographic separation (between group centroids).
#A lack of robust relationship means that significant pairwise Phi between groups are more likely to be driven by real patterning in the data

PhiMatrix_global.dist <- as.dist(PhiMatrix_global) #convert pairwise Phi matrix to distance matrix (larger values = greater differentiation between groups)

PhiMatrix_tb1.dist <- as.dist(PhiMatrix_tb1)
PhiMatrix_tb2.dist <- as.dist(PhiMatrix_tb2)
PhiMatrix_tb3.dist <- as.dist(PhiMatrix_tb3)
PhiMatrix_tb4.dist <- as.dist(PhiMatrix_tb4)
PhiMatrix_tb5.dist <- as.dist(PhiMatrix_tb5)

PhiMatrix_stl.dist <- as.dist(PhiMatrix_stl)
PhiMatrix_ny.dist <- as.dist(PhiMatrix_ny)
PhiMatrix_on.dist <- as.dist(PhiMatrix_on)

PhiMatrix_ny.tb1.dist <- as.dist(PhiMatrix_ny.tb1)
PhiMatrix_ny.tb2.dist <- as.dist(PhiMatrix_ny.tb2)
PhiMatrix_ny.tb3.dist <- as.dist(PhiMatrix_ny.tb3)
PhiMatrix_ny.tb4.dist <- as.dist(PhiMatrix_ny.tb4)
PhiMatrix_ny.tb5.dist <- as.dist(PhiMatrix_ny.tb5)

PhiMatrix_on.tb1.dist <- as.dist(PhiMatrix_on.tb1)
PhiMatrix_on.tb2.dist <- as.dist(PhiMatrix_on.tb2)
PhiMatrix_on.tb3.dist <- as.dist(PhiMatrix_on.tb3)
PhiMatrix_on.tb4.dist <- as.dist(PhiMatrix_on.tb4)
PhiMatrix_on.tb5.dist <- as.dist(PhiMatrix_on.tb5)

PhiMatrix_stl.tb1.dist <- as.dist(PhiMatrix_stl.tb1)
PhiMatrix_stl.tb2.dist <- as.dist(PhiMatrix_stl.tb2)
PhiMatrix_stl.tb3.dist <- as.dist(PhiMatrix_stl.tb3)

spatialDist_global <- distMat(XYT_global)

spatialDist_tb1 <- distMat(XYT_tb1)
spatialDist_tb2 <- distMat(XYT_tb2)
spatialDist_tb3 <- distMat(XYT_tb3)
spatialDist_tb4 <- distMat(XYT_tb4)
spatialDist_tb5 <- distMat(XYT_tb5)

spatialDist_stl <- distMat(XYT_stl)
spatialDist_ny <- distMat(XYT_ny)
spatialDist_on <- distMat(XYT_on)

spatialDist_ny.tb1 <- distMat(XYT_ny.tb1)
spatialDist_ny.tb2 <- distMat(XYT_ny.tb2)
spatialDist_ny.tb3 <- distMat(XYT_ny.tb3)
spatialDist_ny.tb4 <- distMat(XYT_ny.tb4)
spatialDist_ny.tb5 <- distMat(XYT_ny.tb5)

spatialDist_on.tb1 <- distMat(XYT_on.tb1)
spatialDist_on.tb2 <- distMat(XYT_on.tb2)
spatialDist_on.tb3 <- distMat(XYT_on.tb3)
spatialDist_on.tb4 <- distMat(XYT_on.tb4)
spatialDist_on.tb5 <- distMat(XYT_on.tb5)

spatialDist_stl.tb1 <- distMat(XYT_stl.tb1)
spatialDist_stl.tb2 <- distMat(XYT_stl.tb2)
spatialDist_stl.tb3 <- distMat(XYT_stl.tb3)

#All spatial Groups All Timeblocks#
mantel(PhiMatrix_global.dist,spatialDist_global,method = "spearman",permutations = 10000,parallel = 4)

#All spatial Groups by Timeblock#
mantel(PhiMatrix_tb1,spatialDist_tb1,method = "spearman",permutations = 10000,parallel = 4)
mantel(PhiMatrix_tb2,spatialDist_tb2,method = "spearman",permutations = 10000,parallel = 4)
mantel(PhiMatrix_tb3,spatialDist_tb3,method = "spearman",permutations = 10000,parallel = 4)
mantel(PhiMatrix_tb4,spatialDist_tb4,method = "spearman",permutations = 10000,parallel = 4)
mantel(PhiMatrix_tb5,spatialDist_tb5,method = "spearman",permutations = 10000,parallel = 4)

#All Timeblocks All spatial groups by combined group#
mantel(PhiMatrix_ny.dist,spatialDist_ny,method = "spearman",permutations = 10000)
mantel(PhiMatrix_stl.dist,spatialDist_stl,method = "spearman",permutations = 10000)
mantel(PhiMatrix_on.dist,spatialDist_on,method = "spearman",permutations = 10000)

#all spatial groups by Timeblock in NY combined group#
mantel(PhiMatrix_ny.tb1.dist,spatialDist_ny.tb1,method = "spearman",permutations = 10000)
mantel(PhiMatrix_ny.tb2.dist,spatialDist_ny.tb2,method = "spearman",permutations = 10000)
mantel(PhiMatrix_ny.tb3.dist,spatialDist_ny.tb3,method = "spearman",permutations = 10000)
mantel(PhiMatrix_ny.tb4.dist,spatialDist_ny.tb4,method = "spearman",permutations = 10000)
mantel(PhiMatrix_ny.tb5.dist,spatialDist_ny.tb5,method = "spearman",permutations = 10000)

#all spatial groups by Timeblock in ON combined group#
mantel(PhiMatrix_on.tb1.dist,spatialDist_on.tb1,method = "spearman",permutations = 10000)
mantel(PhiMatrix_on.tb2.dist,spatialDist_on.tb2,method = "spearman",permutations = 10000)
mantel(PhiMatrix_on.tb3.dist,spatialDist_on.tb3,method = "spearman",permutations = 10000)
mantel(PhiMatrix_on.tb4.dist,spatialDist_on.tb4,method = "spearman",permutations = 10000)
mantel(PhiMatrix_on.tb5.dist,spatialDist_on.tb5,method = "spearman",permutations = 10000)

#all spatial groups by Timeblock in STL combined group#
mantel(PhiMatrix_stl.tb1.dist,spatialDist_stl.tb1,method = "spearman",permutations = 10000)
mantel(PhiMatrix_stl.tb2.dist,spatialDist_stl.tb2,method = "spearman",permutations = 10000)
mantel(PhiMatrix_stl.tb3.dist,spatialDist_stl.tb3,method = "spearman",permutations = 10000)

###Partial Mantel Test Binary Distance Matrix spatial Group Affiliation ~ BR Holding Geographic Distance###
#purpose of this test is to test the correlation between spatial group affiliation and BR similarity by controlling for the confounding influence of distance between sites
#sites in the same group = 1 sites in other groups = 0
#Create binary distance matrix
#all sites all time periods group by spatial
binaryspgr_global<-as.dist(BinaryCulture(dt$groupcode))#generate binary distance matrix function written by Crema

#all sites all time periods group by regional group
binarycgroup_global<-as.dist(BinaryCulture(dt$cgroup_code))

#all sites by timeblock grouped by spatial
binaryspgr_tb1 <- as.dist(BinaryCulture(dt$groupcode[dt$timeslicecode == "I" | dt$timeslicecode == "II"]))
binaryspgr_tb2 <- as.dist(BinaryCulture(dt$groupcode[dt$timeslicecode == "II" | dt$timeslicecode == "III"]))
binaryspgr_tb3 <- as.dist(BinaryCulture(dt$groupcode[dt$timeslicecode == "III" | dt$timeslicecode == "IV"]))
binaryspgr_tb4 <- as.dist(BinaryCulture(dt$groupcode[dt$timeslicecode == "IV" | dt$timeslicecode == "V"]))
binaryspgr_tb5 <- as.dist(BinaryCulture(dt$groupcode[dt$timeslicecode == "V" | dt$timeslicecode == "VI"]))

#all sites by timeblock grouped by regional group
binaryspgr_tb1.cgroups <- as.dist(BinaryCulture(dt$cgroup_code[dt$timeslicecode == "I" | dt$timeslicecode == "II"]))
binaryspgr_tb2.cgroups <- as.dist(BinaryCulture(dt$cgroup_code[dt$timeslicecode == "II" | dt$timeslicecode == "III"]))
binaryspgr_tb3.cgroups <- as.dist(BinaryCulture(dt$cgroup_code[dt$timeslicecode == "III" | dt$timeslicecode == "IV"]))
binaryspgr_tb4.cgroups <- as.dist(BinaryCulture(dt$cgroup_code[dt$timeslicecode == "IV" | dt$timeslicecode == "V"]))
binaryspgr_tb5.cgroups <- as.dist(BinaryCulture(dt$cgroup_code[dt$timeslicecode == "V" | dt$timeslicecode == "VI"]))

#regional group sites all time periods grouped by spatial
binaryspgr_ny<-as.dist(BinaryCulture(dt$groupcode[dt$cgroup_code=="C"])) #NY combined group
binaryspgr_on<-as.dist(BinaryCulture(dt$groupcode[dt$cgroup_code=="D"])) #ON combined group
binaryspgr_stl<-as.dist(BinaryCulture(dt$groupcode[dt$cgroup_code=="E"])) #STL combined group

#NY regional group by timeblock grouped by spatial group
binaryspgr_ny.tb1 <- as.dist(BinaryCulture(ny.dt$groupcode[ny.dt$timeslicecode == "I" | ny.dt$timeslicecode == "II"]))
binaryspgr_ny.tb2 <- as.dist(BinaryCulture(ny.dt$groupcode[ny.dt$timeslicecode == "II" | ny.dt$timeslicecode == "III"]))
binaryspgr_ny.tb3 <- as.dist(BinaryCulture(ny.dt$groupcode[ny.dt$timeslicecode == "III" | ny.dt$timeslicecode == "IV"]))
binaryspgr_ny.tb4 <- as.dist(BinaryCulture(ny.dt$groupcode[ny.dt$timeslicecode == "IV" | ny.dt$timeslicecode == "V"]))
binaryspgr_ny.tb5 <- as.dist(BinaryCulture(ny.dt$groupcode[ny.dt$timeslicecode == "V" | ny.dt$timeslicecode == "VI"]))

#ON regional group by timeblock grouped by spatial group
binaryspgr_on.tb1 <- as.dist(BinaryCulture(on.dt$groupcode[on.dt$timeslicecode == "I" | on.dt$timeslicecode == "II"]))
binaryspgr_on.tb2 <- as.dist(BinaryCulture(on.dt$groupcode[on.dt$timeslicecode == "II" | on.dt$timeslicecode == "III"]))
binaryspgr_on.tb3 <- as.dist(BinaryCulture(on.dt$groupcode[on.dt$timeslicecode == "III" | on.dt$timeslicecode == "IV"]))
binaryspgr_on.tb4 <- as.dist(BinaryCulture(on.dt$groupcode[on.dt$timeslicecode == "IV" | on.dt$timeslicecode == "V"]))
binaryspgr_on.tb5 <- as.dist(BinaryCulture(on.dt$groupcode[on.dt$timeslicecode == "V" | on.dt$timeslicecode == "VI"]))

#STL regional group by timeblock grouped by spatial group
binaryspgr_stl.tb1 <- as.dist(BinaryCulture(stl.dt$groupcode[stl.dt$timeslicecode == "I" | stl.dt$timeslicecode == "II"]))
binaryspgr_stl.tb2 <- as.dist(BinaryCulture(stl.dt$groupcode[stl.dt$timeslicecode == "II" | stl.dt$timeslicecode == "III"]))
binaryspgr_stl.tb3 <- as.dist(BinaryCulture(stl.dt$groupcode[stl.dt$timeslicecode == "III" | stl.dt$timeslicecode == "IV"]))

#Partial Mantel Test X=BR distance,Y= binary distance matrix of spatial group affiliation, Z=Geographic distance#
#all sites all timeblocks by spatial group
mantel.partial(brdist_global,binaryspgr_global,geodist_global,method = "spearman",permutations = 10000,parallel = 4)

#all sites all timeblocks by regional group
mantel.partial(brdist_global,binarycgroup_global,geodist_global,method = "spearman",permutations = 10000,parallel = 4)

#all sites by timeblock grouped by spatial group
mantel.partial(timeblock1.brdist,binaryspgr_tb1,timeblock1.geodist,method = "spearman",permutations = 10000,parallel = 4)
mantel.partial(timeblock2.brdist,binaryspgr_tb2,timeblock2.geodist,method = "spearman",permutations = 10000,parallel = 4)
mantel.partial(timeblock3.brdist,binaryspgr_tb3,timeblock3.geodist,method = "spearman",permutations = 10000,parallel = 4)
mantel.partial(timeblock4.brdist,binaryspgr_tb4,timeblock4.geodist,method = "spearman",permutations = 10000,parallel = 4)
mantel.partial(timeblock5.brdist,binaryspgr_tb5,timeblock5.geodist,method = "spearman",permutations = 10000,parallel = 4)

#all sites by timeblock grouped by region
mantel.partial(timeblock1.brdist,binaryspgr_tb1.cgroups,timeblock1.geodist,method = "spearman",permutations = 10000,parallel = 4)
mantel.partial(timeblock2.brdist,binaryspgr_tb2.cgroups,timeblock2.geodist,method = "spearman",permutations = 10000,parallel = 4)
mantel.partial(timeblock3.brdist,binaryspgr_tb3.cgroups,timeblock3.geodist,method = "spearman",permutations = 10000,parallel = 4)
mantel.partial(timeblock4.brdist,binaryspgr_tb4.cgroups,timeblock4.geodist,method = "spearman",permutations = 10000,parallel = 4)
mantel.partial(timeblock5.brdist,binaryspgr_tb5.cgroups,timeblock5.geodist,method = "spearman",permutations = 10000,parallel = 4)

#Regional groups all timeblocks grouped by spatial group
mantel.partial(stl_brdist,binaryspgr_stl,stl_geodist,method = "spearman",permutations = 10000,parallel = 4) #NY group
mantel.partial(ny_brdist,binaryspgr_ny,ny_geodist,method = "spearman",permutations = 10000,parallel = 4) #ON group
mantel.partial(on_brdist,binaryspgr_on,on_geodist,method = "spearman",permutations = 10000,parallel = 4) #STL group

#NY regional group by timeblock grouped by spatial group
mantel.partial(ny_brdist.tb1,binaryspgr_ny.tb1,ny_geodist.tb1,method = "spearman",permutations = 10000,parallel = 4)
mantel.partial(ny_brdist.tb2,binaryspgr_ny.tb2,ny_geodist.tb2,method = "spearman",permutations = 10000,parallel = 4)
mantel.partial(ny_brdist.tb3,binaryspgr_ny.tb3,ny_geodist.tb3,method = "spearman",permutations = 10000,parallel = 4)
mantel.partial(ny_brdist.tb4,binaryspgr_ny.tb4,ny_geodist.tb4,method = "spearman",permutations = 10000,parallel = 4)
mantel.partial(ny_brdist.tb5,binaryspgr_ny.tb5,ny_geodist.tb5,method = "spearman",permutations = 10000,parallel = 4)

#ON regional group by timeblock grouped by spatial group
mantel.partial(on_brdist.tb1,binaryspgr_on.tb1,on_geodist.tb1,method = "spearman",permutations = 10000,parallel = 4)
mantel.partial(on_brdist.tb2,binaryspgr_on.tb2,on_geodist.tb2,method = "spearman",permutations = 10000,parallel = 4)
mantel.partial(on_brdist.tb3,binaryspgr_on.tb3,on_geodist.tb3,method = "spearman",permutations = 10000,parallel = 4)
mantel.partial(on_brdist.tb4,binaryspgr_on.tb4,on_geodist.tb4,method = "spearman",permutations = 10000,parallel = 4)
mantel.partial(on_brdist.tb5,binaryspgr_on.tb5,on_geodist.tb5,method = "spearman",permutations = 10000,parallel = 4)

#STL regional group by timeblock grouped by spatial group
mantel.partial(stl_brdist.tb1,binaryspgr_stl.tb1,stl_geodist.tb1,method = "spearman",permutations = 10000,parallel = 4)
mantel.partial(stl_brdist.tb2,binaryspgr_stl.tb2,stl_geodist.tb2,method = "spearman",permutations = 10000,parallel = 4)
mantel.partial(stl_brdist.tb3,binaryspgr_stl.tb3,stl_geodist.tb3,method = "spearman",permutations = 10000,parallel = 4)

### ANOSIM ###
#The purpose of this test is to measure the relative separation in the central tendencies of *ranked* dissimilarities of all sites in each spatial group
#Basically an ANOVA for distance matrices; assumes equal ranges of ranked dissimilarities within groups

#all sites all timeblocks grouped by spatial group
anosim_glob <- anosim(brdist_global,grouping = groups_global,permutations = 10000)
#all sites all timeblocks grouped by regional group
anosim_glob.cgroups <- anosim(brdist_global,grouping = cgroups_global,permutations = 10000)

#all sites by timeblock grouped by spatial group
anosim_tb1 <- anosim(timeblock1.brdist,grouping = groups_tb1,permutations = 10000)
anosim_tb2 <- anosim(timeblock2.brdist,grouping = groups_tb2,permutations = 10000)
anosim_tb3 <- anosim(timeblock3.brdist,grouping = groups_tb3,permutations = 10000)
anosim_tb4 <- anosim(timeblock4.brdist,grouping = groups_tb4,permutations = 10000)
anosim_tb5 <- anosim(timeblock5.brdist,grouping = groups_tb5,permutations = 10000)

#all sites by timeblock grouped by regional group
anosim_tb1.cgroups <- anosim(timeblock1.brdist,grouping = cgroups_tb1,permutations = 10000)
anosim_tb2.cgroups <- anosim(timeblock2.brdist,grouping = cgroups_tb2,permutations = 10000)
anosim_tb3.cgroups <- anosim(timeblock3.brdist,grouping = cgroups_tb3,permutations = 10000)
anosim_tb4.cgroups <- anosim(timeblock4.brdist,grouping = cgroups_tb4,permutations = 10000)
anosim_tb5.cgroups <- anosim(timeblock5.brdist,grouping = cgroups_tb5,permutations = 10000)

#regional groups all timeblocks grouped by spatial group
anosim_ny <- anosim(ny_brdist,grouping = groups_ny,permutations = 10000) #NY group
anosim_stl <- anosim(stl_brdist,grouping = groups_stl,permutations = 10000) #ON group
anosim_on <- anosim(on_brdist,grouping = groups_on,permutations = 10000) #STL group

#NY combined group by timeblocks grouped by spatial groups
anosim_ny.tb1 <- anosim(ny_brdist.tb1,grouping = groups_ny.tb1,permutations = 10000)
anosim_ny.tb2 <- anosim(ny_brdist.tb2,grouping = groups_ny.tb2,permutations = 10000)
anosim_ny.tb3 <- anosim(ny_brdist.tb3,grouping = groups_ny.tb3,permutations = 10000)
anosim_ny.tb4 <- anosim(ny_brdist.tb4,grouping = groups_ny.tb4,permutations = 10000)
anosim_ny.tb5 <- anosim(ny_brdist.tb5,grouping = groups_ny.tb5,permutations = 10000)

#ON combined group by timeblocks grouped by spatial groups
anosim_on.tb1 <- anosim(on_brdist.tb1,grouping = groups_on.tb1,permutations = 10000)
anosim_on.tb2 <- anosim(on_brdist.tb2,grouping = groups_on.tb2,permutations = 10000)
anosim_on.tb3 <- anosim(on_brdist.tb3,grouping = groups_on.tb3,permutations = 10000)
anosim_on.tb4 <- anosim(on_brdist.tb4,grouping = groups_on.tb4,permutations = 10000)
anosim_on.tb5 <- anosim(on_brdist.tb5,grouping = groups_on.tb5,permutations = 10000)

#STL combined group by timeblocks grouped by spatial groups
anosim_stl.tb1 <- anosim(stl_brdist.tb1,grouping = groups_stl.tb1,permutations = 10000)
anosim_stl.tb2 <- anosim(stl_brdist.tb2,grouping = groups_stl.tb2,permutations = 10000)
anosim_stl.tb3 <- anosim(stl_brdist.tb3,grouping = groups_stl.tb3,permutations = 10000)

#generate all plots
source("plots.R",echo = T) #make sure script is in working directory