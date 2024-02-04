# Code to count genera in four guild categories and compare site faunal composition to modern communities
# first set working directory to where data are saved
setwd("C:/Users/abbap/OneDrive - University of Helsinki/BKV SizeDietBias GitHub")
#libraries
library(sf)
library(tidyverse)
library(RColorBrewer)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(RANN)
library(scico)
library(aplot)
library(mixtools)
library(gllvm)
library(geometry)
theme_set(theme_bw())

# Count large/small, herb/nonherb genera at fossil sites ####
# Define large and small mammals

#occurring genera, as defined in 1_SiteOccurrenceMatrix
mmass_occ %>% count(LargeSmall) #how many of each, in occurring fossils? 
# This uses 1.76kg cutoff defined for modern genera in Eurasia (above)
colnames(mmass_occ)
lars=mmass_occ[mmass_occ$LargeSmall=="Large",1] #list of large genera
smas <- mmass_occ[mmass_occ$LargeSmall=="Small",1] # list of small genera
# select these genus names, plus other information in site/occurrence matrix
lar <- c("NAME", lars[!(is.na(lars))], colnames(SO[454:ncol(SO)]))  #colnames for large site/occ matrix
sma <- c("NAME", smas[!(is.na(smas))], colnames(SO[454:ncol(SO)])) #colnames for small site/occ matrix

#Name herbivore, non-herbivore
mmass_occ %>% count(Diet)
her=mmass_occ[mmass_occ$Diet=="Herbivore",1] #list of large genera
non <- mmass_occ[mmass_occ$Diet=="Non-Herbivore",1] # list of small genera
# select these genus names, plus other information in site/occurrence matrix
herb <- c("NAME", her[!(is.na(her))], colnames(SO[454:ncol(SO)]))  #colnames for large site/occ matrix
nonh <- c("NAME", non[!(is.na(non))], colnames(SO[454:ncol(SO)])) #colnames for small site/occ matrix

#Count occs in sites
#filter to 4+ genera sites
SO_use <- SO %>% filter(NAME %in% use_subset$NAME)

#Count large genera per site
SO_large <- SO_use[, which((names(SO_use) %in% lar)==TRUE)] #filters out genera <1.4kg
colnames(SO_large)
jdat <- SO_large[,2:234]
jdat[jdat == 0] <- NA
jdat <- jdat %>% mutate(Large_GenCount=rowSums(!(is.na(.)))) #large genus counts per site
SO_use$Large_GenCount <- jdat$Large_GenCount
#count large herbivores and non-herbivores
fossil_largeherb <- jdat[, which((names(jdat) %in% herb)==TRUE)] #filters out genera <1.4kg
fossil_largeherb <- fossil_largeherb %>% mutate(largeherb_GenCount=rowSums(!(is.na(.)))) #herb genus counts per site
SO_use$LargeHerb_GenCount <- fossil_largeherb$largeherb_GenCount
fossil_largenonherb <- jdat[, which((names(jdat) %in% nonh)==TRUE)] #filters out genera <1.4kg
fossil_largenonherb <- fossil_largenonherb %>% mutate(largenonherb_GenCount=rowSums(!(is.na(.)))) #herb genus counts per site
#Assign diet ratio/ proportion within large genera
SO_use$LargeNonherb_GenCount <- fossil_largenonherb$largenonherb_GenCount
SO_use <- SO_use %>% mutate(Large_Herbprop= LargeHerb_GenCount/Large_GenCount)

# count small genera per site
SO_small <- SO_use[, which((names(SO_use) %in% sma)==TRUE)] #filters out genera <1.4kg
colnames(SO_small)
jdat <- SO_small[,2:219]
jdat[jdat == 0] <- NA
jdat <- jdat %>% mutate(Small_GenCount=rowSums(!(is.na(.)))) #small genus counts per site
SO_use$Small_GenCount <- jdat$Small_GenCount
#count small herbivores and non-herbivores
fossil_smallherb <- jdat[, which((names(jdat) %in% herb)==TRUE)] #filters out genera <1.4kg
fossil_smallherb <- fossil_smallherb %>% mutate(smallherb_GenCount=rowSums(!(is.na(.)))) #herb genus counts per site
SO_use$SmallHerb_GenCount <- fossil_smallherb$smallherb_GenCount
fossil_smallnonherb <- jdat[, which((names(jdat) %in% nonh)==TRUE)] #filters out genera <1.4kg
fossil_smallnonherb <- fossil_smallnonherb %>% mutate(smallnonherb_GenCount=rowSums(!(is.na(.)))) #herb genus counts per site
SO_use$SmallNonherb_GenCount <- fossil_smallnonherb$smallnonherb_GenCount
SO_use <- SO_use %>% mutate(Small_Herbprop= SmallHerb_GenCount/Small_GenCount)

#Assign diet ratio/ proportion within small genera
SO_use <- SO_use %>% mutate(smallperlarge= Small_GenCount/Large_GenCount)
SO_use <- SO_use %>% mutate(smallprop= Small_GenCount/Gen_Count)

#Count overall Herbivores and Non-Herbivores per site
#Count Herb and Non-Herb per site
SO_herb <- SO_use[, which((names(SO_use) %in% herb)==TRUE)] #filters out genera <1.4kg
colnames(SO_herb)
jdat <- SO_herb[,2:293]
jdat[jdat == 0] <- NA
jdat <- jdat %>% mutate(Herb_GenCount=rowSums(!(is.na(.)))) #large genus counts per site
SO_use$Herb_GenCount <- jdat$Herb_GenCount

SO_nonherb <- SO_use[, which((names(SO_use) %in% nonh)==TRUE)] #filters out genera <1.4kg
colnames(SO_nonherb)
jdat <- SO_nonherb[,2:145]
jdat[jdat == 0] <- NA
jdat <- jdat %>% mutate(Nonherb_GenCount=rowSums(!(is.na(.)))) #nonherb genus counts per site
SO_use$Nonherb_GenCount <- jdat$Nonherb_GenCount
#Assign diet ratio/ proportion for whole community
SO_use <- SO_use %>% mutate(DietRatio= Herb_GenCount/Nonherb_GenCount)
SO_use <- SO_use %>% mutate(HerbProp= Herb_GenCount/Gen_Count)
SO_use <- SO_use %>% mutate(mid_age= ((MAX_AGE+MIN_AGE)/2))
summary(SO_use[,450:480])
sd(SO_use$Gen_Count)
sd(SO_use$Nonherb_GenCount)
sd(SO_use$smallprop)
sd(SO_use$HerbProp)
sd(SO_use$age_range)
sd(SO_use$mid_age)

#write.csv(SO_use, "Siteoccmatrix_DietSizeCounts.csv")
#SO_use <- read.csv("Siteoccmatrix_DietSizeCounts.csv", header=T, row.names = 1)


# then check number of each guild at large/small/mixed sites
# Input matrices for large and small mammal sites, excluding genera from opposite group
largeonly_mat <- SO_use %>%  filter(NAME %in% largeonly$NAME)
medonly_mat <- SO_use %>% filter(NAME %in% medonly$NAME)
smallonly_mat <- SO_use %>%  filter (NAME %in% smallonly$NAME)
summary(smallsite_mat[,460:474])
summary(largesite_mat[,460:474])
summary(smallonly_mat[,460:474])
summary(medonly_mat[,460:474])
summary(largeonly_mat[,460:473])
sd(largeonly_mat$age_range)
sd(largeonly_mat$mid_age)



# Count large/small, herb/nonherb genera in modern Eco-ISEA3H community samples
#Compare to modern-- small/large and herb/nonherb genus counts ####
#Define large and small genus names
largegen <- phy_gen_av_terr %>% filter(LogMass >3.25) 
smallgen <- phy_gen_av_terr %>% filter(LogMass <3.25) 
nrow(largegen)
nrow(smallgen)
colnames(Genus_terr_eurasia)
#colnames for l
largeM <- largegen$Genus.1.2
smallM <- smallgen$Genus.1.2

#Count large and small genera at each sites (3.14 division)
# save in hid_pts with N_Gen total
Mod_large <- Genus_terr_eurasia[, which((names(Genus_terr_eurasia) %in% largeM)==TRUE)] #filters out genera <1.4kg
colnames(Mod_large)
Mod_large[Mod_large == 0] <- NA
Mod_large <- Mod_large %>% mutate(Large_GenCount=rowSums(!(is.na(.)))) #large genus counts per site
hid_pts$Large_GenCount <- Mod_large$Large_GenCount

Mod_small <- Genus_terr_eurasia[, which((names(Genus_terr_eurasia) %in% smallM)==TRUE)] #filters out genera <1.4kg
colnames(Mod_small)
Mod_small[Mod_small == 0] <- NA
Mod_small <- Mod_small %>% mutate(Small_GenCount=rowSums(!(is.na(.)))) #large genus counts per site
hid_pts$Small_GenCount <- Mod_small$Small_GenCount

#Find ratios
hid_pts <- hid_pts %>% mutate(smallperlarge= Small_GenCount/Large_GenCount)
hid_pts <- hid_pts %>% mutate(smallprop= Small_GenCount/Gen_Count)

#Define herbivore and non genus names
herbgen <- phy_gen_av_terr %>% filter(Plant >=50) 
nonherbgen <- phy_gen_av_terr %>% filter(Plant<50) 
nrow(herbgen)
nrow(nonherbgen)
#colnames for diet groups
herbM <- herbgen$Genus.1.2
nonherbM <- nonherbgen$Genus.1.2

#Count genera by 2 diet types at each sites, save in hid_pts with N_Gen total
Mod_herb <- Genus_terr_eurasia[, which((names(Genus_terr_eurasia) %in% herbM)==TRUE)] #filters out genera <1.4kg
colnames(Mod_herb)
Mod_herb[Mod_herb == 0] <- NA
Mod_herb <- Mod_herb %>% mutate(herb_GenCount=rowSums(!(is.na(.)))) #herb genus counts per site
hid_pts$Herb_GenCount <- Mod_herb$herb_GenCount

Mod_nonherb <- Genus_terr_eurasia[, which((names(Genus_terr_eurasia) %in% nonherbM)==TRUE)] #filters out genera <1.4kg
colnames(Mod_nonherb)
Mod_nonherb[Mod_nonherb == 0] <- NA
Mod_nonherb <- Mod_nonherb %>% mutate(nonherb_GenCount=rowSums(!(is.na(.)))) #herb genus counts per site
hid_pts$Nonherb_GenCount <- Mod_nonherb$nonherb_GenCount

#Find diet ratios
hid_pts <- hid_pts %>% mutate(DietRatio= Herb_GenCount/Nonherb_GenCount)
hid_pts <- hid_pts %>% mutate(Herbprop= Herb_GenCount/Gen_Count)

# looking at only large mammal modern communities, what is diet breakdown?
Mod_largeherb <- Mod_large[, which((names(Mod_large) %in% herbM)==TRUE)] #filters out genera <1.4kg
Mod_largeherb <- Mod_largeherb %>% mutate(largeherb_GenCount=rowSums(!(is.na(.)))) #herb genus counts per site
hid_pts$LargeHerb_GenCount <- Mod_largeherb$largeherb_GenCount

Mod_largenonherb <- Mod_large[, which((names(Mod_large) %in% nonherbM)==TRUE)] #filters out genera <1.4kg
Mod_largenonherb <- Mod_largenonherb %>% mutate(largenonherb_GenCount=rowSums(!(is.na(.)))) #herb genus counts per site
hid_pts$LargeNonherb_GenCount <- Mod_largenonherb$largenonherb_GenCount

#Find ratios
hid_pts <- hid_pts %>% mutate(LargeDietRatio= LargeHerb_GenCount/LargeNonherb_GenCount)
hid_pts <- hid_pts %>% mutate(Large_Herbprop= LargeHerb_GenCount/Large_GenCount)

#looking only at small mammal communities, what is diet breakdown?
Mod_smallherb <- Mod_small[, which((names(Mod_small) %in% herbM)==TRUE)] #filters out genera <1.4kg
Mod_smallherb <- Mod_smallherb %>% mutate(smallherb_GenCount=rowSums(!(is.na(.)))) #herb genus counts per site
hid_pts$SmallHerb_GenCount <- Mod_smallherb$smallherb_GenCount

Mod_smallnonherb <- Mod_small[, which((names(Mod_small) %in% nonherbM)==TRUE)] #filters out genera <1.4kg
Mod_smallnonherb <- Mod_smallnonherb %>% mutate(smallnonherb_GenCount=rowSums(!(is.na(.)))) #herb genus counts per site
hid_pts$SmallNonherb_GenCount <- Mod_smallnonherb$smallnonherb_GenCount

#Find ratios
hid_pts <- hid_pts %>% mutate(smallDietRatio= SmallHerb_GenCount/SmallNonherb_GenCount)
hid_pts <- hid_pts %>% mutate(Small_Herbprop= SmallHerb_GenCount/Small_GenCount)

#How many large mammal communities, and their diet ratios
hid_pts_lpresent <- hid_pts %>% filter(Large_GenCount>0)
summary(hid_pts_lpresent[15:31])
sd(hid_pts_lpresent$LargeNonherb_GenCount)
sd(hid_pts_lpresent$LargeHerb_GenCount)
sd(hid_pts_lpresent$Large_Herbprop)
sd(hid_pts_lpresent$Large_GenCount)

#How many small mammal communities, and their diet ratios
hid_pts_spresent <- hid_pts %>% filter(Small_GenCount>0)
summary(hid_pts_spresent[15:31])
sd(hid_pts_spresent$SmallNonherb_GenCount)
sd(hid_pts_spresent$SmallHerb_GenCount)
sd(hid_pts_spresent$Small_Herbprop)
sd(hid_pts_spresent$Small_GenCount)

# total numbers of each genus type
# find genus columns with eurasian occurrences
Mod_largeherb_occ <- Mod_largeherb[,colSums(is.na(Mod_largeherb))< nrow(Mod_largeherb)]
Mod_largenonherb_occ <- Mod_largenonherb[,colSums(is.na(Mod_largenonherb))< nrow(Mod_largenonherb)]
Mod_smallherb_occ <- Mod_smallherb[,colSums(is.na(Mod_smallherb))< nrow(Mod_smallherb)]
Mod_smallnonherb_occ <- Mod_smallnonherb[,colSums(is.na(Mod_smallnonherb))< nrow(Mod_smallnonherb)]
colnames(Mod_smallnonherb_occ)

# Count how many modern genera are in each guild in total in Eurasia
ncol(Mod_largeherb_occ)-1
ncol(Mod_largenonherb_occ)-1
ncol(Mod_smallherb_occ)-1
ncol(Mod_smallnonherb_occ)-1

#write.csv(hid_pts, "ModernEurasia_DietSizeCounts.csv")
#hid_pts <- read.csv("ModernEurasia_DietSizeCounts.csv", row.names = 1, header=T)

#Diet ratios for modern sites >4 genera
hid_pts_use <- hid_pts %>% filter(Gen_Count>3)
summary(hid_pts_use[,15:31])
sd(hid_pts_use$Nonherb_GenCount)
sd(hid_pts_use$Herb_GenCount)
sd(hid_pts_use$Herbprop)
sd(hid_pts_use$Gen_Count)

# For large mammal communities, where 3 or more large genera (lowest count in fossil large sites)
hid_pts_lpresent3 <- hid_pts_use %>% filter(Large_GenCount>3)
summary(hid_pts_lpresent3[15:31])
sd(hid_pts_lpresent3$LargeNonherb_GenCount)
sd(hid_pts_lpresent3$LargeHerb_GenCount)
sd(hid_pts_lpresent3$Large_Herbprop)
sd(hid_pts_lpresent3$Large_GenCount)
#Same for small mammal communities, where 3 or more small genera
hid_pts_spresent3 <- hid_pts_use %>% filter(Small_GenCount>3)
summary(hid_pts_spresent3[15:31])
sd(hid_pts_spresent3$SmallNonherb_GenCount)
sd(hid_pts_spresent3$SmallHerb_GenCount)
sd(hid_pts_spresent3$Small_Herbprop)
sd(hid_pts_spresent3$Small_GenCount)

# Combine moderns and fossil community guilds counts, to analyze ####
#modern points
colnames(hid_pts_use)
hid_pts_use$Dataset <- rep("ModernEurasia", nrow(hid_pts_use))
Mcounts <- hid_pts_use %>% dplyr::select(c("HID", "LargeHerb_GenCount","LargeNonherb_GenCount",
                                           "SmallHerb_GenCount", "SmallNonherb_GenCount", "Dataset",
                                           "smallprop", "Herbprop", "X", "Y"))

#fossil points-- find combos
SO_use$Dataset <- rep(NA, nrow(SO_use))
SO_use$Dataset[SO_use$NAME %in% largeonly$NAME] <- "LargeSites"
SO_use$Dataset[SO_use$NAME %in% medonly$NAME] <- "MixedSites"
SO_use$Dataset[SO_use$NAME %in% smallonly$NAME] <- "SmallSites"
Fcount <-  SO_use %>% dplyr::select(c("NAME", "LargeHerb_GenCount","LargeNonherb_GenCount",
                                      "SmallHerb_GenCount", "SmallNonherb_GenCount", "Dataset",
                                      "smallprop", "HerbProp", "LAT", "LONG"))
colnames(Fcount)
#match colnames
colnames(Mcounts)[1] <- "NAME"
colnames(Mcounts)[8:10] <- c("HerbProp", "LONG", "LAT")
guilds <- rbind(Mcounts, Fcount)   #combine 4-guild counts for modern and fossil samples
guilds %>% count(Dataset)

#for GLLVM-- downsample modern sites to 8000 pts
Msample <- sample_n(Mcounts, 8000, replace = FALSE)
guildsS <- rbind(Msample, Fcount)

# Plot Fig 6A 
cols <- c("LargeSites"="#0000FF80",  "MixedSites"="#912CEE98",
          "SmallSites" ="#FF000080", "ModernEurasia"="#77889990")

ggplot(guilds)+geom_point(aes(x=smallprop, y=HerbProp, color=Dataset), size=3)+ 
  scale_color_manual(values=cols)+xlab("Proprtion small mammals")+ylab("Proportion herbivores")

#Add convex hull
# sample modern over 0.1 small, under 0.8 herb
Msampmiddle <- Mcounts %>% filter(smallprop>0.1) %>% filter(HerbProp<0.8)
hull_mod <- Msampmiddle %>% slice(chull(smallprop, HerbProp))

ggplot(guilds)+geom_point(aes(x=smallprop, y=HerbProp, color=Dataset), size=3)+ 
  geom_polygon(data = hull_mod, alpha = 0.2, aes(x=smallprop, y=HerbProp, color= Dataset), fill="grey", 
               lwd=1, lty=4)+
  scale_color_manual(values=cols)+xlab("Proportion small mammals")+ylab("Proportion herbivores")

#find proportion of sites within modern hull, including all modern sites
mod_hulln <- convhulln(Mcounts[,7:8])
insidehull <- inhulln(mod_hulln,as.matrix(guilds[,7:8]))  #true are points within hull of p/pn comms
guilds$InsideModhull <- insidehull  #add column

#summaries of site counts inside/outside 
guilds %>% filter(Dataset=="LargeSites") %>% count(InsideModhull)
guilds %>% filter(Dataset=="MixedSites") %>% count(InsideModhull)
guilds %>% filter(Dataset=="SmallSites") %>% count(InsideModhull)
94/483
97/119
31/116

#Fig. 6B and C:
#Geog/climate gradients within modern hull
colnames(hid_pts_use)
tmod <- ggplot(hid_pts)+geom_point(aes(x=smallprop, y=Herbprop, color=BIO01_Mean), size=3)+ 
  geom_polygon(data = hull_mod, alpha = 0.2, aes(x=smallprop, y=HerbProp), fill="grey", 
               lwd=1, lty=4)+
  scale_color_scico(palette="lajolla", name="MAT")+xlab("Proportion small mammals")+ylab("Proportion herbivores")

lmod <- ggplot(hid_pts)+geom_point(aes(x=smallprop, y=Herbprop, color=Y), size=3)+ 
  geom_polygon(data = hull_mod, alpha = 0.2, aes(x=smallprop, y=HerbProp), fill="grey", 
               lwd=1, lty=4)+ ylab(element_blank())+
  scale_color_scico(palette="roma", name="Latitude")+xlab("Proportion small mammals")+ylab("Proportion herbivores")

tmod  %>% insert_right(lmod) #plots

#Correlation plot between guild variables, in modern sample
MODcomp <- hid_pts_use %>% dplyr::select(LargeHerb_GenCount,LargeNonherb_GenCount,
                                         SmallHerb_GenCount, SmallNonherb_GenCount)
col2 <- colorRampPalette(brewer.pal(9,"RdBu"))

correlation_df<-cor(MODcomp, method="pearson", use="pairwise.complete.obs")
corrplot(correlation_df, method = "square", # this is better matrix, indicates few high correlation values
         col=col2(20),tl.col="black",
         addCoef.col = "black", diag=F) #type="lower")

# Run GLLVM
#simple GLLVM for 4 categories: l/s, herb/nonherb
#install.packages("gllvm")
library(gllvm) #error in this function
LVres <- gllvm(y = guildsS[,c(2:5)], family = "negative.binomial", num.lv = 1,   #runs model, with negative binomial fit for variables
               row.eff="random", maxit=1000)  #add random row effect because of spatial structure ; fit only 1 latent variable
guildsS$LV1 <- LVres$lvs[,1]
#write.csv(guilds, "BKV_modern8000fossil_LV.csv")
LVres
plot(LVres)
coef(LVres)

guildsS <- guildsS %>% mutate(LVscale=100000*LV1)  #scale up LV scores
hist(guildsS$LV1)

cols <- c("LargeSites"="#0000FF80",  "MixedSites"="#912CEE98",
          "SmallSites" ="#FF000080", "ModernEurasia"="#778899")

# Plot Fig. 7, ridgelines to show distributions of LV scores by dataset
ggplot(data=lvres) + geom_density_ridges(aes(x=LVscale,y=Dataset,
                                             fill=Dataset), linewidth=0.8) +
  theme_bw() + scale_fill_manual(values=cols)+ xlab("LV1") +ylab("Density")+
  theme(legend.position = "none")

# Test for differences in means and distributions
guildsS %>% count(Dataset)   
# separate datasets
large <- guildsS %>% filter(Dataset=="LargeSites")
mixed <- guildsS %>% filter(Dataset=="MixedSites")
small <- guildsS %>% filter(Dataset=="SmallSites")
fossil <- guildsS %>% filter(Dataset %in% c("LargeSites", "MixedSites", "SmallSites"))
modern <- guildsS %>% filter(Dataset=="ModernEurasia")

#Run t.tests, significance 0.05
t.test(large$LV1, modern$LV1)
t.test(mixed$LV1, modern$LV1) #mixed sites not significantly different from modern
t.test(small$LV1, modern$LV)
t.test(fossil$LV1, modern$LV1) #all other significantly different

#K Run Kolmogorov-Smirnov test for equality of distributions
ks.test(large$LV1, modern$LV1)
ks.test(mixed$LV1, modern$LV1) #mixed sites not significantly different from modern
ks.test(small$LV1, modern$LV1)
ks.test(fossil$LV1, modern$LV1) #all other significantly different
