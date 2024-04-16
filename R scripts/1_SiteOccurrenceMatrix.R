# Code to read in NOW occurrence data, create site/occurrence matrix, and plot genus counts
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

# Load genus mass data
mmass_fin <- read.csv("MammalGenera_MassDiet.csv", header=T, row.names=1)

#Read in all occurrence data (NOW database format)
df_Pleist <- read.csv('NOWdata_supplemented_filtered.csv')

sitenames <- unique(df_Pleist$NAME)  # list of sites
gen <- unique(df_Pleist$GENUS)
gen <- gen[! (gen %in% c("indet.", "gen.", "X", "Indet."))]  #list of genera

#Looking at dietary categories
mmass_occ <- mmass_fin %>% filter(Genus %in% gen)   #find traits of occurring genera
mmass_occ %>% count(MassSource) #how many from each mass data source
mmass_occ %>% count(Diet) #how many herb/nonherb genera
mmass_diet <- mmass_occ %>%  filter(!(is.na(Diet)))
mmass_diet %>% count(DietSource)  #recorded diet sources 
famav <- mmass_diet %>% filter(DietSource=="GroupAv") #which use family averages for diets

# show proportion of occurrences identified to genus/ species level
df_Pleist %>% count(FAMILY=="indet.")
df_Pleist %>% count(GENUS=="indet.")   
df_Pleist %>% count(SPECIES=="indet.") 

#show count of constrained Holocene sites
holocene <- df_Pleist %>% filter(MAX_AGE < 0.0117)
length(unique(holocene$NAME))
# show count of sites with age range crossing Pleistocene/Holocene boundary
phoverlap <- df_Pleist %>% filter(MIN_AGE < 0.0117)
length(unique(phoverlap$NAME))

#plot age_range vs. genus count ####

#summarize genus count
gencounts <- df_Pleist %>% count(NAME, GENUS)  #number of each genus at each site
siterich <- gencounts %>% count(NAME) #number of genera at each site
summary(siterich$n)   #how many genera per site
xy_sites <- unique(df_Pleist[,c(2, 5:6)]) 
siterich_latlong <- merge(siterich, xy_sites, by="NAME", all.x=T, all.y=F)
colnames(siterich_latlong) <- c("NAME", "Gen_Count", "Y", "X")

# Create site/occurrence matrix
mat <- matrix(nrow=length(sitenames), ncol=length(gen))
row.names(mat) <- sitenames
colnames(mat)<- c(gen)

for (p in 1:length(sitenames)){
  for (l in 1:length(gen)){
    s <- sitenames[p]
    g <- gen[l]
    n <- sum(with(df_Pleist,NAME==s & GENUS==g))
    #mat[p,l+1] <- 
    if (n==0){
      mat[p,l] <- 0
    }
    else if (n>=1){
      mat[p,l] <- 1
    }
  }
} 
colnames(mat) <- gen
mat <- cbind(sitenames, mat)
colnames(mat)[1] <- "NAME"

#append locality information to site/occ matrix
v <- rep(NA,(ncol(mat)+13))
for (u in 1:length(sitenames)){
  site <- sitenames[u]
  info <- df_Pleist %>% filter(NAME==site)
  print(u)
  inf <- info[1, c(2, 5, 6, 7, 8, 9, 10, 12, 13, 14, 20, 103)]   #check these, site info columns in df
  drow <- merge(mat,inf,by="NAME")  #this is slow
  v <- rbind(v, drow)
}
v <- v[2:nrow(v),]
SO <- as.data.frame(v)  #this site/occurrence matrix has locality info

colnames(SO)
jdat <- SO[,2:453] #select genus occurrence columns
jdat[jdat == 0] <- NA
jdat <- jdat %>% mutate(GenCount=rowSums(!(is.na(.))))  #count how many occurrences per row
SO$Gen_Count <- jdat$GenCount   #add genus count to df
summary(SO$Gen_Count)
#write.csv(SO, "Pleist_AllEurasia_SiteOccMatrix.csv")

#Alternative site/occurrence matrix, with entries varied if multiple species in genus
mat2 <- matrix(nrow=length(sitenames), ncol=length(gen))
row.names(mat2) <- sitenames
colnames(mat2)<- c(gen)

for (p in 1:length(sitenames)){
  for (l in 1:length(gen)){
    s <- sitenames[p]
    g <- gen[l]
    n <- sum(with(df_Pleist,NAME==s & GENUS==g))
    if (n==0){
      mat2[p,l] <- 0
    }
    else if (n==1){
      mat2[p,l] <- 1
    }
    else if (n==2){
      mat2[p,l] <- 2
    }
    else if (n==3){
      mat2[p,l] <- 3
    }
    else if (n==4){
      mat2[p,l] <- 4
    }
    else if (n>2){
      mat2[p,l] <- 5
    }
  }
} 
colnames(mat2) <- gen
table(unlist(mat2)) #shows how many genus occurrences have 1-5+ species per genus at site

#plot from saved version
SO <- read.csv("Pleist_AllEurasia_SiteOccMatrix.csv", row.names=1)

#Appendix Figure 1
ggplot(data=SO) + geom_point(aes(x=age_range, y=Gen_Count), size=2) + 
  xlab("Site Age Range (Myr)") +ylab("Fossil Genus Count")
# histograms of components
ggplot(data=SO) + geom_histogram(aes(x=age_range))
ggplot()+ geom_histogram(data=SO, aes(x=Gen_Count))+ xlab("Genus count at fossil sites")+
  ylab("Site count")+ theme_bw()
