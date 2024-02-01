# Total code used for figures and analysis in Parker et al. Bjorn Kurten centary volume

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
# Load genus mass data
setwd("C:/Users/abbap/OneDrive - University of Helsinki")
setwd("C:/Users/abigpark/OneDrive - University of Helsinki")
mmass_fin <- read.csv("MammalMassMaster_Jan24_BKVsplit.csv", header=T, row.names=1)

#Read in site data
df_Pleist <- read.csv('NOWdata_Pleist_Sep23_siteedits3.csv')

sites <- unique(df_Pleist$NAME)
gen <- unique(df_Pleist$GENUS)
gen <- gen[! (gen %in% c("indet.", "gen.", "X", "Indet."))]

#Looking at dietary categories
mmass_occ <- mmass_fin %>% filter(Genus %in% gen)
mmass_occ %>% count(MassSource) #how many from each mass source
mmass_occ %>% count(Diet) #how many herb/nonherb
mmass_diet <- mmass_occ %>%  filter(!(is.na(Diet)))
mmass_diet %>% count(GuildDietSource)  #recorded diet sources 
famav <- mmass_diet %>% filter(GuildDietSource=="GroupAv") #which use family averages

# show proportion of occurrences identified to genus/ species level
df_Pleist %>% count(FAMILY=="indet.")
df_Pleist %>% count(GENUS=="indet.")   #2.9% of occurrences without genus id
df_Pleist %>% count(SPECIES=="indet.")  #23.7% of occurrences without species id

#plot age_range vs. genus count ####

#summarize genus count
gencounts <- df_Pleist %>% count(NAME, GENUS)  #number of each genus at each site
siterich <- gencounts %>% count(NAME) #number of genera at each site
colnames(siterich)[1] <- "loc_name"
summary(siterich$n)   #how many genera per site

# create site/occurrence matrix
mat <- matrix(nrow=length(sites), ncol=length(gen))
row.names(mat) <- sites
colnames(mat)<- c(gen)

for (p in 1:length(sites)){
  for (l in 1:length(gen)){
    s <- sites[p]
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
mat <- cbind(sites, mat)
colnames(mat)[1] <- "NAME"

#append locality information to site/occ matrix
v <- rep(NA,(ncol(mat)+13))
for (u in 1:length(sites)){
  site <- sites[u]
  info <- df_Pleist %>% filter(NAME==site)
  print(u)
  inf <- info[1, c(2, 5, 6, 7, 8, 9, 10, 12, 13, 14, 20, 103)]   #check these
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
write.csv(SO, "BKV_AllEurasia_SiteOccMatrix_final1.24.csv")

#Alternative site/occurrence matrix, with entries varied if multiple species in genus
mat2 <- matrix(nrow=length(sites), ncol=length(gen))
row.names(mat2) <- sites
colnames(mat2)<- c(gen)

for (p in 1:length(sites)){
  for (l in 1:length(gen)){
    s <- sites[p]
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
SO <- read.csv("BKV_AllEurasia_SiteOccMatrix_final1.24.csv", row.names=1)

ggplot(data=SO) + geom_point(aes(x=age_range, y=Gen_Count), size=2) + #geom_smooth(aes(x=age_range, y=Gen_Count))
  xlab("Site Age Range (Myr)") +ylab("Fossil Genus Count")
ggplot(data=SO) + geom_histogram(aes(x=age_range))
ggplot()+ geom_histogram(data=SO, aes(x=Gen_Count))+ xlab("Genus count at fossil sites")+
  ylab("Site count")+ theme_bw()

# compare genus count to genus count at closest modern Eco-ISEA3H sample point ####

# read in geography of hexagon ID's ####
setwd("C:/Users/abigpark/OneDrive - University of Helsinki")
setwd("C:/Users/abbap/OneDrive - University of Helsinki")
HIDs <- read.csv("AllTerrestrialHIDs_continents.csv", header=T, row.names=1)   #hexagon info and climate data
colnames(HIDs)
EUHIDs <- HIDs %>% filter(Eurasia=="TRUE") %>% dplyr::select(HID, X, Y, Realm_Mode, Tree_Mean, Elevation_Mean,
                                                             BIO01_Mean, BIO04_Mean, BIO05_Mean,
                                                             BIO06_Mean, BIO12_Mean, 
                                                             BIO15_Mean, NPP, whichNPP)

#read in genus occurrences
setwd("C:/Users/abigpark/OneDrive - University of Helsinki/Phylacine data/GENUS 2024")
setwd("C:/Users/abbap/OneDrive - University of Helsinki/Phylacine data/GENUS 2024")
# GenusMatrix <- read.csv("Phylacine_PresentNatural_AllGenusOccurrences_Res9.csv", header=T, row.names=1)
GenusMatrix <- read.csv("Phylacine_Present_TerrestrialGenusOccurrences_Res9_Jan24.csv", header=T, row.names=1)
GenusCountedMatrix <- read.csv("Phylacine_Present_TerrestrialGenusOccurrences_spcounts_Res9_Jan24.csv", header=T, row.names=1)

#find genus count per point
sum <- rowSums(GenusMatrix[,2:1334])
summary(sum) 
GenusMatrix$N_Gen <- as.numeric(sum)

#Select Eurasian points
Genus_terr_eurasia <- GenusMatrix %>% filter(HID %in% EUHIDs$HID)
Genus_sp_eurasia <- GenusCountedMatrix %>% filter(HID %in% EUHIDs$HID)
summary(Genus_terr_eurasia$N_Gen) #up to 111 co-occurring genera at these sites

f2 <- function(data, levels) { 
  tt <- function(x, levels) tabulate(factor(x, levels), length(levels))
  cc <- vapply(data, tt, integer(2L), levels, USE.NAMES = FALSE)
  res <- as.integer(.rowSums(cc, 2L, length(data))) 
  names(res) <- levels
  res
}
genpre <- f2(Genus_terr_eurasia[,2:1334], c(0,1))

f3 <- function(data, levels) { 
  tt <- function(x, levels) tabulate(factor(x, levels), length(levels))
  cc <- vapply(data, tt, integer(19L), levels, USE.NAMES = FALSE)
  res <- as.integer(.rowSums(cc, 19L, length(data))) 
  names(res) <- levels
  res
}
countedsp <- f3(Genus_sp_eurasia[,2:1334], c(0:18)) #output counts within genera for all occurrences
#check proportions of 2+ species per genus
genpre #total genus occurrence count
countedsp #how many species per point
countedsp[2]/genpre[2] #proportion 1-species occurrences
countedsp[3]/genpre[2] #proportion 2-species in genus
countedsp[4]/genpre[2]  #proportion 3-species in genus
countedsp[5]/genpre[2]  #proportion 4-species in genus

#find closest modern HID for fossil site and compare their genus counts

#make site/occ matrix into sf file
eurasia <- ne_countries(scale="large", continent = c('europe', 'asia'), returnclass = "sf")
desiredCRS <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
sites_sf <- SO %>%
  st_as_sf(
    coords = c("LONG", "LAT"),
    crs = desiredCRS)
site_coords <- do.call(rbind, st_geometry(sites_sf))
site_coords <- cbind(SO$NAME, site_coords)

#make hex ID matrix into sf file
hex_sf <- EUHIDs %>%
  st_as_sf(
    coords = c("X", "Y"),
    crs = desiredCRS)
hex_coords <- do.call(rbind, st_geometry(hex_sf))
hex_coords <- cbind(hex_sf$HID, hex_coords)

#match closest points
neighbours <- RANN::nn2(hex_coords[,-1], site_coords[,-1], 1) # 1 stands for 1 nearest neighbours
closest <- sapply(neighbours, cbind) %>% as_tibble #for each fossil site, closest hex pt as rownumber
sites_sf$Nearest_HID <- rep(NA,1018)
sites_sf$HID_geom <- rep(NA,1018)
colnames(sites_sf)
for (u in 1:1018) {
  rnum= closest[u,1]
  hid <- hex_sf[rnum[[1]],][[1]]  
  coords <- hex_sf[rnum[[1]],][[13]]
  sites_sf[u,465] <- hid
  sites_sf[u,466] <- paste0(coords)   #site_sf now has column for which hex each point is closest to
}

#merge in genus count for modern communities
GenC <- GenusMatrix[,c(1,1335)] #select HID and genus count
colnames(GenC)[2] <- "Modern_GenCount"
sites <- merge(sites_sf, GenC, by.x= "Nearest_HID", by.y="HID")  #merge with fossil sites mammal data
# and save for modern points
hid_pts <- cbind(EUHIDs, Genus_terr_eurasia[,1335]) #Add Genus count to hex ids
colnames(hid_pts)[15] <-  "Gen_Count"

setwd("C:/Users/abigpark/OneDrive - University of Helsinki") 
setwd("C:/Users/abbap/OneDrive - University of Helsinki") 
sites_2 <- read.csv("BKV_SiteOccs_PresentComparison.csv", row.names=1) #X and Y for site long/lat
colnames(sites_2)
sites <- merge(sites, sites_2[,c(2,499:500)], by="NAME") #save geometry column for sf plots
#write.csv(sites,"BKV_SiteOccMatrix_final_ModernGenusCount.csv", row.names = F) #save to read in later


#merge in modern MAT and MAP (to do)

#plot fossil vs. nearest modern genus count
ggplot(sites)+ geom_point(aes(x=Gen_Count, y= Modern_GenCount, color=Y), size=2) + 
  scale_color_scico(palette="turku", name="Latitude")+ xlab("Fossil Site Genus Count")+ 
  ylab ("Nearest Modern Sample Genus Count")+ xlim(c(0,55)) + ylim(c(0,115))+
  geom_abline(intercept = c(0,0), slope = 1, color="slategrey")+
 scale_x_continuous(breaks=c(0,10,20,30,40,50))

ggplot(sites)+ geom_point(aes(x=Gen_Count, y= Modern_GenCount, color=X), size=3) + 
  scale_color_scico(palette="turku", name="Longitude")+ xlab("Fossil Genus Count")+ 
  ylab ("Nearest Modern Sample Genus Count")+ xlim(c(0,60)) + ylim(c(0,150))

#Plotting geographic patterns of genus count ####
#genus count vs. latitude
# for fossil occurrences
fossilgenlat <- ggplot(sites)+ geom_point(aes(x=Gen_Count, y= Y, color=X), size=2) + 
  scale_color_scico(palette="corkO", name="Longitude")+
  ylab ("Latitude")+ xlim(c(0,115)) +xlab(element_blank())
# for modern nearest points
ggplot(sites)+ geom_point(aes(x=Modern_GenCount, y= Y, color=X), size=3) + 
 scale_color_scico(palette="corkO", name="Long")+
  ylab ("Latitude")+ xlim(c(0,115))
# for modern all eurasia
Eumat <- hid_pts %>% filter(X>-10.1) #remove points past longitude 180, for consistent color scale
moderngenlat <- ggplot(Eumat)+ geom_point(aes(x=Gen_Count, y= Y, color=X), size=2) + 
  scale_color_scico(palette="corkO", name="Longitude")+
  ylab ("Latitude")+ xlim(c(0,115)) +xlab("Number of Genera")

fossilgenlat %>% 
  insert_bottom(moderngenlat)

setwd("C:/Users/abbap/OneDrive - University of Helsinki/Figures/BK vol biases")
#ggsave("NGenLat2.eps", width = 10, height = 14, units = "cm")

#genus count vs. longitude
# for fossil occurrences
fossilgenlong <- ggplot(sites)+ geom_point(aes(x=X, y= Gen_Count, color=Y), size=2) + 
  scale_color_scico(palette="buda", name="Latitude")+ xlab(element_blank())+
  ylab ("Fossil Site Genus Count")
#for modern nearest points
ggplot(sites)+ geom_point(aes(x=X, y=Modern_GenCount), size=3) + 
  # scale_color_scico(palette="corkO", name="Min Site Age")+
  ylab ("Latitude")+ xlim(c(0,180))
# for modern all eurasia
moderngenlong<- ggplot(hid_pts)+ geom_point(aes(x=X, y= Gen_Count, color=Y), size=2) + 
  scale_color_scico(palette="buda", name="Latitude")+
  ylab ("Modern Sample Genus Count")+ xlim(c(0,180))

fossilgenlong %>% 
  insert_bottom(moderngenlong)

#3-part graph: map of sites with genus count
#define color scale
palette_gen <- colorRampPalette(colors = c("lightskyblue1",
                                           "darkorchid4", "midnightblue"))(20)
scales::show_col(palette_gen)
palette_gen #use first half of color scale for less-genus rich fossil plot, full scale for modern


sitemap <- ggplot()+ geom_sf(data=eurasia)+  #plot genus count of RS output
  geom_sf(data=sites, aes(color=Gen_Count), shape=18, size=3)+
  scale_color_gradient(low="#B0E2FF", high="#6B2C91")+
  coord_sf(xlim=c(-12, 170), ylim= c(-10,80), expand=F)+
 # guides(color=guide_legend(title="Mammal Genus Count")) +
  #theme(legend.position="left")+ 
  #guides(shape = guide_legend(override.aes = list(size = 0.3)))
  theme(legend.position="none")

#genus count in 5-degree longitude bins
binsum_long <- sites %>%    #finding mean values by bin to color histogram by them
  mutate(bins = cut(x= sites$X, 
                    breaks = hist(sites$X, seq(-10,170, by=5), plot=F)$breaks, 
                    labels = hist(sites$X, seq(-10,170, by=5), plot=F)$mids)) %>%
  group_by(bins) %>%
  summarize(n = n(), mean = mean(Gen_Count), breaks=sites$breaks) 

#binsum_long <- binsum_long %>% mutate(bin_num= as.numeric(as.character((bins)))) %>% 
 # mutate(bin_low= round_any(bin_num, 5, floor))  # shows minimum long for each bin (x axis)

plt2 <- ggplot(binsum_long)+geom_point(aes(x=n, y=mean, color=mean), size=4)+ #define color scale for mean genus count
  scale_color_gradient(low="#B0E2FF", high="#6B2C91")
ld2 <- layer_data(plt2)
my_colA <- ld2$colour
meancols <- c(my_colA[1:28], "white", "white", my_colA[29], "white", my_colA[30:31], "white", "white")

longB <- ggplot(sites)+ geom_histogram(aes(x=X), breaks = seq(-10, 170, 5),
                              color="black", fill=meancols)+
  scale_fill_gradient(low="#B0E2FF", high="#6B2C91") +
  xlab(element_blank())+ ylab("n")

# genus count in 5 degree latitude bins
binsum_lat <- sites %>%    #finding mean values by bin to color histogram by them
  mutate(bins = cut(x= sites$Y, 
                    breaks = hist(sites$Y, seq(-10,80, by=5), plot=F)$breaks,
                    labels = hist(sites$Y, seq(-10,80, by=5), plot=F)$mids)) %>%
  group_by(bins) %>%
  summarize(n = n(), mean = mean(Gen_Count)) 

plt1 <- ggplot(binsum_lat)+geom_point(aes(x=n, y=mean, color=mean), size=4)+ #define color scale for mean genus count
  scale_color_gradient(low="#B0E2FF", high="#6B2C91")
ld1 <- layer_data(plt1)
my_colB <- ld1$colour

latB <- ggplot(sites)+ geom_histogram(aes(x=Y), breaks = seq(-10, 80, 5),
                                       color="black", fill=my_colB)+
  scale_fill_gradient(low="#B0E2FF", high="#6B2C91") +
  xlab(element_blank())+ ylab("n") + coord_flip()
latB

# alternate bar graph: summarizing average genus count by 5 degree long bins
longB_mean <- ggplot() + 
  stat_summary_bin(data= sites, mapping=aes(x=X, y=Gen_Count),
                   fun = "mean", geom="bar", breaks= seq(-10,170, by=5),
                   color= "black", fill=my_colA) +
  scale_fill_gradient(low="#B0E2FF", high="#6B2C91")+
  ylab("Mean Genus Count") + xlab("Longitude Bins")+
  xlim(-10,170)
longB_mean

# same for lat bins
latB_mean <- ggplot() + 
  stat_summary_bin(data= sites, mapping=aes(x=Y, y=Gen_Count),
                   fun = "mean", geom="bar", breaks= seq(-10,180, by=5),
                   color= "black", fill=my_colB) +
  scale_fill_gradient("#B0E2FF", high="#6B2C91")+
  ylab("Mean Genus Count") + xlab("Latitude Bins")+
  xlim(-10, 80) + coord_flip()
#latB_mean


#combine plots-- aligned axes by hand here
ggplot()+
  coord_equal(xlim = c(-4, 50), ylim = c(0, 45), expand = FALSE) +
  annotation_custom(ggplotGrob(sitemap), xmin = 0, xmax = 39, ymin = 3, 
                    ymax = 27) +
  annotation_custom(ggplotGrob(longB), xmin = -3, xmax = 40, ymin = 27, 
                    ymax = 34)+
  annotation_custom(ggplotGrob(latB), xmin = 39, xmax = 46, ymin = 1, 
                    ymax = 28) +
 theme(panel.spacing = unit(0, "lines"))

#same for modern
hex_sf <- cbind(hex_sf, hid_pts[,c(2:3,15)]) #Add Genus count to hex ids
#colnames(hex_sf)[15] <- "Gen_Count"
modmap <- ggplot()+ geom_sf(data=eurasia)+  #plot genus count of RS output
  geom_sf(data=hex_sf, aes(color=Gen_Count), shape=18, size=1.8)+
  scale_color_gradientn(colors=palette_gen)+
  coord_sf(xlim=c(-12, 170), ylim= c(-10,80), expand=F)+
  theme(legend.position="none")
modmap

#genus count in 5-degree longitude bins
hex_long <- hex_sf %>% filter(X > -10.01) %>% filter(X < 170)
Mbinsum_long <- hex_long %>%    #finding mean values by bin to color histogram by them
  mutate(bins = cut(x= hex_long$X, 
                    breaks = hist(hex_long$X, seq(-10,170, by=5), plot=F)$breaks, 
                    labels = hist(hex_long$X, seq(-10,170, by=5), plot=F)$mids)) %>%
  group_by(bins) %>%
  summarize(n = n(), mean = mean(Gen_Count), breaks=hex_sf$breaks) 

plt3 <- ggplot(Mbinsum_long)+geom_point(aes(x=n, y=mean, color=mean), size=4)+ #define color scale for mean genus count
  scale_color_gradientn(colors=palette_gen)
ld3 <- layer_data(plt3)
my_colC <- ld3$colour

MlongB <- ggplot(hex_long)+ geom_histogram(aes(x=X), breaks = seq(-10, 170, 5),
                                       color="black", fill=my_colC)+
  scale_fill_gradientn(colors=palette_gen) +
  xlab(element_blank())+ ylab("n")
MlongB

# genus count in 5 degree latitude bins
hex_lat <- hex_sf %>% filter(Y > -10.001) %>% filter(Y <80.01)  #for consistent plotting range
Mbinsum_lat <- hex_lat %>%    #finding mean values by bin to color histogram by them
  mutate(bins = cut(x= hex_lat$Y, 
                    breaks = hist(hex_lat$Y, seq(-10,80, by=5), plot=F)$breaks,
                    labels = hist(hex_lat$Y, seq(-10,80, by=5), plot=F)$mids)) %>%
  group_by(bins) %>%
  summarize(n = n(), mean = mean(Gen_Count)) 

plt4 <- ggplot(Mbinsum_lat)+geom_point(aes(x=n, y=mean, color=mean), size=4)+ #define color scale for mean genus count
  scale_color_gradientn(colors=palette_gen)
ld4 <- layer_data(plt4)
my_colD <- ld4$colour

MlatB <- ggplot(hex_long)+ geom_histogram(aes(x=Y), breaks = seq(-10, 80, 5),
                                      color="black", fill=my_colD)+
  scale_fill_gradientn(colors=palette_gen) +
  xlab(element_blank())+ ylab("n") + coord_flip()
MlatB

#combined
ggplot()+
  coord_equal(xlim = c(-4, 50), ylim = c(0, 45), expand = FALSE) +
  annotation_custom(ggplotGrob(modmap), xmin = 0, xmax = 39, ymin = 3, 
                    ymax = 27) +
  annotation_custom(ggplotGrob(MlongB), xmin = -3, xmax = 40, ymin = 27, 
                    ymax = 34)+
  annotation_custom(ggplotGrob(MlatB), xmin = 39, xmax = 46, ymin = 1, 
                    ymax = 28) +
  theme(panel.spacing = unit(0, "lines"))




#now, compare to climate



# Divide small/large mammal genera ####

EM <- matrix(nrow=1,ncol=(1+max(SO$Gen_Count))) #creating empty data frame, for HID plus up to 198 species names
MassFrame <- as.data.frame(EM)

colnames(mmass_fin)  #logmass of genus is col 27

jdat <- SO[,2:453] #select genus occurrence columns of site/occurrence matrix
jdat[jdat == 0] <- NA

for(k in 1:nrow(SO)){   #for each locality
  V <- rep(NA, len=(1+max(SO$Gen_Count)))      #empty vector with length of max number of genera +1
  V[1] <- SO[k,1] #loc_name as first entry
  indx <- which(!(is.na(jdat[k,])), arr.ind=TRUE)  #which columns in row k have 1 
  genu <- colnames(jdat)[indx[,2]]     #list of genera at point
  print(genu[1])
  for (p in (1:length(genu))) {
    genus <- as.character(genu[p])  #pulls genus name
    ma <- mmass_fin[mmass_fin$Genus==genus,27]  #find logmass
    V[p+1] <- ma[1]
  }
  MassFrame <- rbind(MassFrame, V)
}

InputMass <- MassFrame[2:nrow(MassFrame),]
colnames(InputMass)[1] <- "loc_name"
write.csv(InputMass,"BKV_AllEurasia_InputMassDist_1.24.csv")
InputMass <- read.csv("BKV_AllEurasia_InputMassDist_1.24.csv", row.names = 1)

# Find distribution statistics for masses at each input site
clim_pts <- SO %>% dplyr::select(NAME, LAT, LONG, ALTITUDE, MAX_AGE, MIN_AGE, BFA_MAX, BFA_MAX_ABS, BFA_MIN, BFA_MIN_ABS,
                                 COUNTRY, age_range, Gen_Count) %>% 
                              mutate(mid_age= ((MAX_AGE+MIN_AGE)/2))

clim_pts$mean <- rowMeans(InputMass[,2:ncol(InputMass)], na.rm=T)

sd <- rep(NA, nrow(InputMass))
for (p in 1:nrow(InputMass)){      
  sd[p] <- sd(InputMass[p,2:ncol(InputMass)],na.rm=T)
}
clim_pts$SD <- sd

mx <- rep(NA, nrow(InputMass))
for (p in 1:nrow(InputMass)){            
  mx[p] <- max(InputMass[p,2:ncol(InputMass)],na.rm=T)
}
mx[is.infinite(mx)] <- NA
clim_pts$max <- mx

mn <- rep(NA, nrow(InputMass))
for (p in 1:nrow(InputMass)){            
  mn[p] <- min(InputMass[p,2:ncol(InputMass)],na.rm=T)
}
mn[is.infinite(mn)] <- NA
clim_pts$min <- mn

# Look at large/small mammal site divisions
ggplot()+ geom_histogram(data=clim_pts , aes(x=mean))+ xlab("Site Mean of Log(Body Mass g)")+
  ylab("Site count") +
  ggtitle("Histogram of mean body masses for all NOW Pleistocene sites")
ggsave("AllSiteMeanHistogram.eps", width = 7.5, units = "cm")

#2 or more species for Jan/May 23 runs
use_subset <- clim_pts %>% filter(Gen_Count >3)

ggplot()+ geom_histogram(data=use_subset , aes(x=mean))+ xlab("Mean of Log(Body Mass g)")+
  ylab("Site count") +
  ggtitle("Histogram of mean body masses \nfor NOW Pleistocene sites with 4 or more gen.")

# From these genera occurring in NOW, look at mass distribution
avfree <- mmass_fin %>% filter(!(is.na(LogMass))) %>%  filter(!(MassSource=="Family average")) %>% 
  filter(Genus %in% gen)  #genera occurring in NOW, excluding those with family average masses

ggplot()+ geom_histogram(data=avfree, aes(x=LogMass), binwidth=0.2)

#Bimodal distribution: fit model
model <- normalmixEM(avfree$LogMass,
                     mu = c(1, 5), 
                     sigma = 1,
                     maxit = 200)
summary(model)

# function to plot two distributions on histogram
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

mixmdl.df = data.frame(x = model$x)

f.plot = ggplot(mixmdl.df) +
  geom_histogram(aes(x, ..density..),
                 binwidth = 0.2, colour = "black",
                 fill = "azure2",
                 alpha = 0.9
  ) +
  stat_function(
    geom = "line", fun = plot_mix_comps,
    args = list(model$mu[1], model$sigma[1], lam = model$lambda[1]),
    colour = "#FF000080", lwd = 1.5
  ) +
  stat_function(
    geom = "line", fun = plot_mix_comps,
    args = list(model$mu[2], model$sigma[2], lam = model$lambda[2]),
    colour = "#0000FF80", lwd = 1.5
  ) +
  ylab(" Fossil \nGenus Density") + xlab(element_blank())

#extract coordinates of both curves from plot
line.df = data.frame(x = ggplot_build(f.plot)$data[[2]]$x,
                     red = ggplot_build(f.plot)$data[[2]]$y,
                     blue = ggplot_build(f.plot)$data[[3]]$y)

#find the minimal distance between lines along y axis
line.df$delta = line.df$red - line.df$blue

#find x coordinate for the minimal delta y
f.x_coord = line.df$x[which(diff(sign(diff((abs(line.df$delta))))) == 2)+1]
f.x_coord  #this is cutoff between two groups
10^f.x_coord # kg value of cutoff

f.plot <- f.plot+ geom_vline(xintercept = f.x_coord, lty=2)

#If two distributions have same lambda, find crossover point
avfree$Small <- pnorm(avfree$LogMass, mean=model$mu[1], sd=model$sigma[1], lower.tail = F)
avfree$Large <- pnorm(avfree$LogMass, mean=model$mu[2], sd=model$sigma[2], lower.tail = T)
#Find where likelihoods are similar (green and red lines cross)
avfree <- avfree %>% mutate(Diff= Large - Small)
which.min(abs(avfree$Diff))
avfree[141,27]  #find cutoff among genus masses
10^3.41

#Modern mammals instead (Phylacine database)

#first by species ####
setwd("C:/Users/abbap/OneDrive - University of Helsinki/Phylacine data")
setwd("C:/Users/abigpark/OneDrive - University of Helsinki/Phylacine data")
phy <- read.csv("Trait_data.csv", header=T)
# removed species with phylogenetically imputed masses
phy <- phy %>% mutate(LogMass= log10(Mass.g)) %>% filter(!(Mass.Method == "Imputed")) #%>% filter(Aerial==0) 
#remove aquatic species 
marine <- phy %>% filter(Marine==1)
marine <- marine %>% filter(!(Family.1.2 %in% c("Hippopotamidae", "Ursidae","Muridae","Canidae", 
                                                "Otariidae", "Mustelidae")))
# add freshwater fully-aquatic animals
fwat <- phy %>% filter(Family.1.2 %in% c("Delphiniidae", "Iniidae", "Platanistidae",
                                            "Phocoenidae", "Trichechidae"))
aq <- rbind(marine, fwat) #add marine, freshwater families, plus these 3 awquatic species
aq_sp <- rbind(marine, fwat, "Enhydra_lutris", "Lontra_felina", "Neovison_macrodon")
#save species names
aqs_sp <- aq_sp$Binomial.1.2
phy_terr <- phy %>% filter(!(Binomial.1.2 %in% aqs_sp)) #%>% filter(Aerial==0) 
#save genus names 
aq_gen <- rbind(marine, fwat, "Enhydra","Neovison") #Lontra more terrestrial species
aq_gen <- aq_gen$Genus.1.2

#without imputed species masses
ggplot()+ geom_histogram(data=phy_terr, aes(x=LogMass))+ xlab("Log(Body Mass g)")+
  ylab("Species count") +
  ggtitle("Histogram of mean body masses \nfor modern terrestrial mammal species")
modelp <- normalmixEM(phy_terr$LogMass, 
                       mu = c(1, 5), 
                       sigma = 1, 
                       maxit=200)

summary(modelp)

# function to plot two distributions on histogram
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

mixmdl.df = data.frame(x = modelp$x)

p.plot = ggplot(mixmdl.df) +
  geom_histogram(aes(x, ..density..),
                 binwidth = 0.2, colour = "black",
                 fill = "azure2",
                 alpha = 0.9
  ) +
  stat_function(
    geom = "line", fun = plot_mix_comps,
    args = list(modelp$mu[1], modelp$sigma[1], lam = modelp$lambda[1]),
    colour = "#FF000080", lwd = 1.5
  ) +
  stat_function(
    geom = "line", fun = plot_mix_comps,
    args = list(modelp$mu[2], modelp$sigma[2], lam = modelp$lambda[2]),
    colour = "#0000FF80", lwd = 1.5
  ) +
  ylab("Modern Global \nSpecies Density") + xlab(element_blank())

#extract coordinates of both curves from plot
line.df = data.frame(x = ggplot_build(p.plot)$data[[2]]$x,
                     red = ggplot_build(p.plot)$data[[2]]$y,
                     blue = ggplot_build(p.plot)$data[[3]]$y)

#find the minimal distance between lines along y axis
line.df$delta = line.df$red - line.df$blue

#find x coordinate for the minimal delta y
p.x_coord = line.df$x[which(diff(sign(diff((abs(line.df$delta))))) == 2)+1]
p.x_coord  #this is cutoff between two groups
10^p.x_coord # kg value of cutoff

p.plot <- p.plot+ geom_vline(xintercept = p.x_coord, lty=2)

#assign likelihood of each genus belonging to large or small group, with equal lambdas 
phy_terr$Small <- pnorm(phy_terr$LogMass, mean=modelp$mu[1], sd=modelp$sigma[1], lower.tail = F)
phy_terr$Large <- pnorm(phy_terr$LogMass, mean=modelp$mu[2], sd=modelp$sigma[2], lower.tail = T)
#Find where likelihoods are similar (green and red lines cross)
phy_terr <- phy_terr %>% mutate(Diff= Large - Small)
which.min(abs(phy_terr$Diff))  #most similar probabilities of large/small
phy_terr[3645,12] #kg value of crossover
phy_terr[3645,25] # log value of crossover

#then by genus ####
phy_gen_av <- phy %>% group_by(Genus.1.2) %>% summarize(Plant = mean(Diet.Plant, na.rm = TRUE),  #genus average %plant
                LogMass = mean(log10(Mass.g),na.rm=T)) #save genus average mass (not including imputed species)

phy_gen_av_terr <- phy_gen_av %>% filter(!(Genus.1.2 %in% aq_gen)) #remove aquatic genera

ggplot()+ geom_histogram(data=phy_gen_av_terr, aes(x=LogMass))+ xlab("Log(Body Mass g)")+
  ylab("Species count") +
  ggtitle("Histogram of mean body masses \nfor modern mammal genera")

modelpg <- normalmixEM(phy_gen_av_terr$LogMass,
                       mu = c(1, 5), 
                       sigma= 1,
                       maxit=200)

summary(modelpg)

# function to plot two distributions on histogram
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}
mixmdl.df = data.frame(x = modelpg$x)

pg.plot = ggplot(mixmdl.df) +
  geom_histogram(aes(x, ..density..),
                 binwidth = 0.2, colour = "black",
                 fill = "azure2",
                 alpha = 0.9
  ) +
  stat_function(
    geom = "line", fun = plot_mix_comps,
    args = list(modelpg$mu[1], modelpg$sigma[1], lam = modelpg$lambda[1]),
    colour = "#FF000080", lwd = 1.5
  ) +
  stat_function(
    geom = "line", fun = plot_mix_comps,
    args = list(modelpg$mu[2], modelpg$sigma[2], lam = modelpg$lambda[2]),
    colour = "#0000FF80", lwd = 1.5
  ) +
  ylab(" Modern Global \nGenus Density") + xlab(element_blank()) 

#extract coordinates of both curves from plot
line.df = data.frame(x = ggplot_build(pg.plot)$data[[2]]$x,
                     red = ggplot_build(pg.plot)$data[[2]]$y,
                     blue = ggplot_build(pg.plot)$data[[3]]$y)

#find the minimal distance between lines along y axis
line.df$delta = line.df$red - line.df$blue

#find x coordinate for the minimal delta y
pg.x_coord = line.df$x[which(diff(sign(diff((abs(line.df$delta))))) == 2)+1]
pg.x_coord  #this is cutoff between two groups
10^pg.x_coord # kg value of cutoff

pg.plot <- pg.plot+ geom_vline(xintercept = pg.x_coord, lty=2)

f.plot %>% 
  insert_bottom(p.plot) %>% 
  insert_bottom(pg.plot) 

#assign likelihood of each genus belonging to large or small group (done in 12_12 MammalMassMaster)
phy_gen_av_terr$Small <- pnorm(phy_gen_av_terr$LogMass, mean=modelpg$mu[1], sd=modelpg$sigma[1], lower.tail = F)
phy_gen_av_terr$Large <- pnorm(phy_gen_av_terr$LogMass, mean=modelpg$mu[2], sd=modelpg$sigma[2], lower.tail = T)
#Find where likelihoods are similar (green and red lines cross)
phy_gen_av_terr <- phy_gen_av_terr %>% mutate(Diff= Large - Small)
which.min(abs(phy_gen_av_terr$Diff))  #most similar probabilities of large/small
phy_gen_av_terr[966,3] # log value of crossover
10^3.28

# Phylacine genera in Eurasia ####
#genera occurring in Eurasia 
colnames(Genus_terr_eurasia)[1:10]
colnames(Genus_terr_eurasia)[1330:1335]
EUmat_occs <- Genus_terr_eurasia[,2:1334]
Eurasia_occ <- EUmat_occs[which (colSums(EUmat_occs)>0)] #genus columns with occurrences
Eurasia_genera <- colnames(Eurasia_occ)

gen_av_eurasia <- phy_gen_av %>% filter(Genus.1.2 %in% Eurasia_genera)

ggplot()+ geom_histogram(data=gen_av_eurasia, aes(x=LogMass))+ xlab("Log(Body Mass g)")+
  ylab("Genus count") +
  ggtitle("Histogram of mean body masses \nfor modern mammal genera")

modelpge <- normalmixEM(gen_av_eurasia$LogMass,
                       mu = c(1, 5), 
                       sigma= 1,
                       maxit=200)

summary(modelpge)

# function to plot two distributions on histogram
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}
mixmdl.df = data.frame(x = modelpge$x)

pge.plot = ggplot(mixmdl.df) +
  geom_histogram(aes(x, ..density..),
                 binwidth = 0.2, colour = "black",
                 fill = "azure2",
                 alpha = 0.9
  ) +
  stat_function(
    geom = "line", fun = plot_mix_comps,
    args = list(modelpge$mu[1], modelpge$sigma[1], lam = modelpge$lambda[1]),
    colour = "#FF000080", lwd = 1.5
  ) +
  stat_function(
    geom = "line", fun = plot_mix_comps,
    args = list(modelpge$mu[2], modelpge$sigma[2], lam = modelpge$lambda[2]),
    colour = "#0000FF80", lwd = 1.5
  ) +
  ylab("Modern Eurasian \nGenus Density")+ xlab("Log Mass (kg)") 

#extract coordinates of both curves from plot
line.df = data.frame(x = ggplot_build(pge.plot)$data[[2]]$x,
                     red = ggplot_build(pge.plot)$data[[2]]$y,
                     blue = ggplot_build(pge.plot)$data[[3]]$y)

#find the minimal distance between lines along y axis
line.df$delta = line.df$red - line.df$blue

#find x coordinate for the minimal delta y
pge.x_coord = line.df$x[which(diff(sign(diff((abs(line.df$delta))))) == 2)+1]
pge.x_coord  #this is cutoff between two groups
10^pge.x_coord # kg value of cutoff

pge.plot <- pge.plot+ geom_vline(xintercept = pge.x_coord, lty=2)

f.plot %>% 
  insert_bottom(p.plot) %>% 
  insert_bottom(pg.plot) %>% 
  insert_bottom(pge.plot)

#assign likelihood of each genus belonging to large or small group (done in 12_12 MammalMassMaster)
gen_av_eurasia$Small <- pnorm(gen_av_eurasia$LogMass, mean=modelpge$mu[1], sd=modelpge$sigma[1], lower.tail = F)
gen_av_eurasia$Large <- pnorm(gen_av_eurasia$LogMass, mean=modelpge$mu[2], sd=modelpge$sigma[2], lower.tail = T)
#Find where likelihoods are similar (green and red lines cross)
gen_av_eurasia <- gen_av_eurasia %>% mutate(Diff= Large - Small)
which.min(abs(gen_av_eurasia$Diff))  #most similar probabilities of large/small
gen_av_eurasia[303,3] # log value of crossover
10^3.03


#Eurasia species
setwd("C:/Users/abigpark/OneDrive - University of Helsinki/Phylacine Data/Occurrence Matrices")
setwd("C:/Users/abbap/OneDrive - University of Helsinki/Phylacine Data/Occurrence Matrices")
A1samples <- read.csv("PA1_SpeciesList_GlobalHex9Terrestrial_small_nonherb_bats.csv", header=T)
A1EuraSpList <- A1samples %>%  filter(HID %in% EUHIDs$HID)
colnames(A1EuraSpList)
A1eur_spec <- A1EuraSpList[,3:115]
A1eur_species <- unique(unlist(A1eur_spec))

A2EuraSpList <- read.csv("PA2_Eurasia_SpeciesList_GlobalHex9Terrestrial_small_herb_bats.csv", header=T)
colnames(A2EuraSpList) 
A2eur_spec <- A2EuraSpList[,3:64]
A2eur_species <- unique(unlist(A2eur_spec))

#exclude these (option to check excluding bats, instead of A1 and A2)
setwd("C:/Users/abigpark/OneDrive - University of Helsinki/Phylacine Data/Occurrence Matrices/New 2023")
setwd("C:/Users/abbap/OneDrive - University of Helsinki/Phylacine Data/Occurrence Matrices/New 2023")
A3samples <- read.csv("PA3_SpeciesList_GlobalHex9Terrestrial_small_nonherb_nobats.csv", header=T)
A3EuraSpList <- A3samples %>%  filter(HID %in% EUHIDs$HID)
A3eur_spec <- A3EuraSpList[,3:52]
A3eur_species <- unique(unlist(A3eur_spec))

A4samples <- read.csv("PA4_SpeciesList_GlobalHex9Terrestrial_small_herb_nobats.csv", header=T)
colnames(A4samples)
A4EuraSpList <- A4samples %>%  filter(HID %in% EUHIDs$HID)
A4eur_spec <- A4EuraSpList[,3:56]
A4eur_species <- unique(unlist(A4eur_spec))

setwd("C:/Users/abbap/OneDrive - University of Helsinki/Phylacine Data/Occurrence Matrices/New 2023")
setwd("C:/Users/abigpark/OneDrive - University of Helsinki/Phylacine Data/Occurrence Matrices/New 2023")
A5samples <- read.csv("PA5_SpeciesList_GlobalHex9Terrestrial_medium_nonherb_nobats.csv", header=T)
colnames(A5samples)
A5EuraSpList <- A5samples %>%  filter(HID %in% EUHIDs$HID)
A5eur_spec <- A5EuraSpList[,3:29]
A5eur_species <- unique(unlist(A5eur_spec))

A6samples <- read.csv("PA6_SpeciesList_GlobalHex9Terrestrial_medium_herb_nobats.csv", header=T)
colnames(A6samples)
A6EuraSpList <- A6samples %>%  filter(HID %in% EUHIDs$HID)
A6eur_spec <- A6EuraSpList[,3:37]
A6eur_species <- unique(unlist(A6eur_spec))

A7samples <- read.csv("PA7_SpeciesList_GlobalHex9Terrestrial_large_nonherb_nobats.csv", header=T)
colnames(A7samples)
A7EuraSpList <- A7samples %>%  filter(HID %in% EUHIDs$HID)
A7eur_spec <- A7EuraSpList[,3:7]
A7eur_species <- unique(unlist(A7eur_spec))

A8samples <- read.csv("PA8_SpeciesList_GlobalHex9Terrestrial_large_herb_nobats.csv", header=T)
colnames(A8samples)
A8EuraSpList <- A8samples %>%  filter(HID %in% EUHIDs$HID)
A8eur_spec <- A8EuraSpList[,3:27]
A8eur_species <- unique(unlist(A8eur_spec))

Eurasia_species <- c(A1eur_species, A2eur_species, A5eur_species, A6eur_species,
                     A7eur_species, A8eur_species)
Eurasia_species <- unique(Eurasia_species)
Eurasia_species <- Eurasia_species[!is.na(Eurasia_species)]  #identified 760 species occurring in present Eurasia

phy_eurasia <- phy_terr %>% filter(Binomial.1.2 %in% Eurasia_species) #only species in Eurasia, 1877
#without imputed species masses
ggplot()+ geom_histogram(data=phy_eurasia, aes(x=LogMass))+ xlab("Log(Body Mass g)")+
  ylab("Species count") +
  ggtitle("Histogram of mean body masses \nfor modern mammal species")
modelpe <- normalmixEM(phy_eurasia$LogMass, 
                      mu = c(1, 5), 
                      sigma = 1, 
                      maxit=200)

summary(modelpe)

# function to plot two distributions on histogram
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

mixmdl.df = data.frame(x = modelpe$x)

pe.plot = ggplot(mixmdl.df) +
  geom_histogram(aes(x, ..density..),
                 binwidth = 0.2, colour = "black",
                 fill = "azure2",
                 alpha = 0.9
  ) +
  stat_function(
    geom = "line", fun = plot_mix_comps,
    args = list(modelpe$mu[1], modelpe$sigma[1], lam = modelpe$lambda[1]),
    colour = "#FF000080", lwd = 1.5
  ) +
  stat_function(
    geom = "line", fun = plot_mix_comps,
    args = list(modelpe$mu[2], modelpe$sigma[2], lam = modelpe$lambda[2]),
    colour = "#0000FF80", lwd = 1.5
  ) +
  ylab("Modern Eurasian \nSpecies Density") + xlab(element_blank())

#extract coordinates of both curves from plot
line.df = data.frame(x = ggplot_build(pe.plot)$data[[2]]$x,
                     red = ggplot_build(pe.plot)$data[[2]]$y,
                     blue = ggplot_build(pe.plot)$data[[3]]$y)

#find the minimal distance between lines along y axis
line.df$delta = line.df$red - line.df$blue

#find x coordinate for the minimal delta y
pe.x_coord = line.df$x[which(diff(sign(diff((abs(line.df$delta))))) == 2)+1]
pe.x_coord  #this is cutoff between two groups
10^pe.x_coord # kg value of cutoff

pe.plot <- pe.plot+ geom_vline(xintercept = pe.x_coord, lty=2) #+
pe.plot

#assign likelihood of each genus belonging to large or small group , with equal distr
phy_eurasia$Small <- pnorm(phy_eurasia$LogMass, mean=modelpg$mu[1], sd=modelpg$sigma[1], lower.tail = F)
phy_eurasia$Large <- pnorm(phy_eurasia$LogMass, mean=modelpg$mu[2], sd=modelpg$sigma[2], lower.tail = T)
#Find where likelihoods are similar (green and red lines cross)
phy_eurasia <- phy_eurasia %>% mutate(Diff= Large - Small)
which.min(abs(phy_eurasia$Diff))  #most similar probabilities of large/small
phy_eurasia[1367,25] # log value of crossover
10^3.28

# by genus plots together
f.plot %>% 
  #insert_bottom(p.plot) %>% 
  insert_bottom(pg.plot) %>% 
  insert_bottom(pge.plot) 

# by species plots together
p.plot %>%  insert_bottom(pe.plot)


#Splitting fossil sites ####
ggplot()+ geom_histogram(data=use_subset , aes(x=mean))+ xlab("Mean of Log(Body Mass g)")+
  ylab("Site count") +
  ggtitle("Histogram of mean body masses \nfor NOW Pleistocene sites with 4 or more gen.")
colnames(use_subset)

#bimodal distribution
model_fossil <- normalmixEM(use_subset$mean, 
                       mu = c(1, 5), 
                       sigma = 1, 
                       maxit=200)

summary(model_fossil)

mixmdl.df = data.frame(x = use_subset$mean)

fossil.plot2 = ggplot(mixmdl.df) +
  geom_histogram(aes(x, ..density..),
                 binwidth = 0.2, colour = "black",
                 fill = "azure2",
                 alpha = 0.9
  ) +
  stat_function(
    geom = "line", fun = plot_mix_comps,
    args = list(model_fossil$mu[1], model_fossil$sigma[1], lam = model_fossil$lambda[1]),
    colour = "#FF000080", lwd = 1.5
  ) +
  stat_function(
    geom = "line", fun = plot_mix_comps,
    args = list(model_fossil$mu[2], model_fossil$sigma[2], lam = model_fossil$lambda[2]),
    colour = "#0000FF80", lwd = 1.5
  ) +
  ylab("Site Density") + xlab(element_blank())
fossil.plot2

#extract coordinates of both curves from plot
line.df = data.frame(x = ggplot_build(fossil.plot2)$data[[2]]$x,
                     red = ggplot_build(fossil.plot2)$data[[2]]$y,
                     blue = ggplot_build(fossil.plot2)$data[[3]]$y)

#find the minimal distance between lines along y axis
line.df$delta = line.df$red - line.df$blue

#find x coordinate for the minimal delta y
fossil.x_coord = line.df$x[which(diff(sign(diff((abs(line.df$delta))))) == 2)+1]
fossil.x_coord  #this is cutoff between two groups
10^fossil.x_coord # kg value of cutoff

fossil.plot2 <- fossil.plot2 + geom_vline(xintercept = fossil.x_coord, lty=2)

#trimodal distribution
model_fossil3 <- normalmixEM(use_subset$mean, 
                            k=3, 
                            mu=c(1, 3, 5),
                            sigma=1,
                            maxit=200)

summary(model_fossil3)

mixmdl.df = data.frame(x = use_subset$mean)

fossil.plot = ggplot(mixmdl.df) +
  geom_histogram(aes(x, ..density..),
                 binwidth = 0.2, colour = "black",
                 fill = "azure2",
                 alpha = 0.9
  ) +
  stat_function(
    geom = "line", fun = plot_mix_comps,
    args = list(model_fossil3$mu[1], model_fossil3$sigma[1], lam = model_fossil3$lambda[1]),
    colour = "#FF000080", lwd = 1.5
  ) +
  stat_function(
    geom = "line", fun = plot_mix_comps,
    args = list(model_fossil3$mu[2], model_fossil3$sigma[2], lam = model_fossil3$lambda[2]),
    colour = "purple2",alpha=0.8,  lwd = 1.5
  ) +
  stat_function(
    geom = "line", fun = plot_mix_comps,
    args = list(model_fossil3$mu[3], model_fossil3$sigma[3], lam = model_fossil3$lambda[3]),
    colour = "#0000FF80", lwd = 1.5
  ) +
  ylab("Site Density") + xlab(element_blank())
fossil.plot

#extract coordinates of both curves from plot
line.df = data.frame(x = ggplot_build(fossil.plot)$data[[2]]$x,
                     red = ggplot_build(fossil.plot)$data[[2]]$y,
                     blue = ggplot_build(fossil.plot)$data[[3]]$y)

#find the minimal distance between lines along y axis
line.df$delta = line.df$red - line.df$blue

#find x coordinate for the minimal delta y
fossil.x_coord.1 = line.df$x[which(diff(sign(diff((abs(line.df$delta))))) == 2)+1]
fossil.x_coord.1  #this is cutoff between two groups
10^fossil.x_coord.1 # kg value of cutoff

#extract coordinates of both curves from plot
line.df = data.frame(x = ggplot_build(fossil.plot)$data[[3]]$x,
                     red = ggplot_build(fossil.plot)$data[[3]]$y,
                     blue = ggplot_build(fossil.plot)$data[[4]]$y)

#find the minimal distance between lines along y axis
line.df$delta = line.df$red - line.df$blue

#find x coordinate for the minimal delta y
fossil.x_coord.2 = line.df$x[which(diff(sign(diff((abs(line.df$delta))))) == 2)+1]
fossil.x_coord.2  #this is cutoff between two groups
10^fossil.x_coord.2 # kg value of cutoff

fossil.plot3 <- fossil.plot+ geom_vline(xintercept = fossil.x_coord.1, lty=2) +
  geom_vline(xintercept = fossil.x_coord.2, lty=2)

fossil.plot2 %>% 
  insert_bottom(fossil.plot3)

#check which genera the difference in these cutoffs changes s/l category
org <- mmass_fin %>% filter(LogMass >3.14) 
ne <- mmass_fin %>% filter(LogMass >3.25)
spc <- mmass_fin %>% filter(LogMass >3.49)
org$Genus[!(org$Genus %in% ne$Genus)] 
ne$Genus[!(ne$Genus %in% spc$Genus)] 

# Look at community composition of site types ####

#Count large and small genera at each sites 

#name small and large genera (3.25 division)
#colnames(mmass_occ)
#mmass_nas <- mmass_occ[is.na(mmass_occ$LogMass),]
#mmass_occs <- mmass_occ[!(is.na(mmass_occ$LogMass)),]
#mmass_occs[mmass_occs$LogMass >= 3.25, 28] <- "Large"
#mmass_occs[mmass_occs$LogMass < 3.25, 28] <- "Small"
#mmass_occ <- rbind(mmass_occs, mmass_nas)
#mmass_occs %>% count(LargeSmall)

mmass_occ %>% count(LargeSmall)

lars=mmass_occ[mmass_occ$LargeSmall=="Large",3] #list of large genera
smas <- mmass_occ[mmass_occ$LargeSmall=="Small",3] # list of small genera
# select these genus names, plus other information in site/occurrence matrix
lar <- c("NAME", lars[!(is.na(lars))], colnames(SO[454:ncol(SO)]))  #colnames for large site/occ matrix
sma <- c("NAME", smas[!(is.na(smas))], colnames(SO[454:ncol(SO)])) #colnames for small site/occ matrix

#Name herbivore, non-herbivore
mmass_occ %>% count(Diet)
her=mmass_occ[mmass_occ$Diet=="Herbivore",3] #list of large genera
non <- mmass_occ[mmass_occ$Diet=="Non-Herbivore",3] # list of small genera
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


#Assign size ratio/ proportion
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
#Assign diet ratio/ proportion
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

write.csv(SO_use, "BKV_siteoccmatrix_DietSizeCounts_1.26.csv")
SO_use <- read.csv("BKV_siteoccmatrix_DietSizeCounts_1.26.csv", header=T, row.names = 1)

#3 group division
largeonly <- use_subset %>% filter(mean > fossil.x_coord.2) 
medonly <- use_subset %>% filter(mean >= fossil.x_coord.1) %>% filter(mean <=fossil.x_coord.2)  
smallonly <- use_subset %>% filter(mean < fossil.x_coord.1)
# 2 group division
largesites <- use_subset %>% filter(mean >= fossil.x_coord)
smallsites <- use_subset %>% filter(mean < fossil.x_coord)

# then check number of large mammal genera at sites vs. small 
# Input matrices for large and small mammal sites, excluding genera from opposite group
largesite_mat <- SO_use %>% filter(NAME %in% largesites$NAME)
smallsite_mat <- SO_use %>% filter(NAME %in% smallsites$NAME)
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

sd(largeonly_mat$Gen_Count)
sd(largeonly_mat$Small_GenCount)
sd(largeonly_mat$smallprop)
sd(largeonly_mat$HerbProp)



#For datasci project, write matrices
#largesite_mat_L <- largesite_mat[, which((names(largesite_mat) %in% lar)==TRUE)]
#largesite_mat_L$Large_GenCount <- rowSums(largesite_mat_L[,2:238])
#setwd("C:/Users/abigpark/OneDrive - University of Helsinki/Data Science RS")
#write.csv(largesite_mat_L, "LargeSites_SiteOccurrences_LargeGenera_26.1.24.csv")
#write.csv(SO_use, "AllSites_SiteOccurrences_AllGenera_26.1.24.csv")

#How these subsets compare in age range
Mage <- ggplot(medonly_mat)+ geom_histogram(aes(x=age_range))+ ylim(0,60)+xlab(element_blank())+ ylab("Mixed Site Count")
Sage <- ggplot(smallonly_mat)+ geom_histogram(aes(x=age_range))+ ylim(0,60)+ ylab("Small Site Count")
Lage <- ggplot(largeonly_mat)+ geom_histogram(aes(x=age_range))+ ylim(0,60)+xlab(element_blank())+ ylab("Large Site Count")
LSage <- ggplot(largesite_mat)+ geom_histogram(aes(x=age_range))+ ylim(0,60)+xlab(element_blank())+ ylab("Large Site Count")
SSage <- ggplot(smallsite_mat)+ geom_histogram(aes(x=age_range))+ ylim(0,60)+ ylab("Small Site Count")

allage <- ggplot(use_subset)+ geom_histogram(aes(x=age_range), fill="blue3") + ylab("Total Site Count")+xlab(element_blank())
allage %>% 
insert_bottom(Lage) %>% insert_bottom(Mage) %>% insert_bottom(Sage)

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

#Find ratios
hid_pts <- hid_pts %>% mutate(DietRatio= Herb_GenCount/Nonherb_GenCount)
hid_pts <- hid_pts %>% mutate(Herbprop= Herb_GenCount/Gen_Count)

summary(hid_pts[15:23])
sd(hid_pts$Nonherb_GenCount)
sd(hid_pts$smallprop)
sd(hid_pts$Herbprop)
sd(hid_pts$Gen_Count)

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

ncol(Mod_largeherb_occ)-1
ncol(Mod_largenonherb_occ)-1
ncol(Mod_smallherb_occ)-1
ncol(Mod_smallnonherb_occ)-1

write.csv(hid_pts, "ModernEurasia_DietSizeCounts.csv")
hid_pts <- read.csv("ModernEurasia_DietSizeCounts.csv", row.names = 1, header=T)

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


# simple PCA visualizing these 4 categories
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

# run PCA
library(psych)
pc <- prcomp(guilds[,2:5],   #here, uses singular value decomposition (usually preferred)
             center = TRUE,  #T makes all variables centered at 1
             scale. = TRUE) #scale is for normalization of PC scores 
              
print(pc)  #has yielded 14 PC's
summary(pc)

guilds$PC1 <- pc$x[,1]
guilds$PC2 <- pc$x[,2]

#plot
cols <- c("LargeSites"="#0000FF80",  "MixedSites"="#912CEE98",
          "SmallSites" ="#FF000080", "ModernEurasia"="#77889990")

ggplot(guilds)+geom_point(aes(x=smallprop, y=HerbProp, color=Dataset), size=3)+ 
  scale_color_manual(values=cols)+xlab("Proprtion small mammals")+ylab("Proportion herbivores")

ggplot(guilds)+geom_point(aes(x=PC1, y=PC2, color=Dataset), size=3)+ ggtitle("PCA of guild composition")+
  scale_color_manual(values=cols)+xlab("PC1")+ylab("PC2")

#Looking at geography of PF points outside of convex hull of modern PC1 and PC2-- or try PC3 instead if less
      # related to total genus count
library(geometry)
find_hull_bydataset <- function(guilds) guilds[chull(guilds$smallprop, guilds$HerbProp),]
hulls_Dataset <- ddply(guilds, "Dataset", find_hull_bydataset)

hull_dt <- guilds %>%
  group_by(Dataset) %>%
  slice(chull(smallprop, HerbProp))
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
96/119
31/116

#Geog/cliamte gradients within modern hull
colnames(hid_pts_use)
tmod <- ggplot(hid_pts)+geom_point(aes(x=smallprop, y=Herbprop, color=BIO01_Mean), size=3)+ 
  geom_polygon(data = hull_mod, alpha = 0.2, aes(x=smallprop, y=HerbProp), fill="grey", 
               lwd=1, lty=4)+
  scale_color_scico(palette="lajolla", name="MAT")+xlab("Proportion small mammals")+ylab("Proportion herbivores")

lmod <- ggplot(hid_pts)+geom_point(aes(x=smallprop, y=Herbprop, color=Y), size=3)+ 
  geom_polygon(data = hull_mod, alpha = 0.2, aes(x=smallprop, y=HerbProp), fill="grey", 
               lwd=1, lty=4)+ ylab(element_blank())+
  scale_color_scico(palette="roma", name="Latitude")+xlab("Proportion small mammals")+ylab("Proportion herbivores")

tmod  %>% insert_right(lmod)

vcomp <- guilds %>% dplyr::select(PC1, PC2, LargeHerb_GenCount,LargeNonherb_GenCount,
                                           SmallHerb_GenCount, SmallNonherb_GenCount)
ucomp <- guilds %>% dplyr::select(LargeHerb_GenCount,LargeNonherb_GenCount,
                                  SmallHerb_GenCount, SmallNonherb_GenCount)
MODcomp <- hid_pts_use %>% dplyr::select(LargeHerb_GenCount,LargeNonherb_GenCount,
                                  SmallHerb_GenCount, SmallNonherb_GenCount)
library(corrplot)
col2 <- colorRampPalette(brewer.pal(9,"RdBu"))

correlation_df<-cor(MODcomp, method="pearson", use="pairwise.complete.obs")
corrplot(correlation_df, method = "square", # this is better matrix, indicates few high correlation values
         col=col2(20),tl.col="black",
         addCoef.col = "black", diag=F) #type="lower")

#simple GLLVM for 4 categories: l/s, herb/nonherb
#install.packages("gllvm")
library(gllvm) #error in this function
LVres <- gllvm(y = guildsS[,c(2:5)], family = "negative.binomial", num.lv = 1, 
               row.eff="random", maxit=1000)
guilds$LV1 <- LVres$lvs[,1]
#write.csv(guilds, "BKV_modern8000fossil_LV.csv")
LVres
plot(LVres)
coef(LVres)
ordiplot(LVres, biplot=T)
getResidualCov(LVres)$var.q
lv_scale <- scale(as.matrix(LVres$lvs))

lvres <- read.csv("BKV_modern8000fossil_LV.csv", header=T, row.names = 1)
lvres <- lvres %>% mutate(LVbig=100000*LV1)
summary(lvres$LVbig)
hist(lvres$LV1)

cols <- c("LargeSites"="#0000FF80",  "MixedSites"="#912CEE98",
          "SmallSites" ="#FF000080", "ModernEurasia"="#778899")

lvres$DS <- factor(lvres$Dataset, levels = c( "ModernEurasia", "LargeSites",
                                              "MixedSites", "SmallSites"), ordered = TRUE)

ggplot(data=lvres) + geom_density_ridges(aes(x=LVbig,y=Dataset,
                                             fill=DS), linewidth=0.8) +
  theme_bw() + scale_fill_manual(values=cols)+ xlab("LV1") +ylab("Density")+
  theme(legend.position = "none")

lvres %>% count(Dataset)   

large <- lvres %>% filter(Dataset=="LargeSites")
mixed <- lvres %>% filter(Dataset=="MixedSites")
small <- lvres %>% filter(Dataset=="SmallSites")
fossil <- lvres %>% filter(Dataset %in% c("LargeSites", "MixedSites", "SmallSites"))
modern <- lvres %>% filter(Dataset=="ModernEurasia")

t.test(large$LVbig, modern$LVbig)
t.test(mixed$LVbig, modern$LVbig) #mixed sites not significantly different from modern
t.test(small$LVbig, modern$LVbig)
t.test(fossil$LVbig, modern$LVbig) #all other significantly different

#Kolmogorov-Smirnov test for equality of distributions
ks.test(large$LV1, modern$LV1)
ks.test(mixed$LV1, modern$LV1) #mixed sites not significantly different from modern
ks.test(small$LV1, modern$LV1)
ks.test(fossil$LV1, modern$LV1) #all other significantly different


#Compare large/small and diet counts from sites and modern samples ####
SO_dat <- SO_use[,c(1,454:ncol(SO_use))]
sl_comp <- merge(sites[,1:2], SO_dat, by="NAME") #save nearest modern HID from earlier
mod_comp <- hid_pts_use
colnames(mod_comp) <- paste(colnames(hid_pts_use), "mod", sep="_")
sl_comp <- merge(sl_comp, mod_comp, by.x="Nearest_HID", by.y="HID_mod")
lfm <- ggplot(sl_comp)+ geom_point(aes(x=Large_GenCount, y= Large_GenCount_mod, color=LAT), size=2) +
  scale_color_scico(palette="turku", name="Latitude")+ xlab("Fossil Site Large Genus Count")+ 
  ylab ("Nearest Modern Sample \nLarge Genus Count")+ ylim(0,30)+ xlim(0,30)+  theme_bw() +
  geom_abline(intercept = c(0,0), slope = 1, color="slategrey") + theme(legend.position="none")
  
sfm <- ggplot(sl_comp)+ geom_point(aes(x=Small_GenCount, y= Small_GenCount_mod, color=LAT), size=2) +
  scale_color_scico(palette="turku", name="Latitude")+ xlab("Fossil Site Small Genus Count")+ 
  ylab ("Nearest Modern Sample \nSmall Genus Count")+ ylim(0,30)+ xlim(0,30)+ theme_bw() +
  geom_abline(intercept = c(0,0), slope = 1, color="slategrey") + theme(legend.position="none")

hfm <- ggplot(sl_comp)+ geom_point(aes(x=Herb_GenCount, y= Herb_GenCount_mod, color=LAT), size=2) +
  scale_color_scico(palette="turku", name="Latitude")+ xlab("Fossil Site Herbivore Genus Count")+
  theme_bw() +
  ylab ("Nearest Modern Sample \nHerbivore Genus Count")+ ylim(0,30)+ xlim(0,30)+
  geom_abline(intercept = c(0,0), slope = 1, color="slategrey") + theme(legend.position="none")

nfm <- ggplot(sl_comp)+ geom_point(aes(x=Nonherb_GenCount, y= Nonherb_GenCount_mod, color=LAT), size=2) +
  theme_bw()+ theme(legend.position="none")+ ylim(0,30)+ xlim(0,30)+ 
  scale_color_scico(palette="turku", name="Latitude")+ xlab("Fossil Site Non-Herbivore Genus Count")+ 
  ylab ("Nearest Modern Sample \nNon-Herbivore Genus Count")+
  geom_abline(intercept = c(0,0), slope = 1, color="slategrey") 

library(gridExtra)
grid.arrange(lfm, hfm, sfm, nfm, ncol=2)

#check for other relationships in modern
ggplot(hid_pts_use)+ geom_point(aes(x=BIO01_Mean, y=Herbprop, color=Y))
ggplot(hid_pts_use)+ geom_point(aes(x=BIO12_Mean, y=Herbprop, color=Y))
ggplot(hid_pts_use)+ geom_point(aes(x=NPP, y=smallprop, color=Y))

#Check how many species vs genus count at modern sites

#Modern by species count per site
spl <- cbind(A1EuraSpList,A2eur_spec, A5eur_spec, A6eur_spec, A7eur_spec, A8eur_spec)
jdatS <- spl[,3:ncol(spl)] #select species occurrence columns
jdatS <- jdatS %>% mutate(SpCount=rowSums(!(is.na(.))))  #count how many occurrences per row
spl$SpCount <- jdatS$SpCount   
summary(spl$SpCount) #up to 199 species at sample points
#Merge this number with Genus counts of HID's
hex_sf$SpCount <- spl$SpCount
ggplot(hex_sf)+ geom_point(aes(x=Gen_Count, y= SpCount, color=Y), size=3) + 
  scale_color_scico(palette="corkO", name="Latitude")+
  ylab ("Species Count")+ xlab("Genus Count")+
  geom_abline(intercept = c(0,0), slope = 1, color="slategrey")
hex_sf <- hex_sf %>% mutate(Species_per_genus=SpCount / Gen_Count)
hex_sf$Species_per_genus[is.infinite(hex_sf$Species_per_genus)] <- NA
summary(hex_sf$Species_per_genus)

#Species: genus counts for only large mammals 
splL <- cbind(A5eur_spec, A6eur_spec, A7eur_spec, A8eur_spec)
jdatS_L <- splL[,3:ncol(splL)] #select species occurrence columns
jdatS_L <- jdatS_L %>% mutate(SpCount=rowSums(!(is.na(.))))  #count how many occurrences per row
spl$SpCount_medlarge <- jdatS$SpCount   
summary(spl$SpCount_medlarge) #up to 199 species at sample points
#Merge this number with Genus counts of HID's
hex_sf$SpCount_medlarge <- spl$SpCount_medlarge
ggplot(hex_sf)+ geom_point(aes(x=Gen_Count, y= SpCount_medlarge, color=Y), size=3) + 
  scale_color_scico(palette="corkO", name="Latitude")+
  ylab ("Species Count")+ xlab("Genus Count")
hex_sf <- hex_sf %>% mutate(SpeciesML_per_genus=SpCount / Gen_Count)
hex_sf$SpeciesML_per_genus[is.infinite(hex_sf$SpeciesML_per_genus)] <- NA
summary(hex_sf$SpeciesML_per_genus)
