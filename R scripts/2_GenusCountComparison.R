# Code to compare fossil genus count to closet modern community genus count
# and plot Fig. 1 (maps of faunas with lat/long histograms colored by mean genus diversity

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

# compare genus count to genus count at closest modern Eco-ISEA3H sample point ####

# read in geography of hexagon ID's ####
HIDs <- read.csv("EcoISEA3H_TerrestrialHIDs_continents.csv", header=T, row.names=1)   #hexagon info and climate data
colnames(HIDs)
EUHIDs <- HIDs %>% filter(Eurasia=="TRUE") %>% dplyr::select(HID, X, Y, Realm_Mode, Tree_Mean, Elevation_Mean,
                                                             BIO01_Mean, BIO04_Mean, BIO05_Mean,
                                                             BIO06_Mean, BIO12_Mean, 
                                                             BIO15_Mean, NPP, whichNPP)

#read in modern genus occurrences
#These were previously processed drawing species/genus occurrences from the Phylacine database for each Eco-ISEA3H hexagonal grid cell centroid point
GenusMatrix <- read.csv("Phylacine_Present_TerrestrialGenusOccurrences_Res9_Jan24.csv", header=T, row.names=1)
GenusCountedMatrix <- read.csv("Phylacine_Present_TerrestrialGenusOccurrences_spcounts_Res9_Jan24.csv", header=T, row.names=1)

#find genus count per point
sum <- rowSums(GenusMatrix[,2:1334]) #0/1 occurrence columns
summary(sum) 
GenusMatrix$N_Gen <- as.numeric(sum)

#Select Eurasian points
Genus_terr_eurasia <- GenusMatrix %>% filter(HID %in% EUHIDs$HID)
Genus_sp_eurasia <- GenusCountedMatrix %>% filter(HID %in% EUHIDs$HID)
summary(Genus_terr_eurasia$N_Gen) #up to 111 co-occurring genera at these sites

#Find how many species occur within each genus at each HID

f2 <- function(data, levels) { 
  tt <- function(x, levels) tabulate(factor(x, levels), length(levels))
  cc <- vapply(data, tt, integer(2L), levels, USE.NAMES = FALSE)
  res <- as.integer(.rowSums(cc, 2L, length(data))) 
  names(res) <- levels
  res
}
genpre <- f2(Genus_terr_eurasia[,2:1334], c(0,1))  #this counts 0s and 1s in modern occurrence matrix

# This counts species within each genus at occurrences
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

#make Eco-ISEA3H hexagon ID matrix into sf file
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
# and save genus count with climate/geography data for modern points
hid_pts <- cbind(EUHIDs, Genus_terr_eurasia[,1335]) #Add Genus count to hex ids
colnames(hid_pts)[15] <-  "Gen_Count"

#save lat/long as x/y for plotting below
sites <- merge(sites, siterich_latlong[,c(1,3:4)], by="NAME") 


#plot fossil vs. nearest modern genus count, Fig. 4
ggplot(sites)+ geom_point(aes(x=Gen_Count, y= Modern_GenCount, color=Y), size=2) + 
  scale_color_scico(palette="turku", name="Latitude")+ xlab("Fossil Site Genus Count")+ 
  ylab ("Nearest Modern Sample Genus Count")+ xlim(c(0,55)) + ylim(c(0,115))+
  geom_abline(intercept = c(0,0), slope = 1, color="slategrey")+
  scale_x_continuous(breaks=c(0,10,20,30,40,50))
#alternative colored by longitude
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

# Plot App. Fig. 5
fossilgenlat %>%    
  insert_bottom(moderngenlat)

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

# PLot Fig. 5
fossilgenlong %>% 
  insert_bottom(moderngenlong)

#Plotting Fig. 1 ####
#3-part graph: map of fossil sites with genus count
#define color scale
palette_gen <- colorRampPalette(colors = c("lightskyblue1",
                                           "darkorchid4", "midnightblue"))(20)
scales::show_col(palette_gen)
palette_gen #use first half of color scale for less-genus rich fossil plot, full scale for modern

# Map of fossil sites
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

# alternate bar graph: summarizing average genus count by 5 degree long bins in bar heights
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
latB_mean

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

#3-part graph: map of modern samples with genus count

hex_sf <- cbind(hex_sf, hid_pts[,c(2:3,15)]) #Add Genus count to hex ids sf object
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

#combined plot (Fig. 1B)
ggplot()+
  coord_equal(xlim = c(-4, 50), ylim = c(0, 45), expand = FALSE) +
  annotation_custom(ggplotGrob(modmap), xmin = 0, xmax = 39, ymin = 3, 
                    ymax = 27) +
  annotation_custom(ggplotGrob(MlongB), xmin = -3, xmax = 40, ymin = 27, 
                    ymax = 34)+
  annotation_custom(ggplotGrob(MlatB), xmin = 39, xmax = 46, ymin = 1, 
                    ymax = 28) +
  theme(panel.spacing = unit(0, "lines"))
