# Divide small/large mammal genera ####
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

# Matrix build: for every fossil site, body mass for each occurring genus in row
EM <- matrix(nrow=1,ncol=(1+max(SO$Gen_Count))) #creating empty data frame, for HID plus up to 52 species names
MassFrame <- as.data.frame(EM)

colnames(mmass_fin)  #logmass of genus is col 5

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
    ma <- mmass_fin[mmass_fin$Genus==genus,5]  #find logmass
    V[p+1] <- ma[1]
  }
  MassFrame <- rbind(MassFrame, V)
}

InputMass <- MassFrame[2:nrow(MassFrame),]
colnames(InputMass)[1] <- "loc_name"
write.csv(InputMass,"Pleist_AllEurasia_InputMassDist.csv")
InputMass <- read.csv("Pleist_AllEurasia_InputMassDist.csv", row.names = 1) #read in if previously calculated

# Find distribution statistics for masses at each input site
fossil_pts <- SO %>% dplyr::select(NAME, LAT, LONG, ALTITUDE, MAX_AGE, MIN_AGE, BFA_MAX, BFA_MAX_ABS, BFA_MIN, BFA_MIN_ABS,
                                 COUNTRY, age_range, Gen_Count) %>% 
                                  mutate(mid_age= ((MAX_AGE+MIN_AGE)/2))   #find site info

#Find mean body mass at each site
fossil_pts$mean <- rowMeans(InputMass[,2:ncol(InputMass)], na.rm=T)

#These other distribution statistics can be used for other analysis
sd <- rep(NA, nrow(InputMass))
for (p in 1:nrow(InputMass)){      
  sd[p] <- sd(InputMass[p,2:ncol(InputMass)],na.rm=T)
}
fossil_pts$SD <- sd

mx <- rep(NA, nrow(InputMass))
for (p in 1:nrow(InputMass)){            
  mx[p] <- max(InputMass[p,2:ncol(InputMass)],na.rm=T)
}
mx[is.infinite(mx)] <- NA
fossil_pts$max <- mx

mn <- rep(NA, nrow(InputMass))
for (p in 1:nrow(InputMass)){            
  mn[p] <- min(InputMass[p,2:ncol(InputMass)],na.rm=T)
}
mn[is.infinite(mn)] <- NA
fossil_pts$min <- mn

# Look at mean body mass at all fossil sites, App. Fig. 3
ggplot()+ geom_histogram(data=fossil_pts , aes(x=mean))+ xlab("Site Mean of Log(Body Mass g)")+
  ylab("Site count") 

# select sites with 4 or more genera occurring
use_subset <- fossil_pts %>% filter(Gen_Count >3)

ggplot()+ geom_histogram(data=use_subset , aes(x=mean))+ xlab("Mean of Log(Body Mass g)")+
  ylab("Site count") +
  ggtitle("Histogram of mean body masses \nfor NOW Pleistocene sites with 4 or more gen.")

# For these genera occurring in NOW, look at mass distribution (now by taxa, not site averages)
avfree <- mmass_fin %>% filter(!(is.na(LogMass))) %>%  filter(!(MassSource=="Family average")) %>% 
  filter(Genus %in% gen)  #genera occurring in NOW, excluding those with family average masses

ggplot()+ geom_histogram(data=avfree, aes(x=LogMass), binwidth=0.2)

#Bimodal distribution: fit model
model <- normalmixEM(avfree$LogMass,
                     mu = c(1, 5), 
                     sigma = 1,
                     maxit = 200)
summary(model) #model output: mu are means, sigma are sd's, lambda is proportion in each group (comp)

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

f.plot <- f.plot+ geom_vline(xintercept = f.x_coord, lty=2) # bimodal plot for fossil genera


#Modern mammals instead (Phylacine database)

#first by species ####
#trait data, download directly from Phylacine database: https://megapast2future.github.io/PHYLACINE_1.2/
phy <- read.csv("Trait_data_Phylacine.csv", header=T)  
# removed species with phylogenetically imputed masses
phy <- phy %>% mutate(LogMass= log10(Mass.g)) %>% filter(!(Mass.Method == "Imputed"))  
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
aq_gen <- rbind(marine, fwat, "Enhydra","Neovison") #Lontra more terrestrial species, not listed as aquatic genus
aq_gen <- aq_gen$Genus.1.2

#without imputed species masses
ggplot()+ geom_histogram(data=phy_terr, aes(x=LogMass))+ xlab("Log(Body Mass g)")+
  ylab("Species count") +
  ggtitle("Histogram of mean body masses \nfor modern terrestrial mammal species")

# fit bimodal model
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

#then by genus ####
# find genus averages of mass and diet percent plant
phy_gen_av <- phy %>% group_by(Genus.1.2) %>% summarize(Plant = mean(Diet.Plant, na.rm = TRUE),  #genus average %plant
                                                        LogMass = mean(log10(Mass.g),na.rm=T)) #save genus average mass (not including imputed species)

phy_gen_av_terr <- phy_gen_av %>% filter(!(Genus.1.2 %in% aq_gen)) #remove aquatic genera

ggplot()+ geom_histogram(data=phy_gen_av_terr, aes(x=LogMass))+ xlab("Log(Body Mass g)")+
  ylab("Species count") +
  ggtitle("Histogram of mean body masses \nfor modern mammal genera")

#fit bimodal model
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

# Phylacine genera, but only those occurring in Eurasia ####
#genera occurring in Eurasia 
colnames(Genus_terr_eurasia)[1:10]
colnames(Genus_terr_eurasia)[1330:1335]
EUmat_occs <- Genus_terr_eurasia[,2:1334] #jsut occurrence columns
Eurasia_occ <- EUmat_occs[which (colSums(EUmat_occs)>0)] #genus columns with occurrences
Eurasia_genera <- colnames(Eurasia_occ)
#Find genera occurring in Eurasia
gen_av_eurasia <- phy_gen_av %>% filter(Genus.1.2 %in% Eurasia_genera)

ggplot()+ geom_histogram(data=gen_av_eurasia, aes(x=LogMass))+ xlab("Log(Body Mass g)")+
  ylab("Genus count") +
  ggtitle("Histogram of mean body masses \nfor modern mammal genera in Eurasia")

#fit bimodal model
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

#Divide large and small, only Eurasia species

#list all species occurring in Eurasia
Eurasia_species <- c(A1eur_species, A2eur_species, A5eur_species, A6eur_species,   
                     A7eur_species, A8eur_species)  #these are lists from other works of species within size/diet guilds
Eurasia_species <- unique(Eurasia_species)
Eurasia_species <- Eurasia_species[!is.na(Eurasia_species)]  #identified 760 species occurring in present Eurasia

phy_eurasia <- phy_terr %>% filter(Binomial.1.2 %in% Eurasia_species) #only species in Eurasia, 1877 species

#without imputed species masses
ggplot()+ geom_histogram(data=phy_eurasia, aes(x=LogMass))+ xlab("Log(Body Mass g)")+
  ylab("Species count") +
  ggtitle("Histogram of mean body masses \nfor modern mammal species in Eurasia")

# fit bimodal model 
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

# Fig. 2: by genus plots together
f.plot %>% 
  insert_bottom(pg.plot) %>% 
  insert_bottom(pge.plot) 

# App. Fig. 2: by species plots together
p.plot %>%  insert_bottom(pe.plot)


#Splitting fossil sites into 2/3 size categories ####
ggplot()+ geom_histogram(data=use_subset , aes(x=mean))+ xlab("Mean of Log(Body Mass g)")+
  ylab("Site count") +
  ggtitle("Histogram of mean body masses \nfor NOW Pleistocene sites with 4 or more gen.")

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

# fit trimodal distribution
model_fossil3 <- normalmixEM(use_subset$mean, 
                             k=3, 
                             mu=c(1, 3, 5),
                             sigma=1,
                             maxit=200)

summary(model_fossil3)

mixmdl.df = data.frame(x = use_subset$mean)

fossil.plot3 = ggplot(mixmdl.df) +
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
fossil.plot3

#extract coordinates of both curves from plot
line.df = data.frame(x = ggplot_build(fossil.plot3)$data[[2]]$x,
                     red = ggplot_build(fossil.plot3)$data[[2]]$y,
                     blue = ggplot_build(fossil.plot3)$data[[3]]$y)

#find the minimal distance between lines along y axis
line.df$delta = line.df$red - line.df$blue

#find x coordinate for cutoff between small and mixed groups
fossil.x_coord.1 = line.df$x[which(diff(sign(diff((abs(line.df$delta))))) == 2)+1]
fossil.x_coord.1  #this is cutoff between two groups
10^fossil.x_coord.1 # kg value of cutoff

#extract coordinates of both curves from plot
line.df = data.frame(x = ggplot_build(fossil.plot)$data[[3]]$x,
                     red = ggplot_build(fossil.plot)$data[[3]]$y,
                     blue = ggplot_build(fossil.plot)$data[[4]]$y)

#find the minimal distance between lines along y axis
line.df$delta = line.df$red - line.df$blue

#find x coordinate for cutoff between mixed and large groups
fossil.x_coord.2 = line.df$x[which(diff(sign(diff((abs(line.df$delta))))) == 2)+1]
fossil.x_coord.2  #this is cutoff between two groups
10^fossil.x_coord.2 # kg value of cutoff

fossil.plot3 <- fossil.plot3+ geom_vline(xintercept = fossil.x_coord.1, lty=2) +
  geom_vline(xintercept = fossil.x_coord.2, lty=2)

#Plot Fig. 3
fossil.plot2 %>% 
  insert_bottom(fossil.plot3)

#check which genera the difference in these cutoffs changes s/l category
org <- mmass_fin %>% filter(LogMass >3.14) 
ne <- mmass_fin %>% filter(LogMass >3.25)
spc <- mmass_fin %>% filter(LogMass >3.49)
org$Genus[!(org$Genus %in% ne$Genus)] 
ne$Genus[!(ne$Genus %in% spc$Genus)] 

#Define subsets of sites for further analysis (genus counts by size/diet in 4_GuildComposition)
#3 group division
largeonly <- use_subset %>% filter(mean > fossil.x_coord.2) 
medonly <- use_subset %>% filter(mean >= fossil.x_coord.1) %>% filter(mean <=fossil.x_coord.2)  
smallonly <- use_subset %>% filter(mean < fossil.x_coord.1)
# 2 group division
largesites <- use_subset %>% filter(mean >= fossil.x_coord)
smallsites <- use_subset %>% filter(mean < fossil.x_coord)


#How these site subsets compare in age range
Mage <- ggplot(medonly)+ geom_histogram(aes(x=age_range))+ ylim(0,60)+xlab(element_blank())+ ylab("Mixed Site Count")
Sage <- ggplot(smallonly)+ geom_histogram(aes(x=age_range))+ ylim(0,60)+ ylab("Small Site Count")
Lage <- ggplot(largeonly)+ geom_histogram(aes(x=age_range))+ ylim(0,60)+xlab(element_blank())+ ylab("Large Site Count")

allage <- ggplot(use_subset)+ geom_histogram(aes(x=age_range), fill="blue3") + ylab("Total Site Count")+xlab(element_blank())

allage %>% 
  insert_bottom(Lage) %>% insert_bottom(Mage) %>% insert_bottom(Sage) #similar age range distributions

#Map these 3 site subsets
 # sites sf object defined in 2_GenusCountComparison
sites <- sites %>% mutate(SiteType = rep(NA))
sites$SiteType[sites$NAME %in% largeonly$NAME] <-  "Large"
sites$SiteType[sites$NAME %in% medonly$NAME] <-  "Mixed"
sites$SiteType[sites$NAME %in% smallonly$NAME] <-  "Small"

site_plot <- sites %>% filter(!(is.na(SiteType)))

ggplot()+ geom_sf(data=eurasia)+  #plot genus count of RS output
  geom_sf(data=site_plot, aes(color=SiteType), shape=18, size=3)+
  scale_color_manual(values= c("Large"="#0000FF80", "Mixed"="#912CEE",
                     "Small" ="#FF000080"), name= "Site Type")+
  coord_sf(xlim=c(-12, 170), ylim= c(-10,80), expand=F)
  
