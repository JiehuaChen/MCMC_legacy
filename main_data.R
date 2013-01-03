#data
#MCMC for multilevel model
#library(gstat, pos=4)
library(sp, pos=4)
library(lattice)
library(proj4)
library(RandomFields)
library(geoR)
require(rgdal)
library(arm)

#LS regression

afhori <- read.table("Horizon.csv", header=T, sep=",", stringsAsFactors=FALSE)
afprof <- read.table("Profile.csv", header=T, sep=",")
covprof <- read.table("legacy_profile_withcov.csv", header=T, sep=",")

afprof.withcov <- cbind(afprof, covprof)
aflegacy <- merge(afhori, afprof.withcov, by="SPID")


#aflegacy <- cbind(aflegacy, lambertx = aflegacy.lambertcord[,1], lamberty = aflegacy.lambertcord[,2])
aflegacy <- as.data.frame(aflegacy)
locations <- cbind(x=aflegacy$coords.x1, y=aflegacy$coords.x2)
locations <- as.data.frame(locations)
aflegacy <- aflegacy[!is.na(locations$x+locations$y), ]
locations <- locations[!is.na(locations$x+locations$y), ]

 # aflegacy$bio1 <- ifelse(aflegacy$bio1==0, NA, aflegacy$bio1) 
 # aflegacy$bio12 <- ifelse(aflegacy$bio12==0, NA, aflegacy$bio12) 
 # aflegacy$CTI_1K <- ifelse(aflegacy$CTI_1K==0, NA, aflegacy$CTI_1K) 
 # aflegacy$ELEV_1K <- ifelse(aflegacy$ELEV_1K==0, NA, aflegacy$ELEV_1K) 
 # aflegacy$EVIM_1K <- ifelse(aflegacy$EVIM_1K==0, NA, aflegacy$EVIM_1K) 
 # aflegacy$M13RB1ALT <- ifelse(aflegacy$M13RB1ALT==0, NA, aflegacy$M13RB1ALT) 
 # aflegacy$M13RB2ALT <- ifelse(aflegacy$M13RB2ALT==0, NA, aflegacy$M13RB2ALT) 
 # aflegacy$M13RB3ALT <- ifelse(aflegacy$M13RB3ALT==0, NA, aflegacy$M13RB3ALT) 
 # aflegacy$M13RB7ALT <- ifelse(aflegacy$M13RB7ALT==0, NA, aflegacy$M13RB7ALT) 
 # aflegacy$NPP_Mean_1 <- ifelse(aflegacy$NPP_Mean_1==0, NA, aflegacy$NPP_Mean_1) 
 # aflegacy$RELIEF_1K <- ifelse(aflegacy$RELIEF_1K==0, NA, aflegacy$RELIEF_1K) 
 # aflegacy$SLOPE_1K <- ifelse(aflegacy$SLOPE_1K==0, NA, aflegacy$SLOPE_1K)

aflegacy <-aflegacy[aflegacy$EVI_mask_1==1, ]


index.na <- (1:dim(aflegacy)[1])*is.na(rowSums(aflegacy[, c("SOC", "Bot", "bio1", "bio12", "CTI_1K", "ELEV_1K", "EVIM_1K", "M13RB1ALT", "M13RB2ALT", "M13RB3ALT", "M13RB7ALT", "NPP_Mean_1", "RELIEF_1K", "SLOPE_1K", "lstday", "lstnight")]))
aflegacy <- aflegacy[-index.na, ]
aflegacy <- aflegacy[aflegacy$lstday>400, ]

logccm <- log(aflegacy$SOC)
N <- length(logccm)
log.depth <- log(aflegacy$Bot)

bio1.scaled <- as.vector(scale(aflegacy$bio1)) 
bio12.scaled<- as.vector(scale(aflegacy$bio12 )) 
CTI_1K.scaled <- as.vector(scale(aflegacy$CTI_1K)) 
ELEV_1K.scaled <- as.vector(scale(aflegacy$ELEV_1K)) 
EVIM_1K.scaled <- as.vector(scale(aflegacy$EVIM_1K)) 
M13RB1ALT.scaled <- as.vector(scale(aflegacy$M13RB1ALT)) 
M13RB2ALT.scaled <- as.vector(scale(aflegacy$M13RB2ALT )) 
M13RB3ALT.scaled <- as.vector(scale(aflegacy$M13RB3ALT)) 
M13RB7ALT.scaled <- as.vector(scale(aflegacy$M13RB7ALT)) 
NPP_Mean.scaled <- as.vector(scale(aflegacy$NPP_Mean_1)) 
RELIEF.scaled <- as.vector(scale(aflegacy$RELIEF)) 
SLOPE_1K.scaled <- as.vector(scale(aflegacy$SLOPE_1K)) 
lstday.scaled <- as.vector(scale(aflegacy$lstday))
lstnight.scaled <- as.vector(scale(aflegacy$lstnight))

Y <- log(aflegacy$SOC)
X <- data.frame(int=1, bio1.scaled=bio1.scaled, bio12=bio12.scaled, CTI_1K.scaled=CTI_1K.scaled, 
              ELEV_1K.scaled=ELEV_1K.scaled, EVIM_1K.scaled=EVIM_1K.scaled, M13RB1ALT.scaled=M13RB1ALT.scaled,
              NPP_Mean.scaled=NPP_Mean.scaled, RELIEF.scaled=RELIEF.scaled, lstday.scaled=lstday.scaled, lstnight.scaled=lstnight.scaled, bio1.scaled.depth=bio1.scaled*log.depth, bio12=bio12.scaled*log.depth, CTI_1K.scaled.depth=CTI_1K.scaled*log.depth, 
              ELEV_1K.scaled.depth=ELEV_1K.scaled*log.depth, EVIM_1K.scaled.depth=EVIM_1K.scaled*log.depth, M13RB1ALT.scaled.depth=M13RB1ALT.scaled*log.depth,
              NPP_Mean.scaled.depth=NPP_Mean.scaled*log.depth, RELIEF.scaled.depth=RELIEF.scaled*log.depth,lstday.scaled.depth=lstday.scaled*log.depth, lstnight.scaled.depth=lstnight.scaled*log.depth, log.depth=log(aflegacy$Bot))
X <- as.matrix(X)

aflegacy <- aflegacy[!is.na(rowSums(cbind(Y,X))), ]

aflegacy.lambertcord <- project(cbind(aflegacy$Lon, aflegacy$Lat), "+proj=laea +datum=WGS84 +lat_0=5 +lon_0=20")
aflegacy.lambertcord <- aflegacy.lambertcord/1000

#group ID for unique locations
#clustering of points
#distance matrix
locations_PID <- cbind(unique(aflegacy.lambertcord), 1:dim(unique(aflegacy.lambertcord))[1])

N <- dim(aflegacy)[1]
PIDn <- rep(NA, N)
for(i in 1:N){
  PIDn[i] <- locations_PID[locations_PID[,1]==aflegacy.lambertcord[i,1]&locations_PID[,2]==aflegacy.lambertcord[i,2],3]
 }
 
PIDn <- as.numeric(PIDn)

library(fpc)
dbscan.cl <- dbscan(locations_PID[,1:2], 100, 2)
table.dbscan <- table(dbscan.cl$cluster)
d.site<- vector("list", (dim(table.dbscan)[1]-1))
sph.cor20.site <- vector("list",  (dim(table.dbscan)[1]-1))

for(i in 2:dim(table.dbscan)[1]){
	d.site[[(i-1)]] <- as.matrix(dist(locations_PID[dbscan.cl$cluster==(i-1), 1:2]))
	sph.cor20.site[[(i-1)]] <- ((1-3/2*d.site[[(i-1)]]/100+1/2*(d.site[[(i-1)]]/100)^3)*(d.site[[(i-1)]]<100))
}


cluster.asign <- dbscan.cl$cluster[PIDn]

#sort dataset
aflegacy <- cbind(aflegacy, cluster.asign, PIDn)[order(cluster.asign, PIDn), ]
PIDn.new <- rep(NA, N)

PIDn.new[1] <- 1
for(i in 2:N){
	if(aflegacy$PIDn[i]==aflegacy$PIDn[i-1]){
  		PIDn.new[i] <- PIDn.new[i-1]
  	}
  	if(aflegacy$PIDn[i]!=aflegacy$PIDn[i-1]){
  		PIDn.new[i] <- PIDn.new[i-1]+1
  	}
 }
 
logccm <- log(aflegacy$SOC)
N <- length(logccm)
log.depth <- log(aflegacy$Bot)

bio1.scaled <- as.vector(scale(aflegacy$bio1)) 
bio12.scaled<- as.vector(scale(aflegacy$bio12 )) 
CTI_1K.scaled <- as.vector(scale(aflegacy$CTI_1K)) 
ELEV_1K.scaled <- as.vector(scale(aflegacy$ELEV_1K)) 
EVIM_1K.scaled <- as.vector(scale(aflegacy$EVIM_1K)) 
M13RB1ALT.scaled <- as.vector(scale(aflegacy$M13RB1ALT)) 
M13RB2ALT.scaled <- as.vector(scale(aflegacy$M13RB2ALT )) 
M13RB3ALT.scaled <- as.vector(scale(aflegacy$M13RB3ALT)) 
M13RB7ALT.scaled <- as.vector(scale(aflegacy$M13RB7ALT)) 
NPP_Mean.scaled <- as.vector(scale(aflegacy$NPP_Mean_1)) 
RELIEF.scaled <- as.vector(scale(aflegacy$RELIEF)) 
SLOPE_1K.scaled <- as.vector(scale(aflegacy$SLOPE_1K)) 
lstday.scaled <- as.vector(scale(aflegacy$lstday))
lstnight.scaled <- as.vector(scale(aflegacy$lstnight))

Y <- log(aflegacy$SOC)
X <- data.frame(int=1, bio1.scaled=bio1.scaled, bio12=bio12.scaled, CTI_1K.scaled=CTI_1K.scaled, 
              ELEV_1K.scaled=ELEV_1K.scaled, EVIM_1K.scaled=EVIM_1K.scaled, M13RB1ALT.scaled=M13RB1ALT.scaled,
              NPP_Mean.scaled=NPP_Mean.scaled, RELIEF.scaled=RELIEF.scaled, lstday.scaled=lstday.scaled, lstnight.scaled=lstnight.scaled, bio1.scaled.depth=bio1.scaled*log.depth, bio12=bio12.scaled*log.depth, CTI_1K.scaled.depth=CTI_1K.scaled*log.depth, 
              ELEV_1K.scaled.depth=ELEV_1K.scaled*log.depth, EVIM_1K.scaled.depth=EVIM_1K.scaled*log.depth, M13RB1ALT.scaled.depth=M13RB1ALT.scaled*log.depth,
              NPP_Mean.scaled.depth=NPP_Mean.scaled*log.depth, RELIEF.scaled.depth=RELIEF.scaled*log.depth,lstday.scaled.depth=lstday.scaled*log.depth, lstnight.scaled.depth=lstnight.scaled*log.depth, log.depth=log(aflegacy$Bot))
X <- as.matrix(X)
PIDn <- PIDn.new

aflegacy.lambertcord <- project(cbind(aflegacy$Lon, aflegacy$Lat), "+proj=laea +datum=WGS84 +lat_0=5 +lon_0=20")
aflegacy.lambertcord <- aflegacy.lambertcord/1000

#ordered unique locations
locations.est <- cbind(aggregate(aflegacy.lambertcord[,1], by=list(PID = PIDn), mean)[,2], aggregate(aflegacy.lambertcord[,2], by=list(PID = PIDn), mean)[,2])
locations.est.list  <-  vector("list", (length(d.site)+1))
sizes.noise <- table.dbscan[1]

locations.est.list[[1]] <- locations.est[1:sizes.noise, ]
k <- sizes.noise+1
sizes.site <- lapply(d.site, dim)

for(i in 1:(length(d.site))){
	locations.est.list[[(i+1)]] <-locations.est[(k:(k+sizes.site[[i]][1]-1)), ]
	k <- k+(sizes.site[[i]])[1]
}


