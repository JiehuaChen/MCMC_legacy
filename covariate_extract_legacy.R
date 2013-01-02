library(rgdal)

afhori <- read.table("Horizon.csv", header=T, sep=",", stringsAsFactors=FALSE)
afprof <- read.table("Profile.csv", header=T, sep=",")

aflegacy.lambertcord <- project(cbind(afprof$Lon, afprof$Lat), "+proj=laea +datum=WGS84 +lat_0=5 +lon_0=20")
 mapnames <- c( "bio1", "bio12", "CTI_1K", "ELEV_1K", "EVI_mask_1K", "EVIM_1K", "M13RB1ALT", "M13RB2ALT", "M13RB3ALT", "M13RB7ALT","MASK",  "NPP_Mean_1", "RELIEF_1K", "SLOPE_1K", "lstday", "lstnight")
#covar.selected.map <-  c(1, 2, 3, 4, 5, 6, 7, 12, 13, 15, 16)

tiffile <- "predictgrid.tif"
tif.info <- GDALinfo(tiffile, silent=TRUE)
		totalcols <- tif.info[2]
		totalrows <- tif.info[1]
		origin.x <- tif.info[4]
		res.x <- tif.info[6]
		res.y <- tif.info[7]
		if(attr(tif.info, "ysign")==-1){
		#origin of data array should be the left upper corner
			origin.y <- tif.info[5] + totalrows*res.y
		}else{
		origin.y <- tif.info[5]	
		}


cov.extract <- rep(NA, length(mapnames))
for(i in 1:dim(aflegacy.lambertcord)[1]){
	pixeln.x <- ((aflegacy.lambertcord[i,1]-origin.x)/res.x)
	pixeln.y <- ((origin.y-aflegacy.lambertcord[i,2])/res.y)
	predcov.values <- readGDAL(tiffile, offset=c(pixeln.y, pixeln.x), region.dim=c(1,1), silent=TRUE)@data
	cov.extract <- rbind(cov.extract, predcov.values)
	print(i)
}

cov.extract <- cov.extract[-1,]
cov.extract <- as.data.frame(cov.extract)
names(cov.extract) <- mapnames
write.csv(cov.extract, "legacy_profile_withcov.csv", row.names=FALSE, quote=FALSE)
