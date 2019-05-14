###### Making models 
## Code by Sara Mortara, first version on 29.04.2019

# figures with 
# https://kimura.univ-montp2.fr/~rousset/spaMM/example_raster.html

# loading packages
library(caret)
# library(geosphere)
# library(geoR)
library(ggplot2)
library(gstat)
library(raster)
library(spaMM)
library(rgdal)
library(rasterVis)

# reading data
data.raw <- read.csv("../data/data.csv", as.is=TRUE)
dim(data.raw)

names(data.raw)[2] <- "Ocean"

data <- data.raw[!duplicated(data.raw[,c("species_name", "dec_lon_new", "dec_lat_new", "Mean")]),]

data$SD <- as.numeric(ifelse(data$SD=="n.g.", NA, data$SD))

summary(data$SD)

head(data)

data$NS <- ifelse(stringr::str_detect(data$Ocean, "North"), "N", "S")

## removing missing values
xy <- data[,c("dec_lon_new", "dec_lat_new")]

env <- data[,c("MS_biogeo06_bathy_slope_5m",
               "BO2_lightbotmax_bdmin",
               "BO2_tempmax_bdmin", "BO2_tempmin_bdmin", "BO2_tempmean_bdmin", "BO2_temprange_bdmin")]

names(data)[names(data)%in%names(env)] <- c("slope", "light", "tempmax", "tempmean", "tempmin", "temprange")

data <- data[!is.na(data$tempmax),]

head(data)

# checking distribution for phlorotanin mean
hist(data$Mean)
hist(log(data$Mean)) # will fit models with log()

# creting vector with log to use in the models
data$flor <- log(data$Mean+0.01)

summary(data$flor)

# env.cor <- cor(na.omit(env))
# env.cor
# id.cor <- findCorrelation(env.cor)


# saving object with raster
tempmax <- raster("../data/tif/BO2_tempmax_bdmin_lonlat.tif")
# tempmin <- raster("../data/tif/BO2_tempmin_bdmin_lonlat.tif")
# tempran <- raster("../data/tif/BO2_temprange_bdmin_lonlat.tif")
# slope <- raster("../data/tif/MS_biogeo06_bathy_slope_5m_lonlat.tif") 
# bath <- raster("../data/tif/BO_bathymean_lonlat.tif")
# carb <- raster("../data/tif/BO2_carbonphytomean_bdmin_lonlat.tif")
# light <- raster("../data/tif/BO2_lightbotmax_bdmin_lonlat.tif")

# visualisation of the data points
plot(tempmax)
points(xy)

# temp vs phlorotanin concentration
ggplot(data=data, aes(x=tempmax, y=flor, color=Ocean)) +
  labs(x="Mean temperature at min depth", y="Phlorotannins concentration") +
  geom_point() + 
  facet_grid(.~Ocean) +
  #geom_smooth(method='lm') + 
  theme_classic()


ggplot(data=data, aes(x=slope, y=flor, color=Ocean)) +
  labs(x="Mean carbon concentration at min depth", y="Phlorotannins concentration") +
  geom_point() + 
  facet_grid(.~Ocean) +
  #geom_smooth(method='lm') + 
  theme_classic()

####  Selection of the random structure ####

# finding the best random structure for models
m.rand1 <- fitme(flor ~ tempmax*Ocean +
                  (1|species_name) + (1|site) + (1|Reference) +
                  Matern(1|dec_lon_new + dec_lat_new),
                data=data)


m.rand2 <- fitme(flor ~ tempmax*Ocean +
                  (1|species_name) + (1|site) +
                  Matern(1|dec_lon_new + dec_lat_new),
                  data=data)

m.rand3 <- fitme(flor ~ tempmax*Ocean +
                   (1|species_name) +
                   Matern(1|dec_lon_new + dec_lat_new),
                  data=data)

m.rand4 <- fitme(flor ~ tempmax*Ocean +
                   (1|species_name) + (1|Reference) +
                   Matern(1|dec_lon_new + dec_lat_new),
                  data=data)

m.rand5 <- fitme(flor ~ tempmax*Ocean +
                   (1|Reference) +
                   Matern(1|dec_lon_new + dec_lat_new),
                   data=data)

m.rand6 <- fitme(flor ~ tempmax*Ocean +
                   (1|site) +
                   Matern(1|dec_lon_new + dec_lat_new),
                 data=data)

m.rand7 <- fitme(flor ~ tempmax*Ocean +
                   Matern(1|dec_lon_new + dec_lat_new),
                 data=data)

rand.sel <- list(m.rand1, m.rand2, m.rand3, m.rand4, m.rand5, m.rand6, m.rand7)
lapply(rand.sel, extractAIC)

### making models

m01 <- fitme(flor ~ tempmax*Ocean +
               (1|species_name) + (1|Reference) +
               Matern(1|dec_lon_new + dec_lat_new),
             data=data)
m02 <- fitme(flor ~  tempmax+Ocean +
               (1|species_name) + (1|Reference) +
               Matern(1|dec_lon_new + dec_lat_new),
             data=data)
m03 <- fitme(flor ~ tempmax +
               (1|species_name) + (1|Reference) +
               Matern(1|dec_lon_new + dec_lat_new),
             data=data)


m03b <- corrHLfit(flor ~ tempmax +
               (1|species_name) + (1|Reference) +
               Matern(1|dec_lon_new + dec_lat_new),
             data=data, HLmethod="ML")

ci.b <- confint(m03b, "tempmax")
summary(m03b)

m04 <- fitme(flor ~  Ocean +
               (1|species_name)  + (1|Reference) +
               Matern(1|dec_lon_new + dec_lat_new),
             data=data)
m.null <- fitme(flor ~  1 + (1|species_name) + (1|Reference) +
                  Matern(1|dec_lon_new + dec_lat_new),
                data=data)

m.list <- list(m01, m02, m03, m04, m.null)
names(m.list) <- c("m01", "m02", "m03", "m04", "m.null")

AIC.list <- lapply(m.list, extractAIC)
AIC.list2 <- lapply(m.list, AIC)

AICd.vals <- sapply(AIC.list2, function(x) x[3])
AICc.vals <- sapply(AIC.list2, function(x) x[2])
AIC.vals <- data.frame(AIC=sapply(AIC.list, function(x) x[2]))

AIC.vals
AIC.vals$df <- sapply(AIC.list, function(x) x[1])
AIC.vals$delta <- AIC.vals$AIC-sort(AIC.vals$AIC)[1]

AIC.vals$model <- c("Temp:ocean + Temp + Ocean", "Temp + Ocean", "Temp", "Ocean", "Null")

AICw <- function(x){ exp(-0.5*x)/sum(exp(-0.5*AIC.vals$delta)) }

AIC.vals$weights <- round(AICw(AIC.vals$delta), 3)
AIC.vals$weights[AIC.vals$weights==0] <- "<0.001"

aic.tab <- AIC.vals[order(AIC.vals$delta),c("model", "delta", "df", "weights", "AIC")]

write.table(aic.tab, "aic_tab.csv", row.names=FALSE, col.names=TRUE, sep=",")

# AIC from extractAIC function
sort(AIC.vals)-sort(AIC.vals)[1]
# dispersal AIC
sort(AICd.vals)-sort(AICd.vals)[1]
# conditional AIC
sort(AICc.vals)-sort(AICc.vals)[1]

LRT(m01, m.null)
LRT(m02, m.null)
LRT(m03, m.null)
LRT(m04, m.null)

summary(m03)


# extracting values from model

### using all temperature values to make the prediction map
temp.xy <- rasterToPoints(tempmax)

temp.xy <- as.data.frame(temp.xy)

## extracting all temperature values
temp.vals <- na.omit(temp.xy) 
head(temp.vals)
temp.pred <- data.frame(tempmax=temp.vals$BO2_tempmax_bdmin_lonlat, 
                        dec_lon_new=temp.xy$x, 
                        dec_lat_new=temp.xy$y)
## calculating predictors
pred <- predict(m03, newdata=temp.pred, re.form=~Matern(1|dec_lon_new + dec_lat_new))
## data frame with all cell values
temp.xy$pred <- NA
temp.xy$pred[!is.na(temp.xy$BO2_tempmax_bdmin_lonlat)] <- pred

head(temp.xy)

## calculating intervals
#int <- get_intervals(m03, newdata=temp.pred, re.form=~Matern(1|dec_lon_new + dec_lat_new))

int <- get_intervals(m03, re.form=NA)
pred2 <- predict(m03, re.form=NA)
head(int)

pred.table <- data.frame(pred=pred2, lwr=int[,1], upr=int[,2], temp=data$tempmax)
head(pred.table)

## create raster from predictors 

CreateRaster <- function(
  long,
  lat,
  values,
  proj="+proj=longlat +datum=WGS84",
  save.spatial.files=FALSE,
  filename="data_raster",
  overwrite.spatial.files=TRUE
) {
  # Args:
  #   long: a vector of the longitudes of the raster cells
  #   lat: a vector of the latitudes of the raster cells
  #   values: a vector of the values of the raster cells
  #   proj: the projection system for the raster
  #   save.spatial.files: logical indicating if 
  #          an hard copy of the raster should be saved (as ascii)
  #   filename: name of the file for the hard copy
  #   overwrite.spatial.files: logical indicating if 
  #          an existing hard copy should be overwritten or not
  #
  # Returns:
  #   The raster.
  #
  data <- data.frame(long=long, lat=lat, values=values)  # a dataframe 
  #                 with longitudes, lattitudes, and values is being created
  coordinates(data) <- ~long+lat  # coordinates are being set for the raster
  proj4string(data) <- CRS(proj)  # projection is being set for the raster
  gridded(data) <- TRUE  # a gridded structure is being set for the raster
  data.raster <- raster(data)  # the raster is being created
  if(save.spatial.files) writeRaster(
    data.raster,
    filename=paste(filename, ".asc", sep=""),
    overwrite=overwrite.spatial.files
  )  # if save=TRUE the raster is exported as an ascii file
  return(data.raster) # the raster is being returned
}

raster.predict <- CreateRaster(long=temp.xy$x, 
                               lat=temp.xy$y,
                               values=exp(temp.xy$pred), 
                               save.spatial.files = TRUE,
                               filename = "phlorotannin_predict")
### making map 
data(worldcountries)
projection="+proj=longlat +datum=WGS84"
worldcountries <- spTransform(worldcountries, projection)

country.layer <- layer(sp.polygons(worldcountries, fill=fill),
                       data=list(sp.polygons=sp.polygons, 
                                 worldcountries=worldcountries, 
                                 fill="transparent"))

coordinates(xy) <- ~dec_lon_new+dec_lat_new

palette <- spaMM.colors() #wesanderson::wes_palette("Zissou1", 16, type="continuous")
myTheme <- rasterTheme(region=paleta)

png("predicted_map.png", res=300, width=2100, height=1200)
levelplot(raster.predict, margin=FALSE,  
          cuts=length(palette)-1,
          col.regions=palette) + #, par.settings=myTheme
#latticeExtra::layer(sp.points(xy,cex=1.3, col=1, pch=1)) +
  country.layer
dev.off()

summary(m03)

# making plot

head(temp.xy)

png("temp_vs_phlor.png", res = 300, width=1500, height = 1200)
ggplot(data=data, aes(x=tempmax, y=flor)) +
  labs(x="Maximum temperature at min depth (°C)", y="Phlorotannin concentration (Log)") +
  geom_point(alpha=0.7) + 
  geom_ribbon(aes(x=temp, y=pred, ymin=lwr, ymax=upr), data=pred.table, alpha=0.1) +
  #facet_grid(.~Ocean) +
  #geom_smooth(method='lm') + 
  theme_classic()
dev.off()


