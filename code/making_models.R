###### Making models 
## Code by Sara Mortara, first version on 29.04.2019

## how to deal with spatial autocorrelated data: 
# https://stackoverflow.com/questions/53185550/how-to-deal-with-spatially-autocorrelated-residuals-in-glmm
# https://stats.idre.ucla.edu/r/faq/how-do-i-model-a-spatially-autocorrelated-outcome/

## loading packages
library(caret)
# library(geosphere)
# library(geoR)
library(ggplot2)
library(gstat)
library(raster)
library(spaMM)
library(rgdal)
library(rasterVis)


## reading data
data.raw <- read.csv("../data/data.csv", as.is=TRUE)
dim(data.raw)

names(data.raw)[2] <- "Ocean"

data <- data.raw[!duplicated(data.raw[,c("species_name", "dec_lon_new", "dec_lat_new", "Mean")]),]

data$SD <- as.numeric(ifelse(data$SD=="n.g.", NA, data$SD))

summary(data$SD)

head(data)

data$NS <- ifelse(stringr::str_detect(data$Ocean, "North"), "N", "S")

head(data)

## finding a distribution for phlorotanin mean
hist(data$Mean)
hist(log(data$Mean))

# creting vector with log
data$flor <- log(data$Mean+0.01)

summary(data$flor)

names(data)

head(data)

names(data)

env <- data[,c("MS_biogeo06_bathy_slope_5m",
               "BO2_lightbotmax_bdmin",
               "BO2_tempmax_bdmin", "BO2_tempmin_bdmin", "BO2_tempmean_bdmin", "BO2_temprange_bdmin")]

head(env)

env.cor <- cor(na.omit(env))
env.cor
id.cor <- findCorrelation(env.cor)

names(env)[id.cor]

names(data)[names(data)%in%names(env)] <- c("slope", "light", "tempmax", "tempmean", "tempmin", "temprange")



## building models

## 20 km de distancia

# removing missing values
data <- data[!is.na(data$tempmax),]

xy <- data[,c("dec_lon_new", "dec_lat_new")]

tempmin <- raster("../data/tif/BO2_tempmin_bdmin_lonlat.tif")
tempmax <- raster("../data/tif/BO2_tempmax_bdmin_lonlat.tif")
tempran <- raster("../data/tif/BO2_temprange_bdmin_lonlat.tif")
slope <- raster("../data/tif/MS_biogeo06_bathy_slope_5m_lonlat.tif") 
bath <- raster("../data/tif/BO_bathymean_lonlat.tif")
carb <- raster("../data/tif/BO2_carbonphytomean_bdmin_lonlat.tif")
light <- raster("../data/tif/BO2_lightbotmax_bdmin_lonlat.tif")

# visualisation of the data points
plot(tempmax)
points(xy)

temp.vals <- getValues(tempmax)

plot(tempran)
points(xy)


plot(bath)
points(xy)

plot(carb)
points(xy)

plot(slope)
points(xy)

plot(light)
points(xy)


head(data)

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

# ggplot(data=data, aes(x=BO_bathymean, y=Mean, color=Ocean)) +
#   labs(x="Mean temperature at min depth", y="Phlorotannins concentration") +
#   geom_point() + 
#   facet_grid(.~Ocean) +
#   geom_smooth(method='lm') + 
#   theme_classic()

# calculate geographic distance between coordinates
summary(dist(xy))

dist.xy <- distm(xy, fun = distHaversine)

summary(dist.xy[dist.xy!=0])

# we have several sites with distance < 20 km, so we will use a model to take into account spatial autocorrelation
table(dist.xy<20000)


# finding the best model
### finding the best random structure for models
m.full <- fitme(flor ~ slope + tempmax + light +
                       (1|species_name) + (1|site) + (1|NS) +
                       Matern(1|dec_lon_new + dec_lat_new),
                      data=data)

m.rand1 <- fitme(flor ~ slope + tempmax + light +
                  (1|species_name) + (1|site) +
                  Matern(1|dec_lon_new + dec_lat_new),
                  data=data)

m.rand2 <- fitme(flor ~ slope + tempmax + light +
                   (1|species_name) +
                   Matern(1|dec_lon_new + dec_lat_new),
                  data=data)

m.rand3 <- fitme(flor ~ slope + tempmax + light +
                   (1|species_name) + (1|NS) +
                   Matern(1|dec_lon_new + dec_lat_new),
                  data=data)

m.rand4 <- fitme(flor ~ slope + tempmax + light +
                   (1|NS) +
                   Matern(1|dec_lon_new + dec_lat_new),
                   data=data)


rand.sel <- list(m.full, m.rand1, m.rand2, m.rand3, m.rand4)
lapply(rand.sel, extractAIC)

### making models
temp <- data$tempmin

m.full <- fitme(flor ~ slope + temp + light + 
                (1|species_name) + (1|site) + #(1|NS) +
                Matern(1|dec_lon_new + dec_lat_new),
                data=data)

m01 <- fitme(flor ~ slope + temp + light + 
               (1|species_name) + (1|site) + #(1|NS) +
               Matern(1|dec_lon_new + dec_lat_new),
             data=data)

m02 <- fitme(flor ~  slope + temp + 
               (1|species_name) + (1|site) + #(1|NS) +
               Matern(1|dec_lon_new + dec_lat_new),
             data=data)

m03 <- fitme(flor ~ light + temp + 
               (1|species_name) + (1|site) + #(1|NS) +
               Matern(1|dec_lon_new + dec_lat_new),
             data=data)

m04 <- fitme(flor ~  temp +
               (1|species_name) + (1|site) + #(1|NS) +
               Matern(1|dec_lon_new + dec_lat_new),
             data=data)

m.null <- fitme(flor ~  1 + (1|species_name) + (1|site) + #(1|NS) +
               Matern(1|dec_lon_new + dec_lat_new),
              data=data)


m.list <- list(m01, m02, m03, m04, m.null)
names(m.list) <- c("m01", "m02", "m03", "m04","m.null")

AIC.list <- lapply(m.list, extractAIC)
AIC.list2 <- lapply(m.list, AIC)

AICd.vals <- sapply(AIC.list2, function(x) x[3])
AICc.vals <- sapply(AIC.list2, function(x) x[2])
AIC.vals <- sapply(AIC.list, function(x) x[2])

# AIC from extractAIC function
sort(AIC.vals)-sort(AIC.vals)[1]
# dispersal AIC
sort(AICd.vals)-sort(AICd.vals)[1]
# conditional AIC
sort(AICc.vals)-sort(AICc.vals)[1]

LRT(m.null, m04)


### making plot from best model
png("predicted_map.png", res=300, width=2100, height=1200)
filled.mapMM(m09, add.map = TRUE, decorations=NULL)
points(dec_lat_new ~ dec_lon_new, data=data)
dev.off()


filled.mapMM(m09, add.map = TRUE)

# calculating predicted values
xlat <- seq(-180, 180, by=0.1)
ylat <- seq(-70, 70, length.out = length(xlat))

temp <- raster::extract(tempmax, x=xlat, y=ylat)

new.data <- data.frame(dec_lon_new=xlat, dec_lat_new=ylat)

head(new.data)
dim(new.data)

new.data$pred <- predict(m09, re.form=~ Matern(1|dec_lon_new + dec_lat_new), newdata=new.data)


temp

#plot_effects(m09, focal_var="BO2_tempmean_bdmin", ylim=c(0, max(data$flor)))

### models with lme4

library(lme4)
m10 <- lmer(flor ~ Ocean +
              (1|Species) + (1|site), data=data)

# checking autocorrelation in residuals
E <- as.numeric(residuals(m10))
#data$res <- E

mydata <- cbind(E, xy)
coordinates(mydata) <- c("dec_lon_new", "dec_lat_new")

bubble(mydata, "E", col = c("black","grey"),
       main = "Residuals", xlab = "X-coordinates",
       ylab = "Y-coordinates")

Vario1 <- variogram(E ~ 1, mydata)
plot(Vario1)





##### confint


