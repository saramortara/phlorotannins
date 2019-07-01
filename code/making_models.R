###### Making models 
## Code by Sara Mortara, first version on 29.04.2019

# figures with 
# https://kimura.univ-montp2.fr/~rousset/spaMM/example_raster.html



#### 1. loading packages ####
library(caret)
# library(geosphere)
# library(geoR)
library(ggplot2)
library(gstat)
library(raster)
library(spaMM)
library(rgdal)
library(rasterVis)
library(dplyr)
library(vegan)
library(caret)
source("functions.R")

##### 2. reading and organizing data ####
data.raw <- read.csv("../data/data.csv", as.is=TRUE)
dim(data.raw)

names(data.raw)[2] <- "Ocean"

data <- data.raw[!duplicated(data.raw[,c("species_name", "dec_lon_new", "dec_lat_new", "Mean")]),]

head(data)
data$SD <- as.numeric(ifelse(data$SD=="n.g.", NA, data$SD))
data$CV <- data$SD/data$Mean
summary(data$SD)
summary(data$CV...)
unique(data$CV)
summary(data$CV)

head(data)

data$NS <- ifelse(stringr::str_detect(data$Ocean, "North"), "N", "S")

names(data)

## removing missing values
xy <- data[,c("dec_lon_new", "dec_lat_new")]

env <- data[,c("BO_parmax", "BO_salinity", 
         "BO2_tempmax_bdmin", "BO2_tempmin_bdmin", "BO2_tempmean_bdmin", "BO2_temprange_bdmin")] 

#env <- data[,c("MS_biogeo06_bathy_slope_5m",
#               "BO2_lightbotmax_bdmin",
#               "BO2_tempmax_bdmin", "BO2_tempmin_bdmin", "BO2_tempmean_bdmin", "BO2_temprange_bdmin")]

#env.cor <- cor(na.omit(env))
#env.cor
#id.cor <- findCorrelation(env.cor)

names(data)[names(data)%in%names(env)] <- c("PARmax", "salinity", "tempmax", "tempmean", "tempmin", "temprange")

data <- data[!is.na(data$tempmax),]

# checking distribution for phlorotanin mean
hist(data$Mean) # will fit models with Gamma family
hist(data$CV)

summary(data$Mean)
summary(data$CV)


# removing one zero value in data
data <- data[data$Mean>0,]

# separating data per order
data.ord <- list()
ord <- unique(data$Order) 

for(i in 1:length(ord)){
  data.ord[[i]] <- data[data$Order==ord[i],] 
}

ord

names(data.ord) <- ord

length(data.ord)

n.ord <- sapply(data.ord, nrow)
id.ord <- n.ord>7

data.or <- data.ord[id.ord]
length(data.or)

other <- data.ord[!id.ord] 
data.ot <- dplyr::bind_rows(other)

head(data.ot)

data.all <- c(data.or, Other=list(data.ot))

names(data.all)
length(data.all)


# creating columns with order
other_orders <- setdiff(ord, names(data.all))

names(data)

data$Order2 <- ifelse(data$Order%in%other_orders, "Other", data$Order)

table(data$Order2)
table(data$Order)

##### 3. visualizing data ####

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
ggplot(data=data, aes(x=tempmax, y=Mean, color=Order2)) +
  labs(x="Maximum temperature at min depth", y="Phlorotannins concentration") +
  geom_point() + 
  facet_grid(.~Order2) +
  #geom_smooth(method='lm') + 
  theme_classic()


ggplot(data=data, aes(x=tempmax, y=CV, color=Order2)) +
  labs(x="Mean carbon concentration at min depth", y="Phlorotannins concentration") +
  geom_point() + 
  facet_grid(.~Order2) +
  #geom_smooth(method='lm') + 
  theme_classic()

# boxplot
# temp vs phlorotanin concentration
ggplot(data=data, aes(x=Order2, y=Mean, fill=Ocean)) +
  labs(x="Order", y="Phlorotannins concentration") +
  geom_boxplot() + 
  #facet_grid(.~Order2 + Ocean) +
  #geom_smooth(method='lm') + 
  theme_classic()

data$Order2 <- as.factor(data$Order2)
my_xlab <- paste(levels(data$Order2),"\n(N=",table(data$Order2),")",sep="")

p_mean <- ggplot(data=data, aes(x=Order2, y=Mean)) +
  labs(x="Order", y="Phlorotannins concentration") +
  geom_boxplot(varwidth=TRUE) + 
  theme_classic() +
  scale_x_discrete(labels=my_xlab) +
  scale_y_log10()

p_cv <- ggplot(data=data, aes(x=Order2, y=CV)) +
  labs(x="Order", y="Phlorotannins concentration") +
  geom_boxplot(varwidth=TRUE) + 
  theme_classic() +
  scale_x_discrete(labels=my_xlab) +
  scale_y_log10()

png("../results/concentration_vs_order.png", res=300, 
    width=1400, height=1200)




dev.off()

ggplot(data=data, aes(x=Ocean, y=Mean)) +
  labs(x="Ocean", y="Phlorotannins concentration") +
  geom_boxplot() + 
  theme_classic()


#### 4. Fitting models ####
# excluding NA to run models with CV
data.cv <- lapply(data.all, na.omit)

### one order as a test
m.list <- my.models(data.all[[1]], var.y="Mean", var.t="tempmax")
m.list2 <- my.models(data.all[[1]], var.y="CV", var.t="temprange")

### 4.1 making models separated per order

# models for mean
models.mean.or <- list()
for(i in 1:length(data.all)){
message(paste("running mean models for", names(data.all)[i]))
models.mean.or[[i]] <- my.models(data.all[[i]], var.y="Mean", var.t="tempmax")
#message(paste("running CV models for", names(data.all)[i]))
#models.cv[[i]] <- my.models(data.cv[[i]], "CV") #using data set without NAs
}

# models for CV
models.cv.or <- list()
for(i in 1:length(data.all)){
message(paste("running CV models for", names(data.all)[i]))
models.cv.or[[i]] <- my.models(data.cv[[i]], var.y="CV", var.t="temprange") #using data set without NAs
}

models.cv.or2 <- list()
for(i in 1:length(data.all)){
  message(paste("running CV models for", names(data.all)[i]))
  models.cv.or2[[i]] <- my.models(data.cv[[i]], var.y="CV", var.t="tempmax") #using data set without NAs
}

models.cv.or2

## aic per order 
all.aic <- list()
for(i in 1:length(data.all)){
  all.aic[[i]] <- my.aic(models.mean.or[[i]], or=names(data.all)[i])
 # all.aic.cv[[i]] <- my.aic(models.cv.or[[i]], or=names(data.all)[i])
  }
all.aic <- bind_rows(all.aic)

all.aic.cv <- list()
for(i in 1:3){
  all.aic.cv[[i]] <- my.aic(models.cv.or[[i]], or=names(data.all)[i])
  # all.aic.cv[[i]] <- my.aic(models.cv.or[[i]], or=names(data.all)[i])
}


all.aic.cv <- bind_rows(all.aic.cv)

best.aic <- all.aic[all.aic$dAIC<2,]
best.aic.cv <- all.aic.cv[all.aic.cv$dAIC<2,]

best.aic
best.aic.cv

### 4.2. All orders together

m.full <- fitme(Mean ~ tempmax * Ocean + Order2 +
                  (1|species_name) + (1|Reference) +
                  Matern(1|dec_lon_new + dec_lat_new),
                  family=Gamma(log),
                  data=data)

m.01 <- fitme(Mean ~ tempmax * Ocean + Order2 +
                  (1|species_name) + (1|Reference) +
                  Matern(1|dec_lon_new + dec_lat_new),
                family=Gamma(log),
                data=data)

m.null <- fitme(Mean ~ 1 + 
                  (1|species_name) + (1|Reference) +
                  Matern(1|dec_lon_new + dec_lat_new),
                 family=Gamma(log),
                 data=data)

extractAIC(m.full)
extractAIC(m.null)

## fitting models
models.cv <- my.models.all(data, "CV", "temprange")
models.mean <- my.models.all(data, "Mean", "tempmax")
## aic table
my.aic(models.cv, "all")
my.aic(models.mean, "all")

#AIC.list2 <- lapply(m.list, AIC)
#AICd.vals <- sapply(AIC.list2, function(x) x[3])
#AICc.vals <- sapply(AIC.list2, function(x) x[2])

#write.table(all.tab, "../results/aic_tab_mean.csv", row.names=FALSE, col.names=TRUE, sep=",")

#### 5. Calculating predictors ####

# extracting values from model

### using all temperature values to make the prediction map
temp.xy <- rasterToPoints(tempmax)

temp.xy <- as.data.frame(temp.xy)

head(data.all[[1]])

## extracting all temperature values
temp.vals <- na.omit(temp.xy) 
head(temp.vals)
temp.pred <- data.frame(tempmax=temp.vals$BO2_tempmax_bdmin_lonlat, 
                        dec_lon_new=temp.xy$x, 
                        dec_lat_new=temp.xy$y)


id.temp <- names(data.all)%in%c("Fucales", "Laminariales", "Other") 
data.temp <- data.all[id.temp]
models.temp <- models.mean.or[id.temp]


pred <- matrix(NA, nrow=nrow(temp.pred), ncol=length(data.temp))
for(i in 1:length(data.temp)){
pred[,i] <- predict(models.temp[[i]]$'Temp*Ocean', newdata=temp.pred, re.form=~Matern(1|dec_lon_new + dec_lat_new))
}
## data frame with all cell values
temp.xy$pred <- NA
temp.xy$pred[!is.na(temp.xy$BO2_tempmax_bdmin_lonlat)] <- pred

head(temp.xy)
}


## create raster from predictors 
raster.predict <- CreateRaster(long=temp.xy$x, 
                               lat=temp.xy$y,
                               values=exp(temp.xy$pred), 
                               save.spatial.files = TRUE,
                               filename = "phlorotannin_predict")

#### 6. Making plots ####

### default plot from package

filled.mapMM(models.temp[[1]]$`Temp*Ocean`, add.map = TRUE)




### 6.2 making plot

## calculating predict and intervals 
best.mean <- lapply(models.mean.or, best.mod) 
pred.mean <- list()
for(i in 1:length(data.all)){
pred.mean[[i]] <- pred.table(best.mean[[i]], data.all[[i]])
}

head(temp.xy)

png("temp_vs_phlor.png", res = 300, width=1500, height = 1200)
ggplot(data=data, aes(x=tempmax, y=Mean)) +
  labs(x="Maximum temperature at min depth (°C)", y="Phlorotannin concentration (Log)") +
  geom_point(alpha=0.7) + 
  geom_ribbon(aes(x=temp, y=pred, ymin=lwr, ymax=upr), data=pred.table, alpha=0.1) +
  facet_grid(.~Ocean) +
  #geom_smooth(method='lm') + 
  theme_classic()
dev.off()



### 6.1 making map 
data(worldcountries)
projection="+proj=longlat +datum=WGS84"
worldcountries <- spTransform(worldcountries, projection)

country.layer <- layer(sp.polygons(worldcountries, fill=fill),
                       data=list(sp.polygons=sp.polygons, 
                                 worldcountries=worldcountries, 
                                 fill="transparent"))

coordinates(xy) <- ~dec_lon_new+dec_lat_new

palette <- wesanderson::wes_palette("Zissou1", 16, type="continuous")
myTheme <- rasterTheme(region=palette)

png("predicted_map.png", res=300, width=2100, height=1200)
levelplot(raster.predict, margin=FALSE)+
  #cuts=length(palette)-1,
  #col.regions=palette) + #, par.settings=myTheme
  #latticeExtra::layer(sp.points(xy,cex=1.3, col=1, pch=1)) +
  country.layer
dev.off()

summary(models.mean.all$Temp)



#### exporting S1 table

head(data)

stable <- data[,c("Ocean", "dec_lon_new", "dec_lat_new", 'species_name', "Order", "Mean", "CV", "Reference")]

stable$species_name <- gsub("_", " ", stable$species_name)

names(stable)[2:4] <- c("Longitude", "Latitude", "Scientific name")

head(stable)
dim(stable)

write.table(stable, "../results/S_table_01.csv", sep=",", 
            col.names = TRUE, 
            row.names=FALSE)
