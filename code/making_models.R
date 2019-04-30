###### Making models 
## Code by Sara Mortara, first version on 29.04.2019

## how to deal with spatial autocorrelated data: 
# https://stackoverflow.com/questions/53185550/how-to-deal-with-spatially-autocorrelated-residuals-in-glmm
# https://stats.idre.ucla.edu/r/faq/how-do-i-model-a-spatially-autocorrelated-outcome/

## reading data
data <- read.csv("../data/data.csv", as.is=TRUE)

head(data)

names(data)[2] <- "Ocean"

## loading packages
library(nlme)
library(bbmle)
library(caret)
library(geosphere)
library(geoR)
library(ggplot2)
library(gstat)


## finding a distribution for phlorotanin mean
hist(data$Mean)
hist(log(data$Mean))

names(data)

head(data)

env.cor <- cor(na.omit(data[,c(25:33)]))
env.cor
findCorrelation(env.cor)

## building models

## 20 km de distancia

# removing missing values
data <- data[!is.na(data$BO2_tempmean_bdmin),]

xy <- data[,c("dec_lon_new", "dec_lat_new")]

temp <- raster("../data/tif/BO2_tempmean_bdmin_lonlat.tif")

# visualisation of the data points
plot(temp)
points(xy)

head(data)

# temp vs phlorotanin concentration
ggplot(data=data, aes(x=BO2_tempmean_bdmin, y=Mean, color=Ocean)) +
  labs(x="Mean temperature at min depth", y="Phlorotannins concentration") +
  geom_point() + 
  facet_grid(.~Ocean) +
  geom_smooth(method='lm') + 
  theme_classic()

# calculate geographic distance between coordinates
summary(dist(xy))

dist.xy <- distm(xy, fun = distHaversine)

summary(dist.xy[dist.xy!=0])


# variogram for phlorotannin concentration
variog(coords=xy, data=data)

# we have several sites with distance < 20 km, so we will use a model to take into account spatial autocorrelation
table(dist.xy<20000)

# checking autocorrelation in residuals
m01 <- lme(Mean ~ BO_bathymean * BO2_tempmean_bdmin, 
             random=list(~ 1|site, ~1|species_name, ~1|Ocean), data=data)


data$res <-as.numeric(residuals(m01))
E <- data$res

mydata <- cbind(E, xy)
coordinates(mydata) <- c("dec_lon_new", "dec_lat_new")

bubble(mydata, "E", col = c("black","grey"),
       main = "Residuals", xlab = "X-coordinates",
       ylab = "Y-coordinates")

Vario1 <- variogram(E ~ 1, mydata)
plot(Vario1)

head(data)

m01sp <- update(m01, correlation = corSpher(form = ~ dec_lon_new + dec_lat_new,
                                            nugget = TRUE), data=data)


#### example


corSpher(form = ~ xy[,1] + xy[,1], nugget = TRUE)
corGaus(form = ~ xy[,1] + xy[,1], nugget = TRUE)
