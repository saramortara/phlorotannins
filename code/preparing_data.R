#### Code for analysis of global phlorotannin concentration ####
## First version by Sara Mortara 11.04.2019

# loading packages
library(stringr)
library(measurements)
library(sdmpredictors)
#library("leaflet")
library(raster)


# reading data
data <- read.csv("../data/global_table.csv", sep=";", as.is=TRUE, fill=TRUE, encoding = "UTF-8")[,-10]

head(data)

# cleaning coordinates columns
data$dec_lat <- gsub('N', '', data$Latitude)
data$dec_lat <- gsub('S', '', data$dec_lat)
data$dec_lon <- gsub('W', '', data$Longitude)
data$dec_lon <- gsub('E', '', data$dec_lon)

# converting coordinates to decimal degree
data$dec_lat <- gsub('°', ' ', data$dec_lat)
data$dec_lat <- gsub("'", '', data$dec_lat)
data$dec_lat <- gsub("’", '', data$dec_lat)
data$dec_lat <- gsub('"', '', data$dec_lat)
data$dec_lon <- gsub('°', ' ', data$dec_lon)
data$dec_lon <- gsub("'", '', data$dec_lon)
data$dec_lon <- gsub("’", '', data$dec_lon)
data$dec_lon <- gsub('"', '', data$dec_lon)

data$dec_lon_deg <- sapply(strsplit(data$dec_lon, " "), function(x)x[1])
data$dec_lon_min <- sapply(strsplit(data$dec_lon, " "), function(x)x[2])
data$dec_lat_deg <- sapply(strsplit(data$dec_lat, " "), function(x)x[1])
data$dec_lat_min <- sapply(strsplit(data$dec_lat, " "), function(x)x[2])

data$dec_lon_min[is.na(data$dec_lon_min)] <- "00"
data$dec_lat_min[is.na(data$dec_lat_min)] <- "00"

# creating vector with coordinates in the correct format for conversion
data$dec_lon_paste <- paste(data$dec_lon_deg, data$dec_lon_min)
data$dec_lat_paste <- paste(data$dec_lat_deg, data$dec_lat_min)

head(data)
tail(data)
  
data$dec_lon_new <- as.numeric(measurements::conv_unit(data$dec_lon_paste, from = 'deg_dec_min', to = 'dec_deg'))
data$dec_lat_new <- as.numeric(measurements::conv_unit(data$dec_lat_paste, from = 'deg_dec_min', to = 'dec_deg'))

head(data)
tail(data)

summary(data)

# adding negative signal
data$dec_lat_new <- ifelse(str_detect(data$Latitude, "S"), 
                       data$dec_lat_new*(-1), data$dec_lat_new*1)

data$dec_lon_new <- ifelse(str_detect(data$Longitude, "W"), 
                       data$dec_lon_new*-1, data$dec_lon_new*1)

head(data)
summary(data)

# creating column with species name only 
data$genus <- sapply(str_split(data$Species, "( )"), function(x) x[1])
data$epithet <- sapply(str_split(data$Species, "( )"), function(x) x[2])
data$species_name <- paste(data$genus, data$epithet, sep="_")
head(data)


# creating id column

data$site.id <- paste(data$dec_lon_new, data$dec_lat_new, sep='_')
head(data)

site <- data.frame(site.id=unique(data$site.id), 
                   site=paste("site", sprintf("%03d", 1:length(unique(data$site.id))), sep="_"))

data.site <- merge(data, site, by='site.id')


head(data.site)

# getting environmental variables
df <- list_layers()
mar <- df[df$marine==TRUE,]

# mar$layer_code
# mar$name

# creating object with names of environmental variables
var.names <- c("Bathymetry (mean)",   "Sea water temperature (mean at max depth)" , 
               "Sea water temperature (mean at mean depth)" ,  "Sea water temperature (mean at min depth)", 
               "Sea surface temperature (mean)", "Primary production (mean)",  "Carbon phytoplankton biomass (mean at max depth)", 
               "Carbon phytoplankton biomass (mean at mean depth)", "Carbon phytoplankton biomass (mean at min depth)")    


var.code <- mar$layer_code[mar$name%in%var.names]
var.code

# Download environmental data layers in var.code object
env.layers <- lapply(var.code, load_layers, equalarea=FALSE, rasterstack=TRUE, datadir="../data/tif/")

#env.layers <- load_layers(layercodes = var.code, 
#                           equalarea=FALSE, rasterstack=TRUE, datadir="../data/tif/") 

# creating object with environmental variables for each coordinate
env <- matrix(NA, ncol=length(env.layers), nrow=nrow(data))
colnames(env) <- var.code
for(i in 1:length(env.layers)){
env[,i] <- extract(env.layers[[i]], data[,c("dec_lon_new", "dec_lat_new")])
}

head(env)

data.env <- cbind(data.site, env)

head(data.env)

#### exporting table for analysis
write.table(data.env, "../data/data.csv", sep = ',', 
            row.names = FALSE, col.names=TRUE, dec=".")
