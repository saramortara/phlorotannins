my.models <- function(data, var.y, var.t){
  data$y <- data[,var.y]
  data$temp <- data[,var.t]
  data <- data[data$y>0,]
  mod.nomes <- c("Temp*Ocean", "Ocean", #biogeo
                 "Sal+PAR", #gradient
                 "Temp*Ocean+Sal+PAR", "Ocean+Sal+PAR", #combined
                 "Null")
   # message(paste("running", mod.nomes[1], "..."))
   # # biogeographical hypothesis
   # m01 <- fitme(y ~ temp +
   #                (1|species_name) + (1|Reference) +
   #                Matern(1|dec_lon_new + dec_lat_new),
   #              family=Gamma(log),
   #              data=data)
  message(paste("running", mod.nomes[1], "..."))
  m02 <- fitme(y ~ temp*Ocean +
                 (1|species_name) + (1|Reference) +
                 Matern(1|dec_lon_new + dec_lat_new),
               family=Gamma(log),
               data=data)
  message(paste("running", mod.nomes[2], "..."))
  m03 <- fitme(y ~ Ocean +
                 (1|species_name) + (1|Reference) +
                 Matern(1|dec_lon_new + dec_lat_new),
               family=Gamma(log),
               data=data)
  # environmental gradient hypothesis
  message(paste("running", mod.nomes[3], "..."))
  m04 <- fitme(y ~  salinity+PARmax +
                 (1|species_name) + (1|Reference) +
                 Matern(1|dec_lon_new + dec_lat_new),
               family=Gamma(log),
               data=data)
  # biogeographical and environmental gradient hypothesis
  message(paste("running", mod.nomes[4], "..."))
  m05 <- fitme(y ~ temp*Ocean+salinity+PARmax +
                 (1|species_name) + (1|Reference) +
                 Matern(1|dec_lon_new + dec_lat_new),
               family=Gamma(log),
               data=data)
  # message(paste("running", mod.nomes[6], "..."))
  # m06 <- fitme(y ~ temp+salinity+PARmax +
  #                (1|species_name) + (1|Reference) +
  #                Matern(1|dec_lon_new + dec_lat_new),
  #              family=Gamma(log),
  #              data=data)
  message(paste("running", mod.nomes[5], "..."))
  m07 <- fitme(y ~ Ocean+salinity+PARmax +
                 (1|species_name) + (1|Reference) +
                 Matern(1|dec_lon_new + dec_lat_new),
               family=Gamma(log),
               data=data)
  message(paste("running", mod.nomes[6], "!"))
  m.null <- fitme(y ~  1 + (1|species_name) + (1|Reference) +
                    Matern(1|dec_lon_new + dec_lat_new),
                  family=Gamma(log),
                  data=data)
  m.list <- list(m02, m03, m04, m05, m07, m.null)
  names(m.list) <- mod.nomes
  return(m.list)
}

## function to operform model selection
my.models.all <- function(data, var.y, var.t){
  data$y <- data[,var.y]
  data$temp <- data[,var.t]
  data <- data[data$y>0,]
  mod.nomes <- c("Temp", "Temp+Ocean", "Sal+PAR", "Temp+Ocean+Sal+PAR","Temp+Sal+PAR", "Null")
  message(paste("running", mod.nomes[1], "..."))
  m00 <- fitme(y ~ temp +
                 (1|Order) + (1|species_name) + (1|Reference) +
                 Matern(1|dec_lon_new + dec_lat_new),
               family=Gamma(log),
               data=data)
  message(paste("running", mod.nomes[2], "..."))
  m01 <- fitme(y ~ temp+Ocean +
                 (1|Order) + (1|species_name) + (1|Reference) +
                 Matern(1|dec_lon_new + dec_lat_new),
               family=Gamma(log),
               data=data)
  message(paste("running", mod.nomes[3]), "...")
  m02 <- fitme(y ~  salinity+PARmax +
                 (1|Order) + (1|species_name) + (1|Reference) +
                 Matern(1|dec_lon_new + dec_lat_new),
               family=Gamma(log),
               data=data)
  message(paste("running", mod.nomes[4]), "...")
  m03 <- fitme(y ~ temp + Ocean+salinity+PARmax +
                 (1|Order) + (1|species_name) + (1|Reference) +
                 Matern(1|dec_lon_new + dec_lat_new),
               family=Gamma(log),
               data=data)
  message(paste("running", mod.nomes[5]), "...")
  m04 <- fitme(y ~ temp+salinity+PARmax +
                 (1|Order) + (1|species_name) + (1|Reference) +
                 Matern(1|dec_lon_new + dec_lat_new),
               family=Gamma(log),
               data=data)
  message(paste("running", mod.nomes[6]), "!")
  m.null <- fitme(y ~  1 + (1|species_name) + (1|Reference) +
                    (1|Order) +  Matern(1|dec_lon_new + dec_lat_new),
                  family=Gamma(log),
                  data=data)
  m.list <- list(m00, m01, m02, m03, m04, m.null)
  names(m.list) <- mod.nomes
  return(m.list)
}


my.models.or <- function(data, var){
  data$y <- data[,var]
  data <- data[data$y>0,]
  mod.nomes <- c("Temp+Ocean", "Temp", "Ocean", "Null")
  message(paste("running", mod.nomes[2]), "...")
  m02 <- fitme(y ~  tempmax+Ocean + Order +
                 (1|species_name) + (1|Reference) +
                 Matern(1|dec_lon_new + dec_lat_new),
               family=Gamma(log),
               data=data)
  message(paste("running", mod.nomes[3]), "...")
  m03 <- fitme(y ~ tempmax + Order +
                 (1|species_name) + (1|Reference) +
                 Matern(1|dec_lon_new + dec_lat_new),
               family=Gamma(log),
               data=data)
  message(paste("running", mod.nomes[4]), "...")
  m04 <- fitme(y ~  Ocean + Order +
                 (1|species_name)  + (1|Reference) +
                 Matern(1|dec_lon_new + dec_lat_new),
               family=Gamma(log),
               data=data)
  message(paste("running", mod.nomes[5]), "!")
  m.null <- fitme(y ~  1 + (1|species_name) + (1|Reference) +
                    Matern(1|dec_lon_new + dec_lat_new),
                  family=Gamma(log),
                  data=data)
  m.list <- list(m01, m02, m03, m04, m.null)
  names(m.list) <- mod.nomes
  return(m.list)
}

#### function to calculate AIC
my.aic <- function(m.list, or){
  AIC.list <- lapply(m.list, extractAIC)
  AIC.vals <- data.frame(Order=or, AIC=sapply(AIC.list, function(x) x[2]))
  AIC.vals$df <- sapply(AIC.list, function(x) x[1])
  AIC.vals$dAIC <- AIC.vals$AIC-sort(AIC.vals$AIC)[1]
  AIC.vals$Model <- names(m.list)
  AICw <- function(x){exp(-0.5*x)/sum(exp(-0.5*AIC.vals$dAIC))}
  AIC.vals$Weights <- round(AICw(AIC.vals$dAIC), 3)
  AIC.vals$Weights[AIC.vals$Weights==0] <- "<0.001"
  aic.tab <- AIC.vals[order(AIC.vals$dAIC),c("Order", "Model", "dAIC", "df", "Weights", "AIC")]
  row.names(aic.tab) <- NULL
  return(aic.tab)
}

# Create raster from predictors
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


# function to select best model
best.mod <- function(modlist){
AIC.list <- lapply(modlist, extractAIC)
AICvals <- sapply(AIC.list, function(x) x[2])
best <- modlist[which.min(AICvals)]
return(best)
}  



# function to generate predict table for each model
pred.table <- function(mod, data){
int <- get_intervals(mod, re.form=NA)
pred2 <- predict(mod, re.form=NA)
pred.table <- data.frame(pred=pred2, lwr=int[,1], upr=int[,2], temp=data$tempmax)
return(pred.table)
}

