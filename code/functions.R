## function to operform model selection
my.models <- function(data, var){
  data$y <- data[,var]
  data <- data[data$y>0,]
  mod.nomes <- c("Temp*Ocean", "Temp+Ocean", "Temp", "Ocean", "Null")
  message(paste("running", mod.nomes[1], "..."))
  m01 <- fitme(y ~ tempmax*Ocean +
                 (1|species_name) + (1|Reference) +
                 Matern(1|dec_lon_new + dec_lat_new),
               family=Gamma(log),
               data=data)
  message(paste("running", mod.nomes[2]), "...")
  m02 <- fitme(y ~  tempmax+Ocean +
                 (1|species_name) + (1|Reference) +
                 Matern(1|dec_lon_new + dec_lat_new),
               family=Gamma(log),
               data=data)
  message(paste("running", mod.nomes[3]), "...")
  m03 <- fitme(y ~ tempmax +
                 (1|species_name) + (1|Reference) +
                 Matern(1|dec_lon_new + dec_lat_new),
               family=Gamma(log),
               data=data)
  message(paste("running", mod.nomes[4]), "...")
  m04 <- fitme(y ~  Ocean +
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

