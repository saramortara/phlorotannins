#### Code for analysis of global phlorotannin concentration ####
## First version by Sara Mortara 11.04.2019

# loading packages
library(stringr)

# reading data
data <- read.csv("../data/florotaninos.csv", sep=";", as.is=TRUE, fill=TRUE)

head(data)
names(data)[c(1,4)] <- c("latitude", "TPC_perc_DW") 

# creating column with decimal latitude
data$dec_latitude <- as.numeric(sapply(str_split(data$latitude,"Â"), function(x) x[1]))
data$NS <- sapply(str_split(data$latitude,"Â"), function(x) x[2])

data$dec_latitude <- ifelse(str_detect(data$NS, "S"),  
                            data$dec_latitude*-1, data$dec_latitude*1)

head(data)
tail(data)

# creating column with species name only 
data$genus <- sapply(str_split(data$Species, "( )"), function(x) x[1])
data$epithet <- sapply(str_split(data$Species, "( )"), function(x) x[2])
data$species_name <- paste(data$genus, data$epithet, sep="_")
head(data)

# calculating richness per study site
studies <- aggregate(data$species_name, by=list(Reference=data$Reference, dec_latitude=data$dec_latitude), 
                      function(x) length(unique(x)))
names(studies)[3] <- "Richness"

head(studies)

studies$TPC <- aggregate(data$TPC_perc_DW, by=list(Reference=data$Reference, dec_latitude=data$dec_latitude), 
                         mean)$x

head(studies)

plot(TPC ~ Richness, data=studies, las=1)
plot(Richness ~ dec_latitude, data=studies, las=1)
plot(TPC ~ dec_latitude, data=studies, las=1, log='y')

