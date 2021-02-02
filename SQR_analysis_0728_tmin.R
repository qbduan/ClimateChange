#include the seanatility term in the model.

library("Rcpp", lib.loc = '/home/duanq/ClimateProj/package/')
library("RcppArmadillo", lib.loc = '/home/duanq/ClimateProj/package/')
sourceCpp("SpatialQR_cpp.cpp")
source("SpatialQR_QB_2.R")
runid = 917 #as.numeric(Sys.getenv("PBS_ARRAY_INDEX"));

data.tg_elen<-readRDS("mydata365_19.60_n.rds")
data.tg_elen<-data.tg_elen[which(!is.na(data.tg_elen$elevation)),]  #rm elevation is NA

list.na20<- c(96003, 23013, 8095,  39128, 9500,  9053,  86038, 14626, 22008, 40211, 85279, 47048, 33045) #list of station id that miss data is >20%
id <- which(!data.tg_elen$stationID%in%list.na20)
data.tg_elen <-data.tg_elen[id,]

#sds<-c(66:70)
list.stations<-unique(data.tg_elen[c("stationID","Lat","Lon","elevation")])#[sds,]

id <- which(data.tg_elen$stationID%in%list.stations$stationID)
data.tg_elen <-data.tg_elen[id,]

t.rnum <- nrow(data.tg_elen)
t.cnum <- nrow(list.stations)

library(dplyr)
data.byst<-data.tg_elen %>%
  group_split(stationID) 

data.tg_elen$Lats <- (data.tg_elen$Lat-min(data.tg_elen$Lat))/(max(data.tg_elen$Lat)-min(data.tg_elen$Lat))
data.tg_elen$Lons <- (data.tg_elen$Lon-min(data.tg_elen$Lon))/(max(data.tg_elen$Lon)-min(data.tg_elen$Lon))
data.tg_elen$elevations <- (data.tg_elen$elevation - min(data.tg_elen$elevation))/(max(data.tg_elen$elevation)-min(data.tg_elen$elevation))
data.tg_elen$msois <- (data.tg_elen$msoi - min(data.tg_elen$msoi))/(max(data.tg_elen$msoi)-min(data.tg_elen$msoi))

data.byst_elen_y <- data.tg_elen %>%group_split(stationID) 

library(fields)
library(mvtnorm)

ns<- t.cnum                    #Number of spatial locations
nt<- nrow(data.byst_elen_y[[1]])        #Number of observations

sc01 <- function(x){(x-min(x))/(max(x)-min(x))}
scaled_lat <- sc01(list.stations$Lat)
scaled_long <- sc01(list.stations$Lon)

s<- cbind(scaled_lat,scaled_long)                    #Spatial locations (1-D)

#tau.grid<-seq(0.01,0.99,0.01) #Quantile levels to keep track of
tau.grid<-seq(0.025,0.975,0.025)
ntau<-length(tau.grid)

sigma_st <- readRDS("sigma_tmin.rds")[1:72,]

sigma_st60 <- apply(sigma_st, 1, FUN = function(x){matrix(x,1,365*60,byrow = TRUE)}) #21900    73


FSk <- function(x, k = 4) {
  X <- matrix(0, length(x), 2 * k)
  for (j in 1:k) {
    #if unit is in days
    X[, 2*j-1] <- cos(2*pi*x * j/365)
    X[, 2*j  ] <- sin(2*pi*x * j/365)
    #if unit is year
    #X[, 2 * j - 1] <- cos(2 * pi * x * j)
    #X[, 2 * j] <- sin(2 * pi * x * j)
  }
  dimnames(X) <-
    list(NULL, as.vector(outer(c("c", "s"), 1:k, paste, sep = "")))
  X
  
}

FS4<- FSk(1:nt,k=4)
#  X[t,s,k] is the kth covariate at site number s for year t.
X<-array(1,c(nt,ns,3+8))              #X[,,1] is the intercept

year<-seq(0,1,length=nt)
for(t in 1:nt){
  X[t,,2]<-year[t]
  X[t,,3] <- FS4[t,1]
  X[t,,4] <- FS4[t,2]
  X[t,,5] <- FS4[t,3]
  X[t,,6] <- FS4[t,4]
  X[t,,7] <- FS4[t,5]
  X[t,,8] <- FS4[t,6]
  X[t,,9] <- FS4[t,7]
  X[t,,10] <- FS4[t,8]
  
} #X[,,2] is the year (t) ##X[,,3] is another covariate

for(j in 1:ns)
{
  X[,j,11]<- sigma_st60[,j] # varying sigma_t
  #X[,j,12]<-data.byst_elen_y[[j]]$msois 
 }

Y<- matrix(0,nt,ns)

for(i in 1:ns){Y[,i]<-data.byst_elen_y[[i]]$tmin} #with seasonality

rm(myslctdata,data.tg_elen,data.byst_elen)
iniallist = readRDS("intermlist.rds")
fitourcase1<-QR_Spatial_2(Y,s,X,L=4,tau=tau.grid,iters=1100,burn = 100, iniallist =iniallist )

filenames <- paste("fittingsqr_allsts_tmin_nosoi",runid,".rds",sep = "")
saveRDS(fitourcase1,file =filenames)




