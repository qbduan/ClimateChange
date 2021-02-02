
# per decade
rm(list = ls())
library(raster)
library(sp)
library(rgdal)
library(gstat)
library(raster)
library(maptools)
library(lattice)
library(maps)
library(mapdata)
library(oz)
library(ggplot2)
setwd("C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Rfiles_QB")

data.tg_elen<-readRDS("mydata365_19.60_n.rds")
data.tg_elen<-data.tg_elen[which(!is.na(data.tg_elen$elevation)),]  #rm elevation is NA

list.na20<- c(96003, 23013, 8095,  39128, 9500,  9053,  86038, 14626, 22008, 40211, 85279, 47048, 33045) #list of station id that miss data is >20%
id <- which(!data.tg_elen$stationID%in%list.na20)
data.tg_elen <-data.tg_elen[id,]

list.stations<-unique(data.tg_elen[c("stationID","Lat","Lon","elevation")])

######Get Australian map data##########

xy2df <- function(xy) {
  data.frame(x=xy$x,y=xy$y)
}
closeRing <- function(coords){
  
  if(coords$x[length(coords$x)]!=coords$x[1]){
    coords <- rbind(coords,coords[1,])
  }
  coords
}
oz2sp <- function(oz) {
  oz.reg <- ozRegion()
  ply <- list()
  #WA
  coords=rbind(xy2df(oz.reg$lines[[8]])[nrow(xy2df(oz.reg$lines[[8]])):1,],
               xy2df(oz.reg$lines[[1]]),xy2df(oz.reg$lines[[9]])[nrow(xy2df(oz.reg$lines[[9]])):1,])
  ply[[1]]<-Polygons(list(Polygon(coords)),ID="WA")
  #NT
  coords=rbind(xy2df(oz.reg$lines[[2]]),
               xy2df(oz.reg$lines[[11]])[nrow(xy2df(oz.reg$lines[[11]])):1,],
               xy2df(oz.reg$lines[[10]])[nrow(xy2df(oz.reg$lines[[10]])):1,],
               xy2df(oz.reg$lines[[9]]))
  ply[[2]]<-Polygons(list(Polygon(coords)),ID="NT")
  #QLD
  coords=rbind(xy2df(oz.reg$lines[[3]]),
               xy2df(oz.reg$lines[[13]]),
               xy2df(oz.reg$lines[[12]])[nrow(xy2df(oz.reg$lines[[12]])):1,],
               xy2df(oz.reg$lines[[11]]))
  ply[[3]]<-Polygons(list(Polygon(coords)),ID="QLD")
  #NSW
  coords=rbind(xy2df(oz.reg$lines[[4]]),
               xy2df(oz.reg$lines[[15]]),#[nrow(xy2df(oz.reg$lines[[15]])):1,],
               xy2df(oz.reg$lines[[14]])[nrow(xy2df(oz.reg$lines[[14]])):1,],
               xy2df(oz.reg$lines[[13]])[nrow(xy2df(oz.reg$lines[[13]])):1,])
  ply[[4]]<-Polygons(list(Polygon(coords)),ID="NSW")
  #VIC
  coords=rbind(xy2df(oz.reg$lines[[16]]),
               xy2df(oz.reg$lines[[5]])[nrow(xy2df(oz.reg$lines[[5]])):1,],
               xy2df(oz.reg$lines[[15]]))
  ply[[5]]<-Polygons(list(Polygon(coords)),ID="VIC")
  #TAS
  ply[[6]]<-Polygons(list(Polygon(closeRing(xy2df(oz.reg$lines[[6]])))),ID="TAS")
  #SA
  coords=rbind(xy2df(oz.reg$lines[[7]]),
               xy2df(oz.reg$lines[[8]]),
               xy2df(oz.reg$lines[[10]]),
               xy2df(oz.reg$lines[[12]]),
               xy2df(oz.reg$lines[[14]]),
               xy2df(oz.reg$lines[[16]]))
  ply[[7]]<-Polygons(list(Polygon(coords)),ID="SA")
  SpatialPolygons(ply)
}

coordoz.sp<-oz2sp(oz())


s1 <-list.stations$Lon
s2 <-list.stations$Lat
plys <- coordoz.sp@polygons


###################################################
##########plot section############################

maxmintext<-c("max","min")

###################################
#########Daily maxtimum############
resall_tmax<-readRDS("fittingsqr_allsts_tmax904.rds")
mm <- 1

resultall <- resall_tmax
mapqtf1 <- resultall$qfx.mn[,2,]/6 # per decade
AST_all<-rowMeans(mapqtf1)*10

#mapqtvar1 <- resultall$qfx.var[,2,]
coordinates(list.stations) <- ~Lon + Lat
insidepts <- over(list.stations,coordoz.sp)
region<- names(coordoz.sp)


difmxmn<-apply(mapqtf1,2,max) - apply(mapqtf1,2,min)

myqtf1<-as.data.frame(cbind(s1,s2,t(mapqtf1)))

colnames(myqtf1)[3:41]<-paste0("tau",seq(0.025,0.975,by=0.025),seq="")

dsp <- SpatialPoints(myqtf1[,1:2], proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
dsp <- SpatialPointsDataFrame(dsp,myqtf1)


blues <-col <- c( "#2B00FFFF", "#9F0FF0FF", "#F647B8FF","#FF6C93FF","#FF916EFF","#FFB649FF","#FFDB24FF","#FFFF09FF")#,"#FFFF60FF")

pols <- list("sp.polygons",coordoz.sp, fill = "transparent")#

cuts <- c(-0.4,-0.2,0,0.1,0.2,0.3,0.4,0.5,0.6)

filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/figure2_dotsplot_09_tmax_de.pdf")
pdf(file = filedir,width = 5.6,height = 6.0) # The height of the plot in inches

#spplot(dsp,c("tau0.9"),sp.layout=pols,xlab="",ylab= "",xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
#       scales=list(draw = TRUE),key.space=list(x = 1, y = 0.9, corner = c(0, 1)),col.regions = blues, main= "Quantile function of 0.9 quantile",cuts= cuts, cex=1.5)
spplot(dsp,c("tau0.9"),sp.layout=pols,xlab=NULL,ylab= NULL,xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = FALSE),key.space=list(x = 0.31, y = 0.31, corner = c(0, 1)),col.regions = blues, 
       main= "Trend function of 0.9 quantile (Dmx)",cuts= cuts, cex=1.5,
       par.settings = list(axis.line = list(col = 0)))

dev.off()

filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/figure2_dotsplot_05_tmax_de.pdf")
pdf(file = filedir,width = 5.6,height = 6.0) # The height of the plot in inches

spplot(dsp,c("tau0.5"),sp.layout=pols,xlab=NULL,ylab= NULL,xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = FALSE),key.space=list(x = 0.31, y = 0.31, corner = c(0, 1)),col.regions = blues, 
       main= "Trend function of 0.5 quantile (Dmx)",cuts= cuts, cex=1.5,
       par.settings = list(axis.line = list(col = 0)))

dev.off()


filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/figure2_dotsplot_01_tmax_de.pdf")
pdf(file = filedir,width = 5.6,height = 6.0) # The height of the plot in inches

spplot(dsp,c("tau0.1"),sp.layout=pols,xlab=NULL,ylab= NULL,xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = FALSE),key.space=list(x = 0.31, y = 0.31, corner = c(0, 1)),col.regions = blues, 
       main= "Trend function of 0.1 quantile (Dmx)",cuts= cuts, cex=1.5,
       par.settings = list(axis.line = list(col = 0)))

dev.off()


trellis.par.set(sp.theme())
trellis.pars.save <- trellis.par.get() # for good measure
trellis.pars <- trellis.par.get()
layout.height <- trellis.par.get("layout.heights")
layout.widths <- trellis.par.get("layout.widths")
layout.axis <- trellis.par.get("axis.line")
layout.axis$alpha <- 0
#names(layout.height)
#layout.height$top.padding
layout.height$top.padding <-
  0.1 ### here is where margins etc. can be changed
layout.height$main.key.padding <- 0
layout.height$bottom.padding <- 0
layout.height$between <- 0.1
layout.widths$left.padding <- 0
layout.widths$right.padding <- 0
trellis.par.set(layout.heights = layout.height, layout.widths = layout.widths,layout.axis =layout.axis)

labels <- list("sp.text",coordinates(list.stations),seq(1,73,1),1)
spplot(dsp1,c("tau0.1"),sp.layout=list(pols,labels),xlab="",ylab= "",xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = TRUE),col.regions = blues, main= NA,cuts=cuts,auto.key = FALSE,cex=1.2)

WGS84<- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")#
dsp1<-spTransform(dsp,WGS84)
proj4string(coordoz.sp)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
bound1<-spTransform(coordoz.sp,WGS84)


library(dismo)
library(rgeos)
v <- voronoi(dsp,ext = c(113.1939, 153.6692, -43.59316, -10.93156))
plot(v)
centroids <- coordinates(coordoz.sp)
x <- centroids[,1]
y <- centroids[,2]
au <- SpatialPolygonsDataFrame(coordoz.sp,
                                   data=data.frame(x=x, y=y, row.names=row.names(coordoz.sp)))
auta <- spTransform(au, WGS84)
auta <- gSimplify(auta, tol = 0.00001)
#auta <- gBuffer(auta, byid=TRUE, width=0)

au <- aggregate(auta)

vca <- intersect(au,v)

##################
col <- c( "#00EAFF", "#00FFFF", "#FFEA00",
         "#FFD500", "#FFBF00", "#FFAA00", "#FF9500", "#FF8000", "#FF6A00", "#FF5500",
         "#FF4000", "#FF2B00", "#FF1500", "#FF0000")

spplot(vca, 'tau0.1', col.regions=col,xlim=c(113,154.1),ylim=c(-43.5,-10.5),colorkey = list(
  at = seq(from = -0.5, to = 3, by = 0.25),
  labels = list(at = seq(0, 3, by = 1),
                labels = as.character(seq(
                  from = 0, to = 3, by = 1
                )))
))


#layout(matrix(1:4, ncol=2, byrow=TRUE))

r<-500
blank_raster<-raster(nrow=r,ncol=r,extent(coordoz.sp))
values(blank_raster)<-0
bound_raster<-rasterize(coordoz.sp,blank_raster)
bound_raster[!(is.na(bound_raster))] <- 0
#plot(bound_raster,main=paste("Res: ",r,"*",r))
#plot(coordoz.sp,add=T)
crs(blank_raster) <-WGS84

vr <- rasterize(vca,bound_raster,"tau0.1")
plot(vr)

#0.1 5.58X4.76
gs01<-gstat(formula=tau0.1~1,location=dsp1,nmax=5,set=list(idp=0))
crs(gs01) <- WGS84
nn<-interpolate(blank_raster,gs01) #not the resized one
nnmask<-mask(nn,vr)##????????????
plot(nnmask)

gs01 <- gstat(formula=tau0.1~1, locations=dsp1)
idw01 <- interpolate(blank_raster, gs01)
idwmask01<-mask(idw01,vr)
plot(idwmask01)

gs05<-gstat(formula=tau0.5~1,location=dsp1,nmax=5,set=list(idp=0))
crs(gs05) <- WGS84
nn<-interpolate(blank_raster,gs05) #not the resized one
nnmask<-mask(nn,vr)##????????????
plot(nnmask)

gs05 <- gstat(formula=tau0.5~1, locations=dsp1)
idw05 <- interpolate(blank_raster, gs05)
idwmask05<-mask(idw05,vr)
plot(idwmask05)

#spplot(idwmask)
gs09<-gstat(formula=tau0.9~1,location=dsp1,nmax=5,set=list(idp=0))
crs(gs09) <- WGS84
nn<-interpolate(blank_raster,gs09) #not the resized one
nnmask<-mask(nn,vr)##????????????
plot(nnmask)

gs09 <- gstat(formula=tau0.9~1, locations=dsp1)
idw09 <- interpolate(blank_raster, gs09)
idwmask09<-mask(idw09,vr)
plot(idwmask09)

colss<-c(get_col_regions()[seq(1,100,4)])

pts <- list('sp.points',
            dsp,
            col = "black",
            pch = 16,
            cex = 0.6)

######### 0.1 #############

filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/Figure_map_all_01_tmax_de.pdf")
pdf(file = filedir,width = 5.6,height = 5.0) # The height of the plot in inches
spplot(
  idwmask01,
  sp.layout = list(pts),
  scales = list(draw = TRUE),
  col.regions = colss[14:25],
  cuts = 12,
  colorkey = list(
    at = seq(from = 0, to = 0.6, by = 0.05)
  ),
  main = expression(paste(
    "Trend function of Dmx, ", tau, "=", "0.1"
  ))
)
dev.off()

#was [14:25] when 0-0.05 is #D934CBFF

######### 0.5 #############

filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/Figure_map_all_05_tmax_de.pdf")
pdf(file = filedir,width = 5.6,height = 5.0) # The height of the plot in inches
spplot(
  idwmask05,
  sp.layout = list(pts),
  scales = list(draw = TRUE),
  cuts = 15,
  col.regions = colss[9:23],
  colorkey = list(
    at = seq(from = -0.25, to = 0.5, by = 0.05)
  ),
  main = expression(paste(
    "Trend function of Dmx, ", tau, "=", "0.5"
  ))
)
dev.off()


######### 0.9 #############

filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/Figure_map_all_09_tmax_de.pdf")
pdf(file = filedir,width = 5.6,height = 5.0) # The height of the plot in inches
spplot(
  idwmask09,
  sp.layout = list(pts),
  scales = list(draw = TRUE),
  col.regions = colss[7:22],
  cuts = 16,
  colorkey = list(
    at = seq(from = -0.35, to = 0.45, by = 0.05)
  ),
  main = expression(paste(
    "Trend function of Dmx, ", tau, "=", "0.9"
  ))
)

dev.off()


###################################
#########Daily minimum ############
s1 <-list.stations$Lon
s2 <-list.stations$Lat
plys <- coordoz.sp@polygons
resall_tmin<-readRDS("fittingsqr_allsts_tmin904.rds")

resultall <- resall_tmin
mapqtf1 <- resultall$qfx.mn[,2,]/6

mapqtvar1 <- resultall$qfx.var[,2,]

insidepts <- over(list.stations,coordoz.sp)
region<- names(coordoz.sp)
#stid <- which(insidepts==3)
#mapqtf1_stid <- mapqtf1[,which(insidepts==3)]

difmxmn<-apply(mapqtf1,2,max) - apply(mapqtf1,2,min)

myqtf1<-as.data.frame(cbind(s1,s2,t(mapqtf1)))

colnames(myqtf1)[3:41]<-paste0("tau",seq(0.025,0.975,by=0.025),seq="")

dsp <- SpatialPoints(myqtf1[,1:2], proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
dsp <- SpatialPointsDataFrame(dsp,myqtf1)

max(mapqtf1[c(4,20,36),])

blues <-col <- c( "#2B00FFFF", "#9F0FF0FF", "#F647B8FF","#FF6C93FF","#FF916EFF","#FFB649FF","#FFDB24FF","#FFFF09FF","#FFFF60FF")
pols <- list("sp.polygons",coordoz.sp, fill = "transparent")#
cuts <- c(-0.4,-0.2,0,0.1,0.2,0.3,0.4,0.5,0.6,0.70)


spplot(dsp,c("tau0.1"),sp.layout=pols,xlab=NULL,ylab= NULL,xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = FALSE),key.space=list(x = 0.31, y = 0.31, corner = c(0, 1)),col.regions = blues, 
       main= "Trend function of 0.1 quantile",cuts= cuts, cex=1.5,
       par.settings = list(axis.line = list(col = 0)))

filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/figure2_dotsplot_09_tmin_de.pdf")
pdf(file = filedir,width = 5.6,height = 6.0) # The height of the plot in inches

spplot(dsp,c("tau0.9"),sp.layout=pols,xlab=NULL,ylab= NULL,xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = FALSE),key.space=list(x = 0.31, y = 0.31, corner = c(0, 1)),col.regions = blues, 
       main= "Trend function of 0.9 quantile (Dmn)",cuts= cuts, cex=1.5,
       par.settings = list(axis.line = list(col = 0)))

dev.off()

filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/figure2_dotsplot_05_tmin_de.pdf")
pdf(file = filedir,width = 5.6,height = 6.0) # The height of the plot in inches

spplot(dsp,c("tau0.5"),sp.layout=pols,xlab=NULL,ylab= NULL,xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = FALSE),key.space=list(x = 0.31, y = 0.31, corner = c(0, 1)),col.regions = blues, 
       main= "Trend function of 0.5 quantile (Dmn)",cuts= cuts, cex=1.5,
       par.settings = list(axis.line = list(col = 0)))

dev.off()

filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/figure2_dotsplot_01_tmin_de.pdf")
pdf(file = filedir,width = 5.6,height = 6.0) # The height of the plot in inches

spplot(dsp,c("tau0.1"),sp.layout=pols,xlab=NULL,ylab= NULL,xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = FALSE),key.space=list(x = 0.31, y = 0.31, corner = c(0, 1)),col.regions = blues, 
       main= "Trend function of 0.1 quantile (Dmn)",cuts= cuts, cex=1.5,
       par.settings = list(axis.line = list(col = 0)))

dev.off()


labels <- list("sp.text",coordinates(list.stations),seq(1,72,1),1)
spplot(dsp1,c("tau0.1"),sp.layout=list(pols,labels),xlab="",ylab= "",xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = TRUE),col.regions = blues, main= NA,cuts=cuts,auto.key = FALSE,cex=1.2)

WGS84<- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")#
dsp1<-spTransform(dsp,WGS84)#
proj4string(coordoz.sp)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
bound1<-spTransform(coordoz.sp,WGS84)


r<-500
blank_raster<-raster(nrow=r,ncol=r,extent(coordoz.sp))
values(blank_raster)<-0
bound_raster<-rasterize(coordoz.sp,blank_raster)
bound_raster[!(is.na(bound_raster))] <- 0
#plot(bound_raster,main=paste("Res: ",r,"*",r))
#plot(coordoz.sp,add=T)
crs(blank_raster) <-WGS84

vr1 <- rasterize(vca,bound_raster,"tau0.1")
plot(vr1)

#0.1 5.58X4.76
gs01<-gstat(formula=tau0.1~1,location=dsp1,nmax=5,set=list(idp=0))
crs(gs01) <- WGS84
nn<-interpolate(blank_raster,gs01) #not the resized one
nnmask<-mask(nn,vr1)##????????????
plot(nnmask)

gs01 <- gstat(formula=tau0.1~1, locations=dsp1)
idw01 <- interpolate(blank_raster, gs01)
idwmask01<-mask(idw01,vr1)
plot(idwmask01)

#spplot(idwmask)

#0.5
vr5 <- rasterize(vca,bound_raster,"tau0.5")
plot(vr5)
gs05<-gstat(formula=tau0.5~1,location=dsp1,nmax=5,set=list(idp=0))
crs(gs05) <- WGS84
nn<-interpolate(blank_raster,gs05) #not the resized one
nnmask<-mask(nn,vr5)##????????????
plot(nnmask)

gs05 <- gstat(formula=tau0.5~1, locations=dsp1)
idw05 <- interpolate(blank_raster, gs05)
idwmask05<-mask(idw05,vr5)
plot(idwmask05)

#spplot(idwmask)

#0.9

vr9 <- rasterize(vca,bound_raster,"tau0.9")
plot(vr9)

gs09<-gstat(formula=tau0.9~1,location=dsp1,nmax=5,set=list(idp=0))
crs(gs09) <- WGS84
nn<-interpolate(blank_raster,gs09) #not the resized one
nnmask<-mask(nn,vr9)##????????????
plot(nnmask)

gs09 <- gstat(formula=tau0.9~1, locations=dsp1)
idw09 <- interpolate(blank_raster, gs09)
idwmask09<-mask(idw09,vr9)
plot(idwmask09)

#spplot(idwmask)

colss<-c(get_col_regions()[seq(1,100,4)],"#FFFF26FF","#FFFF98FF")
#cols<-colorRampPalette(rev(c('red','yellow', 'blue')))
bds <-
  list(
    "sp.polygons",
    coordoz.sp,
    fill = "transparent",
    col = "black",
    border = "green"
  )#
pts <- list('sp.points',
            dsp,
            col = "black",
            pch = 16,
            cex = 0.6)

filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/Figure_map_all_01_tmin_de.pdf")
pdf(file = filedir,width = 5.6,height = 5.0) # The height of the plot in inches
spplot(
  idwmask01,
  sp.layout = list(pts),
  scales = list(draw = TRUE),
  cuts = 20,
  col.regions = colss[7:26],
  colorkey = list(
    at = seq(from = -0.35, to = 0.65, by = 0.05)
  ),
  main = expression(paste(
    "Trend function of Dmn, ", tau, "=", "0.1"
  ))
)
dev.off()

filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/Figure_map_all_05_tmin_de.pdf")
pdf(file = filedir,width = 5.6,height = 5.0) # The height of the plot in inches
spplot(
  idwmask05,
  sp.layout = list(pts),
  scales = list(draw = TRUE),
  cuts = 15,
  col.regions = colss[9:23],
  colorkey = list(
    at = seq(from = -0.25, to = 0.5, by = 0.05)
  ),
  main = expression(paste(
    "Trend function of Dmn, ", tau, "=", "0.5"
  ))
)
dev.off()


filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/Figure_map_all_09_tmin_de.pdf")
pdf(file = filedir,width = 5.6,height = 5.0) # The height of the plot in inches
spplot(
  idwmask09,
  sp.layout = list(pts),
  scales = list(draw = TRUE),
  col.regions = colss[11:21],
  cuts= 10,
  colorkey = list(
    at = seq(from = -0.15, to = 0.40, by = 0.05)
  ),
  main = expression(paste(
    "Trend function of Dmn, ", tau, "=", "0.9"
  ))
)
dev.off()
col.regions = colss[11:21]
at =   seq(from = -0.15, to = 0.4, by = 0.05)
cbind(col.regions,at)

spplot(
  idwmask09,
  sp.layout = list(pts),
  scales = list(draw = TRUE),
  cuts= 10,
  col.regions = colss[11:21],
  main = expression(paste(
    "Trend function of Dmn, ", tau, "=", "0.9"
  ))
)

#################################################
##### plot by two seasons #######################

coordinates(list.stations) <- ~Lon + Lat
insidepts <- over(list.stations,coordoz.sp)
region<- names(coordoz.sp)

blues <-col <- c( "#2B00FFFF", "#9F0FF0FF", "#F647B8FF","#FF6C93FF","#FF916EFF","#FFB649FF","#FFDB24FF","#FFFF09FF","#FFFF60FF")
pols <- list("sp.polygons",coordoz.sp, fill = "transparent")#
cuts <- c(-0.4,-0.2,0,0.1,0.2,0.3,0.4,0.5,0.6,0.70)

ss <-1
myses <- c("Summer","Winter")

###########Summer maximum #######################
resall<-readRDS("fittingsqr_allsts_tmax_ss1.rds") 

resultall <- resall
mapqtf1 <- resultall$qfx.mn[,2,]/6
mapqtvar1 <- resultall$qfx.var[,2,]

difmxmn<-apply(mapqtf1,2,max) - apply(mapqtf1,2,min)
myqtf1<-as.data.frame(cbind(s1,s2,t(mapqtf1)))
colnames(myqtf1)[3:41]<-paste0("tau",seq(0.025,0.975,by=0.025),seq="")
dsp <- SpatialPoints(myqtf1[,1:2], proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
dsp <- SpatialPointsDataFrame(dsp,myqtf1)


filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/figure2_dotsplot_09_tmax_",myses[ss],"_de.pdf")
pdf(file = filedir,width = 5.6,height = 6.0) # The height of the plot in inches

spplot(dsp,c("tau0.9"),sp.layout=pols,xlab=NULL,ylab= NULL,xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = FALSE),key.space=list(x = 0.31, y = 0.31, corner = c(0, 1)),col.regions = blues, 
       main= "Trend function of 0.9 quantile (summer, Dmx)",cuts= cuts, cex=1.5,
       par.settings = list(axis.line = list(col = 0)))

dev.off()

filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/figure2_dotsplot_05_tmax_",myses[ss],"_de.pdf")
pdf(file = filedir,width = 5.6,height = 6.0) # The height of the plot in inches

spplot(dsp,c("tau0.5"),sp.layout=pols,xlab=NULL,ylab= NULL,xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = FALSE),key.space=list(x = 0.31, y = 0.31, corner = c(0, 1)),col.regions = blues, 
       main= "Trend function of 0.5 quantile (summer, Dmx)",cuts= cuts, cex=1.5,
       par.settings = list(axis.line = list(col = 0)))

dev.off()


filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/figure2_dotsplot_01_tmax",myses[ss],"_de.pdf")
pdf(file = filedir,width = 5.6,height = 6.0) # The height of the plot in inches

spplot(dsp,c("tau0.1"),sp.layout=pols,xlab=NULL,ylab= NULL,xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = FALSE),key.space=list(x = 0.31, y = 0.31, corner = c(0, 1)),col.regions = blues, 
       main= "Trend function of 0.1 quantile (summer, Dmx)",cuts= cuts, cex=1.5,
       par.settings = list(axis.line = list(col = 0)))

dev.off()

###########Winter maximum #######################
ss <-2
myses <- c("Summer","Winter")

resall<-readRDS("fittingsqr_allsts_tmax_ss3.rds")
resultall <- resall
mapqtf1 <- resultall$qfx.mn[,2,]/6
mapqtvar1 <- resultall$qfx.var[,2,]
myqtf1<-as.data.frame(cbind(s1,s2,t(mapqtf1)))
colnames(myqtf1)[3:41]<-paste0("tau",seq(0.025,0.975,by=0.025),seq="")
#myqtf1[31,"tau0.9"] <- -2
dsp <- SpatialPoints(myqtf1[,1:2], proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
dsp <- SpatialPointsDataFrame(dsp,myqtf1)



filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/figure2_dotsplot_09_tmax_",myses[ss],"_de.pdf")
pdf(file = filedir,width = 5.6,height = 6.0) # The height of the plot in inches

spplot(dsp,c("tau0.9"),sp.layout=pols,xlab=NULL,ylab= NULL,xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = FALSE),key.space=list(x = 0.31, y = 0.31, corner = c(0, 1)),col.regions = blues, 
       main= "Trend function of 0.9 quantile (winter, Dmx)",cuts= cuts, cex=1.5,
       par.settings = list(axis.line = list(col = 0)))

dev.off()

filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/figure2_dotsplot_05_tmax_",myses[ss],"_de.pdf")
pdf(file = filedir,width = 5.6,height = 6.0) # The height of the plot in inches

spplot(dsp,c("tau0.5"),sp.layout=pols,xlab=NULL,ylab= NULL,xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = FALSE),key.space=list(x = 0.31, y = 0.31, corner = c(0, 1)),col.regions = blues, 
       main= "Trend function of 0.5 quantile (winter, Dmx)",cuts= cuts, cex=1.5,
       par.settings = list(axis.line = list(col = 0)))

dev.off()


filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/figure2_dotsplot_01_tmax",myses[ss],"_de.pdf")
pdf(file = filedir,width = 5.6,height = 6.0) # The height of the plot in inches

spplot(dsp,c("tau0.1"),sp.layout=pols,xlab=NULL,ylab= NULL,xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = FALSE),key.space=list(x = 0.31, y = 0.31, corner = c(0, 1)),col.regions = blues, 
       main= "Trend function of 0.1 quantile (winter, Dmx)",cuts= cuts, cex=1.5,
       par.settings = list(axis.line = list(col = 0)))

dev.off()

###########Summer minimum #######################
ss <-1
myses <- c("Summer","Winter")

resall<-readRDS("fittingsqr_allsts_tmin_ss1.rds")
#resall<-readRDS("fittingsqr_allsts_tmin812.rds")
resultall <- resall
mapqtf1 <- resultall$qfx.mn[,2,]/6
mapqtvar1 <- resultall$qfx.var[,2,]

difmxmn<-apply(mapqtf1,2,max) - apply(mapqtf1,2,min)


myqtf1<-as.data.frame(cbind(s1,s2,t(mapqtf1)))

colnames(myqtf1)[3:41]<-paste0("tau",seq(0.025,0.975,by=0.025),seq="")

dsp <- SpatialPoints(myqtf1[,1:2], proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
dsp <- SpatialPointsDataFrame(dsp,myqtf1)


filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/figure2_dotsplot_09_tmin_",myses[ss],"_de.pdf")
pdf(file = filedir,width = 5.6,height = 6.0) # The height of the plot in inches

spplot(dsp,c("tau0.9"),sp.layout=pols,xlab=NULL,ylab= NULL,xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = FALSE),key.space=list(x = 0.31, y = 0.31, corner = c(0, 1)),col.regions = blues, 
       main= "Trend function of 0.9 quantile (summer, Dmn)",cuts= cuts, cex=1.5,
       par.settings = list(axis.line = list(col = 0)))

dev.off()

filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/figure2_dotsplot_05_tmin_",myses[ss],"_de.pdf")
pdf(file = filedir,width = 5.6,height = 6.0) # The height of the plot in inches

spplot(dsp,c("tau0.5"),sp.layout=pols,xlab=NULL,ylab= NULL,xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = FALSE),key.space=list(x = 0.31, y = 0.31, corner = c(0, 1)),col.regions = blues, 
       main= "Trend function of 0.5 quantile (summer, Dmn)",cuts= cuts, cex=1.5,
       par.settings = list(axis.line = list(col = 0)))

dev.off()


filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/figure2_dotsplot_01_tmin",myses[ss],"_de.pdf")
pdf(file = filedir,width = 5.6,height = 6.0) # The height of the plot in inches

spplot(dsp,c("tau0.1"),sp.layout=pols,xlab=NULL,ylab= NULL,xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = FALSE),key.space=list(x = 0.31, y = 0.31, corner = c(0, 1)),col.regions = blues, 
       main= "Trend function of 0.1 quantile (summer, Dmn)",cuts= cuts, cex=1.5,
       par.settings = list(axis.line = list(col = 0)))

dev.off()

###########Winter minimum #######################

ss <-2
myses <- c("Summer","Winter")

resall<-readRDS("fittingsqr_allsts_tmin_ss3.rds")
resultall <- resall
mapqtf1 <- resultall$qfx.mn[,2,]/6
mapqtvar1 <- resultall$qfx.var[,2,]

myqtf1<-as.data.frame(cbind(s1,s2,t(mapqtf1)))

colnames(myqtf1)[3:41]<-paste0("tau",seq(0.025,0.975,by=0.025),seq="")
#myqtf1[55,"tau0.1"] <- 4 # was 4.000496


dsp <- SpatialPoints(myqtf1[,1:2], proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
dsp <- SpatialPointsDataFrame(dsp,myqtf1)


filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/figure2_dotsplot_09_tmin_",myses[ss],"_de.pdf")
pdf(file = filedir,width = 5.6,height = 6.0) # The height of the plot in inches

spplot(dsp,c("tau0.9"),sp.layout=pols,xlab=NULL,ylab= NULL,xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = FALSE),key.space=list(x = 0.31, y = 0.31, corner = c(0, 1)),col.regions = blues, 
       main= "Trend function of 0.9 quantile (winter, Dmn)",cuts= cuts, cex=1.5,
       par.settings = list(axis.line = list(col = 0)))
dev.off()

filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/figure2_dotsplot_05_tmin_",myses[ss],"_de.pdf")
pdf(file = filedir,width = 5.6,height = 6.0) # The height of the plot in inches

spplot(dsp,c("tau0.5"),sp.layout=pols,xlab=NULL,ylab= NULL,xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = FALSE),key.space=list(x = 0.31, y = 0.31, corner = c(0, 1)),col.regions = blues, 
       main= "Trend function of 0.5 quantile (winter, Dmn)",cuts= cuts, cex=1.5,
       par.settings = list(axis.line = list(col = 0)))
dev.off()


filedir <-
  paste0(
    "C:/Users/duanq/Queensland University of Technology/You-Gan Wang - Climate(Glen)/Latex/Figures/figure2_dotsplot_01_tmin",myses[ss],"_de.pdf")
pdf(file = filedir,width = 5.6,height = 6.0) # The height of the plot in inches

spplot(dsp,c("tau0.1"),sp.layout=pols,xlab=NULL,ylab= NULL,xlim=c(113,154.1),ylim=c(-43.5,-10.5), 
       scales=list(draw = FALSE),key.space=list(x = 0.31, y = 0.31, corner = c(0, 1)),col.regions = blues, 
       main= "Trend function of 0.1 quantile (winter, Dmn)",cuts= cuts, cex=1.5,
       par.settings = list(axis.line = list(col = 0)))
dev.off()



#####################################
