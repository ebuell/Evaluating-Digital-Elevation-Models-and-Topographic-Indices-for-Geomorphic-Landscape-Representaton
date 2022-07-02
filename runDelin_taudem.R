#Run R delineation for !!!!!!!!!!!one specific DEM and DEM process!!!!!!!!!!!!!!
#Save the outputs from the delineation (slope, TIC, TIV, and SCA)
#to google drive
#Next step for regression analysis (aka not what we are writing
#a paper about - extract spatial characteristics for the soil 
#sample locations
#######please note that the user will have to change the structure of the readin and readout files ######
#######please note that this code relies heavily on correct installation of TauDEM (https://hydrology.usu.edu/taudem/taudem5/) - install before running code

if (!require("pacman")) install.packages("pacman")
pacman::p_load(raster,rgdal,googlesheets4,shapefiles)
library(beepr);library('RColorBrewer');library(corrplot)

#directory
basedir="G:/My Drive/Which DEM Paper/DEMFiles"
diroptions = c("1arcsec","1_3arcsec","LIDAR/2010","LIDAR/2018")
fileoptions = c("/1arcsec_arcproj_bilin.tif","/1_3arcsec_arcproj_bilin.tif","/2010LIDAR_arcproj_bilin.tif","/2018LIDAR_arcproj_bilin.tif")

#decide which DEM
demorigin = readline('1 for usgs, 2 for LIDAR ')
if(demorigin==1){
  aschoice = readline('1 for 1arcsec, 2 for 1/3arcsec ')
  if(aschoice==1){DEMchoice = 1;name = "1arcsec"}
  if(aschoice==2){DEMchoice = 2;name = "1_3arcsec"}
}
if(demorigin==2){
  aschoice = readline('1 for 2010, 2 for 2018 ')
  if(aschoice==1){DEMchoice = 3;name = "2010LIDAR"}
  if(aschoice==2){DEMchoice = 4;name = "2018LIDAR"}
}

d8ordinf = readline('1 for D8, 2 for Dinf ')
if(d8ordinf==1){d8ordinfstr = "d8"}
if(d8ordinf==2){d8ordinfstr = "dinf"}

#Get and process DEM desired
folder = paste(diroptions[DEMchoice],sep = "")
url=paste0(basedir,"/",folder,fileoptions[DEMchoice])
dem_utm=raster(url)


#Alter LIDAR DEM
if(DEMchoice==3||DEMchoice==4){
  #subtracting out the difference from the overpass on 460 and culvert at outlet
  url="G:/My Drive/Which DEM Paper/SWATHelperFiles/pathsunderoverpass/pathsunderoverpass.shp"
  overpass = readOGR(url)
  r = raster(ext = extent(dem_utm),resolution = res(dem_utm))
  overpassrast = rasterize(overpass,r,field = as.numeric(overpass$Id),background = 0)
  dem_utm = dem_utm - overpassrast
  #z coordinates to meters
  dem_utm = calc(dem_utm,fun = function(x)x/3.281)
}
###End DEM manipulations


#Import outlet point
url="G:/My Drive/Which DEM Paper/SWATHelperFiles/Monitoring_Point.shp"
outlet = readOGR(url)
crs_project = crs(outlet) #check that we are in the correct coordinate system
points(outlet,col = 'blue')

#pit remove
setwd("C:/Users/ebuell/Documents/RDump_TauDEM")
writeRaster(dem_utm,"logan",format = "GTiff",overwrite = TRUE)
system("mpiexec -n 8 pitremove -z logan.tif -fel loganfel.tif")
dem_utm_pitrm=raster("loganfel.tif")


#slope and contributing area calculation
if(d8ordinf=="1"){
  #D8 flow directions
  system("mpiexec -n 8 D8Flowdir -p loganp.tif -sd8 logans.tif -fel loganfel.tif",show.output.on.console=F,invisible=F)
  p=raster("loganp.tif")
  slp=raster("logans.tif")
  # Contributing area
  system("mpiexec -n 8 AreaD8 -p loganp.tif -ad8 logana.tif")
  ca=raster("logana.tif")
  sca = ca*res(dem_utm)[1]
}
if(d8ordinf=="2"){
  # DInf flow directions
  system("mpiexec -n 8 DinfFlowdir -ang loganang.tif -slp loganslp.tif -fel loganfel.tif",show.output.on.console=F,invisible=F)
  ang=raster("loganang.tif")
  slp=raster("loganslp.tif")
  # Contributing area
  system("mpiexec -n 8 AreaDinf -ang loganang.tif -sca logansca.tif")
  sca=raster("logansca.tif")
}

#define threshhold for picking outlet point
threshhold = ceiling(max(values(sca),na.rm = TRUE)/20/res(dem_utm)[1])
system(paste0("mpiexec -n 8 Threshold -ssa logana.tif -src logansrc.tif -thresh ",threshhold))
src=raster("logansrc.tif")
plot(src)

#write outlet
shapefile(outlet, "approxoutlet.shp",overwrite = TRUE)

# Move Outlet
system("mpiexec -n 8 moveoutletstostreams -p loganp.tif -src logansrc.tif -o approxoutlet.shp -om Outlet.shp")
outpt=read.shp("Outlet.shp")

# Contributing area upstream of outlet
system("mpiexec -n 8 Aread8 -p loganp.tif -o Outlet.shp -ad8 loganssa.tif")
ssa=raster("loganssa.tif")
plot(ssa) 

#clip the catchment area and slp to the watershed
sca = mask(sca,ssa)
slp = mask(slp,ssa)
riverrast = mask(src,ssa); values(riverrast)[which(values(riverrast)==0)] = NA
if(DEMchoice==3|DEMchoice==4){
  if(DEMchoice==3){thinby = 3}
  if(DEMchoice==4){thinby = 5}
  ind = which(values(riverrast)==1)
  for(i in 1:length(ind)){
    if(ind[i]%%thinby==0){ind[i] = NA}
  }
  ind = ind[!is.na(ind)]
  values(riverrast)[ind] = NA
}
river = SpatialPointsDataFrame(rasterToPoints(riverrast),data.frame(hold = matrix(1,length(rasterToPoints(riverrast)[,1]),1)))
plot(river,pch = 15)
if(d8ordinf==2){p = ang}
aspect = mask(p,ssa)


#calculated tiv
TI = log((sca+1)/(slp+0.00001))
plot(TI)

#split into TICs
pacman::p_load(classInt)
nTIclass=10 #number of TI classes, currently equal area, can adjust method various ways e.g., classIntervals(v, n = nTIclass, style = "jenks")
v=values(TI); v=v[!is.na(v)]
brks.qt = classIntervals(v, n = nTIclass, style = "quantile")$brks #length nTIclass+1 of just the numeric breakpoints
TIC = cut(TI, breaks=brks.qt, include.lowest = T, right=T)
plot(TIC)#,col = rainbow(10))

if(d8ordinf=="1"){d8ordinfstr2 = "_D8"}
if(d8ordinf=="2"){d8ordinfstr2 = "_Dinf"}

name2 = name
if(DEMchoice==3){name2 = "LIDAR_2010"}
if(DEMchoice==4){name2 = "LIDAR_2018"}

#### UNITS of things!!! #####
#Slope: length/length
#SCA: meters
#TIC: unitless
#TI: something complicated that relates to slope and SCA units
#Note about SCA: SCA_D8 = catchment area(m2)/length of pixel(m); SCA_Dinf calculated as seen in Tarboton 1997

setwd(paste0("G:/My Drive/Which DEM Paper/TiffFiles/TIClass_",name2,"/TIClass_",name2,d8ordinfstr2))
writeRaster(TIC,paste0("TIClass_",name,d8ordinfstr2),format = "GTiff",overwrite = TRUE)
setwd(paste0("G:/My Drive/Which DEM Paper/TiffFiles/TIValue_",name2,"/TIValue_",name2,d8ordinfstr2))
writeRaster(TI,paste0("TIValue_",name,d8ordinfstr2),format = "GTiff",overwrite = TRUE)
setwd(paste0("G:/My Drive/Which DEM Paper/TiffFiles/Slope_",name2,"/Slope_",name2,d8ordinfstr2))
writeRaster(slp,paste0("Slope_",name,d8ordinfstr2),format = "GTiff",overwrite = TRUE)
setwd(paste0("G:/My Drive/Which DEM Paper/TiffFiles/SCA_",name2,"/SCA_",name2,d8ordinfstr2))
writeRaster(sca,paste0("SCA_",name,d8ordinfstr2),format = "GTiff",overwrite = TRUE)
setwd(paste0("G:/My Drive/Which DEM Paper/TiffFiles/Aspect_",name2,"/Aspect_",name2,d8ordinfstr2))
writeRaster(p,paste0("Aspect_",name,d8ordinfstr2),format = "GTiff",overwrite = TRUE)
setwd(paste0("G:/My Drive/Which DEM Paper/SWATHelperFiles/Rivers/",name2))
writeOGR(river,d8ordinfstr, layer = 'hold',driver="ESRI Shapefile",overwrite_layer = TRUE)