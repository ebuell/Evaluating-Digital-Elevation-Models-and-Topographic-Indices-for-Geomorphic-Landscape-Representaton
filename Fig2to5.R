#The purpose of this code is to generate figures 2-5 and S2-S5 in occompanying paper

if (!require("pacman")) install.packages("pacman")
pacman::p_load(raster,rgdal,googlesheets4,shapefiles)
library('RColorBrewer');library(corrplot); library("hydroGOF"); library("ggpubr")

#setting up for plotting
whichgraphic = 2 #=2 when you want to compare D8 spatial data; =1 when you want to compare Dinf dinf spatial data
whichspatial = "Slope_" #TIValue_ or TIClass_ or SCA_ or Slope_
zlimit = .3 #TIC=10; slp=0.3; ln(SCA)=10; TIV= 15
zlimithist = 2.5 #=zlimit for not Slope, =2.5 for Slope
filename = "Supplemental/FigS2_slpd8.pdf" #Fig5_TICdinf.pdf; Fig4_TIVdinf.pdf; Fig3_lnscadinf.pdf; Fig2_slpdinf.pdf
#Supplemental/FigS5_TICd8.pdf; Supplemental/FigS4_TIVd8.pdf; Supplemental/FigS3_lnscad8.pdf; Supplemental/FigS2_slpd8.pdf
#^change this
if(whichgraphic==2){d8ordinf = 1}
if(whichgraphic==1){d8ordinf = 2}
whichspatial2 = whichspatial
if(whichspatial=="Slope_"){whichspatialnum = 1}
if(whichspatial=="SCA_"){whichspatialnum = 2; whichspatial2 = "ln(SCA)"}
if(whichspatial=="TIValue_"){whichspatialnum = 3}
if(whichspatial=="TIClass_"){whichspatialnum = 4}


# Elyce's directory
#please note user will have to rearrange this
name3 = c("1arcsec","1_3arcsec","LIDAR_2010","LIDAR_2018")
name4 = c("1arcsec","1_3arcsec","2010LIDAR","2018LIDAR")
d8ordinfstr3 = c("_D8","_Dinf")



#create comparison rasters
a = NULL; b = NULL; c = NULL; d = NULL; e = NULL; f = NULL;
for (i in 1:6){
  #lower resolution is DEMchoice2
  if(i == 1){DEMchoice = 2; DEMchoice2 = 1}
  if(i == 2){DEMchoice = 3; DEMchoice2 = 1}
  if(i == 3){DEMchoice = 4; DEMchoice2 = 1}
  if(i == 4){DEMchoice = 3; DEMchoice2 = 2}
  if(i == 5){DEMchoice = 4; DEMchoice2 = 2}
  if(i == 6){DEMchoice = 4; DEMchoice2 = 3}
  if(DEMchoice<DEMchoice2){stop()}
  url = paste0("G:/My Drive/Which DEM Paper/TiffFiles/",whichspatial,name3[DEMchoice],"/",whichspatial,name3[DEMchoice],d8ordinfstr3[d8ordinf],"/",whichspatial,name4[DEMchoice],d8ordinfstr3[d8ordinf],".tif")
  url2 = paste0("G:/My Drive/Which DEM Paper/TiffFiles/",whichspatial,name3[DEMchoice2],"/",whichspatial,name3[DEMchoice2],d8ordinfstr3[d8ordinf],"/",whichspatial,name4[DEMchoice2],d8ordinfstr3[d8ordinf],".tif")
  rast=raster(url); rast2 = raster(url2);
  x = SpatialPoints(rast)
  rast2altered = raster::extract(rast2,x)
  rast2alteredrast = rast
  values(rast2alteredrast) = rast2altered
  if(whichspatial=="SCA_"){rast2alteredrast = log(rast2alteredrast); rast=log(rast)}
  if(i == 1){a = rast2alteredrast-rast}#;values(a) = round(values(a),1)}
  if(i == 2){b = rast2alteredrast-rast}#;values(b) = round(values(b),1)}
  if(i == 3){c = rast2alteredrast-rast}#;values(c) = round(values(c),1)}
  if(i == 4){d = rast2alteredrast-rast}#;values(d) = round(values(d),1)}
  if(i == 5){e = rast2alteredrast-rast}#;values(e) = round(values(e),1)}
  if(i == 6){f = rast2alteredrast-rast}#;values(f) = round(values(f),1)}
  print(i)
}

#create watershed boundaries
as1 = NULL; as13 = NULL; LI10 = NULL; LI18 = NULL;
for (i in 1:4){
  url = paste0("G:/My Drive/Which DEM Paper/TiffFiles/SCA_",name3[i],"/SCA_",name3[i],d8ordinfstr3[1],"/SCA_",name4[i],d8ordinfstr3[1],".tif")
  rast = raster(url); values(rast)[!is.na(values(rast))] = 0
  if(i ==1){as1 = rasterToPolygons(rast, dissolve=TRUE,na.rm = TRUE)}
  if(i ==2){rast = aggregate(rast,2);as13 = rasterToPolygons(rast, dissolve=TRUE,na.rm = TRUE)}
  if(i ==3){rast = aggregate(rast,10);LI10 = rasterToPolygons(rast, dissolve=TRUE,na.rm = TRUE)}
  if(i ==4){rast = aggregate(rast,20);LI18 = rasterToPolygons(rast, dissolve=TRUE,na.rm = TRUE)}
  print(i)
}

#set up colors for plotting
exty = c(4117200,4118590); extx =  c(549150,551470)
colorforTIC = c(brewer.pal(9,"Oranges")[9:1],"grey100",brewer.pal(9,"BuGn"))
color = brewer.pal(5,"Set1")
color = color[2:5]
colorforDEM = matrix("grey",10,1)
colorforDEM = (paste0(colorforDEM,c(ceiling(seq(0,90,length.out = 10)))))


#import slope data
usgs1as = raster(paste0("G:/My Drive/Which DEM Paper/TiffFiles/",whichspatial,name3[1],"/",whichspatial,name3[1],d8ordinfstr3[d8ordinf],"/",whichspatial,name4[1],d8ordinfstr3[d8ordinf],".tif"))
usgs13as = raster(paste0("G:/My Drive/Which DEM Paper/TiffFiles/",whichspatial,name3[2],"/",whichspatial,name3[2],d8ordinfstr3[d8ordinf],"/",whichspatial,name4[2],d8ordinfstr3[d8ordinf],".tif"))
LI2010 = raster(paste0("G:/My Drive/Which DEM Paper/TiffFiles/",whichspatial,name3[3],"/",whichspatial,name3[3],d8ordinfstr3[d8ordinf],"/",whichspatial,name4[3],d8ordinfstr3[d8ordinf],".tif"))
LI2018 = raster(paste0("G:/My Drive/Which DEM Paper/TiffFiles/",whichspatial,name3[4],"/",whichspatial,name3[4],d8ordinfstr3[d8ordinf],"/",whichspatial,name4[4],d8ordinfstr3[d8ordinf],".tif"))
if(whichspatialnum==2){
  values(usgs1as) = log(values(usgs1as))
  values(usgs13as) = log(values(usgs13as))
  values(LI2018) = log(values(LI2018))
  values(LI2010) = log(values(LI2010))
}

values(usgs1as)[which(values(usgs1as)>zlimit)] = zlimit
values(usgs13as)[which(values(usgs13as)>zlimit)] = zlimit
values(LI2010)[which(values(LI2010)>zlimit)] = zlimit
values(LI2018)[which(values(LI2018)>zlimit)] = zlimit


url="G:/My Drive/Which DEM Paper/SWATHelperFiles/Monitoring_Point.shp"
outlet = readOGR(url)
crs_project = crs(outlet) #check that we are in the correct coordinate system

#define units
zlimmin = 0

pdf(file=paste0("G:/My Drive/Which DEM Paper/PaperFigs/Final paper figs/",filename),width = 11.09,height = 8.92)

#create figure
par(mfrow = c(4,4))
par(mar = c(3,1,2,3.75))
#exty = c(4117000,4118800); extx =  c(549050,551600)
  
#1
plot(LI2018, xaxt = "n", yaxt = "n",ylim = exty,zlim = c(zlimmin,zlimit),xlim = extx,col = rev(colorforDEM),legend = FALSE,main =  paste("2018 LI ",gsub("_","",whichspatial2)),cex.main = 2,maxpixels=400000)
par(xpd=TRUE)
legend(x = 551450, y = 4118600,bty='n',
       legend = c(zlimit,"","","","","","","","",zlimmin),
       fill = colorforDEM,border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
par(xpd=FALSE)
  
#2
plot(f, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = c(-zlimit,zlimit),col = colorforTIC,legend = FALSE)
par(xpd=TRUE)
legend(x = 551450, y = 4118600,bty='n',
       legend = c(zlimit,"","","","","","","","","0","","","","","","","","",-zlimit),
       fill = rev(colorforTIC),border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
par(xpd=FALSE)
  
#3
plot(e, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = c(-zlimit,zlimit),col = colorforTIC,legend = FALSE)
par(xpd=TRUE)
legend(x = 551450, y = 4118600,bty='n',
       legend = c(zlimit,"","","","","","","","","0","","","","","","","","",-zlimit),
       fill = rev(colorforTIC),border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
par(xpd=FALSE)
  
#4
plot(c, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = c(-zlimit,zlimit),col = colorforTIC,legend = FALSE)
par(xpd=TRUE)
legend(x = 551450, y = 4118600,bty='n',
       legend = c(zlimit,"","","","","","","","","0","","","","","","","","",-zlimit),
       fill = rev(colorforTIC),border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
par(xpd=FALSE)
  
#5
hist(f, yaxt = "n",main = paste("\n SD:",sprintf("%.2f",round(sd(values(f),na.rm = TRUE),2)),"\n"),xatx = "n",xlim = c(-zlimithist,zlimithist))

#6
plot(LI2010, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = c(zlimmin,zlimit),col = rev(colorforDEM),legend = FALSE,main =  paste("2010 LI ",gsub("_","",whichspatial2)),cex.main = 2)
par(xpd=TRUE)
legend(x = 551450, y = 4118600,bty='n',
       legend = c(zlimit,"","","","","","","","",zlimmin),
       fill = colorforDEM,border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
par(xpd=FALSE)
  
#7
plot(d, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = c(-zlimit,zlimit),col = colorforTIC,legend = FALSE)
par(xpd=TRUE)
legend(x = 551450, y = 4118600,bty='n',
       legend = c(zlimit,"","","","","","","","","0","","","","","","","","",-zlimit),
       fill = rev(colorforTIC),border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
par(xpd=FALSE)
  
#8
plot(b, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = c(-zlimit,zlimit),col = colorforTIC,legend = FALSE)
par(xpd=TRUE)
legend(x = 551450, y = 4118600,bty='n',
      legend = c(zlimit,"","","","","","","","","0","","","","","","","","",-zlimit),
      fill = rev(colorforTIC),border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
par(xpd=FALSE)

#9
hist(e, yaxt = "n",main = paste("\n SD:",sprintf("%.2f",round(sd(values(e),na.rm = TRUE),2)),"\n"),xatx = "n",xlim = c(-zlimithist,zlimithist))
  
#10 
hist(d, yaxt = "n",main = paste("\n SD:",sprintf("%.2f",round(sd(values(d),na.rm = TRUE),2)),"\n"),xatx = "n",xlim = c(-zlimithist,zlimithist))
  
#11
plot(usgs13as, xaxt = "n", yaxt = "n",ylim = exty,zlim = c(zlimmin,zlimit),xlim = extx,col = rev(colorforDEM),legend = FALSE,main =  paste("USGS 1/3as ",gsub("_","",whichspatial2)),cex.main = 2)
par(xpd=TRUE)
legend(x = 551450, y = 4118600,bty='n',
       legend = c(zlimit,"","","","","","","","",zlimmin),
       fill = colorforDEM,border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
par(xpd=FALSE)
  
#12
plot(a, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = c(-zlimit,zlimit),col = colorforTIC,legend = FALSE)
par(xpd=TRUE)
legend(x = 551450, y = 4118600,bty='n',
       legend = c(zlimit,"","","","","","","","","0","","","","","","","","",-zlimit),
       fill = rev(colorforTIC),border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
par(xpd=FALSE)

#13
hist(c, yaxt = "n",main = paste("\n SD:", sprintf("%.2f",round(sd(values(c),na.rm = TRUE),2)),"\n"),xatx = "n",xlim = c(-zlimithist,zlimithist))

#14
hist(b, yaxt = "n",main = paste("\n SD:",sprintf("%.2f",round(sd(values(b),na.rm = TRUE),2)),"\n"),xatx = "n",xlim = c(-zlimithist,zlimithist))
  
#15
hist(a, yaxt = "n",main = paste("\n SD:",sprintf("%.2f",round(sd(values(a),na.rm = TRUE),2)),"\n"),xatx = "n",xlim = c(-zlimithist,zlimithist))
  
#16
plot(usgs1as, xaxt = "n", yaxt = "n",ylim = exty,zlim = c(zlimmin,zlimit),xlim = extx,col = rev(colorforDEM),legend = FALSE,main =  paste("USGS 1as ",gsub("_","",whichspatial2)),cex.main = 2)
par(xpd=TRUE)
legend(x = 551450, y = 4118600,bty='n',
       legend = c(zlimit,"","","","","","","","",zlimmin),
       fill = colorforDEM,border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
par(xpd=FALSE)

dev.off()

stop()
