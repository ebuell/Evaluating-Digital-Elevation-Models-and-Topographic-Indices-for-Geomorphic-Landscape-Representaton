#The purpose of this code is to generate the supplemental figures for the corresponding papaer
#S1 and S6

if (!require("pacman")) install.packages("pacman")
pacman::p_load(raster,rgdal,googlesheets4,shapefiles,rgeos)
library('RColorBrewer');library(corrplot); library("hydroGOF"); library("ggpubr")

# Elyce's directory
#please note user will have to rearrange this
basedir="G:/My Drive/Which DEM Paper/DEMFiles"
diroptions = c("1arcsec","1_3arcsec","LIDAR/2010","LIDAR/2018")
fileoptions = c("/1arcsec_arcproj_bilin.tif","/1_3arcsec_arcproj_bilin.tif","/2010LIDAR_arcproj_bilin.tif","/2018LIDAR_arcproj_bilin.tif")
name3 = c("1arcsec","1_3arcsec","LIDAR_2010","LIDAR_2018")
name4 = c("1arcsec","1_3arcsec","2010LIDAR","2018LIDAR")
d8ordinfstr3 = c("_D8","_Dinf")

#Paper figure: S6
construction = readOGR("G:/My Drive/Which DEM Paper/SWATHelperFiles/Construction_exclusion.shp")
url = paste0("G:/My Drive/Which DEM Paper/TiffFiles/TIClass_",name3[4],"/TIClass_",name3[4],d8ordinfstr3[1],"/TIClass_",name4[4],d8ordinfstr3[1],".tif")
rast=raster(url); rast = mask(rast,construction,inverse = TRUE); x = SpatialPoints(rast)
corspatialTIV = NULL
for (i in 1:8){
  #lower resolution is DEMchoice2
  if(i == 1){DEMchoice = 4; d8ordinf = 2}
  if(i == 2){DEMchoice = 4; d8ordinf = 1}
  if(i == 3){DEMchoice = 3; d8ordinf = 2}
  if(i == 4){DEMchoice = 3; d8ordinf = 1}
  if(i == 5){DEMchoice = 2; d8ordinf = 2}
  if(i == 6){DEMchoice = 2; d8ordinf = 1}
  if(i == 7){DEMchoice = 1; d8ordinf = 2}
  if(i == 8){DEMchoice = 1; d8ordinf = 1}
  url2 = paste0("G:/My Drive/Which DEM Paper/TiffFiles/TIValue_",name3[DEMchoice],"/TIValue_",name3[DEMchoice],d8ordinfstr3[d8ordinf],"/TIValue_",name4[DEMchoice],d8ordinfstr3[d8ordinf],".tif")
  rast2 = raster(url2); rast2 = mask(rast2,construction,inverse = TRUE)
  hold = raster::extract(rast2,x)
  if(i == 1){corspatialTIV$LI18Dinf = hold}
  if(i == 2){corspatialTIV$LI18D8 = hold}
  if(i == 3){corspatialTIV$LI10Dinf = hold}
  if(i == 4){corspatialTIV$LI10D8 = hold}
  if(i == 5){corspatialTIV$as13Dinf = hold}
  if(i == 6){corspatialTIV$as13D8 = hold}
  if(i == 7){corspatialTIV$as1Dinf = hold}
  if(i == 8){corspatialTIV$as1D8 = hold}
  print(i)
}
    
TIVbrs2018LI = c(3.536021,4.035139,4.436527,4.800556,5.152543,5.521609,5.946593,6.513813,7.565863)
TIVbrs2010LI = c(3.5086740,3.9886785,4.3935137,4.7687975,5.1389298,5.5206985,5.9472978,6.4836664,7.4608381)
TIVbrs13as = c(5.664676,6.014089,6.280968,6.518596,6.769484,7.048799,7.398143,7.842223,8.681655)
TIVbrs1as = c(6.102894,6.40374,6.626421,6.880704,7.149188,7.425146,7.786315,8.357537,9.792187)
limmin = 2; limmax = 13
offset=.4
binsize=70
xoff = 0
df = data.frame(x=corspatialTIV$as1Dinf,y=corspatialTIV$LI18Dinf)
hold <- ggplot()+ggtitle("\n TIV density mapping for\n coarsest and finest rasters\n and corresponding\n TIC breaks")+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid =  element_blank(),panel.background =  element_blank(),
        legend.position = c(.2,.6))+
  xlab("")+ylab("")+xlim(0,0)+ylim(0,0)
c<-ggplot(df, aes(x = x, y = y))+theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid =  element_blank(),panel.background =  element_blank(),legend.position = c(.9, .9))+
  geom_hex(bins=binsize,na.rm = TRUE,show.legend = TRUE)+scale_fill_gradient2()+
  geom_hline(yintercept = TIVbrs2018LI,col='cyan4')+geom_vline(xintercept = TIVbrs1as,col='cyan4')+
  xlab("")+ylab("")+xlim(limmin-offset,limmax)+ylim(limmin-offset,limmax)
#c<-c+guides(fill=guide_legend("Number of rasters\nmapped to the\nintersection of\n2018LI Dinf TIV and\n1as Dinf TIV"))
#c<-c+geom_point(x = TIVbrs1as[1], y = TIVbrs2018LI[1],col = 'darkred')
#c<-c+geom_point(x = TIVbrs1as[2], y = TIVbrs2018LI[2],col = 'darkred')
#c<-c+geom_point(x = TIVbrs1as[3], y = TIVbrs2018LI[3],col = 'darkred')
#c<-c+geom_point(x = TIVbrs1as[4], y = TIVbrs2018LI[4],col = 'darkred')
#c<-c+geom_point(x = TIVbrs1as[5], y = TIVbrs2018LI[5],col = 'darkred')
#c<-c+geom_point(x = TIVbrs1as[6], y = TIVbrs2018LI[6],col = 'darkred')
#c<-c+geom_point(x = TIVbrs1as[7], y = TIVbrs2018LI[7],col = 'darkred')
#c<-c+geom_point(x = TIVbrs1as[8], y = TIVbrs2018LI[8],col = 'darkred')
#c<-c+geom_point(x = TIVbrs1as[9], y = TIVbrs2018LI[9],col = 'darkred')

a <- ggplot(df, aes(x=y))+ylab("Density") + coord_flip()+scale_y_reverse()+theme(axis.ticks = element_blank(),panel.grid =  element_blank(),panel.background =  element_blank())+
  geom_density(color="black",size=1.5)+xlab("2018LI Dinf TIV")+xlim(limmin,limmax)+geom_vline(xintercept=TIVbrs2018LI,col='cyan4')

d<-ggplot(df, aes(x=x)) + ylab("Density")+geom_density(color="black",size=1.5)+scale_color_manual(values = c("TIC numeric breaks" = "cyan4", "Locations of equivalence between 2018LI Dinf and 1as Dinf TICs" = "darkred"))+xlab("1as Dinf TIV")+xlim(limmin,limmax)+
  geom_vline(xintercept = TIVbrs1as,col='cyan4')+theme(axis.ticks = element_blank(),panel.grid =  element_blank(),panel.background =  element_blank(),legend.position = c(.95, .95))
ggarrange(hold,d,a,c,ncol = 2,nrow = 2,widths = c(1,3),heights = c(1,3))


#paper figure S1
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
  folder = paste(diroptions[DEMchoice],sep = "")
  url=paste0(basedir,"/",folder,fileoptions[DEMchoice])
  rast=raster(url); if(DEMchoice==3||DEMchoice==4){values(rast) = values(rast)*0.3048}
  folder = paste(diroptions[DEMchoice2],sep = "")
  url=paste0(basedir,"/",folder,fileoptions[DEMchoice2])
  rast2 = raster(url); if(DEMchoice2==3||DEMchoice2==4){values(rast2) = values(rast2)*0.3048}
  x = SpatialPoints(rast)
  rast2altered = raster::extract(rast2,x)
  rast2alteredrast = rast
  values(rast2alteredrast) = rast2altered
  if(i == 1){a = rast2alteredrast-rast}#;values(a) = round(values(a),1)}
  if(i == 2){b = rast2alteredrast-rast}#;values(b) = round(values(b),1)}
  if(i == 3){c = rast2alteredrast-rast}#;values(c) = round(values(c),1)}
  if(i == 4){d = rast2alteredrast-rast}#;values(d) = round(values(d),1)}
  if(i == 5){e = rast2alteredrast-rast}#;values(e) = round(values(e),1)}
  if(i == 6){f = rast2alteredrast-rast}#;values(f) = round(values(f),1)}
  print(i)
}

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
m <- do.call(bind, list(as1,as13)); m = gUnaryUnion(m); a = mask(a,m)
m <- do.call(bind, list(as1,LI10)); m = gUnaryUnion(m); b = mask(b,m)
m <- do.call(bind, list(as1,LI18)); m = gUnaryUnion(m); c = mask(c,m)
m <- do.call(bind, list(as13,LI10)); m = gUnaryUnion(m); d = mask(d,m)
m <- do.call(bind, list(as13,LI18)); m = gUnaryUnion(m); e = mask(e,m)
m <- do.call(bind, list(LI10,LI18)); m = gUnaryUnion(m); f = mask(f,m)

exty = c(4117200,4118590); extx =  c(549150,551470)
colorforDEM = matrix("grey",17,1)
colorforDEM = paste0(colorforDEM,c(ceiling(seq(0,100,length.out = 8)),ceiling(seq(86,0,length.out = 7))))
color = brewer.pal(5,"Set1")
color = color[2:5]

values(a)[which(values(a)>8)] = 8;values(a)[which(values(a)< -8)] = -8;
values(b)[which(values(b)>8)] = 8;values(b)[which(values(b)< -8)] = -8;
values(c)[which(values(c)>8)] = 8;values(c)[which(values(c)< -8)] = -8;
values(d)[which(values(d)>8)] = 8;values(d)[which(values(d)< -8)] = -8;
values(e)[which(values(e)>8)] = 8;values(e)[which(values(e)< -8)] = -8;
values(f)[which(values(f)>8)] = 8;values(f)[which(values(f)< -8)] = -8;
zlimit = c(-8,8)

dem2018 = raster("G:/My Drive/Which DEM Paper/DEMFiles/LIDAR/2018/2018LIDAR_arcproj_bilin.tif")
maskdem = raster("G:/My Drive/Which DEM Paper/TiffFiles/TIClass_LIDAR_2018/TIClass_LIDAR_2018_D8/TIClass_2018LIDAR_D8.tif")
dem2018 = mask(dem2018,maskdem)
dem2018 = dem2018*.3048

dem2010 = raster("G:/My Drive/Which DEM Paper/DEMFiles/LIDAR/2010/2010LIDAR_arcproj_bilin.tif")
maskdem = raster("G:/My Drive/Which DEM Paper/TiffFiles/TIClass_LIDAR_2010/TIClass_LIDAR_2010_D8/TIClass_2010LIDAR_D8.tif")
dem2010 = mask(dem2010,maskdem)
dem2010 = dem2010*.3048

dem13as = raster("G:/My Drive/Which DEM Paper/DEMFiles/1_3arcsec/1_3arcsec_arcproj_bilin.tif")
maskdem = raster("G:/My Drive/Which DEM Paper/TiffFiles/TIClass_1_3arcsec/TIClass_1_3arcsec_D8/TIClass_1_3arcsec_D8.tif")
dem13as = mask(dem13as,maskdem)

dem1as = raster("G:/My Drive/Which DEM Paper/DEMFiles/1arcsec/1arcsec_arcproj_bilin.tif")
maskdem = raster("G:/My Drive/Which DEM Paper/TiffFiles/TIClass_1arcsec/TIClass_1arcsec_D8/TIClass_1arcsec_D8.tif")
dem1as = mask(dem1as,maskdem)

par(mfrow = c(4,4))
par(mar = c(2.5,1,2,4))
  
#1
plot(dem2018, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,col = rev(colorforDEM[1:7]),maxpixels=500000,main = "2018 LiDAR DEM",legend=FALSE)
par(xpd=TRUE)
legend(x = 551450, y = 4118600,bty='n',
       legend = c(paste0(round(max(values(dem2018),na.rm = TRUE),-1),"m"),"","","","","",paste0(round(min(values(dem2018),na.rm = TRUE),-1),"m")),
       fill = (colorforDEM[1:7]),border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
par(xpd=FALSE)

#2
plot(f, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = zlimit,col = colorforDEM,maxpixels=500000,legend=FALSE)
par(xpd=TRUE)
legend(x = 551500, y = 4118600,bty='n',
       legend = c("8 m","","","","","","","","0","","","","","","","","-8 m"),
       fill = colorforDEM,border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
par(xpd=FALSE)
lines(LI18, col = color[1])
lines(LI10, col = color[2])

#3
plot(e, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = zlimit,col = colorforDEM,maxpixels=500000,legend=FALSE)
par(xpd=TRUE)
legend(x = 551500, y = 4118600,bty='n',
       legend = c("8 m","","","","","","","","0","","","","","","","","-8 m"),
       fill = colorforDEM,border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
par(xpd=FALSE)
lines(LI18, col = color[1])
lines(as13, col = color[3])
  
#4
plot(c, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = zlimit,col = colorforDEM,maxpixels=500000,legend=FALSE)
par(xpd=TRUE)
legend(x = 551500, y = 4118600,bty='n',
       legend = c("8 m","","","","","","","","0","","","","","","","","-8 m"),
       fill = colorforDEM,border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
par(xpd=FALSE)
lines(LI18, col = color[1])
lines(as1, col = color[4])
  
#5
hist(f, yaxt = "n",main = "",xlim = c(-7.5,7.5))
  
#6
plot(dem2010, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,col = rev(colorforDEM[1:7]),maxpixels=500000,main = "2010 LiDAR DEM",legend=FALSE)
par(xpd=TRUE)
legend(x = 551450, y = 4118600,bty='n',
       legend = c(paste0(round(max(values(dem2010),na.rm = TRUE),-1),"m"),"","","","","",paste0(round(min(values(dem2010),na.rm = TRUE),-1),"m")),
       fill = (colorforDEM[1:7]),border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
par(xpd=FALSE)
  
#7
plot(d, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = zlimit,col = colorforDEM,maxpixels=500000,legend=FALSE)
par(xpd=TRUE)
legend(x = 551500, y = 4118600,bty='n',
       legend = c("8 m","","","","","","","","0","","","","","","","","-8 m"),
       fill = colorforDEM,border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
par(xpd=FALSE)
lines(as13, col = color[3])
lines(LI10, col = color[2])
  
#8
plot(b, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = zlimit,col = colorforDEM,maxpixels=500000,legend=FALSE)
par(xpd=TRUE)
legend(x = 551500, y = 4118600,bty='n',
   legend = c("8 m","","","","","","","","0","","","","","","","","-8 m"),
   fill = colorforDEM,border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
par(xpd=FALSE)
lines(as1, col = color[4])
lines(LI10, col = color[2])
  
#9
hist(e, yaxt = "n",main = "",xlim = c(-7.5,7.5))
  
#10 
hist(d, yaxt = "n",main = "",xlim = c(-7.5,7.5))
  
#11
plot(dem13as, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,col = rev(colorforDEM[1:7]),maxpixels=500000,main = "USGS 1/3 arcsec DEM",legend=FALSE)
par(xpd=TRUE)
legend(x = 551450, y = 4118600,bty='n',
       legend = c(paste0(round(max(values(dem13as),na.rm = TRUE),-1),"m"),"","","","","",paste0(round(min(values(dem13as),na.rm = TRUE),-1),"m")),
       fill = (colorforDEM[1:7]),border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
par(xpd=FALSE)
  
#12
plot(a, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,zlim = zlimit,col = colorforDEM,maxpixels=500000,legend=FALSE)
par(xpd=TRUE)
legend(x = 551500, y = 4118600,bty='n',
       legend = c("8 m","","","","","","","","0","","","","","","","","-8 m"),
       fill = colorforDEM,border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
par(xpd=FALSE)
lines(as1, col = color[4])
lines(as13, col = color[3])
  
#13
hist(c, yaxt = "n",main = "",xlim = c(-7.5,7.5))

#14
hist(b, yaxt = "n",main = "",xlim = c(-7.5,7.5))

#15
hist(a, yaxt = "n",main = "",xlim = c(-7.5,7.5))

#16
plot(dem1as, xaxt = "n", yaxt = "n",ylim = exty,xlim = extx,col = rev(colorforDEM[1:7]),maxpixels=500000,main = "USGS 1 arcsec DEM",legend=FALSE)
par(xpd=TRUE)
legend(x = 551450, y = 4118600,bty='n',
       legend = c(paste0(round(max(values(dem1as),na.rm = TRUE),-1),"m"),"","","","","",paste0(round(min(values(dem1as),na.rm = TRUE),-1),"m")),
       fill = (colorforDEM[1:7]),border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
par(xpd=FALSE)
