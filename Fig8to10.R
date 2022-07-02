#The purpose of this code is to generate bootstrapped predictions and distributed soil maps for 6
#soil properties and plot ssurgo soils and predicted soils against measured
#generate figures 8 to 10 in occompanying paper


#### bootstrap ####
library('sp')
library("hydroGOF")
library('googlesheets4')
library("RColorBrewer")

##### set up dataframes #####
#soils dataframe
PhysResults_sheetid="1HnsQjb8j9kUlBrlbPuz20czVoQvGtQwy3tmZ3vyULfc"
PhysicalProperties_Meas <- read_sheet(PhysResults_sheetid, sheet = "PhysicalResults",col_types = "?_dd_d_d__d___dd___d______")
names(PhysicalProperties_Meas) = c("SSID", "Lat", "Long", "DepthBa", "DepthBt", "OMA", "ClayA","OMBa", "ClayBa")
#Depth of A, Ba and Bt
PhysicalProperties_Meas$A_thickness = PhysicalProperties_Meas$DepthBa
PhysicalProperties_Meas$A_thickness[is.na(PhysicalProperties_Meas$A_thickness)] = PhysicalProperties_Meas$DepthBt[is.na(PhysicalProperties_Meas$A_thickness)]
PhysicalProperties_Meas$Ba_thickness = PhysicalProperties_Meas$DepthBt - PhysicalProperties_Meas$A_thickness
PhysicalProperties_Meas$Ba_thickness[PhysicalProperties_Meas$Ba_thickness==0] = NA

#Other DFs
TIClasses <- read_sheet(PhysResults_sheetid, sheet = "TIClass")
Slopes <- read_sheet(PhysResults_sheetid, sheet = "Slope")
SCA <- read_sheet(PhysResults_sheetid, sheet = "SCA")
TIValues <- read_sheet(PhysResults_sheetid, sheet = "TIValue")

##### define training and boostrapping criteria #####
N = 10000
############

##### import SSURGO soils #####
if (!require("pacman")) install.packages("pacman")
pacman::p_load(raster,soilDB)
library('googlesheets4'); library('soilDB'); library('rgdal'); library('aqp')

#if this code has not been run in a couple of days go to https://websoilsurvey.sc.egov.usda.gov and get a new download link
url = "https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/szrl5mkohklzlqayp0kuslzw/wss_aoi_2022-07-01_10-00-35.zip"
download.file(url,"mysoil.zip")
unzip("mysoil.zip")
list.files()
list.files("wss_aoi_2021-04-20_12-34-40/spatial/",pattern = "shp")
list.files("wss_aoi_2021-04-20_12-34-40/tabular/")

# Elyce's directory
bboxfilepath ="G:/My Drive/Which DEM Paper/DEMFiles/LIDAR/2010/2010LIDAR_arcproj_bilin.tif"

#Get extension of soils map
rast=raster(bboxfilepath)
#Projection change
proj4_utm = paste0("+proj=utm +zone=", 17, " +datum=WGS84 +units=m +no_defs")
proj4_ll = "+proj=longlat"; crs_ll=CRS(proj4_ll); crs_utm=CRS(proj4_utm)
dem_utm <- projectRaster(rast,crs = crs_ll)
#Define bbox
mybbox = bbox(dem_utm)

# Get those soils!
mysoil = mapunit_geom_by_ll_bbox(mybbox,source = 'sda')
proj4string(mysoil) = proj4_ll

# First associate mukey with cokey from component
mukey_statement = format_SQL_in_statement(unique(mysoil$mukey))
q_mu2co = paste("SELECT mukey,cokey FROM component WHERE mukey IN ", mukey_statement, sep="")
mu2co = SDA_query(q_mu2co)
# Second associate cokey with ksat_r,awc_r,hzdepb_r from chorizon
cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
q_co2ch = paste("SELECT cokey,chkey,hzname,OM_r,claytotal_r,hzdept_r,hzdepb_r  FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
co2ch = SDA_query(q_co2ch)
mu2ch=merge(mu2co,co2ch)
q.2 = paste("SELECT cokey,comppct_r FROM legend INNER JOIN mapunit ON mapunit.lkey = legend.lkey
       LEFT OUTER JOIN component ON component.mukey = mapunit.mukey WHERE cokey IN ", cokey_statement, sep="")
co2dom = SDA_query(q.2)
mu2ch=merge(mu2ch,co2dom)

mu2chdomH1 = data.frame(mukey = unique(mu2ch$mukey[which(mu2ch$hzname=="H1")]),cokey=NA,OM_r=NA,
                        claytotal_r = NA,hzdepb_r = NA,comppct_r=NA)
mu2chdomH2 = data.frame(mukey = unique(mu2ch$mukey[which(mu2ch$hzname=="H2")]),cokey=NA,OM_r=NA,
                        claytotal_r = NA,hzdept_r = NA,hzdepb_r = NA,comppct_r=NA)

for(i in 1:length(mu2chdomH1$mukey)){
  ind = which(mu2ch$mukey==mu2chdomH1$mukey[i])
  ind2 = ind[which(mu2ch$hzname[ind]=="H1")]
  ind3 = ind2[which.max(mu2ch$comppct_r[ind2])]
  mu2chdomH1$cokey[i]=mu2ch$cokey[ind3]
  mu2chdomH1$OM_r[i]=mu2ch$OM_r[ind3]
  mu2chdomH1$claytotal_r[i]=mu2ch$claytotal_r[ind3]
  mu2chdomH1$hzdepb_r[i]=mu2ch$hzdepb_r[ind3]
  mu2chdomH1$comppct_r[i] = mu2ch$comppct_r[ind3]
}
mu2chdomH1$depth = mu2chdomH1$hzdepb_r

for(i in 1:length(mu2chdomH2$mukey)){
  ind = which(mu2ch$mukey==mu2chdomH2$mukey[i])
  ind2 = ind[which(mu2ch$hzname[ind]=="H2")]
  ind3 = ind2[which.max(mu2ch$comppct_r[ind2])]
  mu2chdomH2$cokey[i]=mu2ch$cokey[ind3]
  mu2chdomH2$OM_r[i]=mu2ch$OM_r[ind3]
  mu2chdomH2$claytotal_r[i]=mu2ch$claytotal_r[ind3]
  mu2chdomH2$hzdept_r[i]=mu2ch$hzdept_r[ind3]
  mu2chdomH2$hzdepb_r[i]=mu2ch$hzdepb_r[ind3]
  mu2chdomH2$comppct_r[i] = mu2ch$comppct_r[ind3]
}
mu2chdomH2$depth = mu2chdomH2$hzdepb_r-mu2chdomH2$hzdept_r

######create SSURGO rasters
#generate OM map from SSURGO
TIV13as = raster("G:/My Drive/Which DEM Paper/TiffFiles/TIValue_1_3arcsec/TIValue_1_3arcsec_Dinf/TIValue_1_3arcsec_Dinf.tif")
x = SpatialPoints(TIV13as,proj4string = crs_utm)

#Athickness named a;  
#AOM named b;        
#AClay names c;      
#Bathickness named d;
#BaOM named e;       
#BaClay named f;   

#set up
a = TIV13as; values(a)[!is.na(values(a))] = 0
mysoil2 = spTransform(mysoil,CRSobj = crs_utm)
values(a) = as.numeric(over(x, mysoil2)$mukey)
values(a)[is.na(values(TIV13as))] = NA
b = a; c = a; d = a; e = a; f = a


#Athick
for(i in 1:length(mu2chdomH1$mukey)){
  ind = which(values(a)==mu2chdomH1$mukey[i])
  values(a)[ind] = mu2chdomH1$depth[i]
}

#AOM
for(i in 1:length(mu2chdomH1$mukey)){
  ind = which(values(b)==mu2chdomH1$mukey[i])
  values(b)[ind] = mu2chdomH1$OM_r[i]/100
}

#Aclay
for(i in 1:length(mu2chdomH1$mukey)){
  ind = which(values(c)==mu2chdomH1$mukey[i])
  values(c)[ind] = mu2chdomH1$claytotal_r[i]/100
}

#Bathick
for(i in 1:length(mu2chdomH2$mukey)){
  ind = which(values(d)==mu2chdomH2$mukey[i])
  values(d)[ind] = mu2chdomH2$depth[i]
}

#BaOM
for(i in 1:length(mu2chdomH2$mukey)){
  ind = which(values(e)==mu2chdomH2$mukey[i])
  values(e)[ind] = mu2chdomH2$OM_r[i]/100
}

#Baclay
for(i in 1:length(mu2chdomH2$mukey)){
  ind = which(values(f)==mu2chdomH2$mukey[i])
  values(f)[ind] = mu2chdomH2$claytotal_r[i]/100
}

######create relevant rasters for bootstrapping fig creation (figures 9 and 10)
#import relevant rasters
Slope2010d8 = raster("G:/My Drive/Which DEM Paper/TiffFiles/Slope_LIDAR_2010/Slope_LIDAR_2010_D8/Slope_2010LIDAR_D8.tif")
Slope13asd8 = raster("G:/My Drive/Which DEM Paper/TiffFiles/Slope_1_3arcsec/Slope_1_3arcsec_D8/Slope_1_3arcsec_D8.tif")
Slope1asdinf = raster("G:/My Drive/Which DEM Paper/TiffFiles/Slope_1arcsec/Slope_1arcsec_Dinf/Slope_1arcsec_Dinf.tif")
Slope1asd8 = raster("G:/My Drive/Which DEM Paper/TiffFiles/Slope_1arcsec/Slope_1arcsec_D8/Slope_1arcsec_D8.tif")

lnSCA2018d8 = log(raster("G:/My Drive/Which DEM Paper/TiffFiles/SCA_LIDAR_2018/SCA_LIDAR_2018_D8/SCA_2018LIDAR_D8.tif"))
lnSCA2010dinf = log(raster("G:/My Drive/Which DEM Paper/TiffFiles/SCA_LIDAR_2010/SCA_LIDAR_2010_Dinf/SCA_2010LIDAR_Dinf.tif"))
lnSCA13asd8 = log(raster("G:/My Drive/Which DEM Paper/TiffFiles/SCA_1_3arcsec/SCA_1_3arcsec_D8/SCA_1_3arcsec_D8.tif"))
lnSCA1asdinf = log(raster("G:/My Drive/Which DEM Paper/TiffFiles/SCA_1arcsec/SCA_1arcsec_Dinf/SCA_1arcsec_Dinf.tif"))

#Athickness named a;  a1: slp1asdinfas2010li, a2: ln(sca2010lidinf)
#AOM named b;         b1: slp1asd8as2010li;   b2: ln(sca2010lidinf)
#AClay names c;       c1: slp1asdinf,         c2: ln(sca1asdinf)
#Bathickness named d; d1: slp13asd8as2018li,  d2: ln(sca2018lid8)
#BaOM named e;        e1: slp1asd8as2018li,   e2: ln(sca2018lid8)
#BaClay named f;      f1: slp2010lid8,        f2: ln(sca13asd8)as2010li

a2 = lnSCA2010dinf; b2 = a2
c1 = Slope1asdinf; c2 = lnSCA1asdinf
d2 = lnSCA2018d8; e2 = d2
f1 = Slope2010d8

sp2010li = SpatialPoints(lnSCA2010dinf,proj4string = crs_utm)
sp2018li = SpatialPoints(lnSCA2018d8,proj4string = crs_utm)

#Athickness named a
a1 = a2
values(a1)[which(!is.na(values(a1)))] = 0
values(a1) = extract(Slope1asdinf,sp2010li)

#AOM named b
b1 = b2
values(b1)[which(!is.na(values(b1)))] = 0
values(b1) = extract(Slope1asd8,sp2010li)

#AClay named c -> same resolutions skip

#Bathickness named d
d1 = d2
values(d1)[which(!is.na(values(d1)))] = 0
values(d1) = extract(Slope13asd8,sp2018li)

#BaOM named e
e1 = e2
values(e1)[which(!is.na(values(e1)))] = 0
values(e1) = extract(Slope1asd8,sp2018li)

#BaClay named f
f2 = f1
values(f2)[which(!is.na(values(f2)))] = 0
values(f2) = extract(lnSCA13asd8,sp2010li)

#seq for letters
letterseq = c("a","b","c","d","e","f")

#rast zlims
zlimits = rbind(c(0,40),
                c(0,.1),
                c(0,.2),
                c(0,40),
                c(0,.05),
                c(0,.3))
colorforplots = brewer.pal(9,"YlGnBu")[2:9]
colorforplots2 = c(brewer.pal(9,"YlGnBu")[c(2,4,6,8:9)],brewer.pal(9,"Reds")[2:9])
colorforplots4 = c(brewer.pal(9,"YlGnBu")[c(2,6,9)],brewer.pal(9,"Reds")[1:9])
xlabtext = c("A horizon thickness",
             "A horizon organic matter",
             "A horizon clay content",
             "Ba horizon thickness",
             "Ba horizon organic matter",
             "Ba horizon clay content")
unitstext = c("cm","kg/kg","kg/kg","cm","kg/kg","kg/kg")

#stats holder
fitstats = data.frame(iter = 1:N, trainx = NA, trainsig = NA, testx = NA,testsig = NA,
                      finestresdinf = NA,finestresd8 = NA,
                      bestguesstrainxbar = NA, bestguesstrainsd = NA, bestguesstestxbar = NA, bestguesstestsd = NA,
                      bestguesscoefp = NA, bestguesscorrtrain=NA, bestguesscorrtest=NA)
forplottcoefvals = data.frame(iter = 1:N,athickc = NA, athickb1 = NA, athickb2 = NA,
                              aOMc = NA, aOMb1 = NA, aOMb2 = NA,
                              aclayc = NA, aclayb1 = NA, aclayb2 = NA,
                              bathickc = NA, bathickb1 = NA, bathickb2 = NA,
                              baOMc = NA, baOMb1 = NA, baOMb2 = NA,
                              baclayc = NA, baclayb1 = NA, baclayb2 = NA)
forplottcoefvalsind = rbind(c(2:4),c(5:7),c(8:10),c(11:13),c(14:16),c(17:19))

spatial = data.frame(slp2010d8 = Slopes$Slope_LIDAR_2010_D8,slp13asd8 = Slopes$Slope_1_3arcsec_D8,slp1asd8 = Slopes$Slope_1arcsec_D8, slp1asdinf = Slopes$Slope_1arcsec_Dinf,
                     lnsca2018d8 = log(SCA$SCA_LIDAR_2018_D8),lnsca2018dinf = log(SCA$SCA_LIDAR_2018_Dinf),lnsca2010dinf = log(SCA$SCA_LIDAR_2010_Dinf),
                     lnsca13asd8 = log(SCA$SCA_1_3arcsec_D8),lnsca1asdinf = log(SCA$SCA_1arcsec_Dinf))
spatialind = rbind(c(4,7),
                   c(3,7),
                   c(4,9),
                   c(2,5),
                   c(4,6),
                   c(1,8))
rownames(spatialind) = c("athick","oma","claya","bathick","omba","clayba")
colnames(spatialind) = c("Slope","lnSCA")
spatialind = data.frame(spatialind)

#set up soils df of interest
soils = data.frame(SSID = unlist(PhysicalProperties_Meas$SSID),
                   Athick = 0,OMA = PhysicalProperties_Meas$OMA,ClayA = PhysicalProperties_Meas$ClayA,
                   Bathick = 0,OMBa = PhysicalProperties_Meas$OMBa,ClayBa = PhysicalProperties_Meas$ClayBa)
soils$Athick = as.numeric(PhysicalProperties_Meas$DepthBa)
soils$Athick[is.na(soils$Athick)] = PhysicalProperties_Meas$DepthBt[is.na(soils$Athick)]
soils$Bathick = PhysicalProperties_Meas$DepthBt - as.numeric(unlist(PhysicalProperties_Meas$DepthBa))



##### Run bootstrap ####
#set up bootstraph dfs
dataforbootstrap = data.frame(TIV1as = TIValues$TIValue_1arcsec_Dinf,TIC1as = TIClasses$TIClass_1arcsec_Dinf,
                              TIV2010 = TIValues$TIValue_LIDAR_2010_Dinf, TIC2010 = TIClasses$TIClass_LIDAR_2010_Dinf,
                              TIVmixed = log((SCA$SCA_LIDAR_2010_Dinf+1)/(Slopes$Slope_1arcsec_Dinf+0.00001)),
                              bestSCA = log(SCA$SCA_LIDAR_2010_Dinf),bestslp = Slopes$Slope_1arcsec_Dinf)


#generate figure 9
pdf(file = "G:/My Drive/Which DEM Paper/PaperFigs/Final paper figs/Fig9_Adistsoil.pdf",width = 11.09,height = 8.92)
par(mfrow = c(3,3),oma=c(0,0,0,0))
#bootstrap
for(j in 1:3){
  indnna = !is.na(soils[,j+1])
  soilsdata = soils[indnna,j+1]
  trainnum = length(soilsdata)/2
  for(i in 1:N){
    #define subsamples
    subsample = sample(1:length(soilsdata),trainnum,replace=FALSE)
    trainingsample = soils[subsample,j+1]
    subtest = !(1:length(soilsdata) %in% subsample)
    testsample = soilsdata[subtest]
    
    #define statistics for subsamples
    fitstats$trainx[i] = mean(trainingsample)
    fitstats$testx[i] = mean(testsample)
    fitstats$trainsig[i] = sd(trainingsample)
    fitstats$testsig[i] = sd(testsample)
    
    #select data for bootstrap
    dataforbootstrap = data.frame(bestSCA = spatial[indnna,spatialind$lnSCA[j]],bestslp = spatial[indnna,spatialind$Slope[j]])
     
    
    #for finest res D8
    df = data.frame(slplid8 = Slopes$Slope_LIDAR_2018_D8[indnna],lnscalid8 = log(SCA$SCA_LIDAR_2018_D8)[indnna])
    linearmodel = lm(trainingsample~.,data = df[subtest,])
    #predict values
    testvals = predict(linearmodel,newdata = df[subsample,])
    fitstats$finestresd8[i] = cor(testvals,testsample)
    
    #for finest res Dinf
    df = data.frame(slplidinf = Slopes$Slope_LIDAR_2018_Dinf[indnna],lnscalidinf = log(SCA$SCA_LIDAR_2018_Dinf)[indnna])
    linearmodel = lm(trainingsample~.,data = df[subtest,])
    #predict values
    testvals = predict(linearmodel,newdata = df[subsample,])
    fitstats$finestresdinf[i] = cor(testvals,testsample)
    
    #for best guess
    linearmodelbest = lm(trainingsample~.,data = dataforbootstrap[subsample,])
    fitstats$bestguesstrainxbar[i] = mean(linearmodelbest$fitted.values)
    fitstats$bestguesstrainsd[i] = sd(linearmodelbest$fitted.values)
    fitstats$bestguessinterp[i] = summary(linearmodelbest)$coefficients[1,4]
    fitstats$bestguessp[i] = pf(summary(linearmodelbest)$fstatistic[1],summary(linearmodelbest)$fstatistic[2],summary(linearmodelbest)$fstatistic[3],lower.tail = FALSE)
    fitstats$bestguesscorrtrain[i] = sqrt(summary(linearmodelbest)$r.squared)
    #predict values
    testvals = predict(linearmodelbest,newdata = dataforbootstrap[subtest,])
    fitstats$bestguesstestxbar[i] = mean(testvals)
    fitstats$bestguesstestsd[i] = sd(testvals)
    fitstats$bestguesscorrtest[i] = cor(testvals,testsample)
    
    forplottcoefvals[i,forplottcoefvalsind[j,]] = linearmodelbest$coefficients
    
  }
  colorforplots3 = colorforplots
  if(j==3){values(c)[which(values(c)>.2)] = .2}
  
  #plot densities 
  par(mar = c(4,4,4,1))
  plot(density(fitstats$bestguesscorrtest),col = 'darkgoldenrod3',xlab = "Correlation",ylab="Density",bty='n',main=paste0("\n(A",j,") ",xlabtext[j]),xlim = c(-1,1),lwd=3,ylim = c(0,2.75))
   lines(density(fitstats$finestresd8),col = 'coral3',lwd=3)
   lines(density(fitstats$finestresdinf),col = 'darkcyan',lwd=3)
   legend('topleft',legend = c("Multivariate Regression","Finest res D8 slp&sca","Finest res Dinf slp&sca"),
          col = c("darkgoldenrod3","coral3","darkcyan"),lwd=3,bty = 'n')
  
   
   #plot proposed soil dist
   par(mar=c(4,2.5,4,5))
   ylimit = c(4117200,4118540); xlimit = c(549150,551500)
   slope = get(paste0(letterseq[j],"1"))
   lnsca = get(paste0(letterseq[j],"2"))
   rast = mean(forplottcoefvals[,forplottcoefvalsind[j,][1]])+mean(forplottcoefvals[,forplottcoefvalsind[j,][2]])*lnsca+mean(forplottcoefvals[,forplottcoefvalsind[j,][3]])*slope
   values(rast)[which(values(rast)<0)] = 0
   if(j==1){values(rast)[which(values(rast)>40)]=40}
   plot(rast,main = paste0("\n (B",j,") ",xlabtext[j]),bty='n',xaxt='n',yaxt='n',xlim=xlimit,ylim=ylimit,zlim=zlimits[j,],col=colorforplots3,legend=FALSE)
   if(j==1){
     par(xpd=TRUE)
     text(550340,4119040,"Multivariate TIV Regression",cex=2)
     par(xpd=FALSE)
   }
   #plot ssurgo soil dist
   par(mar=c(4,0,4,8.5))
   plot(get(paste0(letterseq[j])),main = paste0("\n (C",j,") ",xlabtext[j]),bty='n',xaxt='n',yaxt='n',xlim=xlimit,ylim=ylimit,zlim=zlimits[j,],col=colorforplots3,legend=FALSE)
   
   
   #if(j==1){textunitoffsetx = 700}else{textunitoffsetx = 850}
   par(xpd=TRUE)
   #text(xlimit[2]+textunitoffsetx,ylimit[2]-260,unitstext[j])
   #text(xlimit[2]+textunitoffsetx,ylimit[1]+245,unitstext[j])
   legend(x = 551500, y = 4118200,bty='n',legend = c(paste(">",zlimits[j,2],unitstext[j]),"","","","","","",paste(zlimits[j,1],unitstext[j])),fill = rev(colorforplots3),
          border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
   par(xpd=FALSE)
   if(j==1){
     par(xpd=TRUE)
     text(550340,4119120,"SSURGO Estimated",cex=2)
     par(xpd=FALSE)
   }
   
}
dev.off()

#generate figure 10
pdf(file = "G:/My Drive/Which DEM Paper/PaperFigs/Final paper figs/Fig10_Badistsoil.pdf",width = 11.09,height = 8.92)
par(mfrow = c(3,3),oma=c(0,0,0,0))
#bootstrap
for(j in 4:6){
  indnna = !is.na(soils[,j+1])
  soilsdata = soils[indnna,j+1]
  trainnum = length(soilsdata)/2
  for(i in 1:N){
    #define subsamples
    subsample = sample(1:length(soilsdata),trainnum,replace=FALSE)
    trainingsample = soils[subsample,j+1]
    subtest = !(1:length(soilsdata) %in% subsample)
    testsample = soilsdata[subtest]
    
    #define statistics for subsamples
    fitstats$trainx[i] = mean(trainingsample)
    fitstats$testx[i] = mean(testsample)
    fitstats$trainsig[i] = sd(trainingsample)
    fitstats$testsig[i] = sd(testsample)
    
    #select data for bootstrap
    dataforbootstrap = data.frame(bestSCA = spatial[indnna,spatialind$lnSCA[j]],bestslp = spatial[indnna,spatialind$Slope[j]])
    
    
    #for finest res D8
    df = data.frame(slplid8 = Slopes$Slope_LIDAR_2018_D8[indnna],lnscalid8 = log(SCA$SCA_LIDAR_2018_D8)[indnna])
    linearmodel = lm(trainingsample~.,data = df[subtest,])
    #predict values
    testvals = predict(linearmodel,newdata = df[subsample,])
    fitstats$finestresd8[i] = cor(testvals,testsample)
    
    #for finest res Dinf
    df = data.frame(slplidinf = Slopes$Slope_LIDAR_2018_Dinf[indnna],lnscalidinf = log(SCA$SCA_LIDAR_2018_Dinf)[indnna])
    linearmodel = lm(trainingsample~.,data = df[subtest,])
    #predict values
    testvals = predict(linearmodel,newdata = df[subsample,])
    fitstats$finestresdinf[i] = cor(testvals,testsample)
    
    #for best guess
    linearmodelbest = lm(trainingsample~.,data = dataforbootstrap[subsample,])
    fitstats$bestguesstrainxbar[i] = mean(linearmodelbest$fitted.values)
    fitstats$bestguesstrainsd[i] = sd(linearmodelbest$fitted.values)
    fitstats$bestguessinterp[i] = summary(linearmodelbest)$coefficients[1,4]
    fitstats$bestguessp[i] = pf(summary(linearmodelbest)$fstatistic[1],summary(linearmodelbest)$fstatistic[2],summary(linearmodelbest)$fstatistic[3],lower.tail = FALSE)
    fitstats$bestguesscorrtrain[i] = sqrt(summary(linearmodelbest)$r.squared)
    #predict values
    testvals = predict(linearmodelbest,newdata = dataforbootstrap[subtest,])
    fitstats$bestguesstestxbar[i] = mean(testvals)
    fitstats$bestguesstestsd[i] = sd(testvals)
    fitstats$bestguesscorrtest[i] = cor(testvals,testsample)
    
    forplottcoefvals[i,forplottcoefvalsind[j,]] = linearmodelbest$coefficients
    
  }
  
  colorforplots3 = colorforplots
  if(j==4){values(d)[which(values(d)>39)] = 39}
  if(j==6){values(f)[which(values(f)>.3)] = .3}
  
  #plot densities 
  par(mar = c(4,4,4,1))
  plot(density(fitstats$bestguesscorrtest),col = 'darkgoldenrod3',xlab = "Correlation",ylab="Density",bty='n',main=paste0("\n(A",j,") ",xlabtext[j]),xlim = c(-1,1),lwd=3,ylim = c(0,2.75))
  lines(density(fitstats$finestresd8),col = 'coral3',lwd=3)
  lines(density(fitstats$finestresdinf),col = 'darkcyan',lwd=3)
  legend('topleft',legend = c("Multivariate Regression","Finest res D8 slp&sca","Finest res Dinf slp&sca"),
         col = c("darkgoldenrod3","coral3","darkcyan"),lwd=3,bty = 'n')
  
  #plot proposed soil dist
  par(mar=c(4,2.5,4,5))
  ylimit = c(4117200,4118540); xlimit = c(549150,551500)
  slope = get(paste0(letterseq[j],"1"))
  lnsca = get(paste0(letterseq[j],"2"))
  rast = mean(forplottcoefvals[,forplottcoefvalsind[j,][1]])+mean(forplottcoefvals[,forplottcoefvalsind[j,][2]])*lnsca+mean(forplottcoefvals[,forplottcoefvalsind[j,][3]])*slope
  values(rast)[which(values(rast)<0)] = 0
  if(j==1){values(rast)[which(values(rast)>40)]=40}
  plot(rast,main = paste0("\n (B",j,") ",xlabtext[j]),bty='n',xaxt='n',yaxt='n',xlim=xlimit,ylim=ylimit,zlim=zlimits[j,],col=colorforplots3,legend=FALSE)
  if(j==4){
    par(xpd=TRUE)
    text(550340,4119040,"Multivariate TIV Regression",cex=2)
    par(xpd=FALSE)
  }
  
  #plot ssurgo soil dist
  par(mar=c(4,0,4,8.5))
  plot(get(paste0(letterseq[j])),main = paste0("\n (C",j,") ",xlabtext[j]),bty='n',xaxt='n',yaxt='n',xlim=xlimit,ylim=ylimit,zlim=zlimits[j,],col=colorforplots3,legend=FALSE)
  
  
  #if(j==1){textunitoffsetx = 700}else{textunitoffsetx = 850}
  par(xpd=TRUE)
  #text(xlimit[2]+textunitoffsetx,ylimit[2]-260,unitstext[j])
  #text(xlimit[2]+textunitoffsetx,ylimit[1]+245,unitstext[j])
  legend(x = 551500, y = 4118200,bty='n',legend = c(paste(">",zlimits[j,2],unitstext[j]),"","","","","","",paste(zlimits[j,1],unitstext[j])),fill = rev(colorforplots3),
         border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
  par(xpd=FALSE)
  if(j==4){
    par(xpd=TRUE)
    text(550340,4119120,"SSURGO Estimated",cex=2)
    par(xpd=FALSE)
  }
  
}
dev.off()
#####

#Do the soildb relate well with measured?: figure 8
PhysResults_sheetid="1HnsQjb8j9kUlBrlbPuz20czVoQvGtQwy3tmZ3vyULfc"
PhysicalProperties_Meas <- read_sheet(PhysResults_sheetid, sheet = "PhysicalResults",col_types = "?_dd_d_d__d___dd___d______")
names(PhysicalProperties_Meas) = c("SSID", "Lat", "Long", "DepthBa", "DepthBt", "OMA", "ClayA","OMBa", "ClayBa")
#Depth of A, Ba and Bt
PhysicalProperties_Meas$A_thickness = PhysicalProperties_Meas$DepthBa
PhysicalProperties_Meas$A_thickness[is.na(PhysicalProperties_Meas$A_thickness)] = PhysicalProperties_Meas$DepthBt[is.na(PhysicalProperties_Meas$A_thickness)]
PhysicalProperties_Meas$Ba_thickness = PhysicalProperties_Meas$DepthBt - PhysicalProperties_Meas$A_thickness
PhysicalProperties_Meas$Ba_thickness[PhysicalProperties_Meas$Ba_thickness==0] = NA

SP_ll=SpatialPoints(matrix(c(PhysicalProperties_Meas$Long,PhysicalProperties_Meas$Lat), 
                           ncol = 2, byrow = FALSE),proj4string =crs_ll)
extractmukeysoilsamp = over(SP_ll, mysoil)
extractmukeysoilsamp$SSID = PhysicalProperties_Meas$SSID

####Measured Properties ####
compareprops = data.frame(SSID = unlist(PhysicalProperties_Meas$SSID),
                          thickH1_meas = PhysicalProperties_Meas$A_thickness,
                          OMH1_meas = PhysicalProperties_Meas$OMA,
                          clayH1_meas = PhysicalProperties_Meas$ClayA,
                          thickH2_meas = PhysicalProperties_Meas$Ba_thickness,
                          OMH2_meas = PhysicalProperties_Meas$OMBa,
                          clayH2_meas = PhysicalProperties_Meas$ClayBa,
                          thickH1_ssurgo_r = NA,
                          OMH1_ssurgo_r = NA,
                          clayH1_ssurgo_r = NA,
                          thickH2_ssurgo_r = NA,
                          OMH2_ssurgo_r = NA,
                          clayH2_ssurgo_r = NA)

for(i in 1:length(compareprops$SSID)){
  ind = which(mu2chdomH1$mukey==extractmukeysoilsamp$mukey[i])
  compareprops$thickH1_ssurgo_r[i] = mu2chdomH1$hzdepb_r[ind]
  compareprops$OMH1_ssurgo_r[i] = mu2chdomH1$OM_r[ind]/100
  compareprops$clayH1_ssurgo_r[i] = mu2chdomH1$claytotal_r[ind]/100
  
  ind = which(mu2chdomH2$mukey==extractmukeysoilsamp$mukey[i])
  compareprops$thickH2_ssurgo_r[i] = mu2chdomH2$hzdepb_r[ind]
  compareprops$OMH2_ssurgo_r[i] = mu2chdomH2$OM_r[ind]/100
  compareprops$clayH2_ssurgo_r[i] = mu2chdomH2$claytotal_r[ind]/100
  
}

par(mfrow=c(1,3))
par(mar = c(6,6,4,4))

linm1 = lm(soils$Athick~spatial$slp1asdinf+spatial$lnsca2010dinf)
plot(compareprops$thickH1_meas,compareprops$thickH1_ssurgo_r,xlim = c(0,125),ylim = c(0,125),ylab = "Predicted horizon thickness (cm)",xlab = "Measured horizon thickness (cm)",main = "Predicted vs Measured horizon thickenss",pch=3,col = 'red',cex.main = 1.5,cex.lab = 1.5,cex.axis = 1.25)
points(compareprops$thickH1_meas,linm1$fitted.values,pch=16)
linm2 = lm(soils$Bathick~spatial$slp13asd8+spatial$lnsca2018d8)
points(compareprops$thickH2_meas,compareprops$thickH2_ssurgo_r,pch=3,col = 'blue')
points(compareprops$thickH2_meas[!is.na(compareprops$thickH2_meas)],linm2$fitted.values,pch=16,col="darkgreen")
segments(0,0,65,65,lty="dotted")
SEline = (compareprops$thickH1_meas-compareprops$thickH1_ssurgo_r)^2
SEybar = (compareprops$thickH1_meas-mean(compareprops$thickH1_meas))^2
r21 = round(1-sum(SEline)/sum(SEybar),2)
SEline = (compareprops$thickH1_meas-linm1$fitted.values)^2
SEybar = (compareprops$thickH1_meas-mean(compareprops$thickH1_meas))^2
r22 = round(1-sum(SEline)/sum(SEybar),2)

SEline = (compareprops$thickH2_meas-compareprops$thickH2_ssurgo_r)^2
SEybar = (compareprops$thickH2_meas-mean(compareprops$thickH2_meas,na.rm=TRUE))^2
r23 = round(1-sum(SEline,na.rm=TRUE)/sum(SEybar,na.rm=TRUE),0)
SEline = (compareprops$thickH2_meas[!is.na(compareprops$thickH2_meas)]-linm2$fitted.values)^2
SEybar = (compareprops$thickH2_meas-mean(compareprops$thickH2_meas,na.rm=TRUE))^2
r24 = round(1-sum(SEline,na.rm=TRUE)/sum(SEybar,na.rm=TRUE),2)
legend(51,129.5,box.lty=0,legend = c(paste("A horizon SSURGO\nR2=",r21),paste("A horizon multivariate\nregression; R2=",r22),paste("Ba horizon SSURGO\nR2=",r23),paste("Ba horizon multivariate\nregression; R2=",r24)),pch=c(3,16,3,16),col = c('red','black',"blue","darkgreen"),border='white',cex=1.25)


linm1 = lm(soils$OMA~spatial$slp1asd8+spatial$lnsca2010dinf)
plot(compareprops$OMH1_meas,compareprops$OMH1_ssurgo_r,xlim = c(0,.1),ylim = c(0,.1),ylab = "Predicted organic matter (kg/kg)",xlab = "Measured organic matter (kg/kg)",main = "Predicted vs Measured organic matter",pch=3,col = 'red',cex.main = 1.5,cex.lab = 1.5,cex.axis = 1.25)
points(compareprops$OMH1_meas,linm1$fitted.values,pch=16)
linm2 = lm(soils$OMBa~spatial$slp1asdinf+spatial$lnsca2018dinf)
points(compareprops$OMH2_meas,compareprops$OMH2_ssurgo_r,pch=3,col = 'blue')
points(compareprops$OMH2_meas[!is.na(compareprops$OMH2_meas)],linm2$fitted.values,pch=16,col="darkgreen")
segments(0,0,.085,.085,lty="dotted")
SEline = (compareprops$OMH1_meas-compareprops$OMH1_ssurgo_r)^2
SEybar = (compareprops$OMH1_meas-mean(compareprops$OMH1_meas))^2
r21 = round(1-sum(SEline)/sum(SEybar),2)
SEline = (compareprops$OMH1_meas-linm1$fitted.values)^2
SEybar = (compareprops$OMH1_meas-mean(compareprops$OMH1_meas))^2
r22 = round(1-sum(SEline)/sum(SEybar),2)

SEline = (compareprops$OMH2_meas-compareprops$OMH2_ssurgo_r)^2
SEybar = (compareprops$OMH2_meas-mean(compareprops$OMH2_meas,na.rm=TRUE))^2
r23 = round(1-sum(SEline,na.rm=TRUE)/sum(SEybar,na.rm=TRUE),2)
SEline = (compareprops$OMH2_meas[!is.na(compareprops$OMH2_meas)]-linm2$fitted.values)^2
SEybar = (compareprops$OMH2_meas-mean(compareprops$OMH2_meas,na.rm=TRUE))^2
r24 = round(1-sum(SEline,na.rm=TRUE)/sum(SEybar,na.rm=TRUE),2)
legend(-.00385,.1038,box.lty=0,legend = c(paste("A horizon SSURGO; R2=",r21),paste("A horizon multivariate regression; R2=",r22)),pch=c(3,16),col = c('red','black'),border='white',cex=1.25)
legend(-.00385,.095,bty="n",legend = c(paste("Ba horizon SSURGO; R2=",r23),paste("Ba horizon multivariate\nregression; R2=",r24)),pch=c(3,16),col = c("blue","darkgreen"),border='white',cex=1.25)

linm1 = lm(soils$ClayA~spatial$slp1asdinf+spatial$lnsca1asdinf)
plot(compareprops$clayH1_meas,compareprops$clayH1_ssurgo_r,xlim = c(0,.6),ylim = c(0,.6),ylab = "Predicted clay content (kg/kg)",xlab = "Measured clay content (kg/kg)",main = "Predicted vs Measured clay",pch=3,col = 'red',cex.main = 1.5,cex.lab = 1.5,cex.axis = 1.25)
points(compareprops$clayH1_meas,linm1$fitted.values,pch=16)
linm2 = lm(soils$ClayBa~spatial$slp2010d8+spatial$lnsca13asd8)
points(compareprops$clayH2_meas,compareprops$clayH2_ssurgo_r,pch=3,col = 'blue')
points(compareprops$clayH2_meas[!is.na(compareprops$clayH2_meas)],linm2$fitted.values,pch=16,col="darkgreen")
segments(0,0,.3,.3,lty="dotted")
SEline = (compareprops$clayH1_meas-compareprops$clayH1_ssurgo_r)^2
SEybar = (compareprops$clayH1_meas-mean(compareprops$clayH1_meas))^2
r21 = round(1-sum(SEline)/sum(SEybar),2)
SEline = (compareprops$clayH1_meas-linm1$fitted.values)^2
SEybar = (compareprops$clayH1_meas-mean(compareprops$clayH1_meas))^2
r22 = round(1-sum(SEline)/sum(SEybar),2)

SEline = (compareprops$clayH2_meas-compareprops$clayH2_ssurgo_r)^2
SEybar = (compareprops$clayH2_meas-mean(compareprops$clayH2_meas,na.rm=TRUE))^2
r23 = round(1-sum(SEline,na.rm=TRUE)/sum(SEybar,na.rm=TRUE),0)
SEline = (compareprops$clayH2_meas[!is.na(compareprops$clayH2_meas)]-linm2$fitted.values)^2
SEybar = (compareprops$clayH2_meas-mean(compareprops$clayH2_meas,na.rm=TRUE))^2
r24 = round(1-sum(SEline,na.rm=TRUE)/sum(SEybar,na.rm=TRUE),2)
legend(.25,.62,bty = 'n',legend = c(paste("A horizon SSURGO\nR2=",r21),paste("A horizon multivariate\nregression; R2=",r22),paste("Ba horizon SSURGO\nR2=",r23),paste("Ba horizon multivariate\nregression; R2=",r24)),pch=c(3,16,3,16),col = c('red','black',"blue","darkgreen"),border='white',cex=1.25)



