#the purpose of this code is to generate figures 6 and 7 for occompanying paper

# The purpose of this code create correlation plots for paired regressions for to resulting data

##### Step 0: add in libraries #####
pacman::p_load(devtools)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rgdal,elevatr,raster,devtools,rgdal,readxl,plot.matrix,grid,corrplot,RColorBrewer)
colorforplots = c(brewer.pal(9,"Reds")[c(9:2)],"grey100",brewer.pal(9,"Blues")[c(2:9)])
#######

##### Step 1: Load in data ####
#please note user will need to reorganized input data
#soils data
PhysResults_sheetid="1HnsQjb8j9kUlBrlbPuz20czVoQvGtQwy3tmZ3vyULfc"
PhysicalProperties <- read_sheet(PhysResults_sheetid, sheet = "PhysicalResults")
names(PhysicalProperties) = c("SSID", "GISID", "Lat", "Long", "K_doublering", "DepthBa", "A_BaBoundary", "DepthBt", "Ba_BtBoundary", "DepthRestricting", "OMA", "GravelA", "SandA", "SiltA", "ClayA","OMBa", "GravelBa", "SandBa", "SiltBa", "ClayBa","OMBt", "GravelBt", "SandBt", "SiltBt", "ClayBt","Notes")

TIClasses <- read_sheet(PhysResults_sheetid, sheet = "TIClass")
Slopes <- read_sheet(PhysResults_sheetid, sheet = "Slope")
SCA <- read_sheet(PhysResults_sheetid, sheet = "SCA")
TIValues <- read_sheet(PhysResults_sheetid, sheet = "TIValue")
########

##### Step 2: clean soils data #####
soils = data.frame(SSID = unlist(PhysicalProperties$SSID),
                   Athick = NA,OMA = PhysicalProperties$OMA,ClayA = PhysicalProperties$ClayA,
                   Bathick = NA,OMBa = PhysicalProperties$OMBa, ClayBa = PhysicalProperties$ClayBa)
soils$Athick = as.numeric(PhysicalProperties$DepthBa)
soils$Athick[is.na(soils$Athick)] = PhysicalProperties$DepthBt[is.na(soils$Athick)]
soils$Bathick = PhysicalProperties$DepthBt - as.numeric(unlist(PhysicalProperties$DepthBa))
##############

###### Step 3: calculate and plot corplot ####
#Creating a cross correalation matrix for phys properties
cor4plot = data.frame(A_thick = matrix(0,32,1), OMA = 0, ClayA = 0,
                      Ba_thick = 0, OMBa = 0, ClayBa = 0)

rownames(cor4plot) = c("Slp 2018LIDinf","ln(SCA) 2018LIDinf","TIV 2018LIDinf","TIC 2018LIDinf","Slp 2018LID8","ln(SCA) 2018LID8","TIV 2018LID8","TIC 2018LID8",
                       "Slp 2010LIDinf","ln(SCA) 2010LIDinf","TIV 2010LIDinf","TIC 2010LIDinf","Slp 2010LID8","ln(SCA) 2010LID8","TIV 2010LID8","TIC 2010LID8",
                       "Slp 1/3asDinf","ln(SCA) 1/3asDinf","TIV 1/3asDinf","TIC 1/3asDinf","Slp 1/3asD8","ln(SCA) 1/3asD8","TIV 1/3asD8","TIC 1/3asD8",
                       "Slp 1asDinf","ln(SCA) 1asDinf","TIV 1asDinf","TIC 1asDinf","Slp 1asD8","ln(SCA) 1asD8","TIV 1asD8","TIC 1asD8")
spatial = cbind(Slopes$Slope_LIDAR_2018_Dinf,log(SCA$SCA_LIDAR_2018_Dinf),TIValues$TIValue_LIDAR_2018_Dinf,TIClasses$TIClass_LIDAR_2018_Dinf,
                Slopes$Slope_LIDAR_2018_D8,log(SCA$SCA_LIDAR_2018_D8),TIValues$TIValue_LIDAR_2018_D8,TIClasses$TIClass_LIDAR_2018_D8,
                
                Slopes$Slope_LIDAR_2010_Dinf,log(SCA$SCA_LIDAR_2010_Dinf),TIValues$TIValue_LIDAR_2010_Dinf,TIClasses$TIClass_LIDAR_2010_Dinf,
                Slopes$Slope_LIDAR_2010_D8,log(SCA$SCA_LIDAR_2010_D8),TIValues$TIValue_LIDAR_2010_D8,TIClasses$TIClass_LIDAR_2010_D8,
                
                Slopes$Slope_1_3arcsec_Dinf,log(SCA$SCA_1_3arcsec_Dinf),TIValues$TIValue_1_3arcsec_Dinf,TIClasses$TIClass_1_3arcsec_Dinf,
                Slopes$Slope_1_3arcsec_D8,log(SCA$SCA_1_3arcsec_D8),TIValues$TIValue_1_3arcsec_D8,TIClasses$TIClass_1_3arcsec_D8,
                
                Slopes$Slope_1arcsec_Dinf,log(SCA$SCA_1arcsec_Dinf),TIValues$TIValue_1arcsec_Dinf,TIClasses$TIClass_1arcsec_Dinf,
                Slopes$Slope_1arcsec_D8,log(SCA$SCA_1arcsec_D8),TIValues$TIValue_1arcsec_D8,TIClasses$TIClass_1arcsec_D8)

physical = soils[,2:7]

#correlation for spatial data
for(i in 1:32){
  #i is rows of correlation matrix (dems)
  for(k in 1:length(names(physical))){
    cor4plot[i,k] = cor(spatial[,i],physical[,k],use = "complete.obs")
  }
}

par(xpd=FALSE)
rownames(cor4plot) = c("(1) Slope","(2) ln(SCA)","(3) TIV","(4) TIC","(5) Slope","(6) ln(SCA)","(7) TIV","(8) TIC","(9) Slope","(10) ln(SCA)","(11) TIV","(12) TIC","(13) Slope","(14) ln(SCA)","(15) TIV","(16) TIC","(17) Slope","(18) ln(SCA)","(19) TIV","(20) TIC","(21) Slope","(22) ln(SCA)","(23) TIV","(24) TIC","(25) Slope","(26) ln(SCA)","(27) TIV","(28) TIC","(29) Slope","(30) ln(SCA)","(31) TIV","(32) TIC")
corrplot(t(as.matrix(cor4plot)),tl.col = "black",mar=c(0,0,3,0),col = colorforplots,tl.cex=.75,main = 'Measured Soil Correlations',cl.length = 5,number.cex=.5,number.digits = 2)
x = 8;y=3.5
text(c(x-y,x*2-y,x*3-y,x*4-y),14,c("2018 LIDAR\n0.76m","2010 LIDAR\n1.5m","1/3as USGS\n9m","1as USGS\n28m"),cex=.75)
y=5.5
text(c(x-y,x*2-y,x*3-y,x*4-y),12,c("Dinf"),cex=.75)
y=1.5
text(c(x-y,x*2-y,x*3-y,x*4-y),12,c("D8"),cex=.75)
segheight = 17
segments(.54,.5,.54,segheight); segments(4.54,.5,4.54,segheight-4,lty = 4); segments(8.54,.5,8.54,segheight); segments(12.54,.5,12.54,segheight-4,lty = 4); segments(16.54,.5,16.54,segheight); 
segments(20.54,.5,20.54,segheight-4,lty = 4); segments(24.54,.5,24.54,segheight); segments(28.54,.5,28.54,segheight-4,lty = 4); segments(32.54,.5,32.54,segheight); 
segments(.5,3.5,32.5,3.5)
#####


##### Step 4: Multiple (2) regression #####
#models for regressions
#please note that this portion of analysis was not automated by code, AICs and correlations etc were calculated manually via R
step = matrix(NA, 6,1)
rownames(step) = colnames(cor4plot)
colnames(step) = c("Regression")
coefval = matrix(NA,6,8)
rownames(coefval) = rownames(step)
colnames(coefval) = c("2018slp","2010slp","13asslp","1asslp","2018sca","2010sca","13assca","1assca")

#a horizon
model = lm(soils$Athick~.,data = slpsca[,c(7,11)]) 
coefval[1,4] = round(model$coefficients[2],2)
coefval[1,6] = round(model$coefficients[3],2)
step[1,1] = summary(model)$r.squared

model = lm(soils$OMA~.,data = slpsca[,c(8,11)])
coefval[2,4] = round(model$coefficients[2],2)
coefval[2,6] = round(model$coefficients[3],3)
step[2,1] = summary(model)$r.squared

model = lm(soils$ClayA~.,data = slpsca[,c(7,15)]) #insignificant (0.166)
coefval[3,4] = round(model$coefficients[2],3)
coefval[3,8] = round(model$coefficients[3],3)
step[3,1] = summary(model)$r.squared

model = lm(soils$Bathick~.,data=slpsca[,c(6,10)])
coefval[4,3] = round(model$coefficients[2],2)
coefval[4,5] = round(model$coefficients[3],2)
step[4,1] = summary(model)$r.squared

model = lm(soils$OMBa~.,slpsca[,c(7,9)])
coefval[5,4] = round(model$coefficients[2],2)
coefval[5,5] = round(model$coefficients[3],3)
step[5,1] = summary(model)$r.squared

model = lm(soils$ClayBa~.,data=slpsca[,c(4,14)]) 
coefval[6,2] = round(model$coefficients[2],2)
coefval[6,7] = round(model$coefficients[3],2)
step[6,1] = summary(model)$r.squared

combinedtext = matrix("",6,8)
for(i in 1:6){
  for(j in 1:8){
    if(!is.na(coefval[i,j])){combinedtext[i,j] = coefval[i,j]}
  }
}

legendmatrix = matrix(NA,2,2)
rownames(legendmatrix) = c("Dinf","D8")
colnames(legendmatrix) = c("p < 0.05","p > 0.05")
legendmatrix[] = cbind(c(6,3),c(4,1))

TFmatrixd8dinf = rbind(c(0,0,0,6, 0,6,0,0), #athick
                       c(0,0,0,3, 0,6,0,0), #oma
                       c(0,0,0,4, 0,0,0,4), #claya
                       c(0,0,3,0, 3,0,0,0), #bathick
                       c(0,0,0,6, 6,0,0,0), #omba
                       c(0,3,0,0, 0,0,3,0)) #clayba

layout(t(matrix(c(1:4))),widths = c(.95,2,2,1),heights = 1)
colnames(TFmatrixd8dinf) = c("Slope LIDAR 2018","Slope LIDAR 2010","Slope 1/3as","Slope 1as","SCA LIDAR 2018","SCA LIDAR 2010","SCA 1/3as","SCA 1as")
colnames(step)=""
corrplot(step,tl.col = "black",mar=c(6,0,0,0),col = colorforplots,method = 'number',tl.cex=1.5,cl.pos='n',number.digits = 2,number.cex = 3)
segments(.5,3.5,1.5,3.5,lwd=2)

colorTF = c('white',brewer.pal(3,"OrRd")[2:3],brewer.pal(3,"BuGn")[2:3])
par(mar = c(16,0.25,3.75,0.25))
plot((TFmatrixd8dinf[,1:4]),col = colorTF,las=2,main = "",xlab="",ylab="",key=NULL,cex.axis=1.45,axis.row = NULL,axis.col=list(side=1))
text(c(1,2,3,4),6,combinedtext[1,1:4],cex=1.5)
text(c(1,2,3,4),5,combinedtext[2,1:4],cex=1.5)
text(c(1,2,3,4),4,combinedtext[3,1:4],cex=1.5)
text(c(1,2,3,4),3,combinedtext[4,1:4],cex=1.5)
text(c(1,2,3,4),2,combinedtext[5,1:4],cex=1.5)
text(c(1,2,3,4),1,combinedtext[6,1:4],cex=1.5)
segments(0,3.5,30,3.5,lwd=2)
#par(mar = c(12,0.25,1,6))

plot((TFmatrixd8dinf[,5:8]),col = colorTF,las=2,main = "",xlab="",ylab="",key=NULL,cex.axis=1.45,axis.row = NULL,axis.col=list(side=1))
text(c(1,2,3,4),6,combinedtext[1,5:8],cex=1.5)
text(c(1,2,3,4),5,combinedtext[2,5:8],cex=1.5)
text(c(1,2,3,4),4,combinedtext[3,5:8],cex=1.5)
text(c(1,2,3,4),3,combinedtext[4,5:8],cex=1.5)
text(c(1,2,3,4),2,combinedtext[5,5:8],cex=1.5)
text(c(1,2,3,4),1,combinedtext[6,5:8],cex=1.5)
segments(0,3.5,30,3.5,lwd=2)
par(mar = c(27,3,10,3))
plot(legendmatrix,col = colorTF[2:5],las=2,main = "Legend",xlab="",ylab="",key = NULL,cex.axis=1.25)
