#The purpose of this code is to find the texture given set of hydrometer
# measurements using the evalysis package
#hydrometer method in accordance with ASTM D4222

#Summary of hydrometer corrections
#Density correction calculations from controlled experiment
#ignore miniscus correction because it only complicates things for this calculation (both readings include miniscus so the miniscus correction cancels out)
#correct temperture correction readings: temp corr = raw reading - miniscus
#Truth = meas - (dencorr+tempcorr+minisuscorr)
#dencorr ~ +5.7; tempcorr ~ -1.5 to -3; miniscuscorr ~ -0.5

#Texture method:
#Use texture function to find clay with the following corrections
#Temp correction - report maximum and minimum (for uncertainty) 
#   Temperature correction is done this way because raw measured corrections did not relate to temperature
#                 It is likely that the thermometer was faulty
#do analysis on the mean of the max and min correction of them 

library('googlesheets4')
library('envalysis')

#import data
#please note user will need to change this in order to run this code
#please import texture and QAQCtexture found in the SoilTexture excel under the sheet names TextureSummary and QAQCTexctureSummary
PhysResults_sheetid="1HnsQjb8j9kUlBrlbPuz20czVoQvGtQwy3tmZ3vyULfc"
Texture <- read_sheet(PhysResults_sheetid, sheet = "TextureSummary")
QAQCTexture <- read_sheet(PhysResults_sheetid, sheet = "QAQCTextureSummary")

#Find average density correction
dencorrmeas = c(5.5, 6.5, 6.5, 5.5, 5.5, 5.5, 5.5, 5.5, 6.5, 5.5, 5.5, 6) #direct density correction measurements
dencorr = mean(dencorrmeas)

#set up dataframe
ClaySand = data.frame(SSID = Texture$SSID, claytmin = NA, claytmax = NA,
                      sandtmin = NA, sandtmax = NA, sand = NA, clay = NA)

#Texture analysis with R envanalysis library function
#Temperature correction added to blank for bounds of correction
miniscus = -.5
blank = matrix(dencorr,7,1)
temp = matrix(20,7,1)
for (i in 1:length(ClaySand$SSID)){
  time = c(unlist(as.numeric(Texture[i,which(grepl("Time",names(Texture)))])))
  reading = c(unlist(as.numeric(Texture[i,which(grepl("Hyd",names(Texture)))])))-(-1.5+miniscus)
  weight = matrix(as.numeric(Texture[i,which(grepl("Weight",names(Texture)))]),7,1)
  x = texture(reading ~ blank + time + temp, con = weight, hydrometer = "152H")
  ClaySand$claytmin[i] = round(x$din[1,1],3)
  ClaySand$sandtmin[i] = round(x$din[1,3],3)
  reading = c(unlist(as.numeric(Texture[i,which(grepl("Hyd",names(Texture)))])))-(-3+miniscus)
  x = texture(reading ~ blank + time + temp, con = weight, hydrometer = "152H")
  ClaySand$claytmax[i] = round(x$din[1,1],3)
  ClaySand$sandtmax[i] = round(x$din[1,3],3)
}
ClaySand$sand = rowMeans(ClaySand[,4:5])
ClaySand$clay = rowMeans(ClaySand[,2:3])
hist(ClaySand$claytmax-ClaySand$claytmin,30)

#please note user may also need to change this
export = data.frame(SSID = ClaySand$SSID, Clay = ClaySand$clay, Sand = ClaySand$sand)
write_sheet(export,PhysResults_sheetid,sheet = "TextureFinal")


#set up dataframe
QAQCClaySand = data.frame(SSID = QAQCTexture$SSID, claytmin = NA, claytmax = NA,
                      sandtmin = NA, sandtmax = NA, sand = NA, clay = NA)

#Texture analysis with R envanalysis library function
#Temperature correction added to blank for bounds of correction
for (i in 1:length(QAQCTexture$SSID)){
  time = c(unlist(as.numeric(QAQCTexture[i,which(grepl("Time",names(Texture)))])))
  reading = c(unlist(as.numeric(QAQCTexture[i,which(grepl("Hyd",names(Texture)))])))-(-1.5+miniscus)
  weight = matrix(as.numeric(QAQCTexture[i,which(grepl("Weight",names(Texture)))]),7,1)
  x = texture(reading ~ blank + time + temp, con = weight, hydrometer = "152H")
  QAQCClaySand$claytmin[i] = round(x$din[1,1],3)
  QAQCClaySand$sandtmin[i] = round(x$din[1,3],3)
  reading = c(unlist(as.numeric(QAQCTexture[i,which(grepl("Hyd",names(Texture)))])))-(-3+miniscus)
  x = texture(reading ~ blank + time + temp, con = weight, hydrometer = "152H")
  QAQCClaySand$claytmax[i] = round(x$din[1,1],3)
  QAQCClaySand$sandtmax[i] = round(x$din[1,3],3)
}
QAQCClaySand$sand = rowMeans(QAQCClaySand[,4:5])
QAQCClaySand$clay = rowMeans(QAQCClaySand[,2:3])

#please note user may also need to change this
export = data.frame(SSID = QAQCClaySand$SSID, Clay = QAQCClaySand$clay, Sand = QAQCClaySand$sand)
write_sheet(export,PhysResults_sheetid,sheet = "QAQCTextureFinal")
