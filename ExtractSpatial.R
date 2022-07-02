# The purpose of this code is to extract calculated spatial properties for 
#rasters (generated via the "runDelin_taudem.R") for each location of soil sample
#please note the user will have to change the input and export data locations

#directory
basedir="G:/My Drive/Which DEM Paper/"

if (!require("pacman")) install.packages("pacman")
pacman::p_load(googlesheets4,raster)

#Options for file names
sheetname = c("Slope", "TIClass", "TIValue", "SCA")
resolution = c("_LIDAR_2010","_LIDAR_2018","_1_3arcsec","_1arcsec")
resolution2 = c("_2010LIDAR","_2018LIDAR","_1_3arcsec","_1arcsec")
methodoptions = c("_D8","_Dinf")

#Load Sheet
sheetid="1HnsQjb8j9kUlBrlbPuz20czVoQvGtQwy3tmZ3vyULfc"
PhysResults=read_sheet(sheetid,sheet = "PhysicalResults")
TIClass=read_sheet(sheetid,sheet = "TIClass") #unitless
SCA=read_sheet(sheetid,sheet = "SCA") #meters
Slope=read_sheet(sheetid,sheet = "Slope") #length/length
TIValues=read_sheet(sheetid,sheet = "TIValue") #unitless

for (i in 1:length(sheetname)){
  for (j in 1:length(resolution)){
    for (k in 1:length(methodoptions)){
      #Get extraction raster
      folder = paste(sheetname[i],resolution[j],sep = "")
      Extraction = paste(folder,methodoptions[k],sep = "") 
      dir=paste("TiffFiles/",folder,"/",Extraction,"/",sep = "")
      file=paste(sheetname[i],resolution2[j],methodoptions[k],".tif",sep = "")
      url=paste0(basedir,dir,file)
      rast=raster(url)
      proj4_utm = proj4string(rast)
      proj4_ll = "+proj=longlat"
      
      # Now we will build our proj4strings which define our “Coordinate 
      # Reference Systems” or CRS in future geographic manipulations. 
      crs_ll=CRS(proj4_ll)
      crs_utm=CRS(proj4_utm)
      
      
      #extraction
      SP_ll=SpatialPoints(matrix(c(PhysResults$Long,PhysResults$Lat), 
            ncol = 2, byrow = FALSE),proj4string =crs_ll)
      SP_utm <- spTransform(SP_ll, crs_utm)
      if(sheetname[i] == "TIClass"){
        TIClass[paste(Extraction)]=raster::extract(rast,SP_utm)
      }else if(sheetname[i] == "Slope"){
        Slope[paste(Extraction)]=raster::extract(rast,SP_utm)
      }else if (sheetname[i] == "SCA"){
        SCA[paste(Extraction)]=raster::extract(rast,SP_utm)
      }else if (sheetname[i] == "Aspect"){
        Aspect[paste(Extraction)]=raster::extract(rast,SP_utm)
      }else {
        TIValues[paste(Extraction)]=raster::extract(rast,SP_utm)
      }
    }
  }
}

#Write to Google Sheets so can use them in the regression analysis code
write_sheet(Slope,sheetid,sheet = sheetname[1])
write_sheet(TIClass,sheetid,sheet = sheetname[2])
write_sheet(TIValues,sheetid,sheet = sheetname[3])
write_sheet(SCA,sheetid,sheet = sheetname[4])
