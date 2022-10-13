Evaluating Digital Elevation Models and Topographic Indices for Geomorphic Landscape Representaton [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6791484.svg)](https://doi.org/10.5281/zenodo.6791484)
=================

This repository contains all the data, spatial files, and code used to execute the analysis in Evaluating Digital Elevation Models and Topographic Indices for Geomorphic Landscape Representation.

If you have any questions regarding this publication please contact Elyce Buell (<enb46@cornell.edu> or <ebuell@vt.edu>).

## Links
See the following links for more information on  `R` and `RStudio` download and installation:

- An introduction to `R`: <https://cran.r-project.org/doc/manuals/r-release/R-intro.pdf>
- `R` download: <https://www.r-project.org/>
- `RStudio` download: <https://www.rstudio.com/>

There is also a cloud-based `RStudio` sever at the following location:

- Cloud-based `RStudio` server: <https://rstudio.cloud/>

See the following likes for information regarding `DEM` data

- USGS DEM download <https://apps.nationalmap.gov/downloader/#/>
- For information regarding LiDAR data please contact Elyce Buell (enb46@cornell.edu or ebuell@vt.edu) or Johnathan Resop (resop@vt.edu) 

## Excels (excel and sheet names):

### PhysicalProperties_datainpaper.xlsx
Summarizes physical properties of soil cores taken at 36 locations in southwest virginia
#### PhysicalResults	
Final texture, organic matter, and horizon depth information found here. All categories included in this excel are used in the paper analysis. 
Details on direct measurements and QAQC can be found in KSatMeas, OrgMatter, and SoilTexture.
#### WSView 		      
Locations of soil cores (numbers correspond to WSView ID)

PhysicalProperties_alldata
	PhysicalResults	- 	Summarizes physical properties of soil cores taken at 36 locations in southwest virginia
				              Final texture, organic matter, and horizon depth information found here. Not all categories included in this excel are used in the paper        
                      analysis. Details on direct measurements and QAQC can be found in KSatMeas, OrgMatter, and SoilTexture.
	WSView 		      -	  Locations of soil cores (numbers correspond to WSView ID)

OrganicMatter
	OM	            -	  Summarizes OM for all horizons. Loss on ignition tests are conducted on each horizon for all soil cores. Soil samples ranged from 50g to 250g. 
                      Samples were dried for 2 days in a Fisher Scientific Isotherm oven set to 105 degC and are subsequently weighed. These samples are incinerated 
                      and before and after weights are compared to yield percent organic matter by weight. The samples are incinerated at 425 degC. The length of 
                      incineration was dictated by the incinerator used. For the Barnstead International Model 1730 12 sample incinerator, samples are incinerated for 
                      a minimum of 6 hours. For the Thermolyne Furnace Type 10500 1 sample incinerator, samples are incinerated for a minimum of 4 hours. A 4 hour 
                      incineration in the Barnstead International Model 1730 was shown to have magnitude errors of up to 1% organic matter and thus 2 additional hours 
                      were added on to the incineration time.
	WSView 		      - 	Locations of soil cores (numbers correspond to WSView ID)
	OrgMatter	      -	  Direct weight measurements and calculations for organic matter
	QAQC_NewSoil	  -	  New soil subsamples from the same core/horizon are analyzed for organic matter
	QAQC_NewSoilAnalysis-	Analysis to see if the errors seen in these QAQC runs is tied to horizon thickness or the amount of time that the soil has been left unfrozen 
                      in between runs
	QAQC_SampleRerun-	  Already incinerated samples were rerun in order to check that the current methodology was suffcient accurate within reason.

SoilTexture
	TextureFinal	  -	  Final texture (sand and clay) values decided by R code. For more information on texture descision making for this calculation is in the R code
	WSView 		      -	  Locations of soil cores (numbers correspond to WSView ID)
	SoilGravel	    -	  Measured gravel content using a 2mm sieve
	TextureSummary	-	  Relevant measurement values for calculating soil texture
	QAQCTextureSummary-	Summarizes QAQC and orgininal runs; this includes two types of QAQC: reruns (i.e. shake the soil column again with same soil) and new soil from 
                      soil column

HydrometerMeas		Direct hydrometer measurements following the method described in ASTM D4222. In order to directly calculate the hydrometer corrects (see ASTM D4222), 
                  several corrections need to be measured. This ultimately did not pan out and texture was calculated using envalysis package.
	TextureSummary	-	  Relevant measurement values for calculating soil texture
	QAQCTextureSummary-	Relevant measurement values for calculating QAQC soil texture
	WSView 		      -	  Locations of soil cores (numbers correspond to WSView ID)	
	SoilTextureA	  - 	All measurements for hydrometer trials for A horizon
	SoilTextureBa	  -	  "	"	"	"	"	"  Ba horizon
	SoilTextureBt 	-	  "	"	"	"	"	"  Bt horizon
	TempCorrMeas	  -	  Measured hydrometer correction values in distilled water at a variety of temperatures. A very poor relationship between temperature and 
                      hydrometer correction is found implying the poor reliability of the thermometers being used in the hydrometer method.
	CorrectionSummary-	Direct measurement corrections.

KSatMeas
	KSAT		        - 	Ksat values and lats and longs
	WSView 		      -	  Locations of soil cores (numbers correspond to WSView ID)
	SatHydCondTestx	-	  Sheets with direct field readings for the Eijkelkamp Soil & Water double-ring infiltrometer test following (Bouwer, 1986;  DOI: 
                      10.2136/sssabookser5.1.2ed.c32)

DEMderivedSpatialAtt	All numbers in this document are derived using TauDEM: available at: https://github.com/dtarb/TauDEM. Both D8 and Dinf flow direction algorithms 
                      (discussed in	Tarboton DG. 1997; DOI: 10.1029/96WR03137) are used as well as the four DEMs (LIDAR 2018, LIDAR 2010, 1/3as, and 1as).
	SCA		          -	  Speicfic catchment area (expressed in m)
	Slope		        -	  Slope (expressed in length/length)
	TIValue		      - 	TIV = ln(SCA/slope)
	TIClass		      -	  TIC bins the TIVs into 10 equally sized categories ranging from 1-10 (where the highest TIVs correspond to TIC = 10 etc.)
	WSView 		      -	  Locations of soil cores (numbers correspond to WSView ID)
	Aspect		      -	  Aspects range from 1-8. Refer to TauDEM documentation to understand these values.

MeasSpatProp
	SpatialProp	    -	  Measured spatial properties (slope and aspect); a phone app clinometer 
                      (https://play.google.com/store/apps/details?id=com.plaincode.clinometer&hl=en_US&gl=US) and phone compass are used to measure these values. The 
                      measurement was taken with the assistance of a 4-wheel vehicle; the vehicle was oriented perpendicularly to the slope and the slope that the 
                      utility vehicle is experiancing is recorded. The distance between the two front wheels of the utility vehicle is 1.2m. The aspect is measured by 
                      pointing the compass perpendicularly and downslope of how the utility vehicle is oriented for the slope measurement.
	WSView 		      -	  Locations of soil cores (numbers correspond to WSView ID)
---------------------------------------------------------------------------------------


R Code
---------------------------------------------------------------------------------------
Texture		        -	  Calculates texture of samples using envalysis package
runDelin_taudem	  -	  Delineates watershed and export tiff files for slope, catchment area, aspect, topographic index value, and topographic index class
                      Outlet for the watershed studied (in lat long): 37.204526, -80.445175
ExtractSpatial	  -	  Extracts spatial data (rasters resulting from runDelin_taudem) from all soil sampling locations
Fig2to5		        -	  Figure generation for figures 2-5 and figures S2-S5
Fig6and7	        -	  Figure generation for figures 6 and 7
Fig8to10	        -	  Figure generation for figures 8-10
Figsupp		        - 	Supplemental figure generation
--------------------------------------------------------------------------------------

Spatial data
--------------------------------------------------------------------------------------
Construction_exclusion      -   A box that surround highway 460 throughout the entire watershed so analysis can exclude this area
Monitoring_point            -	  Outlet location
pathsunderoverpass          -	  A path that is created to mimic storm drainage culverts that run below highway 460. Only used for LIDAR analysis; USGS DEM did not have 
                                any highway interfereance.
SoilSampleLocations         -	  Soil sampling locations
1as_arcproj_bilin           -	  1 as DEM bilinear projection executed in ArcMap
1_3as_arcproj_bilin         -	  1/3 as DEM bilinear projection executed in ArcMap
2010LIDAR_arcproj_bilin     -   2010LIDAR DEM bilinear projection executed in ArcMap
2018LIDAR_arcproj_bilin     -   2018LIDAR DEM bilinear projection executed in ArcMap; dataset too big to upload to github, please find via google drive here: 
                                https://drive.google.com/file/d/1fc_7hbiMM8hfNFi38inOW2mjg_jeq-K6/view?usp=sharing
--------------------------------------------------------------------------------------
