# Shorebirds in the Amazon Basin
![R](https://img.shields.io/badge/language-R-blue)
[![DOI](https://img.shields.io/badge/DOI-10.1093/ornithapp/duae034%2Fxxxxx-green)](https://doi.org/10.1093/ornithapp/duae034)

## Overview
This repository contains sample code and data for our manuscript, "The Amazon Basin's Rivers and Lakes Support Nearctic-breeding Shorebirds During Southward Migration", published in Ornithological Applications (https://doi.org/10.1093/ornithapp/duae034). Included is a demonstration of stopover site identification, resource selection functions, and step selection functions, using location data from Hudsonian Godwits (Limosa haemastica). 

<br/>

## Repository contents
- 'data/'
  Tabular data and GIS vector files used in the sample scripts. 
- 'scripts/'
  R scripts used for data processing.

<br/>

## Data
Data for the Resource Selection Function (RSF) and Step Selection Function (SSF) analyses are derived from the Joint Research Commission Global Surface Water occurrence layer, the annual land use and land cover (‘LULC’) data at available sites, from MapBiomas Amazonía (Colección 4.0), and NASA's Shuttle Radar Topography Mission. All data layers were obtained from Google Earth Engine. The shapefile depicting the subbasin of the Amazon Basin (SNAPP_AB.shp) is from the Amazon Waters Initiative. Hudsonian Godwit movement data is stored on Movebank.

<br/>

## Scripts
- '01 identify stopovers.R' identifies stopovers in movement tracks, using example tracks from Hudsonian Godwits (HUGO)   
- '02 RSF habitat selection.R' demonstrates a resource selection function (RSF) analysis for shorebird stopovers  
- '03 simulate alternative steps for SSF.R' demonstrates how to generate plausible alternative steps for a step selection function (SSF) analysis.   
- '04 SSF habitat selection.R' demonstrates an SSF analysis  
- '05 functions.R' contains necessary functions for the scripts above

<br/>
  
## Contact
For questions, contact me at linscotj@email.sc.edu.
