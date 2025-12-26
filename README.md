# Shorebirds in the Amazon Basin
![R](https://img.shields.io/badge/language-R-blue)
[![DOI](https://img.shields.io/badge/DOI-10.1093/ornithapp/duae034%2Fxxxxx-green)](https://doi.org/10.1093/ornithapp/duae034)

## Overview
This repository contains sample code and data for our manuscript, "The Amazon Basin's Rivers and Lakes Support Nearctic-breeding Shorebirds During Southward Migration", published in Ornithological Applications (https://doi.org/10.1093/ornithapp/duae034). Included is a demonstration of stopover site identification, resource selection functions, and step selection functions, using location data from Hudsonian Godwits (Limosa haemastica). how we identified stopover events and how we constructed resource and step selection functions.

## Repository contents
- 'data/'
  Tabular data and GIS vector files used in the sample scripts. 
- 'scripts/'
  R scripts used for data processing.
  
## Data
├── HUGO movement/
│   ├── HUGO preprocessed tracks.csv
│   ├── HUGO stopover locations.csv
│   ├── HUGO stopover centroids.csv
│   ├── KCH_FIT.RDS
│   ├── KCH_SIM.csv
│   ├── KCL_FIT.RDS
│   ├── KCL_SIM.csv
│   ├── KCV_FIT.RDS
│   ├── KCV_SIM.csv
│   ├── MKK_FIT.RDS
│   ├── KCV_FIT.RDS
│   └── MKK_SIM.csv
├── RSF/
│   ├── RSF_DEM_5km_AV.csv
│   ├── RSF_DEM_5km_USED.csv
│   ├── RSF_JRC_AV.csv
│   ├── RSF_JRC_USED.csv
│   ├── RSF_MB_AV.csv
│   ├── RSF_DEM_5km_USED.csv
│   └── RSF_MB_USED.csv
├── SSF/
│    ├── SSF_DEM.csv
│    ├── SSF_JRC.csv
│    ├── SSF_MB_Combined.csv
│    └── SSF_all_steps.cs
└── Amazon Basin/
    └── SNAPP_AB.shp

## Scripts
- '01 identify stopovers.R'
  How to identify stopovers in movement tracks, using Hudsonian Godwit (HUGO) tracks as an example.
- '02 RSF habitat selection.R'
  Example of a Resource Selection Function (RSF) analysis.
- '03 simulate alternative steps for SSF.R'
  Demonstration of how to generate plausible alternative steps for a Step Selection Function (SSF) analysis. 
- '04 SSF habitat selection.R'
  Example of a Resource Selection Function (RSF) analysis.
- '05 functions.R'
  Necessary functions for other scripts.

## Data sources
Data for the Resource Selection Function (RSF) and Step Selection Function (SSF) analyses are derived from the Joint Research Commission Global Surface Water occurrence layer (see Pekel et al. 2016), the annual land use and land cover (‘LULC’) data at available sites, from MapBiomas Amazonía (Colección 4.0, https://amazonia.mapbiomas.org/), and NASA's Shuttle Radar Topography Mission. All data layers were obtained from Google Earth Engine. The shapefile depicting the subbasin of the Amazon Basin (SNAPP_AB.shp) is derived from the Amazon Waters Initiative (see Venticinque et al. 2016). Hudsonian Godwit movement data is stored on Movebank. 

## Contact
For questions, contact me at linscotj@email.sc.edu.
