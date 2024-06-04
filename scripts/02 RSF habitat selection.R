#############################################X
#-- Large-Scale (Lev II) Habitat Selection --X
#############################################X
#------------ Prepare  & Model --------------X
#############################################X

library(sf)
library(mapview)
library(ggplot2)
library(pROC)
library(geojsonsf)
library(MuMIn)
library(lubridate)
library(spdep)
library(broom.mixed)
library(lmtest)
library(sandwich)
source('scripts/05 functions.R')


# Load Geographic Data  --------------------------------------------------------

# SNAPP Amazon Basin shapefiles (Venticinque et al. 2016)
# see: https://snappartnership.net/teams/amazon-waters/
ab <- st_geometry(st_read('data/SNAPP_AB.shp'))
basins <- st_read('data/SNAPP_Subbasins.shp')

# Load stopover centroids 
centroids <- read.csv('data/HUGO stopover centroids.csv')

# Create 'Available' Locations for Comparison ----------------------------------

# Set projection
proj <- 'PROJCS["Custom_Cylindrical_Equal_Area",
               GEOGCS["GCS_WGS_1984",
                      DATUM["D_WGS_1984",
                            SPHEROID["WGS_1984",6378137.0,298.257223563]],
                      PRIMEM["Greenwich",0.0],
                      UNIT["Degree",0.0174532925199433]],
               PROJECTION["Cylindrical_Equal_Area"],
               PARAMETER["False_Easting",0.0],
               PARAMETER["False_Northing",0.0],
               PARAMETER["Central_Meridian",-59.5019531],
               PARAMETER["Standard_Parallel_1",11.4589615],
               UNIT["Meter",1.0]]'

# Negative buffer for the Amazon Basin
ab_buff <- ab %>%
  st_transform(crs = st_crs(proj)) %>%
  st_buffer(dist = -20000) %>%
  st_transform(crs = 4326)

# Inspect
mapview(list(ab,ab_buff),
        col.regions=list("red","blue"),
        col=list("red","blue"))

# Generate random locs
nlocs = 1000 # 700000 used in analysis
random <- st_sample(ab_buff, size = nlocs, type = "random", exact = TRUE)

# Convert to data frame
random_df <- as.data.frame(st_coordinates(random)) %>%
  rename(lon = 'X', lat = 'Y') %>%
  mutate(use = 0, id = 1:n())

# Convert to sf
random_sf <- st_as_sf(random_df, coords = c('lon', 'lat'), crs = 4326)

# Append sub-basin information
random_sf <- st_join(random_sf, basins, join = st_within) 

get_basins <- random_sf %>%
  st_drop_geometry() %>% 
  dplyr::select(id, BL3)

# Final data frame
random_final <- left_join(random_df, get_basins) %>%
  dplyr::select(id, use, lon, lat, BL3) %>%
  group_by(id) %>%
  filter(n() == 1) %>%
  ungroup() 

# Inspect
head(random_final)


# Identify visited subbasins ---------------------------------------------------

# Load tracking data
sp_df <- read.csv('data/HUGO processed tracks.csv') %>% 
  dplyr::select(yearly_id, timestamp, location.long, location.lat)

# Generate straight-line movement paths
sp_sf_polylines <- sp_df %>%
    mutate(timestamp = ymd_hms(timestamp)) %>% 
    st_as_sf(coords = c('location.long', 'location.lat'), crs = 4326) %>%
    group_by(yearly_id) %>% 
    dplyr::summarize(do_union=FALSE) %>% 
    st_cast("LINESTRING")

# Identify basins crossed  
sp_basins <- st_join(sp_sf_polylines, basins, join = st_intersects) %>%
    st_drop_geometry() %>%
    distinct(yearly_id, BL3) %>%
    na.omit()


# Extract geographic data -----------------------------------------------------
# done in GEE

# Habitat data: Mapbiomas (downloaded from GEE; RAISG 2023)
# see: https://amazonia.mapbiomas.org/ 
mb_used_sp <- read.csv('data/RSF_MB_USED.csv') 
mb_av_sp <- read.csv('data/RSF_MB_AV.csv') 

# Habitat data: JRC global surface water mapping (downloaded from GEE; Pekel et al. 2016) 
# see: https://global-surface-water.appspot.com/
jrc_used_sp <- read.csv('data/RSF_JRC_USED.csv') 
jrc_av_sp <- read.csv('data/RSF_JRC_AV.csv') 

# Habitat data: elevation 
# see: NASA Shuttle Radar Topography Mission (SRTM; Farr et al. 2007)
dem_used <- read.csv('data/RSF_DEM_5km_USED.csv') 
dem_av <- read.csv('data/RSF_DEM_5km_AV.csv') 


# Prepare data frame -----------------------------------------------------------

# Vector of sacles
scales = c(2, 3, 4, 5, 10, 15, 20)

# Loop to create data frames for each scale
for(i in 1:length(scales)){
  
  # used DF -- filter to scale of selection 
  mb_used <- mb_used_sp %>% 
    filter(scale == scales[i]) %>% 
    dplyr::select(-scale) 
  jrc_used <- jrc_used_sp %>% 
    filter(scale == scales[i]) %>% 
    dplyr::select(-scale)
  used_vars <- mb_used %>% 
    left_join(jrc_used) %>% 
    left_join(dem_used)
  
  # used DF -- filter to species
  used_vars_sp <- used_vars %>% 
    left_join(centroids) %>% 
    mutate(yearly_id = paste0(id, year(start)))
  
  # available DF -- filter to scale of selection   
  mb_av <- mb_av_sp %>% 
    filter(scale == scales[i]) %>% 
    dplyr::select(-scale)
  jrc_av <- jrc_av_sp %>% 
    filter(scale == scales[i]) %>% 
    dplyr::select(-scale)
  av_vars <- left_join(mb_av, jrc_av) %>% 
    left_join(dem_av) %>% 
    left_join(av)

  # available DF -- filter to basins visited by individual
  indivs <- used_vars_sp %>% 
    pull(yearly_id) %>% 
    unique()
  
  if (scales[i] == 2){ 
    av_vec = c()
    id_vec = c()
    
    for(i in 1:length(indivs)){
      # basins visited
      b <- sp_basins %>% 
        filter(yearly_id == indivs[i])
      # number of used 
      n <- used_vars_sp %>% 
        filter(yearly_id == indivs[i]) 
      # size of used:available ratio (here, 1:500)
      sz <- 500*nrow(n)  
      v <- rep(indivs[i], sz)
      # sample available locations
      a <- av_vars %>% 
        filter(BL3 %in% b$BL3) %>% 
        slice_sample(n = sz) %>% 
        pull(loc_id)

    }
    
    av_vars_sp <- av_vars %>% 
      filter(loc_id %in% a) %>% 
      dplyr::select(-loc_id) 
    av_vars_sp$yearly_id <- v
    habitat_df_2km <- bind_rows(used_vars_sp, av_vars_sp) %>% 
      mutate(ipp_w = ifelse(used == 0, 10000, 1))
  
    } else {
    av_vars_sp <- av_vars %>% 
      filter(loc_id %in% a) %>% 
      dplyr::select(-loc_id)
    av_vars_sp$yearly_id <- v
    w <- bind_rows(used_vars_sp, av_vars_sp) %>% 
      mutate(ipp_w = ifelse(used == 0, 10000, 1))
    assign(paste0('habitat_df_', scales[i], 'km'), w)
  }
}



# Identify characteristic Scale ------------------------------------------------

# Scale-varying variables 
sv_vars <- c('forcat_prop', 'forest_prop', 'savanna_prop', 
             'grassland_prop', 'ag_prop', 'pasture_prop', 
             'wetland_prop', 'farming_prop', 'RLO_prop', 
             'jrc_permanent_prop', 'jrc_semi_prop', 'jrc_temporary_prop', 
             'jrc_dry_prop', 'elevation')

# Create a blank DF
df <- data.frame(matrix(ncol = 9, nrow = 1))
colnames(df) <- c('species', 'covariate', 'AICc_2km', 
                  'AICc_3km', 'AICc_4km', 'AICc_5km', 
                  'AICc_10km', 'AICc_15km', 'AICc_20km')

# Loop to fill data frame
for(i in 1:length(sv_vars)){
  
  df[i,1] <- "HUGO"
  df[i,2] <- sv_vars[i]
  
  rsf_2km <- glm(paste('used ~', sv_vars[i], sep = " "), 
                 family = binomial(link=logit), 
                 data = habitat_df_2km, 
                 weight = ipp_w)
  rsf_3km <- glm(paste('used ~', sv_vars[i], sep = " "), 
                 family = binomial(link=logit), 
                 data = habitat_df_3km, 
                 weight = ipp_w)
  rsf_4km <- glm(paste('used ~', sv_vars[i], sep = " "), 
                 family = binomial(link=logit), 
                 data = habitat_df_4km, weight = ipp_w)
  rsf_5km <- glm(paste('used ~', sv_vars[i], sep = " "), 
                 family = binomial(link=logit), 
                 data = habitat_df_5km, weight = ipp_w)
  rsf_10km <- glm(paste('used ~', sv_vars[i], sep = " "), 
                  family = binomial(link=logit), 
                  data = habitat_df_10km, weight = ipp_w)
  rsf_15km <- glm(paste('used ~', sv_vars[i], sep = " "), 
                  family = binomial(link=logit), 
                  data = habitat_df_15km, weight = ipp_w)
  rsf_20km <- glm(paste('used ~', sv_vars[i], sep = " "), 
                  family = binomial(link=logit), 
                  data = habitat_df_20km, weight = ipp_w)
  

  df[i,3] <- AICc(rsf_2km)
  df[i,4] <- AICc(rsf_3km)
  df[i,5] <- AICc(rsf_4km)
  df[i,6] <- AICc(rsf_5km)
  df[i,7] <- AICc(rsf_10km)
  df[i,8] <- AICc(rsf_15km)
  df[i,9] <- AICc(rsf_20km)
  
  # round to 2 decimal places
  df <- df %>% 
    mutate(across(where(is.numeric), round, digits=2))
  
}

# Inspect; warnings may indicate singularity
df
#df$lowest <- names(df[3:ncol(df)])[apply(df[3:ncol(df)], MARGIN = 1, FUN = which.min)]



# Create final data frame ------------------------------------------------------

# Select variables at their characteristic scales
rsf_df <- habitat_df_2km
rsf_df$pasture_prop <- habitat_df_20km$pasture_prop

# Remove unnecessary columns
rsf_df <- rsf_df %>% 
  dplyr::select(-scale, -BL3)

# Saved as 'data/HUGO_RSF_habitat_df.csv'


# Logistic RSF ----------------------------------------------------------------- 

# Load manuscript results
rsf_df <- read.csv('data/HUGO_RSF_habitat_df.csv')

# Standardize
rsf_df_stz <- rsf_df %>% 
  mutate_at(sv_vars,  ~(scale(.) %>% as.vector))

# Check for correlation (>0.70)
cor(rsf_df_stz %>% dplyr::select(all_of(sv_vars)))

# Generate vector of scales
ch_scale = rep(2, length(sv_vars))

# Inspect coefficients
sv_vars %>%
  map(~ get_betas(.x, 'glm')) %>%
  map_dfr(~bind_rows(.x)) %>% 
  mutate(scale = ch_scale, species = "HUGO") %>%
  relocate(species, scale)

# Inspect model fits 
glm(used ~ RLO_prop, family = binomial(link=logit), 
    data = rsf_df_stz, weights = ipp_w)
glm(used ~ farming_prop, family = binomial(link=logit), 
    data = rsf_df_stz, weights = ipp_w)
glm(used ~ forcat_prop, family = binomial(link=logit), 
    data = rsf_df_stz, weights = ipp_w)
glm(used ~ elevation, family = binomial(link=logit), 
    data = rsf_df_stz, weights = ipp_w)



