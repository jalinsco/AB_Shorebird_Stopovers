#############################################X
#-- Large-Scale (Lev II) Habitat Selection --X
#############################################X
#------------ Prepare  & Model --------------X
#############################################X

library(sf)
library(spatstat)
library(ggplot2)
library(pROC)
library(geojsonsf)
library(MuMIn)
library(lubridate)
library(logistf)
library(spdep)
library(broom.mixed)
library(DHARMa)
library(lmtest)
library(sandwich)
source('scripts/05 functions.R')


# Load Data --------------------------------------------------------------------

# SNAPP Amazon Basin shapefiles (Venticinque et al. 2016)
# see: https://snappartnership.net/teams/amazon-waters/
ab <- st_geometry(st_read('data/Amazon Basin Shapefiles/SNAPP_AB.shp'))
basins <- st_read('data/Amazon Basin Shapefiles/SNAPP_Subbasins.shp')

# Available points
av_basins <- read.csv('data/available_locations.csv') 

# Habitat data: Mapbiomas (downloaded from GEE; RAISG 2023)
# see: https://amazonia.mapbiomas.org/ 
mb_used_all <- read.csv('data/RSF_MB_USED.csv') 
mb_av_all <- read.csv('data/RSF_MB_AV.csv') 

# Habitat data: JRC global surface water mapping (downloaded from GEE; Pekel et al. 2016) 
# see: https://global-surface-water.appspot.com/
jrc_used_all <- read.csv('data/RSF_JRC_USED.csv') 
jrc_av_all <- read.csv('data/RSF_JRC_AV.csv') 

# Habitat data: elevation 
# see: NASA Shuttle Radar Topography Mission (SRTM; Farr et al. 2007)
dem_used <- read.csv('data/RSF_DEM_5km_USED.csv') 
dem_av <- read.csv('data/RSF_DEM_5km_AV.csv') 


# PREPARE ----------------------------------------------------------------------

# Identify visited subbasins
filepathR <- 'data/HUGO processed tracks.csv'
sp_df <- read.csv(filepathR) %>% 
  dplyr::select(yearly_id, timestamp, location.long, location.lat)

sp_sf_polylines <- sp_df %>%
    mutate(timestamp = ymd_hms(timestamp)) %>% 
    st_as_sf(coords = c('location.long', 'location.lat'), crs = 4326) %>%
    group_by(yearly_id) %>% 
    dplyr::summarize(do_union=FALSE) %>% 
    st_cast("LINESTRING")
  
sp_basins <- st_join(sp_sf_polylines, basins, join = st_intersects) %>%
    st_drop_geometry() %>%
    distinct(yearly_id, BL3) %>%
    na.omit()
  

# Prepare data frame -----------------------------------------------------------

# Vector of sacles
scales = c(2, 3, 4, 5, 10, 15, 20)

# Loop to create data frames for each scale
for(i in 1:length(scales)){
  
  # used DF -- filter to scale of selection 
  mb_used <- mb_used_all %>% 
    filter(scale == scales[i]) %>% 
    dplyr::select(-scale) 
  jrc_used <- jrc_used_all %>% 
    filter(scale == scales[i]) %>% 
    dplyr::select(-scale)
  used_vars <- left_join(mb_used, jrc_used) %>% 
    left_join(dem_used)
  
  # used DF -- filter to species
  used_vars_sp <- used_vars %>% 
    left_join(centroids) %>% 
    filter(species == "HUGO") %>% 
    mutate(yearly_id = paste0(id, year(start)))
  
  # available DF -- filter to scale of selection   
  mb_av <- mb_av_all %>% 
    filter(scale == scales[i]) %>% 
    dplyr::select(-scale)
  jrc_av <- jrc_av_all %>% 
    filter(scale == scales[i]) %>% 
    dplyr::select(-scale)
  av_vars <- left_join(mb_av, jrc_av) %>% 
    left_join(dem_av) %>% left_join(av_basins)

  # available DF -- filter to basins visited by individual
  indivs <- used_vars_sp %>% 
    pull(yearly_id) %>% 
    unique()
  
  if (scales[i] == 2){ 
    av_vec = c()
    id_vec = c()
    
    for(i in 1:length(indivs)){
      b <- sp_basins %>% filter(yearly_id == indivs[i])
      n <- used_vars_sp %>% filter(yearly_id == indivs[i]) # number of used 
      sz <- 500*nrow(n) # select 500 av for each used 
      v <- rep(indivs[i], sz)
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



# Identify characteristic Scale ----------------------

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
  df <- df %>% mutate(across(where(is.numeric), round, digits=2))
  
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


# Logistic RSF ----------------------------------------------------------------- 

# Standardize
rsf_df_stz <- rsf_df %>% mutate_at(sv_vars,  ~(scale(.) %>% as.vector))

# Check for correlation (>0.70)
cor(rsf_df_stz %>% dplyr::select(all_of(sv_vars)))

# Generate vector of scales
ch_scale = rep(2, length(sv_vars))

# Fit
betas <- sv_vars %>%
  map(~ get_betas(.x, 'glm')) %>%
  map_dfr(~bind_rows(.x)) %>% 
  mutate(scale = ch_scale, species = "HUGO") %>%
  relocate(species, scale)

# Inspect coefficients
betas

# Inspect individual fits 
glm(used ~ RLO_prop, family = binomial(link=logit), 
    data = rsf_df_stz, weights = ipp_w)
glm(used ~ farming_prop, family = binomial(link=logit), 
    data = rsf_df_stz, weights = ipp_w)
glm(used ~ forcat_prop, family = binomial(link=logit), 
    data = rsf_df_stz, weights = ipp_w)
glm(used ~ elevation, family = binomial(link=logit), 
    data = rsf_df_stz, weights = ipp_w)



