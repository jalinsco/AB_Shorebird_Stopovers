#############################################X
#------ SIMULATE TRACKS for SSF -------------X
#############################################X

library(dplyr)
library(tidyr)
library(lubridate)
library(terra)
library(sf)
library(ggplot2)
library(mapview)
library(mvtnorm)
library(ctmm)
library(units)


# Load data --------------------------------------------------------------------

# Load stopovers list

stop_locs <- read.csv("data/HUGO stopover locations.csv") %>%
  rename(individual.local.identifier = id,
         location.long = lon,
         location.lat = lat) %>%
  mutate(sensor.type = ifelse(lc == "G", "gps", "argos-doppler-shift"),
         timestamp = ymd_hms(timestamp))


# Prepare data frame ------------------------------------------------------------

# Exclude individuals/data with missing argos error ellipse info
stop_locs <- stop_locs %>% 
  filter(!(stop_id == 66 & sensor.type == "argos-doppler-shift"))

# create a list of distinct individuals
indivs_list <- stop_locs %>% 
  distinct(individual.local.identifier) %>% 
  pull()


# Simulate alternative steps ---------------------------------------------------

# used to evaluate the robustness of conclusions to location error 
# Workflow amended from Karagacheva et al. 2023

for (i in 1:length(indivs_list)){
  get_indiv <- indivs_list[i]
  
  indiv <- stop_locs %>% dplyr::filter(individual.local.identifier == get_indiv)
  
  # get projection
  #z <- indiv %>% slice(1) %>% dplyr::select(UTMzone)
  
  # select gps locs
  indiv_gps <- indiv %>% 
    filter(sensor.type == "gps") 
  
  # select argos locs
  indiv_argos <- indiv %>% 
    filter(sensor.type == "argos-doppler-shift")
  
  # check if track contains argos & gps locations
  if (nrow(indiv_gps) > 0 & nrow(indiv_argos) > 0) {
    indiv_gps <- indiv_gps %>%
      dplyr::select(individual.local.identifier, timestamp, 
                    location.long, location.lat, gps.fix.type.raw)
    indiv_gps_telemetry <- as.telemetry(indiv_gps)
    
    # set priors for uncalibrated error -- if 2D and 3D 
    if (length(unique(indiv_gps$gps.fix.type.raw)) == 2) {
      uere(indiv_gps_telemetry) <- c(20,10)
    } else {
      # or if only 3D
      uere(indiv_gps_telemetry) <- 10
    }
    uere(indiv_gps_telemetry)$DOF[] <- 2
    
    indiv_argos <- indiv_argos %>% 
      dplyr::select(-(gps.fix.type.raw))
    indiv_argos_telemetry <- as.telemetry(indiv_argos) 
    
    indiv_telemetry <- tbind(indiv_gps_telemetry, indiv_argos_telemetry)
    
  # if only gps locations  
  } else if(nrow(indiv_gps) > 0){
    indiv_gps <- indiv_gps %>%
      dplyr::select(individual.local.identifier, timestamp, 
                    location.long, location.lat, gps.fix.type.raw)
    indiv_telemetry <- as.telemetry(indiv_gps)
    if (length(unique(indiv_gps$gps.fix.type.raw)) == 2) {
      uere(indiv_telemetry) <- c(20,10)
    } else {
      uere(indiv_telemetry) <- 10
    }
    uere(indiv_telemetry)$DOF[] <- 2 
    
  # if only argos locations  
  } else {
    indiv_argos <- indiv_argos %>% 
      dplyr::select(-(gps.fix.type.raw))
    indiv_telemetry <- as.telemetry(indiv_argos) 
  }
  
  GUESS <- ctmm.guess(indiv_telemetry, CTMM=ctmm(error=TRUE), interactive=FALSE)
  FIT <- ctmm.select(indiv_telemetry,GUESS,trace=TRUE,cores=6)
  
  # save FIT results
  #filepathW = paste0("data/simulations/", get_indiv, "_FIT.RDS")
  #saveRDS(FIT, file=filepathW)

  
  # simulate possible locations by scattering locs in error ellipse
  n.simulations = 100  
  
  all_out <- c()
  for (i in 1:n.simulations) {
    cat('\r', i)
    sim_new <- ctmm::simulate(FIT, data=indiv_telemetry, precompute=TRUE, t=indiv_telemetry$t)
    sim_out <- st_as_sf(x=data.frame(ID=1:length(sim_new@.Data[[2]]), 
                                     x=sim_new@.Data[[2]], y=sim_new@.Data[[3]]), 
                        coords=c('x', 'y'),  crs = projection(sim_new))  %>%
      st_transform(st_crs('EPSG:4326')) %>% 
      st_coordinates()
    all_out <- rbind(all_out, sim_out)
  }
  
  sim_out_df <- all_out %>% 
    as.data.frame() %>%
    mutate(loc = rep(as.character(1:nrow(indiv_df)), 100)) %>%
    rename(location.long = "X", location.lat = "Y") 
  
  sim_final <- left_join(indiv_df, sim_out_df) %>% 
    dplyr::select(-loc)  
  
  # save results
  #filepathW = paste0("data/", get_indiv, "_SIM.csv")
  #write.csv(sim_final, filepathW, row.names=F)
}



# Simulation files
simFiles = list.files(path = 'data/', 
                      pattern = '*_SIM.csv', 
                      full.names = T)
sims_raw <- lapply(simFiles, read_csv, col_types = c('c', 'c', 'd', 'T', 'c', 'd', 'd')) %>% 
  bind_rows()

# Create data frame of 'alternative' steps
sims <- sims_raw %>% 
  mutate(used = 0, sim = 1) %>%
  rename(id = individual.local.identifier, lon = 
           location.long, lat = 
           location.lat) %>%
  mutate(timestamp = ymd_hms(timestamp))

# Check for NA values in simulations
any(is.na(sims$lon))

# Load 'used' steps
stop_locs <- read.csv('data/HUGO stopover locations.csv') %>%
  filter(!(stop_id == 66 & lc != 'G')) %>%
  mutate(timestamp = ymd_hms(timestamp), 
         used = 1, 
         sim = 0)

# Combine 'Used' (empirical) and 'Alternative' (simulated) steps  
ssf_df <- bind_rows(stop_locs, sims)