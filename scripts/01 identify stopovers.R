#############################################X
#-------- Identify Stopover Events ----------X
#############################################X
#-------- from movement tracks --------------X
#############################################X

library(dplyr)
library(tidyr)
library(lubridate)
library(terra)
library(sf)
library(mapview)
library(EMbC)
library(geosphere)
library(purrr)
source('scripts/05 functions.R')


# Load movement tracks  ----------------------------------------------------------

# Load file
# note: includes id (individual identifier) & a yearly_id (unique identifier for each southward track)
filepath <- ('data/HUGO processed tracks.csv') 
sp_df <- read.csv(filepath)



# Identify stationary locations  -----------------------------------------------

# Format timestamps
sp_df$timestamp <- ymd_hms(sp_df$timestamp)

# Stack tracks 
# (performs a species-level clustering) 
path_stk <- sp_df %>% 
  dplyr::select(timestamp, location.long, location.lat, individual.local.identifier) %>%
  arrange(timestamp) %>%
  rename(dTm = timestamp, 
         lon = location.long,
         lat = location.lat) %>%
  group_split(individual.local.identifier)

n_indivs = length(path_stk)

indivs_list <- list()
for (i in seq_along(1:n_indivs)) {
  indivs_list[[i]] = as.data.frame(path_stk[[i]])
}

# Clustering
# note: uses default speed threshold of 40 ms-1
mystk <- stbc(indivs_list)   


# Inspect overall results
#slotNames(mystk)
#stts(mystk) # general statistics: mean and sd; kn = number data points
#mystk@bC@R # delimiters
#mystk@bC@U # likelihood weights
#mystk@bCS[[1]]@pth # original data
#mystk@bCS[[1]]@dst # distance
#mystk@bCS[[1]]@spn # timespan
#mystk@bCS[[1]]@hdg # heading 


# Extract results for each individual
for (i in seq_along(1:n_indivs)){
  d <- slct(mystk, i)
  assign(paste0("indiv", i), d) 
}

#sctr(indiv1) 
#indiv2@R # delimiters (min and max values that limit each region -- same for all indivs)
#cnfm(stbc(expth1, info=-1), slct(mystck, 1)) # side-by-side
#options(scipen=999)
#stats8 <- indiv8@pth 
#indiv8@X

# Create data frame for each individual
for (i in seq_along(1:n_indivs)){
  e <- get(paste0("indiv", i))
  f <- e@pth
  clus <- e@A
  g <- f %>% 
    mutate(clust_id = as.factor(clus))
  assign(paste0("indiv", i, "_df"), g) 
}

# View individual results
#indiv4_sf <- st_as_sf(indiv4_df, coords = c("lon", "lat"), crs = 4326)
#mapview(indiv4_sf, zcol = "clust_id")
#rm(indiv4_sf)


# Produce a final combined data frame
all_df <- indiv1_df
for (i in seq(from = 2, to = n_indivs)){
  h <- get(paste0("indiv", i, "_df"))
  all_df <- all_df %>% 
    bind_rows(h)
}

all_clust <- all_df %>%
  rename(timestamp = dTm,
         location.long = lon,
         location.lat = lat) %>%
  left_join(x = sp_df)



# Group into stopovers----------------------------------------------------------

all_stops <- all_clust %>%
  # select only low-speed locations 
  filter(clust_id == 1 | clust_id == 2) %>%
  # group by annual track
  group_by(yearly_id) %>%
  # calculate distance from the previous point
  mutate(dist_from = c(NA, geosphere::distHaversine(cbind(location.long, location.lat))),
         # if greater than 10km, call this a "jump" point -- i.e., jump to new stopover
         jump = if_else(dist_from > 10000, 1, 0),
         # replace NA with 1 -- this allows for upcoming tidyr fill function
         jump = replace_na(jump, 1)) %>%
  # number the stopovers
  mutate(stop = as.factor(ifelse(jump == 1, cumsum(jump), NA))) %>%
  # fill missing values based on previous value
  tidyr::fill(stop) %>%
  ungroup() %>% 
  group_by(yearly_id, stop) %>%
  # select only stops that contain >1 location
  filter(n()>1) %>%
  ungroup() %>%
  # order by yearly_id
  arrange(yearly_id)


# FILTER: by location ----------------------------------------------------------

# Convert to sf
all_sf <- st_as_sf(all_stops, 
                   coords = c("location.long", "location.lat"), crs = 4326)

# Load AOI
ab <- st_read('data/SNAPP_AB.shp') %>%
  st_set_crs(4326) %>%
  st_make_valid()

# Extract only groups in AOI
all_stopsAB <- st_filter(all_sf, ab)


# Merge Spatially Overlapping Groups  ------------------------------------------

# Merge
all_stopsAB_merged <- merge_overlaps(all_stopsAB)

# Inspect results; update manually if necessary 
#head(all_stopsAB_merged)
#mapview(all_stopsAB, zcol = "stop")
#mapview(all_stopsAB_merged, zcol = "stop")


# FILTER: by duration ----------------------------------------------------------

# Calculate duration
all_stops <- all_stopsAB_merged %>% 
  group_by(yearly_id, stop) %>%
  mutate(start = first(timestamp),
         end = last(timestamp),
         duration = as.numeric(difftime(end, start, units = "days")),
         location.long = unlist(map(geometry,1)),
         location.lat = unlist(map(geometry,2))) %>%
  st_drop_geometry() %>%
  ungroup() 

# Filter 
all_stops <- all_stops %>%
  # remove groups that are too brief (less than 12 hrs) 
  filter(duration >= 0.5) %>%
  # remove groups that begin after Jan 1 (too late in the season)
  filter(!(month(start) < 6)) 


# ADD BACK: last locations -----------------------------------------------------

# the last location in a stopover is sometimes marked non-stationary & not included
# inspect to determine if these should be re-grouped with the stopover 

# List final locations for each stopover
last_locs_sf <- all_stops %>%
  group_by(yearly_id, stop) %>%
  arrange(timestamp) %>%
  slice(n()) %>%
  ungroup() %>%
  st_as_sf(coords = c("location.long", "location.lat"), crs = 4326)

# Extract high-velocity locations 
fast_locs_sf <- all_clust %>%
  filter(clust_id == 3 | clust_id == 4) %>%
  st_as_sf(coords = c("location.long", "location.lat"), crs = 4326)

# Loop to check if any high-velocity locations are w/in 5km of a last location 
add_back <- c()
for (i in 1:nrow(last_locs_sf)) {
  r <- last_locs_sf[i,]
  # t1 = time of last location currently in the data
  t1 <- r$timestamp
  # id1 = id of individual 
  id1 <- r$yearly_id
  # get distances from t1 to all high-velocity locations
  dists <- st_distance(r$geometry, fast_locs_sf$geometry)
  # select only if w/in 5 km
  d <- which(as.vector(dists) < 5000)
  # may return more than 1; choose the last to maximize stopover duration
  if (length(d) > 1) {
    d <- tail(d, 1)
  } else if (length(d) == 1) {
    g <- fast_locs_sf[d,]
    # t2 = time of high-velocity location
    t2 <- fast_locs_sf[g,]$timestamp
    id2 <- g$yearly_id
    # is this location after the currently included last location?
    if (id1 == id2) {
      if (t2 > t1) {
        add_back <- bind_rows(add_back, g)
      }  
    }
  }
}


# Inspect -- disregard if null
add_back_df <- add_back %>%
  mutate(location.long = unlist(map(geometry,1)),
         location.lat = unlist(map(geometry,2))) %>%
  st_drop_geometry() %>%
  filter(!(month(timestamp) < 6)) %>%
  distinct()

# Merge back into the stopovers if present
all_stops <- all_stops %>%
  dplyr::select(-start, -end, -duration, -dist_from) %>%
  bind_rows(add_back_df) %>%
  group_by(yearly_id) %>%
  arrange(timestamp, .by_group=TRUE) %>%
  mutate(stop = ifelse(is.na(stop), lag(stop), stop)) %>%
  ungroup() %>%
  group_by(yearly_id, stop) %>%
  mutate(dist_from = c(NA, geosphere::distHaversine(cbind(location.long, location.lat))),
         n_locs = n(),
         start = first(timestamp),
         end = last(timestamp),
         duration = as.numeric(difftime(end, start, units = "days"))) %>%
  rename(id = individual.local.identifier,
         lon = location.long,
         lat = location.lat) %>%
  ungroup() %>%
  dplyr::select(id, species, yearly_id, timestamp, lon, lat, lc, 
                argos.semi.major, argos.semi.minor, argos.error.radius, argos.orientation,
                gps.fix.type.raw, clust_id, dist_from, jump, stop, n_locs, start, end, duration)
  

# Label stops -----------------------------------------------------------------

# Apply unique & sequential stopover id
all_stops <- all_stops %>%
  group_by(yearly_id, stop) %>%
  dplyr::mutate(stop_id = cur_group_id()) %>%
  dplyr::relocate(stop_id, .after = species) %>%
  ungroup() %>%
  dplyr::select(-stop, -yearly_id, -clust_id, -jump) %>%
  arrange(species, stop_id)


# FILTER: location quality ------------------------------------------------------

# Create stopvoer event table
stops_table <- all_stops %>% 
  group_by(stop_id) %>%
  count(lc) %>%
  pivot_wider(names_from = lc, values_from = n) %>%
  relocate('1', .after = '2') 

# Create location class (LC) table
locs_table <- all_stops %>%
  group_by(stop_id) %>%
  slice(1) %>%
  left_join(stops_table) %>%
  dplyr::select(species, id, stop_id, start, end, duration, G, '3', '2', '1', A, B, '0')

# List stops that have only low-LC locations
low_qual_vec <- locs_table %>%
  rename('LC_G' = 'G',
         'LC_3' = '3',
         'LC_2' = '2',
         'LC_1' = '1') %>%
  filter(is.na(LC_G)) %>%
  filter(is.na(LC_3)) %>%
  filter(is.na(LC_2)) %>%
  filter(is.na(LC_1)) %>%
  dplyr::select(stop_id) %>%
  unlist(use.names = FALSE)

# Filter stops that only have poor-quality locations
all_stops <- all_stops %>% 
  filter(!(stop_id %in% low_qual_vec))


# FILTER: transmitter failures & overwintering ---------------------------------

# identify & remove via manual inspection 
# no instances in HUGO

# via manual inspection; no instances in HUGO


# Save results -----------------------------------------------------------------

# saved as 'data/HUGO stopover locations.csv'


# Create stopover centroids ----------------------------------------------------

centroids <- all_stops %>%
  st_as_sf(all_stops, 
           coords = c('lon', 'lat'), 
           crs = 4326)
  group_by(stop_id) %>% 
  # exclude low-LC locations from centroid creation
  filter(!(lc %in% c('A', 'B', '0'))) %>%
  summarize(geometry = st_union(geometry)) %>%
  st_centroid() %>%
  mutate(lon = unlist(map(geometry,1)),
         lat = unlist(map(geometry,2))) %>%
  st_drop_geometry() 

# Join stopover information
stop_info <- all_stops %>% 
    dplyr::select(id, species, stop_id, start, end, duration) %>% 
    distinct()  
centroids <- left_join(stop_info, centroids) 

# saved as 'data/HUGO stopover centroids.csv'




