#############################################X
#-- Small-Scale (Lev III) Habitat Selection -X
#############################################X
#------------ Prepare  & Model --------------X
#############################################X


library(tidyverse)
library(sf)
library(purrr)
library(mapview)



# Load STEP data ---------------------------------------------------------------

# data frame includes:
# (1) 'used' steps (transmitter-reported locations) -- (used = 1, sim = 0)
# (2) 'available' steps (randomly generated locations) -- (used = 0, sim = 0)
# (3) 'alternative' steps (simulated by ctmm) -- (used = NA, sim = 1)

all_steps <- read.csv('data/HUGO SSF steps.csv') %>%
  mutate(timestamp = parse_date_time(timestamp, orders = "ymd HMS", truncated = 3),
         year = year(timestamp),
         yearly_id = paste0(id, year))

# no MB data currently available for 2022 -- change the year to 2021
all_steps <- all_steps %>% 
  mutate(year = ifelse(year == 2022, 2021, year))

# check
all_steps[which(is.na(all_steps$timestamp)),] 

# Load GEE data -----------------------------------------------------------------

# JRC: load & tidy
jrc <- read.csv('data/SSF_JRC.csv') %>% 
  dplyr::select(loc_id, occurrence, recurrence, seasonality, 
                contains('_Obsv'), contains('_Rec'))
jrc_env <- left_join(all_steps, jrc, by = 'loc_id')

# MAPBIOMAS: load & tidy
mb <- read.csv('data/SSF_MB_Combined.csv') %>% 
  dplyr::select(loc_id, forest:year)
mb_env <- left_join(all_steps, mb, by = c('loc_id', 'year'))


# DEM: load & tidy
dem <- read.csv('data/SSF_DEM.csv') %>% 
  dplyr::select(loc_id, elevation, slope)
dem_env <- left_join(all_steps, dem, by = 'loc_id')

# combine
ssf_env <- left_join(jrc_env, mb_env) %>% left_join(., dem_env) 

# Scale continuous variables
ssf_env$unscaled_occurrence <- ssf_env$occurrence
ssf_env$occurrence <- scale(ssf_env$occurrence)
ssf_env$elevation <- scale(ssf_env$elevation)
ssf_env$slope <- scale(ssf_env$slope)

# Check results
head(ssf_env)
nrow(all_steps) == nrow(ssf_env) # check for missing data
unique(ssf_env$stop_id) # check for stop edits


# Select Ratio -----------------------

# to 1:5
used <- ssf_env %>% 
  filter(used == 1)
av <- ssf_env %>% 
  filter(used == 0) %>% 
  group_by(pair) %>% 
  slice_sample(n = 5) 
reduced_ssf <- bind_rows(used, av)


# Create data frame for modeling ------------------------------------------------

mod_ssf <- reduced_ssf %>% 
  #dplyr::select(-(recurrence:Sept_Rec)) %>%
  #dplyr::select(forest:not_Obs) %>%
  #rowwise() %>%
  mutate(across(forest:not_Obs, function (x) if_else(x > 0, cur_column(), NA))) %>%
  mutate(LULC = coacross(forest:not_Obs)) %>%
  dplyr::select(loc_id, id, species, stop_id, lon, lat, timestamp, lc, used, pair, sim, yearly_id,
                LULC, occurrence, elevation, slope, unscaled_occurrence) 



# Create Used:Available Tables -------------------------------------------------

lc_tbl <- reduced_ssf %>%
  pivot_longer(cols = forest:not_Obs, names_to='LC', values_to = 'LC_use') %>%
  mutate(LC = case_when(LC == 'flooded_Forest' ~ 'forest',
                        LC == 'forest_Formation' ~ 'forest', 
                        LC == 'savanna' ~ 'forest',
                        LC == 'urban' ~ 'other',
                        LC == 'other_Nonveg' ~ 'other',
                        LC == 'other_NF_Natural' ~ 'other',
                        .default = LC)) %>%
  filter(LC_use == 1) %>%
  mutate(used = as.factor(used)) %>%
  group_by(LC, used) %>%
  dplyr::summarize(count = n()) %>%
  ungroup()

totals <- lc_tbl %>%
  group_by(LC) %>%
  dplyr::summarize(total = sum(count)) %>%
  ungroup()

lc_tbl <- left_join(lc_tbl, totals) %>% 
  mutate(pct = count/total) 

lc_list <- unique(lc_tbl$LC)

# Create empty used category for LULCs that individuals never used 
never_used <- vector()
tpix <- vector()

for(i in 1:length(lc_list)){
  
  lc_name <- lc_list[i]
  
  f <- lc_tbl %>% filter(LC == lc_name)
  t <- unique(f$total)
  
  if(nrow(f) == 1){
    never_used <- c(never_used, lc_name)
    tpix <- c(tpix, t)
  } 
}

if(length(never_used) > 0){
  lc_tbl <- lc_tbl %>% 
    bind_rows(data.frame(LC = never_used, 
                         used = as.factor(1), 
                         count = 0, 
                         total = tpix, 
                         pct = 0))
}


# Inspect
#lc_tbl %>% print(n = 50)

# Logistic SSF  -----------------------------------------------------------------

# Check plot (see script 15)
hugo_lc_plot

# Check df (should be correct species)
mod_ssf 

# Simplify LULC
mod_ssf_simple <- mod_ssf %>% 
  mutate(LULC = if_else(LULC == 'RLO', LULC, 'other')) %>% 
  mutate(LULC = factor(LULC, levels = c('other', 'RLO')))

# Fit mods
m.LULC<- clogit(used ~ LULC + strata(pair), 
                data = mod_ssf_simple, 
                method = 'efron', 
                robust=T, 
                cluster = yearly_id)
m.linear <- clogit(used ~ occurrence + LULC + strata(pair), 
                   data = mod_ssf_simple, 
                   method = 'efron', 
                   robust=T, 
                   cluster = yearly_id)
m.quad <- clogit(used ~ occurrence + I(occurrence^2) + LULC + strata(pair), 
                 data = mod_ssf_simple, 
                 method = 'efron', 
                 robust=T, 
                 cluster = yearly_id)


# Make an AICc table
sp.aic <- c(AICc(m.linear), AICc(m.quad), AICc(m.LULC))
mods <- c('LULC+linear_occ', 'LULC+quadratic_occ', 'LULC')
delAIC <- sp.aic - min(sp.aic)
relLik <- exp(-0.5 * delAIC)
aicweight <- relLik/sum(relLik)
sp.aic.table <- data.frame(AICc = sp.aic, delAIC = delAIC, relLik = relLik, weight = aicweight)
round(sp.aic.table, digits = 3) %>% mutate(models = mods) %>% relocate(models) %>% arrange(AICc) 

# Review
summary(m.linear) 
summary(m.LULC)
m.hugo <- m.linear

# Make table of model coefficients
hugo.summary <- broom::tidy(m.hugo) %>% 
  mutate(confint.25 = estimate-1.96*robust.se,
         confint.975 = estimate+1.96*robust.se,
         exp.estimate = exp(estimate),
         exp.confint.25 = exp(confint.25),
         exp.confint.975 = exp(confint.975))

# Check for non-convergence/quasi-separation
hugo.summary



