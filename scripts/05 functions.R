############################################X
#---------------- FUNCTIONS ----------------X
############################################X
#----- A Collection of the Very Largest ----X
############################################X



## Preparing Tracks  -----------------------------------

# removes duplicate locations (same individual, same timestamp)
remove_dupes <- function(data) {
  data = data %>%
    group_by(individual.local.identifier, timestamp) %>%
    mutate(dupe = n()>1) %>%
    mutate(quality = case_when(lc == 'G' ~ 7,
                               lc == 3 ~ 6,
                               lc == 2 ~ 5,
                               lc == 1 ~ 4,
                               lc == 0 ~ 3,
                               lc == 'A' ~ 2,
                               lc == 'B' ~ 1,
                               lc == 'Z' ~ 0)) %>%
    ungroup()
  isdistinct = data %>% 
    filter(dupe != TRUE)
  best_quality_dupes = data %>% 
    filter(dupe == TRUE) %>%
    group_by(individual.local.identifier, timestamp) %>% 
    slice_max(quality, with_ties = FALSE) %>%
    ungroup()
  return(bind_rows(isdistinct, best_quality_dupes))
}

# prepares an analysis-ready data frame 
arrange_dfs <- function(data, species_name) {
  data = data %>%
    mutate(year = year(timestamp),
           lc = ifelse(sensor.type == 'gps', 'G', argos.lc))
  if (species_name == "hugo"){
    data = data %>%
      mutate(species = "HUGO",
             yearly_id = paste0(individual.local.identifier, year)) %>%
      dplyr::select(individual.local.identifier,
                    timestamp,
                    lc,
                    location.long,
                    location.lat,
                    argos.semi.major,
                    argos.semi.minor,
                    argos.error.radius,
                    argos.orientation,
                    gps.fix.type.raw,
                    yearly_id,
                    species) 
  } else if (species_name == "leye"){
    data = data %>%
      mutate(species = "LEYE") %>%
      mutate(yearly_id = paste0(individual.local.identifier, year)) %>%
      dplyr::select(individual.local.identifier,
                    timestamp,
                    lc,
                    location.long,
                    location.lat,
                    argos.semi.major,
                    argos.semi.minor,
                    argos.error.radius,
                    argos.orientation,
                    gps.fix.type.raw,
                    yearly_id,
                    species) 
  } else if (species_name == "amgo"){
    data = data %>%
      mutate(species = "AMGO",
             individual.local.identifier = paste0('A', individual.local.identifier),
             yearly_id = paste0(individual.local.identifier, year)) %>%
      dplyr::select(individual.local.identifier,
                    timestamp,
                    lc,
                    location.long,
                    location.lat,
                    gps.fix.type.raw,
                    yearly_id,
                    species) 
  } else if (species_name == "pesa"){
    data = data %>%
      mutate(species = "PESA",
             individual.local.identifier = paste0('P', individual.local.identifier),
             yearly_id = paste0(individual.local.identifier, year)) %>%
      dplyr::select(individual.local.identifier,
                    timestamp,
                    lc,
                    location.long,
                    location.lat,
                    gps.fix.type.raw,
                    yearly_id,
                    species) %>%
      filter(lc != "Z") %>%
      mutate(argos.semi.major = NA,
             argos.semi.minor = NA,
             argos.error.radius = NA,
             argos.orientation = NA)
  } else if (species_name == "upsa"){
    data = data %>%
      mutate(species = "UPSA",
             individual.local.identifier = paste0('U', tag.local.identifier),
             yearly_id = paste0(individual.local.identifier, year)) %>%
      dplyr::select(individual.local.identifier,
                    timestamp,
                    lc,
                    location.long,
                    location.lat,
                    argos.semi.major,
                    argos.semi.minor,
                    argos.error.radius,
                    argos.orientation,
                    gps.fix.type.raw,
                    yearly_id,
                    species) 
  } else if (species_name == "bbsa"){
    data = data %>%
      mutate(species = "BBSA",
             yearly_id = paste0(individual.local.identifier, year)) %>%
      dplyr::select(individual.local.identifier,
                    timestamp,
                    lc,
                    location.long,
                    location.lat,
                    argos.semi.major,
                    argos.semi.minor,
                    argos.error.radius,
                    argos.orientation,
                    gps.fix.type.raw,
                    yearly_id,
                    species) %>%
      filter(lc != "Z")
  }
}


# Merges overlapping stopovers 
# note 7/2023: this is now problematic!!! after pkg update, does not remove all overlaps
# must check for departures & returns to the same area (e.g., manually in QGIS) 

merge_overlaps = function(x){
  yearly_id_vec = unique(x$yearly_id)
  overlaps = data.frame()
  for (id in seq_along(yearly_id_vec)){
    y = x %>%
      filter(yearly_id == yearly_id_vec[id]) 
    dist.mat = y %>%
      st_distance()
    n = y %>%
      nrow()
    
    ijd = data.frame(expand.grid(i=1:n, j=1:n)) 
    ijd$distance = as.numeric(c(dist.mat))
    ijd$yearly_id = y$yearly_id[ijd$i]
    ijd$stop_i = y$stop[ijd$i]
    ijd$stop_j = y$stop[ijd$j]
    
    v <- ijd %>%
      filter(stop_i != stop_j) %>%
      filter(distance < 10000) %>%
      dplyr::select(yearly_id, stop_i, stop_j)
    
    if (dim(v)[1] != 0) {
      overlaps <- overlaps %>%
        bind_rows(v) %>%
        distinct()
    }
  }
  if (dim(overlaps)[1] != 0) {
    # this part is the problem
    z = overlaps %>%
      filter(row_number() %% 2 == 1)
    for (i in 1:nrow(z)){
      x <- x %>%
        mutate(stop = if_else(yearly_id == z[i,]$yearly_id & stop == z[i,]$stop_j, z[i,]$stop_i, stop))
    }
    w <- x %>%
      group_by(yearly_id) %>%
      arrange(stop, timestamp, .by_group=TRUE) %>%
      ungroup()
    return(w)
  }
  else {
    return(x)
  }
}

## SUMMARIZE INDIVIDUALS -------------------------------------

# add column if it does not eexist
fncols <- function(data, cname) {
  add <-cname[!cname%in%names(data)]
  
  if(length(add)!=0) data[add] <- NA
  data
}


# fix names
fix_names <- function(data, species_name) {
  if (species_name == "leye"){
    data = data %>% mutate(individual.local.identifier = word(individual.local.identifier,1))
  } else if (species_name == "amgo"){
    data = data %>% mutate(individual.local.identifier = paste0('A', individual.local.identifier)) 
  } else if (species_name == "pesa"){
    data = data %>% mutate(individual.local.identifier = paste0('P', individual.local.identifier)) 
  } else if (species_name == "upsa"){
    data = data %>% mutate(individual.local.identifier = paste0('U', tag.local.identifier)) 
  } else {
    data = data 
  }
}

# add device data
add_devices <- function(data, species_name) {
  if (species_name == "leye"){
    data = data %>% mutate(device.type = 'Lotek')
  } else if (species_name == "amgo"){
    data = data %>% mutate(device.type = 'Lotek') 
  } else if (species_name == "pesa"){
    data = data %>% mutate(device.type = ifelse(study == 'Mihai', 'PTT', 'Lotek')) 
  } else if (species_name == "upsa"){
    data = data %>% 
      mutate(device.type = case_when(individual.local.identifier %in% c('U207917', 'U207916', 'U207913', 'U207912', 'U207910', 'U207911', 'U233872',
                                                                        'U233868', 'U233870', 'U233886', 'U233888') ~ 'Lotek', # Jim
                                     individual.local.identifier %in% c('U233890', 'U233889', 'U233882', 'U233881', 'U221492', 'U221493') ~ 'PTT', # Jim
                                     individual.local.identifier %in% c('U158666', 'U158667') ~ 'PTT', # Bret
                                     individual.local.identifier %in% c('U146682', 'U146694', 'U146702', 'U158662', 'U146690') ~ 'Lotek', # Bret
                                     .default = 'unknown')) 
  } else if (species_name == "bbsa"){
    bbsa_devices <- read.csv('raw_data/BBSA from Lee/buffBreastedSandpiper_USGS_ASC_argosGPS_deploymentAttributes.csv') %>%
      dplyr::select(Animal_ID, PTT_Manufacturer) %>%
      distinct() %>%
      rename(individual.local.identifier = Animal_ID,
             device.type = PTT_Manufacturer)
    data = data %>% left_join(bbsa_devices)
  } else if (species_name == 'hugo'){
    data = data %>% mutate(device.type = ifelse(individual.local.identifier == 'KCH', 'PTT', 'Lotek'))
  }
}
## RESOURCE-SELECTION FUNCTIONS ------------------------------

get_betas <- function(x, type){
  
  if(type == 'glmm'){
    unv <- tryCatch((glmmTMB(as.formula(paste0('used ~', x, '+ (1|yearly_id) + (0+ ', x, '|yearly_id)')), family = binomial(link=logit), data = habitat_df_stz, weights = ipp_w)),
                    error = function(e) {unv <- tibble(term = x)})
    rse.25 <- NA
    rse.975 <- NA
    
  } else if(type == 'glm'){
    unv <- tryCatch((glm(as.formula(paste0('used ~', x)), family = binomial(link=logit), data = habitat_df_stz, weights = ipp_w)),
                    error = function(e) {unv <- tibble(term = x)})
    rse.25 <- tryCatch(coefci(unv, vcov = vcovCL, cluster = ~yearly_id, level = 0.95)[2], error = function(e) NA)
    rse.975 <- tryCatch(coefci(unv, vcov = vcovCL, cluster = ~yearly_id, level = 0.95)[4], error = function(e) NA)
  }
  
  
  if(length(unv) > 1){
    unv <- unv %>%
      broom.mixed::tidy() %>%
      dplyr::select(term, estimate, std.error) %>% 
      slice(2) %>%
      mutate(confint.25 = tryCatch(confint(unv)[2,1], error = function(e) NA),
             confint.975 = tryCatch(confint(unv)[2,2], error = function(e) NA),
             confint.25.robust = rse.25,
             confint.975.robust = rse.975)
  }
  
  return(unv)
  
}


# truncated pareto
rpareto_t <- function(n, range, shape, scale = 1) {
  
  # distribution function -- apply CDF to lower and upper bounds to find quantiles of edges
  F.a <- ppareto(min(range), shape = shape, scale = scale)
  F.b <- ppareto(max(range), shape = shape, scale = scale)
  
  # sample from a uniform distribution between those quantiles
  u <- runif(n, min = F.a, max = F.b)
  
  # apply the inverse CDF -- quantile function
  qpareto(u, shape = shape, scale = scale)
}


## STEP-SELECTION FUNCTIONS ------------------------------

simulate_error <- function(df1, df2){
  
  coeff.1 <- vector()
  coeff.2 <- vector()
  coeff.3 <- vector()
  coeff.4 <- vector()
  
  rse.1 <- vector()
  rse.2 <- vector()
  rse.3 <- vector()
  rse.4 <- vector()
  
  s <- df1$species[1]
  
  for(i in seq_along(1:99)){
    
    df <- df1 %>% group_by(pair) %>% slice(0+i) %>% bind_rows(df2) 
    
    if(s == 'HUGO'){
      m <- clogit(used ~ RLO + occurrence + strata(pair), data = df, method = 'efron', robust=T, cluster = yearly_id)
    } else if (s == 'LEYE'){
      m <- clogit(used ~ RLO + grassland + strata(pair), data = df, method = 'efron', robust=T, cluster = yearly_id)
    } else if (s == 'AMGO'){
      m <- clogit(used ~ RLO + grassland + agriculture + strata(pair), data = df, method = 'efron', robust=T, cluster = yearly_id)
    } else if (s == 'UPSA'){
      m <- clogit(used ~ grassland + mosaic_Use + strata(pair), data = df, method = 'efron', robust=T, cluster = yearly_id)
    } else if (s == 'PESA'){
      m <- clogit(used ~ RLO + strata(pair), data = df, method = 'efron', robust=T, cluster = yearly_id)
    } else if (s == 'BBSA'){
      m <- clogit(used ~ RLO + grassland + pasture + agriculture + strata(pair), data = df, method = 'efron', robust=T, cluster = yearly_id)
    }
    
    m <- broom::tidy(m)
    n <- length(m$term)
    
    coeff.1 <- c(coeff.1, m$estimate[1])
    rse.1 <- c(rse.1, m$robust.se[1])
    
    if(n > 1) {
      coeff.2 <- c(coeff.2, m$estimate[2])
      rse.2 <- c(rse.2, m$robust.se[2])
    } 
    
    if(n > 2) {
      coeff.3 <- c(coeff.3, m$estimate[3])
      rse.3 <- c(rse.3, m$robust.se[3])
    }
    
    if(n > 3){
      coeff.4 <- c(coeff.4, m$estimate[4])
      rse.4 <- c(rse.4, m$robust.se[4])
    }
    
  }
  
  e <- data.frame(sim = 1:99, estimate.1 = coeff.1, robust_se.1 = rse.1)
  e$confint.lower.1 = e$estimate.1-1.96*e$robust_se.1
  e$confint.higher.1 = e$estimate.1+1.96*e$robust_se.1
  
  if(n > 1){
    e$estimate.2 = coeff.2
    e$robust_se.2 = rse.2
    e$confint.lower.2 = e$estimate.2-1.96*e$robust_se.2
    e$confint.higher.2 = e$estimate.2+1.96*e$robust_se.2
  }
  
  if(n > 2){
    e$estimate.3 = coeff.3
    e$robust_se.3 = rse.3
    e$confint.lower.3 = e$estimate.3-1.96*e$robust_se.3
    e$confint.higher.3 = e$estimate.3+1.96*e$robust_se.3
  }
  
  if(n > 3){
    e$estimate.4 = coeff.4
    e$robust_se.4 = rse.4
    e$confint.lower.4 = e$estimate.4-1.96*e$robust_se.4
    e$confint.higher.5 = e$estimate.4+1.96*e$robust_se.4
  }
  
  return(e)
  
}

error_influence <- function(df1, df2, s){
  
  df2 <- df2 %>% filter(species == s)
  n <- nrow(df2)
  
  coef <- df2$term
  changed_signs <- vector()
  lost_conf <- vector()
  

  if(n == 1){
    
    changed_signs <- nrow(subset(df1, df1$estimate.1 < 0))
    lost_conf <- ifelse(df2$estimate[1] > 0,  nrow(subset(df1, df1$estimate.1 <= 0)), nrow(subset(df1, df1$estimate.1 >= 0)))  
  
    } else {
      
      for(i in seq_along(1:n)){
        
        z <- df1[ , grepl("estimate", names(df1))]
        l <- df1[ , grepl("confint.lower", names(df1))]
        h <- df1[ , grepl("confint.higher", names(df1))]
      
        if(df2$sig[i] == 'no'){
          changed_signs[i] <- NA
          lost_conf[i] <- NA
          } else if (df2$estimate[i] > 0){
            changed_signs[i] <- nrow(subset(z, z[,i] <= 0))
            lost_conf[i] <- nrow(subset(l, l[,i] <= 0))
            } else if (df2$estimate < 0){
              changed_signs[i] <- nrow(subset(z, z[,i] >= 0))
              lost_conf[i] <- nrow(subset(h, h[,i] >= 0))
      }
    }
  }
  
  return(data.frame(coef = coef, changed_signs = changed_signs, lost_conf = lost_conf))

}
