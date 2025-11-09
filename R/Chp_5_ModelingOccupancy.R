# Chapter 5 modweling occupancy

# Sec 5.1 ----
library(auk)
library(lubridate)
library(sf)
library(dggridR)
library(unmarked)
library(raster)
library(ebirdst)
library(MuMIn)
library(AICcmodavg)
library(fields)
library(tidyverse)
# resolve namespace conflicts
select <- dplyr::select
projection <- raster::projection

# set random number seed to insure fully repeatable results
set.seed(1)

# setup output directory for saved results
if (!dir.exists("output")) {
  dir.create("output")
}

# ebird data
ebird <- read_csv("data/ebd_woothr_june_bcr27_zf.csv") %>% 
  mutate(year = year(observation_date),
         # occupancy modeling requires an integer response
         species_observed = as.integer(species_observed))

# modis land cover covariates
habitat <- read_csv("data/pland-elev_location-year.csv") %>% 
  mutate(year = as.integer(year))

# combine ebird and modis data
ebird_habitat <- inner_join(ebird, habitat, by = c("locality_id", "year"))

# prediction surface
pred_surface <- read_csv("data/pland-elev_prediction-surface.csv")
# latest year of landcover data
max_lc_year <- pred_surface$year[1]
r <- raster("data/prediction-surface.tif")

# load gis data for making maps
map_proj <- st_crs(102003)
ne_land <- read_sf("data/gis-data.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
bcr <- read_sf("data/gis-data.gpkg", "bcr") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_country_lines <- read_sf("data/gis-data.gpkg", "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_state_lines <- read_sf("data/gis-data.gpkg", "ne_state_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()


#Sec 5.2 Data Prep ----

#filter prior to creating occupancy model data
ebird_filtered <- filter(ebird_habitat, 
                         number_observers <= 5,
                         year == max(year))

# 5.2.1 Data formatting ----
occ <- filter_repeat_visits(ebird_filtered, 
                            min_obs = 2, max_obs = 10,
                            annual_closure = TRUE,
                            date_var = "observation_date",
                            site_vars = c("locality_id", "observer_id"))
# entire data set
nrow(ebird_habitat)
#> [1] 48450
# reduced data set
nrow(occ)
#> [1] 3724
# number of individual sites
n_distinct(occ$site)
#> [1] 988

# format for unmarked
occ_wide <- format_unmarked_occu(occ, 
                                 site_id = "site", 
                                 response = "species_observed",
                                 site_covs = c("n_observations", 
                                               "latitude", "longitude", 
                                               "pland_04_deciduous_broadleaf", 
                                               "pland_05_mixed_forest",
                                               "pland_12_cropland",
                                               "pland_13_urban"),
                                 obs_covs = c("time_observations_started", 
                                              "duration_minutes", 
                                              "effort_distance_km", 
                                              "number_observers", 
                                              "protocol_type",
                                              "pland_04_deciduous_broadleaf", 
                                              "pland_05_mixed_forest"))

#Sec. 5.2 Spatial subsampling ----

# generate hexagonal grid with ~ 5 km between cells
dggs <- dgconstruct(spacing = 5)
# get hexagonal cell id for each site
occ_wide_cell <- occ_wide %>% 
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum)
# sample one site per grid cell
occ_ss <- occ_wide_cell %>% 
  group_by(cell) %>% 
  sample_n(size = 1) %>% 
  ungroup() %>% 
  select(-cell)
# calculate the percent decrease in the number of sites
1 - nrow(occ_ss) / nrow(occ_wide)

#5.2.3 {unmarked} object ----
occ_um <- formatWide(occ_ss, type = "unmarkedFrameOccu")
summary(occ_um)


# Section 5.3 Occupancy modeling -----

# fit model
occ_model <- occu(~ time_observations_started + 
                    duration_minutes + 
                    effort_distance_km + 
                    number_observers + 
                    protocol_type +
                    pland_04_deciduous_broadleaf + 
                    pland_05_mixed_forest
                  ~ pland_04_deciduous_broadleaf + 
                    pland_05_mixed_forest + 
                    pland_12_cropland + 
                    pland_13_urban, 
                  data = occ_um)
# look at the regression coefficients from the model
summary(occ_model)


# 5.3.1 Assessment ----

occ_gof <- mb.gof.test(occ_model, nsim = 10, plot.hist = FALSE)
# hide the chisq table to give simpler output
occ_gof$chisq.table <- NULL
print(occ_gof)

# 5.3.2 Model selection ----

# dredge all possible combinations of the occupancy covariates
occ_dredge <- dredge(updateable(occ_model))

# model comparison to explore the results for occupancy
mc <- as.data.frame(occ_dredge) %>% 
  select(starts_with("psi(p"), df, AICc, delta, weight)
# shorten names for printing
names(mc) <- names(mc) %>% 
  str_extract("(?<=psi\\(pland_[0-9]{2}_)[a-z_]+") %>% 
  coalesce(names(mc))

# take a quick peak at the model selection table
mutate_all(mc, ~ round(., 3)) %>% 
  head(18) %>% 
  knitr::kable()
