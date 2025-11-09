# Abundance_estimations from:
# https://cornelllabofornithology.github.io/ebird-best-practices/abundance.html

# 6.2 Data Prep ----

library(lubridate)
library(sf)
library(raster)
library(dggridR)
library(pdp)
library(mgcv)
library(fitdistrplus)
library(viridis)
library(fields)
library(tidyverse)
# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

# set random number seed to insure fully repeatable results
set.seed(1)

# setup output directory for saved results
if (!dir.exists("output")) {
  dir.create("output")
}

# ebird data
ebird <- read_csv("data/ebd_woothr_june_bcr27_zf.csv") %>% 
  mutate(protocol_type = factor(protocol_type, 
                                levels = c("Stationary" , "Traveling"))) %>%
  # remove observations with no count
  filter(!is.na(observation_count))

# modis habitat covariates
habitat <- read_csv("data/pland-elev_location-year.csv") %>% 
  mutate(year = as.integer(year))

# combine ebird and habitat data
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

# 6.2.1 Spatiotemporal subsampling ----

# generate hexagonal grid with ~ 5 km betweeen cells
dggs <- dgconstruct(spacing = 5)
# get hexagonal cell id and week number for each checklist
checklist_cell <- ebird_habitat %>% 
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
         week = week(observation_date))
# sample one checklist per grid cell per week
# sample detection/non-detection independently 
ebird_ss <- checklist_cell %>% 
  group_by(species_observed, year, week, cell) %>% 
  sample_n(size = 1) %>% 
  ungroup() %>% 
  select(-cell, -week)

# 6.2.2 Test train split ----
hab_covs <- c("pland_04_deciduous_broadleaf", 
              "pland_05_mixed_forest",
              "pland_12_cropland",
              "pland_13_urban")
ebird_split <- ebird_ss %>% 
  # select only the columns to be used in the model
  select(observation_count,
         # effort covariates
         day_of_year, time_observations_started, duration_minutes,
         effort_distance_km, number_observers, protocol_type,
         # habitat covariates
         hab_covs) 
# split 80/20
ebird_split <- ebird_split %>% 
  split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))
map_int(ebird_split, nrow)
#>  test train 
#>  3941 16145


# 6.3 Exploratory data analysis ----
p <- par(mfrow = c(1, 2))
# counts with zeros
hist(ebird_ss$observation_count, main = "Histogram of counts", 
     xlab = "Observed count")
# counts without zeros
pos_counts <- keep(ebird_ss$observation_count, ~ . > 0)
hist(pos_counts, main = "Histogram of counts > 0", 
     xlab = "Observed non-zero count")
par(p)

# 6.4 Abundance models ----

# gam parameters
# degrees of freedom for smoothing
k <- 5
# degrees of freedom for cyclic time of day smooth
k_time <- 7 

# continuous predictors
# hold out time to treat seperately since it's cyclic
continuous_covs <- ebird_split$train %>% 
  select(-observation_count, -protocol_type, -time_observations_started) %>% 
  names()

# create model formula for predictors
gam_formula_rhs <- str_glue("s({var}, k = {k})", 
                            var = continuous_covs, k = k) %>% 
  str_flatten(collapse = " + ") %>% 
  str_glue(" ~ ", .,
           " + protocol_type + ",
           "s(time_observations_started, bs = \"cc\", k = {k})", 
           k = k_time) %>% 
  as.formula()

# model formula including response
gam_formula <- update.formula(observation_count ~ ., gam_formula_rhs)
gam_formula
#> observation_count ~ s(day_of_year, k = 5) + s(duration_minutes, 
#>     k = 5) + s(effort_distance_km, k = 5) + s(number_observers, 
#>     k = 5) + s(pland_04_deciduous_broadleaf, k = 5) + s(pland_05_mixed_forest, 
#>     k = 5) + s(pland_12_cropland, k = 5) + s(pland_13_urban, 
#>     k = 5) + protocol_type + s(time_observations_started, bs = "cc", 
#>     k = 7)

# explicitly specify where the knots should occur for time_observations_started
# this ensures that the cyclic spline joins the variable at midnight
# this won't happen by default if there are no data near midnight
time_knots <- list(time_observations_started = seq(0, 24, length.out = k_time))

# zero-inflated poisson
m_ziplss <- gam(list(gam_formula,      # count model
                     gam_formula[-2]), # presence model
                data = ebird_split$train, 
                family = "ziplss", 
                knots = time_knots)

# negative binomial
m_nb <- gam(gam_formula,
            data = ebird_split$train, 
            family = "nb",
            knots = time_knots)

# tweedie distribution
m_tw <- gam(gam_formula,
            data = ebird_split$train, 
            family = "tw",
            knots = time_knots)


# 6.5 Assessmwent ----

obs_count <- select(ebird_split$test, obs = observation_count)

# presence probability is on the complimentary log-log scale
# we can get the inverse link function with
inv_link <- binomial(link = "cloglog")$linkinv
# combine ziplss presence and count predictions
m_ziplss_pred <- predict(m_ziplss, ebird_split$test, type = "link") %>% 
  as.data.frame() %>%
  transmute(family = "Zero-inflated Poisson",
            pred = inv_link(V2) * exp(V1)) %>% 
  bind_cols(obs_count)

m_nb_pred <- predict(m_nb, ebird_split$test, type = "response") %>% 
  tibble(family = "Negative Binomial", pred = .) %>% 
  bind_cols(obs_count)

m_tw_pred <- predict(m_tw, ebird_split$test, type = "response") %>% 
  tibble(family = "Tweedie", pred = .) %>% 
  bind_cols(obs_count)

# combine predictions from all three models
test_pred <- bind_rows(m_ziplss_pred, m_nb_pred, m_tw_pred) %>% 
  mutate(family = as_factor(family))

# 6.5.1 Predictive performance: ranking ----

# spearman’s rank correlation
test_pred %>% 
  group_by(family) %>% 
  summarise(rank_cor = cor.test(obs, pred, 
                                method = "spearman", 
                                exact = FALSE)$estimate) %>% 
  ungroup()

# 6.5.2 Predictive performance: Magnatude ----

#  plot predicted vs. observed
ticks <- c(0, 1, 10, 100, 1000)
mx <- round(max(test_pred$obs))
ggplot(test_pred) +
  aes(x = log10(obs + 1), 
      y = log10(pred + 1)) +
  geom_jitter(alpha = 0.2, height = 0) +
  # y = x line
  geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
  # area where counts off by a factor of 10
  geom_area(data = tibble(x = log10(seq(0, mx - 1) + 1), 
                          y = log10(seq(0, mx - 1) / 10 + 1)),
            mapping = aes(x = x, y = y),
            fill = "red", alpha = 0.2) +
  # loess fit
  geom_smooth(method = "loess", 
              method.args = list(span = 2 / 3, degree = 1)) +
  scale_x_continuous(breaks = log10(ticks + 1), labels = ticks) +
  scale_y_continuous(breaks = log10(ticks + 1), labels = ticks) +
  labs(x = "Observed count",
       y = "Predicted count") +
  facet_wrap(~ family, nrow = 1)

# How many underestimates by an order of magnitude? 
test_pred %>% 
  group_by(family) %>% 
  summarize(n = sum(obs / pred > 10),
            pct = mean(obs / pred > 10))

# mean absolute deviation
test_pred %>% 
  group_by(family) %>% 
  summarise(mad = mean(abs(obs - pred), na.rm = TRUE)) %>% 
  ungroup()


# 6.5.3 Assessing covariates ----

# ggplot function
plot_gam <- function(m, title = NULL, ziplss = c("presence", "abundance")) {
  # capture plot
  tmp <- tempfile()
  png(tmp)
  p <- plot(m, pages = 1)
  dev.off()
  unlink(tmp)
  
  # drop addition models in ziplss
  if (m$family$family == "ziplss") {
    is_presence <- map_lgl(p, ~ str_detect(.$ylab, "^s\\.1"))
    if (ziplss == "presence") {
      p <- p[is_presence]  
    } else {
      p <- p[!is_presence]
    }
  }
  
  # extract data
  p_df <- map_df(p, ~ tibble(cov = rep(.$xlab, length(.$x)),
                             x = .$x, fit = .$fit, se = .$se))
  
  # plot
  g <- ggplot(p_df) +
    aes(x = x, y = fit,
        ymin = fit - se, ymax = fit + se) +
    geom_ribbon(fill = "grey80") +
    geom_line(col = "blue") +
    facet_wrap(~ cov, scales = "free_x") +
    labs(x = NULL,
         y = "Smooth function",
         title = title)
  print(g)
  invisible(p_df)
}

plot_gam(m_nb, title = "Negative Binomial GAM")
#plot_gam(m_zipliss, title = "Zero-inflated Poisson GAM")
plot_gam(m_tw, title = "Tweedie GAM")


# set negative binomial model to be used for predictions
pred_model <- m_nb

# 6.6Prediction ----

# create a dataframe of covariates with a range of start times
seq_tod <- seq(0, 24, length.out = 300)
tod_df <- ebird_split$train %>% 
  # find average pland habitat covariates
  select(starts_with("pland")) %>% 
  summarize_all(mean, na.rm = TRUE) %>% 
  ungroup() %>% 
  # use standard checklist
  mutate(day_of_year = yday(ymd(str_glue("{max_lc_year}-06-15"))),
         duration_minutes = 60,
         effort_distance_km = 1,
         number_observers = 1,
         protocol_type = "Traveling") %>% 
  cbind(time_observations_started = seq_tod)

# predict at different start times
pred_tod <- predict(pred_model, newdata = tod_df, 
                    type = "link", 
                    se.fit = TRUE) %>% 
  as_tibble() %>% 
  # calculate backtransformed confidence limits
  transmute(time_observations_started = seq_tod,
            pred = pred_model$family$linkinv(fit),
            pred_lcl = pred_model$family$linkinv(fit - 1.96 * se.fit),
            pred_ucl = pred_model$family$linkinv(fit + 1.96 * se.fit))

# find optimal time of day
t_peak <- pred_tod$time_observations_started[which.max(pred_tod$pred_lcl)]

# plot the partial dependence plot
ggplot(pred_tod) +
  aes(x = time_observations_started, y = pred,
      ymin = pred_lcl, ymax = pred_ucl) +
  geom_ribbon(fill = "grey80", alpha = 0.5) +
  geom_line() +
  geom_vline(xintercept = t_peak, color = "blue", linetype = "dashed") +
  labs(x = "Hours since midnight",
       y = "Predicted relative abundance",
       title = "Effect of observation start time on Wood Thrush reporting",
       subtitle = "Peak detectability shown as dashed blue line")


# add effort covariates to prediction surface
pred_surface_eff <- pred_surface %>% 
  mutate(day_of_year = yday(ymd(str_glue("{max_lc_year}-06-15"))),
         time_observations_started = t_peak,
         duration_minutes = 60,
         effort_distance_km = 1,
         number_observers = 1,
         protocol_type = "Traveling")

# predict
pred <- predict(pred_model, newdata = pred_surface_eff, 
                type = "link", 
                se.fit = TRUE) %>% 
  as_tibble() %>% 
  # calculate confidence limits and back transform
  transmute(abd = pred_model$family$linkinv(fit),
            abd_se = pred_model$family$linkinv(se.fit),
            abd_lcl = pred_model$family$linkinv(fit - 1.96 * se.fit),
            abd_ucl = pred_model$family$linkinv(fit + 1.96 * se.fit)) %>%
  # add to prediction surface
  bind_cols(pred_surface_eff, .) %>% 
  select(latitude, longitude, abd, abd_se, abd_lcl, abd_ucl)


#convert to spatial features with sf
r_pred <- pred %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  select(abd, abd_se) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize
  rasterize(r)
r_pred <- r_pred[[-1]]

# save the rasters
tif_dir <- "output"
if (!dir.exists(tif_dir)) {
  dir.create(tif_dir)
}
writeRaster(r_pred[["abd"]], 
            filename = file.path(tif_dir, "abundance-model_abd_woothr.tif"),
            overwrite = TRUE)
writeRaster(r_pred[["abd_se"]], 
            filename = file.path(tif_dir, "abundance-model_se_woothr.tif"), 
            overwrite = TRUE)

# any expected abundances below this threshold are set to zero
zero_threshold <- 0.05

# project predictions
r_pred_proj <- projectRaster(r_pred, crs = map_proj$proj4string, method = "ngb")

par(mfrow = c(2, 1))
for (nm in names(r_pred)) {
  r_plot <- r_pred_proj[[nm]]
  
  par(mar = c(3.5, 0.25, 0.25, 0.25))
  # set up plot area
  plot(bcr, col = NA, border = NA)
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
  
  # modified plasma palette
  plasma_rev <- rev(plasma(25, end = 0.9))
  gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
  pal <- c(gray_int(4)[2], plasma_rev)
  
  # abundance vs. se
  if (nm == "abd") {
    title <- "Wood Thrush Relative Abundance"
    # set very low values to zero
    r_plot[r_plot <= zero_threshold] <- NA
    # log transform
    r_plot <- log10(r_plot)
    # breaks and legend
    mx <- ceiling(100 * cellStats(r_plot, max)) / 100
    mn <- floor(100 * cellStats(r_plot, min)) / 100
    brks <- seq(mn, mx, length.out = length(pal) + 1)
    lbl_brks <- sort(c(-2:2, mn, mx))
    lbls <- round(10^lbl_brks, 2)
  } else {
    title <- "Wood Thrush Abundance Uncertainty (SE)"
    # breaks and legend
    mx <- ceiling(1000 * cellStats(r_plot, max)) / 1000
    mn <- floor(1000 * cellStats(r_plot, min)) / 1000
    brks <- seq(mn, mx, length.out = length(pal) + 1)
    lbl_brks <- seq(mn, mx, length.out = 5)
    lbls <- round(lbl_brks, 2)
  }
  
  # abundance
  plot(r_plot, 
       col = pal, breaks = brks, 
       maxpixels = ncell(r_plot),
       legend = FALSE, add = TRUE)
  
  # borders
  plot(bcr, border = "#000000", col = NA, lwd = 1, add = TRUE)
  plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
  plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
  box()
  
  # legend
  par(new = TRUE, mar = c(0, 0, 0, 0))
  image.plot(zlim = range(brks), legend.only = TRUE, col = pal,
             smallplot = c(0.25, 0.75, 0.06, 0.09),
             horizontal = TRUE,
             axis.args = list(at = lbl_brks, 
                              labels = lbls,
                              fg = "black", col.axis = "black",
                              cex.axis = 0.75, lwd.ticks = 0.5,
                              padj = -1.5),
             legend.args = list(text = title,
                                side = 3, col = "black",
                                cex = 1, line = 0))
}
