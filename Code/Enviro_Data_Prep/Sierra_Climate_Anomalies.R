## -------------------------------------------------------------
##
## Script name: Anomaly Data Formatting
##
## Script purpose: Formatting anomaly data from GEE for DCM
##
## Author: Spencer R Keyser
##
## Date Created: 2025-08-06
##
## Email: srk252@cornell.edu
##
## Github: https://github.com/skeyser
##
## -------------------------------------------------------------
##
## Notes:
##
##
## -------------------------------------------------------------

## Defaults
options(scipen = 10, digits = 10)

## -------------------------------------------------------------

## Package Loading
library(dplyr)
library(stringr)
library(ggplot2)
library(here)
library(terra)
library(sf)
library(exactextractr)

## -------------------------------------------------------------

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Data Loading and organizing
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Long-term baseline climate
btmax <- terra::rast("Data/Spatial_Data/Climate_Anomaly_Rasters/Tmax_Baseline_1980_2010.tif")
btmax.ann <- terra::rast("Data/Spatial_Data/Climate_Anomaly_Rasters/Annual_Tmax_Baseline_1980_2010.tif")
btmin <- terra::rast("Data/Spatial_Data/Climate_Anomaly_Rasters/Tmin_Baseline_1980_2010.tif")
tmin.ann <- terra::rast("Data/Spatial_Data/Climate_Anomaly_Rasters/Annual_Tmin_Baseline_1980_2010.tif")
bpcp <- terra::rast("Data/Spatial_Data/Climate_Anomaly_Rasters/Precip_Baseline_1980_2010.tif")
bpcp.ann <- terra::rast("Data/Spatial_Data/Climate_Anomaly_Rasters/Annual_Prcp_Baseline_1980_2010.tif")



bclim <- c(btmax, btmin, bpcp, btmax.ann, tmin.ann, bpcp.ann)
names(bclim) <- c("Tmax_Baseline_JJA", "Tmin_Baseline_JJA", "Prcp_Baseline_JJA", "Tmax_Baseline_Annual", "Tmin_Baseline_Annual", "Prcp_Baseline_Annual")
plot(bclim)

## Order the .tifs for the anomalies by date
atmax.path <- sort(list.files(here("Data/Spatial_Data/Climate_Anomaly_Rasters/"), pattern = "^Tmax_Anomaly*", full.names = T))
atmin.path <- sort(list.files(here("Data/Spatial_Data/Climate_Anomaly_Rasters/"), pattern = "^Tmin_Anomaly*", full.names = T))
aprcp.path <- sort(list.files(here("Data/Spatial_Data/Climate_Anomaly_Rasters/"), pattern = "^Precip_Anomaly*", full.names = T))
antmax.path <- sort(list.files(here("Data/Spatial_Data/Climate_Anomaly_Rasters/"), pattern = "Annual_Tmax_Anomaly*", full.names = T))
antmin.path <- sort(list.files(here("Data/Spatial_Data/Climate_Anomaly_Rasters/"), pattern = "Annual_Tmin_Anomaly*", full.names = T))
anprcp.path <- sort(list.files(here("Data/Spatial_Data/Climate_Anomaly_Rasters/"), pattern = "Annual_Prcp_Anomaly*", full.names = T))

## Trends
trtmax <- terra::rast("Data/Spatial_Data/Climate_Trend_Rasters/Tmax_Trend_Slope_1980_2020.tif")
trtmaxmam <- terra::rast("Data/Spatial_Data/Climate_Trend_Rasters/Spring_Tmax_Trend_Slope_1980_2020.tif")
trtmaxjja <- terra::rast("Data/Spatial_Data/Climate_Trend_Rasters/Summer_Tmax_Trend_Slope_1980_2020.tif")
trtmin <- terra::rast("Data/Spatial_Data/Climate_Trend_Rasters/Tmin_Trend_Slope_1980_2020.tif")
trtminmam <- terra::rast("Data/Spatial_Data/Climate_Trend_Rasters/Spring_Tmax_Trend_Slope_1980_2020.tif")
trtminjja <- terra::rast("Data/Spatial_Data/Climate_Trend_Rasters/Summer_Tmax_Trend_Slope_1980_2020.tif")
trprcp <- terra::rast("Data/Spatial_Data/Climate_Trend_Rasters/Ppt_Trend_Slope_1980_2020.tif")
trprcpmam <- terra::rast("Data/Spatial_Data/Climate_Trend_Rasters/Spring_Ppt_Trend_Slope_1980_2020.tif")
trprcpjja <- terra::rast("Data/Spatial_Data/Climate_Trend_Rasters/Summer_Ppt_Trend_Slope_1980_2020.tif")

trends <- c(trtmax, trtmaxmam, trtmaxjja,
            trtmin, trtminmam, trtminjja,
            trprcp, trprcpmam, trprcpjja)
names(trends) <- c("Trend_Tmax", "Trend_Tmax_MAM", "Trend_Tmax_JJA", 
                   "Trend_Tmin", "Trend_Tmin_MAM", "Trend_Tmin_JJA",
                   "Trend_Prcp", "Trend_Prcp_MAM", "Trend_Prcp_JJA")
plot(trends)

## read in the tifs
tanom.r <- rast(lapply(atmax.path, rast))
manom.r <- rast(lapply(atmin.path, rast))
panom.r <- rast(lapply(aprcp.path, rast))
antanom.r <- rast(lapply(antmax.path, rast))
anmanom.r <- rast(lapply(antmin.path, rast))
anpanom.r <- rast(lapply(anprcp.path, rast))

anom.r <- c(tanom.r, manom.r, panom.r, antanom.r, anmanom.r, anpanom.r)

## Give informative names
names(anom.r) <- stringr::str_extract(sources(anom.r), pattern = "Tmax_Anomaly_\\d{4}|Tmin_Anomaly_\\d{4}|Precip_Anomaly_\\d{4}|Annual_Tmax_Anomaly_\\d{4}|Annual_Tmin_Anomaly_\\d{4}|Annual_Prcp_Anomaly_\\d{4}")

## Load in the ARU locations
locs <- st_read(here("Data/Spatial_Data/ARU_Locs_2021_2025.shp"))
glimpse(locs)
head(unique(sort(locs$Cll_Unt)))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Buffering points and extraction 
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## generate the buffer - starting with 120 but can set variable buffers if needed
locs.buff <- locs |> st_buffer(dist = units::set_units(500, "m"))

## terra extraction
anom.ex <- exact_extract(anom.r, locs.buff, fun = "mean")
base.ex <- exact_extract(bclim, locs.buff, fun = "mean")
trend.ex <- exact_extract(trends, locs.buff, fun = "mean")

## Bind with the data
locs <- cbind(locs, anom.ex, base.ex, trend.ex)

## Clean up the names for the columns
colnames(locs) <- gsub("mean.", "", colnames(locs))
colnames(locs) <- gsub("Precip", "Prcp", colnames(locs))

glimpse(locs)
length(unique(locs$Cll_Unt))/nrow(locs)

locs <- locs |> 
  arrange(Cll_Unt)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Organizing and formatting for DCM
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Create a data sets for static predictors
static_clim_df <- locs |> 
  select(Cell_Unit = Cll_Unt,
         Long, Lat,
         contains("Trend"),
         contains("Baseline")) |>
  st_drop_geometry() |> 
  group_by(Cell_Unit) |> 
  summarise(Long = mean(Long),
            Lat = mean(Lat),
            Trend_Tmax = mean(Trend_Tmax),
            Trend_Tmin = mean(Trend_Tmin),
            Trend_Prcp = mean(Trend_Prcp),
            Tmax_Base_JJA = mean(Tmax_Baseline_JJA),
            Tmax_Base_Annual = mean(Tmax_Baseline_Annual),
            Tmin_Base_JJA = mean(Tmin_Baseline_JJA),
            Tmin_Base_Annual = mean(Tmin_Baseline_Annual),
            Prcp_Base_JJA = mean(Prcp_Baseline_JJA),
            Prcp_Base_Annual = mean(Prcp_Baseline_Annual)) |> 
  ungroup() |> 
  mutate(across(where(is.numeric), 
                ~as.numeric(scale(.)),
                .names = "{.col}_scaled"))

  

hist(static_clim_df$Trend_Tmax)
hist(static_clim_df$Trend_Tmin)
hist(static_clim_df$Trend_Prcp)
hist(static_clim_df$Tmax_Base_JJA)
hist(static_clim_df$Tmin_Base_JJA)
hist(static_clim_df$Tmax_Base_Annual)
hist(static_clim_df$Tmin_Base_Annual)
hist(static_clim_df$Prcp_Base_JJA)
hist(static_clim_df$Prcp_Base_Annual)

## Dynamic climate dataset
## Locs in long format
dyn_clim_df <- locs |>
  select(id:Lat,
         contains("Anomaly")) |> 
  tidyr::pivot_longer(cols = matches("Tmax_Anomaly|Tmin_Anomaly|Prcp_Anomaly|Precip_Anomaly"), names_to = "Var_Year", values_to = "Anomaly") |>
  mutate(Year = as.numeric(stringr::str_extract(Var_Year, "\\d+")),
         ClimVar = stringr::str_extract(Var_Year, "[A-z]+")) |> 
  mutate(ClimVar = gsub("_Anomaly_", "", ClimVar)) |> 
  filter(Year >= 2021) |> 
  st_drop_geometry()

## Make a list to store the datasets
cvars <- unique(dyn_clim_df$ClimVar)
dyn_clim <- vector(mode = "list", length = length(cvars))

for(i in 1:length(cvars)){
  cvar.tmp <- cvars[i]
  tmp <- dyn_clim_df |> 
    select(Cell_Unit = Cll_Unt,
           ClimVar,
           Year,
           Anomaly) |> 
    filter(ClimVar == cvar.tmp) |> 
    tidyr::pivot_wider(names_from = Year,
                       values_from = Anomaly,
                       values_fn = mean) |> 
    select(-ClimVar) |> 
    tibble::column_to_rownames("Cell_Unit") |> 
    as.matrix()
  
  dyn_clim[[i]] <- tmp
  names(dyn_clim)[i] <- cvar.tmp
  
  
}

str(dyn_clim)
str(static_clim_df)


## Write to csv
# locs_long <- locs_long |> st_drop_geometry()
# data.table::fwrite(locs_long, file = here("Data/Tmax_Anomaly_ARU_21_24.csv"))

## Write some model objects
static_clim_df <- as.data.frame(static_clim_df)

clim_dat <- list(static_clim_df, dyn_clim)
str(clim_dat)

saveRDS(clim_dat, file = here("Data/Climate_DCM_Covs_2021_2024.RDS"))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: carpentry and plotting with full extracted data
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ## Check out the extractions
# glimpse(locs)
# 
# ## Look at the overall trend in anomaly
# anom |> 
#   group_by(anomaly_year) |> 
#   summarise(anom_mn = mean(anomaly_value)) |> 
#   ggplot(aes(x = anomaly_year, y = anom_mn)) +
#   geom_point() + 
#   geom_line() + 
#   theme_bw() + 
#   xlab("Year") + 
#   ylab("Mean T. Anomaly from 1980-2010 (Deg. C)")
# 
# ## Original figure year x anom x lat
# anom_plot_points <- locs_long |> 
#   ggplot(aes(x = Year, y = Anomaly, color = Lat)) +
#   geom_point(position = "jitter") + 
#   geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
#   scale_color_viridis_c(option = "C", direction = -1) +
#   scale_x_continuous(breaks = c(2011, 2013, 2015, 2017, 2019, 2021, 2023)) +
#   xlab("Year") + 
#   ylab("Max T. Anomaly from 1980-2010 (Deg. C)") + 
#   labs(color = "Latitude") + 
#   theme_bw()
# 
# ## Some alpha coding
# anom_plot_pointsmooth <- locs_long |> 
#   ggplot(aes(x = Year, y = Anomaly, color = Lat)) +
#   geom_point(alpha = 0.3, position = position_jitter(width = 0.2)) + 
#   geom_smooth(aes(group = 1), color = "black", se = TRUE) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
#   scale_color_viridis_c(option = "D", direction = -1) +
#   scale_x_continuous(breaks = seq(2011, 2023, by = 2)) +
#   xlab("Year") + 
#   ylab("Max T. Anomaly from 1980-2010 (Deg. C)") + 
#   labs(color = "Latitude") + 
#   theme_bw()
# 
# ## Smoothing
# anom_plot_pointsmooth <- locs_long |> 
#   mutate(Lat_category = case_when(
#            Lat >= 40 ~ "High",
#            Lat >= 38 & Lat < 40 ~ "Mid",
#            Lat < 38 ~ "Low"
#          ),
#          Lat_category = factor(Lat_category, levels = c("Low", "Mid", "High"))) |> 
#   ggplot(aes(x = Year, y = Anomaly, color = Lat_category)) +
#   geom_point(alpha = 0.3, position = position_jitter(width = 0.2)) + 
#   geom_smooth(aes(group = Lat_category), method = "loess", se = F) +  # Changed this line
#   geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
#   scale_color_manual(values = c("#FF9E44", "#1E88E5", "#004D40")) +
#   scale_x_continuous(breaks = seq(2011, 2023, by = 2)) +
#   xlab("Year") + 
#   ylab("Max T. Anomaly from 1980-2010 (Deg. C)") + 
#   labs(color = "Latitude Category") + 
#   theme_bw()
# 
# ## Plot anomaly example data
# ## Add in spatial layer for the Sierra and California
# sierra_roi <- st_read(here("Data/Spatial_Data/Sierra_ROI.shp")) |> 
#   st_transform(crs = "EPSG:3310")
# 
# ## Read in a few anomalies
# anom.file <- list.files(here("Data/Spatial_Data/Climate_Anomaly_Rasters/"), pattern = "Tmax_Anomaly", full.names = T)
# 
# 
# anom <- rast(here("Data/Spatial_Data/Climate_Anomaly_Rasters/Tmax_Mean_Anomaly_2011_2023.tif"))
# anom <- project(anom, st_crs(sierra_roi)$wkt)
# anom <- mask(anom, ca)
# plot(anom)
# 
# anom.df <- as.data.frame(anom, xy = T)
# colnames(anom.df)
# 
# anom.rplot <- ggplot() +
#   geom_sf(data = ca, color = "black", fill = "gray", lwd = 1) +
#   geom_raster(data = anom.df, aes(x=x, y=y, fill=tmax)) +
#   #geom_sf(data = sierra_roi, color = "gray", fill = NA, lwd = 0.5) +
#   scale_fill_viridis_c(option = "F",
#                        guide = guide_colorbar(direction = "horizontal",
#                                               title.position = "top",
#                                               title.hjust = 0.5,
#                                               barwidth = 10)) +
#   theme_void() + 
#   labs(fill = expression(Delta~"Mean Max Temp"~degree*"C from (1980-2010)")) + 
#   theme(legend.position = "bottom",
#         legend.text = element_text(family = "serif",
#                                    size = 12))
# 
# # First get the overall min and max values across all files
# get_raster_extremes <- function(file_paths) {
#   # Initialize min and max
#   overall_min <- Inf
#   overall_max <- -Inf
#   
#   # Loop through each file
#   for(file in file_paths) {
#     r <- rast(file)
#     r <- project(r, st_crs(sierra_roi)$wkt)
#     r <- mask(r, ca)
#     
#     # Update min and max
#     overall_min <- min(overall_min, minmax(r)[1], na.rm = TRUE)
#     overall_max <- max(overall_max, minmax(r)[2], na.rm = TRUE)
#   }
#   
#   return(c(overall_min, overall_max))
# }
# 
# # Create plotting function with standardized limits
# create_anomaly_plot <- function(file_path, limits, 
#                                 viridis_option = "F",
#                                 base_fill = "gray",
#                                 text_size = 14,
#                                 title_size = 14) {
#   # Read and process raster
#   anom <- rast(file_path)
#   anom <- project(anom, st_crs(sierra_roi)$wkt)
#   anom <- mask(anom, ca)
#   
#   # Convert to dataframe
#   anom.df <- as.data.frame(anom, xy = T)
#   
#   # Extract year from filename
#   year <- str_extract(basename(file_path), "\\d{4}")
#   
#   # Create plot
#   p <- ggplot() +
#     geom_sf(data = ca, color = "black", fill = base_fill, lwd = 1) +
#     geom_raster(data = anom.df, aes(x=x, y=y, fill=tmax)) +
#     scale_fill_gradient2(low = "blue",
#                          mid = "white",
#                          high = "red",
#                          midpoint = 0,
#                          limits = limits,  # Set standardized limits
#                          guide = guide_colorbar(direction = "horizontal",
#                                                 title.position = "top",
#                                                 title.hjust = 0.5,
#                                                 barwidth = 15)) +
#     theme_void() +
#     labs(fill = expression(Delta~"Summer Max Temp"~degree*"C from (1980-2010)"),
#          title = paste("Temperature Anomaly", year)) + 
#     theme(legend.position = "bottom",
#           legend.text = element_text(family = "serif", size = text_size),
#           legend.title = element_text(family = "serif", size = text_size),
#           plot.title = element_text(hjust = 0.5, family = "serif", size = title_size))
#   
#   return(p)
# }
# 
# # Get all files
# anom.files <- list.files(here("Data/Spatial_Data/Climate_Anomaly_Rasters/"), 
#                          pattern = "Tmax_Anomaly", 
#                          full.names = TRUE)
# 
# # Get overall min and max
# limits <- get_raster_extremes(anom.files)
# 
# # Create list of plots with standardized limits
# plot_list <- lapply(anom.files, function(x) {
#   create_anomaly_plot(x, 
#                       limits = limits,
#                       viridis_option = "F",
#                       base_fill = "lightgray",
#                       text_size = 20,
#                       title_size = 20)
# })
# 
# # Create gif
# gifski::save_gif(
#   expr = {
#     for(i in seq_along(plot_list)) {
#       print(plot_list[[i]])
#     }
#   },
#   gif_file = here("Figures/temperature_anomalies.gif"),
#   width = 1200,
#   height = 1200,
#   delay = 1.5,
#   loop = TRUE
# )
