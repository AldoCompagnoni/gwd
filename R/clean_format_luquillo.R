# luquillo data preparation script
# Data: Kambach

# Clear environment ------------------------------

rm(list=ls())

# Load packages ----------------------------------
library(tidyverse)
# devtools::install_github("idiv-biodiversity/LCVP")
library(LCVP)
# devtools::install_github("idiv-biodiversity/lcvplants")
library(lcvplants)
library(leaflet)

# Read in and filter data ------------------------

kam_means <- read.table( 'data/Demographies_with_dbh_means_Kambach.txt', header = TRUE )
kam_medians <- read.table( 'data/Demographies_with_dbh_medians_Kambach.txt', header = TRUE )
luquillo_means <- dplyr::filter( kam_means, site == "luquillo")
luquillo_medians <- dplyr::filter( kam_medians, site == "luquillo")

# Prepare site table -----------------------------

# Do this by hand
site_out <- data.frame( Site_name = 'Luquillo',
                        Latitude  = 18.326200000000,
                        Longitude = -65.816000000000 )

# Check that the location is sensible sense
leaflet( data = site_out) %>% 
  addTiles() %>% 
  addCircleMarkers(~Longitude, ~Latitude)

# Store site information
write.csv( site_out, 'results/luquillo_site.csv',
           row.names = F )
