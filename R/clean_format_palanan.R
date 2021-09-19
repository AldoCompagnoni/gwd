# palanan data preparation script
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
palanan_means <- dplyr::filter( kam_means, site == "palanan")
palanan_medians <- dplyr::filter( kam_medians, site == "palanan")

# Prepare site table -----------------------------

# Do this by hand
site_out <- data.frame( Site_name = 'Palanan',
                        Latitude  = 17.040200000000,
                        Longitude = 122.388000000000 )

# Check that the location is sensible sense
leaflet( data = site_out) %>% 
  addTiles() %>% 
  addCircleMarkers(~Longitude, ~Latitude)

# Store site information
write.csv( site_out, 'results/palanan_site.csv',
           row.names = F )
