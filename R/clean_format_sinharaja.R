# sinharaja data preparation script
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
sinharaja_means <- dplyr::filter( kam_means, site == "sinharaja")
sinharaja_medians <- dplyr::filter( kam_medians, site == "sinharaja")

# Prepare site table -----------------------------

# Do this by hand
site_out <- data.frame( Site_name = 'Sinharaja',
                        Latitude  = 6.402300000000,
                        Longitude = 80.402300000000 )

# Check that the location is sensible sense
leaflet( data = site_out) %>% 
  addTiles() %>% 
  addCircleMarkers(~Longitude, ~Latitude)

# Store site information
write.csv( site_out, 'results/sinharaja_site.csv',
           row.names = F )
