# BCI data preparation script
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
library(ggplot2)

# Read in and filter data ------------------------

kam_means <- read.table( 'data/Demographies_with_dbh_means_Kambach.txt', header = TRUE )
kam_medians <- read.table( 'data/Demographies_with_dbh_medians_Kambach.txt', header = TRUE )
bci_means <- dplyr::filter( kam_means, site == "bci")
bci_medians <- dplyr::filter( kam_medians, site == "bci")

# Prepare site table -----------------------------

# Do this by hand
site_out <- data.frame( Site_name = 'BCI',
                        Latitude  = 9.154300000000,
                        Longitude = -79.846100000000 )

# Check that the location is sensible sense
leaflet( data = site_out) %>% 
  addTiles() %>% 
  addCircleMarkers(~Longitude, ~Latitude)

# Store site information
write.csv( site_out, 'results/bci_site.csv',
           row.names = F )

# Prepare taxonomic table ------------------------

# Produce the binomial used for checking
taxa_df_bci         <- dplyr::select( bci_means, latin) %>%
  rename( Submitted_Name = latin )

# Create function: get "cleaned" names
get_clean_names <- function( nam, fuzzy = 0 ) lcvp_search( nam, max.distance = fuzzy )

# Clean names from the Leipzig's list of plants
clean_l_bci         <- lapply( taxa_df_bci$Submitted_Name , get_clean_names )
clean_df_bci        <- clean_l_bci %>% bind_rows
clean_df_bci

# check whether there are issues
notfound_df_bci     <- clean_df_bci %>% subset( PL.comparison != 'identical' )

# Rerun Leipzig list with fuzzy matching
reclean_l_bci       <- lapply( notfound_df_bci$Submitted_Name, get_clean_names, 5 )
reclean_df_bci      <- reclean_l_bci %>% bind_rows

# Examine final file
reclean_df_bci

# Prepare demographic table ----------------------

# Check for missing/NULL/Inf values?
# Check for sample size 0
ggplot(data = bci_means, aes(site, growth_layer1_mean)) + geom_boxplot()
#find way to repeat for all parameters

# [Make boxplots to explore the data]

