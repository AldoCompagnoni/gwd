# HKK data preparation script
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
hkk_means <- dplyr::filter( kam_means, site == "hkk")
hkk_medians <- dplyr::filter( kam_medians, site == "hkk")

# Prepare site table -----------------------------

# Do this by hand
site_out <- data.frame( Site_name = 'HKK',
                        Latitude  = 15.632400000000,
                        Longitude = 99.217000000000 )

# Check that the location is sensible sense
leaflet( data = site_out) %>% 
  addTiles() %>% 
  addCircleMarkers(~Longitude, ~Latitude)

# Store site information
write.csv( site_out, 'results/hkk_site.csv',
           row.names = F )

# Prepare taxonomic table ------------------------

# Produce the binomial used for checking
taxa_df         <- dplyr::select( hkk_means, latin ) %>%
  rename( Submitted_Name = latin )

# Separate NAs - genus names only are unresolved
taxa_genus_df <- subset( taxa_df,  is.na( Submitted_Name ) )
taxa_df       <- subset( taxa_df, !is.na( Submitted_Name ) )

# Create function: get "cleaned" names
get_clean_names   <- function( nam, fuzzy = 0.1 ) lcvp_search( nam, max.distance = fuzzy )

# Clean names from the Leipzig's list of plants
clean_l         <- lapply( taxa_df$Submitted_Name, get_clean_names )
clean_df        <- clean_l %>% 
  bind_rows %>% 
  rename( Submitted_Name      = Search,
          First_matched_Name  = Input.Taxon,
          LCVP_Accepted_Taxon = Output.Taxon ) %>% 
  # check for accepted names
  mutate( mismatch_test = str_detect( First_matched_Name, 
                                      Submitted_Name ) )

# visual check of mismatches and alternative spellings
mismatch_df <- clean_df %>% subset( !mismatch_test )
check_mismatches <- lapply( mismatch_df$Submitted_Name, lcvp_fuzzy_search )
# 7 plausible typos which remain in clean data frame

# remove 5 unresolved species containing "cf." and "sp.1" from clean data frame
mismatch_unresvd <- data.frame("Submitted_Name" = grep( 'sp.1|cf.', mismatch_df$Submitted_Name, value = T ))
clean_df_final <- data.frame("Submitted_Name" = grep( 'sp.1|cf.', clean_df$Submitted_Name, value = T, invert = T ))

# Check species without matches
no_match_v <- data.frame( "Submitted_Name" = setdiff( taxa_df$Submitted_Name, 
                       clean_df$Submitted_Name ))

# Rerun Leipzig list with fuzzy matching
reclean_l       <- lapply( no_match_v, lcvp_fuzzy_search )
reclean_df      <- reclean_l %>% bind_rows
# visually select species identified with lcvp_fuzzy_search
# No fuzzy match for 3 species, fuzzy matches for 2 taxa are not reliable (only identified to genus level)

# Final taxonomy files 
# Clean taxa should have LCVP search results
taxa_out        <- lapply( clean_df_final$Submitted_Name, get_clean_names ) %>% 
  bind_rows %>% 
  rename( Submitted_Name      = Search,
          First_matched_Name  = Input.Taxon,
          LCVP_Accepted_Taxon = Output.Taxon ) %>% 
  # check for accepted names
  mutate( mismatch_test = str_detect( First_matched_Name, 
                                      Submitted_Name ) )

# Do "taxa unresolved" by hand (taxa with no matches found)
taxa_unresvd    <- bind_rows( taxa_genus_df, mismatch_unresvd, no_match_v ) %>%
  mutate( site = 'hkk' )
  

# store resolved AND unresolved taxa
write.csv( taxa_out, 'results/hkk_taxa.csv',
           row.names = F )
write.csv( taxa_unresvd, 'results/hkk_taxa_unresvd.csv',
           row.names = F )
