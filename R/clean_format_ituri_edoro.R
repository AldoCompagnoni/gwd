# ituri_edoro data preparation script
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
ituri_edoro_means <- dplyr::filter( kam_means, site == "ituri_edoro")
ituri_edoro_medians <- dplyr::filter( kam_medians, site == "ituri_edoro")

# Prepare site table -----------------------------

# Do this by hand
site_out <- data.frame( Site_name = 'Ituri (Edoro)',
                        Latitude  = 1.436800000000,
                        Longitude = 28.582600000000 )

# Check that the location is sensible sense
leaflet( data = site_out) %>% 
  addTiles() %>% 
  addCircleMarkers(~Longitude, ~Latitude)

# Store site information
write.csv( site_out, 'results/ituri_edoro_site.csv',
           row.names = F )

# Prepare taxonomic table ------------------------

# Produce the binomial used for checking
taxa_df         <- dplyr::select( ituri_edoro_means, latin) %>%
  rename( Submitted_Name = latin )

# Separate NAs
taxa_genus_df <- subset( taxa_df,  is.na( Submitted_Name ) )
taxa_df       <- subset( taxa_df, !is.na( Submitted_Name ) )

# Create function: get "cleaned" names
get_clean_names   <- function( nam, fuzzy = 0.1 ) lcvp_search( nam, max.distance = fuzzy )

# Clean names from the Leipzig's list of plants
clean_l         <- lapply( taxa_df$Submitted_Name , get_clean_names )
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
# 4 plausible typos which remain in clean data frame

# Check species without matches
no_match_v <- setdiff( taxa_df$Submitted_Name, 
                       clean_df$Submitted_Name )

# Rerun Leipzig list with fuzzy matching
reclean_l       <- lapply( no_match_v, lcvp_fuzzy_search )
reclean_df      <- reclean_l %>% bind_rows
# visually select species identified with lcvp_fuzzy_search
# Only genus specified, 2 taxa remain unresolved

# Final taxonomy files 
taxa_nofuzzy    <- clean_df
taxa_fuzzy      <- reclean_df
# Combine 
taxa_out        <- bind_rows( taxa_nofuzzy, taxa_fuzzy )
# Do "taxa unresolved" by hand (taxa with no matches found)
taxa_unresvd    <- data.frame( Submitted_Name = no_match_v,
                               site           = 'ituri_edoro' )

# store resolved AND unresolved taxa
write.csv( taxa_out, 'results/ituri_edoro_taxa.csv',
           row.names = F )
write.csv( taxa_unresvd, 'results/ituri_edoro_taxa_unresvd.csv',
           row.names = F )
