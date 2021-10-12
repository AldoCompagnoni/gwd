# Fushan data preparation script
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
fushan_means <- dplyr::filter( kam_means, site == "fushan")
fushan_medians <- dplyr::filter( kam_medians, site == "fushan")

# Prepare site table -----------------------------

# Do this by hand
site_out <- data.frame( Site_name = 'Fushan',
                        Latitude  = 24.761400000000,
                        Longitude = 121.555000000000 )

# Check that the location is sensible sense
leaflet( data = site_out) %>% 
  addTiles() %>% 
  addCircleMarkers(~Longitude, ~Latitude)

# Store site information
write.csv( site_out, 'results/fushan_site.csv',
           row.names = F )

# Prepare taxonomic table ------------------------

# Produce the binomial used for checking
taxa_df         <- dplyr::select( fushan_means, latin, sp, genus, family, IDlevel ) %>%
  rename( Submitted_Name = latin, Sp_Code = sp, Submitted_Genus = genus, Submitted_Family = family )

# Separate NAs - no NAs
taxa_na_df      <- subset( taxa_df,  is.na( Submitted_Name ) )
taxa_na_rm_df   <- subset( taxa_df, !is.na( Submitted_Name ) )

# Create function: get "cleaned" names
get_clean_names <- function( nam, fuzzy = 0.1 ) lcvp_search( nam, max.distance = fuzzy )

# Clean names from the Leipzig's list of plants
clean_l         <- lapply( taxa_na_rm_df$Submitted_Name , get_clean_names )
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
no_match_v <- data.frame( "Submitted_Name" = setdiff( taxa_na_rm_df$Submitted_Name, 
                       clean_df$Submitted_Name ))
# 0 species without matches

# Rerun Leipzig list with fuzzy matching
reclean_l       <- lapply( no_match_v, lcvp_fuzzy_search )
# visually select species identified with lcvp_fuzzy_search (only genus specified, all remain unresolved)
reclean_df      <- reclean_l %>% bind_rows

# Final taxonomy files 
taxa_nofuzzy    <- clean_df
taxa_fuzzy      <- reclean_df
# Combine non-fuzzy (including typos) and fuzzy matches
taxa_out        <- bind_rows( taxa_nofuzzy, taxa_fuzzy) %>%
                    mutate( site = 'fushan' )
# No unresolved taxa


# store resolved taxa
write.csv( taxa_out, 'results/fushan_taxa.csv',
           row.names = F )

