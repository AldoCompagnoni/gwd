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

# Prepare taxonomic table ------------------------

# Produce the binomial used for checking
taxa_df         <- dplyr::select( palanan_means, latin, sp, genus, family, IDlevel ) %>%
  rename( Submitted_Name = latin, Sp_Code = sp, Submitted_Genus = genus, Submitted_Family = family )

#Capitalise taxa so mismatch test works
upper_case_genus <- function( x ){
  
  x %>% 
    separate( Submitted_Name, c('gen','spp'), 
              sep = ' ', extra = 'merge' ) %>% 
    mutate( gen   = str_to_title(gen) ) %>% 
    mutate( Submitted_Name = paste0(gen,' ',spp) ) %>% 
    dplyr::select( -gen, -spp)
  
}
taxa_df <- upper_case_genus( taxa_df )

# Separate NAs - no NAs
taxa_na_df      <- subset( taxa_df,  is.na( Submitted_Name ) )
taxa_na_rm_df   <- subset( taxa_df, !is.na( Submitted_Name ) )

# Create function: get "cleaned" names
get_clean_names   <- function( nam, fuzzy = 0.1 ) lcvp_search( nam, max.distance = fuzzy )

# Clean names from the Leipzig's list of plants
clean_l         <- lapply( taxa_na_rm_df$Submitted_Name, get_clean_names )
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
# 24 plausible typos which remain in clean data frame

# remove 10 unresolved species containing "cf." from clean data frame
mismatch_unresvd <- data.frame("Submitted_Name" = grep( 'cf.', mismatch_df$Submitted_Name, value = T ))
clean_df_final <- data.frame("Submitted_Name" = grep( 'cf.', clean_df$Submitted_Name, value = T, invert = T )) 

# Check species without matches
no_match_v <- data.frame( "Submitted_Name" = setdiff( taxa_na_rm_df$Submitted_Name, 
                                                      clean_df$Submitted_Name ))

# Rerun Leipzig list with fuzzy matching
reclean_l       <- lapply( no_match_v$Submitted_Name, lcvp_fuzzy_search )
reclean_df      <- reclean_l %>% bind_rows
# visually select species identified with lcvp_fuzzy_search
# Fuzzy matches for these 20 taxa are either not reliable (only identified to genus level) or not available

# Final taxonomy files 
# Clean taxa should have LCVP search results
taxa_out        <- lapply( clean_df_final$Submitted_Name, get_clean_names ) %>% 
  bind_rows %>% 
  rename( Submitted_Name      = Search,
          First_matched_Name  = Input.Taxon,
          LCVP_Accepted_Taxon = Output.Taxon ) %>% 
  # check for accepted names
  mutate( mismatch_test = str_detect( First_matched_Name, 
                                      Submitted_Name ), site = 'palanan' )

# Do "taxa unresolved" by hand (taxa with no matches found), and add back in the submitted genus, family and IDlevel to enable future identification
taxa_unresvd    <- bind_rows(  mismatch_unresvd, no_match_v ) %>%
  inner_join( taxa_df ) %>%
  mutate( site = 'palanan' ) 


# store resolved AND unresolved taxa
write.csv( taxa_out, 'results/palanan_taxa.csv',
           row.names = F )
write.csv( taxa_unresvd, 'results/palanan_taxa_unresvd.csv',
           row.names = F )

