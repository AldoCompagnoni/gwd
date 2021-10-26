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
library(ggplot2)
library(gridExtra)


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

# Prepare taxonomic table ------------------------

# Produce the binomial used for checking
taxa_df         <- dplyr::select( luquillo_means, latin, sp, genus, family, IDlevel ) %>%
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

# Separate NAs 
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
# 2 plausible typos which remain in clean data frame

# remove 5 unresolved species containing "X", "spp" "Â" and "NA NA" from clean data frame
mismatch_unresvd <- data.frame("Submitted_Name" = grep( 'X|spp|Â|NA NA', mismatch_df$Submitted_Name, value = T ))
matched_df <- clean_df %>% subset( mismatch_test )
clean_df_final <- data.frame("Submitted_Name" = grep( 'X|spp|Â|NA NA', mismatch_df$Submitted_Name, value = T, invert = T )) %>% full_join( matched_df )

# Check species without matches
no_match_v <- data.frame( "Submitted_Name" = setdiff( taxa_na_rm_df$Submitted_Name, 
                                                      clean_df$Submitted_Name ))
# from visual inspection, 1 unmatched taxon is unresolved

# Final taxonomy files 
# Clean taxa should have LCVP search results
taxa_out        <- lapply( clean_df_final$Submitted_Name, get_clean_names ) %>% 
  bind_rows %>% 
  rename( Submitted_Name      = Search,
          First_matched_Name  = Input.Taxon,
          LCVP_Accepted_Taxon = Output.Taxon ) %>% 
  # check for accepted names
  mutate( mismatch_test = str_detect( First_matched_Name, 
                                      Submitted_Name ), site = 'luquillo' )

# Do "taxa unresolved" by hand (taxa with no matches found), and add back in the submitted genus, family and IDlevel to enable future identification
taxa_unresvd    <- bind_rows( mismatch_unresvd, no_match_v ) %>%
  inner_join( taxa_df ) %>%
  mutate( site = 'luquillo' ) %>%
  distinct()


# store resolved AND unresolved taxa
write.csv( taxa_out, 'results/luquillo_taxa.csv',
           row.names = F )
write.csv( taxa_unresvd, 'results/luquillo_taxa_unresvd.csv',
           row.names = F )

# Prepare demographic table --------------------------------------

# distiguish taxa with sample size 0 for each growth layer (1-4) and survival layer (1-4)
luquillo_means <- luquillo_means %>% mutate( "growth_layer_1_imputed"   = grepl( "^0", growth_layer_1_obs ),
                                         "growth_layer_2_imputed"   = grepl( "^0", growth_layer_2_obs ), 
                                         "growth_layer_3_imputed"   = grepl( "^0", growth_layer_3_obs ), 
                                         "growth_layer_4_imputed"   = grepl( "^0", growth_layer_4_obs ), 
                                         "survival_layer_1_imputed" = grepl( "^0", survival_layer_1_obs ), 
                                         "survival_layer_2_imputed" = grepl( "^0", survival_layer_2_obs ), 
                                         "survival_layer_3_imputed" = grepl( "^0", survival_layer_3_obs ), 
                                         "survival_layer_4_imputed" = grepl( "^0", survival_layer_4_obs ) )

luquillo_medians <- luquillo_medians %>% mutate( "growth_layer_1_imputed"   = grepl( "^0", growth_layer_1_obs ),
                                             "growth_layer_2_imputed"   = grepl( "^0", growth_layer_2_obs ), 
                                             "growth_layer_3_imputed"   = grepl( "^0", growth_layer_3_obs ), 
                                             "growth_layer_4_imputed"   = grepl( "^0", growth_layer_4_obs ), 
                                             "survival_layer_1_imputed" = grepl( "^0", survival_layer_1_obs ), 
                                             "survival_layer_2_imputed" = grepl( "^0", survival_layer_2_obs ), 
                                             "survival_layer_3_imputed" = grepl( "^0", survival_layer_3_obs ), 
                                             "survival_layer_4_imputed" = grepl( "^0", survival_layer_4_obs ) )

# Analyse relationship between sample size and CI width
# Create variable for CI width
luquillo_means <- luquillo_means %>% mutate( "growth_layer_1_CI90_width"   = growth_layer1_CI.95 - growth_layer1_CI.05,
                                   "growth_layer_2_CI90_width"   = growth_layer2_CI.95 - growth_layer2_CI.05, 
                                   "growth_layer_3_CI90_width"   = growth_layer3_CI.95 - growth_layer3_CI.05, 
                                   "growth_layer_4_CI90_width"   = growth_layer4_CI.95 - growth_layer4_CI.05, 
                                   "survival_layer_1_CI90_width" = survival_layer1_CI.95 - survival_layer1_CI.05, 
                                   "survival_layer_2_CI90_width" = survival_layer2_CI.95 - survival_layer2_CI.05, 
                                   "survival_layer_3_CI90_width" = survival_layer3_CI.95 - survival_layer3_CI.05, 
                                   "survival_layer_4_CI90_width" = survival_layer4_CI.95 - survival_layer4_CI.05 )

# Load plots in a 2x2 grid 
growth_layer_1_graph <- ggplot( data = subset( luquillo_means, growth_layer_1_imputed != TRUE), aes( x = growth_layer_1_obs, y = growth_layer_1_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Growth Layer 1")
growth_layer_2_graph <- ggplot( data = subset( luquillo_means, growth_layer_2_imputed != TRUE), aes( x = growth_layer_2_obs, y = growth_layer_2_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Growth Layer 2")
growth_layer_3_graph <- ggplot( data = subset( luquillo_means, growth_layer_3_imputed != TRUE), aes( x = growth_layer_3_obs, y = growth_layer_3_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Growth Layer 3")
growth_layer_4_graph <- ggplot( data = subset( luquillo_means, growth_layer_4_imputed != TRUE), aes( x = growth_layer_4_obs, y = growth_layer_4_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Growth Layer 4")
grid.arrange( growth_layer_1_graph,
              growth_layer_2_graph,
              growth_layer_3_graph,
              growth_layer_4_graph )

survival_layer_1_graph <- ggplot( data = subset( luquillo_means, survival_layer_1_imputed != TRUE), aes( x = survival_layer_1_obs, y = survival_layer_1_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Survival Layer 1")
survival_layer_2_graph <- ggplot( data = subset( luquillo_means, survival_layer_2_imputed != TRUE), aes( x = survival_layer_2_obs, y = survival_layer_2_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Survival Layer 2")
survival_layer_3_graph <- ggplot( data = subset( luquillo_means, survival_layer_3_imputed != TRUE), aes( x = survival_layer_3_obs, y = survival_layer_3_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Survival Layer 3")
survival_layer_4_graph <- ggplot( data = subset( luquillo_means, survival_layer_4_imputed != TRUE), aes( x = survival_layer_4_obs, y = survival_layer_4_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Survival Layer 4")
grid.arrange( survival_layer_1_graph,
              survival_layer_2_graph,
              survival_layer_3_graph,
              survival_layer_4_graph )
