# nanjenshan data preparation script
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
nanjenshan_means <- dplyr::filter( kam_means, site == "nanjenshan")
nanjenshan_medians <- dplyr::filter( kam_medians, site == "nanjenshan")

# Prepare site table -----------------------------

# Do this by hand
site_out <- data.frame( Site_name = 'Nanjenshan',
                        Latitude  = 22.056000000000,
                        Longitude = 120.860000000000 )

# Check that the location is sensible sense
leaflet( data = site_out) %>% 
  addTiles() %>% 
  addCircleMarkers(~Longitude, ~Latitude)

# Store site information
write.csv( site_out, 'results/nanjenshan_site.csv',
           row.names = F )

# Prepare taxonomic table ------------------------

# Produce the binomial used for checking
taxa_df         <- dplyr::select( nanjenshan_means, latin, sp, genus, family, IDlevel ) %>%
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
# Separate taxa named "Na NA"
taxa_Na_NA_df      <- subset( taxa_df,  Submitted_Name == 'Na NA' )
taxa_na_rm_df   <- subset( taxa_df, !is.na( Submitted_Name ) & Submitted_Name != 'Na NA' )


# Create function: get "cleaned" names
get_clean_names <- function( nam, fuzzy = 0.1 ){
  print( nam )
  lcvp_search( nam, max.distance = fuzzy )
}

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

# remove 13 unresolved species containing "_var_" and " NA" from clean data frame
mismatch_unresvd <- data.frame("Submitted_Name" = grep( '_var_| NA', mismatch_df$Submitted_Name, value = T )) %>%
 inner_join( taxa_df ) %>% 
  distinct()
matched_df <- clean_df %>% subset( mismatch_test ) 
clean_df_final <- data.frame("Submitted_Name" = grep( '_var_| NA', mismatch_df$Submitted_Name, value = T, invert = T )) %>% full_join( matched_df )

# Check species without matches
no_match_v <- data.frame( "Submitted_Name" = setdiff( taxa_na_rm_df$Submitted_Name, 
                                                      clean_df$Submitted_Name ))
# No unmatched taxa

# Final taxonomy files 
# Clean taxa should have LCVP search results
taxa_out        <- lapply( clean_df_final$Submitted_Name, get_clean_names ) %>% 
  bind_rows %>% 
  rename( Submitted_Name      = Search,
          First_matched_Name  = Input.Taxon,
          LCVP_Accepted_Taxon = Output.Taxon ) %>% 
  # check for accepted names
  mutate( mismatch_test = str_detect( First_matched_Name, 
                                      Submitted_Name ), site = 'nanjenshan' )

# Do "taxa unresolved" by hand (taxa with no matches found), and add back in the submitted genus, family and IDlevel to enable future identification
taxa_unresvd    <- bind_rows( taxa_Na_NA_df, mismatch_unresvd ) %>%
  mutate( site = 'nanjenshan' ) 
# Note there are two taxa "Na NA" both with the Sp_Code "celtfo"


# store resolved AND unresolved taxa
write.csv( taxa_out, 'results/nanjenshan_taxa.csv',
           row.names = F )
write.csv( taxa_unresvd, 'results/nanjenshan_taxa_unresvd.csv',
           row.names = F )

# Prepare demographic table --------------------------------------

taxa_out <- read.csv( 'results/nanjenshan_taxa.csv' , header = TRUE)
taxa_unresvd <- read.csv( 'results/nanjenshan_taxa_unresvd.csv' , header = TRUE)

# Select variables for demographic means tables - one for clean taxa and another for unresolved taxa
# Get demographic variables from nanjenshan_means
demog_means_df         <- select( nanjenshan_means, -c( genus, family, IDlevel ) ) %>%
  rename( Submitted_Name = latin, Sp_Code = sp )

#Capitalise taxa in demog_means_df
upper_case_genus <- function( x ){
  
  x %>% 
    separate( Submitted_Name, c('gen','spp'), 
              sep = ' ', extra = 'merge' ) %>% 
    mutate( gen   = str_to_title(gen) ) %>% 
    mutate( Submitted_Name = paste0(gen,' ',spp) ) %>% 
    dplyr::select( -gen, -spp)
  
}
demog_means_df <- upper_case_genus( demog_means_df )


# Join clean taxa to rest of schema via accepted names
demog_means_df_clean   <- data.frame( "Submitted_Name" = taxa_out$Submitted_Name,
                                      "LCVP_Accepted_Taxon" = taxa_out$LCVP_Accepted_Taxon ) %>%
  inner_join( demog_means_df ) %>%
  select( -c( Submitted_Name, Sp_Code ) )

# Join unresolved taxa to rest of schema via species codes and remove two extra taxa "Na NA" both with the Sp_Code "celtfo"
demog_means_df_unresvd <- data.frame( "Submitted_Name" = taxa_unresvd$Submitted_Name,
                                      "Sp_Code" = taxa_unresvd$Sp_Code ) %>%
  inner_join( demog_means_df ) %>%
  slice( -c( 9, 10 ) )

# Distiguish taxa with sample size 0 for each growth layer (1-4) and survival layer (1-4)
demog_means_df_clean <- demog_means_df_clean %>% 
  mutate( "growth_layer_1_imputed"   = grepl( "^0", growth_layer_1_obs ),
          "growth_layer_2_imputed"   = grepl( "^0", growth_layer_2_obs ), 
          "growth_layer_3_imputed"   = grepl( "^0", growth_layer_3_obs ), 
          "growth_layer_4_imputed"   = grepl( "^0", growth_layer_4_obs ), 
          "survival_layer_1_imputed" = grepl( "^0", survival_layer_1_obs ), 
          "survival_layer_2_imputed" = grepl( "^0", survival_layer_2_obs ), 
          "survival_layer_3_imputed" = grepl( "^0", survival_layer_3_obs ), 
          "survival_layer_4_imputed" = grepl( "^0", survival_layer_4_obs ) )

demog_means_df_unresvd <- demog_means_df_unresvd %>% 
  mutate( "growth_layer_1_imputed"   = grepl( "^0", growth_layer_1_obs ),
          "growth_layer_2_imputed"   = grepl( "^0", growth_layer_2_obs ), 
          "growth_layer_3_imputed"   = grepl( "^0", growth_layer_3_obs ), 
          "growth_layer_4_imputed"   = grepl( "^0", growth_layer_4_obs ), 
          "survival_layer_1_imputed" = grepl( "^0", survival_layer_1_obs ), 
          "survival_layer_2_imputed" = grepl( "^0", survival_layer_2_obs ), 
          "survival_layer_3_imputed" = grepl( "^0", survival_layer_3_obs ), 
          "survival_layer_4_imputed" = grepl( "^0", survival_layer_4_obs ) )

# store demographic means table for resolved AND unresolved taxa
write.csv( demog_means_df_clean, 'results/nanjenshan_demog_means.csv',
           row.names = F )
write.csv( demog_means_df_unresvd, 'results/nanjenshan_demog_means_unresvd.csv',
           row.names = F )

# Select variables for median tables - one for clean taxa and another for unresolved taxa
# Get demographic variables from nanjenshan_medians
demog_medians_df         <- select( nanjenshan_medians, -c( genus, family, IDlevel ) ) %>%
  rename( Submitted_Name = latin, Sp_Code = sp )

#Capitalise taxa in demog_medians_df
demog_medians_df <- upper_case_genus( demog_medians_df )

# Join clean taxa to rest of schema via accepted names
demog_medians_df_clean   <- data.frame( "Submitted_Name" = taxa_out$Submitted_Name,
                                        "LCVP_Accepted_Taxon" = taxa_out$LCVP_Accepted_Taxon ) %>%
  inner_join( demog_medians_df ) %>%
  select( -c( Submitted_Name, Sp_Code ) )

# Join unresolved taxa to rest of schema via species codes and remove two extra taxa "Na NA" both with the Sp_Code "celtfo"
demog_medians_df_unresvd <- data.frame( "Submitted_Name" = taxa_unresvd$Submitted_Name,
                                        "Sp_Code" = taxa_unresvd$Sp_Code ) %>%
  inner_join( demog_medians_df ) %>%
  slice( -c( 9, 10 ) )

# Distiguish taxa with sample size 0 for each growth layer (1-4) and survival layer (1-4)
demog_medians_df_clean <- demog_medians_df_clean %>% 
  mutate( "growth_layer_1_imputed"   = grepl( "^0", growth_layer_1_obs ),
          "growth_layer_2_imputed"   = grepl( "^0", growth_layer_2_obs ), 
          "growth_layer_3_imputed"   = grepl( "^0", growth_layer_3_obs ), 
          "growth_layer_4_imputed"   = grepl( "^0", growth_layer_4_obs ), 
          "survival_layer_1_imputed" = grepl( "^0", survival_layer_1_obs ), 
          "survival_layer_2_imputed" = grepl( "^0", survival_layer_2_obs ), 
          "survival_layer_3_imputed" = grepl( "^0", survival_layer_3_obs ), 
          "survival_layer_4_imputed" = grepl( "^0", survival_layer_4_obs ) )

demog_medians_df_unresvd <- demog_medians_df_unresvd %>% 
  mutate( "growth_layer_1_imputed"   = grepl( "^0", growth_layer_1_obs ),
          "growth_layer_2_imputed"   = grepl( "^0", growth_layer_2_obs ), 
          "growth_layer_3_imputed"   = grepl( "^0", growth_layer_3_obs ), 
          "growth_layer_4_imputed"   = grepl( "^0", growth_layer_4_obs ), 
          "survival_layer_1_imputed" = grepl( "^0", survival_layer_1_obs ), 
          "survival_layer_2_imputed" = grepl( "^0", survival_layer_2_obs ), 
          "survival_layer_3_imputed" = grepl( "^0", survival_layer_3_obs ), 
          "survival_layer_4_imputed" = grepl( "^0", survival_layer_4_obs ) )

# store demographic medians table for resolved AND unresolved taxa
write.csv( demog_medians_df_clean, 'results/nanjenshan_demog_medians.csv',
           row.names = F )
write.csv( demog_medians_df_unresvd, 'results/nanjenshan_demog_medians_unresvd.csv',
           row.names = F )


# Analyse relationship between sample size and CI width
# Create variable for CI width
nanjenshan_means <- nanjenshan_means %>% mutate( "growth_layer_1_CI90_width"   = growth_layer1_CI.95 - growth_layer1_CI.05,
                                   "growth_layer_2_CI90_width"   = growth_layer2_CI.95 - growth_layer2_CI.05, 
                                   "growth_layer_3_CI90_width"   = growth_layer3_CI.95 - growth_layer3_CI.05, 
                                   "growth_layer_4_CI90_width"   = growth_layer4_CI.95 - growth_layer4_CI.05, 
                                   "survival_layer_1_CI90_width" = survival_layer1_CI.95 - survival_layer1_CI.05, 
                                   "survival_layer_2_CI90_width" = survival_layer2_CI.95 - survival_layer2_CI.05, 
                                   "survival_layer_3_CI90_width" = survival_layer3_CI.95 - survival_layer3_CI.05, 
                                   "survival_layer_4_CI90_width" = survival_layer4_CI.95 - survival_layer4_CI.05 )



# Load plots in a 2x2 grid 
growth_layer_1_graph <- ggplot( data = subset( nanjenshan_means, growth_layer_1_imputed != TRUE), aes( x = growth_layer_1_obs, y = growth_layer_1_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Growth Layer 1")
growth_layer_2_graph <- ggplot( data = subset( nanjenshan_means, growth_layer_2_imputed != TRUE), aes( x = growth_layer_2_obs, y = growth_layer_2_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Growth Layer 2")
growth_layer_3_graph <- ggplot( data = subset( nanjenshan_means, growth_layer_3_imputed != TRUE), aes( x = growth_layer_3_obs, y = growth_layer_3_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Growth Layer 3")
growth_layer_4_graph <- ggplot( data = subset( nanjenshan_means, growth_layer_4_imputed != TRUE), aes( x = growth_layer_4_obs, y = growth_layer_4_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Growth Layer 4")
grid.arrange( growth_layer_1_graph,
              growth_layer_2_graph,
              growth_layer_3_graph,
              growth_layer_4_graph )

survival_layer_1_graph <- ggplot( data = subset( nanjenshan_means, survival_layer_1_imputed != TRUE), aes( x = survival_layer_1_obs, y = survival_layer_1_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Survival Layer 1")
survival_layer_2_graph <- ggplot( data = subset( nanjenshan_means, survival_layer_2_imputed != TRUE), aes( x = survival_layer_2_obs, y = survival_layer_2_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Survival Layer 2")
survival_layer_3_graph <- ggplot( data = subset( nanjenshan_means, survival_layer_3_imputed != TRUE), aes( x = survival_layer_3_obs, y = survival_layer_3_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Survival Layer 3")
survival_layer_4_graph <- ggplot( data = subset( nanjenshan_means, survival_layer_4_imputed != TRUE), aes( x = survival_layer_4_obs, y = survival_layer_4_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Survival Layer 4")
grid.arrange( survival_layer_1_graph,
              survival_layer_2_graph,
              survival_layer_3_graph,
              survival_layer_4_graph )
