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
library(ggplot2)
library(gridExtra)

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
taxa_df         <- dplyr::select( ituri_edoro_means, latin, sp, genus, family, IDlevel ) %>%
  rename( Submitted_Name = latin, Sp_Code = sp, Submitted_Genus = genus, Submitted_Family = family )

# Separate NAs
taxa_na_df      <- subset( taxa_df,  is.na( Submitted_Name ) )
taxa_na_rm_df   <- subset( taxa_df, !is.na( Submitted_Name ) )

# Create function: get "cleaned" names
get_clean_names <- function( nam, fuzzy = 0.1 ) lcvp_search( nam, max.distance = fuzzy )

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
mismatch_df      <- clean_df %>% subset( !mismatch_test )
check_mismatches <- lapply( mismatch_df$Submitted_Name, lcvp_fuzzy_search )
# 42 plausible typos which remain in clean data frame

# remove 4 unresolved species containing "sp1" "sp2" "sp3" and "xx" from clean data frame
mismatch_unresvd <- data.frame( "Submitted_Name" = grep( 'sp1|sp2|sp3|xx', mismatch_df$Submitted_Name, value = T ) )
clean_df         <- data.frame( "Submitted_Name" = grep( 'sp1|sp2|sp3|xx', clean_df$Submitted_Name, value = T, invert = T ) )

# Check species without matches
no_match_v       <- data.frame( "Submitted_Name" = setdiff( taxa_na_rm_df$Submitted_Name, 
                                                            clean_df$Submitted_Name ))

# Rerun Leipzig list with fuzzy matching
reclean_l       <- lapply( no_match_v, lcvp_fuzzy_search )
reclean_df      <- reclean_l %>% bind_rows
# visually select species identified with lcvp_fuzzy_search
# No fuzzy match for 13 species

# Clean taxa should have LCVP search results
clean_df_final   <- lapply( clean_df$Submitted_Name, get_clean_names ) %>% 
                            bind_rows %>% 
                            rename( Submitted_Name      = Search,
                            First_matched_Name  = Input.Taxon,
                            LCVP_Accepted_Taxon = Output.Taxon ) %>% 
                    # check for accepted names
                    mutate( mismatch_test = str_detect( First_matched_Name, Submitted_Name ), 
                            site = 'ituri_edoro' )

# Check clean dataframe for duplications in accepted taxa
duplicates      <- clean_df_final$LCVP_Accepted_Taxon[ duplicated( clean_df_final$LCVP_Accepted_Taxon ) ]
duplicates_df   <- clean_df_final[ clean_df_final$LCVP_Accepted_Taxon %in% duplicates, ]
# "Cuviera nigrescens (Elliott ex Oliv.) Wernham" is counted twice and must be removed, label issue as "double count"
double_counts   <- c("Cuviera nigrescens (Elliott ex Oliv.) Wernham")
double_count_df <- duplicates_df[ duplicates_df$LCVP_Accepted_Taxon %in% double_counts, ]
double_count_df <- mutate( double_count_df, issue = 'double count' )
# Several taxa duplicated as synonyms must be removed, label issue as "synonym"
synonyms        <- c("Allophylus africanus P.Beauv.", 
                     "Englerophytum iturense (Engl.) L.Gaut.", 
                     "Phyllocosmus africanus (Hook.f.) Klotzsch", 
                     "Vepris afzelii (Engl.) Mziray")
synonyms_df     <- duplicates_df[ duplicates_df$LCVP_Accepted_Taxon %in% synonyms, ]
synonyms_df     <- mutate( synonyms_df, issue = 'synonym' )

# Final resolved taxonomic file with duplicates removed
taxa_out        <- clean_df_final[ !clean_df_final$LCVP_Accepted_Taxon %in% duplicates, ]


# Do "taxa unresolved" by hand (taxa with no matches found), and add back in the submitted genus, family and IDlevel to enable future identification, labelling issue as "not in LCVP" for unresolved taxa
not_in_LCVP     <-  bind_rows( mismatch_unresvd, no_match_v ) %>%
                    inner_join( taxa_df ) %>%
                    bind_rows( taxa_na_df ) %>%
                    mutate( issue = 'not in LCVP' )
taxa_unresvd    <-  bind_rows( not_in_LCVP, double_count_df, synonyms_df ) %>%
                    mutate( site = 'ituri_edoro' )

# store resolved AND unresolved taxa
write.csv( taxa_out, 'results/ituri_edoro_taxa.csv',
           row.names = F )
write.csv( taxa_unresvd, 'results/ituri_edoro_taxa_unresvd.csv',
           row.names = F )

# Prepare demographic table --------------------------------------

taxa_out <- read.csv( 'results/ituri_edoro_taxa.csv' , header = TRUE)
taxa_unresvd <- read.csv( 'results/ituri_edoro_taxa_unresvd.csv' , header = TRUE)

# Select variables for demographic means tables - one for clean taxa and another for unresolved taxa
# Get demographic variables from ituri_edoro_means
demog_means_df         <- select( ituri_edoro_means, -c( genus, family, IDlevel ) ) %>%
  rename( Submitted_Name = latin, Sp_Code = sp )

# Join clean taxa to rest of schema via accepted names
demog_means_df_clean   <- data.frame( "Submitted_Name" = taxa_out$Submitted_Name,
                                      "LCVP_Accepted_Taxon" = taxa_out$LCVP_Accepted_Taxon ) %>%
  inner_join( demog_means_df ) %>%
  select( -c( Submitted_Name, Sp_Code ) )

# Join unresolved taxa to rest of schema via species codes
demog_means_df_unresvd <- data.frame( "Submitted_Name" = taxa_unresvd$Submitted_Name,
                                      "Sp_Code" = taxa_unresvd$Sp_Code ) %>%
  inner_join( demog_means_df, by = "Submitted_Name" ) %>%
  select( -c( Sp_Code.x ) ) %>%
  rename( Sp_Code = Sp_Code.y ) %>%
  distinct( .keep_all = TRUE )

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
write.csv( demog_means_df_clean, 'results/ituri_edoro_demog_means.csv',
           row.names = F )
write.csv( demog_means_df_unresvd, 'results/ituri_edoro_demog_means_unresvd.csv',
           row.names = F )

# Select variables for median tables - one for clean taxa and another for unresolved taxa
# Get demographic variables from ituri_edoro_medians
demog_medians_df         <- select( ituri_edoro_medians, -c( genus, family, IDlevel ) ) %>%
  rename( Submitted_Name = latin, Sp_Code = sp )

# Join clean taxa to rest of schema via accepted names
demog_medians_df_clean   <- data.frame( "Submitted_Name" = taxa_out$Submitted_Name,
                                        "LCVP_Accepted_Taxon" = taxa_out$LCVP_Accepted_Taxon ) %>%
  inner_join( demog_medians_df ) %>%
  select( -c( Submitted_Name, Sp_Code ) )

# Join unresolved taxa to rest of schema via species codes
demog_medians_df_unresvd <- data.frame( "Submitted_Name" = taxa_unresvd$Submitted_Name,
                                      "Sp_Code" = taxa_unresvd$Sp_Code ) %>%
  inner_join( demog_medians_df, by = "Submitted_Name" ) %>%
  select( -c( Sp_Code.x ) ) %>%
  rename( Sp_Code = Sp_Code.y ) %>%
  distinct( .keep_all = TRUE )

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
write.csv( demog_medians_df_clean, 'results/ituri_edoro_demog_medians.csv',
           row.names = F )
write.csv( demog_medians_df_unresvd, 'results/ituri_edoro_demog_medians_unresvd.csv',
           row.names = F )


# Analyse relationship between sample size and CI width
# Create variable for CI width
ituri_edoro_means <- ituri_edoro_means %>% mutate( "growth_layer_1_CI90_width"   = growth_layer1_CI.95 - growth_layer1_CI.05,
                                                 "growth_layer_2_CI90_width"   = growth_layer2_CI.95 - growth_layer2_CI.05, 
                                                 "growth_layer_3_CI90_width"   = growth_layer3_CI.95 - growth_layer3_CI.05, 
                                                 "growth_layer_4_CI90_width"   = growth_layer4_CI.95 - growth_layer4_CI.05, 
                                                 "survival_layer_1_CI90_width" = survival_layer1_CI.95 - survival_layer1_CI.05, 
                                                 "survival_layer_2_CI90_width" = survival_layer2_CI.95 - survival_layer2_CI.05, 
                                                 "survival_layer_3_CI90_width" = survival_layer3_CI.95 - survival_layer3_CI.05, 
                                                 "survival_layer_4_CI90_width" = survival_layer4_CI.95 - survival_layer4_CI.05 )

# Load plots in a 2x2 grid 
growth_layer_1_graph <- ggplot( data = subset( ituri_edoro_means, growth_layer_1_imputed != TRUE), aes( x = growth_layer_1_obs, y = growth_layer_1_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Growth Layer 1")
growth_layer_2_graph <- ggplot( data = subset( ituri_edoro_means, growth_layer_2_imputed != TRUE), aes( x = growth_layer_2_obs, y = growth_layer_2_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Growth Layer 2")
growth_layer_3_graph <- ggplot( data = subset( ituri_edoro_means, growth_layer_3_imputed != TRUE), aes( x = growth_layer_3_obs, y = growth_layer_3_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Growth Layer 3")
growth_layer_4_graph <- ggplot( data = subset( ituri_edoro_means, growth_layer_4_imputed != TRUE), aes( x = growth_layer_4_obs, y = growth_layer_4_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Growth Layer 4")
grid.arrange( growth_layer_1_graph,
              growth_layer_2_graph,
              growth_layer_3_graph,
              growth_layer_4_graph )

survival_layer_1_graph <- ggplot( data = subset( ituri_edoro_means, survival_layer_1_imputed != TRUE), aes( x = survival_layer_1_obs, y = survival_layer_1_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Survival Layer 1")
survival_layer_2_graph <- ggplot( data = subset( ituri_edoro_means, survival_layer_2_imputed != TRUE), aes( x = survival_layer_2_obs, y = survival_layer_2_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Survival Layer 2")
survival_layer_3_graph <- ggplot( data = subset( ituri_edoro_means, survival_layer_3_imputed != TRUE), aes( x = survival_layer_3_obs, y = survival_layer_3_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Survival Layer 3")
survival_layer_4_graph <- ggplot( data = subset( ituri_edoro_means, survival_layer_4_imputed != TRUE), aes( x = survival_layer_4_obs, y = survival_layer_4_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Survival Layer 4")
grid.arrange( survival_layer_1_graph,
              survival_layer_2_graph,
              survival_layer_3_graph,
              survival_layer_4_graph )
