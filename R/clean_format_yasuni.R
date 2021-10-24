# yasuni data preparation script
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
yasuni_means <- dplyr::filter( kam_means, site == "yasuni")
yasuni_medians <- dplyr::filter( kam_medians, site == "yasuni")

# Prepare site table -----------------------------

# Do this by hand
site_out <- data.frame( Site_name = 'Yasuni',
                        Latitude  = -0.685900000000,
                        Longitude = -76.397000000000 )

# Check that the location is sensible sense
leaflet( data = site_out) %>% 
  addTiles() %>% 
  addCircleMarkers(~Longitude, ~Latitude)

# Store site information
write.csv( site_out, 'results/yasuni_site.csv',
           row.names = F )

# Prepare taxonomic table ------------------------

# Produce the binomial used for checking
taxa_df         <- dplyr::select( yasuni_means, latin, sp, genus, family, IDlevel ) %>%
  rename( Submitted_Name = latin, Sp_Code = sp, Submitted_Genus = genus, Submitted_Family = family )

# Separate NAs
taxa_na_df      <- subset( taxa_df,  is.na( Submitted_Name ) )
taxa_na_rm_df   <- subset( taxa_df, !is.na( Submitted_Name ) )

# Create function: get "cleaned" names
get_clean_names <- function( nam, fuzzy = 0.1 ){
  print( nam )
  lcvp_search( nam, max.distance = fuzzy )
}

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
# 25 plausible typos which remain in clean data frame

# several potential typos are possible for 6 taxa; leave them in unresolved table
# no mismatch is plausible for 60 taxa; leave them in unresolved table
unresvd_list <- c('Bunchosia myrt',
                  'Calyptranthes punteada',
                  'Eugenia margot',
                  'Miconia corine',
                  'Piper obglab',
                  'Piper obtomen',
                  'Acalypha pub',
                  'Annona mosaic',
                  'Brownea lore',
                  'Calyptranthes grancauli',
                  'Caryodaphnopsis chica',
                  'Chrysochlamys hugo',
                  'Chrysophyllum minor',
                  'Chrysophyllum tremi',
                  'Coccoloba jill',
                  'Coccoloba ninfi',
                  'Coccoloba papel',
                  'Colubrina arbol',
                  'Dacryodes gorda',
                  'Endlicheria dori',
                  'Eschweilera giga',
                  'Heisteria grande',
                  'Inga 3crasa',
                  'Inga 3oscura',
                  'Inga 6cuadra',
                  'Lacistema med',
                  'Licania filo',
                  'Licania opaca',
                  'Margaritaria roja',
                  'Marila comun',
                  'Matayba ocho',
                  'Maytenus ala',
                  'Meliosma doly',
                  'Miconia karina',
                  'Miconia tipica',
                  'Neea aniboid',
                  'Neea claudia',
                  'Neea garci',
                  'Neea gigante',
                  'Neea micro',
                  'Ocotea luis',
                  'Pera duguet',
                  'Picramnia mini',
                  'Piper bulada',
                  'Piper obchic',
                  'Piper obvil',
                  'Piper renato',
                  'Piper sesivil',
                  'Plinia unop',
                  'Pouteria redondita',
                  'Psychotria dracula',
                  'Quiina mediana',
                  'Randia manolo',
                  'Rudgea fina',
                  'Sloanea guia',
                  'Sloanea oak',
                  'Sloanea oppd',
                  'Solanum granmini',
                  'Talisia 2-retic',
                  'Terminalia ob',
                  'Tovomita tyana',
                  'Guatteria hiena',
                  'Zanthoxylum nervi',
                  'Zanthoxylum perp',
                  'Zanthoxylum suave',
                  'Alibertia lance'
)
mismatch_unresvd <- data.frame( "Submitted_Name" = unresvd_list ) %>% 
  inner_join( taxa_df )
clean_df_final <- data.frame( "Submitted_Name" = setdiff( clean_df$Submitted_Name, 
                                                          mismatch_unresvd$Submitted_Name ) )

# Check species without matches
no_match_v <- data.frame( "Submitted_Name" = setdiff( taxa_na_rm_df$Submitted_Name, 
                                                        clean_df$Submitted_Name ) ) %>% 
  inner_join( taxa_df )

# Rerun Leipzig list with fuzzy matching
reclean_l       <- lapply( no_match_v$Submitted_Name, lcvp_fuzzy_search )
reclean_df      <- reclean_l %>% bind_rows
# Fuzzy matches for these 218 taxa are not available

# Final taxonomy files 
taxa_out        <- lapply( clean_df_final$Submitted_Name, get_clean_names ) %>% 
  bind_rows %>% 
  rename( Submitted_Name      = Search,
          First_matched_Name  = Input.Taxon,
          LCVP_Accepted_Taxon = Output.Taxon ) %>% 
  # check for accepted names
  mutate( mismatch_test = str_detect( First_matched_Name, 
                                      Submitted_Name ), site = 'yasuni' )

# Do "taxa unresolved" by hand (taxa with no matches found), and add back in the submitted genus, family and IDlevel to enable future identification
taxa_unresvd    <- bind_rows( taxa_na_df, mismatch_unresvd, no_match_v ) %>%
  mutate( site = 'yasuni' ) %>%
  distinct()

# store resolved AND unresolved taxa
write.csv( taxa_out, 'results/yasuni_taxa.csv',
           row.names = F )
write.csv( taxa_unresvd, 'results/yasuni_taxa_unresvd.csv',
           row.names = F )

# Prepare demographic table --------------------------------------

# distiguish taxa with sample size 0 for each growth layer (1-4) and survival layer (1-4)
yasuni_means <- yasuni_means %>% mutate( "growth_layer_1_imputed"   = grepl( "^0", growth_layer_1_obs ),
                                               "growth_layer_2_imputed"   = grepl( "^0", growth_layer_2_obs ), 
                                               "growth_layer_3_imputed"   = grepl( "^0", growth_layer_3_obs ), 
                                               "growth_layer_4_imputed"   = grepl( "^0", growth_layer_4_obs ), 
                                               "survival_layer_1_imputed" = grepl( "^0", survival_layer_1_obs ), 
                                               "survival_layer_2_imputed" = grepl( "^0", survival_layer_2_obs ), 
                                               "survival_layer_3_imputed" = grepl( "^0", survival_layer_3_obs ), 
                                               "survival_layer_4_imputed" = grepl( "^0", survival_layer_4_obs ) )

yasuni_medians <- yasuni_medians %>% mutate( "growth_layer_1_imputed"   = grepl( "^0", growth_layer_1_obs ),
                                                   "growth_layer_2_imputed"   = grepl( "^0", growth_layer_2_obs ), 
                                                   "growth_layer_3_imputed"   = grepl( "^0", growth_layer_3_obs ), 
                                                   "growth_layer_4_imputed"   = grepl( "^0", growth_layer_4_obs ), 
                                                   "survival_layer_1_imputed" = grepl( "^0", survival_layer_1_obs ), 
                                                   "survival_layer_2_imputed" = grepl( "^0", survival_layer_2_obs ), 
                                                   "survival_layer_3_imputed" = grepl( "^0", survival_layer_3_obs ), 
                                                   "survival_layer_4_imputed" = grepl( "^0", survival_layer_4_obs ) )