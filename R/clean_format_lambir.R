# lambir data preparation script
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
lambir_means <- dplyr::filter( kam_means, site == "lambir")
lambir_medians <- dplyr::filter( kam_medians, site == "lambir")

# Prepare site table -----------------------------

# Do this by hand
site_out <- data.frame( Site_name = 'Lambir',
                        Latitude  = 4.186500000000,
                        Longitude = 114.017000000000 )

# Check that the location is sensible sense
leaflet( data = site_out) %>% 
  addTiles() %>% 
  addCircleMarkers(~Longitude, ~Latitude)

# Store site information
write.csv( site_out, 'results/lambir_site.csv',
           row.names = F )

# Prepare taxonomic table ------------------------

# Produce the binomial used for checking
taxa_df         <- dplyr::select( lambir_means, latin, sp, genus, family, IDlevel ) %>%
  rename( Submitted_Name = latin, Sp_Code = sp, Submitted_Genus = genus, Submitted_Family = family )

# Separate NAs - no NAs
taxa_na_df      <- subset( taxa_df,  is.na( Submitted_Name ) )
taxa_na_rm_df   <- subset( taxa_df, !is.na( Submitted_Name ) )

# 'Shorea macroptera subsp. macroptera' prevents successful lcvp_search; move this taxon to unreserved
taxa_error      <- subset( taxa_na_rm_df, Submitted_Name == 'Shorea macroptera subsp. macroptera' )
taxa_na_rm_df   <- subset( taxa_na_rm_df, Submitted_Name != 'Shorea macroptera subsp. macroptera' )


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
# 967 plausible typos which remain in clean data frame

# Cinnamomum tahyanum is likely C. tahijanum, not C. bahianum - remove from clean data frame and replace
C_tahyanum_resvd <- check_mismatches[[54]] %>%
                      mutate( Submitted_Name = 'Cinnamomum tahyanum' ) %>%
                      rename( First_matched_Name  = Input.Taxon, 
                              LCVP_Accepted_Taxon = Output.Taxon ) %>%
                      subset( First_matched_Name == 'Cinnamomum tahijanum Kosterm.' )
clean_df <- full_join(clean_df, C_tahyanum_resvd) %>%
              subset( First_matched_Name != 'Cinnamomum bahianum Lukman.')

# remove ... unresolved species containing "aff.", "cf.", "Sp 39", "nom. ined.", "jvl", "nom. incert.", "sp. nov.", "af." and "hairy teeth" from clean data frame
# Archidendron casai and Castanopsis nivea are not plausible typos - move to unresolved
# Unclear whether Kayea longipedicellata is Genus Calea or Genus Mabea - move to unresolved
mismatch_unresvd <- data.frame( "Submitted_Name" = grep(  'aff.|cf.|Sp 39|ined.|jvl|incert.|sp. nov.|af.|hairy teeth|Archidendron casai|Castanopsis nivea|Kayea longipedicellata', mismatch_df$Submitted_Name, value = T )) %>%
  inner_join( taxa_df ) %>%
  distinct()
matched_df <- clean_df %>% subset( mismatch_test )
clean_df_final <- data.frame( "Submitted_Name" = grep( 'aff.|cf.|Sp 39|ined.|jvl|incert.|sp. nov.|af.|hairy teeth|Archidendron casai|Castanopsis nivea|Kayea longipedicellata', mismatch_df$Submitted_Name, value = T, invert = T )) %>% full_join( matched_df )

# Check species without matches; to avoid duplicates, use Sp_Code as a place holder
add_sp_code <- select(taxa_na_rm_df, Sp_Code, Submitted_Name)
clean_df_coded <- inner_join( clean_df, add_sp_code )
no_match_v <- data.frame( "Sp_Code" = setdiff( taxa_na_rm_df$Sp_Code, 
                                                       clean_df_coded$Sp_Code ))
no_match_v <- right_join( taxa_na_rm_df, no_match_v )

# Rerun Leipzig list with fuzzy matching
reclean_l       <- lapply( no_match_v$Submitted_Name, lcvp_fuzzy_search )
reclean_df      <- reclean_l %>% bind_rows
# from visual inspection, no unmatched species can be resolved via fuzzy search
# these 151 unmatched species remain unresolved

# Final taxonomy files 
# All clean taxa should have LCVP search results
taxa_out        <- lapply( clean_df_final$Submitted_Name, get_clean_names ) %>% 
  bind_rows %>% 
  rename( Submitted_Name      = Search,
          First_matched_Name  = Input.Taxon,
          LCVP_Accepted_Taxon = Output.Taxon ) %>% 
  # check for accepted names
  mutate( mismatch_test = str_detect( First_matched_Name, 
                                      Submitted_Name ), site = 'lambir' )

# Do "taxa unresolved" by hand (taxa with no matches found and unresolved mismatches); add in 'Shorea macroptera subsp. macroptera'

taxa_unresvd    <- bind_rows( taxa_error, mismatch_unresvd, no_match_v ) %>%
  inner_join( taxa_df ) %>%
  mutate( site = 'lambir' ) 

# store resolved AND unresolved taxa
write.csv( taxa_out, 'results/lambir_taxa.csv',
           row.names = F )
write.csv( taxa_unresvd, 'results/lambir_taxa_unresvd.csv',
           row.names = F )

# Prepare demographic table --------------------------------------

# distiguish taxa with sample size 0 for each growth layer (1-4) and survival layer (1-4)
lambir_means <- lambir_means %>% mutate( "growth_layer_1_imputed"   = grepl( "^0", growth_layer_1_obs ),
                                       "growth_layer_2_imputed"   = grepl( "^0", growth_layer_2_obs ), 
                                       "growth_layer_3_imputed"   = grepl( "^0", growth_layer_3_obs ), 
                                       "growth_layer_4_imputed"   = grepl( "^0", growth_layer_4_obs ), 
                                       "survival_layer_1_imputed" = grepl( "^0", survival_layer_1_obs ), 
                                       "survival_layer_2_imputed" = grepl( "^0", survival_layer_2_obs ), 
                                       "survival_layer_3_imputed" = grepl( "^0", survival_layer_3_obs ), 
                                       "survival_layer_4_imputed" = grepl( "^0", survival_layer_4_obs ) )

lambir_medians <- lambir_medians %>% mutate( "growth_layer_1_imputed"   = grepl( "^0", growth_layer_1_obs ),
                                           "growth_layer_2_imputed"   = grepl( "^0", growth_layer_2_obs ), 
                                           "growth_layer_3_imputed"   = grepl( "^0", growth_layer_3_obs ), 
                                           "growth_layer_4_imputed"   = grepl( "^0", growth_layer_4_obs ), 
                                           "survival_layer_1_imputed" = grepl( "^0", survival_layer_1_obs ), 
                                           "survival_layer_2_imputed" = grepl( "^0", survival_layer_2_obs ), 
                                           "survival_layer_3_imputed" = grepl( "^0", survival_layer_3_obs ), 
                                           "survival_layer_4_imputed" = grepl( "^0", survival_layer_4_obs ) )