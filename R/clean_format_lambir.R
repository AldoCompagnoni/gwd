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
library(ggplot2)
library(gridExtra)

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
clean_df_matched <- data.frame( "Submitted_Name" = grep( 'aff.|cf.|Sp 39|ined.|jvl|incert.|sp. nov.|af.|hairy teeth|Archidendron casai|Castanopsis nivea|Kayea longipedicellata', mismatch_df$Submitted_Name, value = T, invert = T )) %>% full_join( matched_df )

# Check species without matches; to avoid duplicates, use Sp_Code as a place holder
add_sp_code    <- select( taxa_na_rm_df, Sp_Code, Submitted_Name  )
clean_df_coded <- inner_join( clean_df, add_sp_code )
no_match_v     <- data.frame( "Sp_Code" = setdiff( taxa_na_rm_df$Sp_Code, 
                                                       clean_df_coded$Sp_Code ))
no_match_v     <- right_join( taxa_na_rm_df, no_match_v )

# Rerun Leipzig list with fuzzy matching
reclean_l       <- lapply( no_match_v$Submitted_Name, lcvp_fuzzy_search )
reclean_df      <- reclean_l %>% bind_rows
# from visual inspection, no unmatched species can be resolved via fuzzy search
# these 151 unmatched species remain unresolved

# Clean taxa should have LCVP search results
clean_df_final   <- lapply( clean_df_matched$Submitted_Name, get_clean_names ) %>% 
  bind_rows %>% 
  rename( Submitted_Name      = Search,
          First_matched_Name  = Input.Taxon,
          LCVP_Accepted_Taxon = Output.Taxon ) %>% 
  # check for accepted names
  mutate( mismatch_test = str_detect( First_matched_Name, Submitted_Name ), 
          site = 'lambir' )

# Check clean dataframe for duplications in accepted taxa
duplicates      <- clean_df_final$LCVP_Accepted_Taxon[ duplicated( clean_df_final$LCVP_Accepted_Taxon ) ]
duplicates_df   <- clean_df_final[ clean_df_final$LCVP_Accepted_Taxon %in% duplicates, ]
# Several taxa counted twice and must be removed, label issue as "double count"
double_counts   <- c( "Dacryodes incurvata (Engl.) H.J.Lam",
                      "Litsea oppositifolia Gibbs",
                      "Santiria griffithii (Hook.f.) Engl.",
                      "Shorea macroptera Dyer",
                      "Shorea scaberrima Burck",
                      "Syzygium garciniifolium (King) Merr. & L.M.Perry",
                      "Syzygium oligomyrum Diels",
                      "Vatica oblongifolia Hook.f.")
double_count_df <- duplicates_df[ duplicates_df$LCVP_Accepted_Taxon %in% double_counts, ]
double_count_df <- mutate( double_count_df, issue = 'double count' )
# Duplicated axa which are not double counted are duplicated as a taxonomic synonym and must be removed, label issue as "synonym"
synonyms_df     <- duplicates_df[ !duplicates_df$LCVP_Accepted_Taxon %in% double_counts, ]
synonyms_df     <- mutate( synonyms_df, issue = 'synonym' )

# Final resolved taxonomic file with duplicates removed
taxa_out        <- clean_df_final[ !clean_df_final$LCVP_Accepted_Taxon %in% duplicates, ]

# Do "taxa unresolved" by hand (taxa with no matches found), and add back in the submitted genus, family and IDlevel to enable future identification, labelling issue as "not in LCVP" for unresolved taxa
not_in_LCVP     <-  bind_rows( mismatch_unresvd, no_match_v ) %>%
                    inner_join( taxa_df ) %>%
                    bind_rows( taxa_na_df, taxa_error ) %>%
                    mutate( issue = 'not in LCVP' )
taxa_unresvd    <- bind_rows( not_in_LCVP, double_count_df, synonyms_df ) %>%
                   mutate( site = 'korup' )

# store resolved AND unresolved taxa
write.csv( taxa_out, 'results/lambir_taxa.csv',
           row.names = F )
write.csv( taxa_unresvd, 'results/lambir_taxa_unresvd.csv',
           row.names = F )

# Prepare demographic table --------------------------------------

taxa_out <- read.csv( 'results/lambir_taxa.csv' , header = TRUE)
taxa_unresvd <- read.csv( 'results/lambir_taxa_unresvd.csv' , header = TRUE)

# Select variables for demographic means tables - one for clean taxa and another for unresolved taxa
# Get demographic variables from lambir_means
demog_means_df         <- select( lambir_means, -c( genus, family, IDlevel ) ) %>%
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
write.csv( demog_means_df_clean, 'results/lambir_demog_means.csv',
           row.names = F )
write.csv( demog_means_df_unresvd, 'results/lambir_demog_means_unresvd.csv',
           row.names = F )

# Select variables for median tables - one for clean taxa and another for unresolved taxa
# Get demographic variables from lambir_medians
demog_medians_df         <- select( lambir_medians, -c( genus, family, IDlevel ) ) %>%
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
write.csv( demog_medians_df_clean, 'results/lambir_demog_medians.csv',
           row.names = F )
write.csv( demog_medians_df_unresvd, 'results/lambir_demog_medians_unresvd.csv',
           row.names = F )

# Analyse relationship between sample size and CI width
# Create variable for CI width
lambir_means <- lambir_means %>% mutate( "growth_layer_1_CI90_width"   = growth_layer1_CI.95 - growth_layer1_CI.05,
                                                 "growth_layer_2_CI90_width"   = growth_layer2_CI.95 - growth_layer2_CI.05, 
                                                 "growth_layer_3_CI90_width"   = growth_layer3_CI.95 - growth_layer3_CI.05, 
                                                 "growth_layer_4_CI90_width"   = growth_layer4_CI.95 - growth_layer4_CI.05, 
                                                 "survival_layer_1_CI90_width" = survival_layer1_CI.95 - survival_layer1_CI.05, 
                                                 "survival_layer_2_CI90_width" = survival_layer2_CI.95 - survival_layer2_CI.05, 
                                                 "survival_layer_3_CI90_width" = survival_layer3_CI.95 - survival_layer3_CI.05, 
                                                 "survival_layer_4_CI90_width" = survival_layer4_CI.95 - survival_layer4_CI.05 )

# Load plots in a 2x2 grid 
growth_layer_1_graph <- ggplot( data = subset( lambir_means, growth_layer_1_imputed != TRUE), aes( x = growth_layer_1_obs, y = growth_layer_1_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Growth Layer 1")
growth_layer_2_graph <- ggplot( data = subset( lambir_means, growth_layer_2_imputed != TRUE), aes( x = growth_layer_2_obs, y = growth_layer_2_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Growth Layer 2")
growth_layer_3_graph <- ggplot( data = subset( lambir_means, growth_layer_3_imputed != TRUE), aes( x = growth_layer_3_obs, y = growth_layer_3_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Growth Layer 3")
growth_layer_4_graph <- ggplot( data = subset( lambir_means, growth_layer_4_imputed != TRUE), aes( x = growth_layer_4_obs, y = growth_layer_4_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Growth Layer 4")
grid.arrange( growth_layer_1_graph,
              growth_layer_2_graph,
              growth_layer_3_graph,
              growth_layer_4_graph )

survival_layer_1_graph <- ggplot( data = subset( lambir_means, survival_layer_1_imputed != TRUE), aes( x = survival_layer_1_obs, y = survival_layer_1_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Survival Layer 1")
survival_layer_2_graph <- ggplot( data = subset( lambir_means, survival_layer_2_imputed != TRUE), aes( x = survival_layer_2_obs, y = survival_layer_2_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Survival Layer 2")
survival_layer_3_graph <- ggplot( data = subset( lambir_means, survival_layer_3_imputed != TRUE), aes( x = survival_layer_3_obs, y = survival_layer_3_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Survival Layer 3")
survival_layer_4_graph <- ggplot( data = subset( lambir_means, survival_layer_4_imputed != TRUE), aes( x = survival_layer_4_obs, y = survival_layer_4_CI90_width ) ) + geom_point() + labs( x = "sample size", y = "90% CI width", title = "Survival Layer 4")
grid.arrange( survival_layer_1_graph,
              survival_layer_2_graph,
              survival_layer_3_graph,
              survival_layer_4_graph )
