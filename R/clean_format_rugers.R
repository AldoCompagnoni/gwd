library(tidyverse)
# devtools::install_github("idiv-biodiversity/LCVP")
# devtools::install_github("idiv-biodiversity/lcvplants")
library(LCVP)
library(lcvplants)
library(leaflet)

<<<<<<< HEAD
# read in data
=======
# read in data (Appendix from dx.doi.org/10.1126/science.aaz4797 )
>>>>>>> d4ef71d195a46b3147d9f62a1f9cdc2fd3a74238
rug1 <- readxl::read_xlsx( 'data/aaz4797_Ruger_Data_S1.xlsx', 
                           sheet = 2 )


# Prepare Site table -------------------------------

# Do this by hand (I could not find this data in tabular form yet)
site_out <- data.frame( Site_name = 'BCI',
                        Latitude  = 9.154300000000,
                        Longitude = -79.846100000000 )

# Check that the location is sensible sense
leaflet( data = site_out ) %>% 
  addTiles() %>% 
  addCircleMarkers( ~Longitude, ~Latitude )
  
# store site information
write.csv( site_out, 'results/bci_site.csv',
           row.names = F )


# Prepare taxonomic table ------------------------------

# Produce the binomial used for checking
taxa_df_rugers         <- dplyr::select( rug1, Genus, Species ) %>% 
                     mutate( Submitted_Name = paste0( Genus, ' ', Species) ) %>% 
                     rename( Genus_author   = Genus,
                             Species_author = Species )


# function: get "cleaned" names
get_clean_names <- function( nam, fuzzy = 0 ) lcvp_search( nam, max.distance = fuzzy )

# Clean names from the Leipzig's list of plants
clean_l_rugers         <- lapply( taxa_df_rugers$Submitted_Name , get_clean_names )
clean_df_rugers        <- clean_l_rugers %>% bind_rows

# check whether there are issues
notfound_df_rugers     <- clean_df_rugers %>% subset( PL.comparison != 'identical' )

# Rerun Leipzig list with fuzzy matching
reclean_l_rugers       <- lapply( notfound_df_rugers$Submitted_Name, get_clean_names, 5 )
reclean_df_rugers      <- reclean_l_rugers %>% bind_rows

# Examine final file
reclean_df_rugers      

# Upon scrituny
# 1. There are several spelling mistakes
# 2. We cannot determine the genus for only two species:
# Nectandra and Sapium.

# Final taxonomy files 
taxa_nofuzzy    <- clean_df %>% subset( Score == 'matched' ) 
taxa_fuzzy      <- count( reclean_df, Submitted_Name ) %>% 
                     subset( n > 1 ) %>% 
                     dplyr::select( -n )
taxa_unresvd    <- count( reclean_df, Submitted_Name ) %>% 
                     subset( n > 2 ) %>% 
                     dplyr::select( -n )
taxa_out        <- bind_rows( taxa_nofuzzy, taxa_fuzzy) %>% 
                    left_join( taxa_df )

# store resolved AND unresolved taxa
write.csv( taxa_out, 'results/bci_taxa.csv',
           row.names = F )
write.csv( taxa_unresvd, 'results/bci_taxa_unresvd.csv',
           row.names = F )


# Prepare demographic table --------------------------------------

# Are there missing/NULL/Inf values?
# Make some box plots to explore the data
# Check for sample size 0

# check survival
dplyr::select(rug1, mu1_1:mu4_1 ) %>% 
  lapply( range )

# check growth
dplyr::select(rug1, G1_1:G4_1 ) %>% 
  lapply( range )

# Check fecundity
rug1$F_1 %>% range

# Check fecundity
rug1$SampleSize %>% range
rug1$SampleSize %>% hist(breaks=30)

# Taxa information for merging
taxa_merge_df <- dplyr::select( taxa_out, Submitted_Name, Genus, Species )

# Final "demography table"
demo_out <- dplyr::select( rug1, 
                           Genus, Species,
                           mu1_1:mu4_1, 
                           G1_1:G4_1, F_1, 
                           SampleSize ) %>% 
              # service column for merging with taxonomic information
              mutate( Submitted_Name = paste0( Genus, ' ', Species) ) %>% 
              dplyr::select( -Genus, -Species ) %>% 
              # right join to include only taxa that we identified with LCVP
              right_join( taxa_merge_df ) %>% 
              # add site name (BCI: Barro Colorado Island)
              mutate( Site_name = 'BCI' )
  
# store demographic information
write.csv( demo_out, 'results/bci_demog.csv', row.names = F )
