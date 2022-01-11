library(tidyverse)

# taxa table
taxa_all    <- list.files('results/') %>% 
                  grep( 'taxa', ., value = T )
taxa_ures   <- grep( 'unresvd.csv', taxa_all, value = T )
taxa_tabs   <- setdiff( taxa_all, taxa_ures )
taxa_df     <- lapply( paste0('results/',taxa_tabs), 
                       read.csv ) %>% 
                  bind_rows %>% 
                  rename( Site_name = site )

# site table 
site_tabs   <- list.files('results/') %>% 
                  grep( '_site', ., value = T )
site_df     <- lapply( paste0('results/',site_tabs), 
                       read.csv ) %>% 
                  bind_rows %>% 
                  mutate( Site_name = tolower(Site_name) )

# means table
means_all   <- list.files('results/') %>% 
                  grep( '_means', ., value = T )
means_unres <- grep( 'unresvd.csv', means_all, value = T )
means_tabs  <- setdiff( means_all, means_unres )
means_df    <- lapply( paste0('results/',means_tabs), 
                       read.csv ) %>% 
                  bind_rows %>% 
                  rename( Site_name = site )


write.csv( site_df,  'results/full_db/site_table.csv', row.names = F )
write.csv( taxa_df,  'results/full_db/taxa_table.csv', row.names = F )
write.csv( means_df, 'results/full_db/means_table.csv', row.names = F )

# Metadata files ------------------------------------------------------------------

# store variables tables
taxa_meta <- data.frame( variable = taxa_df %>% names,
                         meaning  = 
                         c("Original name contained in Kambach's file",
                           "Closes matches of 'Submitted_Name' in the LCVP obtained via fuzzy matching",
                           "Nomenclature status of the species name contained in First_matched_Name. Possible values: 'accepted', 'synonym', 'unresolved' or 'external'",
                           "This field provides a direct comparison with ‘The Plant List’. 'The Plant list' is an alternative taxonomy to the Leipzig Catalogue of Plants. This taxonomy is well known, and it is found at http://www.theplantlist.org/. Possible values: ‘identical', 'synonym', 'other synonym', 'different authors', 'missing', 'misspelling' or 'unresolved'",
                           "This field provides a possible alternative name from ‘The Plant List’",
                           "The accepted plant taxa name according to the LCVP",
                           "Self-explanatory",
                           "Self-explanatory",
                           "Was 'Submitted_Name' equal to 'First_matched_Name'? 'FALSE' 1. certifies a mismatch, and 2. indicates that we left this taxon in the DB because we subjectively determined that the mismatch was a typo which could be resolved by the LCVP R package",
                           "Code referring to the ForestGeo site. This is a foreign key that links to the site table" )
                         )


# site metadata
site_meta <- data.frame( variable = site_df %>% names,
                         meaning  = c('Site name',
                                      'Self-explanatory',
                                      'Self-explanatory')
)

# means metadata are done by hand!


# write it down
write.csv( site_meta, 'results/full_db/metadata/site_tab_metadata.csv', row.names = F )
write.csv( taxa_meta, 'results/full_db/metadata/taxa_tab_metadata.csv', row.names = F )


# MERGE ------------------------------------------------------------------

# merge all dataset
all_df <- site_df %>% 
            left_join( taxa_df ) %>% 
            left_join( means_df) 

# write the whole database in a single file
write.csv( all_df, 'results/full_db/merged_gwd.csv',
           row.names = F )
