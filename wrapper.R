# script to clean all data for the GWD database



# launch one site-specific script at a time -----------

# remove the directory
file_v <- grep( '.R', list.files('R'), value = T)

# read files (this will take a lof of time!)
lapply( file_v, function(x) source(x) )


# Combine data from all sites -------------------------

source( 'R/combine_tables/combine_tables.R' )

