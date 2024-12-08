
library(terra)

# function to process and calculate CO2 from multiple Rasters -------------------------------

compute_CO2_flux <- function(input_GWL_raster, input_land_use_raster, computed_number, 
                             output_CO2_directory, output_CO2_name, clipping_extent) {
  # Helper function to compute CO2 flux for a given GWL and land-use raster
  compute_CO2_flux_single <- function(gwl_raster, land_use_raster) {
    # CO2 flux calculations based on land use and GWL
    flux_pf <- exp((3.61 / (1 + exp(30.7 * ((gwl_raster / 100) - 0.01)))) + 
                     (0.33 / (1 + exp(103 * ((gwl_raster / 100) + 0.06)))) - 1.34)
    
    # Agriculture flux calculation, with valid GWL range condition applied
    flux_agriculture <- 68.64 - (196.80 * gwl_raster) - (75.97 * (gwl_raster^2))
    flux_agriculture <- mask(flux_agriculture, gwl_raster >= -90 & gwl_raster <= 0)
    
    # Plantation flux calculation
    flux_plantation <- 48.8052145 * exp(-exp(-0.0002404 * (gwl_raster - 0.9002186)))
    flux_plantation <- mask(flux_plantation, gwl_raster >= -80 & gwl_raster <= 0)
    
    # IPF flux calculation
    flux_IPF <- 86.21 - 53.22 * gwl_raster + 142.65 * (gwl_raster^2)
    flux_IPF <- mask(flux_IPF, gwl_raster >= -110 & gwl_raster <= -50)
    
    # SF flux calculation
    flux_SF <- 50.16 - 58.95 * gwl_raster + 42.05 * (gwl_raster^2)
    flux_SF <- mask(flux_SF, gwl_raster >= -80 & gwl_raster <= -30)
    
    # Global flux formula
    flux_global <- 117.06 + (-70.91) / (1 + exp((-0.10) * (gwl_raster + 75.78)))
    
    # Default flux (global) and then apply based on land-use type
    flux <- flux_global  # Default flux if no condition matches
    
    # Apply land-use conditions
    flux <- mask(flux, land_use_raster == 2001 | land_use_raster == 2005)  # For Primary Forest (PF)
    flux <- flux + flux_pf * (land_use_raster == 2001 | land_use_raster == 2005)
    
    flux <- mask(flux, land_use_raster == 2007 | land_use_raster == 3000 | 
                           land_use_raster == 20071 | land_use_raster == 20091 | 
                           land_use_raster == 20092 | land_use_raster == 20093)  # Agriculture (Agr)
    flux <- flux + flux_agriculture * (land_use_raster == 2007 | land_use_raster == 3000 | 
                                        land_use_raster == 20071 | land_use_raster == 20091 | 
                                        land_use_raster == 20092 | land_use_raster == 20093)
    
    flux <- mask(flux, land_use_raster == 2010)  # Plantation (PL)
    flux <- flux + flux_plantation * (land_use_raster == 2010)
    
    flux <- mask(flux, land_use_raster == 2006)  # IPF (Industrial Plantation Forest)
    flux <- flux + flux_IPF * (land_use_raster == 2006)
    
    flux <- mask(flux, land_use_raster == 2002 | land_use_raster == 20051)  # Secondary Forest (SF)
    flux <- flux + flux_SF * (land_use_raster == 2002 | land_use_raster == 20051)
    
    return(flux)
  }
  
  # Loop over the computed number of GWL rasters
  for (i in computed_number) {
    # Determine the year based on the index
    if (i <= 24) {
      year <- 2021
    } else if (i <= 48) {
      year <- 2022
    } else {
      year <- 2023
    }
  
    # Construct file paths for the current GWL and land-use raster
    gwl_raster_path <- file.path(input_GWL_raster, paste0("predicted_GWL_input_GWL_", i, ".tif"))
    land_use_raster_path <- file.path(input_land_use_raster, paste0("Sumsel_LU_", year, "_input.tif"))
    
    # Load GWL and land-use rasters
    gwl_raster <- rast(gwl_raster_path)
    land_use_raster <- rast(land_use_raster_path)
    
    # Resample rasters to match clipping_extent
    gwl_raster_resampled <- resample(gwl_raster, clipping_extent, method = "bilinear")
    land_use_raster_resampled <- resample(land_use_raster, clipping_extent, method = "near")
    
    # Crop and mask rasters based on clipping_extent
    gwl_raster_clipped <- mask(crop(gwl_raster_resampled, clipping_extent), clipping_extent)
    land_use_raster_clipped <- mask(crop(land_use_raster_resampled, clipping_extent), clipping_extent)
    
    # Handle edge effects: exclude pixels where either raster has NA
    valid_mask <- !is.na(gwl_raster_clipped) & !is.na(land_use_raster_clipped)
    gwl_raster_valid <- mask(gwl_raster_clipped, valid_mask)
    land_use_raster_valid <- mask(land_use_raster_clipped, valid_mask)
    
    # Compute CO2 flux
    flux_raster <- compute_CO2_flux_single(gwl_raster_valid, land_use_raster_valid)
    
    # Save the computed CO2 flux raster directly to a file
    output_file <- file.path(output_CO2_directory, paste0(output_CO2_name, "_", i, ".tif"))
    writeRaster(flux_raster, output_file, overwrite = TRUE)
    
    # Explicitly remove large objects to free up memory
    rm(gwl_raster, land_use_raster, gwl_raster_resampled, land_use_raster_resampled, 
       gwl_raster_clipped, land_use_raster_clipped, gwl_raster_valid, land_use_raster_valid, flux_raster)
    
    # Invoke garbage collection to free memory
    gc()
    
    # Optionally print progress
    cat("Processed and saved CO2 flux raster for GWL_", i, "\n")
  }
}
