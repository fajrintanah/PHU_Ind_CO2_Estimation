
library(grid) # combining plots
library(gridExtra) # combining plots
library(ggpubr) # combining plots
library(patchwork) # combining plots
library(ggfortify) # nice extension for ggplot
library(tidyverse)
library(data.table)
library(timeplyr)
library(lubridate)
library(zoo)
library(slider)
library(writexl)
library(readxl)
library(hrbrthemes)
library(lubridate)
library(anytime)
library(tidyr)
library(splines)
library(hms)
library(mgcv)
library(ggplot2)
library(viridis)
library(scales)
library(outliers)
library(changepoint)
library(FNN)
library(progress)
library(raster)
library(sf)
library(whitebox)
library(tmap)
library(terra)
library(exactextractr)
library(tools)
library(doParallel)
library(data.table)



# # Sequential resampling function with memory management and flexible resolution

# Load necessary libraries
library(terra)
# Set the temporary directory for R's file usage
tempdir_custom <- "E:/KLHK/KHG/Sumsel/R"
# Ensure the temporary directory exists
if (!dir.exists(tempdir_custom)) {
  dir.create(tempdir_custom, recursive = TRUE)
}
message("Temporary directory set to:", tempdir_custom)

# Sequential resampling function with memory management and flexible resolution
sequential_resample <- function(raster, final_resolution = 10, method = "cubicspline") {
	  current_res <- res(raster)
	  resampled_raster <- raster  
	  # Halve the resolution until the largest side reaches <= 2x final resolution
	  while (max(current_res) > 2 * final_resolution) {
		target_res <- current_res / 2
		temp_target <- rast(ext(raster), res = target_res, crs = crs(raster))
		resampled_raster <- resample(resampled_raster, temp_target, method = method)
		current_res <- res(resampled_raster)
	  }    
	  # Final resampling to the exact final_resolution value
	  final_target <- rast(ext(raster), res = c(final_resolution, final_resolution), crs = crs(raster))
	  resampled_raster <- resample(resampled_raster, final_target, method = method)    
	  return(resampled_raster)
}

# Function to resample, clip, stack, and write raster data
resample_clip_stack_write <- function(
	  raster_data, 
	  output_directory, 
	  tempdir_custom, 
	  final_resolution = 10, 
	  extent = NULL, 
	  method = "cubicspline",
	  stack_name = NULL
	) {
	  # Ensure the output directory exists
	  if (!dir.exists(output_directory)) {
		dir.create(output_directory, recursive = TRUE)
		message("Created directory:", output_directory)
	  } else {
		message("Directory already exists:", output_directory)
	  }  
	  # Set up a temporary directory for intermediate operations
	  tempdir(tempdir_custom)
	  message("Temporary directory set to:", tempdir_custom)  
	  # Initialize an empty list to store resampled and clipped raster layers
	  resampled_layers <- list()  
	  # Iterate over each layer in the raster stack
	  for (i in 1:nlyr(raster_data)) {
		message(paste("Processing layer", i, "of", nlyr(raster_data)))    
		# Access the layer by index
		raster_layer <- raster_data[[i]]    
		# Apply sequential resampling
		resampled_layer <- sequential_resample(raster_layer, final_resolution, method)    
		# Clip the resampled layer if an extent is provided
		if (!is.null(extent)) {
		  message("Clipping layer to extent...")
		  if (inherits(extent, "SpatVector")) {
			resampled_layer <- crop(resampled_layer, extent)
		  } else if (inherits(extent, "SpatRaster")) {
			resampled_layer <- mask(resampled_layer, extent)
		  }
		  message("Clipping complete.")
		}    
		# Store the processed layer
		resampled_layers[[i]] <- resampled_layer
	  }  
	  # Combine all processed layers into a single SpatRaster stack
	  stacked_raster <- rast(resampled_layers)
	  message("All layers combined into a single SpatRaster stack.")  
	  # Determine the output file name
	  output_file <- file.path(
		output_directory, 
		ifelse(is.null(stack_name), "Processed_Stack.tif", paste0(stack_name, ".tif"))
	  )  
	  # Write the resulting stack to disk
	  writeRaster(stacked_raster, output_file, overwrite = TRUE)
	  message("Stacked raster written to:", output_file)  
	  # Return the stacked raster as an in-memory SpatRaster object
	  return(stacked_raster)
	}


# ET 2022 have 22 layers (one month gap: october) which algorithm above mistakenly renamed from 
# "17" "18" "21" "22" "23" "24" to "17" "18" "19" "20" "21" "22" 
# so we will fix it --------------------------------------------------------

# Function to rename raster layers from 19-22 to 21-24
rename_raster_layers_correctly <- function(raster) {
  # Get current names
  original_names <- names(raster)  
  # Define the new names for layers 19 to 22
  new_names <- original_names
  new_names[19] <- "21"
  new_names[20] <- "22"
  new_names[21] <- "23"
  new_names[22] <- "24"  
  # Assign the new names to the raster layers
  names(raster) <- new_names  
  return(raster)
}


# Function to calculate Euclidean distance from the nearest river with flexible input for both river_shp and boundary_extent
calc_euclid_dist <- function(river_shp, boundary_extent, buffer_distance, output_directory, output_raster_name) {  
  # Step 1: Check the type of river_shp and load the river shapefile or raster
  if (inherits(river_shp, "SpatVector")) {
    river_poly <- river_shp
  } else if (dir.exists(river_shp)) {
    shapefiles <- list.files(river_shp, pattern = "\\.shp$", full.names = TRUE)
    if (length(shapefiles) == 0) {
      stop("No shapefile found in the specified river directory.")
    }
    river_poly <- vect(shapefiles[1])
  } else if (file.exists(river_shp) && grepl("\\.shp$", river_shp)) {
    river_poly <- vect(river_shp)
  } else {
    stop("Invalid input for river_shp. It must be a SpatVector, .shp file, or directory containing a shapefile.")
  }  
  
  # Step 2: Check the type of boundary_extent and load the boundary shapefile or raster
  if (inherits(boundary_extent, "SpatVector")) {
    boundary_shp <- boundary_extent
  } else if (inherits(boundary_extent, "SpatRaster")) {
    boundary_shp <- as(boundary_extent, "SpatVector")  # Convert raster to vector
  } else if (dir.exists(boundary_extent)) {
    files <- list.files(boundary_extent, pattern = "\\.shp$|\\.tif$", full.names = TRUE)
    if (length(files) == 0) {
      stop("No shapefile or raster found in the specified boundary directory.")
    }
    if (grepl("\\.shp$", files[1])) {
      boundary_shp <- vect(files[1])
    } else if (grepl("\\.tif$", files[1])) {
      boundary_shp <- rast(files[1])  # If it's a raster
    }
  } else {
    stop("Invalid input for boundary_extent. It must be a SpatVector, SpatRaster, or directory containing a shapefile or raster.")
  }  

  # Step 3: Create a buffer around the river (outside of it)
  river_buffer <- buffer(river_poly, width = buffer_distance)  
  
  # Step 4: Create a raster template based on the boundary extent
  if (inherits(boundary_shp, "SpatVector")) {
    r_template <- rast(ext(boundary_shp), resolution = 10, crs = crs(boundary_shp))
  } else if (inherits(boundary_shp, "SpatRaster")) {
    r_template <- boundary_shp  # Use the boundary raster as template
  }
  
  # Step 5: Rasterize the river polygon (river presence = 1)
  river_raster <- rasterize(river_poly, r_template, field = 1)  
  
  # Step 6: Calculate the Euclidean distance from the river
  distance_raster <- distance(river_raster)  
  
  # Step 7: Mask the distance raster to the boundary extent
  distance_raster_masked <- mask(distance_raster, boundary_shp)  
  
  # Step 8: Save the final distance raster
  output_path <- file.path(output_directory, paste0(output_raster_name, ".tif"))
  writeRaster(distance_raster_masked, output_path, overwrite = TRUE)  
  
  # Return the final distance raster for further use
  return(distance_raster_masked)
}


# Function to calculate aspect in degrees
aspect_deg <- function(dem) {
  # Check input
  if (!inherits(dem, "SpatRaster")) stop("Input must be a 'SpatRaster' object.")  
  # Progress bar
  pb <- progress_bar$new(
    format = "  Calculating aspect [:bar] :percent in :elapsed",
    total = 1, clear = FALSE, width = 60
  )  
  pb$tick(0) # Initialize progress  
  # Calculate aspect in radians
  aspect_radians <- terrain(dem, v = "aspect", unit = "radians")  
  # Convert aspect from radians to degrees
  aspect_deg <- (aspect_radians * 180 / pi) %% 360  
  pb$tick() # Finish progress  
  return(aspect_deg)
}



# Function to convert raster to multipoints (SpatVector) 
raster_to_multipoints <- function(raster, output_crs = NULL) {
  # Ensure the raster is a terra SpatRaster
  if (!inherits(raster, "SpatRaster")) {
    stop("Input must be a terra SpatRaster.")
  }
  
  # Convert raster cells to points, excluding NA values
  points <- terra::as.points(raster, na.rm = TRUE)
  
  # Transform CRS if output_crs is provided
  if (!is.null(output_crs)) {
    points <- terra::project(points, output_crs)
  }
  
  return(points)  # Returns a SpatVector
}


# Function to extract raster values and add them as attributes to a SpatVector
extract_raster_values <- function(raster, points_sv, attribute_name = "value") {
  # Ensure the input raster is a terra SpatRaster
  if (!inherits(raster, "SpatRaster")) {
    stop("Input must be a terra SpatRaster.")
  }
  
  # Ensure the input points are a terra SpatVector
  if (!inherits(points_sv, "SpatVector")) {
    stop("Input points must be a terra SpatVector.")
  }
  
  # Check if the CRS of the points matches the raster
  if (!terra::same.crs(raster, points_sv)) {
    points_sv <- terra::project(points_sv, terra::crs(raster))
    message("Points have been reprojected to match the raster CRS.")
  }
  
  # Extract raster values at the point locations
  values <- terra::extract(raster, points_sv)
  
  # Add the extracted values as a new attribute to the SpatVector
  points_sv[[attribute_name]] <- values[, 2]  # The second column contains the extracted values

  return(points_sv)  # Returns the updated SpatVector
}


# Function to assign values from points to a raster
assign_values_to_raster <- function(fishnet_rast, point_sf, column_name) {
  # Ensure the input raster is a terra SpatRaster
  if (!inherits(fishnet_rast, "SpatRaster")) {
    stop("Input must be a terra SpatRaster.")
  }  
  # Ensure the points_sf is a simple features (sf) object
  if (!inherits(point_sf, "sf")) {
    stop("Points input must be a simple features (sf) object.")
  }  
  # Ensure the column exists in the sf object
  if (!column_name %in% colnames(point_sf)) {
    stop(paste("Column", column_name, "does not exist in the points object."))
  }  
  # Convert the sf object to a SpatVector
  point_vect <- vect(point_sf)  
  # Rasterize the points, assigning the specified column's values to the raster
  raster_with_values <- rasterize(point_vect, fishnet_rast, field = column_name)  
  return(raster_with_values)
}


# function to clip all resulted rasters iteratively in the directory to eliminate pixel outside of base raster extent
# build wrapper function
clip_raster <- function(input_directory, clipping_extent, output_directory, pattern = "\\.tif$") {
  # Ensure output directory exists
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }    
  # List raster files in the input directory
  raster_files <- list.files(input_directory, pattern = pattern, full.names = TRUE)    
  # Check if files are found
  if (length(raster_files) == 0) {
    stop("No raster files found in the input directory with the given pattern!")
  }    
  # Extract EPSG code of the clipping_extent
  clipping_epsg <- crs(clipping_extent, proj = TRUE)
  clipping_epsg <- sub(".*EPSG:", "", clipping_epsg)  # Extract the EPSG code    
  # Initialize progress bar
  pb <- progress_bar$new(
    format = "[:bar] :percent :current/:total (:elapsed | ETA: :eta)",
    total = length(raster_files),
    clear = FALSE,
    width = 60
  )    
  # Loop through each raster file
  for (raster_path in raster_files) {
    # Load the raster
    input_raster <- rast(raster_path)      
    # Extract EPSG code of the input raster
    input_epsg <- crs(input_raster, proj = TRUE)
    input_epsg <- sub(".*EPSG:", "", input_epsg)  # Extract the EPSG code        
    # Check if EPSG codes match
    if (input_epsg != clipping_epsg) {
      warning(
        paste0(
          "EPSG mismatch detected for raster: ", raster_path, 
          "\nInput EPSG: ", input_epsg, 
          "\nClipping Extent EPSG: ", clipping_epsg,
          "\nConsider reprojecting the raster to match the clipping extent EPSG."
        )
      )
    }       
    # Perform mask and crop
    processed_raster <- mask(crop(input_raster, clipping_extent), clipping_extent)        
    # Generate output file path
    raster_name <- basename(raster_path)  # Extract the file name
    output_path <- file.path(output_directory, raster_name)        
    # Save the processed raster
    writeRaster(processed_raster, output_path, overwrite = TRUE)        
    # Update progress bar
    pb$tick()
  }    
  message("All files processed successfully!")
}


# Function to extract values from a SpatRaster stack and add them as attributes to an sf object
extract_raster_stack_values <- function(raster_stack, points_sf, attribute_name) {
  # Ensure the input is a SpatRaster with multiple layers
  if (!inherits(raster_stack, "SpatRaster")) {
    stop("Input must be a terra SpatRaster.")
  }
  if (terra::nlyr(raster_stack) < 1) {
    stop("The SpatRaster stack must have at least one layer.")
  }
  # Ensure the CRS of points matches the raster stack
  if (!st_crs(points_sf)$wkt == crs(raster_stack)) {
    stop("CRS of points and raster stack do not match. Reproject points to match raster stack.")
  }
  # Extract raster values for all layers
  point_coords <- st_coordinates(points_sf)  # Extract coordinates from sf object
  values <- terra::extract(raster_stack, point_coords)  # Extract values  
  # Ensure the extraction returned results
  if (is.null(values)) {
    stop("Extraction failed. Check raster extent and point coordinates.")
  }
  # Remove the first column if it contains point IDs (not relevant here, but safe to check)
  if (ncol(values) == terra::nlyr(raster_stack) + 1) {
    values <- values[, -1, drop = FALSE]
  }
  # Set column names for the values data frame
  layer_names <- names(raster_stack)  # Get layer names from the raster stack
  if (is.null(layer_names) || all(layer_names == "")) {
    # Default to numeric indices if layer names are missing
    layer_names <- seq_len(terra::nlyr(raster_stack))
  }
  column_names <- paste0(attribute_name, layer_names)  # Combine attribute name with layer names
  colnames(values) <- column_names
  # Add extracted values to the points_sf object
  for (i in seq_along(column_names)) {
    points_sf[[column_names[i]]] <- values[[i]]
  }  
  return(points_sf)
}

# function to clip all  rasters iteratively in the directory to eliminate pixel outside of base raster extent
# build wrapper function
clip_raster <- function(input_directory, clipping_extent, output_directory, pattern = "\\.tif$") {
  # Ensure output directory exists
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }    
  # List raster files in the input directory
  raster_files <- list.files(input_directory, pattern = pattern, full.names = TRUE)    
  # Check if files are found
  if (length(raster_files) == 0) {
    stop("No raster files found in the input directory with the given pattern!")
  }    
  # Extract EPSG code of the clipping_extent
  clipping_epsg <- crs(clipping_extent, proj = TRUE)
  clipping_epsg <- sub(".*EPSG:", "", clipping_epsg)  # Extract the EPSG code    
  # Initialize progress bar
  pb <- progress_bar$new(
    format = "[:bar] :percent :current/:total (:elapsed | ETA: :eta)",
    total = length(raster_files),
    clear = FALSE,
    width = 60
  )    
  # Loop through each raster file
  for (raster_path in raster_files) {
    # Load the raster
    input_raster <- rast(raster_path)      
    # Extract EPSG code of the input raster
    input_epsg <- crs(input_raster, proj = TRUE)
    input_epsg <- sub(".*EPSG:", "", input_epsg)  # Extract the EPSG code        
    # Check if EPSG codes match
    if (input_epsg != clipping_epsg) {
      warning(
        paste0(
          "EPSG mismatch detected for raster: ", raster_path, 
          "\nInput EPSG: ", input_epsg, 
          "\nClipping Extent EPSG: ", clipping_epsg,
          "\nConsider reprojecting the raster to match the clipping extent EPSG."
        )
      )
    }       
    # Perform mask and crop
    processed_raster <- mask(crop(input_raster, clipping_extent), clipping_extent)        
    # Generate output file path
    raster_name <- basename(raster_path)  # Extract the file name
    output_path <- file.path(output_directory, raster_name)        
    # Save the processed raster
    writeRaster(processed_raster, output_path, overwrite = TRUE)        
    # Update progress bar
    pb$tick()
  }    
  message("All files processed successfully!")
}


# Function to extract values from a SpatRaster stack and add them as attributes to eihter spatvector or sf objects
extract_raster_stack_values <- function(raster_stack, points, attribute_name) {
  # Ensure the input is a SpatRaster with multiple layers
  if (!inherits(raster_stack, "SpatRaster")) {
    stop("Input must be a terra SpatRaster.")
  }
  
  if (terra::nlyr(raster_stack) < 1) {
    stop("The SpatRaster stack must have at least one layer.")
  }

  # Ensure the input is a SpatVector
  if (!inherits(points, "SpatVector")) {
    stop("Points input must be a terra SpatVector.")
  }

  # Check CRS compatibility and reproject points if necessary
  if (!terra::same.crs(raster_stack, points)) {
    message("CRS mismatch detected. Reprojecting points to match raster CRS.")
    points <- terra::project(points, crs(raster_stack))
  }

  # Align extents by cropping the raster stack to the vector's extent
  raster_stack_cropped <- terra::crop(raster_stack, terra::ext(points))

  # Check if raster and vector extents intersect
  if (is.null(terra::intersect(terra::ext(raster_stack_cropped), terra::ext(points)))) {
    stop("The extents of the raster and vector do not overlap.")
  }

  # Extract raster values for all layers at the point locations
  values <- terra::extract(raster_stack_cropped, points)

  # Ensure the extraction returned results
  if (is.null(values)) {
    stop("Extraction failed. Check raster extent and point coordinates.")
  }

  # Remove the first column if it contains point IDs (depends on terra behavior)
  if (ncol(values) == terra::nlyr(raster_stack) + 1) {
    values <- values[, -1, drop = FALSE]
  }

  # Set column names for the extracted values
  layer_names <- names(raster_stack_cropped)
  if (is.null(layer_names) || all(layer_names == "")) {
    layer_names <- seq_len(terra::nlyr(raster_stack_cropped))
  }

  column_names <- paste0(attribute_name, layer_names)
  colnames(values) <- column_names

  # Add extracted values as new attributes in the original points SpatVector
  for (i in seq_along(column_names)) {
    points[[column_names[i]]] <- values[[i]]
  }

  return(points)
}



# Function to re-align the stacks by creating new input raster from fishnet point vector --------------- 
# a wrapper function to assign multiple columns to rasters and create a raster stack with sequential names 
# and also crop all resulting rasters within the stacks according to the base raster
# consume huge memory. Do it if you have at least 32GB RAM with around 60-70% free 

assign_values_to_temp_stack <- function(fishnet_rast, points, column_pattern, temp_dir, stack_name, output_file = NULL, input_extent = NULL) {
  # Ensure the input raster is a terra SpatRaster
  if (!inherits(fishnet_rast, "SpatRaster")) {
    stop("Input raster must be a terra SpatRaster.")
  }
  
  # Check if the points are a SpatVector or sf object
  if (inherits(points, "sf")) {
    # Convert sf to SpatVector
    points <- vect(points)
  } else if (!inherits(points, "SpatVector")) {
    stop("Points input must be either an sf or SpatVector object.")
  }
  
  # Identify all columns in points matching the column pattern
	matching_columns <- grep(column_pattern, terra::names(points), value = TRUE)
	if (length(matching_columns) == 0) {
	stop(paste("No columns match the pattern:", column_pattern))
	}
  
  # Create a temporary directory for storing raster layers
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE)
  }
  
  # Check CRS and reproject if necessary
  if (!terra::same.crs(fishnet_rast, points)) {
    cat("CRS mismatch detected. Reprojecting points to match raster CRS...\n")
    points <- terra::project(points, crs(fishnet_rast))
  }
  
  # Apply cropping and masking if input_extent is provided
  if (!is.null(input_extent)) {
    fishnet_rast <- terra::mask(terra::crop(fishnet_rast, input_extent), input_extent)
  }
  
  # Iterate over each matching column
  temp_files <- character(length(matching_columns))  # Store file paths for later stacking
  for (i in seq_along(matching_columns)) {
    col_name <- matching_columns[i]
    
    # Rasterize the points using the current column
    raster_layer <- terra::rasterize(points, fishnet_rast, field = col_name)
    
    # Save the raster layer to a temporary file
    temp_file <- file.path(temp_dir, paste0(stack_name, "_", i, ".tif"))
    terra::writeRaster(raster_layer, temp_file, overwrite = TRUE)
    temp_files[i] <- temp_file  # Store the file path
    
    # Print progress
    cat("Processed and saved raster layer:", col_name, "to", temp_file, "\n")
  }
  
  # Load all temporary files into a raster stack
  raster_stack <- terra::rast(temp_files)
  
  # Assign numeric names (1:24) to the layers
  names(raster_stack) <- as.character(seq_len(terra::nlyr(raster_stack)))
  cat("All layers processed and stacked.\n")
  
  # Apply cropping and masking to the stack if input_extent is provided
  if (!is.null(input_extent)) {
    cat("Cropping and masking the raster stack to the provided extent...\n")
    raster_stack <- terra::mask(terra::crop(raster_stack, input_extent), input_extent)
    cat("Cropping and masking completed.\n")
  }
  
  # Write the final stack to the specified output file if provided
  if (!is.null(output_file)) {
    terra::writeRaster(raster_stack, output_file, overwrite = TRUE)
    cat("Final raster stack saved to:", output_file, "\n")
  }
  
  return(raster_stack)
}

# Wrapper function to compute D-Infinity Flow Accumulation after filling depressions
comp_d_inf_facc <- function(
  Elev_raster,                   # Input DEM as a SpatRaster or file path
  output_dir,                    # Output directory
  save_fill = TRUE,              # Save the filled DEM or not
  output_fill_name = "filled_dem.tif",      # Output name for filled DEM
  output_facc_name = "d_inf_flow_acc.tif"   # Output name for flow accumulation raster
) {
  # Ensure the output directory ends with a slash
  output_dir <- ifelse(substr(output_dir, nchar(output_dir), nchar(output_dir)) == "/", 
                       output_dir, paste0(output_dir, "/"))

  # Internal output file paths
  internal_filled_path <- file.path(output_dir, output_fill_name)
  internal_facc_path <- file.path(output_dir, output_facc_name)
  # Handle input raster (SpatRaster or file path)
  if (inherits(Elev_raster, "SpatRaster")) {
    temp_elev_path <- file.path(tempdir(), "temp_elev.tif")
    writeRaster(Elev_raster, temp_elev_path, overwrite = TRUE)
  } else {
    temp_elev_path <- Elev_raster
  }
  # Step 1: Fill depressions
  cat("Filling depressions...\n")
  whitebox::wbt_fill_depressions(
    dem = temp_elev_path,
    output = internal_filled_path
  )
  # Step 2: Compute D-Infinity Flow Accumulation
  cat("Computing D-Infinity Flow Accumulation...\n")
  whitebox::wbt_d_inf_flow_accumulation(
    input = internal_filled_path,
    output = internal_facc_path
  )
  # Optionally, delete the filled DEM if save_fill is FALSE
  if (!save_fill) {
    file.remove(internal_filled_path)
    cat("Filled DEM file removed.\n")
  } else {
    cat("Filled DEM saved as:", internal_filled_path, "\n")
  }
  # Return SpatRaster objects as a list with dynamic names
  result <- list()
  if (save_fill) {
    result[[output_fill_name]] <- rast(internal_filled_path)  # Assign filled DEM
  }
  result[[output_facc_name]] <- rast(internal_facc_path)  # Assign flow accumulation
  cat("Flow accumulation file saved as:", internal_facc_path, "\n")
  return(result)
}

# function to compute TWI

compute_twi <- function(sca, slope, output_dir, output_twi_name) {   # compute TWI using custom wrapper function
  # Define the output file path
  output_twi_path <- file.path(output_dir, output_twi_name)  
  # Run the WhiteboxTools wetness index function
  wbt_wetness_index(
    sca = sca,           # Flow Accumulation raster
    slope = slope,       # Slope raster
    output = output_twi_path,  # Output file path
    wd = NULL,
    verbose_mode = NULL,
    compress_rasters = NULL,
    command_only = FALSE
  )  
  # Load and return the output TWI raster as a SpatRaster object
  return(rast(output_twi_path))
}



#  function to make subsurface elevation ---------------------

subtract_rasters <- function(input_raster, subtraction_raster, output_raster_name, output_directory) {  
  # Step 1: Extract EPSG codes from both input rasters
  input_epsg <- crs(input_raster, proj = TRUE)
  input_epsg <- sub(".*EPSG:", "", input_epsg)  # Extract the EPSG code for input raster  
  subtraction_epsg <- crs(subtraction_raster, proj = TRUE)
  subtraction_epsg <- sub(".*EPSG:", "", subtraction_epsg)  # Extract the EPSG code for subtraction raster  
  # Step 2: Check if EPSG codes of both rasters match
  if (input_epsg != subtraction_epsg) {
    warning(
      paste0(
        "EPSG mismatch detected:\n",
        "Input raster EPSG: ", input_epsg, "\n",
        "Subtraction raster EPSG: ", subtraction_epsg, "\n",
        "Consider reprojecting the rasters to match the EPSG codes."
      )
    )
  }  
  # Step 3: Mask the rasters to exclude NAs in the overlapping pixels
  input_raster_masked <- mask(input_raster, subtraction_raster, maskvalues = NA)
  subtraction_raster_masked <- mask(subtraction_raster, input_raster, maskvalues = NA)  
  # Step 4: Perform the pixel-wise subtraction
  result_raster <- input_raster_masked - subtraction_raster_masked  
  # Step 5: Write the resulting raster to the specified output directory
  output_path <- file.path(output_directory, paste0(output_raster_name, ".tif"))
  writeRaster(result_raster, output_path, overwrite = TRUE)  
  # Return the result raster object
  return(result_raster)
}

# Calculate slope as percent rise
slope_pct <- function(dem) {
  # Check input
  if (!inherits(dem, "SpatRaster")) stop("Input must be a 'SpatRaster' object.")  
  # Progress bar
  pb <- progress_bar$new(
    format = "  Calculating slope [:bar] :percent in :elapsed",
    total = 1, clear = FALSE, width = 60
  )  
  pb$tick(0) # Initialize progress  
  # Calculate slope in radians
  slope_radians <- terrain(dem, v = "slope", unit = "radians")  
  # Convert radians to percent rise
  slope_percent <- tan(slope_radians) * 100  
  pb$tick() # Finish progress  
  return(slope_percent)
}


# Function to calculate aspect in degrees
aspect_deg <- function(dem) {
  # Check input
  if (!inherits(dem, "SpatRaster")) stop("Input must be a 'SpatRaster' object.")  
  # Progress bar
  pb <- progress_bar$new(
    format = "  Calculating aspect [:bar] :percent in :elapsed",
    total = 1, clear = FALSE, width = 60
  )  
  pb$tick(0) # Initialize progress  
  # Calculate aspect in radians
  aspect_radians <- terrain(dem, v = "aspect", unit = "radians")  
  # Convert aspect from radians to degrees
  aspect_deg <- (aspect_radians * 180 / pi) %% 360  
  pb$tick() # Finish progress  
  return(aspect_deg)
}

# Function to rename and select columns from a SpatVector
rename_and_select_column <- function(spatvector, old_name, new_name, selected_column) {  
  # Rename the column
  names(spatvector)[names(spatvector) == old_name] <- new_name  
  # Select the specified column
  spatvector <- spatvector[, selected_column, drop = FALSE]  
  return(spatvector)
}

# Function to combine multiple stacks into a single stack 
# All layers are sequentially named 1:72 to match Rst_ID of GWL obs dataframe
combine_and_rename <- function(..., output_directory, output_raster_name) {
  # Capture all input SpatRaster objects into a list
  rasters <- list(...)  
  # Total number of input rasters and layers for progress
  total_layers <- sum(sapply(rasters, nlyr))
  
  # Create a progress bar
  pb <- progress_bar$new(
    format = "  Combining Rasters [:bar] :percent in :elapsed",
    total = total_layers,
    clear = FALSE,
    width = 60
  )  
  # Initialize an empty SpatRaster
  combined_raster <- rast()  
  # Loop through each raster to combine and track progress
  for (i in seq_along(rasters)) {
    combined_raster <- c(combined_raster, rasters[[i]])
    
    # Update progress bar after processing each layer
    for (layer in 1:nlyr(rasters[[i]])) {
      pb$tick()
    }
  }  
  # Rename layers sequentially
  names(combined_raster) <- seq(1, total_layers)  
  # Create the output file path
  output_path <- file.path(output_directory, paste0(output_raster_name, ".tif"))  
  # Save the combined raster to the specified directory
  writeRaster(combined_raster, output_path, overwrite = TRUE)  
  # Return the combined raster stack
  return(combined_raster)
}


# function to import SpatRaster objects. I used spatraster list not stack
import_raster_list <- function(directory, file_extension = "\\.tif$") {
  # List all raster files in the specified directory with the given extension
  raster_files <- list.files(path = directory, pattern = file_extension, full.names = TRUE)
  
  # Check if any raster files were found
  if (length(raster_files) == 0) {
    stop("No raster files found in the specified directory with the given extension.")
  }
  
  # Read the rasters into a list of SpatRaster objects
  raster_list <- lapply(raster_files, rast)
  
  # Return the list of SpatRaster objects
  return(raster_list)
}


Sumsel_static_var <- import_raster_list("D:/Sumsel/3_olah5/Static_var")

# Loop through each raster in the list
for (i in 1:length(Sumsel_static_var)) {
  
  # Get the raster layer
  raster_layer <- Sumsel_static_var[[i]]
  
  # Extract the file path using sources() method
  file_path <- sources(raster_layer)
  
  # Get the file name without the extension
  raster_name <- file_path_sans_ext(basename(file_path))
  
  # Set the raster's name to match the file name
  names(raster_layer) <- raster_name
  
  # Update the list with the modified raster layer
  Sumsel_static_var[[i]] <- raster_layer
  
  # Optionally, print each raster to check its details
  print(raster_layer)
}


# function to extract raster values, clean column names, and remove "ID" columns (except point_ID)
extract_rename_static <- function(raster_list, point_data) {
  # Extract values from rasters for each point
  extracted_values_list <- lapply(raster_list, function(raster) {
    extract(raster, point_data)
  })  
  # Convert the list to a data frame
  extracted_values_df <- as.data.frame(extracted_values_list)  
  # Add the point_ID from the point data
  extracted_values_df$point_ID <- point_data$point_ID  
  # Rename columns by removing "Sumsel_" from the start and "_Sumsel" from the end
  colnames(extracted_values_df) <- gsub("^Sumsel_", "", colnames(extracted_values_df))   # Remove "Sumsel_" from the start
  colnames(extracted_values_df) <- gsub("_Sumsel$", "", colnames(extracted_values_df))   # Remove "_Sumsel" from the end  
  # Remove columns that contain "ID", except for the "point_ID" column
  extracted_values_df_clean <- extracted_values_df[, !grepl("ID", colnames(extracted_values_df)) | colnames(extracted_values_df) == "point_ID"]
  # Return the cleaned data frame
  return(extracted_values_df_clean)
}


