
library(tidyverse) # plotting and manipulation
library(grid) # combining plots
library(gridExtra) # combining plots
library(ggpubr) # combining plots
library(patchwork) # combining plots
library(ggfortify) # nice extension for ggplot
library(mgcv) #fitting gam models
library(GGally) # displaying pairs panel
library(caret)
library(caTools) # split dataset
library(readxl)
library(randomForest)
library(e1071)
library(gbm)          # basic implementation
library(xgboost)      # a faster implementation of gbm
library(caret)        # an aggregator package for performing many machine learning models
library(pdp)          # model visualization
library(lime)         # model visualization
library(neuralnet)
library(rpart)     #rpart for computing decision tree models
library(rsample)     # data splitting 
library(dplyr)       # data wrangling
library(rpart.plot)  # plotting regression trees
library(ipred)       # bagging
library(broom)
library(ranger) 	#efficient RF
library(NeuralNetTools)
library(tidymodels)
library(earth) 		#MARS model
library(iml)		#most robust and efficient relative importance 
library(xgboost)	#extreeme gradient boosting
library(ModelMetrics) #get model metrics
library(Metrics) 	#get ML model metrics
library(Cubist) #Cubist modell
library(iBreakDown)
library(DALEX)
library(viridis)
library(ICEbox)
library(hrbrthemes)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(plotly)
library(vip)
library(fastDummies)
library(fastDummies)

#  function to evaluate predictions
evaluate_and_bind <- function(predictions, y_test) {
  models <- names(predictions)
  
  eval_metrics <- map(models, function(model) {
    pred <- predictions[[model]]
    
    # Check for NA in the predictions and skip faulty models if needed
    if (anyNA(pred)) {
      warning(paste("NA values detected in", model, "predictions. Skipping."))
      return(NULL)
    }
    
    data.frame(
      Rsq = R2(pred, y_test),              # R-squared
      RMSE = RMSE(pred, y_test),           # Root Mean Square Error
      MAE = MAE(pred, y_test),             # Mean Absolute Error
      BIAS = Metrics::bias(y_test, pred),  # Bias (Mean Error)
      Model = model                        # Model name correctly assigned
    )
  }) %>% 
    bind_rows()  # Combine all evaluation results into a single dataframe
}


# function to automatically generate vi plot
generate_vi_plot <- function(model, model_name, x_test, y_test, theme) {
  set.seed(42)  # For reproducibility
  vip(
    model,
    method = "permute",
    train = x_GWL_train,
    target = y_GWL_train,
    metric = "rsq",
    nsim = 1000,
	num_features = 20,
    pred_wrapper = predict,
    geom = "boxplot",
    mapping = aes_string(fill = "Variable"),
    aesthetics = list(color = "grey35")
  ) + 
    ggtitle(paste("GWL-", model_name, sep = "")) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) + theme
}


# function to automatically plot ICE and PDP
ICE_GWL <- function(model, pred.var, ...) {
  partial(
    model, 
    pred.var = pred.var, 
    ice = TRUE,          # Generate ICE plots
    center = TRUE,       # Center the ICE curves
    plot = TRUE,         # Plot the results
    rug = TRUE,          # Add rug plot
    alpha = 0.1,         # Transparency for ICE lines
    plot.engine = "ggplot2", 
    train = x_GWL_train, 
    type = "regression", 
    ...
  )
}


# Custom SHAP wrapper function for GWL models
shapviz_GWL <- function(model, ...) {
  shapviz(
    fastshap::explain(
      model,
      X = x_GWL_train,  # Training features
      nsim = 1000,      # Monte Carlo simulations
      pred_wrapper = predict,
      shap_only = FALSE,
      adjust = TRUE,
      ...
    )
  )
}


# Custom Beeswarm Plot Wrapper
bee_GWL <- function(model, color_bar_title = NULL, ...) {
  sv_importance(
    model, 
    kind = "beeswarm", 
    show_numbers = FALSE, 
    ...
  )
}


# Wrapper function for mean variable importance plot
vi_mean_GWL <- function(model, ...) {
  sv_importance(
    model, 
    kind = "bar", 
    show_numbers = TRUE, 
    ...
  )
}


# applying the GWL model to PHU -------------------------------------------------

#load list of static variables -------------

# Define the function to import SpatRaster objects. I used spatraster list not stack
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


# function to rename each raster inside the list to match the variables used in the GWL model
rename_rasters <- function(raster_list) {
  # Loop through each raster in the list
  for (i in seq_along(raster_list)) {
    
    # Get the raster layer
    raster_layer <- raster_list[[i]]
    
    # Extract the file path using sources() method
    file_path <- sources(raster_layer)
    
    # Get the file name without the extension
    raster_name <- sub("\\.[^.]*$", "", basename(file_path))
    
    # Remove "Sumsel_" or "_Sumsel" from the raster name
    raster_name <- sub("^Sumsel_", "", raster_name)   # Remove prefix "Sumsel_"
    raster_name <- sub("_Sumsel$", "", raster_name)   # Remove suffix "_Sumsel"
    
    # Set the raster's name to the modified name
    names(raster_layer) <- raster_name
    
    # Update the list with the modified raster layer
    raster_list[[i]] <- raster_layer
  }
  
  # Return the updated list
  return(raster_list)
}

# Function to convert a list of rasters to a raster stack
convert_to_stack <- function(raster_list) {
  # Check if all elements in the list are SpatRaster objects
  if (!all(sapply(raster_list, inherits, "SpatRaster"))) {
    stop("All elements in the list must be SpatRaster objects.")
  }
  
  # Convert the list to a raster stack
  raster_stack <- rast(raster_list)
  
  # Return the raster stack
  return(raster_stack)
}


# Function to convert a list of rasters to a raster stack and save it
convert_and_save_stack <- function(raster_list, output_path, file_name) {
  # Check if all elements in the list are SpatRaster objects
  if (!all(sapply(raster_list, inherits, "SpatRaster"))) {
    stop("All elements in the list must be SpatRaster objects.")
  }
  
  # Convert the list to a raster stack
  raster_stack <- rast(raster_list)
  
  # Construct the full file path
  output_file <- file.path(output_path, paste0(file_name, ".tif"))
  
  # Save the raster stack to the specified file path
  writeRaster(raster_stack, output_file, overwrite = TRUE)
  
  # Return a message indicating the file has been saved
  message("Raster stack saved as: ", output_file)
  
  # Return the raster stack (optional)
  return(raster_stack)
}


#----- INPUT GWL ----#
input_gwl_stack <- function(
  Precip_input_stack,   # Dynamic Precip stack (SpatRaster object or directory)
  ET_input_stack,       # Dynamic ET stack (SpatRaster object or directory)
  SM_input_stack,       # Dynamic SM stack (SpatRaster object or directory)
  computed_name,        # The raster name(s) to compute (e.g., 1:72 or "1":"72")
  static_input_stack,   # Static stack (SpatRaster object or directory)
  output_stack_name,    # Base name for the output stack files
  output_stack_directory # Directory to save the output stacks
) {
  # Function to load SpatRaster if input is a directory
  load_stack <- function(input) {
    if (is.character(input)) { # If input is a directory
      cat("Loading raster from:", input, "\n")
      return(rast(input)) # Load the raster stack
    } else if (inherits(input, "SpatRaster")) { # If input is already a SpatRaster
      return(input)
    } else {
      stop("Invalid input: must be either a SpatRaster object or a directory.")
    }
  }

  # Load all input stacks
  Precip_stack <- load_stack(Precip_input_stack)
  ET_stack <- load_stack(ET_input_stack)
  SM_stack <- load_stack(SM_input_stack)
  Static_stack <- load_stack(static_input_stack)

  # Ensure the output directory exists
  if (!dir.exists(output_stack_directory)) {
    dir.create(output_stack_directory, recursive = TRUE)
  }

  # Initialize a list to store results for each computed name
  result_list <- list()

  # Handle `computed_name` as numeric or character
  for (name in computed_name) {
    # Convert to character if `computed_name` is numeric
    name <- as.character(name)
    cat("Processing computed_name:", name, "\n")

    # Extract the raster by name
    precip_raster <- Precip_stack[[name]]
    et_raster <- ET_stack[[name]]
    sm_raster <- SM_stack[[name]]

    # Combine the extracted rasters into a single dynamic stack
    dynamic_stack <- c(precip_raster, et_raster, sm_raster)
    names(dynamic_stack) <- c("Precip", "ET", "SM")  # Rename the dynamic layers

    # Combine with the static stack
    final_stack <- c(dynamic_stack, Static_stack)

    # Define the output file name
    output_file <- file.path(output_stack_directory, paste0(output_stack_name, "_", name, ".tif"))

    # Save the combined stack
    cat("Saving stack to:", output_file, "\n")
    writeRaster(final_stack, output_file, overwrite = TRUE)

    # Store the final stack in the result list
    result_list[[name]] <- final_stack
  }

  cat("Stacks successfully saved to", output_stack_directory, "\n")

  # Return the list of SpatRaster stacks if multiple, or a single SpatRaster if only one
  if (length(result_list) == 1) {
    return(result_list[[1]])
  } else if (length(result_list) > 1) {
    return(result_list)
  } else {
    return(NULL)
  }
}


# function to add dummry years (2021:2023) -------------------------------------------------
add_dummy_year <- function(input_stack, dummyvars_number = 3, GWL_file_computed = 1:72, output_stack) {
  # Ensure GWL_file_computed is numeric
  GWL_file_computed <- as.numeric(GWL_file_computed)

  # Check if output directory exists, if not, create it
  if (!dir.exists(output_stack)) {
    dir.create(output_stack, recursive = TRUE)
  }

  # Generate the file names to process
  stack_files <- file.path(input_stack, paste0("input_GWL_", GWL_file_computed, ".tif"))

  # Check if all specified files exist
  missing_files <- stack_files[!file.exists(stack_files)]
  if (length(missing_files) > 0) {
    stop("The following files are missing: ", paste(missing_files, collapse = ", "))
  }

  # Initialize the progress bar
  pb <- progress_bar$new(total = length(stack_files), format = "  Processing [:bar] :percent eta: :eta")

  # Loop through each specified file and process the stack
  for (i in seq_along(stack_files)) {
    stack <- rast(stack_files[i])

    # Check if the stack has at least 25 layers
    if (nlyr(stack) < 25) {
      warning(paste("Stack", basename(stack_files[i]), "has less than 25 layers. Skipping."))
      next
    }

    # Get the reference layer for masking and extent
    reference_layer <- stack[[25]]
    last_layer <- stack[[nlyr(stack)]]

    # Create dummy layers with all values initially set to 1
    dummy_layers <- vector("list", dummyvars_number)
    for (j in seq_len(dummyvars_number)) {
      dummy_layer <- last_layer
      values(dummy_layer) <- 1  # Set all values to 1
      dummy_layer <- mask(dummy_layer, reference_layer)  # Apply mask
      dummy_layer <- crop(dummy_layer, reference_layer)  # Crop to reference extent
      dummy_layers[[j]] <- dummy_layer
    }

    # Generate names for the dummy layers
    start_year <- 2021
    names(dummy_layers) <- paste0("year_", start_year:(start_year + dummyvars_number - 1))

    # Assign 0 values based on the file index
    if (i <= 24) {
      values(dummy_layers[[2]]) <- 0  # year_2022 = 0
      values(dummy_layers[[3]]) <- 0  # year_2023 = 0
    } else if (i <= 48) {
      values(dummy_layers[[1]]) <- 0  # year_2021 = 0
      values(dummy_layers[[3]]) <- 0  # year_2023 = 0
    } else {
      values(dummy_layers[[1]]) <- 0  # year_2021 = 0
      values(dummy_layers[[2]]) <- 0  # year_2022 = 0
    }

    # Combine the original stack with the new dummy layers
    new_stack <- c(stack, rast(dummy_layers))

    # Save the modified raster stack to the specified output directory
    output_file <- file.path(output_stack, basename(stack_files[i]))
    writeRaster(new_stack, output_file, overwrite = TRUE)

    # Update the progress bar
    pb$tick()
  }
}


# Main function to predict the GWL raster and save the output
predict_GWL_raster <- function(model, input_dataset, input_number = NULL, 
                               output_directory = ".", output_name = "predicted_GWL.tif", 
                               batch_size = 1000, use_parallel = FALSE) {

  # Helper function to prepare data for prediction
  data_prep <- function(input_raster, model_features) {
    raster_features <- names(input_raster)
    
    # Check for missing features
    missing_features <- setdiff(model_features, raster_features)
    if (length(missing_features) > 0) {
      stop(paste("Missing features in input raster:", paste(missing_features, collapse = ", ")))
    }
    
    # Convert raster to data frame (ignoring NA cells) and subset to model features
    raster_matrix <- as.data.frame(input_raster, xy = TRUE, na.rm = TRUE)
    prediction_data <- raster_matrix[, model_features, drop = FALSE]
    
    list(raster_matrix = raster_matrix, prediction_data = prediction_data)
  }

  # Function to predict raster values in batches and return the predicted raster
  predict_GWL_raster_internal <- function(model, input_raster, model_features, batch_size = 1000) {
    # Prepare data for prediction
    data <- data_prep(input_raster, model_features)
    raster_matrix <- data$raster_matrix
    prediction_data <- data$prediction_data
    
    # Initialize progress bar
    n_batches <- ceiling(nrow(prediction_data) / batch_size)
    pb <- progress_bar$new(
      format = "  Predicting [:bar] :percent in :elapsed",
      total = n_batches,
      clear = FALSE,
      width = 60
    )
    
    # Make predictions in batches
    predictions <- numeric(nrow(prediction_data))  # Preallocate predictions vector
    for (i in seq_len(n_batches)) {
      start_idx <- ((i - 1) * batch_size) + 1
      end_idx <- min(i * batch_size, nrow(prediction_data))
      batch_data <- prediction_data[start_idx:end_idx, , drop = FALSE]
      
      # Predict for the current batch
      predictions[start_idx:end_idx] <- predict(model, newdata = batch_data)
      
      # Update progress bar
      pb$tick()
    }
    
    # Attach predictions back to the raster matrix and create a raster
    raster_matrix$prediction <- predictions
    predicted_raster <- terra::rast(raster_matrix[, c("x", "y", "prediction")], type = "xyz")
    
    # Set CRS to match the input raster
    crs(predicted_raster) <- crs(input_raster)
    
    return(predicted_raster)
  }

  # Extract model features
  model_features <- model$finalModel$feature_names

  # Case 1: If input_dataset is a single raster
  if (inherits(input_dataset, "SpatRaster")) {
    predicted_raster <- predict_GWL_raster_internal(model, input_dataset, model_features, batch_size)
    
    # Generate output file name dynamically
    if (grepl("input_GWL_", output_name)) {
      output_name <- paste0(output_name, "_", basename(input_dataset))
    }
    output_file <- file.path(output_directory, output_name)
    writeRaster(predicted_raster, output_file, overwrite = TRUE)
    cat("Prediction complete! Saved as:", output_file, "\n")

  # Case 2: If input_dataset is a list of rasters
  } else if (is.list(input_dataset)) {
    lapply(input_dataset, function(input_raster) {
      predicted_raster <- predict_GWL_raster_internal(model, input_raster, model_features, batch_size)
      
      # Generate output file name dynamically
      output_file <- file.path(output_directory, paste0(output_name, "_", basename(input_raster)))
      writeRaster(predicted_raster, output_file, overwrite = TRUE)
      cat("Prediction complete for", basename(input_raster), "! Saved as:", output_file, "\n")
    })

  # Case 3: If input_dataset is a directory, process each raster file
  } else if (is.character(input_dataset) && dir.exists(input_dataset)) {
    # List all input files in the directory matching the pattern 'input_GWL_*.tif'
    input_files <- list.files(input_dataset, pattern = "^input_GWL_.*\\.tif$", full.names = TRUE)
    
    # Check if input_number is provided and subset the files accordingly
    if (!is.null(input_number)) {
      # Ensure input_number is within the bounds of available files
      if (any(input_number > length(input_files))) {
        stop("input_number exceeds the number of available raster files.")
      }
      input_files <- input_files[input_number]
    }
    
    # Handle parallel processing (optional)
    if (use_parallel) {
      # Detect the number of cores and use parallel processing
      no_cores <- detectCores() - 1
      cl <- makeCluster(no_cores)
      clusterExport(cl, list("predict_GWL_raster_internal", "model", "model_features", "batch_size"))
      
      # Parallel processing of raster files
      parLapply(cl, input_files, function(input_file) {
        input_raster <- rast(input_file)
        predicted_raster <- predict_GWL_raster_internal(model, input_raster, model_features, batch_size)
        
        # Generate output file name dynamically
        output_file <- file.path(output_directory, paste0(output_name, "_", basename(input_file)))
        writeRaster(predicted_raster, output_file, overwrite = TRUE)
        cat("Prediction complete for", basename(input_file), "! Saved as:", output_file, "\n")
      })
      stopCluster(cl)
      
    } else {
      # Process each raster serially
      lapply(input_files, function(input_file) {
        input_raster <- rast(input_file)
        predicted_raster <- predict_GWL_raster_internal(model, input_raster, model_features, batch_size)
        
        # Generate output file name dynamically
        output_file <- file.path(output_directory, paste0(output_name, "_", basename(input_file)))
        writeRaster(predicted_raster, output_file, overwrite = TRUE)
        cat("Prediction complete for", basename(input_file), "! Saved as:", output_file, "\n")
      })
    }
    
  } else {
    stop("Invalid input dataset. Must be a single SpatRaster, list of rasters, or a directory path.")
  }
}









