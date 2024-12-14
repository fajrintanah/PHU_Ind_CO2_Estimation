// download Sentinel-1 soil moisture index ----------------------------------------

// Define the range of years and months to process
var years = [2021, 2022, 2023];
var months = ee.List.sequence(1, 12);

// Load your AOI (Area of Interest)
var aoi = ee.FeatureCollection('projects/ee-fajrintanah/assets/AOI_Rupat');

// Function to get the last day of the month dynamically
function getLastDayOfMonth(year, month) {
  return ee.Date.fromYMD(year, month, 1).advance(1, 'month').advance(-1, 'day');
}

// Function to get the start and end date for each bi-weekly period
function getBiweeklyPeriods(year, month) {
  var startOfMonth = ee.Date.fromYMD(year, month, 1);
  var midOfMonth = startOfMonth.advance(15, 'day');
  var endOfMonth = getLastDayOfMonth(year, month);

  return [
    { start: startOfMonth, end: midOfMonth.advance(-1, 'day') }, // First half of the month
    { start: midOfMonth, end: endOfMonth } // Second half of the month
  ];
}

// Function to convert month number to a zero-padded string (numeric)
function getMonthNumeric(month) {
  return month < 10 ? '0' + month : String(month); // Zero-pad months below 10
}

// Preprocessing function for Sentinel-1 data with dynamic min-max normalization using LiDAR data
function preprocessSentinel1(image) {
  // Apply speckle filtering
  var speckleFiltered = image.focal_mean(3, 'circle', 'pixels');
  
  // Radiometric Terrain Flattening: Remove terrain-induced distortions using LiDAR data
  var lidar = ee.Image('projects/ee-fajrintanah/assets/Lidar_Rupat'); // LiDAR data
  var terrainFlattened = speckleFiltered.subtract(lidar.multiply(0.0001)); // Adjust scaling as needed

  // Dynamically calculate min and max values over the AOI
  var minMax = terrainFlattened.reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry: aoi,
    scale: 10,
    maxPixels: 1e9
  });

  var minVal = ee.Number(minMax.get('VV_min'));
  var maxVal = ee.Number(minMax.get('VV_max'));

  // Normalize to the range 0-100
  var normalized = terrainFlattened
    .subtract(minVal)              // (value - min)
    .divide(maxVal.subtract(minVal)) // (max - min)
    .multiply(100);                 // Scale to 0-100

  return normalized.select(['VV']).clip(aoi); // Select only the VV band and clip to AOI
}

// Function to export bi-weekly soil moisture estimates using Sentinel-1 VV data
function exportAndVisualizeBiweeklySoilMoisture(year, month) {
  // Get the bi-weekly periods for the given year and month
  var periods = getBiweeklyPeriods(year, month);

  periods.forEach(function(period, index) {
    var startDate = period.start;
    var endDate = period.end;

    // Load Sentinel-1 data and apply preprocessing within the period
    var sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD')
      .filterBounds(aoi)
      .filterDate(startDate, endDate)
      .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
      .filter(ee.Filter.eq('instrumentMode', 'IW'))
      .map(preprocessSentinel1);  // Apply preprocessing with normalization

    // Mosaicking step to handle overlap between tiles (combine images)
    var mosaicImage = sentinel1.mosaic();  // Apply mosaic to combine the collection

    // Check image statistics for the exported region
    var stats = mosaicImage.reduceRegion({
      reducer: ee.Reducer.minMax(),
      geometry: aoi.geometry(),
      scale: 10,
      maxPixels: 1e9
    });
    print('Image stats:', stats);

    // Add the mosaicked layer to the map
    var visParams = {
      min: 0,     // Normalized minimum value
      max: 100,   // Normalized maximum value
      palette: ['blue', 'green', 'yellow', 'red'] // Visualization color palette
    };

    Map.addLayer(mosaicImage.clip(aoi), visParams, year + '' + getMonthNumeric(month) + '' + (index + 1)); // Numeric month

    // Define export parameters and export the image to Google Drive
    var exportRegion = aoi.geometry().bounds();
    print('Export Region:', exportRegion);

    Export.image.toDrive({
      image: mosaicImage.clip(aoi),  // Clip the image to the AOI for export
      description: year + '_' + getMonthNumeric(month)  + '_' + (index + 1),
      scale: 10,
      folder: 'SM_COP_RUPAT',
      region: exportRegion,
      crs: 'EPSG:4326',
      maxPixels: 1e13
    });
  });
}

// Center the map on the AOI
Map.centerObject(aoi, 8); // Adjust zoom level as needed

// Loop through each year and month and call the export function
years.forEach(function(year) {
  months.getInfo().forEach(function(month) {
    exportAndVisualizeBiweeklySoilMoisture(year, month);
  });
});
