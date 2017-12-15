#--------------------------------------------------
# R Script for creating mortality rates for California based on annual U.S. Forest Service Aerial Detection Summaries
# Ben Sleeter, U.S. Geological Survey
# For more information contact: bsleeter@usgs.gov; (253) 343-3363
# Created December 4, 2017
#--------------------------------------------------

#--------------------------------------------------
# Process each year of ADS data
# Steps are:
# 1. Import shapefile of mortality areas.
# 2. Convert to raster using the "TPA1" field (dead trees per acre), convert to extent and reference system of LUCAS State Class Map.
# 3. Reclassify ADS raster into 3 severity classes based on 0-10 TPA (low), 10-20 (medium), and 20+ (high).
# 4. Summarize state wide area of mortality and divide by forest area to derive mortality rate.
# 5. Merge into a single data frame and plot results.

#--------------------------------------------------
# Load packages
library(raster)
library(maptools)
library(rgdal)
library(tidyverse)
library(rasterVis)
#--------------------------------------------------



#--------------------------------------------------
# Read in a raster of land cover
lulc = raster("Rasters/stateClass.tif")
levelplot(lulc)
forest = reclassify(lulc, c(-Inf, 5.4, NA, 5.5, 6.5, 1, 6.6, Inf, NA))
forestArea = cellStats(forest, sum)
levelplot(forest)

# Define coordinate reference system and extent
crs = crs(forest)
extent = extent(forest)



#--------------------------------------------------
# Read in a raster of ecoregions
ecoregions = raster("Rasters/ecoregions.tif")
ecoregions[ecoregions > 85] = NA
levelplot(ecoregions)


# Calculate forest area by ecoregion
zonalForest = data.frame(zonal(forest, ecoregions, sum)) %>%
    rename(forest = value) %>%
    mutate(km2 = forest * 0.0729)
head(zonalForest)





#---------------------------
# Change the year in the following line to run for additonal years
# Script write raster output and a data frame summary to disk

year = 2002
# Read in the vector data
v = readOGR(paste("Shapefiles/ads", year, "_mortality.shp", sep = ""))
v = spTransform(v, CRS('+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'))

# Convert to raster
r = raster(extent(lulc), res = res(lulc), crs = crs(lulc))
r = rasterize(v, r, "TPA1", update = TRUE, updateValue = "NA")

# Create low, medium, and high severity class rasters
low = reclassify(r, c(-Inf, 10, 1, 10, Inf, NA))
med = reclassify(r, c(-Inf, 10, NA, 10, 20, 1, 20, Inf, NA))
high = reclassify(r, c(-Inf, 20.0, NA, 20.1, Inf, 1))

levelplot(high, maxpixels = 2e5)

# Calculate zonal statistics data frame
zonalLow = data.frame(zonal(low, ecoregions, sum)) %>%
    mutate(timestep = year, severity = "low") %>%
    left_join(zonalForest, by = "zone")

zonalMed = data.frame(zonal(med, ecoregions, sum)) %>%
    mutate(timestep = year, severity = "medium") %>%
    left_join(zonalForest, by = "zone")

zonalHigh = data.frame(zonal(high, ecoregions, sum)) %>%
    mutate(timestep = year, severity = "high") %>%
    left_join(zonalForest, by = "zone")

zonal = rbind(zonalLow, zonalMed, zonalHigh) %>%
    mutate(mortality = (value*0.0729) / km2)
head(zonal)

write.csv(zonal, paste("R Outputs/zonal", year, ".csv", sep = ""), row.names = F)
writeRaster(low, paste("R Outputs/low", year, "tif", sep = "."), format = "GTiff", overwrite = TRUE)
writeRaster(med, paste("R Outputs/med", year, "tif", sep = "."), format = "GTiff", overwrite = TRUE)
writeRaster(high, paste("R Outputs/high", year, "tif", sep = "."), format = "GTiff", overwrite = TRUE)



save.image()




#---------------------------
# Read in an ecoregion id/name lookup
ecoLookup = read.csv("ecoregion_lookup.csv", header = T)

# Merge the dataframes on disk into a single dataframe
list = list.files(path = "R Outputs", pattern = "*.csv")
ecoSeverity = lapply(paste("R Outputs/", list, sep = ""), read_csv) %>%
    bind_rows() %>%
    left_join(ecoLookup)
head(ecoSeverity)


# Plot the results by ecoregion
ggplot(data = ecoSeverity, aes(x = timestep, y = mortality, fill = severity)) +
    geom_bar(stat = "identity") +
    facet_wrap(~zone) +
    theme_bw()

# Write the combined dataframe to disk
write.csv(ecoSeverity, "R Outputs/ecoregion-severity.csv", row.names = FALSE)

distributions = data.frame(StratumID = ecoSeverity$ecoregion,
                           SecondaryStratumID = "",
                           DistributionTypeID = "",
                           E)
