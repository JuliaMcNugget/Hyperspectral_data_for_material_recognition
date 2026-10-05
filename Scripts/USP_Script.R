
#=======================================================================================================

# Libraries

library(RNetCDF) # netcdf files
library(raster)
library(ncdf4)
library(ggplot2)
library(sp) # for spatial data
library(sf)
library(dplyr) # for consolidations
library(caret)
library(reshape)
library(corrplot)
library(agricolae)
library(klaR)
library(terra)
#=======================================================================================================
setwd("C:/Users/Robot Mistress/Documents/R_Programs/current")
#=======================================================================================================

# Functions

#-------------------------------------------------------------------------------------------------------
# Merge all rfl files in a list to a matrix
# Might not actually need this
all_data_matrix <- function(rfl_list){
  
  data_full <- NULL
  
  for (rfl in rfl_list){
    print("Getting data from: ", rfl)
    data <- rast(rfl)
    
    print("Getting the values")
    vals <- values(data)
    
    print("Getting the coordinates...")
    coord <- crds(data)
    
    print("Merging coordinates and values")
    data_seperated <- cbind(coord, vals)
    
    print("Getting names...")
    
    # Get the wave length names
    waves <- as.numeric(gsub("X|\\.Nanometers", "", names(data)))
    
    colnames(data_seperated) <- c("x", "y", waves)
    
    if(is.null(data_full)){
      
      data_full <- data_seperated
      
    } else{
      print("Merging files...")
      data_full <- rbind(data_full, data_seperated, by = c("x", "y"))
      
    }
    
    
  }
  
  return(data_full)
  
}
#-------------------------------------------------------------------------------------------------------
process_dataset <- function(specdat, dada){
  fulldat <- merge(dada, specdat, by = c("x","y"))
  fulldat <- fulldat %>%
    dplyr::select(Code,
                  MapUnit,
                  x,
                  y,
                  everything(),
                  VM)
  fulldat <- fulldat %>%
    filter(Code != "Unknown")
  return(fulldat)
}

#-------------------------------------------------------------------------------------------------------

#Set the threshold
sample_dataset <- function(dat){  
  print("Reorganize Data...")
  fd <- dada %>% dplyr::select(-x, -y)
  
  print("Counting by classes...")
  unit_counts <- fd %>% count(Code, name = "Pixels")
  
  print("Removing NAs, Water and Unkown...")
  classCode <- as.factor(unit_counts$Code) # define classes
  keep <- !is.na(classCode) & classCode != "Unknown" & classCode != "water"
  unit_counts <- unit_counts[keep, ]
  classCode <- droplevels(classCode[keep])
  
  print("Establishing apropriate Threshold...")
  threshold <- min(unit_counts$Pixels, na.rm = FALSE)
  
  set.seed(1080)
  
  
  balanced <-
    fd %>%
    group_by(Code) %>%
    slice_sample(
      n = threshold,
      replace = FALSE
    ) %>%
    ungroup()
  
  return(balanced)
}
#=======================================================================================================
# Variables

#samples <- 5000 If Needed
nsims <- 3
k <- 10
#=======================================================================================================
# Files and Paths

# Input path for the EMIT tile
setwd("C:/Users/Robot Mistress/Documents/R_Programs/current")
rfl <- file.path("./R_input/Moab/f110609t01p00r13rdn_c_sc01_ort_img")

# Adjust coordinates data

print("Inputting coordinate data and adjusting")

#=======================================================================================================

#                                                       Step 1: Data Conversion
#                                                         & Shape Crop
print("Converting Data")
#=======================================================================================================

# Create empty variables
rfl_dat <- NULL
spec_all<- NULL

print("New Raster")
rfl_new <- rast(rfl)

if(is.rotated(rfl_new)){
  print("Data is rotated, rotating to adjust.")
  rfl_new <- rectify(rfl_new)
}

rast_max <- max(minmax(rfl_new))

#data_projected <- project(rfl_rectified, coord_vector)

print("printing image of each wave")
############ Verification Primarily ##########################
for(i in seq_along(rfl_new[1])){

  wave <- as.data.frame(rfl_new[[i]], xy = TRUE)
  
  file <- paste0("img_",colnames(wave)[3],".png")
  
  ggplot(wave, aes(x = x, y = y, fill = .data[[colnames(wave)[3]]])) + geom_point(aes(color = .data[[colnames(wave)[3]]]), size = 0.5) + scale_color_gradient(low = "darkorchid4", high = "yellow", limits = c(1, rast_max))
  ggsave(file)
}


#=======================================================================================================