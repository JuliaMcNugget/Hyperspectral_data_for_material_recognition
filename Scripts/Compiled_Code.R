

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

#=======================================================================================================

#=======================================================================================================

# Functions

#-------------------------------------------------------------------------------------------------------

process_dataset <- function(specdat, dada){
  fulldat <- merge(dada, specdat, by = c("lon","lat"))
  fulldat <- fulldat %>%
    dplyr::select(Code,
                  MapUnit,
                  lon,
                  lat,
                  everything(),
                  VM)
  fulldat <- fulldat %>%
    filter(Code != "Unknown")
  return(fulldat)
}

#-------------------------------------------------------------------------------------------------------
#Balance out the datasets

#Set the threshold
sample_dataset <- function(dat, dataset_name){
  
  fd <- dat %>%
    dplyr::select(-lon, -lat)
  unit_counts <- fd %>%
    count(Code, MapUnit, name = "Pixels")
  keep_classes <-
    unit_counts %>%
    filter(Pixels >= threshold)
  removed_classes <-
    unit_counts %>%
    filter(Pixels < threshold)
  fd <- fd %>%
    filter(Code %in% keep_classes$Code)
  final_counts <- fd %>%
    count(Code, MapUnit, name = "Pixels")
  
  write.csv(
    final_counts,
    file.path(
      "R_output",
      paste0("TrainingSummary_", dataset_name, ".csv")
    ),
    row.names = FALSE
  )
  
  min_pixels <- min(final_counts$Pixels)
  
  limiting_class <-
    final_counts %>%
    filter(Pixels == min_pixels)
  set.seed(1080)
  balanced <-
    fd %>%
    group_by(Code) %>%
    slice_sample(
      n = min_pixels,
      replace = FALSE
    ) %>%
    ungroup()
  
  return(balanced)
}
#=======================================================================================================


#=======================================================================================================
# Files and Paths

# Input path for the EMIT tile
fp        <- file.path("R_input/EMIT_L2A_RFL_001_20250531T200206_2515113_018_ortho_subset_aid0000.nc")
#=======================================================================================================

#   Step 1: Data Conversion

#=======================================================================================================
# Preparing data 
dat <- ncdf4::nc_open(fp) # open datacube 

# Extract the reflectance values from cube
refl <- ncdf4::ncvar_get(dat, "reflectance")

# Extract essential positioning and wavelength data
wvl <- dat$dim$wavelengths$vals # extract wavelengths
lat <- dat$dim$lat$vals # extract latitude
lon <- dat$dim$lon$vals # extract longitude


# Alternative to ncdf4
EMIT <- brick(fp, subds = "reflectance")
plot(EMIT[[1]])

#--------------------------------
# Uses extent of EMIT to plot
mmini <- extent(EMIT)
plot(mmini)

EMIT_crop <- raster::crop(EMIT, mmini)

# Confirm wavelength metadata matches raster layers
stopifnot(nlayers(EMIT_crop) == length(wvl))
#-------------------------------------------------------------------------------
# this is probably the place to use bad_bands 
# removal once we get everything working





#-------------------------------------------------------------------------------

rfl_values <- raster::getValues(EMIT_crop)

# Raster cell numbers
cells <- seq_len(ncell(EMIT_crop))

# Coordinates of each pixel
xy <- raster::xyFromCell(EMIT_crop, cells)

# Row and column indices in the cropped raster
dtrk <- raster::rowFromCell(EMIT_crop, cells)
ctrk <- raster::colFromCell(EMIT_crop, cells)

# Combine coordinates and spectra
spec <- cbind(
  dtrk = dtrk,
  ctrk = ctrk,
  lat  = xy[, 2],
  lon  = xy[, 1],
  rfl_values
)

# Apply wavelength names to reflectance columns
colnames(spec) <- c(
  "dtrk",
  "ctrk",
  "lat",
  "lon",
  as.character(wvl)
)

spec <- data.frame(spec, check.names = FALSE)
# Save file #
points <- subset(spec, select = c("lon","lat")) # extract the coordinates only
points <- data.frame(points)

  # Cleanup! 
nc_close(dat)


rm(
   EMIT,
   mmini,
   dat,
   EMIT_crop,
   cells,
   xy,
   dtrk,
   ctrk,
   rfl_values,
   wvl,
   fp
)
#=======================================================================================================

#   Step 2: NDVI
print("NVDI")
#=======================================================================================================

  # Identify bands to use for NDVI
red_band <- which.min(abs(as.numeric(colnames(spec)[-(1:4)]) - 660))
nir_band <- which.min(abs(as.numeric(colnames(spec)[-(1:4)]) - 860))

red <- spec[, 4 + red_band]
nir <- spec[, 4 + nir_band]
ndvi <- (nir - red) / (nir + red)

  # Mask non-veg
spec_nonveg <- spec[ndvi <= 0.35, ]
spec_veg    <- spec[ndvi > 0.35, ]

rm(red_band,
   nir_band,
   red,
   nir,
   ndvi
)
#=======================================================================================================

#   Step 3: Extractor
print("Extracting Data")
#=======================================================================================================
  # Read shapefile
shp_fp    <- file.path("R_input/Shapes/Final_map.shp") #*load the shapefile for the Henry map **yellowstone shape**

units <- raster::shapefile(shp_fp) #i think there is an issue with the crs of this object for tomorrow

units <- spTransform(units, CRSobj =CRS("EPSG:32612"))

  # verify shapefile visually
plot(units) 

  #Create the buffer
buffer_dist <- -30
units_sf <- st_as_sf(units)
units_sf <- st_buffer(units_sf, dist = buffer_dist)
units_sf <- units_sf[!sf::st_is_empty(units_sf), ]
units <- as(units_sf, "Spatial")

  # Labeling
uni <- units["MapUnit"]
names(uni) <- "Code"
uni$Code <- as.factor(uni$Code)

poi_sp <- SpatialPointsDataFrame(coords = cbind(points$lon,points$lat), # GPS points are in WGS84
                                 data=points,proj4string=CRS("EPSG:4326"))

  # Transform coordinates of points to match geologic map
poi_UTM <- spTransform(poi_sp,CRSobj =CRS("EPSG:32612")) 

  # check if there's overlap, the red is showing where the emit pixel is
plot(uni)
plot(poi_UTM, cex=2, col="red", pch=16, add = T)
dev.off()

  # Rasterize shapefile
exi <- extent(uni) # bounds of shapefile

empty_ras <- raster(exi, resolution=60, crs="EPSG:32612")# create empty raster based on extent
uni_ras <- raster::rasterize(uni,empty_ras, field="Code")  

  # Create lookup table for code names: check if ok
look <- as.data.frame(matrix(levels(uni$Code), nrow=length(levels(uni$Code)))) # code names from shapefile
names(look)<- "Code"
look$ID <- seq(1:nrow(look)) # assign id to each code

plot(uni_ras)
plot(poi_UTM, add=T)

  # Extract data
poi_classes <- raster::extract(uni_ras, poi_UTM, sp = TRUE)
dat <- as.data.frame(poi_classes)

  # Data wrangling

names(dat)[3:5] <- c("ID", "X_EPSG_32612", "Y_EPSG_32612")
dada <- merge(dat,look, by="ID")

  # Count number of EMIT pixels per unit and save it
pixel_counts <- dada%>%
  group_by(Code) %>%
  summarise(Pixels = n())

print(pixel_counts)

  # Clean up 
rm(dat,empty_ras,exi,look,poi_classes,poi_sp,poi_UTM,points,uni,uni_ras,units, pixel_counts, units_sf, buffer_dist)
dada<-dada[,-1]
dada<-dada[,c(5,2,1,3,4)]

#=======================================================================================================

#   Step 4: Consolidation of Units
print("Consolidating Units and Data")
#=======================================================================================================

output_fp <- file.path("R_output", "UnitCoord_FB_CHM.csv")

# Read data

# Geologic table look-up

unit_lookup <- data.frame(
  Code = c(
    "Qa",
    "Qpce",
    "Qpcw",
    "Qpuc",
    "Qpus",
    "Qyl2",
    "Qyl3",
    "Qyl4",
    "Qyl5",
    "water",
    "Unknown"
  ),
  
  MapUnit = c(
    "Alluvium",
    "Placeholder",
    "Placeholder",
    "Placeholder",
    "Placeholder",
    "Placeholder",
    "Placeholder",
    "Placeholder",
    "Placeholder",
    "Water",
    "Unknown"
  ),
  
  stringsAsFactors = FALSE
  
)

#Join the names to the table
dada <- left_join(dada, unit_lookup, by = "Code")

# Check labels
print(table(dada$Code))
print(unique(dada[, c("Code", "MapUnit")]))

print(names(dada))

print(unique(dada[, c("Code", "MapUnit")]))

rm(
   unit_lookup,
   output_fp
)
# Save #
#write.csv(dada, file = output_fp, row.names = FALSE)

##################### Works up until here ###################################################
#-------------------------------------------------------------------------------------------------------
#=======================================================================================================

# # Step 5: Merger Equal and Threshold
print("Merging Data and establishing thresholds")
# #=======================================================================================================
#               ###Script 5, Merger weighs all training units equally

spec_nonveg$VM <- "nonveg"
spec_veg$VM    <- "veg"

spec_all <- bind_rows(spec_nonveg, spec_veg)

# Merge data
threshold <- 50
rownames(dada)<-NULL
#Create datasets
fulldat_nonveg <- process_dataset(spec_nonveg, dada)
fulldat_veg    <- process_dataset(spec_veg, dada)
fulldat_all    <- process_dataset(spec_all, dada)


#Create new balanced datasets
balanced_nonveg <- sample_dataset(fulldat_nonveg, "nonveg_thresh")
balanced_veg <- sample_dataset(fulldat_veg, "veg_thresh")
balanced_all <- sample_dataset(fulldat_all, "all_thresh")

#Save new datasets

#=======================================================================================================

# Step 6: PLSDA

#=======================================================================================================
#     Formula
#-------------------------------------------------------------------------------------------------------
print("PLSDA")
PLSDA_1 <- function(balanced, training_summary, pdf_fp, nsims = 3, mods = list(), k=15){
  ## Preparation ##
  nzv <- nearZeroVar(balanced)
  F_balanced <- balanced[, -nzv]
  classi <- as.factor(balanced$Code) # define classes
  is.wave <- suppressWarnings(!is.na(as.numeric(colnames(balanced))))
  wvl  <- colnames(balanced)[is.wave]
  spec <- balanced[, wvl]
  spec <- data.frame(spec)
  
  set.seed(1840) ### create random number for obs in each Code/class 
  rndid <- with(balanced, ave(1:nrow(balanced), Code, FUN=function(x) {sample.int(length(x))}))
  
  ## Partition  ##
  set.seed(1080) # for consistency
  inTrain <- createDataPartition(y =classi, p = .7, list = FALSE) # e.g. 70% spectra per class
  training <- spec[inTrain,seq(1,ncol(spec), by=2)] # resampling optional
  testing <- spec[-inTrain,seq(1,ncol(spec), by=2)] # pick all that was not picked in inTrain
  trainclass <- classi[inTrain] # define units for training
  testclass <- classi[-inTrain] # define units for testing
  

  
  table(testclass) # verify
  table(trainclass)
  table(F_balanced$Code)
  ## Finding components  ## 
  
  for (nsim in seq(nsims)){
    inTrain <- createDataPartition(y =classi, p = 0.7, list = FALSE) ## or create partitioning using fractions, e.g. 70:30
    print(nsim)
    flush.console()
    set.seed(nsim)
    traini <- spec[inTrain,seq(1,length(wvl), by=2)] 
    testi <- spec[-inTrain,seq(1,length(wvl), by=2)]
    trainclass <- classi[inTrain];testclass <- classi[!(inTrain)]
    plsFit <- train(traini, trainclass, method = "simpls", tuneLength = k,
                    probMethod="Bayes", trControl = trainControl(method="cv")) 
    mods[[nsim]] <- plsFit
    # , preProcess = "nzv"
  }
  
  #head(mods[[2]]$results) # verify
  
  ncomps <- vector(length = nsims) # create vector
  for (i in 1:nsims){ 
    ncomps[i]<-mods[[i]]$finalModel$ncomp # determining number of components in each model
  }

  ## Kappa statistic ##
  
  kappas <- data.frame(ncomps= 1:k,matrix(NA, nrow = k, ncol = length(mods))) # creating kappa matrix
  for (i in 1:length(mods)){
    kappas[,i+1] <- mods[[i]]$results$Kappa # defines the coefficenets
  }
  
  ## Tukey test ##
  
  kapp <- as.data.frame(as.numeric(t(kappas[,-1]))) # combines all kappas
  kapp <- cbind(kapp, rep(1:k, each=length(mods))) # labels kappas per model
  names(kapp) <- c("Kappa", "ncomps") # nicer labeling
  
  kapp$ncomps <- as.factor(kapp$ncomps) # change to factor
  
  modi <- lm(Kappa~ncomps, kapp) # combine and reformat data
  tuk <- HSD.test(modi,"ncomps") # perform tukey test
  
  tuk_dat <- as.data.frame(tuk$groups) # convert tukey test into dataframe
  tuk_dat$var <- as.numeric(row.names(tuk_dat)) # add the component column
  tuk_dat <- tuk_dat[order(tuk_dat$var,decreasing = F),] # reorder from small to large
  
  letters <- as.character(tuk_dat$groups) # constrain order of latter for graph
  
  ## Kappa Tukey plot ## 
  
  pdf(pdf_fp, width = 5, height = 4)
  par(bty = "l")
  boxplot(kapp$Kappa ~ kapp$ncomps, ylim = c(0,1),
          xlab = "Number of components", ylab = "Kappa")
  text(x = 1:k, y = rep(1, k), letters)
  
  return(tuk_dat)
  
}

#----------# All #-----------#

# PDF Output
pdf_fp    <- file.path("R_output/All", "PLSDA_kappas_all_thresh.pdf")

#Testing pixel count
print(table(balanced_all$Code))
print(nrow(balanced_all))

training_summary <- as.data.frame(table(balanced_all$Code))
names(training_summary) <- c("Code", "TrainingPixels")
training_summary$Dataset <- "RandUnitCoord_all_thresh.csv"

tuk_dat_all <- PLSDA_1(balanced_all, training_summary, pdf_fp)

dev.off()

#----------# Veg #-----------#
# PDF Output
pdf_fp    <- file.path("R_output/Veg/", "PLSDA_kappas_veg_thresh.pdf")

#Testing pixel count
print(table(balanced_veg$Code))
print(nrow(balanced_veg))

training_summary <- as.data.frame(table(balanced_veg$Code))
names(training_summary) <- c("Code", "TrainingPixels")
training_summary$Dataset <- "RandUnitCoord_veg_thresh.csv"

tuk_dat_veg <- PLSDA_1(balanced_veg, training_summary, pdf_fp)

dev.off()

#----------# Non Veg #-----------#
# PDF Output
pdf_fp    <- file.path("R_output/Nonveg", "PLSDA_kappas_nonveg_thresh.pdf")

#Testing pixel count
print(table(balanced_nonveg$Code))
print(nrow(balanced_nonveg))

training_summary <- as.data.frame(table(balanced_nonveg$Code))
names(training_summary) <- c("Code", "TrainingPixels")
training_summary$Dataset <- "RandUnitCoord_nonveg_thresh.csv"
tuk_dat_nonveg <- PLSDA_1(balanced_nonveg, training_summary, pdf_fp)

dev.off()

#=======================================================================================================
# Step 7: Extractor
#=======================================================================================================
#     Function
#-------------------------------------------------------------------------------------------------------

print("Extracting Data and Finalizing PLSDA")
PLSDA_2 <- function(balanced, out_finmods, out_scores, out_probmean, out_probplot, 
                    out_confumean,out_confuperc, corrplot, out_loadings){
  
  ## Preperation ##
  nzv <- nearZeroVar(balanced)
  F_balanced <- balanced[, -nzv]
  nsims <- 10 # iterations, as it increase time increases but better average
  classi <- as.factor(balanced$Code)# define classes
  is.wave <- suppressWarnings(!is.na(as.numeric(colnames(balanced))))
  wvl <- colnames(balanced)[is.wave] # define wvl range for spectral matrix, check your column names
  spec <- balanced[, wvl]
  nwaves <- 1
  compi <- 8 ### select number of components based on Kappa plot / Tukey test from pt 1
  finmods <- list() #
  
  for (nsim in 1:nsims){
    print(nsim)
    flush.console()
    set.seed(nsim)
    #inTrain <- rndid[[nsim]]<= ??? # number train per code
    inTrain <- createDataPartition(y =classi, p = 0.7, list = FALSE)
    training <- spec[inTrain,seq(1,length(wvl), by=nwaves)]
    testing <- spec[-inTrain,seq(1,length(wvl), by=nwaves)]
    trainclass <- as.factor(classi[inTrain]); testclass <- as.factor(classi[!(inTrain)])
    
    finalModel <- plsda(training,trainclass, ncomp=compi, probMethod = "Bayes", method = "simpls")
    finmods[[nsim]] <- finalModel
  }
  
  saveRDS(finmods, out_finmods)
  
  
  # Scores
  scores_list <- lapply(seq_along(finmods), function(i) {
    sc <- as.data.frame(finmods[[i]]$scores)
    sc$Code <- as.character(trainclass)
    sc$Sim <- i
    sc
  })
  
  scores_df <- do.call(rbind, scores_list)
  write.csv(scores_df, out_scores, row.names = FALSE)
  
  ## Probabilities and confusion matrix ##
  probis <- list()
  confus <- list()
  
  for (nsim in seq(nsims)){
    print(nsim)
    flush.console()
    set.seed(nsim)
    #inTrain <- rndid[[nsim]]<= ???
    inTrain <- createDataPartition(y = classi, p = 0.7, list = FALSE)
    testing <- spec[-inTrain,seq(1,length(wvl), by=nwaves)]
    testclass <- as.factor(classi[-inTrain])
    
    plsProbs <- predict(finmods[[nsim]], newdata = testing, type = "prob")
    plsClasses <- predict(finmods[[nsim]], newdata = testing)
    confus[[nsim]] <- confusionMatrix(data = plsClasses, testclass)
    
    probs <- as.data.frame(plsProbs)
    names(probs) <- sapply(strsplit(names(probs),split = ".n"),"[",1)
    probs <- cbind(testclass, probs)
    probis[[nsim]] <- probs 
  }
  
  ## Probability mean ##
  
  arr <- array(unlist(probis), dim = c(dim(probis[[1]]), nsims))
  prob_mean <- apply(arr, 1:2, mean)
  prob_mean <- as.data.frame(prob_mean)
  colnames(prob_mean) <- colnames(probis[[1]]) #was 98
  prob_mean$testclass <- probis[[1]]$testclass #was 97
  write.csv(prob_mean, out_probmean, row.names = FALSE)
  
  ## Probability Plot ##
  pp <- melt(prob_mean, id = "testclass")
  pp$position <- ifelse(pp$testclass == pp$variable, 2, 1)
  pp$testclass <- factor(pp$testclass, levels = rev(levels(pp$testclass)))
  
  pdf(out_probplot, width = 6, height = 4)
  ggplot(pp, aes(x = testclass, y = value, fill = variable, group = position)) +
    geom_bar(stat = "identity", position = "fill") +
    theme(legend.title = element_blank()) +
    labs(x = "Code", y = "") +
    coord_flip()
  dev.off()
  
  ## Confusion matrix! ##
  tabs <- lapply(confus, function(x) x$table)
  tabsi <- Reduce('+', tabs)
  tab_mean <- as.data.frame.matrix(tabsi / length(confus))
  write.csv(tab_mean, out_confumean, row.names = TRUE)
  
  sums <- colSums(tab_mean)
  tabs_perc <- sweep(tab_mean, 2, sums, FUN = "/")
  write.csv(tabs_perc, out_confuperc, row.names = TRUE)
  
  # Confu plot
  pdf(out_corrplot, width = 7, height = 6, pointsize = 9)
  corrplot(as.matrix(tabs_perc), addCoef.col = 1, tl.srt = 70,
           col = colorRampPalette(c("black","brown","gold","forestgreen"))(20),
           tl.col = 1, tl.offset = 1.5, cl.ratio = 0.2,
           cl.align.text = "l", cl.cex = 0.9, mar = c(1,3,3,3))
  mtext("Prediction", 2, line = 2.5, cex = 1.2)
  mtext("Reference", at = 2, line = 1, cex = 1.2)
  dev.off()
  
  ## Importance of bands ##
  lls <- lapply(finmods, function(m) abs(loadings(m)[, 1:compi]))
  sumis <- lapply(lls, rowSums)
  mm <- apply(simplify2array(sumis), 1, mean)
  ss <- apply(simplify2array(sumis), 1, sd)
  mm_df <- data.frame(mean = mm, sd = ss)
  
  pdf(out_loadings, width = 6, height = 4)
  plot(mm_df$mean, type = "l", bty = "l", ylab = "abs(loadings)", xlab = "Bands")
  polygon(x = c(1:length(mm), length(mm):1),
          y = c(mm_df$mean - mm_df$sd, rev(mm_df$mean + mm_df$sd)),
          col = "grey", border = "grey")
  lines(mm_df$mean, type = "l")
  dev.off()
  
  ## Model statistics ##
  
  accu   <- sapply(confus, function(x) sum(diag(x$table)) / sum(x$table))
  kappas <- sapply(confus, function(x) x$overall["Kappa"])
  sensi  <- sapply(confus, function(x) mean(x$byClass[, "Sensitivity"], na.rm = TRUE))
  specifi<- sapply(confus, function(x) mean(x$byClass[, "Specificity"], na.rm = TRUE))
  preci  <- sapply(confus, function(x) mean(x$byClass[, "Precision"], na.rm = TRUE))
  
  ssum <- data.frame(
    Stat = c("Accuracy","Kappa","Sensitivity","Specificity","Precision"),
    Mean = c(mean(accu), mean(kappas), mean(sensi), mean(specifi), mean(preci)),
    SD   = c(sd(accu), sd(kappas), sd(sensi), sd(specifi), sd(preci))
  )
}
##----------- ALL -------------##
# Input/Output paths

out_finmods   <- file.path("R_output/All/", "PLSDA_finmods_all_thresh.rds")
out_scores    <- file.path("R_output/All/", "PLSDA_Scores_all_thresh.csv")
out_probmean  <- file.path("R_output/All/", "PLSDA_probmean_all_thresh.csv")
out_probplot  <- file.path("R_output/All/", "probability_plot_all_thresh.pdf")
out_confumean <- file.path("R_output/All/", "PLSDA_confumean_all_thresh.csv")
out_confuperc <- file.path("R_output/All/", "PLSDA_confuperc_all_thresh.csv")
out_corrplot  <- file.path("R_output/All/", "PLSDA_corrplot_all_thresh.pdf")
out_loadings  <- file.path("R_output/All/", "PLSDA_loadings_all_thresh.pdf")
out_summary   <- file.path("R_output/All/", "PLSDA_summary_all_thresh.csv")

ssum_all <-  PLSDA_2(balanced_all, out_finmods, out_scores, out_probmean, out_probplot, 
                     out_confumean,out_confuperc, corrplot, out_loadings)

write.csv(ssum_all, out_summary, row.names = FALSE)

##----------- VEG -------------##
# Input/Output paths

out_finmods   <- file.path("R_output/Veg/", "PLSDA_finmods_veg_thresh.rds")
out_scores    <- file.path("R_output/Veg/", "PLSDA_Scores_veg_thresh.csv")
out_probmean  <- file.path("R_output/Veg/", "PLSDA_probmean_veg_thresh.csv")
out_probplot  <- file.path("R_output/Veg/", "probability_plot_veg_thresh.pdf")
out_confumean <- file.path("R_output/Veg/", "PLSDA_confumean_veg_thresh.csv")
out_confuperc <- file.path("R_output/Veg/", "PLSDA_confuperc_veg_thresh.csv")
out_corrplot  <- file.path("R_output/Veg/", "PLSDA_corrplot_veg_thresh.pdf")
out_loadings  <- file.path("R_output/Veg/", "PLSDA_loadings_veg_thresh.pdf")
out_summary   <- file.path("R_output/Veg/", "PLSDA_summary_veg_thresh.csv")

ssum_veg <-  PLSDA_2(balanced_all, out_finmods, out_scores, out_probmean, out_probplot, 
                     out_confumean,out_confuperc, corrplot, out_loadings)

write.csv(ssum_veg, out_summary, row.names = FALSE)

##----------- NON VEG -------------##
# Input/Output paths

out_finmods   <- file.path("R_output/Nonveg/", "PLSDA_finmods_nonveg_thresh.rds")
out_scores    <- file.path("R_output/Nonveg/", "PLSDA_Scores_nonveg_thresh.csv")
out_probmean  <- file.path("R_output/Nonveg/", "PLSDA_probmean_nonveg_thresh.csv")
out_probplot  <- file.path("R_output/Nonveg/", "probability_plot_nonveg_thresh.pdf")
out_confumean <- file.path("R_output/Nonveg/", "PLSDA_confumean_nonveg_thresh.csv")
out_confuperc <- file.path("R_output/Nonveg/", "PLSDA_confuperc_nonveg_thresh.csv")
out_corrplot  <- file.path("R_output/Nonveg/", "PLSDA_corrplot_nonveg_thresh.pdf")
out_loadings  <- file.path("R_output/Nonveg/", "PLSDA_loadings_nonveg_thresh.pdf")
out_summary   <- file.path("R_output/Nonveg/", "PLSDA_summary_nonveg_thresh.csv")

ssum_nonveg <-  PLSDA_2(balanced_all, out_finmods, out_scores, out_probmean, out_probplot, 
                     out_confumean,out_confuperc, corrplot, out_loadings)

write.csv(ssum_nonveg, out_summary, row.names = FALSE)


#=======================================================================================================
# Step 8: Extractor
# #=======================================================================================================
###Script 8: Compare the models: This is a complete rewrite, so good luck to us all

#Input files

#input_fp <- file.path("R_output", "UnitCoordSpec_all.csv")
model_all_fp <- file.path("R_output/All/", "PLSDA_finmods_all_thresh.rds")
model_veg_fp <- file.path("R_output/Veg/", "PLSDA_finmods_veg_thresh.rds")
model_nonveg_fp <- file.path("R_output/Nonveg/", "PLSDA_finmods_nonveg_thresh.rds")

#Set output files
pred_all_fp <- file.path("R_output/Final_Outputs/", "PredUnitCoord_all.csv")
pred_veg_fp <- file.path("R_output/Final_Outputs/", "PredUnitCoord_veg.csv")
pred_nonveg_fp <- file.path("R_output/Final_Outputs/", "PredUnitCoord_nonveg.csv")

map_all_fp <- file.path("R_output/Final_Outputs/", "PredictionMap_all.png")
map_veg_fp <- file.path("R_output/Final_Outputs/", "PredictionMap_veg.png")
map_nonveg_fp <- file.path("R_output/Final_Outputs/", "PredictionMap_nonveg.png")

#Read data
dati <- fulldat_all
coords <- dati[, c("Code", "X_EPSG_32612", "Y_EPSG_32612")]
is.wave <- suppressWarnings(!is.na(as.numeric(colnames(dati))))
spec <- dati[, is.wave]

#Load models
mods_all <- readRDS(model_all_fp)
mods_veg <- readRDS(model_veg_fp)
mods_nonveg <- readRDS(model_nonveg_fp)


##Just bug fixing
ncol(spec)
ncol(mods_all[[1]]$scores)
dim(mods_all[[1]]$projection)
class(mods_all[[1]])

print("Making Predictions:")

#Predicting Using Majority Vote (Admitidely, maybe not the most efficient way to do this, but I was curious how this would work)
predict_majority <- function(model_list, spec){
  nsims <- length(model_list)
  pred_matrix <- matrix(
    NA_character_,
    nrow = nrow(spec),
    ncol = nsims
  )
  for(i in seq_len(nsims)){
    pred_matrix[, i] <-
      as.character(
        predict(
          model_list[[i]],
          spec
        )
      )
  }
  majority <- apply(
    pred_matrix,
    1,
    function(x){
      tab <- table(x)
      names(tab)[which.max(tab)]
    }
  )
  return(as.factor(majority))
}

print("Predicting All")
#Make prediction
pred_all <- predict_majority(mods_all, spec)

print("Predicting Vegitation")
pred_veg <- predict_majority(mods_veg, spec)

print("Predicting Non Vegitation")
pred_nonveg <- predict_majority(mods_nonveg, spec)


nrow(spec)
class(pred_all)
length(pred_all)
table(pred_all)

#Save csv files
print("Finalizing results")
pred_all_df <- coords
pred_all_df$Prediction <- pred_all

pred_veg_df <- coords
pred_veg_df$Prediction <- pred_veg

pred_nonveg_df <- coords
pred_nonveg_df$Prediction <- pred_nonveg

print("Saving Data...")
write.csv(pred_all_df, pred_all_fp, row.names = FALSE)
write.csv(pred_veg_df, pred_veg_fp, row.names = FALSE)
write.csv(pred_nonveg_df, pred_nonveg_fp, row.names = FALSE)

#Map Making!!
#Set up how we want them to look
print("Plotting Data")
plot_prediction <- function(df, outfile, title){
  p <- ggplot(df,aes(x = X_EPSG_32612, y = Y_EPSG_32612, color = Prediction)) +
    geom_point(size = 0.35) +
    coord_equal() +
    theme_bw() +
    labs(title = title, x = "UTM Easting (m)", y = "UTM Northing (m)", color = "Prediction")
  ggsave(filename = outfile, plot = p, width = 8, height = 6, dpi = 300)
}

print("Creating Maps")
#Create the maps
plot_prediction(pred_all_df, map_all_fp, "PLSDA_Prediction_All_Traing")
plot_prediction(pred_veg_df, map_veg_fp, "PLSDA_Prediction_Veg_Training")
plot_prediction(pred_nonveg_df, map_nonveg_fp, "PLSDA_Prediction_NonVeg_Training")

print("Done!")
