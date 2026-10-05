
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
rfl_list <- file.path("./R_input/Moab/f110512t01p00r11rdn_c_sc01_ort_img")
shp_fp    <- file.path("R_input/Shapes/Final_map.shp") #*load the shapefile for the Henry map **yellowstone shape**
#=======================================================================================================

#                                                       Step 1: Data Conversion
#                                                         & Shape Crop
print("Converting Data")
#=======================================================================================================
  # Create empty variables
rfl_dat <- NULL

  # add all rasters to a list
  # Add values to another list
for(rfl in rfl_list){
  print("New Raster")
  rfl_new <- rast(rfl)
  
  print("Adding new Raster to collection")
  rfl_dat <- c(rfl_dat, rfl_new)
  
  print("Converting to Spatial Point")
  crs(rfl_new) <- "EPSG:32612"
  
}

#   Shape Details

print("Extracting shape and cropping")
#=======================================================================================================
print("Creating the shape from file")
units <- vect(shp_fp) 

units <- project(units, "EPSG:32612")

uni <- units["MapUnit"]

print("Shape converted to EPSG:32612")
print("Plotting to verify")

# verify shapefile visually
plot(uni) 

# add in the codes while removing the " "
uni$Code <- as.factor(uni$MapUnit)

#Create the buffer
buffer_dist <- buffer(units, width = -30)

spec = NULL

# set training Area
print("Extracting data within established bounds")

for(i in seq_along(rfl_dat)){
  
  spec <- extract(rfl_dat[[i]], uni, cells = TRUE, xy=TRUE)

  xy <- xyFromCell(rfl_dat[[i]], spec$cell)

  spec_all<- cbind(xy, spec)
  
}

# check if there's overlap, the red is showing where the emit pixel is
print("Plotting for verification")
plot(uni)

print("Creating spacial points for graphing and plotting overlap")

points <- subset(spec, select = c("x","y")) # extract the coordinates only
poi_UTM <- vect(points, geom = c("x", "y"), crs = "EPSG:32612")
plot(poi_UTM, cex=2, col="red", pch=16, add = T)

print("Adding Codes to raster")
empty_rast <- rast(ext(uni), resolution = 14.8, crs="EPSG:32612")
uni_ras <- rasterize(uni, empty_rast, field="Code")

print("Plotting to verify")
plot(uni_ras)

print("Extracting and compiling data")
# Extract data
dat <- extract(uni_ras, poi_UTM, ID = FALSE)

# Create lookup table for code names: check if ok
look <- as.data.frame(matrix(levels(uni$Code), nrow=length(levels(uni$Code)))) # code names from shapefile
names(look)<- "Code"
look$ID <- seq(1:nrow(look)) # assign id to each code

# Data wrangling

dat$ID <- match(dat$Code, look$Code)
names(dat$Code) <- "Code"
dada <- cbind(dat$Code, spec)
names(dada)[names(dada)=="dat$Code"]<-"Code"
dada <- dada %>% dplyr::select(ID,x,y,Code,everything())
names(dada) <- gsub(" Nanometers", "",names(dada))
# Count number of EMIT pixels per unit and save it
pixel_counts <- dada%>% group_by(Code) %>% summarise(Pixels = n())

print(pixel_counts)

# Clean up 
rm(dat, look, poi_UTM, points, uni, uni_ras, units, pixel_counts, buffer_dist)
dada<-dada[,-1]

#=======================================================================================================
#=======================================================================================================

##################### Works up until here ###################################################

#   Step 3: NDVI
print("NVDI")
#=======================================================================================================
print("Establishing the Vegetation and Non tables")
# Identify bands to use for NDVI
red_band <- which.min(abs(as.numeric(gsub(" Nanometers", "",names(dada[5:(length(spec)-1)]))) - 660))
nir_band <- which.min(abs(as.numeric(gsub(" Nanometers", "",names(dada[5:(length(spec)-1)]))) - 860))

red <- dada[, 4 + red_band]
nir <- dada[, 4 + nir_band]
ndvi <- (nir - red) / (nir + red)

# Mask non-veg
spec_nonveg <- dada[ndvi <= 0.35, ]
spec_veg    <- dada[ndvi > 0.35, ]

rm(red_band,
   nir_band,
   red,
   nir,
   ndvi
)
#=======================================================================================================

#-------------------------------------------------------------------------------------------------------
#=======================================================================================================

# # Step 5: Merger Equal and Threshold
print("Merging Data and establishing thresholds")
# #=======================================================================================================
#               ###Script 5, Merger weighs all training units equally
print("Vegetation column created and set to appropriate rows")
# Create Data Sets
spec_nonveg$VM <- "nonveg"
spec_veg$VM    <- "veg"

print("Data re-merged")
# Merge data

spec_all <- bind_rows(spec_nonveg, spec_veg)

#   Step 2: Consolidation of Units
print("Consolidating Units and Data")
#=======================================================================================================

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
print("Map Units being added to the data table")
#Join the names to the table
spec_nonveg <- left_join(spec_nonveg, unit_lookup, by = "Code")
spec_veg <- left_join(spec_veg, unit_lookup, by = "Code")
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

#=======================================================================================================

# Step 6: PLSDA

#=======================================================================================================
#     Formula
#-------------------------------------------------------------------------------------------------------
print("PLSDA")

#----------# All #-----------#

# PDF Output
pdf_fp    <- file.path("R_output/All", "PLSDA_kappas_all_thresh.pdf")

#Testing pixel count
print(table(spec_all$Code))
print(nrow(spec_all))
balanced <- sample_dataset(dada)

training_summary <- as.data.frame(table(balanced$Code))
names(training_summary) <- c("Code", "TrainingPixels")
training_summary$Dataset <- "RandUnitCoord_all_thresh.csv"

print("PLSDA training starting!")

print("Defining Classes")
classi <- as.factor(balanced$Code) # define classes
is.wave <- suppressWarnings(!is.na(as.numeric(colnames(balanced))))

print("Removing uneccesary data...")
keep <- !is.na(classi) & classi != "Unknown" & classi != "water"

spec_temp <- balanced[keep, ]
classi <- droplevels(classi[keep])

## Preparation ##

print("Getting Wavelengths")
wvl  <- colnames(spec_temp)[is.wave]
spec <- spec_temp[, wvl]

## Partition  ##
set.seed(1080) # for consistency
inTrain <- createDataPartition(y =classi, p = .7, list = FALSE) # e.g. 70% spectra per class
training <- spec[inTrain,seq(1,ncol(spec), by=5)] # resampling optional
testing <- spec[-inTrain,seq(1,ncol(spec), by=5)] # pick all that was not picked in inTrain
trainclass <- classi[inTrain] # define units for training
testclass <- classi[-inTrain] # define units for testing

set.seed(1840) ### create random number for obs in each Code/class 

table(testclass) # verify
table(trainclass)
table(spec_temp$Code)
## Finding components  ## 

print("Training...")
mods <- list()
for (nsim in seq(nsims)){
  inTrain <- createDataPartition(y =classi, p = 0.7, list = FALSE) ## or create partitioning using fractions, e.g. 70:30
  print(nsim)
  flush.console()
  set.seed(nsim)
  traini <- spec[inTrain,seq(1,length(wvl), by=5)] 
  testi <- spec[-inTrain,seq(1,length(wvl), by=5)]
  trainclass <- classi[inTrain]
  testclass <- classi[!(inTrain)]
  plsFit <- train(traini, trainclass, method = "simpls", tuneLength = k, 
                  trControl = trainControl(method="cv", number = 5, 
                                           classProbs = TRUE, 
                                           savePredictions = "final")
                  ) 
  mods[[nsim]] <- plsFit
}

table(testclass) # verify
table(trainclass)
table(balanced$Code)
## Finding components  ## 

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
print("Training completed. Compiling Data...")
kapp <- as.data.frame(as.numeric(t(kappas[,-1]))) # combines all kappas
kapp <- cbind(kapp, rep(1:k, each=length(mods))) # labels kappas per model
names(kapp) <- c("Kappa", "ncomps") # nicer labeling

kapp$ncomps <- as.factor(kapp$ncomps) # change to factor

modi <- lm(Kappa~ncomps, kapp) # combine and reformat data
tuk <- HSD.test(modi,"ncomps") # perform tukey test

tuk_dat <- as.data.frame(tuk$groups) # convert tukey test into dataframe
tuk_dat$var <- as.numeric(row.names(tuk_dat)) # add the component column
tuk_dat_all <- tuk_dat[order(tuk_dat$var,decreasing = F),] # reorder from small to large

letters <- as.character(tuk_dat$groups) # constrain order of latter for graph

## Kappa Tukey plot ## 

pdf(pdf_fp, width = 5, height = 4)
par(bty = "l")
boxplot(kapp$Kappa ~ kapp$ncomps, ylim = c(0,1),
        xlab = "Number of components", ylab = "Kappa")
text(x = 1:k, y = rep(1, k), letters)


dev.off()

#=======================================================================================================
# Step 7: Extractor
#=======================================================================================================
#     Function
#-------------------------------------------------------------------------------------------------------

print("Extracting Data and Finalizing PLSDA")
PLSDA_2 <- function(balanced, out_finmods, out_scores, out_probmean, out_probplot, 
                    out_confumean,out_confuperc, corrplot, out_loadings, nsims){

  print("Optimizing raster...")
  is.wave <- suppressWarnings(!is.na(as.numeric(colnames(balanced))))
  wvl <- colnames(balanced)[is.wave] # define wvl range for spectral matrix, check your column names
  spec <- balanced[, wvl]
  
  ## Preperation ##
  # define classes
  
  classi <- as.factor(balanced$Code)
  nwaves <- 1
  #compi <- 8 ### select number of components based on Kappa plot / Tukey test from pt 1
  compi <- 2:20
  finmods <- list() #
  
  print("Removing uneccesary data...")
  # Remove unlabeled pixels
  keep <- !(is.na(classi))
  spec <- spec[keep, ]
  classi <- classi[keep]

  spec <- data.frame(spec)
  
  for (nsim in 1:nsims){
    print(nsim)
    flush.console()
    set.seed(nsim)
    #inTrain <- rndid[[nsim]]<= ??? # number train per code
    inTrain <- createDataPartition(y =classi, p = 0.7, list = FALSE)
    training <- spec[inTrain,seq(1,length(wvl), by=nwaves)]
    trainclass <- as.factor(classi[inTrain])
    
    finalModel <- plsda(training, trainclass, ncomp=compi, method = "simpls")
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
    
    print("Training: Pt2...")
    print("Creating Partitions...")
    inTrain <- createDataPartition(y = classi, p = 0.7, list = FALSE)
    testing <- spec[-inTrain,seq(1,length(wvl), by=nwaves)]
    testclass <- as.factor(classi[-inTrain])
    
    print("Establishing Probabilities...")
    plsProbs <- predict(finmods[[nsim]], newdata = testing, type = "prob")
    
    print("Establishing PLS Classes")
    plsClasses <- predict(finmods[[nsim]], newdata = testing)
    
    print("Creating Confusion Matrix")
    confus[[nsim]] <- confusionMatrix(data = plsClasses, testclass)
    
    print("Converting Probabilities to Data Frame")
    probs <- as.data.frame(plsProbs)
    
    print("Converting Names into proper strings")
    names(probs) <- sapply(strsplit(names(probs),split = ".n"),"[",1)
    
    print("Binding Data")
    probs <- cbind(testclass, probs)
    
    print("Compiling Data")
    probis[[nsim]] <- probs 
  }
  
  ## Probability mean ##
  print("Finding Mean...")
  arr <- array(unlist(probis), dim = c(dim(probis[[1]]), nsims))
  prob_mean <- apply(arr, 1:2, mean)
  prob_mean <- as.data.frame(prob_mean)
  colnames(prob_mean) <- colnames(probis[[1]]) #was 98
  prob_mean$testclass <- probis[[1]]$testclass #was 97
  write.csv(prob_mean, out_probmean, row.names = FALSE)
  
  ## Probability Plot ##
  print("Plotting...")
  pp <- melt(prob_mean, id = "testclass")
  pp$position <- ifelse(pp$testclass == gsub(".8 comps$","",pp$variable), 2, 1)
  pp$testclass <- factor(pp$testclass, levels = rev(levels(pp$testclass)))
  
  pdf(out_probplot, width = 6, height = 4)
  ggplot(pp, aes(x = testclass, y = value, fill = variable, group = position)) +
    geom_bar(stat = "identity", position = "fill") +
    theme(legend.title = element_blank()) +
    labs(x = "Code", y = "") +
    coord_flip()
  dev.off()
  
  print("Creating Confusion Matrix...")
  ## Confusion matrix! ##
  tabs <- lapply(confus, function(x) x$table)
  tabsi <- Reduce('+', tabs)
  tab_mean <- as.data.frame.matrix(tabsi / length(confus))
  write.csv(tab_mean, out_confumean, row.names = TRUE)
  
  sums <- colSums(tab_mean)
  tabs_perc <- sweep(tab_mean, 2, sums, FUN = "/")
  write.csv(tabs_perc, out_confuperc, row.names = TRUE)
  
  print("Plotting the matrix")
  # Confu plot
  pdf(out_corrplot, width = 7, height = 6, pointsize = 9)
  corrplot::corrplot(as.matrix(tabs_perc), addCoef.col = 1, tl.srt = 70,
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
  print("Compiling Statistics")
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

ssum_all <-  PLSDA_2(balanced, out_finmods, out_scores, out_probmean, out_probplot, 
                     out_confumean,out_confuperc, corrplot, out_loadings, nsims)


write.csv(ssum_all, out_summary, row.names = FALSE)

#=======================================================================================================
# Step 8: Extractor
# #=======================================================================================================
###Script 8: Compare the models: This is a complete rewrite, so good luck to us all

#Input files

#input_fp <- file.path("R_output", "UnitCoordSpec_all.csv")
model_all_fp <- file.path("R_output/All/", "PLSDA_finmods_all_thresh.rds")
#model_veg_fp <- file.path("R_output/Veg/", "PLSDA_finmods_veg_thresh.rds")
#model_nonveg_fp <- file.path("R_output/Nonveg/", "PLSDA_finmods_nonveg_thresh.rds")

#Set output files
pred_all_fp <- file.path("R_output/Final_Outputs/", "PredUnitCoord_all.csv")
#pred_veg_fp <- file.path("R_output/Final_Outputs/", "PredUnitCoord_veg.csv")
#pred_nonveg_fp <- file.path("R_output/Final_Outputs/", "PredUnitCoord_nonveg.csv")

map_all_fp <- file.path("R_output/Final_Outputs/", "PredictionMap_all.png")
#map_veg_fp <- file.path("R_output/Final_Outputs/", "PredictionMap_veg.png")
#map_nonveg_fp <- file.path("R_output/Final_Outputs/", "PredictionMap_nonveg.png")

#Read data
dati <- dada
coords <- dati[, c("Code", "x", "y")]
is.wave <- suppressWarnings(!is.na(as.numeric(colnames(dati))))
spec <- dati[, is.wave]

#Load models
mods_all <- readRDS(model_all_fp)
#mods_veg <- readRDS(model_veg_fp)
#mods_nonveg <- readRDS(model_nonveg_fp)


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

#print("Predicting Vegitation")
#pred_veg <- predict_majority(mods_veg, spec)

#print("Predicting Non Vegitation")
#pred_nonveg <- predict_majority(mods_nonveg, spec)


nrow(spec)
class(pred_all)
length(pred_all)
table(pred_all)

#Save csv files
print("Finalizing results")
pred_all_df <- coords
pred_all_df$Prediction <- pred_all

# pred_veg_df <- coords
# pred_veg_df$Prediction <- pred_veg
# 
# pred_nonveg_df <- coords
# pred_nonveg_df$Prediction <- pred_nonveg

print("Saving Data...")
write.csv(pred_all_df, pred_all_fp, row.names = FALSE)
# write.csv(pred_veg_df, pred_veg_fp, row.names = FALSE)
# write.csv(pred_nonveg_df, pred_nonveg_fp, row.names = FALSE)

#Map Making!!
#Set up how we want them to look
print("Plotting Data")
plot_prediction <- function(df, outfile, title){
  p <- ggplot(df,aes(x = x, y = y, color = Prediction)) +
    geom_point(size = 0.35) +
    coord_equal() +
    theme_bw() +
    labs(title = title, x = "UTM Easting (m)", y = "UTM Northing (m)", color = "Prediction")
  ggsave(filename = outfile, plot = p, width = 8, height = 6, dpi = 300)
}

print("Creating Maps")
#Create the maps
plot_prediction(pred_all_df, map_all_fp, "PLSDA_Prediction_All_Traing")
# plot_prediction(pred_veg_df, map_veg_fp, "PLSDA_Prediction_Veg_Training")
# plot_prediction(pred_nonveg_df, map_nonveg_fp, "PLSDA_Prediction_NonVeg_Training")

print("Done!")
