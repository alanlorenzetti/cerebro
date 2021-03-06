# alorenzetti 20190128
# version 0.3

# declaring package requirements
libraries = c("EBImage", "ggplot2")

# loading packages 
invisible(lapply(libraries, function(x){library(x, character.only = T)}))

# getting parameters
#args=c("images", "processedImages", "results", "100", "0.8", "0.5")
args=commandArgs(trailingOnly=T)

# control input
helpMessage="Please, check README.md file for instructions.\n\nUsage:\nRscript cerebro.R <imageInputDir> <imageOutputDir> <resultsDir> <area> <circularity> <eccentricity>\n\nExample:\nRscript cerebro.R images processedImages results 100 0.8 0.5\n\nCerebro: Spontaneous Mutant Finder\nhttps://github.com/alanlorenzetti/cerebro\n\n"
if(args[1] == "--help" | args[1] == "-help" | args[1] == "-h"){stop(helpMessage)}
if(length(args) != 6){stop("You must provide all the six arguments with the correct order. Please, check README.md for more information.")}

## setting variables
# input dir
inputDir=args[1]

# processed image dir
imgOutputDir=args[2]

# results dir
resultsDir=args[3]

# area of colony must be at least area px
area = as.numeric(args[4])

# circularity must be greater equal than circ
circ = as.numeric(args[5])

# eccentricity must be lesser than eccen
eccen = as.numeric(args[6])

# creating output directories
if(!dir.exists(imgOutputDir)){
  dir.create(imgOutputDir)
}else{stop(paste(imgOutputDir, "directory already exists. Aborting."))}

if(!dir.exists(resultsDir)){
  dir.create(resultsDir)
}else{stop(paste(resultsDir, "directory already exists. Aborting."))}

# defining alternative display function
display2 = function(x){ display(x, method="raster") }

# getting list of files
files = list.files(inputDir, "*.tif", full.names = T)

# defining function to get features for every plate and then compute parameters
getFeatureTable = function(imgIN = imgIN){
  # setting prefix for filenames
  prefix = sub(".tif$", "", imgIN)
  prefix = sub("^.*/", "", prefix)
  
  # reading img
  imgOriginal = readImage(imgIN)
  
  # image rgb
  imgRGB = channel(imgOriginal, mode = "rgb")
  
  # split channel RED (mutants are almost vanished in this one)
  imgRed = channel(imgRGB, mode = "red")
  
  # grayscale image
  imgGray = channel(imgOriginal, "gray")
  
  ## removing borders
  # length of vector corresponding to dark/blank field
  # in this way the first pixel will tell us which value correspond to borders
  zeros = which(imgGray == imageData(imgGray)[[1]])
  l = length(zeros)
  
  # creating rnormal values to replace dark field
  bgcenter = median(imgGray[-zeros])
  bgsd = sd(imgGray[-zeros]) * 0.05
  bgvec = rnorm(l, bgcenter, bgsd)
  imgGray[zeros] = bgvec
  
  # computing max of kde to make background more homogeneous
  x = density(imgGray)
  thr = x$x[which(x$y == max(x$y))]
  
  # apply median only to values equal or above the mean
  imgGray[which(imgGray >= thr)] = thr
  
  # inverting img
  imgGray = 1 - imgGray
  
  # apply adaptive thresholding to binarize img
  imgBin = thresh(imgGray, w = 25, h = 25, offset = 0.01)
  
  # applying segmentation and creating mask
  imgSeg = bwlabel(imgBin)
  
  # filling holes to avoid doughnut shapes
  imgSeg = fillHull(imgSeg)
  
  # computing features for segmented objects
  fts = computeFeatures.shape(imgSeg)
  fts2 = computeFeatures.basic(imgSeg, ref = imgRed)
  fts3 = computeFeatures.moment(imgSeg, ref = imgGray)
  circularity = (fts[,"s.area"] / fts[,"s.perimeter"]^2)*4*pi
  
  # creating a new dataframe with features
  features = as.data.frame(cbind(fts, fts2, fts3, circularity), stringsAsFactors = F)
  features$plate = rep(prefix, dim(features)[1])
  
  # filtering features by area, circularity and eccentricity
  filter = features[,"s.area"] >= area & features[,"circularity"] >= circ & features[,"m.eccentricity"] < eccen
  featuresBad = rownames(features[!filter,])
  features = features[filter,]
  features$ID = paste0("colony_", 1:dim(features)[1])
  
  return(features)
}

# defining function for image processing
workflow = function(imgIN = imgIN, dir = dir, intensityThr = intensityThr, intensitySdThr = intensitySdThr){
  # setting prefix for filenames
  prefix = sub(".tif$", "", imgIN)
  prefix = sub("^.*/", "", prefix)
  
  # echoing progress
  print(paste("Processing", prefix))  

  # reading img
  imgOriginal = readImage(imgIN)
  
  # image rgb
  imgRGB = channel(imgOriginal, mode = "rgb")
  
  # split channel RED (mutants are almost vanished in this one)
  imgRed = channel(imgRGB, mode = "red")
  
  # grayscale image
  imgGray = channel(imgOriginal, "gray")
  
  ## removing borders
  # length of vector corresponding to dark/blank field
  # in this way the first pixel will tell us which value correspond to borders
  zeros = which(imgGray == imageData(imgGray)[[1]])
  l = length(zeros)
  
  # creating rnormal values to replace dark field
  bgcenter = median(imgGray[-zeros])
  bgsd = sd(imgGray[-zeros]) * 0.05
  bgvec = rnorm(l, bgcenter, bgsd)
  imgGray[zeros] = bgvec

  # computing max of kde to make background more homogeneous
  x = density(imgGray)
  thr = x$x[which(x$y == max(x$y))]

  # apply median only to values equal or above the mean
  imgGray[which(imgGray >= thr)] = thr
  
  # inverting img
  imgGray = 1 - imgGray
  
  # apply adaptive thresholding to binarize img
  imgBin = thresh(imgGray, w = 25, h = 25, offset = 0.01)
  
  # applying segmentation and creating mask
  imgSeg = bwlabel(imgBin)
  
  # filling holes to avoid doughnut shapes
  imgSeg = fillHull(imgSeg)
  
  # computing features for segmented objects
  fts = computeFeatures.shape(imgSeg)
  fts2 = computeFeatures.basic(imgSeg, ref = imgRed)
  fts3 = computeFeatures.moment(imgSeg, ref = imgGray)
  circularity = (fts[,"s.area"] / fts[,"s.perimeter"]^2)*4*pi
  
  # creating a new dataframe with features
  features = as.data.frame(cbind(fts, fts2, fts3, circularity), stringsAsFactors = F)
  features$plate = rep(prefix, dim(features)[1])
  
  # filtering features by area, circularity and eccentricity
  filter = features[,"s.area"] >= area & features[,"circularity"] >= circ & features[,"m.eccentricity"] < eccen
  featuresBad = rownames(features[!filter,])
  features = features[filter,]
  features$ID = paste0("colony_", 1:dim(features)[1])
  
  # removing objects not satisfying filter
  imgFinal = rmObjects(imgSeg, featuresBad, reenumerate = T)
  
  # showing valid objects in gray scale processed image
  imgPainted = paintObjects(imgFinal, channel(imgGray, mode="rgb"), opac = 1, thick = T, col = "green")
  
  # counting final image objects
  ncolonies = max(imgFinal)
  
  # finding mutants
  RedIntensity = intensityThr
  RedIntensitySd = intensitySdThr
  mutants = which(features$b.mean >= RedIntensity & features$b.sd <= RedIntensitySd)
  nonMutants = which(features$b.mean < RedIntensity & features$b.sd > RedIntensitySd)
  imgMutants = rmObjects(imgFinal, nonMutants)
  
  # showing valid mutants
  imgPainted = paintObjects(imgMutants, imgPainted, opac = 1, thick = T, col = "red")
  
  # counting mutants
  nmutants = max(imgMutants)
  
  # saving images to one file
  outfile = paste0(dir, "/", prefix, ".tif")
  
  png(outfile, width = dim(imgBin)[2] * 2, height = dim(imgBin)[1] * 1, res = 600, units = "px")
  par(mfrow = c(1,2), bg = "black")
  display2(imgOriginal)
  display2(imgPainted)
  dev.off()
  
  # printing results
  return(c(ncolonies, nmutants, imgIN))
}

# running function to get parameters to classify mutants
print("Computing Parameters for Classification...")
featureTable = NULL
for(i in files){
  featureTable = rbind(featureTable, getFeatureTable(i))
}
print("Done!")

write.table(featureTable, paste0(resultsDir, "/", "features.txt"), sep="\t", col.names=T, row.names=F, quote=F)
intensityThr = quantile(featureTable$b.mean, probs = 0.99)
intensitySdThr = quantile(featureTable$b.sd, probs = 0.01)
mutantIndex = featureTable$b.mean >= intensityThr & featureTable$b.sd <= intensitySdThr
featureTable$class = NA
featureTable$class[mutantIndex] = "Mutant"
featureTable$class[!mutantIndex] = "Normal"

# running analysis to classify mutants and
# write images
results = NULL
for(i in files){
    line = workflow(i, imgOutputDir, intensityThr, intensitySdThr)
    results = rbind(results, line)
}
print("Done!")

# making a copy of results
df = results

# adjusting output matrix
results = data.frame(plateID = sub("^.*/(.*)\\..*$", "\\1", as.character(results[,3])),
                     counts = as.numeric(results[,1]),
                     mutcounts = as.numeric(results[,2]),
                     row.names=NULL)

write.table(results, paste0(resultsDir, "/", "generalAndMutCounts.txt"),
            row.names = F, col.names = T, quote = F, sep="\t")

# Creating classification chart
print("Plotting classification chart")
theme_set(theme_bw())
svg(paste0(resultsDir, "/", "colony_redIntensity_scatter.svg"), height = 4, width = 5)
ggplot(featureTable, aes(x=b.mean, y=b.sd, colour=class)) +
  scale_colour_manual(values = c("Mutant" = "red", "Normal" = "green"),
                      guide = guide_legend(title = "Class")) +
  geom_point(alpha=0.3) +
  geom_vline(xintercept = intensityThr) +
  geom_hline(yintercept = intensitySdThr) +
  xlab("Colony Mean Red Intensity") +
  ylab("Colony Standard Deviation of Red Intensity")
dev.off()
print("Done!")

