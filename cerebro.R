# alorenzetti 20190125
# version 0.2

# declaring package requirements
libraries = c("EBImage", "plyr", "ggplot2")

# loading packages 
invisible(lapply(libraries, function(x){library(x, character.only = T)}))

# getting parameters
args=commandArgs(trailingOnly = T)
args=c("100", "0.8", "0.5")

# area of colony must be at least area px
area = 100
area = args[1]
# circularity must be greater equal than circ
circ = 0.8
circ = args[2]
# eccentricity must be lesser than eccen
eccen = 0.5
eccen = args[3]

# defining alternative display function
display2 = function(x){ display(x, method="raster") }

# getting list of files
files = list.files("croppedPlates", "*.tif", full.names = T)

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
  circularity = fts[,"s.area"] / fts[,"s.perimeter"]^2 *4*pi
  
  # creating a new dataframe with features
  features = as.data.frame(cbind(fts, fts2, fts3, circularity), stringsAsFactors = F)
  features$plate = rep(prefix, dim(features)[1])
  
  # filtering features by area, circularity and eccentricity
  filter = features[,"s.area"] > area & features[,"circularity"] >= circ & features[,"m.eccentricity"] < eccen
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
  circularity = fts[,"s.area"] / fts[,"s.perimeter"]^2 *4*pi
  
  # creating a new dataframe with features
  features = as.data.frame(cbind(fts, fts2, fts3, circularity), stringsAsFactors = F)
  features$plate = rep(prefix, dim(features)[1])
  
  # filtering features by area, circularity and eccentricity
  filter = features[,"s.area"] > area & features[,"circularity"] >= circ & features[,"m.eccentricity"] < eccen
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
featureTable = NULL
for(i in files){
  featureTable = rbind(featureTable, getFeatureTable(i))
}

write.table(featureTable, "features.txt", sep="\t", col.names=T, row.names=F, quote=F)
intensityThr = quantile(featureTable$b.mean, probs = 0.99)
intensitySdThr = quantile(featureTable$b.sd, probs = 0.01)
mutantIndex = featureTable$b.mean >= intensityThr & featureTable$b.sd <= intensitySdThr
featureTable$class = NA
featureTable$class[mutantIndex] = "Mutant"
featureTable$class[!mutantIndex] = "Normal"

# running analysis to classify mutants and
# write images
if(!dir.exists("resultsR")){
  dir.create("resultsR")
  
  results = NULL
  for(i in files){
    line = workflow(i, "resultsR", intensityThr, intensitySdThr)
    results = rbind(results, line)
  }
}

# making a copy of results
df = results

# adjusting ouput matrix
results = data.frame(plateID = sub("^.*/(.*)\\..*$", "\\1", as.character(results[,3])),
                     strain = sub("^.*/(.*?)_.*$", "\\1", as.character(results[,3])),
                     replicate = sub("^.*/.*_(.*_.*).tif", "\\1", as.character(results[,3])),
                     counts = as.numeric(results[,1]),
                     mutcounts = as.numeric(results[,2]),
                     row.names=NULL)

write.table(results, "generalAndMutCounts.txt",
            row.names = F, col.names = T, quote = F, sep="\t")

# analyzing the data and creating charts
theme_set(theme_bw())
svg("colony_redIntensity_scatter.svg", height = 4, width = 5)
ggplot(featureTable, aes(x=b.mean, y=b.sd, colour=class)) +
  scale_colour_manual(values = c("Mutant" = "red", "Normal" = "green"),
                      guide = guide_legend(title = "Class")) +
  geom_point(alpha=0.3) +
  geom_vline(xintercept = intensityThr) +
  geom_hline(yintercept = intensitySdThr) +
  xlab("Colony Mean Red Intensity") +
  ylab("Colony Standard Deviation of Red Intensity")
dev.off()