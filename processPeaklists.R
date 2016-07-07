#' This function performs preprocessing on MS peaklists using the MALDIQuant package
#' Note: This script performs data preprocessing on all the MS peaklists from different bacterial samples and their replicates and generates a feature matrix with the m/z list in the first column of the matrix.
#' @title processPeaklists.R
#' @return Writes a feature matrix to a .csv file

## Date 4 July 2016
## Author: Purva Kulkarni

library(MALDIquant)
library(MALDIquantForeign)

setwd("/Users/localadm/PhD_Work/Server/Lab_proceedings_Current_Folder/Work/Paolina_Riya-LESA/mzXML 190416/Finalized_Analysis/")

processPeaklists <- function()
{
  # Read the directory containing the peaklists in .csv format
  peaklistObjects <- importCsv('./MS_Spectra/', centroided = FALSE)
  
  # Create a canvas to display the plots 
  par(mfrow=c(3,3))

  # Display plot for replicate 1
  plot(peaklistObjects[[1]], main = "Raw data - Replicate 1")
  
  # Tranform intensities using the selected method for variance stabolization and display
  peaklistObjectsTransformed <- transformIntensity(peaklistObjects, method = "sqrt")
  plot(peaklistObjectsTransformed[[1]], main = "Sqrt-transformation")
  
  # Smooth intensities using the selected method and display
  peaklistObjectsTransformedSmoothed <- smoothIntensity(peaklistObjectsTransformed, method = "SavitzkyGolay", halfWindowSize=10)
  plot(peaklistObjectsTransformedSmoothed[[1]], main = "Smoothing")
  
  # Estimate baseline, plot it and subtract it from the data and plot again
  baseline <- estimateBaseline(peaklistObjectsTransformedSmoothed[[1]], method = "TopHat", halfWindowSize=75)
  plot(peaklistObjectsTransformedSmoothed[[1]], main = "Baseline estimation")
  lines(baseline, col = "red", lwd = 2)
  peaklistObjectsTransformedSmoothedBaseline <- removeBaseline(peaklistObjectsTransformedSmoothed, method = "TopHat", halfWindowSize=75)
  plot(peaklistObjectsTransformedSmoothedBaseline[[1]], main = "Baseline correction")
  
  # Equalize intensity values using intesity calibration
  peaklistObjectsTransformedSmoothedBaselineCalibrated <- calibrateIntensity(peaklistObjectsTransformedSmoothedBaseline, method = "TIC")

  # Perfrom spectral alignment and plot
  peaklistObjectsTransformedSmoothedBaselineCalibratedAligned <- alignSpectra(peaklistObjectsTransformedSmoothedBaselineCalibrated, halfWindowSize = 20, SNR = 2, tolerance = 0.02, warpingMethod = "lowess")
  plot(peaklistObjectsTransformedSmoothedBaselineCalibratedAligned[[1]], main = "Calibration and alignment")

  # Average mass spectra/ technical replicates and plot
  averageSpectra <- averageMassSpectra(peaklistObjectsTransformedSmoothedBaselineCalibratedAligned, method = "mean")
  plot(averageSpectra, main = "Averaged spectra")
  
  # Estimate noise in average spectra and plot 
  noise <- estimateNoise(averageSpectra)
  #plot(averageSpectra, xlim = c(600,800), ylim = c(0,0.02), main = "Estimated noise")
  #lines(noise, col = "red")
  #lines(noise[,1], noise[,2]*2, col = "blue")
  
  # Detect peaks in mass peak objects and plot one of the mass peak objects
  peaks <- detectPeaks(peaklistObjectsTransformedSmoothedBaselineCalibratedAligned, method = "MAD", halfWindowSize = 20, SNR = 2)
  plot(peaklistObjectsTransformedSmoothedBaselineCalibratedAligned[[1]], xlim = c(600,800), ylim = c(0,0.02), main = "Detected peaks")
  points(peaks[[1]], col = "red", pch = 18)
  
  # Find similar peaks across mass peak obbjects using peak binning and filter them by removing infrequently occuring peaks
  peaks <- binPeaks(peaks, tolerance=0.002)
  peaks <- filterPeaks(peaks, minFrequency=0.25)
  plot(peaks[[1]], main = "Binned and filtered")
  
  # Generate feature matrix
  # Returns a matrix containing intensities of all MassPeaks objects of peaks
  # If a peak is missing the corresponding intensity value of the spectrum is used. If spectra is missing NA is used instead
  # There is an additional attribute "mass" that stores the mass values.
  featureMatrix <- intensityMatrix(peaks, peaklistObjectsTransformedSmoothedBaselineCalibratedAligned)
  
  # Save feature matrix as csv file
  write.csv(t(featureMatrix), file = "FeatureMatrix.csv")
}