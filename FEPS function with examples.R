# These codes were tested in R 4.2.3 and on a Windows 10 computer.
# Examples:
# install.packages("signal")
# install.packages("dplR") 
library(signal)
library(dplR)
data(ca533)
feps(ca533,4,10,stc=c(3,2,1))  # FEPS value of 10-year low-frequency components, using a 4th-order Butterworth filter
fstats(ca533,4,20)  # stats, including FEPS (still named eps), of 20-year low-frequency components 

#### Codes (functions) used for calculate FEPS ####
# The "butter" function from the "signal" package is used instead of "pass.filt" 
# from the "dplR" package because it runs faster and is suitable for handling 
# large amounts of data or Monte Carlo simulations.

feps <- function(rwl,n,cutperiod,stc=NULL){
  # calculate feps from raw ring widths
  # rwl is data.frame, n is filter order of butterworth
  # cutperiod is cutoff period, unit is year
  # stc is a vector of three integral values or character string "auto", 
  # more information see help of "read.ids"
  # extract eps from stats
  stats_lowfreq <- fstats(rwl,n,cutperiod,stc=NULL)
  eps_lowfreq <- stats_lowfreq$eps
  eps_lowfreq
}

fstats <- function(rwl,n,cutperiod,stc=NULL){
  # calculate descriptive statistics of raw ring
  # rwl is a data.frame, n is filter order of butterworth
  # for example:
  # stats of a 10-years low-frequency components, core IDs like "KTS02A", "KTS02B"
  # stats_decade <- fstats(rwl,10,stc = c(0, 5, 1))
  rwi_low <- rwilow(rwl, n, cutperiod)
  ids <- dplR::read.ids(rwi_low,stc=stc)  
  stats <- dplR::rwi.stats(rwi_low,ids)
  stats
}

rwilow <- function(rwl, n, cutperiod){
  # rwl is a data.frame
  # if cutperiod is set to 0, then return unfiltered series with detrend only 
  if (cutperiod != 0){
    # number of cores 
    n_core <- ncol(rwl)
    # detrend use ModNegExp, can be replaced by other methods
    rwi <- dplR::detrend(rwl,method = "ModNegExp")
    # low-frequency variable initialization
    rwi_low <- rwi
    for (i in c(1:n_core)) {
      noEmpty_row <- which(!is.na(rwi[,i]))
      rwi_low[noEmpty_row,i] <- butterLow_signal(rwi[noEmpty_row,i],n,cutperiod)
      # butterLow_signal can be replaced by pass.filt function in "dplR" package
      # rwi_low[noEmpty_row,i] <- pass.filt(rwi[noEmpty_row,i], W=1/cutperiod, type="low", method="Butterworth")
    }
  } else if (cutperiod == 0){
    rwi_low <- dplR::detrend(rwl,method = "ModNegExp")
  }
  rwi_low
}

butterLow_signal <- function(x,n,cutperiod){
  # low-passing filter using "signal" package
  # sampling frequency is set to 1 for annual time series
  sample_frequency <- 1;  
  nyquist_frequency <- 0.5*sample_frequency
  # cutoff frequency
  Wn = 1/cutperiod/nyquist_frequency
  # using butter function in "signal" package, and remove the end effect  
  bf_low <- signal::butter(n, Wn, "low")
  y_low <- endEffect(bf_low,x)
  y_low
}

endEffect <- function(filt,x) {
  # end effect supression by padding series
  signal::filtfilt(filt,c(rev(x),x,rev(x)))[(length(x) + 1):(2 * length(x))]
}