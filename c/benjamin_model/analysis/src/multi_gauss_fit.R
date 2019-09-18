# multi_gauss_fit.R
# An R script to perform multi-modal Gaussian fits

# Load the required packages
if ("mixtools" %in% rownames(installed.packages()) == FALSE) {
   cat("Missing required package 'mixtools'!\n")
   quit()
}
suppressMessages(require(mixtools))
if ("rootSolve" %in% rownames(installed.packages()) == FALSE) {
  cat("Missing required package 'rootSolve'!\n")
  quit()
}
suppressMessages(require(rootSolve))

# Get the arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) {
   cat("usage: time_col value_col tstart tend val_min val_max bin_size plot_hist filename\n")
   quit()
}

time_col <- as.integer(args[1])
value_col <- as.integer(args[2])
tstart <- as.integer(args[3])
tend <- as.integer(args[4])
val_min <- as.double(args[5])
val_max <- as.double(args[6])
bin_size <- as.double(args[7])
plot_hist <- as.integer(args[8])

ncomps <- 2 # Only do bi-modal Gaussian fit (to classify modality easily)

# Load the data files
istart <- 9
x <- NULL
for (i in istart:length(args)) {
  filename <- args[i]
  data <- read.csv(filename, sep=' ', header=FALSE)

  # Filter the data for the specific time frames
  data <- data[(data[,time_col] >= tstart) & (data[,time_col] <= tend),]
  
  # Select the value column only
  if (i == istart) {
    x <- data[[value_col]]    
  } else {
    x <- c(x,data[[value_col]])
  }
}

cat("Done loading data\n")

# Do multi-modal Gaussian fitting
fit <- normalmixEM(x, k=ncomps, lambda=1.0/ncomps)

means <- fit['mu'][[1]]
stdevs <- fit['sigma'][[1]]
weights <- fit['lambda'][[1]]

cat("mean:", means, "\n")
cat("stdev:", stdevs, "\n")
cat("weight:", weights, "\n")

# Determine if the mixture is unimodal or bimodal
# Use the test from the paper by Robertson and Fryer (1969)
m <- abs(means[2]-means[1])/stdevs[1]
s <- stdevs[2]/stdevs[1]

m0 <- (2*(s^4-s^2+1)^(1.5)-(2*s^6-3*s^4-3*s^2+2))^(0.5)/s

if (m <= m0) {
  cat("modality: unimodal\n")
} else {

  a <- (s^2-1)
  b <- -m*(s^2-2)
  c <- -m^2
  d <- m*s^2
  cubic <- function(x) {
    return(a*x^3+b*x^2+c*x+d)
  }
  roots <- uniroot.all(cubic,interval=c(0,m))
  get_p <- function(y) {
    return(1.0/(1.0+y*s^3/(m-y)*exp(-0.5*y^2+0.5*((y-m)/s)^2)))
  }
  p <- sort(get_p(roots))
  if (weights[1] > p[1] && weights[1] < p[2] &&
      weights[2] > p[1] && weights[2] < p[2]) {
    cat("modality: bimodal\n")
  } else {
    cat("modality: unimodal\n")
  }
}

# Plot the histogram and the estimated guassian mixtures
if (plot_hist) {
  suppressMessages(require(tcltk))
  X11()
  h <- hist(x, freq = FALSE)
  xx <- seq(val_min,val_max,bin_size)
  lines(xx, weights[1]*dnorm(xx, mean=means[1], sd=stdevs[1]), type='l')
  lines(xx, weights[2]*dnorm(xx, mean=means[2], sd=stdevs[2]), type='l')
  prompt <- "hit spacebar to close plots"
  capture <- tk_messageBox(message=prompt)
}