# multi_gauss_fit.R
# An R script to perform multi-modal Gaussian fits

# Get the arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
   cat("usage: filename value_col ncomps\n")
   quit()
}

filename <- args[1]
value_col <- as.integer(args[2])
ncomps <- as.integer(args[3])

# Plot histogram or not
plot_hist = FALSE

# Load the package
if ("mixtools" %in% rownames(installed.packages()) == FALSE) {
   cat("Missing required package 'mixtools'!\n")
   quit()
}
suppressMessages(require(mixtools))

# Load the data file
data <- read.csv(filename, sep=' ')
x <- data[[value_col]]

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
  if ("rootSolve" %in% rownames(installed.packages()) == FALSE) {
    cat("Missing required package 'rootSolve'!\n")
    quit()
  }
  suppressMessages(require(rootSolve))

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
  hist(x, freq = FALSE)
  xx <- seq(1.0,1.025,0.0001)
  lines(xx, weights[1]*dnorm(xx, mean=means[1], sd=stdevs[1]), type='l')
  lines(xx, weights[2]*dnorm(xx, mean=means[2], sd=stdevs[2]), type='l')
  prompt <- "hit spacebar to clsoe plots"
  capture <- tk_messageBox(message=prompt)
}
