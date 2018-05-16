#!/usr/local/bin/Rscript
### script to process FC data in R creating all same figures in 1 document
### example of use: Rscript FC_DensPlots_in1.R /Users/mvlkova/Documents/RecA_Promoter/FC/180222/plate1/ 12,8
### example of use2: Rscript FC_DensPlots_in1.R 180222/plate1/ 12,8
### (/Users/mvlkova/Documents/RecA_Promoter/FC/ is folder where the script is located)
### second argument is optional and defines how to partition the final pdf plot (columns,rows)
library("flowCore") ### handling FCS files
library("flowViz")  ### package for plotting
library("hexbin")   ### to plot 2D histograms
library("MASS")     ### to run dov.rob method in gating, which gives more oval shaped gates than previously used covMcd
library("grid")     ### to add text in densityplot (alternative for mtext)

### incorporate function finding all divisors
if(!exists("foo", mode="function")) source("Divisors.R")
### incorporate gating function
if(!exists("foo", mode="function")) source("FCgating.R")

args = commandArgs(trailingOnly=TRUE)     ### read arguments from command line when running the script
root.path <- c("/Users/mvlkova/Documents/RecA_Promoter/FC/")

### stop and print error message if path argument is not provided
if(length(args) == 0) {
  stop("Absolute or relative path to analyse data in must be provided as the first argument!
    Full path example: /Users/mvlkova/Documents/RecA_Promoter/FC/180222/plate1/
    Relative path example: 180222/plate1/\n")
}else{
  if(startsWith(args[1], "/")) {                                 ### if path starts with "/"
    full.path <- args[1]                                         ### set the argument as absolute path
  }else{                                                         ### else
    full.path <- paste(root.path, args[1], sep = "")             ### merge pre-set root path and the first argument as absolute path
  }
  ### stop & print error message if the directory doesn't exist
  if(!dir.exists(full.path)) {
    stop(paste(full.path, "directory doesn't exist!
    If you enter relative path, the root path is set as:", root.path, "
    Relative path cannot start with '/'", sep = " "))
  }
}



setwd(full.path)    ### move to the desired directory
print(paste("Processing files in:", full.path, sep = " "))

fileNames <- Sys.glob("*.fcs")  ### set all files names ending with .fsc in variable

dim <- as.numeric(unlist(strsplit(args[2], split = ",", fixed = T)))

if (length(dim) == 2) {
  ncol = dim[1]
  nrow = dim[2]
} else {
  nfil <- length(fileNames)       ### save number of fcs files
  divfil <- divisors(nfil)        ### find all divisors of number of files

  ### if number of number of files is prime number greater than 2
  if ((nfil > 2) && (length(divfil) == 2)) {
    nfil = nfil + 1               ### increase the number by 1
    divfil <- divisors(nfil)      ### and count divisors again
  }
  ### if there would be to many rows
  if (((divfil[length(divfil) / 2 + 1]) - (divfil[length(divfil) / 2])) > 10) {
    nfil = nfil + 1               ### increase the number by 1
    divfil <- divisors(nfil)      ### and count divisors again
  }
  if ((length(divfil) %% 2) == 0) {         ### when there is even number of divisors
    ncol <- divfil[length(divfil) / 2]      ### set number of columns in final pdf file
    nrow <- divfil[length(divfil) / 2 + 1]  ### set number of rows in final pdf file
  } else {                                  ### in case of odd number of divisors
    ncol <- divfil[length(divfil) / 2 + 1]
    nrow <- ncol
  }
}

pdf(file = "DensPlotAll.pdf", width = (ncol * 5), height = (nrow * 5))
par(mfcol = c(nrow, ncol))        ### densityplot function ignores par function in this case

  for(fcs in fileNames) {         ### loop through all files in the desired directory

    name <- unlist(strsplit(fcs, split = ".", fixed = T))[1]  ### extract filename without file ending
    well <- unlist(strsplit(name, "_"))[3]                ### extract processed well from filename

    data <- read.FCS(fcs, alter.names = T)  ### read FCS file; tranformation argument is "linearize" by default
                                                                              ### to avoid conversion of some measured parameters into "a * 10^(x / R)"
    ### exprs function enables handle flowFrame data as numerics
    fsc.h <- as.vector(exprs(data$FSC.H))
    ssc.h <- as.vector(exprs(data$SSC.H))
    fitc.h <- as.vector(exprs(data$FITC.H))
    time <- as.vector(exprs(data$Time)) / 100             ### the values saved in fcs files are in hundreds of seconds
    fsc.log <- log10(fsc.h)
    ssc.log <- log10(ssc.h)

    ### gating by FCgating.R function
    print(paste("gating", fcs))
    cells <- FCgating(data)

    ### extract fitc values of cells in log10
    fitc.log <- log10(exprs(cells$FITC.H))
    ### replace infinite values by 0s
    fitc.log[is.infinite(fitc.log)] <- 0
    ### count density
    denspl <- density(fitc.log)

    ### count mean value for gated forward scatter
    fsc.m <- mean(log10(exprs(cells$FSC.H)))
    ### count mean value for gated side scatter
    ssc.m <- mean(log10(exprs(cells$SSC.H)))

    ### plotting
    plot(denspl, xlim = c(0.1, 5.5), ylim = c(0.1, 5), main = well,
        xlab = "log10 (FITC.H)", ylab = "Density", col = "red")
    mtext(paste0(" gated: ", length(fitc.log), " ; ",
          format((length(fitc.log) / length(fsc.h) * 100), digits = 4), "%"),
          side = 3, line = -2, adj = 0)
    mtext(paste0(" FITC mean: ", format(mean(fitc.log), digits = 4),
          " ± ", format(sd(fitc.log), digits = 3)),
          side = 3, line = -4, adj = 0)
    mtext(paste0(" FSC mean: ", format(fsc.m, digits = 4)),
          side = 3, line = -6, adj = 0)
    mtext(paste0(" SSC mean: ", format(ssc.m, digits = 4)),
          side = 3, line = -8, adj = 0)
    if ((length(fsc.h) %% 10000) != 0) {
      mtext(paste("Possible overflow", length(fsc.h), "events in",
            format(max(time), digits = 6), "seconds", sep = " "),
            side = 3, line = 0, adj = 0)
    } else {
      mtext(paste("No overflow", sep = " "),
            side = 3, line = 0, adj = 0)
    }


    print(paste(fcs, "processed", sep = " "))

  }

dev.off()
