#!/usr/local/bin/Rscript
### script visualise number of events during time
### example of use: Rscript EventS.R /Users/mvlkova/Documents/RecA_Promoter/FC/180222/plate1/ 12,8
### example of use2: Rscript EventS.R 180222/plate1/ 12,8
### (/Users/mvlkova/Documents/RecA_Promoter/FC/ is folder where the script is located)
### second argument is optional and defines how to partition the final pdf plot (columns,rows)
library("flowCore") ### handling FCS files

### incorporate function finding all divisors
if(!exists("foo", mode="function")) source("Divisors.R")

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
  if (((divfil[length(divfil) / 2 + 1]) - (divfil[length(divfil) / 2])) > 20) {
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

mxt <- 0                                  ### initialize maximum time value
mnt <- 0                                  ### initialize minimum time value

for (fcs in fileNames) {

  data <- read.FCS(fcs, alter.names = T)                ### read FCS file; tranformation argument is "linearize" by default
                                                                            ### to avoid conversion of some measured parameters into "a * 10^(x / R)"
  ### exprs function enables handle flowFrame data as numerics
  time <- as.vector(exprs(data$Time)) / 100   ### the values saved in fcs files are in hundreds of seconds

  ### save maximum time value
  if (mxt < max(time)) {
      mxt <- max(time)
  }
  ### save minimum time value
  if (mnt > min(time)) {
      mnt <- min(time)
  }

}

pdf(file = "EventS.pdf", width = (ncol * 5), height = (nrow * 5))
par(mfcol = c(nrow, ncol))        ### densityplot function ignores par function in this case

  for(fcs in fileNames) {         ### loop through all files in the desired directory

    name <- unlist(strsplit(fcs, split = ".", fixed = T))[1]  ### extract filename without file ending
    well <- unlist(strsplit(name, "_"))[3]                ### extract processed well from filename

    data <- read.FCS(fcs, alter.names = T)                ### read FCS file; tranformation argument is "linearize" by default

    ### exprs function enables handle flowFrame data as numerics
    fsc.h <- as.vector(exprs(data$FSC.H))
    time <- as.vector(exprs(data$Time)) / 100   ### the values saved in fcs files are in hundreds of seconds
    time.r <- round(time, digits = 0) + 1       ### round and add 1 to each value, because first seconds would be 0s

    ### specify number of bins
    bins <- seq(min(time.r), max(time.r), 1)

    ### plottinh
    hist(time.r, main = well, xlim = c(mnt + 1, mxt + 1), ylim = c(0, 20000),
        breaks = bins, xlab = "Time (s)", ylab = "Count", col = "red")

    ### adding info about possible overflows
    if ((length(fsc.h) %% 10000) != 0) {
      mtext(paste("Possible overflow", length(fsc.h), "events in",
            format(max(time), digits = 6), "seconds", sep = " "),
            side = 3, line = 0, adj = 0)
    } else {
      mtext(paste("No overflow", sep = " "),
            side = 3, line = 0, adj = 0)
    }

  }

dev.off()

print(paste0("Max time for ", args[1], "is: ",
      format(mxt, digits = 7), " seconds"))
