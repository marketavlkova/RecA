#!/usr/local/bin/Rscript
### gating function for FC data

library("MASS")     ### to run dov.rob method in gating, which gives more oval shaped gates than previously used covMcd

FCgating <- function(data) {

  ### create an elipse gate
  gate.elipse <- norm2Filter("FSC.H", "SSC.H", method = "cov.rob", scale.factor = 1, n = 50000, filterId = "NiceCells")
  ### make a preliminary filter based on the rectangle gate
  filt1 <- filter(data, gate.elipse)
  ### make a preliminary subset from original data
  f1.data <- Subset(data, filt1)
  ### extract values for rectangle gate
  fsc.min <- min(exprs(f1.data$FSC.H))
  fsc.max <- max(exprs(f1.data$FSC.H))
  ssc.min <- min(exprs(f1.data$SSC.H))
  ssc.max <- max(exprs(f1.data$SSC.H))
  if ((log10(fsc.max) - log10(fsc.min)) > 0.5) {
    print("scatter values for second gating corrected")
    fsc.min <- fsc.min + 10000
    fsc.max <- fsc.max + 100000
    ssc.min <- ssc.min + 10000
    ssc.max <- ssc.max + 100000
  }
  ### create a rectangle gate based on first filtered subset
  gate.rectangle <- rectangleGate(filterId = "RemoveMess", "FSC.H" = c((fsc.min - 450), (fsc.max + 3000)), "SSC.H" = c((ssc.min - 500), (ssc.max + 8000)))
  ### make second filter based on values from preliminary filter & filter the original dataset
  filt2 <- filter(data, gate.rectangle)
  f2.data <- Subset(data, filt2)
  ### make final filter based on rectangle-filtered datapoints
  filt <- filter(f2.data, gate.elipse)
  ### make final filtered subset from preliminary one
  cells <- Subset(f2.data, filt)

}
