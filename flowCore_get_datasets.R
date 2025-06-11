library('flowCore')

# fcs.name <- 'C:/Users/mholmes/Documents UNamur/Lab work/03_Experiments/15_2022-06_DIVERCE/2022-06-17/01 Mono/2022-06-17_at_08-34-54pm_10x.fcs'

# this code has been ripped out of the flowCore functions and reworked - the names weren't
# chosen by me, nor the weird structure. I don't know why it is this way but it just
# needs to be

getdatasets <- function(fcs.name) {
  # get well ids
  con <- file(fcs.name, "rb") # connection to file
  offsets <- flowCore:::readFCSheader(con) # get info from file about where data is stored (byte-offsets)
  offsets <- matrix(offsets, nrow = 1, dimnames = list(NULL, names(offsets))) # as matrix
  txt <- flowCore:::readFCStext(con, offsets[1, ], emptyValue = FALSE) # horrible text data
  
  addOff <- 0 # looping variable (I think add offset? like you have an offset and add to it?)
  
  # if there is only one well, use that, else find where it ends
  if("$NEXTDATA" %in% names(txt)) {
    nd <- as.numeric(txt[["$NEXTDATA"]]) } else {nd <- 0} # nd = next data set position
  
  # empty list to store stuff in
  well.ids <- sample.ids <- date.time <- dil <- dens <- vol <- list()
  
  # addOff <- addOff 
  
  # extending data frame to include info about more wells
  offsets <- rbind(offsets, flowCore:::readFCSheader(con, addOff))
  this.txt <- flowCore:::readFCStext(con, 
                                     offsets[nrow(offsets),], 
                                     emptyValue = FALSE)
  # extract information about the wells looked at so far
  well.ids[[1]] <- this.txt[['$WELLID']] 
  sample.ids[[1]] <- this.txt[['GTI$SAMPLEID']]
  date.time[[1]] <- paste0(this.txt[c("$DATE", "$ETIM")], collapse = '_')
  dil[[1]] <- as.numeric(this.txt[['GTI$DILUTIONFACTOR']]) # get dilution
  # compute cell density:
  n.read <- as.numeric(this.txt[['GTI$TOTALNUMREADINGS']]) # get number of readings 
  timestep <- as.numeric(this.txt[['$TIMESTEP']]) # get timesteps
  flowrateconstant <- as.numeric(this.txt[['GTI$FLOWRATECAL']]) # get flow rate constant
  pumpspeed <-as.numeric(this.txt[['GTI$PUMPSAMPLESPEED']]) # get pump speed
  timetaken <- n.read * timestep # calculate time taken
  flowrate <- flowrateconstant * pumpspeed # calculate flow rate
  vol[[1]] <- timetaken * flowrate # calculate volume
  cells <- as.numeric(this.txt[['$TOT']]) # get n cells
  dens[[1]] <- cells / vol[[1]] * dil[[1]]# compute density
  
  i <- 1
  while(nd != 0) { # loop over datasets
    i <- i + 1
    addOff <- addOff + nd
    offsets <- rbind(offsets, flowCore:::readFCSheader(con, addOff))
    this.txt <- flowCore:::readFCStext(con, 
                                       offsets[nrow(offsets),], 
                                       emptyValue = FALSE)
    nd <- as.numeric(this.txt[["$NEXTDATA"]])
    well.ids[[i]] <- this.txt[['$WELLID']]
    sample.ids[[i]] <- this.txt[['GTI$SAMPLEID']]
    date.time[[i]] <- paste0(this.txt[c("$DATE", "$ETIM")], collapse = '_')
    dil[[i]] <- as.numeric(this.txt[['GTI$DILUTIONFACTOR']])
    n.read <- as.numeric(this.txt[['GTI$TOTALNUMREADINGS']]) # get number of readings
    timestep <- as.numeric(this.txt[['$TIMESTEP']]) # get timesteps
    flowrateconstant <- as.numeric(this.txt[['GTI$FLOWRATECAL']]) # get flow rate constant
    pumpspeed <-as.numeric(this.txt[['GTI$PUMPSAMPLESPEED']]) # get pump speed
    timetaken <- n.read * timestep # calculate time taken
    flowrate <- flowrateconstant * pumpspeed # calculate flow rate
    vol[[i]] <- timetaken * flowrate # calculate volume
    cells <- as.numeric(this.txt[['$TOT']]) # get n cells
    dens[[i]] <- cells / vol[[i]] * dil[[i]]# compute density
  }
  
  close(con)
  
  date.time <- gsub(' ', '', date.time)
  # return data frame with all useful info
  return(data.frame('date.time' = unlist(date.time, recursive = FALSE),
                    'filename' = fcs.name, 
                    'n' = 1:length(well.ids), 
                    'well' = unlist(well.ids, recursive = FALSE),
                    'id' = unlist(sample.ids, recursive = FALSE),
                    'dil' = unlist(dil, recursive = FALSE),
                    'vol' = unlist(vol, recursive = FALSE),
                    'density.uncleaned' = unlist(dens, recursive = FALSE)))
}

removedebris <- function(fcs) {
  # select width channel for margin event detection
  if ('SSC.W' %in% colnames(fcs@exprs)) {
    fcs <- fcs %>%
      cellMargin(Channel = 'SSC.W', 
                 type = 'estimate', 
                 y_toplot = 'SSC.HLin') %>%
      reducedFlowframe
  } else {
    fcs <- fcs %>%
      cellMargin(Channel = 'FSC.W', 
                 type = 'estimate', 
                 y_toplot = 'FSC.HLin') %>%
      reducedFlowframe
  } 
}