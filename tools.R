library(mgcv)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(gratia)
library(tidygam)
library(gridExtra)
library('tidyverse')
#library('RColorBrewer') # not needed anymore
library('flowCore')
library('flowDensity')
library('flowClust')
library('cyanoFilter')
library('PeacoQC')
library('lubridate')
library('EMCluster')

# functions ====
source('flowCore_get_datasets.R') #getdatasets function

#Read and clean an FCS file: remove margin events, remove doublets, bad data, NAs, and log-transform
#n = the well index within the FCS file (e.g. if the file contains 5 observations/measured wells, the valid options for "dataset" are 1, 2, 3, 4, or 5)
#channelALin = the "area" channel supplied in the flow cytometer settings - used to compare against the "height" channel to distinguish between singlets and doublets. See the arguments for PeacoQC::RemoveDoublets
#channelHLin = the "height" channel corresponding to the "area" channel. See the arguments for PeacoQC::RemoveDoublets
clean_data <- function(filename, n=n, channelALin='SSC.ALin', 
                       channelHLin = 'SSC.HLin',...) {
  read.FCS(filename = filename, dataset = n, 
           emptyValue = FALSE, alter.names = TRUE,
           truncate_max_range = FALSE) %>% 
    # I SUGGEST IGNORING THE SIDE SCATTER WIDTH (10TH CHANNEL) SINCE IT IS OVERLY ZEALOUS
    # PROBABLY SAFER TO ACTUALLY CREATE A NAMED VECTOR BUT THIS WORKS
    RemoveMargins(channels = 1:9) %>% # remove margin events and doublets 
    RemoveDoublets(channel1 = channelALin, channel2 = channelHLin, nmad = 2) %>%
    noNeg %>% noNA %>% lnTrans %>% # remove bad data and log-transform
    exprs
}

#Read the FCS file to get metadata, 
# do a bunch of cleaning steps with clean_data, and calculate cell density
# folder = folder with the fcs files
# id = a tibble linking the wells to the treatments and the biota they contain
#Output: a tibble with
# date.time: precise time of sampling of cyto 
# filename: where the data come from
# n: well number; not used but can come in handy
# well: well-id
# dil: dilution of the sample
# vol: volume of the sample
# density.uncleaned: density prior to cleaning of the data (e.g. removal of margin events and debris)
# strains: strains that are in the well
# nr_of_strains: nr of strains in the well
# treat: treatment code
# repl: replicate
# date: date, starting from day 0
# hour_in_seconds: precise hour of sampling of cytometer
# data_cleaned: matrix where every line is a cell and every column is a measurement on the cell
# cells: the number of cells (i.e. the number of lines in data_cleaned)
# density: cells/volume * dilution
read_and_clean_data <- function(folder, id) { 
  tibble(fcs.name=list.files(folder, recursive = TRUE, full.names = TRUE, pattern = '\\.fcs')) %>% #all fcs files within a directory
    mutate(data = pmap(.,getdatasets)) %>% 
    unnest(data) %>% 
    left_join(id, by = 'well') %>% # add the info from the mono.id table
    dplyr::filter(!is.na(treat)) %>% #kick out wells that are not in id
    dplyr::select(-all_of(c("fcs.name","id"))) %>%
    na.omit %>% # remove any bad rows
    distinct %>% # omit any exact copies
    mutate(date=date(dmy_hms(date.time))) %>% #split date from time
    mutate(hour_in_seconds=as.numeric(hour(dmy_hms(date.time))*3600 + #h in seconds on a given day 
                                        minute(dmy_hms(date.time))*60 + 
                                        seconds(dmy_hms(date.time)))) %>%
    group_by(strains, treat, repl) %>%
    mutate(date = difftime(date, min(date), units = 'days')) %>% # count from first in days 
    ungroup() %>%
    mutate(data_cleaned = pmap(.,clean_data)) %>% #clean the data
    mutate(cells = pmap(., function(data_cleaned, ...){nrow(data_cleaned)})) %>% # count cells
    unnest(cells) %>%
    mutate(density = cells / vol * dil) %>%   # calculate cell concentration
    rowwise() %>% #get nr of strains per observation
    mutate(nr_of_strains = length(strains), .after=strains) %>%
    ungroup() 
}

#Collapse fcs data of n different single strains into a single fcs data set
# Only works in a dataframe where single strains have been measured
# Done for selected_date and treatments
# Also replicates within dates and treatments are merged
#Output:
# treat: treatment
# date: date
# data_cleaned: the matrices of the two strains stacked together
collapse_fcs <- function(data, strains_per_combo=2, selected_dates=c(0)){
  data_single <- data %>% #just keep single strain rows, and selected date
    dplyr::filter(nr_of_strains==1, date %in% selected_dates)
  combos <- combn(unique(unlist(data_single$strains)), strains_per_combo, simplify=F)   # all possible combos of strains_per_combo strains
  merged_data <- data_single %>%
    rowwise()%>%
    mutate(strains=unlist(strains)) %>% #can be done as there's only a single entry per row now
    mutate(data_cleaned=list(as_tibble(cbind(data_cleaned, strains)))) %>% #add strain info to fcs data
    group_by(date, treat) %>%
    summarise(data_cleaned = list(bind_rows(data_cleaned))) %>% #merge all fcs data from all strains 
    ungroup()
  #now make all combos and for every combo add merged fcs data for all strains
  data_single %>%
    dplyr::select(treat, date) %>%
    distinct %>% # omit any exact copies (caused by dropping repl)
    tidyr::expand(nesting(treat, date), combo=combos) %>%
    left_join(merged_data, by=c("treat", "date")) %>%
    mutate(data_cleaned = pmap(., function(data_cleaned, combo,...){
      data_cleaned %>% dplyr::filter(strains %in% combo)})) #only keep strains in the combo
}

# Train the clustering method on artificially combining fcs files from monocultures
# A wrapper around init.EM
# data_cleaned = the tibble containing the flow cyto data,...
#..variables of these data are the various channels, plus some 
#..other stuff, including a label indicating 
#..what strain the observation belongs to. 
# combo = a string containing the strain names.
# label = the variable name in data_cleaned indicating the strain.
train_clustering <- function(data_cleaned, channels, combo, label, ...) {
  init.EM(data_cleaned[,channels], nclass = length(combo), 
          lab = as.numeric(as.factor(unlist(data_cleaned[, label]))), 
          min.n.iter = 1, EMC = .EMC.Rnd, method = 'Rnd.EM')
}

#Do clustering. For the first day, use the cluster made from the ...
#... collapsed fcs files as initial estimate. 
#For the following days, use the cluster from ...
#... previous day as initial estimate. 
#data = tibble containing the following variables:
# init_cluster_result_i (with i taking the value of the date): initial estimate from collapsed fcs files
# data_cleaned_i (with i taking the value of the date): fcs data for that date
#dates = the dates for which data were collected.
#Add for every day and treatment and species combo the cluster result under "cluster_result_i" with i = the date  
do_clustering <- function(data, dates){
  for (i in c(1:length(dates))){
    initial_guess <- ifelse(i==1, paste("init_cluster_result_", dates[i], sep=""),
                            paste("cluster_result_", dates[i-1], sep=""))
    actual_data <- paste("data_cleaned_", dates[i], sep="")
    cluster <- paste("cluster_result_", dates[i], sep="")
    data <- data %>%
      mutate(actual_data := !!ensym(actual_data)) %>%
      mutate(initial_guess := !!ensym(initial_guess)) %>%
      mutate({{cluster}}:=pmap(., function(actual_data, initial_guess,...){
        emcluster(as_tibble(actual_data) %>% dplyr::select(all_of(channels)), 
                  emobj=initial_guess)}))
  }
  data %>%
    rename_with(~ gsub("cluster_result_", "day ", .x, fixed = TRUE)) %>%
    dplyr::select(all_of(c("strains", "treat", "repl", paste("day", dates)))) %>%
    pivot_longer(cols=all_of(paste("day", dates)), 
                 names_to = "date", 
                 values_to="cluster_result")
  
}
