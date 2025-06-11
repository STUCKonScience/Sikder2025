#load required packages====
source("tools.R")

#Define the channel columns to read .FCS data====
channels <- c("FSC.HLin", "SSC.HLin", "GRN.B.HLin", "YEL.B.HLin", 
              "RED.B.HLin", "NIR.B.HLin", "RED.R.HLin", "NIR.R.HLin")

# Read IDs and data====
idJan1 <- read_delim("idJan1.csv")
idJan2 <- read_delim("idJan2.csv")
idOct1 <- read_delim("idOct1.csv")
idOct2 <- read_delim("idOct2.csv")
idNov1 <- read_delim("idNov1.csv")
idNov2 <- read_delim("idNov2.csv")

#Define the pathway for each output file====
outJan1 <- read_and_clean_data(folder= "../DATA/January_24/phase2/plate1",
                               id = idJan1)
outJan2<- read_and_clean_data(folder= "../DATA/January_24/phase2/plate2",
                              id = idJan2)
outOct1 <- read_and_clean_data(folder= "../DATA/October_23/phase2/plate1",
                               id = idOct1)
outOct2<- read_and_clean_data(folder= "../DATA/October_23/phase2/plate2",
                              id = idOct2)
outNov1 <- read_and_clean_data(folder= "../DATA/November_23/phase2/plate1",
                               id = idNov1)
outNov2<- read_and_clean_data(folder= "../DATA/November_23/phase2/plate2",
                              id = idNov2)

#Combine all outputs, define phenotypes and treatment conditions====
out_final <- bind_rows(outJan1, outJan2, outOct1, outOct2, 
                       outNov1, outNov2)|>   
  group_by(date, treat, repl, strains) |>  
##to keep only the last observations for each day since more than one dilution was used===
  dplyr::filter(hour_in_seconds == max(hour_in_seconds)) |>
  separate_wider_delim(cols=treat, delim="_", names=c("treat", "acclimation", "setup")) |> 
  dplyr::filter(setup=="DN")|>
  mutate(Phenotype_Pigment =case_when(
    strains %in% c(2375,2524,556) ~ "2",
    strains %in% c(2383,2555) ~ "1",
    strains==2385 ~ "3a")) %>% ungroup()
#Final output with Scaled trait values and modified columns====
out_final <- out_final %>% 
  mutate(meanTraits = map(data_cleaned, ~ colMeans(.x[,channels]))) %>%
  unnest_wider(meanTraits, names_repair = "unique") %>% 
  mutate(fsc = FSC.HLin,
         ssc = SSC.HLin/FSC.HLin,
         redb = RED.B.HLin/FSC.HLin,
         redr = RED.R.HLin/FSC.HLin,
         yelb = YEL.B.HLin/FSC.HLin,
         acclimation = as.factor(acclimation),
         repl = as.factor(repl),
         date.num = as.numeric(date),
         density = as.numeric(density),
         Habitat = case_when(strains == "2375" ~ "Atlantic WSahara",
                             strains == "2385" ~ "Atlantic Norway",
                             strains == "2555" ~ "Atlantic USA",
                             strains == "2383" ~ "Red Sea1",
                             strains == "556" ~ "Red Sea2",
                             strains == "2524" ~ "Mediterranean Sea"))

#CHECK THE DILUTION FACTORS====
ggplot(out_final) + 
  aes(x=date, y=log10(density), col=as.factor(dil), pch = as.factor(repl)) + 
  geom_point() + 
  facet_grid(cols=vars(treat), rows=vars(strains))+
  theme_bw()

#fix bad dilutions===
out_final <- out_final %>%
  mutate(density = if_else(strains == "2524" & date == 8 & acclimation == "T" & treat == "P",
    density / 10, density))

out_final <- out_final %>%
  mutate(density = if_else(strains == "2524" & date %in% 7:10 & repl %in% c("1", "3") &
      acclimation %in% c("P","T") & treat == "T",density * 10, density))

out_final <- out_final %>% mutate(density = if_else(strains == "2375" & date %in% 6:10 &
      repl=="1" & acclimation %in% c("P") & treat == "T",
    density * 10, density))

out_final <- out_final %>% mutate(density = if_else(strains == "2383" & date=="10" &
                                                      repl=="1" & acclimation %in% c("T") 
                                                    & treat == "T",
                                                density * 10, density))
#test for outliers in data====
ggplot(out_final, aes(x = Habitat, y = fsc)) +
  geom_boxplot(outlier.shape = 16, outlier.colour = "red", outlier.size = 3) +
  labs(title = "Boxplot of Density with Outliers Highlighted", x = "Strain", y = "Density") +
  facet_grid(treat~acclimation)+
  theme_minimal()
#remove a bad fsc point====
out_final <- out_final|> dplyr::filter(fsc > 1)

#Save the final output as .RDS file====
saveRDS(out_final, file= "./output.RDS")
#Check if it's saved correctly===
a <- readRDS("output.RDS")
