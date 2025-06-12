library(tidyverse)
library(mgcv)
library(ggplot2)
library (ggpubr)
library (emmeans)
library (car)
library (rstatix)
library (afex)
library(gridExtra)
library (gt)
library(writexl)
cbPalette <- c("Cell.size" = "#1F77B4","Chlorophyll" ="#228B22","Density"= "#FF7F0E",
               "PCGR" = "#8C564B")
env_palette <- c("olivedrab", "orange", "darkred")
PCA_palette <-  c("C_C" = "#2E8B57","T_T" = "#8B0000",
  "P_P" = "#EFC000", "P_C" = "#66C2A5", "T_C" = "#228B22",
  "C_T" = "#CD5C5C","P_T" = "#B22222",
  "C_P" = "#FFD700", "T_P" = "#DAA520")

#1. Loading data====
dat <- readRDS("output.RDS") |>
  dplyr::select(strains, Habitat, Phenotype_Pigment, treat, setup, density,fsc,
                ssc, redb, redr, yelb,
                acclimation, date.num, repl) |>
  dplyr::filter(setup == "DN") |>
  rename(strain = Habitat)|>
  mutate(strain = as_factor(strain)) |>
  group_by(strain, acclimation, treat, repl) |>
  nest() |>
  ungroup() |>
#2. Fit the gams====
#gam for density
  mutate(mod1 = map(data, ~ gam(log10(density) ~ s(date.num, k =7),
                                method="REML",
                                data = .x)))|>
  mutate(data = map2(data, mod1, ~ .x |> 
                       mutate(density_p = (predict)(.y, se.fit = TRUE)$fit)))|>
  dplyr::select(!mod1) |>
#gam for chlorophyll
  mutate(
    mod2 = map(data, ~ gam(log10(redb) ~ s(date.num, k =7),
                           method="REML",
                           data = .x)))|>
  mutate(data = map2(data, mod2, ~ .x |> 
                       mutate(chl_p = 10^predict(.y, se.fit = TRUE)$fit)))|>
  dplyr::select(!mod2) |>
#gam for fsc
mutate(mod3 = map(data, ~ gam(log10(fsc) ~ s(date.num, k = 5),
                              method="REML", data = .x)))|>
mutate(data = map2(data, mod3, ~ .x |> 
                       mutate(fsc_p = 10^predict(.y, se.fit = TRUE)$fit,
                              fsc_p_se = 10^(predict(.y, se.fit = TRUE)$se.fit))))|>
  dplyr::select(!mod3)|>
  unnest(data) |>
  group_by(strain, acclimation, treat, repl)|>
  arrange(strain, acclimation, treat, repl, date.num)|>
#compute pcgr with predicted density values
  mutate(pcgr_p = log10(lead(10^density_p,1)/10^density_p)/lead(date.num,1)-date.num)
  #dplyr::filter(!is.na(pcgr_p))
dat <- dat |> unite("sequence", acclimation, treat, sep = "_")
dat$pcgr_p <- as.numeric(dat$pcgr_p)
##checking gam diagnostics===
gam.check(mod1)
#figure1====
figure <- expand_grid(x = c(2:10), y = c(2:6))
ggplot(figure) + 
  geom_point(aes(x = x, y = y)) + 
  geom_hline(yintercept = 0, color = "gray40", linetype = "dotted") +
  geom_vline(xintercept = 0, color = "gray40", linetype = "dotted") +
  annotate("text", x = 0.5, y = 0.8, label = "ZONE II (Amplified)", size = 4, fontface = "italic") +
  annotate("text", x = -0.5, y = 0.8, label = "ZONE I (Overcompensated)", size = 4, fontface = "italic") +
  annotate("text", x = -0.5, y = -0.5, label = "ZONE IV (Depressed)", size = 4, fontface = "italic") +
  annotate("text", x = 0.5, y = -0.5, label = "ZONE III (Constrained)", size = 4, fontface = "italic") +
  annotate("text", x = 0, y = 0, label = "No Effect", size = 4, fontface = "italic") +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(
    title = "Mechanisms of legacy effects based on sensitivity",
    x = "Sensitivity to X (XX - YY)",
    y = "Legacy of X in Y (XY - YY)") +
  theme_pubr() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

#ANALYSIS====
#CHECK THE RAW DATA
#Get the mean_values per replicate, this is temporal means==
response_values <- dat |> group_by(strain, sequence, repl) |>
  summarise(
   Density = max(density_p, na.rm = TRUE),
   PCGR = max(pcgr_p, na.rm = TRUE),
   Chlorophyll = mean(chl_p, na.rm = TRUE),
  Cell.size = mean(fsc_p, na.rm = TRUE))

##Background stuff,normality of residuals===
lm <- lm(PCGR ~ sequence*strain, data = response_values)
resid <- lm$residuals
qqPlot(resid)
hist(resid)
shapiro.test(resid)
leveneTest(lm)
outlierTest(lm)
#clean the data===
# Remove rows for the specified indices for each variable
# Step 1: Remove specific rows based on indices
response_values_clean <- response_values[-c(159), ]
response_values_clean <- response_values[-c(86, 83, 107, 131), ]
response_values_clean <- response_values[-c(89, 90, 144, 85, 158), ]
response_values_clean <- response_values[-c(83, 77, 80, 72, 71, 98, 104), ]

# Step 2: Replace NA values in specified columns with column means
response_values_clean <- response_values_clean |>
mutate(across(c(Chlorophyll, Density, Cell.size, PCGR), ~ replace_na(., mean(., na.rm = TRUE))))
# Checked data and repeated tests===
#ok for all variables

#5. Chronic sequence analysis====
#filter the data and set reference level===
chronic <-  response_values_clean |>
  ungroup()|> 
  filter(sequence %in% c("C_C", "P_P", "T_T"))|>
  mutate(subject = paste(strain, repl, sep = "_"))|>
  rename(exposure = sequence) |>
  mutate(exposure = relevel(exposure, "C_C", "P_P", "T_T"))
##chronic_results===
## aov_ez is specific for factorial designs, results include effect size of each predictor
#and significantly different pairs. Tukey's correction is applied to reduce error in p value discovery.
chronic_results <- tibble(response = c("Density", "Cell.size", "PCGR", "Chlorophyll")) |>
  mutate(model = map(response, ~ aov_ez(
    id = "subject",
    dv = .x,
    data = chronic,
    between = "strain",
    within= "exposure",
    type = 3)),
  anova_table = map(model, nice),
  emmeans = map(model, ~ emmeans(.x, ~ exposure | strain, adjust = "tukey", level = 0.95)),
  pairs = map(emmeans, ~ contrast(.x, method = list(
    "P_P - C_C" = c(-1, 1, 0),
    "T_T - C_C" = c(-1, 0, 1)))))

anova_chronic <- chronic_results |> select(response, anova_table)|>
  unnest(anova_table)|>
  mutate(
    df_num = sub(",.*", "", df) %>% as.numeric(),
    df_den = sub(".*, ", "", df) %>% as.numeric(),
    F = gsub("\\*+", "", F)) %>%
  select(response, Effect, df_num, df_den, MSE, F, ges, p.value)

library(gt)
anova_chronic|>
  gt() |>
  tab_header(
    title = "Mixed-effect ANOVA results",
    subtitle = "Variability due to Chronic exposures compared to control")|>
  cols_label(
    response = "Response Variable",
    Effect = "Effect",
    df_num = "df (num)",
    df_den = "df (den)",
    MSE = "Mean Sq. Error",
    F = "F-value",
    ges = "Generalized Eta Squared",
    p.value = "p-value")|>
  tab_options(
    table.font.size = px(12),
    heading.title.font.size = px(14),
    heading.subtitle.font.size = px(12),
    data_row.padding = px(5))

pairs_chronic <- chronic_results |>
  mutate(pairs_df = map(pairs, ~ as.data.frame(.x))) |>
  select(response, pairs_df) |>
  unnest(pairs_df)|>
  filter(contrast!= "P_P - T_T")|>
  mutate(
    significance = case_when(
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",   
      p.value < 0.10 ~ ".",   
      TRUE ~ "NA")) 

library(writexl)
write_xlsx(pairs_chronic, "pairs_chronic.xlsx")

#6. Sequential Exposure analysis====
sequence_results <- response_values_clean |>
  mutate(subject = paste(strain, repl, sep = "_")) |>
  separate(sequence, into = c("past", "subsequent")) |>
  expand_grid(response = c("Density", "Cell.size", "PCGR", "Chlorophyll")) |>
  group_by(subsequent, response) |>
  nest() |>
  ungroup() |>
  mutate(data = map2(data, subsequent, ~ .x |>
                       mutate(past = relevel(factor(past), ref = .y)))) |>
  mutate(model = map2(data, response, ~ {
    aov_ez("subject", dv = .y, data = .x,
           between = "strain",
            within = "past",
           type = 3)})) |>
  mutate(
    emmeans = map(model, ~ emmeans(.x, ~ past | strain, adjust = "tukey", level = 0.95)),
    pairs = map(emmeans, ~ pairs(.x)),
    anova_table = map(model, nice)) |>
  select(subsequent, response, anova_table, emmeans, pairs)

anova_sequence <- sequence_results |> 
  select(subsequent, response, anova_table) |>
  unnest(anova_table)|>
  mutate(
    df_num = sub(",.*", "", df) %>% as.numeric(),
    df_den = sub(".*, ", "", df) %>% as.numeric(),
    F = gsub("\\*+", "", F)) |>
  select(subsequent, response, Effect, df_num, df_den, MSE, F, ges, p.value)

anova_sequence |>
  gt() |>
  tab_header(
    title = "Mixed-effect ANOVA Results",
    subtitle = "Sequential exposures compared to chronic counterpart")|>
  cols_label(
    response = "Response Variable",
    subsequent = "Subsequent Environmment",
    Effect = "Effect",
    df_num = "df (num)",
    df_den = "df (den)",
    MSE = "Mean Sq. Error",
    F = "F-value",
    ges = "Generalized Eta Squared",
    p.value = "p-value")|>
  tab_options(
    table.font.size = px(12),
    heading.title.font.size = px(14),
    heading.subtitle.font.size = px(12),
    data_row.padding = px(5))


write_xlsx(anova_sequence, "anova_sequence.xlsx")

pairs_sequence <- sequence_results |>
  mutate(pairs_df = map(pairs, ~ as.data.frame(.x))) |>
  select(response, pairs_df) |>
  unnest(pairs_df) |>
  mutate(
    sequence = str_remove(contrast, " - C_C"),
    sequence = str_remove_all(sequence, " - "),
    significance = case_when(
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      p.value < 0.10 ~ ".",
      TRUE ~ "NA"))

#save in excel===
library(writexl)
write_xlsx(pairs_sequence, "pairs_sequence.xlsx")

#Mean values for plotting, mean of replicates====
sensitivity <- response_values_clean |> 
  group_by(strain, sequence) |> 
  summarise(
    Density_m = mean(Density),
    Density_sd = sd(Density),
    PCGR_m = mean(PCGR),
    PCGR_sd = sd(PCGR),
    Chlorophyll_m = mean(Chlorophyll),
    Chlorophyll_sd = sd(Chlorophyll),
    Cell.size_m = mean(Cell.size),
    Cell.size_sd = sd(Cell.size),
    .groups = "drop") |> 
  mutate(across(ends_with("_sd"), ~ ifelse(is.na(.), 0, .)))|>
  pivot_longer(
    cols = ends_with ("_m") | ends_with("_sd"),
    names_to = c("response", "metric"),
    names_sep = "_",
    values_to = "value")|>
  pivot_wider(names_from = c(sequence),
              values_from = value)

sensitivity_diff <- sensitivity|> dplyr::filter(metric=="m")|>
  group_by(strain, response, metric)|>
  summarise(
   PP = P_P - C_C,
   PC = P_C - C_C,
   TT = T_T - C_C,
   TC = T_C - C_C,
   TP1 = T_T - P_P,
   TP = T_P - P_P,
   PT1 = P_P - T_T,
   PT = P_T - T_T,
   CT1 = C_C - T_T,
   CT = C_T - T_T,
   CP1 = C_C - P_P,
   CP = C_P - P_P) |>
  pivot_longer(cols = 4:15, names_to = "sequence", values_to = "value")

sensitivity_se <- sensitivity|> dplyr::filter(metric=="sd")|>
  group_by(strain, response, metric)|>
  summarise(
    PP = sqrt((P_P^2/3)+ (C_C^2)/3),
    PC = sqrt((P_C^2/3) + (C_C^2/3)),
    TT = sqrt((T_T^2/3) + (C_C^2/3)),
    TC = sqrt((C_T^2/3) + (T_T^2/3)),
    TP1 = sqrt((T_T^2/3) + (P_P^2/3)),
    TP = sqrt((T_P^2/3) + (P_P^2/3)),
    PT1 = sqrt((P_P^2/3) + (T_T^2/3)),
    PT = sqrt((P_T^2/3) + (T_T^2/3)),
    CT1= sqrt((C_C^2/3) + (T_T^2/3)),
    CT = sqrt((C_T^2/3) + (T_T^2/3)),
    CP1 = sqrt((C_C^2/3) + (P_P^2/3)),
    CP = sqrt((C_P^2/3) + (P_P^2/3)))|>
  pivot_longer(cols = 4:15, names_to = "sequence", values_to = "se_value")

sensitivity_final <- merge(sensitivity_diff, sensitivity_se, 
                                by= c("strain", "response", "sequence")) |>
  dplyr::select(-c(metric.x, metric.y))

##PREPARE DATA FOR PLOTTING====
figure_3 <- ggplot(sensitivity_final|> dplyr::filter(sequence%in%c ("PP", "TT"))) +
  geom_point(aes(x = as.factor(sequence), y = value, shape = strain, 
                 colour = response),position = position_dodge(0.5),size = 3)+
  scale_color_manual (values = cbPalette)+
  geom_errorbar(
    aes(
      x = as.factor(sequence),
      ymin = value - (se_value*1.96),
      ymax = value + (se_value*1.96),
      group = strain),
    width = 0.2,
    colour = "grey40",
    position = position_dodge(0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black")+
  facet_wrap(~ response, nrow = 2, scales = "free_y",
             labeller = as_labeller(c(
               "Density" = "maximum Density",
               "PCGR" = "maximum PCGR",
               "Cell.size" = "mean Cell Size",
               "Chlorophyll" = "mean Chlorophyll"
             )))+
  labs(title = "Trait and growth response to chronic exposure",
       x = "Chronic Environment",
       y = "log response ratio")+
  guides(shape = guide_legend(title = "strain"))+
  guides(colour = "none")+
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

##FOR ACUTE RESPONSES ====
#order the groups ===
plot4a <- sensitivity_final |> filter(!sequence %in% c("TT", "PP","CT1", "CP1", "PT1", "TP1"),
                                      response %in% c("Cell.size", "Chlorophyll"))
plot4a <- reorder_levels(plot4a,sequence, order = c("PC", "TC", "CP", "TP", "CT", "PT"))

figure_4a <- ggplot(plot4a) +
  geom_point(aes(x = as.factor(sequence), y = value, pch = strain, colour = response),
             position = position_dodge(0.5),size = 3)+
  scale_color_manual (values = cbPalette)+
  geom_errorbar(
    aes(
      x = as.factor(sequence),
      ymin = value - (se_value*1.96),
      ymax = value + (se_value*1.96),
      group = strain),
    width = 0.2,
    colour = "grey40",
    position = position_dodge(0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black")+
  facet_wrap(~ response, ncol = 2, scales = "free_y",
             labeller = as_labeller(c(
               "Cell.size" = "mean Cell Size",
               "Chlorophyll" = "mean Chlorophyll"
             )))+
  labs(title = "Trait and growth responses during sequential exposures",
       x = "Sequential environment",
       y = "log response ratio")+
  guides(shape = guide_legend(title = "strain"))+
  guides(colour = "none")+
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5))


#for density===
plot4b <- sensitivity_final |> filter(response %in% c("Density","PCGR"),
                                      !sequence %in% c("TT", "PP","CT1", "CP1", "PT1", "TP1"))
plot4b <- reorder_levels(plot4b, sequence, order = c("PC", "TC", "CP", "TP", "CT", "PT"))

figure_4b <- ggplot(plot4b) +
  geom_point(aes(x = as.factor(sequence), y = value, shape = strain, colour = response),
             position = position_dodge(0.5), size = 3)+
  geom_errorbar(aes(x = as.factor(sequence),
                    ymin = value - (se_value * 1.96),
                    ymax = value + (se_value * 1.96),
                    group = strain),
                width = 0.2, colour = "grey40", position = position_dodge(0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  scale_color_manual(values = c("#FF7F0E", "#8C564B")) + 
  facet_wrap(~ response, ncol = 2, scales = "free_y",
             labeller = as_labeller(c(
               "Density" = "maximum Density",
               "PCGR" = "maximum PCGR")))+
  guides(shape = guide_legend(title = "strain", override.aes = list(size = 4)),
         colour = "none") +
  labs(
    x = "Sequential environment",
    y = "log response ratio") +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    strip.text = element_text(size = 11, face = "bold"))

fig4 <- ggarrange(
  figure_4a, figure_4b,
  nrow = 2,
  common.legend = TRUE,
  legend = "bottom" )

#7. Sensitivity plot====
##prepare the data===
final_plot <- sensitivity_final |>
  pivot_wider(
    names_from = "sequence",
    values_from = c(value, se_value)) |>
  rename_with(~ gsub("^value_", "", .x), starts_with("value_")) |>
  rename_with(~ gsub("^se_value_", "se_", .x), starts_with("se_value_"))

build_plot <- function(xval, yval, xerr, yerr) {
  color_scale <- scale_colour_manual(
    values = cbPalette,
    name = "Variable")
  shape_scale <- scale_shape_discrete(name = "strain")
  plot_theme <- theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),                        
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.position = "right",                           
      legend.box = "vertical",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      strip.text = element_text(size = 11, face = "bold"))

  ggplot(final_plot) +
    geom_point(aes(x = {{xval}}, y = {{yval}}, shape = strain, colour = response),
               size = 2.2) +
    geom_errorbar(aes(x = {{xval}}, ymin = {{yval}} - {{yerr}}, ymax = {{yval}} + {{yerr}}, 
                      colour = response),
                  width = 0.01, alpha = 0.6, linewidth = 0.3) +
    geom_errorbarh(aes(y = {{yval}}, xmin = {{xval}} - {{xerr}}, xmax = {{xval}} + {{xerr}}, 
                       colour = response),
                   height = 0.01, alpha = 0.6, linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    color_scale +
    shape_scale +
    plot_theme}

S1 <- build_plot(CT1, CT, se_CT1, se_CT) +
  labs(x = "Sensitivity to control vs warming", y = "Legacy of control in warming")

S2 <- build_plot(CP1, CP, se_CP1, se_CP) +
  labs(x = "Sensitivity to control vs pollution", y = "Legacy of control in pollution")

S3 <- build_plot(PP, PC, se_PP, se_PC) +
  labs(x = "Sensitivity to pollution vs control", y = "Legacy of pollution in control")

S4 <- build_plot(TT, TC, se_TT, se_TC) +
  labs(x = "Sensitivity to warming vs control", y = "Legacy of warming in control")

S5 <- build_plot(PT1, PT, se_PT1, se_PT) +
  labs(x = "Sensitivity to pollution vs warming", y = "Legacy of pollution in warming")

S6 <- build_plot(TP1, TP, se_TP1, se_TP) +
  labs(x = "Sensitvity to warming vs pollution", y = "Legacy of warming in pollution")

figure_5 <- ggarrange(
    S1, S4, S5, S2, S3, S6,
    ncol = 3, nrow = 2,
    common.legend = TRUE,
    legend = "right")
  
ggsave("sensitivity.pdf", plot = figure_5,
       width = 14, height = 10, dpi = 600)
    
#Figure for supplement====
#PCA of response variables====
  # 1. Select only the response variables
  response_vars <- response_values_clean|> group_by(strain, sequence) |>
    select("Density", "PCGR", "Chlorophyll", "Cell.size") |> ungroup()
  
  # 2. Correlation matrix plot
  GGally::ggcorr(response_vars,
                 label = TRUE, 
                 label_round = 2, 
                 label_size = 4, 
                 low = "blue", mid = "white", high = "red", 
                 name = "Correlation")
  
  # 3. Scale only the response variables (not strain/sequence)
  scaled_data <- response_vars |>
    mutate(across(c("Density", "PCGR", "Chlorophyll", "Cell.size"), scale))
  scaled_data <- rename(scaled_data, exposure = sequence)
  
  pca_result <- prcomp(scaled_data |> select("Density", "PCGR", "Chlorophyll", "Cell.size"), 
                       center = TRUE, scale. = TRUE)
  
  # 4. Extract PCA scores
  pca_scores <- as.data.frame(pca_result$x) 
  pca_scores <- pca_scores |> bind_cols(select(scaled_data, strain, exposure))
  loadings <- pca_result$rotation
  loadings <- as.data.frame(loadings)
  loadings$varname <- rownames(loadings)
  
  # 5. Plot PCA====
  #install and load "ggfortify"
  library("ggfortify")
  autoplot(pca_result, 
           data = scaled_data, 
           color = "exposure",
           shape = "strain",
           loadings = TRUE, 
           loadings.label = TRUE)+
    scale_color_manual(values = PCA_palette) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      panel.grid.major = element_line(color = "gray90"))
  
