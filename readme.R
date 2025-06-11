readme_text <- '
# README: R Code for Legacy Effects and Trait Plasticity Analysis

## Overview

This R script accompanies the analysis for the manuscript titled **"[
Trait and Growth Responses to Sequential Environmental Change linked to Sensitivity in Synechococcus populations
]"**. 
The analysis investigates the impact of environmental exposure sequences (e.g., control, warming, pollution) on phytoplankton trait plasticity and legacy effects. 
The code includes data preprocessing, modeling with Generalized Additive Models (GAMs), mixed-effect ANOVA tests, sensitivity/legacy metric calculations, and plotting.

---

## Dependencies

The following R packages are used:

```r
tidyverse, mgcv, ggplot2, ggpubr, emmeans, car, rstatix, afex, gridExtra, writexl, GGally, gt, ggfortify
install.packages(c("tidyverse", "mgcv", "ggplot2", "ggpubr", "emmeans", "car",
                   "rstatix", "afex", "gridExtra", "writexl", "GGally", "gt", "ggfortify"))
Workflow Summary
1. Data Loading and Preprocessing

    Loads output.RDS.

    Selects and renames relevant variables.

    Filters for specific experimental setup.

    Nests data by treatment and strain.

2. GAM Modeling

    Fits smoothed time series (GAMs) for:

        Density

        Chlorophyll

        Cell size (FSC)

    Predicts values and calculates PCGR.

3. Chronic and Sequence Exposure Analysis

    Mixed-effect ANOVAs using aov_ez from afex.

    Evaluates chronic exposure sequences (C_C, P_P, T_T).

    Analyzes sequential exposures by comparing pre/post environmental switches.

4. Sensitivity & Legacy Metrics

    Computes log-response ratios between treatment combinations.

    Quantifies:

        Sensitivity to prior conditions

        Legacy effects from past exposures

5. Visualization

    Publication-quality figures:

        Zone diagrams for conceptual legacy types

        Trait and growth response plots

        Sensitivity vs. legacy scatter plots

        PCA of response variables

Outputs

    pairs_chronic.xlsx, pairs_sequence.xlsx: Pairwise statistical comparisons

    anova_sequence.png: ANOVA summary figure

    sensitivity.pdf: Summary of sensitivity and legacy metrics

    PCA plots and correlation matrix for supplementary material

Reproducibility Notes

    Ensure output.RDS is located in the working directory.

    Use R version ≥ 4.1.

    For reproducibility of plots, consider using set.seed() if randomness is involved.

Citation

If you use this code, please cite:

[Arunima Sikder, Colin T Kremer, Frederik De Laender].
"[Trait and Growth Responses to Sequential Environmental Change linked to Sensitivity in Synechococcus populations]." 

Contact

    Author: [Arunima Sikder]

    Email: [arunima.sikder@unamur.be]

    Institution: [University of Namur]
    '