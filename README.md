# Thermal history and the condition of symbionts lost by heat-stressed corals from Ningaloo Reef, Australia

Data and R scripts for Moody et al., *Frontiers in Ecology and Evolution* (manuscript 2043607, in review).

## Files

| File | Contents |
|---|---|
| `20yr_SST_CRW.csv` | Daily NOAA Coral Reef Watch SST and DHW for the two sites |
| `montiPAMvalues_2.csv` | Extreme heat assay (*Montipora* spp.): Fv/Fm of coral fragments and of released symbionts, colour, necrosis and survival |
| `Adults_Heat_Stress_Ningaloo2024_long_Aug2025.csv` | Moderate heat assay (*A. tenuis* and *Montipora* spp.): daily survival, Fv/Fm and colour |
| `DATA_DICTIONARY.md` | Column definitions for each data file, and corrections made in revision |
| `01_site_temperatures.R` | Methods (site selection), Figure 1B |
| `02_extreme_assay.R` | Results section 1, Figures 2, 3 and S2, Supplementary Tables S1 and S2 |
| `03_moderate_assay.R` | Results section 2, Figure 4, Supplementary Table S3 |
| `04_revision_checks.R` | Analyses added in revision: results without the partial years 2004 and 2024, the confidence interval for the extreme-assay interaction and model diagnostics, losses before day 0, overall heat hazard ratios (Figure 4A, B) and proportional hazards checks |

## Running the scripts

Download the repository, open R or RStudio, set the working directory to the repository folder (Session > Set Working Directory > To Source File Location) and source the scripts in order (01 to 04).

Packages (install once): `install.packages(c("lme4", "lmerTest", "emmeans", "pbkrtest", "coxme"))`. The `survival` package comes with R.

Results are printed in the Console and saved in a `results/` folder: `results/tables/` holds the CSV files behind the Supplementary Tables, each script saves a log of its output, and `sessionInfo` files record the R and package versions.

## Notes

- Analyses were run in R 4.6.1. Figures were produced in Python with matplotlib from the tables written by the R scripts.
- These scripts replace the analysis used in an earlier version of the manuscript.
