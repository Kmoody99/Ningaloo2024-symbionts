# Frontiers 2026 reanalysis: R scripts

Statistical analysis for "Thermal history and the condition of symbionts lost by heat-stressed corals from Ningaloo Reef" (*Frontiers in Ecology and Evolution*, 2026). These scripts replace the analysis in the original MEPS submission.

Put the three scripts in the same folder as the data files already in this repository:

- `20yr_SST_CRW.csv`
- `montiPAMvalues_2.csv`
- `Adults_Heat_Stress_Ningaloo2024_long_Aug2025.csv`

Then, in RStudio, set the working directory to that folder (Session > Set Working Directory > To Source File Location) and Source the scripts in order:

| Script | Paper |
|---|---|
| `01_site_temperatures.R` | Methods (site selection), Figure 1B |
| `02_extreme_assay.R` | Results section 1, Figures 2, 3 and S2, Supplementary Tables S1 and S2 |
| `03_moderate_assay.R` | Results section 2, Figure 4, Supplementary Table S3 |

Packages (install once): `install.packages(c("lme4", "lmerTest", "emmeans", "pbkrtest", "coxme"))`

Results are printed in the Console and saved in a `results/` folder: `results/tables/` holds the CSV files behind Supplementary Tables S1 to S3, each script saves a log of its output, and `sessionInfo` files record the R and package versions.
