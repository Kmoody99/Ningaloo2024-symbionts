# Data dictionary

Column definitions for the three data files used by the analysis scripts. Columns marked "not used" are kept as recorded but are not read by the scripts.

## 20yr_SST_CRW.csv

Daily NOAA Coral Reef Watch (CRW) v3.1 5 km data for the pixels at Tantabiddi (21.925 S, 113.975 E) and Bruboodjoo (22.975 S, 113.775 E), downloaded from the NOAA CoastWatch ERDDAP server (dataset NOAA_DHW). Script 01 uses 1 May 2004 to 30 April 2024 and summarises values by calendar year, so 2004 (May to December) and 2024 (January to April) are partial years.

| Column | Description |
|---|---|
| site | Site name (Tantabiddi or Bruboodjoo) |
| time | Date (daily, UTC) |
| latitude, longitude | Pixel coordinates (decimal degrees) |
| CRW_DHW | Degree heating weeks (°C-weeks) |
| CRW_SST | Sea surface temperature (°C) |

## montiPAMvalues_2.csv (extreme heat assay, *Montipora* spp.)

One row per sample: each colony has one heated and one control fragment, each in its own aquarium. Coral fragments (in hospite) were measured at t0 to t3; symbionts filtered from aquarium water (ex hospite) at t0, t2 and t3.

| Column | Description |
|---|---|
| SampleID | Colony, treatment (c = control, h = heat), timepoint and sample type, e.g. bm1h-t2-sym |
| ID | Colony, sample type and date |
| Type | coral = coral fragment (in hospite); sym = filter of symbionts released into the water (ex hospite) |
| Treatment | control (27 °C) or heat |
| Genotype | Colony ID: bm = Bruboodjoo *Montipora*, tm = Tantabiddi *Montipora*, followed by colony number |
| Geno.Treat, Reef.Treat, Reef.Temp | Combined labels (not used) |
| Reef | Site |
| Temp | Aquarium temperature at sampling (°C) |
| Date.Temp | Label of timepoint, treatment and temperature, e.g. t2-H-33 (not used; timepoints are assigned from Date) |
| Date | Sampling date (dd.mm.yy): 12.04.24 = t0 (27 °C), 14.04.24 = t1 (30 °C), 17.04.24 = t2 (33 °C), 19.04.24 = t3 (36 °C) |
| CoralWatch | CoralWatch colour score (1–6) |
| ChangeCoralWatch | Units lost relative to t0; "dead" once the fragment had died |
| necrosis_notes | Field notes |
| Necrosis | Necrosis score, 0–10 (percentage of necrotic tissue in 10% steps); "dead" once the fragment had died |
| Alive | alive or dead (>50% skeleton exposed). Samples marked dead are excluded from analysis, including any measurements taken on the day the fragment was classed as dead |
| PAR | 0 for measured samples; "dead" where the fragment had died (used with Alive to exclude dead samples) |
| F1–F5 | Minimum fluorescence (F0) in each of five focal areas |
| fm1–fm5 | Maximum fluorescence (Fm) in each of five focal areas |
| fvfm1–fvfm5 | Fv/Fm in each of five focal areas; the sample value is the mean of the five. A value of 0 means no detectable variable fluorescence and is kept in the mean |

Colony BM5 is in the file but excluded from analysis (see Methods).

## Adults_Heat_Stress_Ningaloo2024_long_Aug2025.csv (moderate heat assay, *A. tenuis* and *Montipora* spp.)

One row per fragment per day. Each colony had one fragment in each of six aquaria (three control at 27.5 °C, three heated at 31.1 °C).

| Column | Description |
|---|---|
| Treatment | Control or Heat |
| Tank | Aquarium: C1–C3 (control), H1–H3 (heat) |
| Stick, Position | Location of the fragment within the aquarium (not used) |
| Frag_tank_ID, TreatGeno, ID, Frag_ID | Alternative labels (not used) |
| Genotype | Colony ID: BT/TT = Bruboodjoo/Tantabiddi *A. tenuis*; BM/TM = Bruboodjoo/Tantabiddi *Montipora*, followed by colony number |
| FragID_better | Unique fragment ID used in the analysis |
| Species | Acropora_tenuis or Montipora |
| Site | Bruboodjoo or Tantabiddi |
| Site.Treat, Cor.Site.Treat | Combined labels (not used) |
| Index | Sampling index; day relative to the target temperature = Index − 2 (day 0 = 17 April 2024, when heated aquaria reached 31.1 °C). Index 0 and 1 are 13 and 15 April, during the ramp |
| Date | Sampling date (m/d/yy) |
| FvFm | Field notes on the Fv/Fm measurement, e.g. "low" (not used) |
| FvFm1–FvFm3 | Fv/Fm in each of three focal areas |
| FvFmAvg | Mean Fv/Fm of the three focal areas (used in the analysis) |
| TotalArea, NecrosisArea | Fragment area and necrotic area traced in ImageJ (pixels); "geno_exclude" for colonies excluded from analysis |
| Surv_percent, Surv_prop | Proportion of live tissue (not used) |
| Survival | 1 = alive, 0 = dead (>50% skeleton exposed), blank = not recorded |
| Obs | Field notes |
| BleachingScore_update | CoralWatch colour score (1–6) used in the analysis; x = not scored |
| Last column (unnamed) | Notes from checking necrosis and colour scores against photographs (not used) |

## Corrections made in revision

- montiPAMvalues_2.csv: Date.Temp for the 18 heated 36 °C (t3) rows read "t2-H-36" and has been corrected to "t3-H-36". The scripts assign timepoints from Date, so results are unchanged.
- Adults_Heat_Stress_Ningaloo2024_long_Aug2025.csv: a focal-area reading recorded as "0.0.661" (fragment C3B.TT.16, Index 5) has been corrected to 0.661, and FvFmAvg for that row updated from 0.671 to 0.668. Results are unchanged.
