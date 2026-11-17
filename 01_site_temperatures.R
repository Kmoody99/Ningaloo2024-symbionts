################################################################################
# Ningaloo 2024: thermal history of the two collection sites
# Script 01 of 03  (Methods: "Selection of sites with distinct temperature
# conditions"; Figure 1B)
#
# Input  : data/20yr_SST_CRW.csv
#          NOAA Coral Reef Watch daily 5 km SST and DHW (ERDDAP dataset NOAA_DHW)
#          for the pixels at Tantabiddi (21.925 S, 113.975 E) and
#          Bruboodjoo (22.975 S, 113.775 E).
# Output : results/site_temperatures/annual_MMM_DHW.csv  (used by figure1.py)
#          results/site_temperatures/site_temperature_stats.txt
# Packages: base R only
#
# Run with the working directory set to the folder that contains the data
# (open the .Rproj file, or Session > Set Working Directory > To Source File Location) and click Source.
################################################################################

# Data files are looked for in data/ or in the working directory, so the script
# also works when it sits in the same folder as the CSV files.
if (!dir.exists("data") && dir.exists("../data")) setwd("..")
find_data <- function(f) {
  for (p in c(file.path("data", f), f)) if (file.exists(p)) return(p)
  stop("Cannot find ", f, ". Put it in the same folder as this script (or in data/) ",
       "and set the working directory to that folder (Session > Set Working Directory > To Source File Location).")
}
out_dir <- file.path("results", "site_temperatures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(out_dir, "site_temperature_stats.txt"), split = TRUE)

## ---- 1. Read data -----------------------------------------------------------
sst <- read.csv(find_data("20yr_SST_CRW.csv"))
sst$date <- as.Date(substr(sst$time, 1, 10))
sst$site <- factor(sst$site, levels = c("Tantabiddi", "Bruboodjoo"))

# 20-year window ending at the April 2024 experiments (1 May 2004 to 30 April 2024).
# As in the original analysis, values are summarised by calendar year within
# this window (21 years, 2004 to 2024).
sst <- subset(sst, date >= as.Date("2004-05-01") & date <= as.Date("2024-04-30"))
sst$year  <- as.integer(format(sst$date, "%Y"))
sst$month <- as.integer(format(sst$date, "%m"))

## ---- 2. Annual metrics per site ---------------------------------------------
# Maximum monthly mean (MMM): highest of the 12 monthly mean SSTs in each year.
monthly <- aggregate(CRW_SST ~ site + year + month, data = sst, FUN = mean)
mmm     <- aggregate(CRW_SST ~ site + year, data = monthly, FUN = max)
names(mmm)[3] <- "MMM"
dhw  <- aggregate(CRW_DHW ~ site + year, data = sst, FUN = max)
names(dhw)[3] <- "max_DHW"
vsd  <- aggregate(CRW_SST ~ site + year, data = sst, FUN = sd)
names(vsd)[3] <- "SST_SD"
vrg  <- aggregate(CRW_SST ~ site + year, data = sst, FUN = function(x) diff(range(x)))
names(vrg)[3] <- "SST_range"
annual <- Reduce(function(a, b) merge(a, b, by = c("site", "year")), list(mmm, dhw, vsd, vrg))
annual <- annual[order(annual$site, annual$year), ]
write.csv(annual, file.path(out_dir, "annual_MMM_DHW.csv"), row.names = FALSE)

wide <- function(v) {
  w <- reshape(annual[, c("site", "year", v)], idvar = "year", timevar = "site", direction = "wide")
  names(w) <- sub(paste0(v, "\\."), "", names(w)); w[order(w$year), ]
}

cat("==================== MAXIMUM MONTHLY MEAN (MMM) SST ====================\n")
m <- wide("MMM")
print(round(sapply(m[, c("Tantabiddi", "Bruboodjoo")], function(x) c(mean = mean(x), SD = sd(x), n = length(x))), 2))
print(t.test(m$Tantabiddi, m$Bruboodjoo, paired = TRUE))

cat("\nAssay temperatures relative to each site's mean MMM (deg C above MMM):\n")
mean_mmm <- colMeans(m[, c("Tantabiddi", "Bruboodjoo")])
print(round(rbind(`31.1 C (moderate assay)` = 31.1 - mean_mmm,
                  `33 C (extreme assay, t2)` = 33 - mean_mmm), 1))

cat("\n==================== ANNUAL MAXIMUM DHW ====================\n")
d <- wide("max_DHW")
cat("Years with DHW > 8:\n"); print(colSums(d[, c("Tantabiddi", "Bruboodjoo")] > 8))
cat("Years exceeding 8 DHW:\n"); print(d[d$Tantabiddi > 8 | d$Bruboodjoo > 8, ], row.names = FALSE)
higher <- sum(d$Bruboodjoo > d$Tantabiddi); lower <- sum(d$Bruboodjoo < d$Tantabiddi)
ties <- sum(d$Bruboodjoo == d$Tantabiddi)
cat(sprintf("\nBruboodjoo higher in %d years, lower in %d, equal in %d (both 0 DHW).\n", higher, lower, ties))
cat("Sign test (ties excluded, as is standard):\n")
print(binom.test(higher, higher + lower, p = 0.5))
cat("For reference, the original MEPS analysis kept the tied years in n (binom.test(",
    higher, ", ", nrow(d), ")): p = ", round(binom.test(higher, nrow(d))$p.value, 3), "\n", sep = "")

cat("\n==================== ANNUAL VARIABILITY ====================\n")
cat("Annual SST range (max minus min, deg C):\n")
print(round(tapply(annual$SST_range, annual$site, mean), 2))
print(summary(aov(SST_range ~ site, data = annual)))
cat("Annual SD of daily SST:\n")
print(round(tapply(annual$SST_SD, annual$site, mean), 2))
print(summary(aov(SST_SD ~ site, data = annual)))

sink()
cat("\nDone. Results in", out_dir, "\n")
