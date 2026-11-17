################################################################################
# Ningaloo 2024: additional analyses for the revised manuscript (Reviewer 1)
# Script 04 of 04
#
# Input  : 20yr_SST_CRW.csv, montiPAMvalues_2.csv,
#          Adults_Heat_Stress_Ningaloo2024_long_Aug2025.csv
# Output : results/revision_checks/revision_checks_log.txt
#          results/tables/Table_S1F_extreme_interaction_contrast.csv
#          results/tables/Table_S3F_moderate_cox_main_effects.csv
#          results/tables/Table_S3G_moderate_losses_before_day0.csv
# Packages: lme4, lmerTest, emmeans, pbkrtest, survival
#
# Reproduces the numbers added in revision:
#  1. Site temperatures: results excluding the partial years 2004 and 2024
#     (Methods, site selection).
#  2. Extreme assay: confidence interval for the treatment x site interaction
#     contrast, within-site confidence intervals and model diagnostics
#     (Methods, Statistical analysis; Results section 1).
#  3. Moderate assay: fragments alive at day 0 (losses before day 0), overall
#     heat and site hazard ratios from main-effects Cox models, and checks of
#     the proportional hazards assumption (Methods, Statistical analysis;
#     Results section 2; Figure 4A, B).
################################################################################

if (!dir.exists("data") && dir.exists("../data")) setwd("..")
find_data <- function(f) {
  for (p in c(file.path("data", f), f)) if (file.exists(p)) return(p)
  stop("Cannot find ", f, ". Put it in the same folder as this script (or in data/) ",
       "and set the working directory to that folder (Session > Set Working Directory > To Source File Location).")
}
suppressMessages({ library(lme4); library(lmerTest); library(emmeans); library(survival) })
emm_options(lmer.df = "kenward-roger")
out_dir <- file.path("results", "revision_checks"); tab_dir <- file.path("results", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(out_dir, "revision_checks_log.txt"), split = TRUE)
banner <- function(txt) cat("\n\n==================== ", txt, " ====================\n", sep = "")

## ---- 1. Site temperatures without the partial years -------------------------
# Same window and annual metrics as script 01 (1 May 2004 to 30 April 2024,
# summarised by calendar year), so 2004 (May to Dec) and 2024 (Jan to Apr)
# are partial years.
banner("1. SITE TEMPERATURES: EXCLUDING PARTIAL YEARS")
sst <- read.csv(find_data("20yr_SST_CRW.csv"))
sst$date <- as.Date(substr(sst$time, 1, 10))
sst <- subset(sst, date >= as.Date("2004-05-01") & date <= as.Date("2024-04-30"))
sst$year <- as.integer(format(sst$date, "%Y")); sst$month <- as.integer(format(sst$date, "%m"))
monthly <- aggregate(CRW_SST ~ site + year + month, data = sst, FUN = mean)
mmm <- aggregate(CRW_SST ~ site + year, data = monthly, FUN = max); names(mmm)[3] <- "MMM"
dhw <- aggregate(CRW_DHW ~ site + year, data = sst, FUN = max); names(dhw)[3] <- "max_DHW"
annual <- merge(mmm, dhw, by = c("site", "year"))
wide <- function(a, v) { w <- reshape(a[, c("site", "year", v)], idvar = "year", timevar = "site", direction = "wide")
  names(w) <- sub(paste0(v, "\\."), "", names(w)); w }
for (ex in list(integer(0), 2024L, c(2004L, 2024L))) {
  a <- subset(annual, !(year %in% ex)); m <- wide(a, "MMM"); d <- wide(a, "max_DHW")
  h <- sum(d$Bruboodjoo > d$Tantabiddi); l <- sum(d$Bruboodjoo < d$Tantabiddi)
  cat(sprintf("Excluding %-10s n = %d years; MMM Tantabiddi %.2f vs Bruboodjoo %.2f, paired t-test p = %.2g; years with DHW > 8: Tantabiddi %d, Bruboodjoo %d; sign test %d of %d, p = %.4f\n",
              if (length(ex)) paste(ex, collapse = ", ") else "none", nrow(m), mean(m$Tantabiddi), mean(m$Bruboodjoo),
              t.test(m$Tantabiddi, m$Bruboodjoo, paired = TRUE)$p.value, sum(d$Tantabiddi > 8), sum(d$Bruboodjoo > 8),
              h, h + l, binom.test(h, h + l)$p.value))
}

## ---- 2. Extreme assay: interaction contrast and diagnostics -----------------
banner("2. EXTREME ASSAY: INTERACTION CONTRAST AND DIAGNOSTICS")
raw <- read.csv(find_data("montiPAMvalues_2.csv"), encoding = "latin1", stringsAsFactors = FALSE)
raw[] <- lapply(raw, function(x) if (is.character(x)) trimws(x) else x)
aoi <- paste0("fvfm", 1:5); raw[aoi] <- lapply(raw[aoi], function(x) suppressWarnings(as.numeric(x)))
raw$FvFm <- rowMeans(raw[aoi], na.rm = TRUE); raw$FvFm[is.nan(raw$FvFm)] <- NA
raw$Timepoint <- factor(c("12.04.24" = "t0", "14.04.24" = "t1", "17.04.24" = "t2", "19.04.24" = "t3")[raw$Date], levels = c("t0", "t1", "t2", "t3"))
raw$Colony <- toupper(raw$Genotype); raw$Site <- factor(raw$Reef, levels = c("Bruboodjoo", "Tantabiddi"))
raw$Treatment <- factor(raw$Treatment, levels = c("control", "heat"), labels = c("Control", "Heat"))
raw$Sample <- ifelse(raw$Type == "sym", "ex hospite", "in hospite")
pam <- subset(raw, Colony != "BM5" & PAR != "dead" & Alive == "alive")
ic_tab <- NULL
for (z in list(c("ex hospite", "t2", "33 C"), c("in hospite", "t2", "33 C"), c("in hospite", "t1", "30 C"))) {
  cat("\n---", z[1], "at", z[3], "---\n")
  m <- suppressMessages(lmer(FvFm ~ Treatment * Site + (1 | Colony), data = subset(pam, Sample == z[1] & Timepoint == z[2])))
  e <- emmeans(m, ~ Treatment | Site)
  ws <- as.data.frame(confint(pairs(e, reverse = TRUE)))
  ic <- as.data.frame(summary(contrast(e, interaction = c(tr = "revpairwise", site = "revpairwise"), by = NULL), infer = c(TRUE, TRUE)))
  print(ws); print(ic)
  vc <- as.data.frame(VarCorr(m)); sw <- shapiro.test(resid(m))
  cat(sprintf("Colony variance = %.3g (singular fit: %s); Shapiro-Wilk on residuals W = %.3f, p = %.3f\n",
              vc$vcov[vc$grp == "Colony"], isSingular(m), sw$statistic, sw$p.value))
  ic_tab <- rbind(ic_tab, data.frame(Response = z[1], Temperature = z[3],
    Contrast = "Heat response, Tantabiddi minus Bruboodjoo", Estimate = ic$estimate, SE = ic$SE, df = ic$df,
    CI_low = ic$lower.CL, CI_high = ic$upper.CL, p = ic$p.value,
    Bruboodjoo_heat_minus_control = ws$estimate[1], Bruboodjoo_CI_low = ws$lower.CL[1], Bruboodjoo_CI_high = ws$upper.CL[1],
    Tantabiddi_heat_minus_control = ws$estimate[2], Tantabiddi_CI_low = ws$lower.CL[2], Tantabiddi_CI_high = ws$upper.CL[2],
    Colony_variance = vc$vcov[vc$grp == "Colony"], Shapiro_W = unname(sw$statistic), Shapiro_p = sw$p.value))
}
write.csv(ic_tab, file.path(tab_dir, "Table_S1F_extreme_interaction_contrast.csv"), row.names = FALSE)

## ---- 3. Moderate assay: losses before day 0, main-effects Cox models --------
banner("3. MODERATE ASSAY: LOSSES BEFORE DAY 0 AND MAIN-EFFECTS COX MODELS")
d <- read.csv(find_data("Adults_Heat_Stress_Ningaloo2024_long_Aug2025.csv"), encoding = "latin1", stringsAsFactors = FALSE, check.names = FALSE)
names(d) <- trimws(names(d)); names(d)[grep("Treatment", names(d))[1]] <- "Treatment"
d <- d[, names(d) != "" & !grepl("^X(\\.|$)|^V[0-9]+$|^Unnamed", names(d)), drop = FALSE]
d[] <- lapply(d, function(x) if (is.character(x)) trimws(x) else x)
d$Day <- as.numeric(d$Index) - 2; d$Alive <- suppressWarnings(as.numeric(d$Survival))
d$Treatment <- factor(d$Treatment, levels = c("Control", "Heat")); d$Site <- factor(d$Site, levels = c("Bruboodjoo", "Tantabiddi"))
d$Species <- factor(ifelse(d$Species == "Montipora", "Montipora spp.", "A. tenuis"), levels = c("A. tenuis", "Montipora spp."))
d$Fragment <- d$FragID_better
d0 <- subset(d, Day == 0)
alive0 <- aggregate(Alive ~ Species + Site + Genotype + Treatment, data = d0, FUN = sum)
keep <- aggregate(Alive ~ Species + Genotype, data = alive0, FUN = min); keep_geno <- keep$Genotype[keep$Alive >= 2]
fr <- unique(d[, c("Species", "Site", "Genotype", "Fragment")])
fr$alive_day0 <- fr$Fragment %in% d0$Fragment[d0$Alive %in% 1]
fr$analysed <- fr$alive_day0 & fr$Genotype %in% keep_geno
loss <- do.call(rbind, lapply(split(fr, list(fr$Species, fr$Site), drop = TRUE), function(g) data.frame(
  Species = g$Species[1], Site = g$Site[1], Fragments_set_up = nrow(g), Fragments_alive_day0 = sum(g$alive_day0),
  Fragments_analysed = sum(g$analysed), Colonies_set_up = length(unique(g$Genotype)),
  Colonies_analysed = length(unique(g$Genotype[g$Genotype %in% keep_geno])))))
rownames(loss) <- NULL; print(loss, row.names = FALSE)
write.csv(loss, file.path(tab_dir, "Table_S3G_moderate_losses_before_day0.csv"), row.names = FALSE)

dd <- subset(d, Fragment %in% fr$Fragment[fr$analysed] & Day >= 0)
surv <- do.call(rbind, lapply(split(dd, dd$Fragment), function(f) {
  f <- f[order(f$Day), ]; f <- f[!is.na(f$Alive), ]; died <- which(f$Alive == 0)
  data.frame(Species = f$Species[1], Site = f$Site[1], Treatment = f$Treatment[1], Genotype = f$Genotype[1],
             time = if (length(died)) f$Day[min(died)] else max(f$Day), event = as.integer(length(died) > 0)) }))
me_tab <- NULL
for (sp in levels(surv$Species)) {
  s <- subset(surv, Species == sp)
  cat("\n---", sp, ": Cox model with treatment and site as main effects (robust SE clustered by colony) ---\n")
  m <- coxph(Surv(time, event) ~ Treatment + Site + cluster(Genotype), data = s)
  ci <- summary(m)$conf.int; co <- summary(m)$coefficients; print(round(ci[, c(1, 3, 4)], 3)); print(round(co, 4))
  me_tab <- rbind(me_tab, data.frame(Species = sp, Model = "Cox, main effects, robust SE clustered by colony",
    Term = c("Heat", "Site (Tantabiddi)"), HR = ci[, 1], CI_low = ci[, 3], CI_high = ci[, 4], p = co[, ncol(co)], row.names = NULL))
  cat("Proportional hazards check (Schoenfeld residuals, interaction model):\n")
  print(cox.zph(coxph(Surv(time, event) ~ Treatment * Site, data = s)))
}
write.csv(me_tab, file.path(tab_dir, "Table_S3F_moderate_cox_main_effects.csv"), row.names = FALSE)

sink()
capture.output(sessionInfo(), file = file.path("results", "sessionInfo_revision_checks.txt"))
cat("\nDone. Results in results/revision_checks and results/tables\n")
