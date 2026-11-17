################################################################################
# Ningaloo 2024: MODERATE HEAT ASSAY (A. tenuis and Montipora spp., 31.1 C)
# Script 03 of 03  (Results section 2; Figure 4; Table S3)
#
# Input  : data/Adults_Heat_Stress_Ningaloo2024_long_Aug2025.csv
# Output : results/moderate/  survival and change tables (used by figure4.py),
#          console log
#          results/tables/    Table_S3*.csv
# Packages: survival (installed with R), coxme, lme4, lmerTest, emmeans, pbkrtest
#          install once with
#          install.packages(c("coxme","lme4","lmerTest","emmeans","pbkrtest"))
#
# Design: 3 control aquaria (C1 to C3, 27.5 C) and 3 heated aquaria (H1 to H3,
# 31.1 C). Each colony has one fragment in every aquarium (3 per treatment).
# Day 0 = 17 April 2024, when the heated aquaria reached 31.1 C. Monitoring
# ended at day 5 (A. tenuis) and day 9 (Montipora spp.).
#
# Analysis choices (see Methods, "Statistical analysis"):
#  * Only fragments alive at day 0; colonies kept if at least 2 fragments per
#    treatment were alive at day 0.
#  * Survival: one record per fragment (day of death, or last day of monitoring
#    if it survived). Cox models with Treatment * Site and robust standard
#    errors clustered by colony; heat effects within each site and site effects
#    within each treatment from separate models fitted to each subset; checked
#    with mixed-effects Cox models (colony and aquarium as random intercepts).
#  * Montipora spp. sensitivity analysis: refit without the Tantabiddi colonies
#    in which every control fragment died.
#  * Change in Fv/Fm and CoralWatch score from day 0 to the last day alive:
#    linear mixed models, Treatment * Site + (1 | colony) + (1 | aquarium);
#    within-site contrasts with emmeans (Kenward-Roger df).
################################################################################

# Data files are looked for in data/ or in the working directory, so the script
# also works when it sits in the same folder as the CSV files.
if (!dir.exists("data") && dir.exists("../data")) setwd("..")
find_data <- function(f) {
  for (p in c(file.path("data", f), f)) if (file.exists(p)) return(p)
  stop("Cannot find ", f, ". Put it in the same folder as this script (or in data/) ",
       "and set the working directory to that folder (Session > Set Working Directory > To Source File Location).")
}
suppressMessages({ library(survival); library(coxme); library(lme4); library(lmerTest); library(emmeans) })
emm_options(lmer.df = "kenward-roger")
out_dir <- file.path("results", "moderate"); tab_dir <- file.path("results", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(out_dir, "moderate_assay_log.txt"), split = TRUE)
banner <- function(txt) cat("\n\n==================== ", txt, " ====================\n", sep = "")

## ---- 0. Read and clean ----------------------------------------------------
d <- read.csv(find_data("Adults_Heat_Stress_Ningaloo2024_long_Aug2025.csv"),
              encoding = "latin1", stringsAsFactors = FALSE, check.names = FALSE)
names(d) <- trimws(names(d))
names(d)[grep("Treatment", names(d))[1]] <- "Treatment"   # first column name carries a byte-order mark
d <- d[, names(d) != "" & !grepl("^X(\\.|$)|^V[0-9]+$|^Unnamed", names(d)), drop = FALSE]
d[] <- lapply(d, function(x) if (is.character(x)) trimws(x) else x)
stopifnot(nrow(d) == 2376)
d$Day       <- as.numeric(d$Index) - 2                      # day 0 = 17 April 2024
d$Alive     <- suppressWarnings(as.numeric(d$Survival))     # 1 alive, 0 dead, NA not recorded
d$FvFm      <- suppressWarnings(as.numeric(d$FvFmAvg))      # mean of 3 focal areas
d$Colour    <- suppressWarnings(as.numeric(d$BleachingScore_update))
d$Treatment <- factor(d$Treatment, levels = c("Control", "Heat"))
d$Site      <- factor(d$Site, levels = c("Bruboodjoo", "Tantabiddi"))
d$Species   <- factor(ifelse(d$Species == "Montipora", "Montipora spp.", "A. tenuis"),
                      levels = c("A. tenuis", "Montipora spp."))
d$Fragment  <- d$FragID_better

## ---- 1. Colonies and fragments analysed --------------------------------------
d0 <- subset(d, Day == 0)
alive0 <- aggregate(Alive ~ Species + Site + Genotype + Treatment, data = d0, FUN = sum)
keep <- aggregate(Alive ~ Species + Genotype, data = alive0, FUN = min)
keep_geno <- keep$Genotype[keep$Alive >= 2]
banner("1. COLONIES AND FRAGMENTS ANALYSED")
cat("Colonies excluded (fewer than 2 fragments alive per treatment at day 0, or no data):\n")
print(sort(setdiff(unique(d$Genotype), keep_geno)))
frag_ok <- d0$Fragment[d0$Alive %in% 1 & d0$Genotype %in% keep_geno]
dd <- subset(d, Fragment %in% frag_ok & Day >= 0)
print(aggregate(Genotype ~ Species + Site, data = unique(dd[, c("Species", "Site", "Genotype")]), FUN = length))
print(ftable(with(unique(dd[, c("Species", "Site", "Treatment", "Fragment")]), table(Species, Site, Treatment))))

## ---- 2. Survival ---------------------------------------------------------------
surv_rows <- do.call(rbind, lapply(split(dd, dd$Fragment), function(f) {
  f <- f[order(f$Day), ]; f <- f[!is.na(f$Alive), ]
  died <- which(f$Alive == 0)
  data.frame(Fragment = f$Fragment[1], Species = f$Species[1], Site = f$Site[1],
             Treatment = f$Treatment[1], Genotype = f$Genotype[1], Tank = f$Tank[1],
             time  = if (length(died)) f$Day[min(died)] else max(f$Day),
             event = as.integer(length(died) > 0))
}))
rownames(surv_rows) <- NULL
write.csv(surv_rows, file.path(out_dir, "survival_one_row_per_fragment.csv"), row.names = FALSE)

hr_row <- function(m, term, species, model, label = term) {
  ci <- summary(m)$conf.int; co <- summary(m)$coefficients
  data.frame(Species = species, Model = model, Term = label, HR = ci[term, 1],
             CI_low = ci[term, 3], CI_high = ci[term, 4], p = co[term, ncol(co)], row.names = NULL)
}
coxme_rows <- function(m, species, model) {
  b <- fixef(m); se <- sqrt(diag(as.matrix(vcov(m)))); z <- b / se
  data.frame(Species = species, Model = model, Term = names(b), HR = exp(b),
             CI_low = exp(b - 1.96 * se), CI_high = exp(b + 1.96 * se), p = 2 * pnorm(-abs(z)), row.names = NULL)
}
cox_set <- function(s, species, tag = "") {
  out <- NULL
  m <- coxph(Surv(time, event) ~ Treatment * Site + cluster(Genotype), data = s)
  print(round(summary(m)$conf.int[, c(1, 3, 4)], 3)); print(round(summary(m)$coefficients, 4))
  for (tm in rownames(summary(m)$coefficients))
    out <- rbind(out, hr_row(m, tm, species, paste0("Cox, robust SE clustered by colony", tag)))
  for (st in levels(s$Site)) {
    ms <- coxph(Surv(time, event) ~ Treatment + cluster(Genotype), data = subset(s, Site == st))
    out <- rbind(out, hr_row(ms, "TreatmentHeat", species, paste0("Heat effect within site", tag), paste("Heat vs control,", st)))
  }
  for (tr in levels(s$Treatment)) {
    ms <- coxph(Surv(time, event) ~ Site + cluster(Genotype), data = subset(s, Treatment == tr))
    out <- rbind(out, hr_row(ms, "SiteTantabiddi", species, paste0("Site effect within treatment", tag), paste("Tantabiddi vs Bruboodjoo,", tr)))
  }
  out
}
km_end <- function(s) {
  do.call(rbind, lapply(split(s, list(s$Site, s$Treatment), drop = TRUE), function(g) {
    k <- survfit(Surv(time, event) ~ 1, data = g)
    data.frame(Site = g$Site[1], Treatment = g$Treatment[1], n_fragments = nrow(g),
               n_colonies = length(unique(g$Genotype)), n_alive_end = sum(g$event == 0),
               pct_alive_end = 100 * tail(k$surv, 1), row.names = NULL)
  }))
}

surv_tab <- NULL; cox_tab <- NULL
for (sp in levels(surv_rows$Species)) {
  s <- subset(surv_rows, Species == sp)
  banner(paste0("2. SURVIVAL: ", sp, " (to day ", max(s$time), ")"))
  ke <- km_end(s); print(ke, digits = 3, row.names = FALSE)
  surv_tab <- rbind(surv_tab, cbind(Species = sp, ke))
  hk <- survfit(Surv(time, event) ~ 1, data = subset(s, Treatment == "Heat"))
  cat("Heated fragments, both sites pooled: % alive by day\n")
  print(round(setNames(100 * summary(hk, times = 0:max(s$time))$surv, 0:max(s$time))))
  cox_tab <- rbind(cox_tab, cox_set(s, sp))
  mm <- coxme(Surv(time, event) ~ Treatment * Site + (1 | Genotype) + (1 | Tank), data = s)
  cat("\nMixed-effects Cox model (colony and aquarium random intercepts):\n"); print(mm)
  cox_tab <- rbind(cox_tab, coxme_rows(mm, sp, "Mixed-effects Cox (colony, aquarium)"))
  cat("Colonies by number of control fragments that died (Tantabiddi):\n")
  cc <- aggregate(event ~ Genotype, data = subset(s, Site == "Tantabiddi" & Treatment == "Control"),
                  FUN = function(x) c(died = sum(x), n = length(x)))
  print(do.call(data.frame, cc), row.names = FALSE)
}

banner("2c. Montipora spp. SENSITIVITY: without colonies in which all controls died")
mo <- subset(surv_rows, Species == "Montipora spp.")
allc <- aggregate(event ~ Genotype + Site, data = subset(mo, Treatment == "Control"), FUN = min)
drop_col <- allc$Genotype[allc$event == 1]
cat("Colonies removed:", drop_col, "\n")
ms <- subset(mo, !(Genotype %in% drop_col))
ke <- km_end(ms); print(ke, digits = 3, row.names = FALSE)
surv_tab <- rbind(surv_tab, cbind(Species = paste("Montipora spp., without", paste(drop_col, collapse = ", ")), ke))
cox_tab <- rbind(cox_tab, cox_set(ms, "Montipora spp.", paste0(" (without ", paste(drop_col, collapse = ", "), ")")))
term_names <- c("TreatmentHeat" = "Heat", "SiteTantabiddi" = "Site (Tantabiddi)",
                "TreatmentHeat:SiteTantabiddi" = "Heat x Site")
cox_tab$Term <- ifelse(cox_tab$Term %in% names(term_names), term_names[cox_tab$Term], cox_tab$Term)
print(cox_tab, digits = 3, row.names = FALSE)

## ---- 3. Fv/Fm and colour: change from day 0 to last day alive -----------------
first_last <- do.call(rbind, lapply(split(dd, dd$Fragment), function(f) {
  f <- f[order(f$Day), ]
  a <- f[f$Alive %in% 1 & !is.na(f$FvFm), ]
  if (nrow(a) < 2 || a$Day[1] != 0) return(NULL)        # need day 0 and a later live value
  ac <- f[f$Alive %in% 1 & !is.na(f$Colour), ]
  data.frame(Fragment = f$Fragment[1], Species = f$Species[1], Site = f$Site[1],
             Treatment = f$Treatment[1], Genotype = f$Genotype[1], Tank = f$Tank[1],
             last_day = max(a$Day), FvFm_day0 = a$FvFm[1], FvFm_last = a$FvFm[nrow(a)],
             dFvFm = a$FvFm[nrow(a)] - a$FvFm[1],
             dColour = if (nrow(ac) >= 2 && ac$Day[1] == 0) ac$Colour[nrow(ac)] - ac$Colour[1] else NA)
}))
rownames(first_last) <- NULL
write.csv(first_last, file.path(out_dir, "change_per_fragment.csv"), row.names = FALSE)

lmm_tab <- NULL; ws_tab <- NULL; desc_tab <- NULL
fit_change <- function(dat, resp, species, label) {
  dat <- dat[!is.na(dat[[resp]]), ]
  cat("\n---", species, ":", label, "---\n")
  ds <- do.call(rbind, lapply(split(dat, list(dat$Site, dat$Treatment)), function(g)
    data.frame(Species = species, Response = label, Site = g$Site[1], Treatment = g$Treatment[1],
               n = nrow(g), Mean = mean(g[[resp]]), SE = sd(g[[resp]]) / sqrt(nrow(g)))))
  print(ds, digits = 3, row.names = FALSE)
  m <- lmer(as.formula(paste(resp, "~ Treatment * Site + (1 | Genotype) + (1 | Tank)")), data = dat)
  a <- anova(m); print(a)
  e <- as.data.frame(pairs(emmeans(m, ~ Treatment | Site), reverse = TRUE)); print(e)
  lmm_tab <<- rbind(lmm_tab, data.frame(Species = species, Response = label,
    Term = sub("Treatment:Site", "Treatment x Site", rownames(a)), NumDF = a$NumDF, DenDF = a$DenDF,
    F = a$`F value`, p = a$`Pr(>F)`, row.names = NULL))
  ws_tab <<- rbind(ws_tab, data.frame(Species = species, Response = label, Site = e$Site,
    Contrast = "Heated - control", Estimate = e$estimate, SE = e$SE, df = e$df, t = e$t.ratio, p = e$p.value, row.names = NULL))
  desc_tab <<- rbind(desc_tab, ds)
}
for (sp in levels(first_last$Species)) {
  s <- subset(first_last, Species == sp)
  banner(paste0("3. ", sp, ": CHANGE FROM DAY 0 TO LAST DAY ALIVE"))
  fit_change(s, "dFvFm", sp, "Change in Fv/Fm")
  fit_change(s, "dColour", sp, "Change in CoralWatch score")
}
banner("3c. ROBUSTNESS: Montipora spp. Fv/Fm without the most extreme fragment")
mo2 <- subset(first_last, Species == "Montipora spp.")
worst <- mo2$Fragment[which.min(mo2$dFvFm)]
cat("Most extreme fragment:", worst, "change =", round(min(mo2$dFvFm), 3), "\n")
fit_change(subset(mo2, Fragment != worst), "dFvFm", "Montipora spp.",
           paste0("Change in Fv/Fm, without ", worst))

## ---- 4. Save tables ------------------------------------------------------------
write.csv(surv_tab, file.path(tab_dir, "Table_S3A_moderate_survival_summary.csv"), row.names = FALSE)
write.csv(cox_tab,  file.path(tab_dir, "Table_S3B_moderate_cox_models.csv"), row.names = FALSE)
write.csv(lmm_tab,  file.path(tab_dir, "Table_S3C_moderate_linear_mixed_models.csv"), row.names = FALSE)
write.csv(ws_tab,   file.path(tab_dir, "Table_S3D_moderate_within_site_contrasts.csv"), row.names = FALSE)
write.csv(desc_tab, file.path(tab_dir, "Table_S3E_moderate_change_means.csv"), row.names = FALSE)

sink()
capture.output(sessionInfo(), file = file.path("results", "sessionInfo_moderate.txt"))
cat("\nDone. Results in results/moderate and results/tables\n")
