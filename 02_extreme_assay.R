################################################################################
# Ningaloo 2024: EXTREME HEAT ASSAY (Montipora spp., 27 -> 30 -> 33 -> 36 C)
# Script 02 of 03  (Results section 1; Figures 2, 3 and S2; Tables S1 and S2)
#
# Input  : data/montiPAMvalues_2.csv
# Output : results/extreme/  per-colony and per-fragment tables (used by the
#          Python figure scripts), leave-one-out results, console log
#          results/tables/    Table_S1*.csv and Table_S2*.csv
# Packages: lme4, lmerTest, emmeans, pbkrtest
#          install once with install.packages(c("lme4","lmerTest","emmeans","pbkrtest"))
#
# Design: each colony has ONE heated and ONE control fragment, each in its own
# 40 L aquarium. Timepoints: t0 27 C, t1 30 C, t2 33 C, t3 36 C. In hospite
# Fv/Fm was measured on the coral fragment (t0 to t3); ex hospite Fv/Fm on
# symbionts filtered from the aquarium water (t0, t2, t3; no filtration at t1).
#
# Analysis choices (see Methods, "Statistical analysis"):
#  * Fv/Fm of a sample = mean of its five focal areas (areas of interest).
#  * Colony BM5 is excluded: necrosis was present at the start of the experiment.
#  * Fragments classed as dead (>50% tissue loss) are not analysed from that
#    timepoint onwards.
#  * Control ex hospite Fv/Fm rose over the experiment, so heated fragments are
#    compared with the SAME-DAY control fragment of the same colony (not with t0).
#  * Linear mixed models: Fv/Fm ~ Treatment * Site + (1 | Colony); F tests with
#    Satterthwaite df; within-site contrasts with emmeans, Kenward-Roger df.
#  * Exact permutation test of the site difference in colony-level
#    (heated minus control) differences: all 126 splits of 9 colonies into 4 + 5.
#  * Leave-one-colony-out refits of the main model.
#  * t3 (36 C): too few survivors, descriptive only.
################################################################################

# Data files are looked for in data/ or in the working directory, so the script
# also works when it sits in the same folder as the CSV files.
if (!dir.exists("data") && dir.exists("../data")) setwd("..")
find_data <- function(f) {
  for (p in c(file.path("data", f), f)) if (file.exists(p)) return(p)
  stop("Cannot find ", f, ". Put it in the same folder as this script (or in data/) ",
       "and set the working directory to that folder (Session > Set Working Directory > To Source File Location).")
}
suppressMessages({ library(lme4); library(lmerTest); library(emmeans) })
emm_options(lmer.df = "kenward-roger")
out_dir <- file.path("results", "extreme"); tab_dir <- file.path("results", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(out_dir, "extreme_assay_log.txt"), split = TRUE)

## ---- 0. Read and clean ----------------------------------------------------
raw <- read.csv(find_data("montiPAMvalues_2.csv"), encoding = "latin1", stringsAsFactors = FALSE)
raw[] <- lapply(raw, function(x) if (is.character(x)) trimws(x) else x)
aoi <- paste0("fvfm", 1:5)
raw[aoi] <- lapply(raw[aoi], function(x) suppressWarnings(as.numeric(x)))
raw$FvFm      <- rowMeans(raw[aoi], na.rm = TRUE)
raw$FvFm[is.nan(raw$FvFm)] <- NA
raw$Timepoint <- factor(c("12.04.24" = "t0", "14.04.24" = "t1", "17.04.24" = "t2",
                          "19.04.24" = "t3")[raw$Date], levels = c("t0", "t1", "t2", "t3"))
raw$Colony    <- toupper(raw$Genotype)
raw$Site      <- factor(raw$Reef, levels = c("Bruboodjoo", "Tantabiddi"))
raw$Treatment <- factor(raw$Treatment, levels = c("control", "heat"), labels = c("Control", "Heat"))
raw$Sample    <- ifelse(raw$Type == "sym", "ex hospite", "in hospite")
raw$Fragment  <- paste(raw$Colony, raw$Treatment)

all_frag <- subset(raw, Colony != "BM5")                        # BM5 excluded
pam      <- subset(all_frag, PAR != "dead" & Alive == "alive")   # live, measured samples
temps    <- c(t0 = "27 C", t1 = "30 C", t2 = "33 C", t3 = "36 C")

## ---- helpers --------------------------------------------------------------
banner <- function(txt) cat("\n\n==================== ", txt, " ====================\n", sep = "")
exact_perm <- function(x, g) {   # exact two-sided permutation test, difference in means
  g <- droplevels(factor(g)); lv <- levels(g)
  obs <- mean(x[g == lv[2]]) - mean(x[g == lv[1]])
  cmb <- combn(length(x), sum(g == lv[1]))
  pr  <- apply(cmb, 2, function(i) mean(x[-i]) - mean(x[i]))
  c(diff = obs, p = mean(abs(pr) >= abs(obs) - 1e-12), n_perm = ncol(cmb))
}
wide <- function(sample_type) {  # one row per colony; columns control_t0, heat_t0, ...
  d <- subset(pam, Sample == sample_type, select = c(Colony, Site, Treatment, Timepoint, FvFm))
  d$key <- paste(tolower(d$Treatment), d$Timepoint, sep = "_")
  w <- reshape(d[, c("Colony", "Site", "key", "FvFm")], idvar = c("Colony", "Site"),
               timevar = "key", direction = "wide")
  names(w) <- sub("FvFm\\.", "", names(w)); w[order(w$Site, w$Colony), ]
}
aov_rows <- function(m, model, response, data_label) {
  a <- anova(m)
  data.frame(Model = model, Response = response, Data = data_label, Term = rownames(a),
             NumDF = a$NumDF, DenDF = a$DenDF, F = a$`F value`, p = a$`Pr(>F)`, row.names = NULL)
}
within_site <- function(m, response, temp) {
  e <- as.data.frame(pairs(emmeans(m, ~ Treatment | Site), reverse = TRUE))
  data.frame(Response = response, Temperature = unname(temp), Site = e$Site, Contrast = "Heated - control",
             Estimate = e$estimate, SE = e$SE, df = e$df, t = e$t.ratio, p = e$p.value,
             row.names = NULL)
}
fit <- function(dat) suppressMessages(lmer(FvFm ~ Treatment * Site + (1 | Colony), data = dat))

## ---- 1. Fragments alive and measured ---------------------------------------
banner("1. FRAGMENTS ALIVE AND MEASURED (in hospite)")
print(ftable(with(subset(pam, Sample == "in hospite"), table(Site, Treatment, Timepoint))))

## ---- 2. Controls over time ---------------------------------------------------
banner("2a. IN HOSPITE controls: site and time")
ci <- subset(pam, Treatment == "Control" & Sample == "in hospite")
print(aggregate(FvFm ~ Site, data = ci, FUN = function(x) round(c(mean = mean(x), SD = sd(x)), 3)))
m_ci <- suppressMessages(lmer(FvFm ~ Site * Timepoint + (1 | Colony), data = ci)); print(anova(m_ci))

banner("2b. EX HOSPITE controls: site and time")
ce <- subset(pam, Treatment == "Control" & Sample == "ex hospite"); ce$Timepoint <- droplevels(ce$Timepoint)
print(aggregate(FvFm ~ Site, data = ce, FUN = function(x) round(c(mean = mean(x), SD = sd(x)), 3)))
print(aggregate(FvFm ~ Site + Timepoint, data = ce, FUN = function(x) round(mean(x), 3)))
m_ce <- suppressMessages(lmer(FvFm ~ Site * Timepoint + (1 | Colony), data = ce)); print(anova(m_ce))
cat("Control ex hospite Fv/Fm as % of control in hospite Fv/Fm (all timepoints):",
    round(100 * mean(ce$FvFm) / mean(ci$FvFm)), "%\n")

## ---- 3. Heat vs same-day control: linear mixed models ------------------------
banner("3. HEAT vs SAME-DAY CONTROL (Treatment * Site + (1 | Colony))")
fits <- list(
  list("in hospite", "t1", fit(subset(pam, Sample == "in hospite" & Timepoint == "t1"))),
  list("in hospite", "t2", fit(subset(pam, Sample == "in hospite" & Timepoint == "t2"))),
  list("ex hospite", "t2", fit(subset(pam, Sample == "ex hospite" & Timepoint == "t2"))))
lmm_tab <- rbind(aov_rows(m_ci, "Controls over time", "in hospite", "Control fragments, 27 C (t0 to t3)"),
                 aov_rows(m_ce, "Controls over time", "ex hospite", "Control aquaria, 27 C (t0, t2, t3)"))
ws_tab <- NULL
for (f in fits) {
  cat("\n---", f[[1]], "at", temps[f[[2]]], "---\n"); print(anova(f[[3]]))
  print(pairs(emmeans(f[[3]], ~ Treatment | Site), reverse = TRUE))
  lmm_tab <- rbind(lmm_tab, aov_rows(f[[3]], "Heated vs same-day control", f[[1]],
                                     paste0("Heated and control fragments, ", temps[f[[2]]], " (", f[[2]], ")")))
  ws_tab  <- rbind(ws_tab, within_site(f[[3]], f[[1]], temps[f[[2]]]))
}
lmm_tab$Term <- sub("Treatment:Site", "Treatment x Site", sub("Site:Timepoint", "Site x Timepoint", lmm_tab$Term))

## ---- 4. Colony-level differences and exact permutation tests -----------------
banner("4. COLONY-LEVEL DIFFERENCES (heated minus same-day control)")
we <- wide("ex hospite"); wi <- wide("in hospite")
for (tp in c("t0", "t2")) we[[paste0("diff_", tp)]] <- we[[paste0("heat_", tp)]] - we[[paste0("control_", tp)]]
for (tp in c("t0", "t1", "t2")) wi[[paste0("diff_", tp)]] <- wi[[paste0("heat_", tp)]] - wi[[paste0("control_", tp)]]
cat("\nEx hospite, 33 C (t2):\n")
print(transform(we[, c("Colony", "Site", "control_t2", "heat_t2", "diff_t2")],
                control_t2 = round(control_t2, 3), heat_t2 = round(heat_t2, 3), diff_t2 = round(diff_t2, 3)),
      row.names = FALSE)
desc <- function(x, s) t(sapply(split(x, s), function(v) round(c(n = length(v), mean = mean(v), SE = sd(v) / sqrt(length(v)), n_higher = sum(v > 0)), 3)))
cat("Mean (SE) of colony differences, and number of colonies where heated > control:\n")
print(desc(we$diff_t2, we$Site))
cat("\nIn hospite, 33 C (t2), % difference heated vs control (mean by site):\n")
print(round(tapply(100 * wi$diff_t2 / wi$control_t2, wi$Site, mean), 1))

perm_tab <- do.call(rbind, lapply(list(
  list("in hospite", "t0", wi$diff_t0, wi$Site), list("ex hospite", "t0", we$diff_t0, we$Site),
  list("ex hospite", "t2", we$diff_t2, we$Site)), function(z) {
    p <- exact_perm(z[[3]], z[[4]])
    data.frame(Response = z[[1]], Temperature = temps[z[[2]]],
               Bruboodjoo_mean_diff = mean(z[[3]][z[[4]] == "Bruboodjoo"]),
               Tantabiddi_mean_diff = mean(z[[3]][z[[4]] == "Tantabiddi"]),
               Site_difference = p["diff"], p = p["p"], n_permutations = p["n_perm"], row.names = NULL)
  }))
cat("\nExact permutation tests (site difference in heated minus control):\n"); print(perm_tab, digits = 3)

## ---- 5. Leave one colony out (main model) ------------------------------------
banner("5. EX HOSPITE 33 C: main model refitted without each colony")
ex2 <- subset(pam, Sample == "ex hospite" & Timepoint == "t2")
loo <- do.call(rbind, lapply(sort(unique(ex2$Colony)), function(cn) {
  a <- anova(fit(subset(ex2, Colony != cn)))
  data.frame(Colony_removed = cn, Site = as.character(ex2$Site[ex2$Colony == cn][1]),
             NumDF = a["Treatment:Site", "NumDF"], DenDF = a["Treatment:Site", "DenDF"],
             F = a["Treatment:Site", "F value"], p = a["Treatment:Site", "Pr(>F)"])
}))
loo <- loo[order(loo$Site, loo$Colony_removed), ]
print(loo, digits = 3, row.names = FALSE)
cat("Interaction p < 0.05 in", sum(loo$p < 0.05), "of", nrow(loo), "refits\n")

## ---- 6. Colour, tissue loss and survival --------------------------------------
banner("6. COLOUR (CoralWatch), TISSUE LOSS AND SURVIVAL")
af <- subset(all_frag, Sample == "in hospite")
af$CoralWatch   <- suppressWarnings(as.numeric(af$CoralWatch))
cw0 <- setNames(af$CoralWatch[af$Timepoint == "t0"], af$Fragment[af$Timepoint == "t0"])
af$CW_change    <- af$CoralWatch - cw0[af$Fragment]
af$Necrosis_pct <- 10 * suppressWarnings(as.numeric(af$Necrosis))
af$Necrosis_pct[af$Alive == "dead" & is.na(af$Necrosis_pct)] <- 100
af$Dead         <- af$Alive == "dead"
cat("Control fragments: colour change from t0 =", paste(unique(na.omit(af$CW_change[af$Treatment == "Control"])), collapse = ", "),
    "; fragments with tissue loss =", sum(af$Necrosis_pct[af$Treatment == "Control"] > 0, na.rm = TRUE), "\n")
b2 <- with(subset(af, Treatment == "Heat" & Timepoint == "t2"), table(Site, Paler = CW_change < 0))
d3 <- with(subset(af, Treatment == "Heat" & Timepoint == "t3"), table(Site, Dead))
cat("\nHeated fragments paler than at t0, 33 C:\n"); print(b2)
cat("\nHeated fragments dead by 36 C:\n"); print(d3)
fisher_tab <- data.frame(
  Outcome = c("Heated fragments paler than at t0 at 33 C", "Heated fragments dead at 36 C"),
  Bruboodjoo = c(paste(b2["Bruboodjoo", "TRUE"], "of", sum(b2["Bruboodjoo", ])), paste(d3["Bruboodjoo", "TRUE"], "of", sum(d3["Bruboodjoo", ]))),
  Tantabiddi = c(paste(b2["Tantabiddi", "TRUE"], "of", sum(b2["Tantabiddi", ])), paste(d3["Tantabiddi", "TRUE"], "of", sum(d3["Tantabiddi", ]))),
  p_Fisher = c(fisher.test(b2)$p.value, fisher.test(d3)$p.value))
print(fisher_tab, digits = 3)
cat("\nHeated fragments alive at 36 C and their Fv/Fm:\n")
print(subset(pam, Timepoint == "t3" & Treatment == "Heat", select = c(Colony, Site, Sample, FvFm)), row.names = FALSE)

## ---- 7. Save -----------------------------------------------------------------
# Per-fragment time series (used for Figure S2 and Table S2)
fr <- merge(af[, c("Colony", "Site", "Treatment", "Timepoint", "CoralWatch", "CW_change", "Necrosis_pct", "Dead")],
            setNames(subset(all_frag, Sample == "in hospite")[, c("Colony", "Treatment", "Timepoint", "FvFm")],
                     c("Colony", "Treatment", "Timepoint", "FvFm_in_hospite")), all.x = TRUE)
fr <- merge(fr, setNames(subset(all_frag, Sample == "ex hospite")[, c("Colony", "Treatment", "Timepoint", "FvFm", "Alive")],
                         c("Colony", "Treatment", "Timepoint", "FvFm_ex_hospite", "ex_alive")), all.x = TRUE)
fr$FvFm_in_hospite[fr$Dead] <- NA                  # not analysed once dead
fr$FvFm_ex_hospite[fr$Dead | fr$ex_alive %in% "dead"] <- NA
fr$ex_alive <- NULL
fr <- fr[order(fr$Site, fr$Colony, fr$Treatment, fr$Timepoint), ]
write.csv(fr, file.path(out_dir, "fragment_timeseries.csv"), row.names = FALSE)
write.csv(we, file.path(out_dir, "exhospite_per_colony.csv"), row.names = FALSE)
write.csv(wi, file.path(out_dir, "inhospite_per_colony.csv"), row.names = FALSE)
write.csv(loo, file.path(out_dir, "exhospite_leave_one_out.csv"), row.names = FALSE)

write.csv(lmm_tab,    file.path(tab_dir, "Table_S1A_extreme_linear_mixed_models.csv"), row.names = FALSE)
write.csv(ws_tab,     file.path(tab_dir, "Table_S1B_extreme_within_site_contrasts.csv"), row.names = FALSE)
write.csv(perm_tab,   file.path(tab_dir, "Table_S1C_extreme_permutation_tests.csv"), row.names = FALSE)
write.csv(loo,        file.path(tab_dir, "Table_S1D_extreme_leave_one_out.csv"), row.names = FALSE)
write.csv(fisher_tab, file.path(tab_dir, "Table_S1E_extreme_bleaching_mortality.csv"), row.names = FALSE)
s2 <- reshape(fr[, c("Colony", "Site", "Treatment", "Timepoint", "FvFm_in_hospite", "FvFm_ex_hospite", "Dead")],
              idvar = c("Colony", "Site", "Treatment"), timevar = "Timepoint", direction = "wide")
write.csv(s2[order(s2$Site, s2$Colony, s2$Treatment), ], file.path(tab_dir, "Table_S2_extreme_per_fragment_FvFm.csv"), row.names = FALSE)

sink()
capture.output(sessionInfo(), file = file.path("results", "sessionInfo_extreme.txt"))
cat("\nDone. Results in results/extreme and results/tables\n")
