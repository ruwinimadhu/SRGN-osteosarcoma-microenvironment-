# =============================================================================
# SRGN and recurrence-free survival in human osteosarcoma (GSE39055)
#
# Reproduces the human survival analysis (Fig. 3) from:
#   "Serglycin is a conserved osteoclast-enriched marker of the osteosarcoma
#    microenvironment with prognostic relevance across dogs and humans"
#
# Data:    GSE39055 (Kelly et al., 2013), Illumina HumanHT-12 WG-6 v3 (GPL14951)
# Method:  SRGN expression averaged across two probes, optimal-cutpoint
#          stratification, Kaplan-Meier + log-rank, Cox proportional hazards.
# Author:  Ruwini Madhushika Herath
# License: MIT
# =============================================================================

# --- Dependencies ------------------------------------------------------------
# install.packages(c("survival", "survminer"))
# BiocManager::install("GEOquery")
library(GEOquery)
library(survival)
library(survminer)

set.seed(1)

# --- 1. Download series matrix + platform annotation -------------------------
gse   <- getGEO("GSE39055", GSEMatrix = TRUE, getGPL = TRUE)[[1]]
expr  <- exprs(gse)
pheno <- pData(gse)

# --- 2. SRGN expression (mean of the two probes on GPL14951) -----------------
srgn_probes <- c("ILMN_1760347", "ILMN_2169152")
srgn_probes <- srgn_probes[srgn_probes %in% rownames(expr)]
stopifnot(length(srgn_probes) > 0)
srgn <- colMeans(expr[srgn_probes, , drop = FALSE])

# --- 3. Locate the recurrence-free survival time + event columns -------------
# GSE39055's phenotype column names are not stable across GEO revisions, so we
# detect them rather than hardcode. Inspect the candidates printed below; if the
# automatic pick is wrong, set time_col / event_col manually.
candidates <- grep("time|surviv|event|recur|relapse|status|rfs",
                   colnames(pheno), ignore.case = TRUE, value = TRUE)
message("Candidate survival columns:\n  ", paste(candidates, collapse = "\n  "))

time_col  <- grep("time|rfs.*time|relapse.*time|recur.*time",
                  candidates, ignore.case = TRUE, value = TRUE)[1]
event_col <- grep("event|status|recur|relapse",
                  candidates, ignore.case = TRUE, value = TRUE)[1]

# Manual override (uncomment and edit if the auto-pick above is wrong):
# time_col  <- "rfs time:ch1"
# event_col <- "rfs event:ch1"

stopifnot(!is.na(time_col), !is.na(event_col))
message("Using time column:  ", time_col)
message("Using event column: ", event_col)

# --- 4. Assemble the analysis frame ------------------------------------------
df <- data.frame(
  time  = as.numeric(pheno[[time_col]]),
  event = as.numeric(pheno[[event_col]]),
  srgn  = srgn
)
df <- df[complete.cases(df), ]

# Self-check against the manuscript (Results: "n=37 patients, 18 events").
message("n = ", nrow(df), " (manuscript: 37)")
message("events = ", sum(df$event), " (manuscript: 18)")

# --- 5a. Median split  -> Results: "did not reach significance (p=0.72)" ------
# Reported as the null contrast to the optimal cutpoint; coded here so every
# number in the Results paragraph has a script behind it.
df$group_median <- factor(ifelse(df$srgn > median(df$srgn), "high", "low"),
                          levels = c("low", "high"))
fit_median <- survfit(Surv(time, event) ~ group_median, data = df)
print(surv_pvalue(fit_median, data = df))     # manuscript: p = 0.72

# --- 5b. Optimal expression cutpoint  -> Methods; Fig. 3 ----------------------
cut <- surv_cutpoint(df, time = "time", event = "event", variables = "srgn")
print(cut)                                    # reported threshold: SRGN = 10.42
df$group <- surv_categorize(cut)$srgn         # "high" / "low"
df$group <- factor(df$group, levels = c("low", "high"))

# Fig. 3 reports SRGN-low n=27, SRGN-high n=10. If this table disagrees, the
# cutpoint or the input expression has drifted from the published version.
print(table(df$group))                        # manuscript: low = 27, high = 10

# --- 6. Kaplan-Meier + log-rank  -> Fig. 3 (log-rank p = 0.04) ----------------
# UNIT CHECK: Fig. 3 states time is in MONTHS. Confirm the detected time column
# is in months; if it is in days, change xlab and re-check the caption.
fit <- survfit(Surv(time, event) ~ group, data = df)

km_plot <- ggsurvplot(
  fit, data = df,
  pval = TRUE, risk.table = TRUE,
  legend.labs = c("SRGN low", "SRGN high"),
  xlab = "Time (months)",
  ylab = "Recurrence-free survival"
)
print(km_plot)

# --- 7. Cox proportional hazards (HR + 95% CI for the figure caption) --------
cox <- coxph(Surv(time, event) ~ group, data = df)
print(summary(cox))                           # reported: HR = 2.66 (1.01-7.01)

# --- 8. Reproducibility ------------------------------------------------------
sessionInfo()
