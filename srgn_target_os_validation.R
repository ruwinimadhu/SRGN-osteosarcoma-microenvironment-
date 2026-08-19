# =============================================================================
#  SRGN prognostic validation in human osteosarcoma (TARGET-OS)
#  Independent-cohort validation for the SRGN / osteosarcoma manuscript.
#
#  Expression : GDC TARGET-OS STAR-FPKM, via UCSC Xena (UCSCXenaTools)
#  Survival   : GDC TARGET-OS harmonised clinical, via TCGAbiolinks
#               (the Xena flat survival file is NOT used — it lacks
#                follow-up time for censored patients)
#  Endpoint   : overall survival (OS). NOTE: the discovery cohort
#               (GSE39055) used recurrence-free survival — a related but
#               NOT identical endpoint. Recurrence times are not recorded
#               in TARGET-OS, so RFS cannot be evaluated here.
#
#  Primary analysis : SRGN as a continuous (z-scored), pre-specified
#                     predictor in Cox regression — no cutpoint search.
#  Author : Ruwini M. Herath
#  License: MIT
# =============================================================================

# ---- 0. Packages ------------------------------------------------------------
# install once:
# install.packages(c("UCSCXenaTools", "survival", "survminer", "dplyr", "tidyr"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("TCGAbiolinks")

suppressPackageStartupMessages({
  library(UCSCXenaTools)
  library(TCGAbiolinks)
  library(survival)
  library(survminer)
  library(dplyr)
  library(tidyr)
})

OUTDIR <- "results"
dir.create(OUTDIR, showWarnings = FALSE)

SRGN_ENSEMBL <- "ENSG00000122862"   # Serglycin
PATIENT_RX   <- "^(TARGET-[0-9]+-[^-]+).*"  # aliquot barcode -> patient stem

# ---- 1. Expression: locate + download STAR-FPKM from the GDC Xena hub --------
os_sets <- XenaData %>% filter(XenaCohorts == "GDC TARGET-OS")

# STAR-FPKM is the current GDC release (HTSeq was retired); fall back if needed.
ds <- os_sets$XenaDatasets
expr_name <- grep("star_fpkm\\.tsv$", ds, value = TRUE, ignore.case = TRUE)
if (length(expr_name) == 0) expr_name <- grep("htseq_fpkm", ds, value = TRUE, ignore.case = TRUE)
if (length(expr_name) == 0) expr_name <- grep("fpkm",       ds, value = TRUE, ignore.case = TRUE)
expr_name <- expr_name[1]
stopifnot(!is.na(expr_name))
message("Expression dataset: ", expr_name)

query <- os_sets %>%
  filter(XenaDatasets == expr_name) %>%
  XenaGenerate() %>%
  XenaQuery()

dl   <- XenaDownload(query, destdir = "xena_target_os", trans_slash = TRUE)
prep <- XenaPrepare(dl)
expr <- if (is.data.frame(prep)) prep else prep[[grep("fpkm", names(prep), ignore.case = TRUE)]]

# ---- 2. Extract SRGN (values are ALREADY log2(fpkm+1) — do NOT re-log) -------
id_col <- names(expr)[1]
srgn <- expr %>%
  filter(grepl(SRGN_ENSEMBL, .data[[id_col]])) %>%
  pivot_longer(-all_of(id_col), names_to = "sample", values_to = "SRGN") %>%
  select(sample, SRGN) %>%
  mutate(patient = sub(PATIENT_RX, "\\1", sample))

stopifnot(nrow(srgn) > 0)
message("SRGN expression rows: ", nrow(srgn))

# ---- 3. Clinical / overall survival from GDC (TCGAbiolinks) ------------------
gclin <- GDCquery_clinic(project = "TARGET-OS", type = "clinical")
gclin <- gclin[, !duplicated(names(gclin))]   # GDCquery_clinic can return dup cols

clin <- gclin %>%
  transmute(
    patient    = submitter_id,
    vital      = vital_status,
    d_death    = suppressWarnings(as.numeric(days_to_death)),
    d_fu       = suppressWarnings(as.numeric(days_to_last_follow_up)),
    metastasis = metastasis_at_diagnosis
  ) %>%
  filter(vital %in% c("Alive", "Dead")) %>%       # drop "Not Reported"/"Unknown"
  mutate(
    OS        = ifelse(vital == "Dead", 1L, 0L),
    # coalesce time: death time for the dead, follow-up time for the living
    OS_days   = ifelse(vital == "Dead", d_death, d_fu),
    OS_months = OS_days / 30.44
  ) %>%
  filter(!is.na(OS_months), OS_months >= 0) %>%
  distinct(patient, .keep_all = TRUE)

message(sprintf("Clinical (full registry): %d patients, %d deaths",
                nrow(clin), sum(clin$OS)))

# ---- 4. Merge expression + survival on the patient stem ---------------------
df <- srgn %>% inner_join(clin, by = "patient")
df$SRGN_z <- as.numeric(scale(df$SRGN))
message(sprintf("Analysis cohort: %d patients, %d events (%.0f%% mortality)",
                nrow(df), sum(df$OS), 100 * mean(df$OS)))

# ---- 5. PRIMARY — continuous Cox (HR per 1 SD of SRGN) ----------------------
cox_cont <- coxph(Surv(OS_months, OS) ~ SRGN_z, data = df)
sc <- summary(cox_cont)
print(sc)

res <- data.frame(
  model       = "continuous (per SD)",
  n           = sc$n,
  events      = sc$nevent,
  HR          = round(sc$conf.int[1, "exp(coef)"], 3),
  CI_low      = round(sc$conf.int[1, "lower .95"], 3),
  CI_high     = round(sc$conf.int[1, "upper .95"], 3),
  p           = signif(sc$coefficients[1, "Pr(>|z|)"], 3),
  concordance = round(sc$concordance[1], 3)
)

# ---- 6. Kaplan-Meier figure (median split — visual only) --------------------
df$SRGN_grp <- factor(ifelse(df$SRGN >= median(df$SRGN), "SRGN-high", "SRGN-low"),
                      levels = c("SRGN-low", "SRGN-high"))
fit <- survfit(Surv(OS_months, OS) ~ SRGN_grp, data = df)
km <- ggsurvplot(fit, data = df, pval = TRUE, risk.table = TRUE,
                 xlab = "Months", ylab = "Overall survival",
                 legend.title = "SRGN", palette = c("#377EB8", "#E41A1C"))

pdf(file.path(OUTDIR, "km_target_os_overall_survival.pdf"), width = 6, height = 7)
print(km, newpage = FALSE)
dev.off()

# ---- 7. Multivariable Cox — adjust for metastasis at diagnosis --------------
# TARGET labels: "Metastasis, NOS" / "No Metastasis" (plus NA/unknown -> dropped)
df_mv <- df %>%
  filter(metastasis %in% c("Metastasis, NOS", "No Metastasis")) %>%
  mutate(metastasis = factor(metastasis,
                             levels = c("No Metastasis", "Metastasis, NOS")))
cox_mv <- coxph(Surv(OS_months, OS) ~ SRGN_z + metastasis, data = df_mv)
smv <- summary(cox_mv)
print(smv)

res_mv <- data.frame(
  model       = "multivariable (SRGN + metastasis)",
  n           = smv$n,
  events      = smv$nevent,
  HR          = round(smv$conf.int["SRGN_z", "exp(coef)"], 3),
  CI_low      = round(smv$conf.int["SRGN_z", "lower .95"], 3),
  CI_high     = round(smv$conf.int["SRGN_z", "upper .95"], 3),
  p           = signif(smv$coefficients["SRGN_z", "Pr(>|z|)"], 3),
  concordance = round(smv$concordance[1], 3)
)

# ---- 8. Save results --------------------------------------------------------
out <- rbind(res, res_mv)
write.csv(out, file.path(OUTDIR, "srgn_target_os_cox_results.csv"), row.names = FALSE)
print(out)

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
message("Done. Outputs written to ", normalizePath(OUTDIR))

