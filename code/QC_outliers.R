#### Outliers: identify using robust statistics 
# raise queries to be addressed by TC/DM:
#   1. blood pressure
#   2. demographics (age, pe_height, pe_weight, pe_bmi, bmi)
#   3. all continuous vars in "Baseline" section


# set the threshold for false positive rate: z_0.99 = 1% false positive rate
# Simon: "Maybe Z_0.995 to be 2-sided"
k <- 0.995
qnorm(k)

# clear the plot window
while (!is.null(dev.list())) dev.off()


# helper: collect variable metadata ----------------------------------------
collect_vars <- function(df, section, source_form = NULL) {
  
  if (is.null(source_form)) {
    if (!"source_form" %in% names(df)) {
      df <- df %>% mutate(source_form = NA_character_)
    }
  } else {
    df <- df %>% mutate(source_form = source_form)
  }
  
  if (!"units" %in% names(df)) {
    df <- df %>% mutate(units = NA_character_)
  }
  
  if (!"form_visit" %in% names(df)) {
    df <- df %>% mutate(form_visit = NA_character_)
  }
  
  df %>%
    mutate(
      vars       = as.character(vars),
      units      = as.character(units),
      form_visit = as.character(form_visit),
      units      = suppressWarnings(iconv(units, from = "", to = "ASCII", sub = "u")),
      units      = trimws(units),
      form_visit = trimws(form_visit),
      units      = dplyr::na_if(units, ""),
      form_visit = dplyr::na_if(form_visit, "")
    ) %>%
    mutate(section = section, .before = 1) %>%
    group_by(section, source_form, vars) %>%
    summarise(
      n_items    = sum(!is.na(value)),
      units      = paste(unique(stats::na.omit(units)), collapse = "; "),
      form_visit = paste(unique(stats::na.omit(form_visit)), collapse = "; "),
      .groups = "drop"
    ) %>%
    mutate(
      units      = dplyr::na_if(units, ""),
      form_visit = dplyr::na_if(form_visit, "")
    ) %>%
    select(section, source_form, vars, units, form_visit, n_items)
}

# initialise variable catalogue --------------------------------------------
vars.catalog <- tibble(
  section = character(),
  source_form = character(),
  vars = character(),
  units = character(),
  form_visit = character(),
  n_items = integer()
)

# initialise outlier summary messages --------------------------------------
outlier_msgs <- character()



# 1. blood pressure ---------------------------------------------------------

## all blood pressure measuraments from eForm: "OMRONRQG"
## (SEATED AutomatedOfficeBloodPressure and HR (OMRON))
qc.bp <- extract_form(form = "OMRONRQG", data = data_long, vars_keep = c("subjid", "study")) %>% 
  lab2val %>% 
  select(-dm_screenidder, -dm_subjidder, -omron_pulseyn, -omron_autobpdat)

# nurse in room (long form to wide)
qc.bp1 <- qc.bp %>% 
  pivot_wider(
    id_cols = c("subjid", "form_visit"),
    names_from = "omron_bp",
    values_from = "omron_avr2_3"
  )

# nurse out of room (first line only)
qc.bp2 <- qc.bp %>% 
  select(subjid, omron_avsysbp:omron_avhr, form_visit) %>%
  group_by(subjid, form_visit) %>% 
  slice(1) %>%
  ungroup()

# merge together & pivot to longer form
qc.bp <- full_join(qc.bp1, qc.bp2, by = c("subjid", "form_visit")) %>% 
  clean_names() %>% 
  pivot_longer(cols = -c("subjid", "form_visit"), names_to = "vars", values_to = "value")

# collect variable names
vars.catalog <- bind_rows(
  vars.catalog,
  collect_vars(qc.bp, section = "blood pressure", source_form = "OMRONRQG")
)

# total number of values
qc.bp %>% 
  filter(!is.na(value)) %>% 
  {table(.$vars, .$form_visit)} %>% 
  addmargins()

n <- qc.bp %>% filter(!is.na(value)) %>% nrow()

# visualise distribution
for (i in unique(qc.bp$vars)) {
  print(i)
  x <- qc.bp %>% filter(vars == i) %>% ungroup()
  hist(x$value, main = paste("Histogram of", i), xlab = NULL, breaks = 50)
  abline(v = mean(x$value, na.rm = TRUE), col = "blue")
  abline(v = median(x$value, na.rm = TRUE), col = "red")
}
rm(x)

# detects outliers
outliers1 <- detect_outliers(qc.bp, k, type = "sd")
outliers2 <- detect_outliers(qc.bp, k, type = "iqr")
outliers.bp <- detect_outliers(qc.bp, k, type = "mad")

msg <- paste0(
  "OMRON - outliers detected: ", nrow(outliers.bp),
  " (", round(100 * nrow(outliers.bp) / n, 1), "%)"
)
print(msg)
outlier_msgs <- c(outlier_msgs, msg)

rm(list = ls(pattern = "^qc\\."))



# 2. demographics -----------------------------------------------------------

## continuous covariates (demographics: age, height, weight, bmi)
qc.demo <- covars %>% 
  select(subjid, age, pe_height, pe_weight, pe_bmi, bmi) %>% 
  pivot_longer(cols = -c("subjid"), names_to = "vars", values_to = "value")

qc.demo$form_visit <- "V1"

# collect variable names
vars.catalog <- bind_rows(
  vars.catalog,
  collect_vars(qc.demo, section = "demographics", source_form = "covars")
)

# total number of values
n <- qc.demo %>% filter(!is.na(value)) %>% nrow()

# visualise distribution
for (i in unique(qc.demo$vars)) {
  print(i)
  x <- qc.demo %>% filter(vars == i) %>% ungroup()
  hist(x$value, main = paste("Histogram of", i), xlab = NULL, breaks = 50)
  abline(v = mean(x$value, na.rm = TRUE), col = "blue")
  abline(v = median(x$value, na.rm = TRUE), col = "red")
}
rm(x)

# detects outliers
outliers1 <- detect_outliers(qc.demo, k, type = "sd")
outliers2 <- detect_outliers(qc.demo, k, type = "iqr")
outliers.demo <- detect_outliers(qc.demo, k, type = "mad")

msg <- paste0(
  "Demographics - outliers detected: ", nrow(outliers.demo),
  " (", round(100 * nrow(outliers.demo) / n, 1), "%)"
)
print(msg)
outlier_msgs <- c(outlier_msgs, msg)

rm(list = ls(pattern = "^qc\\."))



# 3. baseline ---------------------------------------------------------------

## Biochemistry Bloods from eForm: "BLDBCHM" (Bloods - Biochemistry)
qc.blds <- extract_form(form = "BLDBCHM", data = data_long, vars_keep = c("subjid", "study")) %>% 
  lab2val %>%  
  select(subjid, study, lbtestcd, lbtestrsl, lbtestunit, form_visit) %>% 
  rename(vars = lbtestcd, value = lbtestrsl, units = lbtestunit) %>%
  mutate(source_form = "BLDBCHM")

## Plasma Renin/Aldosterone from eForm: "PRA" (Plasma Renin and Aldosterone)
qc.pra <- extract_form(form = "PRA", data = data_long, vars_keep = c("subjid", "study")) %>% 
  lab2val %>% 
  select(subjid, study, praloctestcd, pralocrsl, pralocunit, pralocunitoth, form_visit) %>%
  mutate(pralocunit = ifelse(pralocunit == "Other", pralocunitoth, pralocunit)) %>% 
  select(-pralocunitoth) %>% 
  rename(vars = praloctestcd, value = pralocrsl, units = pralocunit) %>%
  mutate(source_form = "PRA")

## 24hr Urine Sodium/Potassium from eForm: "H24URNNAK" (24 Hour Urine Sodium/Potassium)
qc.urk <- extract_form(form = "H24URNNAK", data = data_long, vars_keep = c("subjid", "study")) %>% 
  lab2val %>% 
  select(subjid, study, urn24htestcd, urn24hrsl, urn24hunit, urn24hunitoth, form_visit) %>%
  mutate(urn24hunit = ifelse(urn24hunit == "Other", urn24hunitoth, urn24hunit)) %>% 
  select(-urn24hunitoth) %>% 
  rename(vars = urn24htestcd, value = urn24hrsl, units = urn24hunit) %>%
  mutate(source_form = "H24URNNAK")

## TANITA Body Composition from eForm: "TANITA" (Body Composition measure with TANITA scales)
qc.tanita <- extract_form(form = "TANITA", data = data_long, vars_keep = c("subjid", "study")) %>% 
  lab2val %>% 
  select(subjid, study, tanitatestcd, tanitarsl, tanitaunit, form_visit) %>% 
  rename(vars = tanitatestcd, value = tanitarsl, units = tanitaunit) %>%
  mutate(source_form = "TANITA")

## Pulse Wave Analysis from eForm: "PWA" (Pulse Wave Analysis (PWA))
qc.pwa <- extract_form(form = "PWA", data = data_long, vars_keep = c("subjid", "study")) %>% 
  lab2val %>% 
  select(subjid, study, pwatestcd, pwarslav, pwaunit, form_visit) %>% 
  rename(vars = pwatestcd, value = pwarslav, units = pwaunit) %>%
  mutate(source_form = "PWA")

## Aortic Pulse Wave Velocity from eForm: "APWV" (Aortic Pulse Wave Velocity (APWV))
qc.apwv_raw <- extract_form(form = "APWV", data = data_long, vars_keep = c("subjid", "study")) %>% 
  lab2val %>% 
  select(subjid, study, pltestcd, plrsl, plunit, apwvrslav, apwvunit, form_visit)

# reshape data (mixture of long/wide)
qc.apwv <- bind_rows(
  qc.apwv_raw %>% 
    select(subjid, study, pltestcd, plrsl, plunit, form_visit) %>%
    mutate(source_form = "APWV"),
  qc.apwv_raw %>% 
    select(subjid, study, apwvrslav, apwvunit, form_visit) %>% 
    filter(!is.na(apwvrslav)) %>% 
    rename(plrsl = apwvrslav, plunit = apwvunit) %>% 
    mutate(
      pltestcd = "Carotid to Femoral (average)",
      source_form = "APWV"
    )
) %>% 
  rename(vars = pltestcd, value = plrsl, units = plunit)

## Dundee Step Test from eForm: "DUNDEE" (Dundee Step Test)
qc.dundee <- extract_form(form = "DUNDEE", data = data_long, vars_keep = c("subjid", "study")) %>% 
  lab2val %>% 
  select(subjid, study, dsttestcd, dstrslav, dstrslaft, dstunit, form_visit)

qc.dundee2 <- bind_rows(
  qc.dundee %>% 
    select(subjid, study, dsttestcd, dstrslav, dstunit, form_visit) %>% 
    mutate(
      vars = paste0("DUNDEE_before_", dsttestcd),
      source_form = "DUNDEE"
    ) %>% 
    select(-dsttestcd) %>% 
    rename(value = dstrslav, units = dstunit),
  qc.dundee %>% 
    select(subjid, study, dsttestcd, dstrslaft, dstunit, form_visit) %>% 
    mutate(
      vars = paste0("DUNDEE_after_", dsttestcd),
      source_form = "DUNDEE"
    ) %>% 
    select(-dsttestcd) %>% 
    rename(value = dstrslaft, units = dstunit)
)

## ECHO from eForm "CO" (Echocardiograph (ECHO))
## only records whether test was performed
qc.co <- extract_form(form = "CO", data = data_long, vars_keep = c("subjid", "study")) %>% 
  lab2val %>% 
  select(subjid, study, cocd, corsl, counit, form_visit) %>% 
  rename(vars = cocd, value = corsl, units = counit) %>%
  mutate(source_form = "CO")

qc.df <- bind_rows(
  qc.blds,       # Biochemistry Bloods
  qc.pra,        # Plasma Renin/Aldosterone
  qc.urk,        # 24hr Urine Sodium/Potassium
  qc.tanita,     # TANITA Body Composition
  qc.pwa,        # Pulse Wave Analysis
  qc.apwv,       # Aortic Pulse Wave Velocity
  qc.dundee2,    # Dundee Step Test
  qc.co          # ECHO
)

# collect variable names
vars.catalog <- bind_rows(
  vars.catalog,
  collect_vars(qc.df, section = "baseline")
)

# total number of values
qc.df %>% 
  filter(!is.na(value)) %>% 
  {table(.$vars, .$form_visit)} %>% 
  addmargins()

n <- qc.df %>% filter(!is.na(value)) %>% nrow()

# visualise distributions
for (i in unique(qc.df$vars)) {
  print(i)
  x <- qc.df %>% filter(vars == i) %>% ungroup()
  hist(x$value, main = paste("Histogram of", i), xlab = NULL, breaks = 50)
  abline(v = mean(x$value, na.rm = TRUE), col = "blue")
  abline(v = median(x$value, na.rm = TRUE), col = "red")
}
rm(x)

# detects outliers
outliers1 <- detect_outliers(qc.df, k, type = "sd")
outliers2 <- detect_outliers(qc.df, k, type = "iqr")
outliers.df <- detect_outliers(qc.df, k, type = "mad")

msg <- paste0(
  "Baseline - outliers detected: ", nrow(outliers.df),
  " (", round(100 * nrow(outliers.df) / n, 1), "%)"
)
print(msg)
outlier_msgs <- c(outlier_msgs, msg)

# fix umol
outliers.df %<>% mutate(units = iconv(units, "latin1", "ASCII", sub = "u"))

rm(list = ls(pattern = "^qc\\."))



# 4. biochemistry -----------------------------------------------------------

## bloods from eForm: "SAFEBL" (Bloods - Biochemistry)
# 1: Sodium (mmol/L) - 1 in umol/L
# 2: Potassium (mmol/L) - all good
# 3: Creatinine (umol/L) - 8 in mmol/L
qc.sf <- extract_form(form = "SAFEBL", data = data_long, vars_keep = c("subjid", "study")) %>% 
  lab2val %>% 
  select(subjid, study, sbltestcd, sbtestrsl, sbtestunit, form_visit) %>% 
  rename(vars = sbltestcd, value = sbtestrsl, units = sbtestunit) %>%
  mutate(source_form = "SAFEBL")

# collect variable names
vars.catalog <- bind_rows(
  vars.catalog,
  collect_vars(qc.sf, section = "biochemistry")
)

# total number of values
qc.sf %>% 
  filter(!is.na(value)) %>% 
  {table(.$vars, .$form_visit)} %>% 
  addmargins()

n <- qc.sf %>% filter(!is.na(value)) %>% nrow()

table(qc.sf$vars, qc.sf$units)

# visualise distributions
for (i in unique(qc.sf$vars)) {
  print(i)
  x <- qc.sf %>% filter(vars == i) %>% ungroup()
  hist(x$value, main = paste("Histogram of", i), xlab = NULL, breaks = 50)
  abline(v = mean(x$value, na.rm = TRUE), col = "blue")
  abline(v = median(x$value, na.rm = TRUE), col = "red")
}
rm(x)

# detects outliers
outliers1 <- detect_outliers(qc.sf, k, type = "sd")
outliers2 <- detect_outliers(qc.sf, k, type = "iqr")
outliers.sf <- detect_outliers(qc.sf, k, type = "mad")

msg <- paste0(
  "Biochemistry - outliers detected: ", nrow(outliers.sf),
  " (", round(100 * nrow(outliers.sf) / n, 1), "%)"
)
print(msg)
outlier_msgs <- c(outlier_msgs, msg)

# fix umol
outliers.sf %<>% mutate(units = iconv(units, "latin1", "ASCII", sub = "u"))

rm(list = ls(pattern = "^qc\\."))


# finalise variable catalogue ----------------------------------------------

vars.catalog <- vars.catalog %>%
  distinct() %>%
  arrange(section, source_form, vars)



# write out ----------------------------------------------------------------

write_xlsx(
  list(
    bp       = outliers.bp   %>% arrange(-mad_score),
    demo     = outliers.demo %>% arrange(-mad_score),
    baseline = outliers.df   %>% arrange(-mad_score),
    biochem  = outliers.sf   %>% arrange(-mad_score),
    vars     = vars.catalog
  ),
  path = file.path("Output","Queries",paste0("outliers_",substr(DATA_PREFIX,13,99),".xlsx")),
  format_headers = FALSE
)


# summary of MACRO DD -----------------------------------------------------
outfile <- file.path("Output", "Queries", paste0("outliers_msg_", substr(DATA_PREFIX, 13, 99), ".txt"))

# MACRO details
study_counts <- data_long %>%
  count(study) %>%
  mutate(study = recode(study, D = "DUAL", M = "MONO"))

subj_counts <- data_long %>%
  distinct(study, subjid) %>%
  count(study) %>%
  mutate(study = recode(study, D = "DUAL", M = "MONO"))

N_rows      <- nrow(data_long)
N_subjid    <- sum(subj_counts$n)
N_dual_rows <- study_counts$n[study_counts$study == "DUAL"]
N_mono_rows <- study_counts$n[study_counts$study == "MONO"]
N_dual_subj <- subj_counts$n[subj_counts$study == "DUAL"]
N_mono_subj <- subj_counts$n[subj_counts$study == "MONO"]

writeLines(
  c(
    paste0("## Lock_Final_", DATA_PREFIX),
    paste0(
      "Total ", N_rows, " data rows in ", N_subjid,
      " subjects (DUAL=", N_dual_subj, ", MONO=", N_mono_subj, ")"
    ),
    "distributed across the two studies:",
    paste0("DUAL=", N_dual_rows, " data rows"),
    paste0("MONO=", N_mono_rows, " data rows")
  ),
  con = outfile
)

# details of variables checked for outliers
N_vars     <- nrow(vars.catalog)
N_items    <- sum(vars.catalog$n_items, na.rm = TRUE)
N_outliers <- nrow(outliers.bp) + nrow(outliers.demo) + nrow(outliers.df) + nrow(outliers.sf)

write(
  c(
    "",
    "## Variable Catalogue",
    paste0("Total variables checked=", N_vars),
    paste0("Total data points analysed=", N_items)
  ),
  file = outfile,
  append = TRUE
)

write(
  c(
    "",
    "## Outlier Summary",
    outlier_msgs,
    "",
    paste0(
      "Total - outliers detected: ", N_outliers,
      " (", round(100 * N_outliers / N_items, 1), "%)"
    )
  ),
  file = outfile,
  append = TRUE
)

message("QC of outliers completed successfully for MACRO ", DATA_PREFIX)


# clean-up -----------------------------------------------------------------
rm(list = ls(pattern = "^outliers"))
rm(list = ls(pattern = "_counts$"))
rm(list = ls(pattern = "^N_"))
rm(
  i, k, n, msg, outfile,
  collect_vars,
  vars.catalog,
  outlier_msgs
)

