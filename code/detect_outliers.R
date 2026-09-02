# function to detect outliers in continuous data in long form based on SD, IQR or MAD
# can specify level of dispersion desired 'k'
detect_outliers <- function(df, k = 0.995, type = "mad") {
  
  type <- match.arg(type, c("sd", "mad", "iqr"))
  
  # calculate descriptive stats
  # assumes df is in long form, with column 'vars' to analyse and values in column 'value'
  tmp <- df %>% 
    group_by(vars) %>% 
    summarise(
      mean   = mean(value, na.rm = TRUE),
      sd     = sd(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      mad    = median(abs(value - median(value, na.rm = TRUE)), na.rm = TRUE),
      q1     = quantile(value, 0.25, na.rm = TRUE),
      q3     = quantile(value, 0.75, na.rm = TRUE),
      iqr    = IQR(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  # merge to df and calculate measures of dispersion:
  # - z-score, for normally distributed date (number of sd away from the mean)
  # - interquartile, for skwed data: outside the IQR 
  # - using the MAD score & the median instead of mean/sd
  df <- left_join(df, tmp, by = "vars") %>% 
    mutate(
      z_score   = ifelse(sd  > 0, abs(value - mean)   / sd,             NA_real_),
      mad_score = ifelse(mad > 0, abs(value - median) / (1.4826 * mad), NA_real_),
      iqr_score = ifelse(iqr > 0, abs(value - median) / (iqr / 1.34898), NA_real_)
    )
  # scaling of mad/iqr scores taken from 
  #://real-statistics.com/descriptive-statistics/mad-and-outliers/
  
  
  # select desired level of upper/lower bounds based on type of estimator
  #qnorm(0.99) gives 1% FP rate (0.995 for 2 sided)
  if (type == "sd") {
    out <- df %>% filter(z_score > qnorm(k))
  } else if (type == "mad") {
    out <- df %>% filter(mad_score > qnorm(k))
  } else {
    out <- df %>% filter(iqr_score > qnorm(k))
  }
  
  out %<>% arrange(subjid, form_visit)
  
  return(out)
}