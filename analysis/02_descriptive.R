# ===================================================================
# 02_descriptive.R -- describes dataset.csv and decides, for every variable
# and group, whether a parametric or a non-parametric test is justified
# Run from the repository root, after 01_cleaning.R
# ===================================================================

## CHECK WE ARE IN THE REPOSITORY ROOT
if (!file.exists("analysis/dataset.csv")) {
  stop("Not in the repository root, or analysis/dataset.csv is missing -- run 01_cleaning.R first.")
}

library(dplyr)
library(ggplot2)
library(xtable)


## FOLDERS
dir_dat <- "analysis"
dir_fig <- "report/figures"
dir_tab <- "report/tables"

for (d in c(dir_fig, dir_tab)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
}


## DATA
df <- read.csv(file.path(dir_dat, "dataset.csv"), sep = ",", header = TRUE)
glimpse(df)


## SETTINGS
n_boot     <- 5000  # number of bootstrap resamples (the guidelines ask for at least 5,000)
skew_limit <- 0.5   # DECISION: |skewness| below this counts as roughly symmetric
min_n_clt  <- 30    # rule of thumb for the central limit theorem, per group


## FUNCTIONS
# skewness: 0 for a symmetric distribution, positive for a long right tail
skewness <- function(x) {
  mean((x - mean(x))^3) / sd(x)^3
}

# resample the mean n_boot times, with the seed fixed so every run is the same
bootstrap_means <- function(x) {
  set.seed(123)
  replicate(n_boot, mean(sample(x, length(x), replace = TRUE)))
}


# ============================================================
# 1. ONE LONG TABLE: ONE LINE PER CAFE PER VARIABLE AND GROUP
# ============================================================
# Each line below is one row of the final table, and says which question needs it.
df_long <- bind_rows(
  # A1, A6, A7, A10
  df %>% transmute(variable = "Distance (km)", group = "All cafes", value = distance_km),
  # A2
  df %>% transmute(variable = "Distance (km)", group = ifelse(is_roastery == 1, "Roastery", "Not a roastery"), value = distance_km),
  # A4 (row 136 has no cafe type)
  df %>% filter(!is.na(cafe_type)) %>% transmute(variable = "Distance (km)", group = cafe_type, value = distance_km),
  # A9
  df %>% transmute(variable = "Distance (km)", group = ifelse(is_east == 1, "East", "West"), value = distance_km),
  # A6, A7, A8
  df %>% transmute(variable = "Number of ratings", group = "All cafes", value = n_ratings),
  # A9
  df %>% transmute(variable = "Number of ratings", group = ifelse(is_east == 1, "East", "West"), value = n_ratings),
  # A8
  df %>% transmute(variable = "Seats inside", group = "All cafes", value = seats_inside),
  # A3: the paired differences
  df %>% transmute(variable = "East-west minus north-south", group = "All cafes", value = east_west_abs - north_south_abs)
)

# a label for each row, kept in the order written above
df_long <- df_long %>%
  mutate(row_label = paste(variable, "-", group)) %>%
  mutate(row_label = factor(row_label, levels = unique(row_label)))

# the same table without missing values, for the plots and the checks
df_complete <- df_long %>% filter(!is.na(value))


# ============================================================
# 2. SUMMARY STATISTICS OF NUMERIC VARIABLES
# ============================================================
summary_table <- df_long %>%
  group_by(row_label) %>%
  summarise(
    n         = sum(!is.na(value)),
    n_missing = sum(is.na(value)),
    mean      = mean(value, na.rm = TRUE),
    sd        = sd(value, na.rm = TRUE),
    median    = median(value, na.rm = TRUE),
    q1        = quantile(value, 0.25, na.rm = TRUE),
    q3        = quantile(value, 0.75, na.rm = TRUE),
    iqr       = q3 - q1,
    min       = min(value, na.rm = TRUE),
    max       = max(value, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(
    se       = sd / sqrt(n),
    ci_lower = mean - qt(0.975, n - 1) * se,
    ci_upper = mean + qt(0.975, n - 1) * se
  ) %>%
  select(row_label, n, n_missing, mean, ci_lower, ci_upper, sd, median, iqr, min, max)

summary_table

summary_export <- summary_table
names(summary_export) <- c("Variable - group", "n", "Missing", "Mean", "CI lower", "CI upper",
                           "SD", "Median", "IQR", "Min", "Max")

print(xtable(summary_export),
      file = file.path(dir_tab, "D_summary.tex"),
      include.rownames = FALSE,
      size = "footnotesize",
      booktabs = TRUE,
      floating = FALSE)


# ============================================================
# 3. COUNTS AND PERCENTAGES OF CATEGORICAL VARIABLES
# ============================================================
category_table <- bind_rows(
  df %>% count(category = as.character(cafe_type)) %>% mutate(variable = "Cafe type"),
  df %>% count(category = as.character(district)) %>% mutate(variable = "District"),
  df %>% count(category = as.character(is_roastery)) %>% mutate(variable = "Roastery (1 = yes)"),
  df %>% count(category = as.character(is_prague1)) %>% mutate(variable = "Prague 1 (1 = yes)"),
  df %>% count(category = as.character(is_intl_chain)) %>% mutate(variable = "International chain (1 = yes)"),
  df %>% count(category = as.character(is_east)) %>% mutate(variable = "East of the Orloj (1 = yes)"),
  df %>% count(category = as.character(is_above_median)) %>% mutate(variable = "Above-median ratings (1 = yes)")
) %>%
  mutate(percent = 100 * n / nrow(df)) %>%
  select(variable, category, n, percent)

category_table

category_export <- category_table
names(category_export) <- c("Variable", "Category", "n", "Percent")

print(xtable(category_export),
      file = file.path(dir_tab, "D_categories.tex"),
      include.rownames = FALSE,
      size = "footnotesize",
      booktabs = TRUE,
      floating = FALSE)


# ============================================================
# 4. BOOTSTRAPPED MEANS FOR THE CENTRAL LIMIT THEOREM CHECK
# ============================================================
df_boot <- bind_rows(lapply(levels(df_complete$row_label), function(r) {
  values <- df_complete$value[df_complete$row_label == r]
  data.frame(row_label = r, boot_mean = bootstrap_means(values))
}))

df_boot <- df_boot %>%
  mutate(row_label = factor(row_label, levels = levels(df_complete$row_label)))

boot_skew_table <- df_boot %>%
  group_by(row_label) %>%
  summarise(boot_skew = skewness(boot_mean), .groups = "drop")


# ============================================================
# 5. THE DECISION TABLE
# ============================================================
# Parametric if the variable looks normal OR the central limit theorem applies.
decision_table <- df_complete %>%
  group_by(row_label) %>%
  summarise(n = n(), skew = skewness(value), .groups = "drop") %>%
  left_join(boot_skew_table, by = "row_label") %>%
  mutate(
    looks_normal = ifelse(abs(skew) < skew_limit, "Yes", "No"),
    clt_applies  = ifelse(n >= min_n_clt & abs(boot_skew) < skew_limit, "Yes", "No"),
    conclusion   = ifelse(looks_normal == "Yes" | clt_applies == "Yes", "Parametric", "Non-parametric")
  )

decision_table

decision_export <- decision_table
names(decision_export) <- c("Variable - group", "n", "Skewness", "Skewness of bootstrapped means",
                            "Looks normal", "CLT applies", "Conclusion")

print(xtable(decision_export),
      file = file.path(dir_tab, "D_table.tex"),
      include.rownames = FALSE,
      size = "footnotesize",
      booktabs = TRUE,
      floating = FALSE)


# ============================================================
# 6. PLOTS
# ============================================================
histograms <- ggplot(df_complete, aes(x = value)) +
  geom_histogram(bins = 20, fill = "#2a78d6") +
  facet_wrap(~ row_label, scales = "free") +
  labs(x = "Value", y = "Number of cafes") +
  theme_minimal()

ggsave(file.path(dir_fig, "D_histograms.pdf"), histograms, width = 11, height = 8)

# groups of the same variable side by side
boxplots <- ggplot(df_complete, aes(x = group, y = value)) +
  geom_boxplot() +
  coord_flip() +
  facet_wrap(~ variable, scales = "free") +
  labs(x = NULL, y = "Value") +
  theme_minimal()

ggsave(file.path(dir_fig, "D_boxplots.pdf"), boxplots, width = 11, height = 8)

qq_plots <- ggplot(df_complete, aes(sample = value)) +
  stat_qq() +
  stat_qq_line(color = "#e34948") +
  facet_wrap(~ row_label, scales = "free") +
  labs(x = "Expected if normal", y = "Observed") +
  theme_minimal()

ggsave(file.path(dir_fig, "D_qqplots.pdf"), qq_plots, width = 11, height = 8)

bootstrap_plots <- ggplot(df_boot, aes(x = boot_mean)) +
  geom_histogram(bins = 30, fill = "#2a78d6") +
  facet_wrap(~ row_label, scales = "free") +
  labs(x = "Bootstrapped mean", y = "Number of resamples") +
  theme_minimal()

ggsave(file.path(dir_fig, "D_bootstrap.pdf"), bootstrap_plots, width = 11, height = 8)
