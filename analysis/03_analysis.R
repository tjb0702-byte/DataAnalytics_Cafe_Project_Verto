# ===================================================================
# 03_analysis.R -- Location topic, research questions A1 to A10
# Run from the repository root, after 01_cleaning.R and 02_descriptive.R
#
# Significance level: 0.05. Cafes in the same district are not independent,
# so the p-values below may be too small. Discuss this in the report.
# ===================================================================

## CHECK WE ARE IN THE REPOSITORY ROOT
if (!file.exists("analysis/03_analysis.R")) {
  stop("Not in the repository root -- set the working directory to the repository root, then run again.")
}

library(dplyr)
library(ggplot2)
library(xtable)
library(rstatix)          # levene_test, dunn_test, cramer_v
library(car)              # vif, avPlots
library(lmtest)           # coeftest, bptest
library(sandwich)         # vcovHC (robust standard errors)
library(marginaleffects)  # avg_slopes


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

n_boot <- 5000  # number of bootstrap resamples (the guidelines ask for at least 5,000)


## FUNCTIONS USED BY SEVERAL QUESTIONS

# p-value as text, three decimals and no leading zero.
# 1000 + thousandths always has four digits (0.09 becomes 1090), so replacing
# the first digit with a dot gives .090 and keeps the final zero.
p_text <- function(p) {
  thousandths <- round(p * 1000)
  case_when(
    p < 0.001          ~ "<.001",
    thousandths == 1000 ~ "1.000",
    TRUE               ~ sub("^1", ".", 1000 + thousandths)
  )
}

# floating = FALSE leaves out \begin{table}, because the report adds its own
save_table <- function(table_data, file_name) {
  print(xtable(table_data),
        file = file.path(dir_tab, file_name),
        include.rownames = FALSE,
        size = "footnotesize",
        booktabs = TRUE,
        floating = FALSE)
}

# the median of x with a bootstrapped 95% interval
boot_median_ci <- function(x) {
  set.seed(123)
  b <- replicate(n_boot, median(sample(x, length(x), replace = TRUE)))
  data.frame(med   = median(x),
             lower = unname(quantile(b, 0.025)),
             upper = unname(quantile(b, 0.975)))
}

# bootstrapped differences in medians, x minus y
bootstrap_median_diff <- function(x, y) {
  set.seed(123)
  replicate(n_boot, median(sample(x, length(x), replace = TRUE)) -
              median(sample(y, length(y), replace = TRUE)))
}

# one line per group: the median distance and its bootstrapped interval
group_medians <- function(data) {
  result <- bind_rows(lapply(levels(data$group), function(g) {
    values <- data$distance_km[data$group == g]
    cbind(group = g, boot_median_ci(values))
  }))
  result$group <- factor(result$group, levels = levels(data$group))
  result
}

# the cafes as points, with the median of each group and its interval
plot_group_medians <- function(data, medians) {
  set.seed(1100)
  ggplot(data, aes(x = group, y = distance_km)) +
    geom_jitter(width = 0.2, height = 0, color = "grey60", alpha = 0.6) +
    geom_errorbar(data = medians, aes(y = med, ymin = lower, ymax = upper),
                  width = 0.15, linewidth = 1, color = "#2a78d6") +
    geom_point(data = medians, aes(y = med), size = 3, color = "#2a78d6") +
    labs(x = NULL, y = "Distance from the Orloj (km)",
         title = "Median distance (blue) with bootstrapped 95% CI") +
    theme_minimal()
}

# a regression table with robust standard errors and 95% intervals
regression_table <- function(fit, model_name) {
  robust <- coeftest(fit, vcov = vcovHC(fit, type = "HC3"))
  result <- data.frame(
    model_name,
    rownames(robust),
    robust[, 1],
    robust[, 2],
    robust[, 1] - qnorm(0.975) * robust[, 2],
    robust[, 1] + qnorm(0.975) * robust[, 2],
    robust[, 3],
    p_text(robust[, 4]),
    length(residuals(fit)),
    summary(fit)$r.squared,
    p_text(bptest(fit)$p.value)
  )
  names(result) <- c("Model", "Term", "Estimate", "Robust SE", "CI lower", "CI upper",
                     "t", "p", "n", "R squared", "Breusch-Pagan p")
  result
}

# the four standard diagnostic plots of a regression, saved to one file
save_diagnostics <- function(fit, file_name) {
  pdf(file.path(dir_fig, file_name))
  par(mfrow = c(2, 2))
  plot(fit)
  par(mfrow = c(1, 1))
  invisible(dev.off())
}

# the five cafes with the largest Cook's distance
influential_cafes <- function(data, cooks) {
  data$cooks <- cooks
  result <- data %>%
    arrange(desc(cooks)) %>%
    select(sheet_row, cafe_name, n_ratings, distance_km, cooks) %>%
    head(5)
  names(result) <- c("Sheet row", "Cafe", "Number of ratings", "Distance (km)", "Cook's distance")
  result
}


# ===================================================================
# A1 -- Is the average distance of cafes from the Orloj different from 2 km?
# ===================================================================

# H0: The mean distance of cafes from the Orloj equals 2 km.
# H1: The mean distance of cafes from the Orloj is different from 2 km.
# Tails and why: Two-tailed, because the question asks "different from" and names no direction.

# Method: one-sample t-test. D_table says parametric for distance (all cafes),
# because the sample is large enough for the central limit theorem.

a1_benchmark <- 2

# keep only the cafes with a known distance
a1_distance  <- df$distance_km[!is.na(df$distance_km)]
a1_n         <- length(a1_distance)
a1_n_missing <- sum(is.na(df$distance_km))

a1_test <- t.test(a1_distance, mu = a1_benchmark)
a1_test

a1_mean   <- mean(a1_distance)
a1_lower  <- a1_test$conf.int[1]
a1_upper  <- a1_test$conf.int[2]
a1_median <- median(a1_distance)
a1_iqr    <- quantile(a1_distance, 0.75) - quantile(a1_distance, 0.25)

# Cohen's d: how many standard deviations the mean is from the benchmark
a1_d <- (a1_mean - a1_benchmark) / sd(a1_distance)

## A1 TABLE
a1_table <- data.frame(
  a1_n,
  a1_n_missing,
  a1_mean,
  a1_lower,
  a1_upper,
  a1_median,
  a1_iqr,
  a1_test$statistic,
  as.integer(a1_test$parameter),
  p_text(a1_test$p.value),
  a1_d
)

names(a1_table) <- c("n", "Missing", "Mean (km)", "CI lower (km)", "CI upper (km)",
                     "Median (km)", "IQR (km)", "t", "df", "p", "Cohen's d")

save_table(a1_table, "A1_table.tex")

## A1 FIGURE
a1_points  <- data.frame(distance_km = a1_distance)
a1_summary <- data.frame(mean = a1_mean, lower = a1_lower, upper = a1_upper)

# the seed makes the random sideways spread of the points the same every run
set.seed(1100)

a1_figure <- ggplot(a1_points, aes(x = "Cafes", y = distance_km)) +
  geom_jitter(width = 0.2, height = 0, color = "grey60", alpha = 0.6) +
  geom_hline(yintercept = a1_benchmark, linetype = "dashed", color = "#e34948") +
  geom_errorbar(data = a1_summary, aes(y = mean, ymin = lower, ymax = upper),
                width = 0.1, linewidth = 1, color = "#2a78d6") +
  geom_point(data = a1_summary, aes(y = mean), size = 3, color = "#2a78d6") +
  labs(x = NULL, y = "Distance from the Orloj (km)",
       title = "Mean distance (blue, 95% CI) vs. the 2 km benchmark (dashed)") +
  theme_minimal()

ggsave(file.path(dir_fig, "A1_figure.pdf"), a1_figure, width = 6, height = 4)


# ===================================================================
# A2 -- Are roasteries farther from the Orloj than other cafes?
# ===================================================================

# H0: The mean distance from the Orloj is the same for roasteries and other cafes.
# H1: The mean distance from the Orloj is greater for roasteries than for other cafes.
# Tails and why: One-tailed, because the question says "farther".

# Method: Wilcoxon rank-sum test. D_table says non-parametric for roasteries:
# the group is under 30 cafes and skewed. A rank test compares typical
# values, so read H0 and H1 as being about the typical distance, and the
# estimates below are medians, not means.

# Roastery is the first level, so "greater" means roastery minus other
df_a2 <- df %>%
  mutate(group = ifelse(is_roastery == 1, "Roastery", "Not a roastery")) %>%
  mutate(group = factor(group, levels = c("Roastery", "Not a roastery")))

a2_roastery <- df_a2$distance_km[df_a2$group == "Roastery"]
a2_other    <- df_a2$distance_km[df_a2$group == "Not a roastery"]

# Assumption check: similar spread in the two groups
a2_levene <- df_a2 %>% levene_test(distance_km ~ group)
a2_levene

# One-tailed test for the p-value
a2_test <- wilcox.test(distance_km ~ group, data = df_a2, alternative = "greater")
a2_test

# Two-sided call only for the 95% interval of the difference
a2_interval <- wilcox.test(distance_km ~ group, data = df_a2, conf.int = TRUE)
a2_interval

## A2 TABLE
a2_table <- data.frame(
  length(a2_roastery),
  length(a2_other),
  median(a2_roastery),
  median(a2_other),
  unname(a2_interval$estimate),
  a2_interval$conf.int[1],
  a2_interval$conf.int[2],
  unname(a2_test$statistic),
  p_text(a2_test$p.value),
  p_text(a2_levene$p)
)

names(a2_table) <- c("n roastery", "n other", "Median roastery (km)", "Median other (km)",
                     "Difference (km)", "CI lower (km)", "CI upper (km)",
                     "W", "p (one-tailed)", "Levene p")

save_table(a2_table, "A2_table.tex")

## A2 FIGURE
a2_medians <- group_medians(df_a2)
a2_figure  <- plot_group_medians(df_a2, a2_medians)

ggsave(file.path(dir_fig, "A2_figure.pdf"), a2_figure, width = 6, height = 4)


# ===================================================================
# A3 -- Are cafes spread farther east-west than north-south from the Orloj?
# ===================================================================

# H0: The mean of (absolute east-west km minus absolute north-south km) equals 0.
# H1: The mean of (absolute east-west km minus absolute north-south km) is greater than 0.
# Tails and why: One-tailed, because the question says "farther east-west".

# Method: paired t-test. D_table says parametric for the differences.
# Difference = absolute east-west km minus absolute north-south km.

a3_difference <- df$east_west_abs - df$north_south_abs
a3_n          <- length(a3_difference)

# One-tailed test for the p-value
a3_test <- t.test(df$east_west_abs, df$north_south_abs, paired = TRUE, alternative = "greater")
a3_test

# Two-sided call only for the 95% interval
a3_interval <- t.test(df$east_west_abs, df$north_south_abs, paired = TRUE)

a3_mean  <- mean(a3_difference)
a3_lower <- a3_interval$conf.int[1]
a3_upper <- a3_interval$conf.int[2]

# Cohen's d_z: the mean difference in standard deviations of the differences
a3_dz <- a3_mean / sd(a3_difference)

## A3 TABLE
a3_table <- data.frame(
  a3_n,
  a3_mean,
  a3_lower,
  a3_upper,
  a3_test$statistic,
  as.integer(a3_test$parameter),
  p_text(a3_test$p.value),
  a3_dz
)

names(a3_table) <- c("n", "Mean difference (km)", "CI lower (km)", "CI upper (km)",
                     "t", "df", "p (one-tailed)", "Cohen's dz")

save_table(a3_table, "A3_table.tex")

## A3 FIGURE
a3_points  <- data.frame(difference = a3_difference)
a3_summary <- data.frame(mean = a3_mean, lower = a3_lower, upper = a3_upper)

set.seed(1100)

a3_figure <- ggplot(a3_points, aes(x = "Cafes", y = difference)) +
  geom_jitter(width = 0.2, height = 0, color = "grey60", alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#e34948") +
  geom_errorbar(data = a3_summary, aes(y = mean, ymin = lower, ymax = upper),
                width = 0.1, linewidth = 1, color = "#2a78d6") +
  geom_point(data = a3_summary, aes(y = mean), size = 3, color = "#2a78d6") +
  labs(x = NULL, y = "East-west km minus north-south km",
       title = "Mean difference (blue, 95% CI) vs. no difference (dashed)") +
  theme_minimal()

ggsave(file.path(dir_fig, "A3_figure.pdf"), a3_figure, width = 6, height = 4)


# ===================================================================
# A4 -- Does the distance from the Orloj differ by cafe type, and how?
# ===================================================================

# H0: The mean distance from the Orloj is the same for all three cafe types.
# H1: The mean distance of at least one cafe type differs from another.
# Tails and why: Non-directional, because the question says "differ" and names no order.

# Method: Kruskal-Wallis test, then Dunn's test with Holm adjustment.
# D_table says non-parametric for international chains, the smallest group.
# A rank test compares typical values, so read H0 and H1 as being about the
# typical distance. Read the pairwise results only if the Kruskal-Wallis
# test is significant.

df_a4 <- df %>%
  filter(!is.na(cafe_type)) %>%
  mutate(group = factor(cafe_type, levels = c("Independent", "Local Chain", "International Chain")))

a4_n_missing <- nrow(df) - nrow(df_a4)

# Assumption check: similar spread in the three groups
a4_levene <- df_a4 %>% levene_test(distance_km ~ group)
a4_levene

a4_test <- kruskal.test(distance_km ~ group, data = df_a4)
a4_test

## A4 TABLES
a4_table <- data.frame(
  nrow(df_a4),
  a4_n_missing,
  unname(a4_test$statistic),
  as.integer(a4_test$parameter),
  p_text(a4_test$p.value),
  p_text(a4_levene$p)
)

names(a4_table) <- c("n", "Missing", "H", "df", "p", "Levene p")
save_table(a4_table, "A4_table.tex")

a4_medians <- group_medians(df_a4)
a4_groups  <- a4_medians %>%
  left_join(df_a4 %>% count(group), by = "group") %>%
  select(group, n, med, lower, upper)

names(a4_groups) <- c("Cafe type", "n", "Median (km)", "CI lower (km)", "CI upper (km)")
save_table(a4_groups, "A4_groups.tex")

a4_dunn <- df_a4 %>% dunn_test(distance_km ~ group, p.adjust.method = "holm")
a4_dunn

# for each pair: group2 minus group1, with a bootstrapped interval
a4_differences <- bind_rows(lapply(seq_len(nrow(a4_dunn)), function(i) {
  values_2 <- df_a4$distance_km[df_a4$group == a4_dunn$group2[i]]
  values_1 <- df_a4$distance_km[df_a4$group == a4_dunn$group1[i]]
  boot     <- bootstrap_median_diff(values_2, values_1)
  data.frame(difference = median(values_2) - median(values_1),
             lower      = unname(quantile(boot, 0.025)),
             upper      = unname(quantile(boot, 0.975)))
}))

a4_pairwise <- data.frame(
  paste(a4_dunn$group2, "minus", a4_dunn$group1),
  a4_differences$difference,
  a4_differences$lower,
  a4_differences$upper,
  a4_dunn$statistic,
  p_text(a4_dunn$p.adj)
)

names(a4_pairwise) <- c("Comparison", "Median difference (km)", "CI lower (km)", "CI upper (km)",
                        "z", "Holm-adjusted p")

save_table(a4_pairwise, "A4_pairwise.tex")

## A4 FIGURE
a4_figure <- plot_group_medians(df_a4, a4_medians)

ggsave(file.path(dir_fig, "A4_figure.pdf"), a4_figure, width = 6, height = 4)


# ===================================================================
# A5 -- Is an international chain more likely to be in Prague 1 than
#       other cafe types, and by how many percentage points?
# ===================================================================

# H0: The share of cafes in Prague 1 is the same for international chains and other cafe types.
# H1: The share of cafes in Prague 1 is higher for international chains than for other cafe types.
# Tails and why: One-tailed, because the question says "more likely".

# Method: Fisher's exact test, because the smallest expected count is below 5.

# International chain and Prague 1 come first, so "greater" means more
# international chains in Prague 1
df_a5 <- df %>%
  filter(!is.na(is_intl_chain)) %>%
  mutate(
    chain    = factor(is_intl_chain, levels = c(1, 0), labels = c("International chain", "Other cafe type")),
    location = factor(is_prague1, levels = c(1, 0), labels = c("Prague 1", "Elsewhere"))
  )

a5_n_missing <- nrow(df) - nrow(df_a5)

a5_counts <- with(df_a5, table(chain, location))
a5_counts

a5_row_percent <- 100 * prop.table(a5_counts, margin = 1)
a5_row_percent

# Assumption check: expected counts (warning switched off because they are small)
a5_chi <- suppressWarnings(chisq.test(a5_counts, correct = FALSE))
a5_smallest_expected <- min(a5_chi$expected)
a5_smallest_expected

a5_v <- suppressWarnings(cramer_v(a5_counts))

# standardised residuals: above 2 in absolute value drives the association
round(a5_chi$stdres, 2)

# One-tailed test for the p-value
a5_test <- fisher.test(a5_counts, alternative = "greater")
a5_test

# Two-sided call only for the 95% interval of the odds ratio
a5_interval <- fisher.test(a5_counts)

# The difference in percentage points, international chains minus the rest.
# With only a few international chains this interval is approximate.
a5_prop <- suppressWarnings(prop.test(a5_counts[, "Prague 1"], rowSums(a5_counts), correct = FALSE))

a5_share_intl  <- 100 * a5_prop$estimate[1]
a5_share_other <- 100 * a5_prop$estimate[2]
a5_difference  <- a5_share_intl - a5_share_other
a5_diff_lower  <- 100 * a5_prop$conf.int[1]
a5_diff_upper  <- 100 * a5_prop$conf.int[2]

## A5 TABLES
a5_count_table <- data.frame(
  c("International chain", "Other cafe type"),
  c(a5_counts["International chain", "Prague 1"], a5_counts["Other cafe type", "Prague 1"]),
  c(a5_counts["International chain", "Elsewhere"], a5_counts["Other cafe type", "Elsewhere"]),
  as.integer(rowSums(a5_counts)),
  c(a5_share_intl, a5_share_other)
)

names(a5_count_table) <- c("Cafe type", "In Prague 1", "Elsewhere", "Total", "Percent in Prague 1")
save_table(a5_count_table, "A5_counts.tex")

a5_table <- data.frame(
  nrow(df_a5),
  a5_n_missing,
  a5_difference,
  a5_diff_lower,
  a5_diff_upper,
  unname(a5_interval$estimate),
  a5_interval$conf.int[1],
  a5_interval$conf.int[2],
  p_text(a5_test$p.value),
  a5_smallest_expected,
  a5_v
)

names(a5_table) <- c("N", "Missing", "Difference (pp)", "CI lower (pp)", "CI upper (pp)",
                     "Odds ratio", "OR lower", "OR upper", "p (one-tailed)",
                     "Smallest expected", "Cramer's V")

save_table(a5_table, "A5_table.tex")

## A5 FIGURE
a5_intl_ci  <- suppressWarnings(prop.test(a5_counts["International chain", "Prague 1"],
                                          sum(a5_counts["International chain", ]), correct = FALSE))
a5_other_ci <- suppressWarnings(prop.test(a5_counts["Other cafe type", "Prague 1"],
                                          sum(a5_counts["Other cafe type", ]), correct = FALSE))

a5_shares <- data.frame(
  chain = c("International chain", "Other cafe type"),
  share = c(a5_share_intl, a5_share_other),
  lower = 100 * c(a5_intl_ci$conf.int[1], a5_other_ci$conf.int[1]),
  upper = 100 * c(a5_intl_ci$conf.int[2], a5_other_ci$conf.int[2]),
  label = paste(c(a5_counts["International chain", "Prague 1"], a5_counts["Other cafe type", "Prague 1"]),
                "of", rowSums(a5_counts), "cafes")
)

a5_figure <- ggplot(a5_shares, aes(x = chain, y = share)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, linewidth = 1, color = "#2a78d6") +
  geom_point(size = 3, color = "#2a78d6") +
  geom_text(aes(label = label), hjust = -0.1) +
  labs(x = NULL, y = "Percent of cafes in Prague 1",
       title = "Share of cafes in Prague 1 (95% CI)") +
  theme_minimal()

ggsave(file.path(dir_fig, "A5_figure.pdf"), a5_figure, width = 6, height = 4)


# ===================================================================
# A6 -- Is the number of Google ratings correlated with the distance
#       from the Orloj, in which direction and how strongly?
# ===================================================================

# H0: The population correlation between number of ratings and distance equals 0.
# H1: The population correlation between number of ratings and distance is not 0.
# Tails and why: Two-tailed, because the question asks "in which direction".

# Method: Spearman's rho with a bootstrapped interval. D_table says
# non-parametric for the number of ratings.

df_a6 <- df %>% filter(!is.na(n_ratings))
a6_n_missing <- nrow(df) - nrow(df_a6)

a6_test <- cor.test(df_a6$n_ratings, df_a6$distance_km, method = "spearman", exact = FALSE)
a6_test

a6_rho <- unname(a6_test$estimate)

# resample whole cafes, so each cafe keeps its own pair of values
set.seed(123)
a6_boot <- replicate(n_boot, {
  i <- sample(nrow(df_a6), replace = TRUE)
  cor(df_a6$n_ratings[i], df_a6$distance_km[i], method = "spearman")
})

a6_lower <- unname(quantile(a6_boot, 0.025))
a6_upper <- unname(quantile(a6_boot, 0.975))

## A6 TABLE
a6_table <- data.frame(
  nrow(df_a6),
  a6_n_missing,
  a6_rho,
  a6_lower,
  a6_upper,
  a6_rho^2,
  unname(a6_test$statistic),
  p_text(a6_test$p.value)
)

names(a6_table) <- c("n", "Missing", "Spearman rho", "CI lower", "CI upper",
                     "rho squared", "S", "p")

save_table(a6_table, "A6_table.tex")

## A6 FIGURE
# A log scale cannot show 0, so 1 is added to every count. Spearman uses ranks,
# so the plot is the only thing this changes.
a6_figure <- ggplot(df_a6, aes(x = distance_km, y = n_ratings + 1)) +
  geom_point(color = "grey40", alpha = 0.6) +
  scale_y_log10() +
  labs(x = "Distance from the Orloj (km)", y = "Number of Google ratings + 1 (log scale)",
       title = paste0("Spearman rho = ", round(a6_rho, 2),
                      ", 95% CI [", round(a6_lower, 2), ", ", round(a6_upper, 2), "]")) +
  theme_minimal()

ggsave(file.path(dir_fig, "A6_figure.pdf"), a6_figure, width = 6, height = 4)


# ===================================================================
# A7 -- By how many ratings does the number of Google ratings change
#       per km from the Orloj?
# ===================================================================

# H0: The slope of number of ratings on distance equals 0 ratings per km.
# H1: The slope of number of ratings on distance is not 0.
# Tails and why: Two-tailed, because "by how many ... change" names no direction.

# Method: simple linear regression with robust (HC3) standard errors.
# DECISION: A7, A8 and A9 all use robust standard errors. The Breusch-Pagan
# test rejects constant variance in A8, and using one method for all three
# keeps their standard errors comparable.

df_a7 <- df %>% filter(!is.na(n_ratings))

a7_fit <- lm(n_ratings ~ distance_km, data = df_a7)
summary(a7_fit)

## A7 ASSUMPTION CHECKS
save_diagnostics(a7_fit, "A7_diagnostics.pdf")
bptest(a7_fit)

a7_cooks <- cooks.distance(a7_fit)
a7_influential <- influential_cafes(df_a7, a7_cooks)
a7_influential
save_table(a7_influential, "A7_influential.tex")

# Does the slope depend on the influential cafes?
a7_trimmed <- lm(n_ratings ~ distance_km, data = df_a7[a7_cooks <= 4 / nrow(df_a7), ])

a7_sensitivity <- data.frame(
  c("All cafes", "Without influential cafes"),
  c(coef(a7_fit)["distance_km"], coef(a7_trimmed)["distance_km"]),
  c(length(residuals(a7_fit)), length(residuals(a7_trimmed)))
)

names(a7_sensitivity) <- c("Sample", "Slope (ratings per km)", "n")
save_table(a7_sensitivity, "A7_sensitivity.tex")

## A7 TABLE
a7_table <- regression_table(a7_fit, "Simple")
save_table(a7_table, "A7_table.tex")

## A7 FIGURE
a7_figure <- ggplot(df_a7, aes(x = distance_km, y = n_ratings)) +
  geom_point(color = "grey40", alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#2a78d6") +
  labs(x = "Distance from the Orloj (km)", y = "Number of Google ratings",
       title = "Regression line with 95% confidence band") +
  theme_minimal()

ggsave(file.path(dir_fig, "A7_figure.pdf"), a7_figure, width = 6, height = 4)


# ===================================================================
# A8 -- Does the distance-ratings link hold once you account for the
#       number of seats inside, and how much of it remains?
# ===================================================================

# H0: Holding seats inside constant, the slope of number of ratings on distance equals 0.
# H1: Holding seats inside constant, the slope of number of ratings on distance is not 0.
# Tails and why: Two-tailed, because "hold" names no direction.

# Method: multiple linear regression with robust (HC3) standard errors.

# both models use the same cafes, so their slopes can be compared
df_a8 <- df %>% filter(!is.na(n_ratings), !is.na(seats_inside))
a8_n_missing <- nrow(df) - nrow(df_a8)
a8_n_missing

a8_simple <- lm(n_ratings ~ distance_km, data = df_a8)
a8_fit    <- lm(n_ratings ~ distance_km + seats_inside, data = df_a8)
summary(a8_fit)

## A8 ASSUMPTION CHECKS
save_diagnostics(a8_fit, "A8_diagnostics.pdf")
bptest(a8_fit)

a8_vif <- vif(a8_fit)
a8_vif

a8_cooks <- cooks.distance(a8_fit)
a8_influential <- influential_cafes(df_a8, a8_cooks)
a8_influential
save_table(a8_influential, "A8_influential.tex")

a8_trimmed <- lm(n_ratings ~ distance_km + seats_inside, data = df_a8[a8_cooks <= 4 / nrow(df_a8), ])

a8_sensitivity <- data.frame(
  c("All cafes", "Without influential cafes"),
  c(coef(a8_fit)["distance_km"], coef(a8_trimmed)["distance_km"]),
  c(length(residuals(a8_fit)), length(residuals(a8_trimmed)))
)

names(a8_sensitivity) <- c("Sample", "Slope of distance (ratings per km)", "n")
save_table(a8_sensitivity, "A8_sensitivity.tex")

## A8 TABLES
a8_table <- bind_rows(
  regression_table(a8_simple, "Without seats"),
  regression_table(a8_fit, "With seats")
)

save_table(a8_table, "A8_table.tex")

# how much of the slope remains, and does adding seats improve the model?
a8_slope_simple <- unname(coef(a8_simple)["distance_km"])
a8_slope_full   <- unname(coef(a8_fit)["distance_km"])
a8_percent_left <- 100 * a8_slope_full / a8_slope_simple

a8_f_test <- anova(a8_simple, a8_fit)
a8_f_test

a8_comparison <- data.frame(
  a8_slope_simple,
  a8_slope_full,
  a8_percent_left,
  a8_f_test$F[2],
  as.integer(a8_f_test$Df[2]),
  as.integer(a8_f_test$Res.Df[2]),
  p_text(a8_f_test$`Pr(>F)`[2]),
  unname(a8_vif["distance_km"])
)

names(a8_comparison) <- c("Slope without seats", "Slope with seats", "Percent remaining",
                          "F", "df1", "df2", "p", "VIF")

save_table(a8_comparison, "A8_comparison.tex")

## A8 FIGURES
a8_slopes <- a8_table %>%
  filter(Term == "distance_km") %>%
  mutate(Model = factor(Model, levels = c("Without seats", "With seats")))

a8_figure <- ggplot(a8_slopes, aes(x = Model, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#e34948") +
  geom_errorbar(aes(ymin = `CI lower`, ymax = `CI upper`), width = 0.1, linewidth = 1, color = "#2a78d6") +
  geom_point(size = 3, color = "#2a78d6") +
  labs(x = NULL, y = "Slope of distance (ratings per km)",
       title = "Slope of distance with robust 95% CI, before and after adding seats") +
  theme_minimal()

ggsave(file.path(dir_fig, "A8_figure.pdf"), a8_figure, width = 6, height = 4)

# the data behind each slope, holding the other predictor constant
pdf(file.path(dir_fig, "A8_partial.pdf"))
avPlots(a8_fit)
invisible(dev.off())


# ===================================================================
# A9 -- Does the distance-ratings relationship differ between cafes
#       east and west of the Orloj, and how?
# ===================================================================

# H0: The slope of number of ratings on distance is the same east and west of the Orloj (interaction = 0).
# H1: The slope of number of ratings on distance differs between east and west (interaction is not 0).
# Tails and why: Two-tailed, because the question says "differ".

# Method: linear regression with an interaction and robust (HC3) standard errors.
# is_east is 0 for west, so distance_km is the slope for the west and
# distance_km:is_east is how much steeper the slope is in the east.

df_a9 <- df %>% filter(!is.na(n_ratings))

a9_fit <- lm(n_ratings ~ distance_km * is_east, data = df_a9)
summary(a9_fit)

## A9 ASSUMPTION CHECKS
save_diagnostics(a9_fit, "A9_diagnostics.pdf")
bptest(a9_fit)

a9_cooks <- cooks.distance(a9_fit)
a9_influential <- influential_cafes(df_a9, a9_cooks)
a9_influential
save_table(a9_influential, "A9_influential.tex")

a9_trimmed <- lm(n_ratings ~ distance_km * is_east, data = df_a9[a9_cooks <= 4 / nrow(df_a9), ])

a9_sensitivity <- data.frame(
  c("All cafes", "Without influential cafes"),
  c(coef(a9_fit)["distance_km:is_east"], coef(a9_trimmed)["distance_km:is_east"]),
  c(length(residuals(a9_fit)), length(residuals(a9_trimmed)))
)

names(a9_sensitivity) <- c("Sample", "Interaction (ratings per km)", "n")
save_table(a9_sensitivity, "A9_sensitivity.tex")

## A9 TABLES
a9_table <- regression_table(a9_fit, "Interaction")
save_table(a9_table, "A9_table.tex")

a9_slope_west <- unname(coef(a9_fit)["distance_km"])
a9_slope_east <- a9_slope_west + unname(coef(a9_fit)["distance_km:is_east"])

a9_slopes <- data.frame(
  c("West of the Orloj", "East of the Orloj"),
  c(a9_slope_west, a9_slope_east)
)

names(a9_slopes) <- c("Side", "Slope (ratings per km)")
save_table(a9_slopes, "A9_slopes.tex")

## A9 FIGURE
df_a9_plot <- df_a9 %>% mutate(side = ifelse(is_east == 1, "East", "West"))

a9_figure <- ggplot(df_a9_plot, aes(x = distance_km, y = n_ratings, color = side)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x) +
  scale_color_manual(values = c("East" = "#2a78d6", "West" = "#e34948")) +
  labs(x = "Distance from the Orloj (km)", y = "Number of Google ratings", color = NULL,
       title = "Regression lines with 95% confidence bands, by side of the Orloj") +
  theme_minimal()

ggsave(file.path(dir_fig, "A9_figure.pdf"), a9_figure, width = 6, height = 4)


# ===================================================================
# A10 -- Does the probability of an above-median number of ratings
#        change with each additional km from the Orloj?
# ===================================================================

# H0: The probability of an above-median number of ratings does not change with distance (coefficient = 0).
# H1: The probability of an above-median number of ratings changes with distance (coefficient is not 0).
# Tails and why: Two-tailed, because the question says "change" and names no direction.

# Method: logistic regression, reported as an average marginal effect.

df_a10 <- df %>% filter(!is.na(is_above_median))
a10_n_missing <- nrow(df) - nrow(df_a10)

a10_fit <- glm(is_above_median ~ distance_km, data = df_a10, family = binomial)
summary(a10_fit)

# the average change in probability per extra km
a10_ame <- avg_slopes(a10_fit, variables = "distance_km")
a10_ame

a10_log_odds <- summary(a10_fit)$coefficients["distance_km", "Estimate"]

## A10 TABLE
a10_table <- data.frame(
  nrow(df_a10),
  a10_n_missing,
  a10_log_odds,
  exp(a10_log_odds),
  100 * a10_ame$estimate,
  100 * a10_ame$conf.low,
  100 * a10_ame$conf.high,
  a10_ame$statistic,
  p_text(a10_ame$p.value)
)

names(a10_table) <- c("n", "Missing", "Log-odds per km", "Odds ratio per km",
                      "Marginal effect (pp per km)", "CI lower (pp)", "CI upper (pp)",
                      "z", "p")

save_table(a10_table, "A10_table.tex")

## A10 FIGURE
set.seed(1100)

a10_figure <- ggplot(df_a10, aes(x = distance_km, y = is_above_median)) +
  geom_jitter(height = 0.03, width = 0, color = "grey40", alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial),
              formula = y ~ x, color = "#2a78d6") +
  labs(x = "Distance from the Orloj (km)", y = "Probability of above-median ratings",
       title = "Logistic curve with 95% confidence band") +
  theme_minimal()

ggsave(file.path(dir_fig, "A10_figure.pdf"), a10_figure, width = 6, height = 4)


# ===================================================================
# SUMMARY TABLE OF ALL TEN ANALYSES (for the main text of the report)
# ===================================================================

# one line of the summary table
summary_row <- function(analysis, test, estimate, unit, lower, upper,
                        effect, effect_name, statistic, p, n) {
  data.frame(analysis, test, estimate, unit, lower, upper, effect, effect_name, statistic, p, n)
}

# the slope of interest from each regression table
a7_row <- a7_table %>% filter(Term == "distance_km")
a8_row <- a8_table %>% filter(Model == "With seats", Term == "distance_km")
a9_row <- a9_table %>% filter(Term == "distance_km:is_east")

summary_table <- bind_rows(
  summary_row("A1 Mean distance vs. 2 km", "One-sample t", a1_mean, "km",
              a1_lower, a1_upper, a1_d, "d",
              unname(a1_test$statistic), p_text(a1_test$p.value), a1_n),
  summary_row("A2 Roastery minus other", "Wilcoxon rank-sum (one-tailed)", unname(a2_interval$estimate), "km",
              a2_interval$conf.int[1], a2_interval$conf.int[2], NA, "",
              unname(a2_test$statistic), p_text(a2_test$p.value), nrow(df_a2)),
  summary_row("A3 East-west minus north-south", "Paired t (one-tailed)", a3_mean, "km",
              a3_lower, a3_upper, a3_dz, "dz",
              unname(a3_test$statistic), p_text(a3_test$p.value), a3_n),
  summary_row("A4 Distance by cafe type", "Kruskal-Wallis", NA, "",
              NA, NA, NA, "",
              unname(a4_test$statistic), p_text(a4_test$p.value), nrow(df_a4)),
  summary_row("A5 International chain minus other, in Prague 1", "Fisher exact (one-tailed)", a5_difference, "pp",
              a5_diff_lower, a5_diff_upper, a5_v, "V",
              NA, p_text(a5_test$p.value), nrow(df_a5)),
  summary_row("A6 Ratings and distance", "Spearman", a6_rho, "rho",
              a6_lower, a6_upper, a6_rho^2, "rho squared",
              unname(a6_test$statistic), p_text(a6_test$p.value), nrow(df_a6)),
  summary_row("A7 Ratings per km", "Regression", a7_row$Estimate, "ratings per km",
              a7_row$`CI lower`, a7_row$`CI upper`, a7_row$`R squared`, "R squared",
              a7_row$t, a7_row$p, a7_row$n),
  summary_row("A8 Ratings per km, given seats", "Multiple regression", a8_row$Estimate, "ratings per km",
              a8_row$`CI lower`, a8_row$`CI upper`, a8_row$`R squared`, "R squared",
              a8_row$t, a8_row$p, a8_row$n),
  summary_row("A9 East minus west slope", "Interaction", a9_row$Estimate, "ratings per km",
              a9_row$`CI lower`, a9_row$`CI upper`, a9_row$`R squared`, "R squared",
              a9_row$t, a9_row$p, a9_row$n),
  summary_row("A10 Probability of above-median ratings", "Logistic regression", 100 * a10_ame$estimate, "pp per km",
              100 * a10_ame$conf.low, 100 * a10_ame$conf.high, NA, "",
              a10_ame$statistic, p_text(a10_ame$p.value), nrow(df_a10))
)

summary_table

names(summary_table) <- c("Analysis", "Test", "Estimate", "Unit", "CI lower", "CI upper",
                          "Effect size", "Effect type", "Statistic", "p", "n")

save_table(summary_table, "summary.tex")
