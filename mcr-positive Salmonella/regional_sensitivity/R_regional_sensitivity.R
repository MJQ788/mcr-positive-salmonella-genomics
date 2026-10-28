# Input: 1,514 isolate-level metadata rows in TableS1_for_R.csv.
# Analytical scope: region-standardized Hybrid proportions and BioProject robustness.
# Requires BASE R only (no packages or internet).

input_file <- "TableS1_for_R.csv"
output_dir <- "R_sensitivity_output"
B <- 30000L
seed <- 20261013L

if (!file.exists(input_file)) {
  stop("Cannot find TableS1_for_R.csv. Set the R working directory to its folder.")
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

d <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE,
              fileEncoding = "UTF-8")
names(d) <- sub("^\ufeff", "", names(d))
need <- c("BioProject", "Continent", "Country", "MLST",
          "mcr subtype", "Isolation source", "Risk_group")
if (!all(need %in% names(d))) stop("Required columns are missing.")
if (anyNA(d[, need])) stop("Missing value in a required column.")

for (col in need) d[[col]] <- trimws(as.character(d[[col]]))
allowed <- c("Human", "Chicken", "Environment", "Egg")
if (!all(d[["Isolation source"]] %in% allowed)) {
  stop("Unexpected Isolation source label.")
}
d$Group <- ifelse(d[["Isolation source"]] == "Human",
                  "Human", "Chicken_related")
d$Hybrid <- as.integer(d[["Risk_group"]] == "Hybrid")

# Input checks for the supplied dataset.
stopifnot(nrow(d) == 1514L,
          sum(d$Group == "Human") == 879L,
          sum(d$Group == "Chicken_related") == 635L,
          sum(d$Hybrid) == 223L)

regions <- c("Asia", "North America", "Europe")
groups <- c("Human", "Chicken_related")
strata <- data.frame(
  Region = rep(regions, each = 2),
  Group = rep(groups, times = length(regions)),
  stringsAsFactors = FALSE
)

# Each of the six source-region cells is fixed in the same order.
cell_key <- paste(strata$Region, strata$Group, sep = "|")
cell_id <- match(paste(d$Continent, d$Group, sep = "|"), cell_key)
within <- !is.na(cell_id)
n <- tabulate(cell_id[within], nbins = 6)
events <- tabulate(cell_id[within & d$Hybrid == 1L], nbins = 6)
stopifnot(sum(n) == 1442L, sum(events) == 202L, all(n > 0L))

regional_n <- n[c(1, 3, 5)] + n[c(2, 4, 6)]
pooled_w <- regional_n / sum(regional_n)

standardize <- function(nn, ee, ww) {
  stopifnot(length(nn) == length(ee),
            length(nn) == 2L * length(ww))
  if (any(nn <= 0L)) {
    return(c(Human = NA_real_, Chicken_related = NA_real_,
             Difference_pp = NA_real_))
  }
  rate <- ee / nn
  human <- sum(ww * rate[seq(1, length(nn), by = 2)])
  chicken <- sum(ww * rate[seq(2, length(nn), by = 2)])
  c(Human = human, Chicken_related = chicken,
    Difference_pp = 100 * (human - chicken))
}

# Test the mathematics, including the empty-stratum boundary case.
toy <- standardize(c(10, 10, 20, 20, 5, 5),
                   c(5, 1, 10, 2, 4, 1), rep(1 / 3, 3))
stopifnot(abs(unname(toy["Difference_pp"]) -
              100 * mean(c(0.4, 0.4, 0.6))) < 1e-10)
stopifnot(is.na(standardize(c(10, 0, 20, 20, 5, 5),
                            c(5, 0, 10, 2, 4, 1),
                            rep(1 / 3, 3))["Difference_pp"]))

std <- standardize(n, events, pooled_w)
full_human <- sum(d$Hybrid[d$Group == "Human"]) /
              sum(d$Group == "Human")
full_chicken <- sum(d$Hybrid[d$Group == "Chicken_related"]) /
                sum(d$Group == "Chicken_related")
full_diff <- 100 * (full_human - full_chicken)
stopifnot(abs(unname(std["Difference_pp"]) - 11.8976944) < 0.0001)

regional_table <- data.frame(
  Region = regions,
  Human_N = n[c(1, 3, 5)], Human_Hybrid = events[c(1, 3, 5)],
  Chicken_related_N = n[c(2, 4, 6)],
  Chicken_related_Hybrid = events[c(2, 4, 6)],
  Pooled_weight = pooled_w,
  Human_rate = events[c(1, 3, 5)] / n[c(1, 3, 5)],
  Chicken_related_rate = events[c(2, 4, 6)] / n[c(2, 4, 6)]
)

# Alternative standards for the same observed data.
equal_w <- rep(1 / 3, 3)
two_w <- regional_n[1:2] / sum(regional_n[1:2])
equal_std <- standardize(n, events, equal_w)
two_std <- standardize(n[1:4], events[1:4], two_w)

# Preserve all 129 projects for project-level resampling.
projects <- sort(unique(d$BioProject))
P <- length(projects)
stopifnot(P == 129L)
project_id <- match(d$BioProject, projects)
pn <- matrix(0L, nrow = P, ncol = 6)
pe <- matrix(0L, nrow = P, ncol = 6)
for (j in seq_len(6)) {
  selected <- !is.na(cell_id) & cell_id == j
  pn[, j] <- tabulate(project_id[selected], nbins = P)
  pe[, j] <- tabulate(project_id[selected & d$Hybrid == 1L], nbins = P)
}
stopifnot(identical(as.integer(colSums(pn)), as.integer(n)),
          identical(as.integer(colSums(pe)), as.integer(events)))

# Leave one entire BioProject out, keeping original pooled weights fixed.
loo <- do.call(rbind, lapply(seq_len(P), function(i) {
  est <- standardize(n - pn[i, ], events - pe[i, ], pooled_w)
  data.frame(
    BioProject = projects[i],
    Removed_N_in_three_regions = sum(pn[i, ]),
    Removed_Hybrid_in_three_regions = sum(pe[i, ]),
    Standardized_Human = unname(est["Human"]),
    Standardized_Chicken_related = unname(est["Chicken_related"]),
    Difference_pp = unname(est["Difference_pp"]),
    stringsAsFactors = FALSE
  )
}))
stopifnot(all(is.finite(loo$Difference_pp)))

# Cluster bootstrap: sample 129 whole BioProjects with replacement.
# Empty source-region cells cannot be standardized and are omitted.
set.seed(seed)
boot_project <- rep(NA_real_, B)
for (b in seq_len(B)) {
  times <- tabulate(sample.int(P, size = P, replace = TRUE), nbins = P)
  bn <- as.vector(crossprod(times, pn))
  be <- as.vector(crossprod(times, pe))
  boot_project[b] <- standardize(bn, be, pooled_w)["Difference_pp"]
}
project_valid <- boot_project[is.finite(boot_project)]
if (length(project_valid) < B * 0.90) {
  stop("Too many invalid project bootstrap replicates.")
}
project_ci <- unname(quantile(project_valid, c(0.025, 0.975)))

# Independent-isolate parametric binomial bootstrap.
set.seed(seed + 1L)
boot_isolate <- replicate(B, {
  be <- rbinom(length(n), size = n, prob = events / n)
  standardize(n, be, pooled_w)["Difference_pp"]
})
isolate_ci <- unname(quantile(boot_isolate, c(0.025, 0.975)))

# China ST34 mcr-1 subset.
china <- d[d$Country == "China" &
           d$MLST == "34" & d[["mcr subtype"]] == "mcr-1", ]
china_h_n <- sum(china$Group == "Human")
china_h_h <- sum(china$Group == "Human" & china$Hybrid == 1L)
china_c_n <- sum(china$Group == "Chicken_related")
china_c_h <- sum(china$Group == "Chicken_related" & china$Hybrid == 1L)
stopifnot(china_h_n == 351L, china_h_h == 125L,
          china_c_n == 13L, china_c_h == 4L)

analysis_summary <- data.frame(
  Measure = c("All isolates: human Hybrid rate",
              "All isolates: chicken-related Hybrid rate",
              "All isolates: difference (pp)",
              "Three-region standardized human rate",
              "Three-region standardized chicken-related rate",
              "Three-region standardized difference (pp)",
              "Equal three-region weights: difference (pp)",
              "Asia + North America: difference (pp)",
              "BioProject bootstrap 95% CI lower (pp)",
              "BioProject bootstrap 95% CI upper (pp)",
              "Isolate bootstrap 95% CI lower (pp)",
              "Isolate bootstrap 95% CI upper (pp)",
              "Leave-one-project-out minimum (pp)",
              "Leave-one-project-out maximum (pp)",
              "Valid project bootstrap replicates",
              "China ST34/mcr-1 human: N, Hybrid",
              "China ST34/mcr-1 chicken-related: N, Hybrid"),
  Value = c(full_human, full_chicken, full_diff,
            unname(std["Human"]), unname(std["Chicken_related"]),
            unname(std["Difference_pp"]),
            unname(equal_std["Difference_pp"]),
            unname(two_std["Difference_pp"]),
            project_ci[1], project_ci[2],
            isolate_ci[1], isolate_ci[2],
            min(loo$Difference_pp), max(loo$Difference_pp),
            length(project_valid),
            paste(china_h_n, china_h_h, sep = ", "),
            paste(china_c_n, china_c_h, sep = ", ")),
  stringsAsFactors = FALSE
)

write.csv(regional_table, file.path(output_dir, "regional_counts.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(loo, file.path(output_dir, "leave_one_BioProject_out.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(analysis_summary, file.path(output_dir, "analysis_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "R_session_info.txt"))

cat("\nR analysis completed. Results in:", normalizePath(output_dir), "\n")
cat("All-isolate difference (pp):", round(full_diff, 4), "\n")
cat("Three-region standardized difference (pp):",
    round(unname(std["Difference_pp"]), 4), "\n")
cat("Project bootstrap 95% interval (pp):",
    paste(round(project_ci, 2), collapse = " to "), "\n")
cat("Leave-one-project-out range (pp):",
    paste(round(range(loo$Difference_pp), 2), collapse = " to "), "\n")
cat("Valid project bootstrap replicates:", length(project_valid), "/", B, "\n")

