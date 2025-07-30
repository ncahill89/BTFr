###########
# run_modern took 5 minutes and 37 seconds
# run_core took 8 minutes and 21 seconds
# run_valid took 13 minutes and 52 seconds
# tested on 16-core CPU laptop
###########


# Required packages
required_packages <- c("parallel", "doParallel", "ggplot2")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
}

# Load libraries
library(BTFr)
library(parallel)
library(doParallel)
library(ggplot2)
library(dplyr)

# Set up parallel processing
n_cores <- detectCores() - 1 # Leave one core free
if (n_cores < 1) n_cores <- 1

# Create and register parallel cluster
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Record the start time for run_modern
start_time_run_modern <- Sys.time()

# Running the modern calibration model
modern_mod <- run_modern(
  # modern_elevation,
  # modern_species,
  parallel = TRUE,
  n_cores = n_cores
)

# Save Model Output
saveRDS(modern_mod, file = "MPmodern_test.rds", compress = TRUE)

# Record the end time
end_time_run_modern <- Sys.time()

# Calculate the duration
duration_run_modern <- difftime(
  end_time_run_modern, start_time_run_modern,
  units = "secs"
)

# Record the start time for run_core
start_time_run_core <- Sys.time()

# With prior
core_mod_prior <- run_core(
  modern_mod,
  parallel = TRUE,
  n_cores = n_cores
)
saveRDS(core_mod_prior, "MPcore_test.rds")

# Record the end time
end_time_run_core <- Sys.time()

# Calculate the duration
duration_run_core <- difftime(
  end_time_run_core, start_time_run_core,
  units = "secs"
)

# Record the start time for run_valid
start_time_run_valid <- Sys.time()


# Parallel cross-validation
valid_run_10fold <- run_valid(
  n_folds = 10,
  parallel = TRUE,
  n_cores = n_cores
)

saveRDS(valid_run_10fold, file = "MPvalid_test.rds")

# Record the end time
end_time_run_valid <- Sys.time()

# Calculate the duration
duration_run_valid <- difftime(
  end_time_run_valid, start_time_run_valid,
  units = "secs"
)


# Calculate summary statistics
valid_fold_summary <- valid_run_10fold %>%
  summarise(
    coverage = sum(lower < True & True < upper) * 100 / n(),
    RMSE = sqrt(mean((True - SWLI)^2))
  ) %>%
  round(2)

print(valid_fold_summary)

# Convert duration to minutes and seconds
duration_run_modern_minutes <- as.numeric(duration_run_modern) %/% 60
duration_run_modern_seconds <- round(as.numeric(duration_run_modern) %% 60)
duration_run_core_minutes <- as.numeric(duration_run_core) %/% 60
duration_run_core_seconds <- round(as.numeric(duration_run_core) %% 60)
duration_run_valid_minutes <- as.numeric(duration_run_valid) %/% 60
duration_run_valid_seconds <- round(as.numeric(duration_run_valid) %% 60)

# Output the duration
cat("run_modern took", duration_run_modern_minutes, "minutes and", duration_run_modern_seconds, "seconds\n")
cat("run_core took", duration_run_core_minutes, "minutes and", duration_run_core_seconds, "seconds\n")
cat("run_valid took", duration_run_valid_minutes, "minutes and", duration_run_valid_seconds, "seconds\n")

# Clean up
stopCluster(cl)
gc()
