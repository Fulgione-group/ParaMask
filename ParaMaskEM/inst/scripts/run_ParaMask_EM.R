#!/usr/bin/env Rscript

# Suppress startup messages when loading the package
suppressPackageStartupMessages(library(ParaMaskEM))

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Help message
if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
  cat("Usage: Rscript run_paramaskem_pipeline.R --het <input.het> [options]\n")
  cat("Options:\n")
  cat("  --het         Path to .het file (REQUIRED)\n")
  cat("  --chrom       Chromosome index (default: 0)\n")
  cat("  --startline   Start line in het file (default: 2)\n")
  cat("  --endline     End line (0 = all) (default: 0)\n")
  cat("  --missingness Missingness filter (default: 0.1)\n")
  cat("  --outdir      Output directory (default: current dir)\n")
  cat("  --ID          Sample ID prefix for outputs\n")
  cat("  --verbose     Verbose mode\n")
  cat("  --tolerance   EM convergence tolerance (default: 0.001)\n")
  cat("  --num_it      Max EM iterations (default: 100)\n")
  cat("  --dist_em_rep Repetitions for distance cutoff (default: 1000)\n")
  cat("  --cdist       Use conservative distance estimation\n")
  cat("  --boundary    Lower,Upper bounds for EM fit, e.g. \"-1000,1000\"\n")
  cat("  --noRRD       Disable read ratio deviation classification\n")
  quit(status = 0)
}

# Set defaults
hetpath <- NULL
outpath <- getwd()
ID <- ""
missingness <- 0.1
startline <- 2
endline <- 0
verbose <- FALSE
tolerance <- 0.001
num_iterations <- 100
dist_em_rep <- 1000
chr <- 0
cdist <- FALSE
min_dist <- 50
max_dist <- 5000
boundfit <- FALSE
lboundary <- -1000
uboundary <- 1000
useRRD <- TRUE

#print args
for(i in 1:length(args)){
  print(args[i])
}

# Simple CLI parser
i <- 1
while (i <= length(args)) {
  if (args[i] %in% c("--het", "-h")) {
    hetpath <- args[i + 1]; i <- i + 1
  } else if (args[i] %in% c("--chrom", "-c")) {
    chr <- as.numeric(args[i + 1]); i <- i + 1
  } else if (args[i] %in% c("--startline", "-s")) {
    startline <- as.numeric(args[i + 1]); i <- i + 1
  } else if (args[i] %in% c("--endline", "-e")) {
    endline <- as.numeric(args[i + 1]); i <- i + 1
  } else if (args[i] %in% c("--missingness", "-m")) {
    missingness <- as.numeric(args[i + 1]); i <- i + 1
  } else if (args[i] %in% c("--outdir", "-o")) {
    outpath <- args[i + 1]; i <- i + 1
  } else if (args[i] %in% c("--ID", "-id")) {
    ID <- args[i + 1]; i <- i + 1
  } else if (args[i] %in% c("--verbose", "-v")) {
    verbose <- TRUE
  } else if (args[i] %in% c("--tolerance", "-t")) {
    tolerance <- as.numeric(args[i + 1]); i <- i + 1
  } else if (args[i] %in% c("--num_it", "-ni")) {
    num_iterations <- as.numeric(args[i + 1]); i <- i + 1
  } else if (args[i] %in% c("--dist_em_rep", "-dr")) {
    dist_em_rep <- as.numeric(args[i + 1]); i <- i + 1
  } else if (args[i] %in% c("--boundary", "-b")) {
    boundary <- args[i + 1]
    lboundary <- as.numeric(strsplit(boundary, ",")[[1]][1])
    uboundary <- as.numeric(strsplit(boundary, ",")[[1]][2])
    boundfit <- TRUE
    i <- i + 1
  } else if (args[i] %in% c("--cdist", "-cd")) {
    cdist <- TRUE
  } else if (args[i] == "--noRRD") {
    useRRD <- FALSE
  }
  i <- i + 1
}

# Ensure output path ends with /
if (!endsWith(outpath, "/")) outpath <- paste0(outpath, "/")

if (is.null(hetpath)) stop("You must provide a --het path to the .het file.")

het <- ParaMaskEM::load_het_file(hetpath, verbose = TRUE)

em_input_data <- ParaMaskEM::prepare_em_input(het, subsample = 10000, verbose = verbose)

intial_em_fit <- ParaMaskEM::fit_initial_model(
  regData2 = em_input_data$regData2,
  weight1 = em_input_data$weight1,
  weight2 = em_input_data$weight2,
  boundfit = boundfit,
  lboundary = lboundary,
  uboundary = uboundary,
  verbose = verbose
)

tryCatch({
  ParaMaskEM::plot_EM_it(het, intial_em_fit$fit, em_input_data$weight1, em_input_data$rs, 0, outpath, ID)
}, error = function(e) message("Initial EM plot failed: ", conditionMessage(e)))

EM_results <- ParaMaskEM::run_em_loop(
  regData2 = em_input_data$regData2,
  predDF = em_input_data$predDF,
  fit = intial_em_fit$fit,
  coef_fit = intial_em_fit$coef_fit,
  offsetfit = intial_em_fit$offsetfit,
  em_input_data$weight1,
  em_input_data$weight2,
  maxiter = num_iterations,
  boundfit = boundfit,
  lboundary = lboundary,
  uboundary = uboundary,
  tolerance = tolerance,
  verbose = verbose
)

tryCatch({
  ParaMaskEM::plot_EM_it(het, EM_results$fit, EM_results$weight1, em_input_data$rs, EM_results$iteration, outpath, ID)
}, error = function(e) message("Final EM plot failed: ", conditionMessage(e)))

tryCatch({
  ParaMaskEM::plot_em_loglikelihood(EM_results$log_likelihood, outpath, ID)
}, error = function(e) message("Log-likelihood plot failed: ", conditionMessage(e)))

het <- ParaMaskEM::classify_snps_by_llr(het, EM_results$probs1, EM_results$probs2, cutoff = 0.99)

tryCatch({
  ParaMaskEM::plot_llr_classification(het, em_input_data$rs, cutoff = 0.99, outpath, ID)
}, error = function(e) message("LLR classification plot failed: ", conditionMessage(e)))

#use read ratio deviation
if(useRRD){
  tryCatch({
    results_classify_by_read_ratio_deviation <- classify_by_read_ratio_deviation(het=het, verbose = verbose)
    het <- results_classify_by_read_ratio_deviation$het
  }, error = function(e) {
    message("Classification with read ratio deviations failed: ", conditionMessage(e))
  })
  tryCatch({
    plot_read_ratio_deviation(het=het, em_input_data$rs, outpath = outpath, ID = ID, q1 = results_classify_by_read_ratio_deviation$q1, q2 = results_classify_by_read_ratio_deviation$q2)
  }, error = function(e) {
    message("Plotting read ratio deviations failed: ", conditionMessage(e))
  })
}

het$Position <- formatC(het$Position, format = "f", digits = 0)
write.table(het, file = paste0(outpath, ID, "_EMresults.het"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
het$Position <- as.numeric(het$Position)

if (sum(het$EM_class == 2, na.rm = TRUE) == 0) {
  warning("No sites classified as multicopy")
}

#estimate distance cutoff
tryCatch({
  dist_cutoff_results <- estimate_distance_cutoff(het = het, dist_em_rep = dist_em_rep, min_dist = min_dist, max_dist = max_dist, cdist = cdist, verbose = verbose)
  write.table(x = dist_cutoff_results$dist_cutoff,file =paste(outpath, ID,"_EMresults.dist", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
}, error = function(e) {
  message("Estimating distance cutoff failed: ", conditionMessage(e))
})

#plot distances
tryCatch({
  plot_distance_cutoff(dist_cutoff_results$distances, dist_cutoff_results$dist_cutoff, dist_cutoff_results$p1,dist_cutoff_results$p2, outpath = outpath, ID = ID)
}, error = function(e) {
  message("Plotting distance distributions failed: ", conditionMessage(e))
})

quit(status = 0)
