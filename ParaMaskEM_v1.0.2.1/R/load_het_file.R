#' Load .het file for ParaMaskEM
#'
#' Reads a .het file with optional line range and chromosome filtering.
#'
#' @param hetpath Path to the .het file.
#' @param startline Line number to start reading from (default = 2).
#' @param endline Line number to stop reading (0 means read full file).
#' @param chr Optional chromosome to subset by (0 = all chromosomes).
#' @param missingness Maximum allowed missingness (default = 0.1).
#' @param verbose Whether to print progress messages.
#'
#' @return A filtered `data.frame` with heterozygosity stats.
#' @export
load_het_file <- function(hetpath, startline = 2, endline = 0, chr = 0, missingness = 0.1, verbose = FALSE) {
  if (verbose) {
    message("Reading het file...")
  }

  het <- tryCatch({
    if (endline == 0) {
      het <- data.table::fread(file = hetpath, header = TRUE)
    } else {
      header <- data.table::fread(file = hetpath, nrows = 1, header = FALSE)
      het <- data.table::fread(
        file = hetpath,
        skip = (startline - 1),
        nrows = (endline - startline + 1),
        header = FALSE
      )
      colnames(het) <- as.character(header)
    }

    if (chr != 0 && "Chromosome" %in% colnames(het)) {
      het <- het[het$Chromosome == chr, ]
    }

    het <- as.data.frame(het, stringsAsFactors = FALSE)

    # Convert all columns except 'Chromosome' to numeric
    num_cols <- setdiff(names(het), "Chromosome")
    het[, num_cols] <- lapply(het[, num_cols], function(x) {
      if (is.factor(x)) x <- as.character(x)
      suppressWarnings(as.numeric(x))
    })

    # Apply basic filter
    het <- het[het$Non.missing >= (1 - missingness) * max(het$Non.missing) & het$Minor.allele.freq > 0, ]


  }, warning = function(w) {
    stop("Warning while loading hetfile: ", conditionMessage(w))
  }, error = function(e) {
    stop("Error while loading hetfile: ", conditionMessage(e))
  })

  return(het)
}
