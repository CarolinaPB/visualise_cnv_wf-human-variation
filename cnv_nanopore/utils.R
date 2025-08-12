get_samples <- function(dir) {
    files <- list.files(dir, pattern = "\\.regions\\.bed\\.gz$", full.names = FALSE)
    samples <- sub("\\..*$", "", files)
    unique(samples)
}
