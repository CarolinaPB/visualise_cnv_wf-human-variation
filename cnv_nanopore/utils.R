get_samples <- function(dir) {
    files <- list.files(dir, pattern = "\\.regions\\.bed\\.gz$", full.names = FALSE)
    samples <- sub("\\..*$", "", files)
    unique(samples)
}

# adjust color opacity
fade_color <- function(color, alpha = 0.1) adjustcolor(color, alpha.f = alpha)

is_all_numeric <- function(x) all(!is.na(suppressWarnings(as.numeric(x))))
