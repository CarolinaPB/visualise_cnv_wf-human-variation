mod_circos_circlize_ui <- function(id) {
    ns <- NS(id)
    tagList(
        plotOutput(ns("circosPlot"), height = "800px"),
        verbatimTextOutput(ns("sv_summary"))
    )
}

init_circos_default <- function(start_degree = 90) {
    par(mar = c(1, 1, 1, 1))
    circos.par("start.degree" = start_degree)
    circos.par("track.height" = 0.15)
    circos.par("cell.padding" = c(0.01, 0.01), "track.margin" = c(0.02, 0.02))
    circos.par("canvas.xlim" = c(-1.1, 1.1), "canvas.ylim" = c(-1.1, 1.1))
    
    circos.initializeWithIdeogram(
        species = "hg38",
        plotType = c("ideogram", "labels"),
        labels.cex = 0.8
    )
}

# Helper: fade color
fade_color <- function(color, alpha = 0.1) adjustcolor(color, alpha.f = alpha)
sv_colors <- c(BND = "#1f78b4", DEL = "#e31a1c", DUP = "#33a02c",
               INV = "#ff7f00", INS = "maroon2")

plot_sv_links <- function(sv_df, 
                          filter_tumor, 
                          min_tumor_DV,
                          min_tumor_VAF,
                          min_MAPQ,
                          filter_normal,
                          max_normal_DV,
                          min_normal_DR,
                          sample_name,
                          alpha_pass = 0.7, 
                          alpha_fail = 0.2) {
    if (nrow(sv_df) == 0) {
        message("No SV events found.")
        return(NULL)
    }
    
    sv_df <- sv_df %>% filter(!is.na(SVTYPE))
    
    # Apply filtering if requested
    tumor_DV   <- sym(paste0(sample_name, "_T.DV"))
    tumor_VAF  <- sym(paste0(sample_name, "_T.VAF"))
    normal_DV  <- sym(paste0(sample_name, "_N.DV"))
    normal_DR  <- sym(paste0(sample_name, "_N.DR"))
    if (filter_tumor) {
        if (filter_normal) {
            sv_df <- sv_df %>%
                mutate(
                    pass_filter = (!!tumor_DV) >= min_tumor_DV &
                        (!!tumor_VAF) >= min_tumor_VAF &
                        MAPQ >= min_MAPQ &
                        (!!normal_DV) <= max_normal_DV &
                        (!!normal_DR) >= min_normal_DR
                )
        } else {
            sv_df <- sv_df %>%
                mutate(
                    pass_filter = (!!tumor_DV) >= min_tumor_DV &
                        (!!tumor_VAF) >= min_tumor_VAF &
                        MAPQ >= min_MAPQ
                )
        }
        
        # Summary counts
        n_pass <- sum(sv_df$pass_filter, na.rm = TRUE)
        n_fail <- sum(!sv_df$pass_filter, na.rm = TRUE)
        message("Sample: ", sample_name)
        message("SVs passing filter: ", n_pass)
        message("SVs failing filter: ", n_fail)
    }
    
    # Plot links
    for (sv_type in unique(sv_df$SVTYPE)) {
        subset_df <- sv_df %>% filter(SVTYPE == sv_type)
        if (nrow(subset_df) == 0) next
        
        region1 <- data.frame(
            chr   = subset_df$CHROM1,
            start = subset_df$POS1,
            end   = ifelse(!is.na(subset_df$SVLEN), subset_df$POS1 + abs(subset_df$SVLEN), subset_df$POS1 + 1)
        )
        
        region2 <- data.frame(
            chr   = subset_df$CHROM2,
            start = subset_df$POS2,
            end   = ifelse(!is.na(subset_df$SVLEN), subset_df$POS2 + abs(subset_df$SVLEN), subset_df$POS2 + 1)
        )
        
        if (filter_tumor) {
            link_colors <- ifelse(
                subset_df$pass_filter,
                fade_color(sv_colors[sv_type], alpha = alpha_pass),
                fade_color(sv_colors[sv_type], alpha = alpha_fail)
            )
        } else {
            link_colors <- fade_color(sv_colors[sv_type], alpha = alpha_pass)
        }
        
        circos.genomicLink(region1, region2, col = link_colors)
    }
    
    # Legend
    legend(
        "topleft",
        legend = names(sv_colors),
        fill   = sv_colors,
        border = NA,
        bty    = "n",
        cex    = 0.8,
        title  = "SV Type"
    )
    
    # Title
    if (filter_tumor) {
        if (filter_normal) {
            title(main = paste(sample_name, "- Filtered (with normal)"), cex.main = 1.2)
        } else {
            title(main = paste(sample_name, "- Filtered (tumor only)"), cex.main = 1.2)
        }
    } else {
        title(main = sample_name, cex.main = 1.2)
    }
}

plot_cnv_track <- function(cnv_df, bed_df, max_cov) {
    cnv_cols <- c(DEL = adjustcolor("red", alpha.f = 0.7),
                  DUP = adjustcolor("green", alpha.f = 0.7))
    
    # precompute CNV colors
    cnv_df$col <- cnv_cols[cnv_df$svtype]
    cnv_df$col[is.na(cnv_df$col)] <- adjustcolor("grey", alpha.f = 0.3)
    
    # split by chromosome for fast lookup
    cnv_by_chr <- split(cnv_df, cnv_df$chr)
    
    circos.genomicTrack(
        bed_df[, c("chr","start","end","coverage")],
        ylim = c(0, max_cov),
        panel.fun = function(region, value, ...) {
            this_cnv <- cnv_by_chr[[CELL_META$sector.index]]
            if (!is.null(this_cnv) && nrow(this_cnv) > 0) {
                circos.rect(
                    xleft   = this_cnv$start,
                    ybottom = 0,
                    xright  = this_cnv$end,
                    ytop    = max_cov,
                    col     = this_cnv$col,
                    border  = NA
                )
            }
            value <- pmin(pmax(value, 0), max_cov)
            circos.genomicLines(region, value, col = "steelblue", lwd = 0.5)
            
        }
    )
}



mod_circos_circlize_server <- function(id, inputs, plots_res) {
    moduleServer(id, function(input, output, session) {
        
        # Reactive values to store SV counts
        sv_counts <- reactiveVal(list(pass = 0, fail = 0))
        
        # Cache CNV data to avoid repeated binning
        cached_cnv_data <- reactiveVal(NULL)
        
        output$circosPlot <- renderPlot({
            req(plots_res$main_plots_ready())
            dat_list <- inputs$data_list()
            sv_df <- dat_list$sv
            cnv_df <- dat_list$cnv
            bed_df <- dat_list$bed
            sample_name <- inputs$sample()
            
            withProgress(message = paste("Generating Circos plot for", sample_name), value = 0, {
                
                # Step 1: Initialize Circos
                incProgress(0.1, detail = "Initializing Circos...")
                circos.clear()
                init_circos_default()
                
                # Step 2: CNV + coverage track
                if(inputs$circos_cnv()) {
                    incProgress(0.4, detail = "Plotting CNV and coverage tracks...")
                    
                    if(is.null(cached_cnv_data())) {
                        # Bin coverage for speed
                        bed_binned <- bed_df %>%
                            group_by(chr, bin = floor(start/1000)) %>%
                            summarize(coverage = mean(coverage), start = min(start), end = max(end), .groups = "drop")
                        
                        cached_cnv_data(list(
                            cnv_df = cnv_df,
                            bed_binned = bed_binned,
                            max_cov = inputs$max_coverage()
                        ))
                    }
                    
                    # Draw CNV track
                    cnv_cache <- cached_cnv_data()
                    plot_cnv_track(cnv_cache$cnv_df, cnv_cache$bed_binned, cnv_cache$max_cov)
                }
                
                # Step 3: SV links
                if(inputs$circos_sv() && !is.null(sv_df) && nrow(sv_df) > 0) {
                    incProgress(0.7, detail = "Plotting SV links...")
                    
                    plot_sv_links(
                        sv_df,
                        filter_tumor = inputs$filter_tumor(),
                        filter_normal = inputs$filter_normal(),
                        min_tumor_DV = inputs$min_tumor_DV(),
                        min_tumor_VAF = inputs$min_tumor_VAF(),
                        min_MAPQ = inputs$min_MAPQ(),
                        max_normal_DV = inputs$max_normal_DV(),
                        min_normal_DR = inputs$min_normal_DR(),
                        sample_name = sample_name
                    )
                    
                    # Update reactive counts for SV summary
                    sv_df <- sv_df %>% filter(!is.na(SVTYPE))
                    tumor_DV <- sym(paste0(sample_name, "_T.DV"))
                    tumor_VAF <- sym(paste0(sample_name, "_T.VAF"))
                    normal_DV <- sym(paste0(sample_name, "_N.DV"))
                    normal_DR <- sym(paste0(sample_name, "_N.DR"))
                    
                    if(inputs$filter_tumor()) {
                        if(inputs$filter_normal()) {
                            sv_df <- sv_df %>%
                                mutate(pass_filter = (!!tumor_DV) >= inputs$min_tumor_DV() &
                                           (!!tumor_VAF) >= inputs$min_tumor_VAF() &
                                           MAPQ >= inputs$min_MAPQ() &
                                           (!!normal_DV) <= inputs$max_normal_DV() &
                                           (!!normal_DR) >= inputs$min_normal_DR())
                        } else {
                            sv_df <- sv_df %>%
                                mutate(pass_filter = (!!tumor_DV) >= inputs$min_tumor_DV() &
                                           (!!tumor_VAF) >= inputs$min_tumor_VAF() &
                                           MAPQ >= inputs$min_MAPQ())
                        }
                        
                        sv_counts(list(
                            pass = sum(sv_df$pass_filter, na.rm = TRUE),
                            fail = sum(!sv_df$pass_filter, na.rm = TRUE)
                        ))
                    }
                }
                
                # Step 4: Finalize
                incProgress(1, detail = "Finalizing Circos plot...")
            })
        })
        
        # Render the SV summary text
        output$sv_summary <- renderText({
            counts <- sv_counts()
            paste0("SVs passing filter: ", counts$pass, "\n",
                   "SVs failing filter: ", counts$fail)
        })
        
        # Reset CNV cache if sample or CNV settings change
        observeEvent(list(inputs$sample(), inputs$circos_cnv()), {
            cached_cnv_data(NULL)
        })
    })
}
