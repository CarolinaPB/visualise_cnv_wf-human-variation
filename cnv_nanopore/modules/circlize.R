# ===================== UI Module =====================
mod_circos_circlize_ui <- function(id) {
    ns <- NS(id)
    tagList(
        plotOutput(ns("circosPlot"), height = "800px"),
        verbatimTextOutput(ns("sv_summary")),
        DT::DTOutput(ns("sv_table"))
    )
}

# ===================== Helper Functions =====================
init_circos_default <- function(start_degree = 90) {
    par(mar = c(0, 0, 2, 0))  # leave top margin for sample name
    circos.par(
        "start.degree" = start_degree,
        "track.height" = 0.15,
        "cell.padding" = c(0, 0),
        "track.margin" = c(0.02, 0.02),  # more gap between ring and labels
        "canvas.xlim" = c(-1, 1),
        "canvas.ylim" = c(-1, 1)
    )
    
    circos.initializeWithIdeogram(
        species = "hg38",
        plotType = c("ideogram", "labels"),
        labels.cex = 1.2
    )
}



sv_colors <- c(
    BND = "#1f78b4", DEL = "#e31a1c", DUP = "#33a02c",
    INV = "#ff7f00", INS = "maroon2"
)

plot_sv_links <- function(sv_df, sample_name, alpha_pass = 0.7, alpha_fail = 0.2) {
    if (nrow(sv_df) == 0) return(NULL)
    
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
        
        link_colors <- ifelse(
            subset_df$pass_filter,
            fade_color(sv_colors[sv_type], alpha = alpha_pass),
            fade_color(sv_colors[sv_type], alpha = alpha_fail)
        )
        
        circos.genomicLink(region1, region2, col = link_colors)
    }
    
    # Add legend with xpd = TRUE to allow drawing outside plotting region
    legend(
        x = -1,          # leftmost x in plot coordinates
        y = 1.1,         # slightly above top of circle
        legend = names(sv_colors),
        fill = sv_colors,
        border = NA,
        bty = "n",
        cex = 1.2,
        pt.cex = 1.5,
        title = "SV Type",
        title.cex = 1.2,
        xpd = TRUE
    )
    
    title(main = sample_name, cex.main = 1.2)
}


plot_cnv_track <- function(cnv_df, bed_df, max_cov) {
    cnv_cols <- c(
        DEL = adjustcolor("red", alpha.f = 0.7),
        DUP = adjustcolor("green", alpha.f = 0.7)
    )
    
    cnv_df$col <- cnv_cols[cnv_df$svtype]
    cnv_df$col[is.na(cnv_df$col)] <- adjustcolor("grey", alpha.f = 0.3)
    cnv_by_chr <- split(cnv_df, cnv_df$chr)
    
    circos.genomicTrack(
        bed_df[, c("chr", "start", "end", "coverage")],
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

# ===================== Server Module =====================
mod_circos_circlize_server <- function(id, inputs, plots_res) {
    moduleServer(id, function(input, output, session) {
        
        # ----------------- Reactive: Filtered SVs -----------------
        filtered_svs <- reactive({
            req(inputs$data_list())
            sv_df <- inputs$data_list()$sv
            sample_name <- inputs$sample()
            if (is.null(sv_df) || nrow(sv_df) == 0) return(NULL)
            
            tumor_DV <- sym(paste0(sample_name, "_T.DV"))
            tumor_VAF <- sym(paste0(sample_name, "_T.VAF"))
            normal_DV <- sym(paste0(sample_name, "_N.DV"))
            normal_DR <- sym(paste0(sample_name, "_N.DR"))
            
            sv_df <- sv_df %>% filter(!is.na(SVTYPE))
            
            if (inputs$filter_tumor()) {
                if (inputs$filter_normal()) {
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
            } else {
                sv_df$pass_filter <- TRUE
            }
            
            sv_df
        })
        
        # ----------------- Reactive: CNV Cache -----------------
        cnv_cache <- reactive({
            req(inputs$data_list())
            cnv_df <- inputs$data_list()$cnv
            bed_df <- inputs$data_list()$bed
            
            if (is.null(cnv_df) || is.null(bed_df)) return(NULL)
            
            bed_binned <- bed_df %>%
                group_by(chr, bin = floor(start/1000)) %>%
                summarize(
                    coverage = mean(coverage),
                    start = min(start),
                    end = max(end),
                    .groups = "drop"
                )
            
            list(
                cnv_df = cnv_df,
                bed_binned = bed_binned,
                max_cov = inputs$max_coverage()
            )
        })
        
        # ----------------- SV Summary -----------------
        sv_counts <- reactive({
            sv_df <- filtered_svs()
            if (is.null(sv_df)) return(list(pass = 0, fail = 0))
            list(
                pass = sum(sv_df$pass_filter, na.rm = TRUE),
                fail = sum(!sv_df$pass_filter, na.rm = TRUE)
            )
        })
        
        output$sv_summary <- renderText({
            counts <- sv_counts()
            paste0("SVs passing filter: ", counts$pass, "\n",
                   "SVs failing filter: ", counts$fail)
        })
        
        # ----------------- SV Table -----------------
        output$sv_table <- DT::renderDT({
            sv_df <- filtered_svs() %>% filter(pass_filter == 1)
            if (is.null(sv_df)) return(NULL)
            
            DT::datatable(
                sv_df %>% dplyr::select(
                    CHROM1, POS1, CHROM2, POS2, SVTYPE,
                    Annotation, Annotation_Impact,
                    Gene_Name, Gene_ID,
                    starts_with(inputs$sample())
                ),
                options = list(pageLength = 10)
            )
        })
        
        # ----------------- Circos Plot -----------------
        output$circosPlot <- renderPlot({
            req(plots_res$main_plots_ready())
            dat_list <- inputs$data_list()
            sample_name <- inputs$sample()
            sv_df <- filtered_svs()
            
            withProgress(message = paste("Generating Circos plot for", sample_name), value = 0, {
                
                # Step 1: Initialize Circos
                incProgress(0.1, detail = "Initializing Circos...")
                circos.clear()
                init_circos_default()
                
                # Step 2: CNV + coverage track
                if (inputs$circos_cnv()) {
                    incProgress(0.4, detail = "Plotting CNV and coverage tracks...")
                    cnv_data <- cnv_cache()
                    if (!is.null(cnv_data)) {
                        plot_cnv_track(cnv_data$cnv_df, cnv_data$bed_binned, cnv_data$max_cov)
                    }
                }
                
                # Step 3: SV links
                if (inputs$circos_sv() && !is.null(sv_df) && nrow(sv_df) > 0) {
                    incProgress(0.7, detail = "Plotting SV links...")
                    plot_sv_links(sv_df, sample_name)
                }
                
                # Step 4: Finalize
                incProgress(1, detail = "Finalizing Circos plot...")
            })
        })
        
    })
}
