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
    
  # initialise empty so that SNV gene labels can be plotted outside the ideogram
    circos.initializeWithIdeogram(plotType = NULL, species = "hg38")
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
        DEL = adjustcolor("#DC143C", alpha.f = 0.7),
        DUP = adjustcolor("#00FF7F", alpha.f = 0.7)
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

plot_snv_gene_labels <- function(snv_df){
  snv_df <- snv_df %>%  mutate(
    mid = floor((position + end_position) / 2)
  ) %>% 
    dplyr::select(chr = chromosome, start = position, end = end_position,
                  value = wf_somatic_snv_allele_fraction,
                  gene_symbol, 
                  mid)
  label_input <- snv_df %>%
    transmute(chr = chr, start = mid, end = mid, label = gene_symbol)
  
  circos.genomicLabels(
    label_input,
    labels.column = 4,
    side = "outside",                 # place labels outside ideogram
    labels_height = mm_h(2),          # distance of label text from circle
    connection_height = mm_h(2),      # length of connector line
    cex = 1.2,
    padding = mm_h(2),                # extra padding, helps readability
  )
}

plot_snv_track <- function(snv_df, track_height = 0.08, alpha = 0.5) {
  if (nrow(snv_df) == 0) return(NULL)
  
  # --- Define activity color mapping once ---
  activity_colors <- c(
    loss   = "#e31a1c",
    gain   = "#33a02c",
    amp    = "#ff7f00",
    normal = "#1f78b4",
    other  = "grey70"
  )
  
  snv_df <- snv_df %>%
    dplyr::mutate(
      activity_color = activity_colors[wf_somatic_snv_inferred_activity] %||% activity_colors["other"]
    ) %>%
    dplyr::select(
      chr = chromosome, start = position, end = end_position,
      value = wf_somatic_snv_allele_fraction,
      activity_color, gene_symbol
    )
  
  # --- Draw SNV points ---
  circos.genomicTrack(
    snv_df,
    ylim = c(0, 100),
    track.height = track_height,
    panel.fun = function(region, value, ...) {
      sec <- CELL_META$sector.index
      sec_idx <- which(snv_df$chr == sec)
      circos.genomicPoints(
        region, value,
        col = adjustcolor(snv_df$activity_color[sec_idx], alpha.f = alpha),
        pch = 16
      )
    }
  )
  
  # --- Draw legend using the same color mapping ---
  legend(
    x = 1.1, y = 1.1,
    legend = names(activity_colors),
    fill = activity_colors,
    border = NA,
    bty = "n",
    cex = 1.5,
    pt.cex = 1.5,
    title = "SNV Activity",
    title.cex = 1.2,
    xpd = TRUE
  )
  

}

plot_chr_labels <- function() {
  sectors <- get.all.sector.index()
  
  # Create a dummy track for chromosome IDs
  circos.track(
    ylim = c(0, 1),  # dummy y-axis
    track.height = 0.05,  # thin track
    panel.fun = function(x, y) {
      sec <- CELL_META$sector.index
      x_center <- mean(CELL_META$xlim)
      circos.text(
        x = x_center,
        y = 0.5,  # slightly above track
        labels = sub("^chr", "", sec),
        facing = "bending.inside",
        niceFacing = TRUE,
        cex = 1.2,
        adj = c(0.5, 0)
      )
    },
    bg.border = NA
  )
}









# ===================== Server Module =====================

mod_circos_circlize_server <- function(id, inputs, plots_res) {
  moduleServer(id, function(input, output, session) {

    # ----------------- Reactive: Filtered SVs -----------------
    sv_data_ready <- reactive({
      req(inputs$data_list())
      sv_df <- inputs$data_list()$sv
      sample_name <- inputs$sample()
      if (is.null(sv_df) || nrow(sv_df) == 0) return(NULL)

      tumor_DV <- sym(paste0(sample_name, "_T.DV"))
      tumor_VAF <- sym(paste0(sample_name, "_T.VAF"))
      normal_DV <- sym(paste0(sample_name, "_N.DV"))
      normal_DR <- sym(paste0(sample_name, "_N.DR"))

      sv_df <- sv_df %>%
        filter(!is.na(SVTYPE)) %>%
        mutate(
          CHROM1 = factor(CHROM1),
          CHROM2 = factor(CHROM2),
          SVTYPE = factor(SVTYPE)
        )

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

    
    cnv_data_ready <- reactive({
      req(inputs$data_list())
      cnv_df <- inputs$data_list()$cnv
      bed_df <- inputs$data_list()$bed
      if (is.null(cnv_df) || is.null(bed_df)) return(NULL)

      bed_binned <- bed_df %>%
        group_by(chr, bin = floor(start / 1000)) %>%
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
      sv_df <- sv_data_ready()
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
      sv_df <- sv_data_ready()
      if (is.null(sv_df)) return(NULL)
      sv_df <- sv_df %>% filter(pass_filter == 1)

      DT::datatable(
        sv_df %>% dplyr::select(
          CHROM1, POS1, CHROM2, POS2, SVTYPE, SVLEN,
          Annotation, Annotation_Impact,
          Gene_Name, Gene_ID,
          starts_with(inputs$sample())
        ),
        extensions = 'Buttons',
        filter = "top",
        options = list(
          pageLength = 10,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', "colvis"),
          searching = TRUE
        )
      )
    })

    # ----------------- Circos Plot -----------------
    output$circosPlot <- renderPlot({
      req(plots_res$main_plots_ready())
      sv_df <- sv_data_ready()
      cnv_data <- cnv_data_ready()
      sample_name <- inputs$sample()
      snv_data <- inputs$data_list()$snv 
    

      withProgress(message = paste("Generating Circos plot for", sample_name), value = 0, {

        # Initialize Circos
        incProgress(0.1, detail = "Initializing Circos...")
        circos.clear()
        init_circos_default() # initialises empty plot so that SNV labels can be plotted outside the ideogram
        
        # Process SNV data once at the start
        if (!is.null(snv_data) && nrow(snv_data) > 0) {
          snv_data <- snv_data %>%
            janitor::clean_names() %>%
            mutate(chromosome = paste0("chr", chromosome)) %>%
            rename_with(~ sub(paste0("^", tolower(sample_name), "_"), "", .x))
        }
        
        # Plot SNV gene labels first if requested and data exists
        if (inputs$circos_snv_genes() && !is.null(snv_data) && nrow(snv_data) > 0) {
          plot_snv_gene_labels(snv_data)
        }
        plot_chr_labels()
        
        # Initialise ideogram (after gene labels)
        circos.genomicIdeogram(species = "hg38")
        
        # Plot SNV track if data exists and track is selected
        if (inputs$circos_snv() && !is.null(snv_data) && nrow(snv_data) > 0) {
          incProgress(0.3, detail = "Plotting SNV track...")
          plot_snv_track(snv_data)
        } else if (inputs$circos_snv()) {
          incProgress(0.3, detail = "No SNVs available for plotting")
        }
        
        

        # CNV + coverage track
        if (inputs$circos_cnv() && !is.null(cnv_data)) {
          incProgress(0.4, detail = "Plotting CNV and coverage tracks...")
          plot_cnv_track(cnv_data$cnv_df, cnv_data$bed_binned, cnv_data$max_cov)
        }

        # SV links
        if (inputs$circos_sv() && !is.null(sv_df) && nrow(sv_df) > 0) {
          incProgress(0.7, detail = "Plotting SV links...")
          plot_sv_links(sv_df, sample_name)
        }

        # Finalize
        incProgress(1, detail = "Done.")
      })
    })
  })
}
