mod_plots_ui <- function(id) {
    ns <- NS(id)
    tagList(
        plotlyOutput(ns("facetPlot"), width = "100%", height = "100%"),
        br(),
        plotlyOutput(ns("chrPlot"))
    )
}

mod_plots_server <- function(id, inputs, config_section = "default") {
    moduleServer(id, function(input, output, session) {
        
        config <- config::get(config = config_section, 
                              file = here::here("cnv_nanopore/configs/config.yaml"))
        
        cytoband_file <- here(config$cytoband_file)
        fai_file <- here(config$fai_file)
        
        gie_colors <- c(
            gneg = "ghostwhite", gpos25 = "gray70", gpos50 = "gray50",
            gpos75 = "gray30", gpos100 = "black", gvar = "gray90",
            stalk = "brown", acen = "red"
        )
        # Load cytoband once outside renderPlotly
        cytoband <- read.table(cytoband_file, header = TRUE, sep = "\t",
                               col.names = c("chr", "start", "end", "name", "gieStain"),
                               stringsAsFactors = FALSE) %>%
            mutate(chr = sub("^chr", "", chr),
                   color = gie_colors[gieStain]) %>% as_tibble() 
        
        # Load fai data once
        fai_df <- read.table(fai_file, header = FALSE, sep = "\t",
                             col.names = c("chr", "length", "start", "value1", "value2"),
                             stringsAsFactors = FALSE) %>%
            mutate(chr = sub("^chr", "", chr),
                   middle = start + (length / 2)) %>%
            as_tibble() %>%
            dplyr::slice(1:24)
        
        
        output$facetPlot <- renderPlotly({
            req(inputs$sample())
            req(inputs$data_list())
            
            dat <- inputs$data_list()
            ordered_chrs <- c(paste0("chr", 1:22), "chrX", "chrY")
            
            bed <- dat$bed %>% mutate(chr = factor(chr, levels = ordered_chrs))
            cnv <- dat$cnv %>% mutate(chr = factor(chr, levels = ordered_chrs))
            
            # Calculate ymin/ymax for CNV rects from bed coverage for consistent vertical range
            ymin_rect <- min(bed$coverage, na.rm = TRUE) - 1
            ymax_rect <- max(bed$coverage, na.rm = TRUE) + 1
            
            p_facet <- ggplot() +
                geom_line(data = bed, aes(x = start, y = coverage), color = "steelblue") +
                geom_rect(
                    data = cnv,
                    aes(
                        xmin = start, xmax = end,
                        fill = svtype,
                        text = paste0(
                            "CNV: ", svtype,
                            "<br>CN: ", cn,
                            "<br>Chr: ", chr,
                            "<br>Start: ", start,
                            "<br>End: ", end
                        )
                    ),
                    ymin = ymin_rect,
                    ymax = ymax_rect,
                    alpha = 0.2
                ) +
                geom_point(
                    data = cnv,
                    aes(
                        x = start,
                        y = ymax_rect,
                        color = svtype,
                        text = paste0(
                            "CNV Start: ", start,
                            "<br>Chr: ", chr,
                            "<br>Type: ", svtype,
                            "<br>CN: ", cn
                        )
                    ),
                    shape = 17,
                    size = 2
                ) +
                scale_fill_manual(values = c(DEL = "red", DUP = "green")) +
                scale_color_manual(values = c(DEL = "red", DUP = "green")) +
                facet_wrap(~ chr, scales = "free", ncol = 4) +
                labs(
                    title = inputs$sample(),
                    x = "Genomic Position",
                    y = "Mean Coverage",
                    fill = "CNV Type",
                    color = "CNV Type"
                ) +
                theme_minimal() +
                theme(panel.grid = element_blank(),
                      axis.text.x = element_blank(),
                      axis.ticks.x = element_blank()
                )
            
            ggplotly(p_facet, tooltip = "text")
        })
        
        output$chrPlot <- renderPlotly({
            req(inputs$sample(), inputs$data_list(), inputs$chromosome())
            if (inputs$chromosome() == "All") return(NULL)
            
            dat <- inputs$data_list()
            
            bed_chr <- dat$bed %>%
                filter(chr == inputs$chromosome(),
                       start >= inputs$start_pos(),
                       end <= inputs$end_pos())
            
            cnv_chr <- dat$cnv %>%
                filter(chr == inputs$chromosome()) %>%
                filter((start >= inputs$start_pos() & start <= inputs$end_pos()) |
                           (end   >= inputs$start_pos() & end   <= inputs$end_pos()) |
                           (start <= inputs$start_pos() & end   >= inputs$end_pos()))
            
            cyto_min <- -(inputs$max_coverage() * 0.15)
            cyto_max <- -(inputs$max_coverage() * 0.05)
            
            cytoband_chr <- cytoband %>%
                filter(chr == gsub(x = inputs$chromosome(), pattern = "chr", replacement = "")) %>%
                mutate(
                    start = start,
                    end   = end,
                    ymin  = cyto_min,
                    ymax  = cyto_max
                )
            
            p_chr <- ggplot() +
                # Cytoband first (drawn below coverage)
                geom_rect(data = cytoband_chr,
                          aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = gieStain,
                              text = paste0(
                                  "Band: ", name,
                                  "<br>Chr: ", chr,
                                  "<br>Start: ", start,
                                  "<br>End: ", end
                              )),
                          inherit.aes = FALSE, show.legend = FALSE) +
                
                # Coverage line
                geom_line(data = bed_chr, aes(x = start, y = coverage),
                          color = "steelblue") +
                
                # CNVs if present
                {
                    if (nrow(cnv_chr) > 0) {
                        list(
                            geom_rect(data = cnv_chr,
                                      aes(xmin = start, xmax = end,
                                          ymin = min(bed_chr$coverage, na.rm = TRUE) - 1,
                                          ymax = max(bed_chr$coverage, na.rm = TRUE) + 1,
                                          fill = svtype,
                                          text = paste0(
                                              "CNV: ", svtype,
                                              "<br>CN: ", cn,
                                              "<br>Start: ", start,
                                              "<br>End: ", end
                                          )),
                                      alpha = 0.2),
                            geom_point(data = cnv_chr,
                                       aes(x = start, y = -5, color = svtype,
                                           text = paste0(
                                               "CNV: ", svtype,
                                               "<br>CN: ", cn,
                                               "<br>Chr: ", chr,
                                               "<br>Start: ", start,
                                               "<br>End: ", end
                                           )),
                                       shape = 17, size = 2),
                            scale_color_manual(values = c(DEL = "red", DUP = "green")) #points
                        )
                    }
                } +
                scale_fill_manual(values = c(gie_colors, DEL = "red", DUP = "green")) + # rectangles
                # Make space for ideogram below zero
                coord_cartesian(ylim = c(-(inputs$max_coverage() * 0.2),
                                         inputs$max_coverage())) +
                labs(
                    title = paste(inputs$sample(), "-", inputs$chromosome()),
                    x = "Genomic Position", y = "Mean Coverage"
                ) +
                theme_minimal() +
                guides(fill = "none")
            
            ggplotly(p_chr, tooltip = "text")
        })
        
    })
}
