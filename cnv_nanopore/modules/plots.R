library(plotly)

mod_plots_ui <- function(id) {
    ns <- NS(id)
    tagList(
        plotlyOutput(ns("facetPlot"), width = "100%", height = "100%"),
        br(),
        plotlyOutput(ns("chrPlot"))
    )
}

mod_plots_server <- function(id, inputs) {
    moduleServer(id, function(input, output, session) {
        
        output$facetPlot <- renderPlotly({
            req(inputs$sample())
            req(inputs$data_list())
            
            dat <- inputs$data_list()
            ordered_chrs <- c(paste0("chr", 1:22), "chrX", "chrY")
            
            bed <- dat$bed %>% mutate(chr = factor(chr, levels = ordered_chrs))
            cnv <- dat$cnv %>% mutate(chr = factor(chr, levels = ordered_chrs))
            
            # Ensure coverage ranges are present (as in original code)
            coverage_ranges <- dat$coverage_ranges
            
            
            p_facet <- ggplot() +
                geom_line(data = bed, aes(x = start, y = coverage), color = "steelblue") +
                geom_rect(
                    data = cnv,
                    aes(
                        xmin = start, xmax = end,
                        ymin = ymin, ymax = ymax,
                        fill = svtype,
                        text = paste0(
                            "CNV: ", svtype,
                            "<br>CN: ", cn,
                            "<br>Chr: ", chr,
                            "<br>Start: ", start,
                            "<br>End: ", end
                        )
                    ),
                    alpha = 0.2
                ) +
                geom_point(
                    data = cnv,
                    aes(
                        x = start,
                        y = ymax,
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
                    title = paste("Coverage + CNVs - Sample:", inputs$sample()),
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
            req(inputs$sample())       # ensure sample is chosen
            req(inputs$data_list())    # ensure data is loaded
            req(inputs$chromosome())   # ensure chromosome is set
            
            if (inputs$chromosome() == "All") return(NULL)
            
            dat <- inputs$data_list()
            bed_chr <- dat$bed %>% filter(chr == inputs$chromosome())
            cnv_chr <- dat$cnv %>% filter(chr == inputs$chromosome())
            
            bed_chr <- bed_chr %>% filter(start >= inputs$start_pos(), end <= inputs$end_pos())
            cnv_chr <- cnv_chr %>%
                filter(
                    (start >= inputs$start_pos() & start <= inputs$end_pos()) |
                        (end >= inputs$start_pos() & end <= inputs$end_pos()) |
                        (start <= inputs$start_pos() & end >= inputs$end_pos())
                )
            
            p_chr <- ggplot() +
                geom_line(data = bed_chr, aes(x = start, y = coverage), color = "steelblue") +
                coord_cartesian(ylim = c(0, inputs$max_coverage())) +
                labs(
                    title = paste("Coverage + CNVs -", inputs$chromosome(), "- Sample:", inputs$sample()),
                    x = "Genomic Position", y = "Mean Coverage"
                ) +
                theme_minimal()
            
            if (nrow(cnv_chr) > 0) {
                p_chr <- p_chr +
                    geom_rect(
                        data = cnv_chr,
                        aes(
                            xmin = start, xmax = end,
                            ymin = min(bed_chr$coverage) - 1,
                            ymax = max(bed_chr$coverage) + 1,
                            fill = svtype,
                            text = paste0(
                                "CNV: ", svtype,
                                "<br>CN: ", cn,
                                "<br>Chr: ", chr,
                                "<br>Start: ", start,
                                "<br>End: ", end
                            )
                        ),
                        alpha = 0.2
                    ) +
                    geom_point(
                        data = cnv_chr,
                        aes(x = start, y = 10, color = svtype,
                            text = paste0(
                                "CNV: ", svtype,
                                "<br>CN: ", cn,
                                "<br>Chr: ", chr,
                                "<br>Start: ", start,
                                "<br>End: ", end
                            )),
                        shape = 17, size = 2
                    ) +
                    scale_fill_manual(values = c(DEL = "red", DUP = "green")) +
                    scale_color_manual(values = c(DEL = "red", DUP = "green"))
            }
            
            ggplotly(p_chr, tooltip = "text")
        })
    })
}
