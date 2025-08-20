load_cytoband <- function(cytoband_file, gie_colors) {
    read.table(cytoband_file, header = TRUE, sep = "\t",
               col.names = c("chr", "start", "end", "name", "gieStain"),
               stringsAsFactors = FALSE) %>%
        mutate(chr = sub("^chr", "", chr),
               color = gie_colors[gieStain]) %>%
        as_tibble()
}

compute_chrom_offsets <- function(bed, ordered_chrs) {
    bed %>%
        mutate(chr = factor(chr, levels = ordered_chrs)) %>%
        group_by(chr) %>%
        summarise(chr_len = max(end), .groups = "drop") %>%
        arrange(chr) %>%
        mutate(offset = lag(cumsum(chr_len), default = 0),
               chr_index = row_number())
}

# build_facet_plot <- function(data_list, sample_name, ordered_chrs) {
#     bed <- data_list$bed %>% mutate(chr = factor(chr, levels = ordered_chrs))
#     cnv <- data_list$cnv %>% mutate(chr = factor(chr, levels = ordered_chrs))
#     
#     p <- ggplot() +
#         geom_line(data = bed,
#                   aes(x = start, y = coverage),
#                   color = "steelblue") +
#         geom_rect(
#             data = cnv,
#             aes(
#                 xmin = start,
#                 xmax = end,
#                 ymin = ymin,
#                 ymax = ymax,
#                 fill = svtype,
#                 color = svtype,
#                 text = paste0(
#                     "CNV: ",
#                     svtype,
#                     "<br>CN: ",
#                     cn,
#                     "<br>Chr: ",
#                     chr,
#                     "<br>Start: ",
#                     start,
#                     "<br>End: ",
#                     end
#                 )
#             ),
#             alpha = 0.2
#         ) +
#         geom_point(
#             data = cnv,
#             aes(
#                 x = start,
#                 y = ymax,
#                 color = svtype,
#                 text = paste0(
#                     "CNV Start: ",
#                     start,
#                     "<br>Chr: ",
#                     chr,
#                     "<br>Type: ",
#                     svtype,
#                     "<br>CN: ",
#                     cn
#                 ),
#                 fill = svtype
#             ),
#             shape = 25,
#             size = 2
#         ) +
#         scale_fill_manual(values = c(DEL = "red", DUP = "green")) +
#         scale_color_manual(values = c(DEL = "red", DUP = "green")) +
#         facet_wrap( ~ chr, scales = "free", ncol = 4) +
#         labs(
#             title = sample_name,
#             x = "Genomic Position",
#             y = "Mean Coverage",
#             fill = "CNV Type",
#             color = "CNV Type"
#         ) + theme_minimal() + theme(
#             panel.grid = element_blank(),
#             axis.text.x = element_blank(),
#             axis.ticks.x = element_blank()
#         )
#     
#     ggplotly(p, tooltip = "text")
# }

build_chr_plot <- function(data_list,
                           chromosome,
                           start_pos,
                           end_pos,
                           max_cov,
                           sample_name,
                           cytoband,
                           gie_colors) {
    bed_chr <- data_list$bed %>% filter(chr == chromosome, start >= start_pos, end <= end_pos)
    cnv_chr <- data_list$cnv %>% filter(chr == chromosome) %>% filter(
        (start >= start_pos & start <= end_pos) |
            (end >= start_pos & end <= end_pos) |
            (start <= start_pos & end >= end_pos)
    )
    
    cyto_min <- -(max_cov * 0.15)
    cyto_max <- -(max_cov * 0.05)
    cytoband_chr <- cytoband %>% filter(chr == gsub("chr", "", chromosome)) %>%
        mutate(ymin = cyto_min, ymax = cyto_max)
    
    p <- ggplot() +
        geom_rect(
            data = cytoband_chr,
            aes(
                xmin = start,
                xmax = end,
                ymin = ymin,
                ymax = ymax,
                fill = gieStain,
                text = paste0(
                    "Band: ",
                    name,
                    "<br>Chr: ",
                    chr,
                    "<br>Start: ",
                    start,
                    "<br>End: ",
                    end
                )
            ),
            show.legend = FALSE
        ) +
        {
            if (nrow(cnv_chr) > 0)
                list(
                    geom_rect(
                        data = cnv_chr,
                        aes(
                            xmin = start,
                            xmax = end,
                            ymin = min(bed_chr$coverage) -
                                1,
                            ymax = max(bed_chr$coverage) +
                                1,
                            fill = svtype,
                            color = svtype,
                            text = paste0(
                                "CNV: ", svtype,
                                "<br>CN: ", cn,
                                "<br>Start: ", start,
                                "<br>End: ", end,
                                "<br>Width: ", end - start
                            )
                        ),
                        alpha = 0.2
                    ),
                    geom_point(
                        data = cnv_chr,
                        aes(
                            x = start,
                            y = -(max_cov * 0.02),
                            color = svtype,
                            text = paste0(
                                "CNV: ", svtype,
                                "<br>CN: ", cn,
                                "<br>Start: ", start,
                                "<br>End: ", end,
                                "<br>Width: ", end - start
                            )
                        ),
                        shape = 17,
                        size = 2
                    ),
                    scale_color_manual(values = c(
                        DEL = "red", DUP = "green"
                    ))
                )
        } +
        geom_line(
            data = bed_chr,
            aes(
                x = start,
                y = coverage,
                group = 1,
                text = paste0(
                    "<br>Chr: ",
                    chr,
                    "<br>Start: ",
                    start,
                    "<br>End: ",
                    end,
                    "<br>Coverage: ",
                    coverage
                )
            ),
            color = "steelblue"
        ) +
        scale_fill_manual(values = c(gie_colors, DEL = "red", DUP = "green")) +
        coord_cartesian(ylim = c(-(max_cov * 0.2), max_cov),
                        xlim = c(start_pos, end_pos)) +
        labs(
            title = paste(sample_name, "-", chromosome),
            x = "Genomic Position",
            y = "Mean Coverage"
        ) +
        theme_minimal() + guides(fill = "none")
    ggplotly(p, tooltip = "text")
}

build_manhattan_plot <- function(data_list,
                                 cytoband,
                                 max_cov,
                                 ordered_chrs,
                                 gie_colors,
                                 sample_name) {
    chrom_len <- compute_chrom_offsets(data_list$bed, ordered_chrs)
    
    bed_global <- data_list$bed %>% mutate(chr = factor(chr, levels = ordered_chrs)) %>%
        left_join(chrom_len, "chr") %>% mutate(global_mid = (start + end) /
                                                   2 + offset)
    cnv_global <- data_list$cnv %>% mutate(chr = factor(chr, levels = ordered_chrs)) %>%
        left_join(chrom_len, "chr") %>% mutate(global_start = start + offset,
                                               global_end = end + offset)
    
    cytoband_global <- cytoband %>% mutate(chr = paste0("chr", chr)) %>%
        left_join(chrom_len, "chr") %>% mutate(
            global_start = start + offset,
            global_end = end + offset,
            ymin = -(max_cov * 0.2),
            ymax = -(max_cov * 0.15)
        )
    
    chrom_bg <- chrom_len %>% filter(chr_index %% 2 == 0) %>% mutate(xmin =
                                                                         offset, xmax = offset + chr_len)
    
    p <- ggplot() +
        geom_rect(
            data = chrom_bg,
            aes(
                xmin = xmin,
                xmax = xmax,
                # ymin = -(max_cov * 0.2),
                ymin = -(max_cov * 1.1),
                ymax = max_cov * 1.1
            ),
            fill = "lightgray",
            alpha = 0.4
        ) +
        geom_rect(
            data = cytoband_global,
            aes(
                xmin = global_start,
                xmax = global_end,
                ymin = ymin,
                ymax = ymax,
                text = paste0(
                    "Band: ",
                    name,
                    "<br>Chr: ",
                    chr,
                    "<br>Start: ",
                    start,
                    "<br>End: ",
                    end
                ),
                fill = gieStain
            ),
        ) +
        geom_rect(
            data = cnv_global,
            aes(
                xmin = global_start,
                xmax = global_end,
                fill = svtype,
                color = svtype,
                text = paste0(
                    "CNV: ", svtype,
                    "<br>CN: ", cn,
                    "<br>Start: ", start,
                    "<br>End: ", end,
                    "<br>Width: ", end - start
                )
            ),
            ymin = 0,
            ymax = max_cov*1.1,
            alpha = 0.3
        ) +
        geom_point(
            data = cnv_global,
            aes(
                x = global_start,
                y = -(max_cov * 0.1),
                fill = svtype,
                color = svtype,
                text = paste0(
                    "CNV: ", svtype,
                    "<br>CN: ", cn,
                    "<br>Start: ", start,
                    "<br>End: ", end,
                    "<br>Width: ", end - start
                )
            ),
            shape = 23,
            size = 3,
            alpha = 0.6
        ) +
        geom_line(
            data = bed_global,
            aes(x = global_mid, y = coverage, group = 1,
                text = paste0(
                    "<br>Chr: ",
                    chr,
                    "<br>Start: ",
                    start,
                    "<br>End: ",
                    end,
                    "<br>Coverage: ",
                    coverage
                )),
            color = "steelblue",
            size = 0.3
        ) +
        scale_x_continuous(breaks = chrom_len$offset + chrom_len$chr_len /
                               2,
                           labels = chrom_len$chr) +
        coord_cartesian(ylim = c(-(max_cov * 0.2), max_cov)) + theme_bw() +
        theme(panel.grid = element_blank(),
              axis.text.x =
                  element_text(angle = 45, hjust = 1)) + labs(x = sample_name, y =
                                                                  "Coverage") +
        scale_fill_manual(values = c(gie_colors, DEL = "red", DUP = "green")) +
        scale_color_manual(values = c(DEL = "red", DUP = "green"))+
        guides(fill = "none")
    ggplotly(p,tooltip = "text")
}


mod_plots_ui <- function(id) {
    ns <- NS(id)
    tagList(plotlyOutput(ns("manhattanChrPlot")), 
            br(), 
            plotlyOutput(ns("chrPlot")))
}

mod_plots_server <- function(id, inputs, config_section = "default") {
    moduleServer(id, function(input, output, session) {
        config <- config::get(
            config = config_section,
            file = here::here("cnv_nanopore/configs/config.yaml")
        )
        gie_colors <- c(
            gneg = "ghostwhite",
            gpos25 = "gray70",
            gpos50 = "gray50",
            gpos75 = "gray30",
            gpos100 = "black",
            gvar = "gray90",
            stalk = "brown",
            acen = "red"
        )
        cytoband <- load_cytoband(here(config$cytoband_file), gie_colors)
        ordered_chrs <- c(paste0("chr", 1:22), "chrX", "chrY")
        
        output$manhattanChrPlot <- renderPlotly({
            req(
                inputs$sample(),
                inputs$data_list(),
                !is.na(inputs$max_coverage())
            )
            
            build_manhattan_plot(
                inputs$data_list(),
                cytoband,
                inputs$max_coverage(),
                ordered_chrs,
                gie_colors,
                inputs$sample()
            )
        })
        output$chrPlot <- renderPlotly({
            req(
                inputs$sample(),
                inputs$data_list(),
                inputs$chromosome(),
                !is.na(inputs$start_pos()),
                !is.na(inputs$end_pos()),
                !is.na(inputs$max_coverage())
            )
            
            
            build_chr_plot(
                inputs$data_list(),
                inputs$chromosome(),
                inputs$start_pos(),
                inputs$end_pos(),
                inputs$max_coverage(),
                inputs$sample(),
                cytoband,
                gie_colors
            )
        })
    })
}
