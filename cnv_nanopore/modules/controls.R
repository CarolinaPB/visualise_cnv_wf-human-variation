mod_controls_ui <- function(id) {
    ns <- NS(id)
    tagList(
        selectInput(ns("chromosome"), "Chromosome:", choices = c(paste0("chr", 1:22), "chrX", "chrY"), selected = "chr1"),
        numericInput(ns("start_pos"), "Start position:", value = NA, min = 1),
        numericInput(ns("end_pos"), "End position:", value = NA, min = 1),
        numericInput(ns("max_coverage"), "Max coverage:", value = NA, min = 1),
        actionButton(ns("reset_btn"), "Reset Inputs")
    )
}

mod_controls_server <- function(id, root_dir, sample_info) {
    moduleServer(id, function(input, output, session) {
        
        # Load the data whenever the directory or sample changes
        data_list <- reactive({
            req(root_dir(), sample_info$sample())  # sample_info$sample is reactive
            sample_name <- sample_info$sample()
            
            message("Loading data for sample: ", sample_name, " from ", root_dir())
            
            bed_file <- file.path(root_dir(), paste0(sample_name, ".regions.bed.gz"))
            cnv_file_plain <- file.path(root_dir(), paste0(sample_name, ".wf_cnv.vcf"))
            cnv_file_gz <- paste0(cnv_file_plain, ".gz")
            
            if (file.exists(cnv_file_plain)) {
                cnv_file <- cnv_file_plain
            } else if (file.exists(cnv_file_gz)) {
                cnv_file <- cnv_file_gz
            } else {
                stop("VCF file not found for sample: ", sample_name)
            }
            
            valid_chrs <- c(paste0("chr", 1:22), "chrX", "chrY")
            
            # Load BED
            bed <- read_tsv(bed_file, col_names = c("chr", "start", "end", "coverage"), progress = FALSE) %>%
                filter(chr %in% valid_chrs)
            
            # Load VCF
            vcf <- readVcf(cnv_file)
            cnv_df <- tibble(
                chr = as.character(seqnames(rowRanges(vcf))),
                start = start(rowRanges(vcf)),
                end = as.numeric(info(vcf)$END),
                svtype = as.character(info(vcf)$SVTYPE),
                cn = as.numeric(info(vcf)$CN)
            ) %>% filter(chr %in% valid_chrs)
            
            coverage_ranges <- bed %>%
                group_by(chr) %>%
                summarise(
                    ymin = min(coverage, na.rm = TRUE) - 1,
                    ymax = max(coverage, na.rm = TRUE) + 1,
                    .groups = "drop"
                )
            
            cnv_df_cov <- cnv_df %>%
                left_join(coverage_ranges, by = "chr")
            
            print(cnv_df_cov)
            
            list(
                bed = bed,
                cnv = cnv_df_cov,
                coverage_ranges = coverage_ranges
            )
        })
        
        # Reset inputs on button press
        observeEvent(input$reset_btn, {
            updateSelectInput(session, "chromosome", selected = "chr1")
            
            dat <- data_list()
            bed_chr <- dat$bed %>% filter(chr == "chr1")
            
            start_min <- min(bed_chr$start, na.rm = TRUE)
            end_max <- max(bed_chr$end, na.rm = TRUE)
            cov_min <- floor(min(bed_chr$coverage, na.rm = TRUE))
            cov_max <- ceiling(max(bed_chr$coverage, na.rm = TRUE))
            
            updateNumericInput(session, "start_pos", value = start_min, min = start_min, max = end_max)
            updateNumericInput(session, "end_pos", value = end_max, min = start_min, max = end_max)
            updateNumericInput(session, "max_coverage", value = 250, min = cov_min, max = cov_max)
        })
        
        observeEvent(
            c(input$chromosome, sample_info$sample()), 
            {
                dat <- data_list()
                
                chr_selected <- input$chromosome
                
                bed_chr <- dat$bed %>% filter(chr == chr_selected)
                
                if(nrow(bed_chr) > 0) {
                    start_min <- min(bed_chr$start, na.rm = TRUE)
                    end_max <- max(bed_chr$end, na.rm = TRUE)
                    cov_min <- floor(min(bed_chr$coverage, na.rm = TRUE))
                    cov_max <- ceiling(max(bed_chr$coverage, na.rm = TRUE))
                    
                    updateNumericInput(session, "start_pos", value = start_min, min = start_min, max = end_max)
                    updateNumericInput(session, "end_pos", value = end_max, min = start_min, max = end_max)
                    updateNumericInput(session, "max_coverage", value = 250, min = cov_min, max = cov_max)
                }
            }, 
            ignoreNULL = TRUE, 
            ignoreInit = FALSE
        )
        
        # Return all controls + sample + data_list to plotting module
        return(list(
            data_list    = data_list,
            chromosome   = reactive(input$chromosome),
            start_pos    = reactive(input$start_pos),
            end_pos      = reactive(input$end_pos),
            max_coverage = reactive(input$max_coverage),
            sample       = sample_info$sample
        ))
    })
}
