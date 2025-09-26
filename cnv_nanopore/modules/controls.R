mod_controls_ui <- function(id) {
    ns <- NS(id)
    tagList(
        selectInput(ns("chromosome"), "Chromosome:", choices = c(paste0("chr", 1:22), "chrX", "chrY"), selected = "chr1"),
        numericInput(ns("start_pos"), "Start position:", value = NA, min = 1),
        numericInput(ns("end_pos"), "End position:", value = NA, min = 1),
        numericInput(ns("max_coverage"), "Max coverage:", value = 250, min = 1),
        br(),
        h4("Circos settings"),
        h5("Include"),
        fluidRow(
            column(6,
                   checkboxInput(ns("circos_cnv"), "CNVs", value = FALSE)
                   ),
            column(6,
                   checkboxInput(ns("circos_sv"), "SVs", value = TRUE)
            ),
        ),
        h5("SV filters"),
        checkboxInput(ns("filter_tumor"), "Somatic filtering", value = TRUE),
        fluidRow(
            column(6,
            numericInput(ns("T_DV"), "T DV >= x:", value = 6, min = 0),
                   ),
            column(6,
            numericInput(ns("T_vaf"), "T VAF >= x:", value = 0.1, min = 0, max = 1),
            )
        ),        
        numericInput(ns("mapq"), "MAPQ >= x:", value = 50, min = 0),
        checkboxInput(ns("filter_normal"), "Germline Filtering:", value = FALSE),
        fluidRow(
            column(6,
            numericInput(ns("N_DV"), "N DV <= x:", value = 2, min = 0),
                   ),
            column(6,
            numericInput(ns("N_DR"), "N DR >= x:", value = 5, min = 0),
            ),
        ),
        actionButton(ns("reset_btn"), "Reset Inputs")
    )
}

is_all_numeric <- function(x) all(!is.na(suppressWarnings(as.numeric(x))))


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
            
            sv_somatic_file <- file.path(root_dir(), paste0(sample_name, ".SV_raw.tsv"))
            sv_all_file <- file.path(root_dir(), paste0("severus_all_", sample_name, ".vcf"))
            
            if (file.exists(cnv_file_plain)) {
                cnv_file <- cnv_file_plain
            } else if (file.exists(cnv_file_gz)) {
                cnv_file <- cnv_file_gz
            } else {
                stop("VCF file not found for sample: ", sample_name)
            }
            
            # SVs (optional)
            column_names <- c(
                "Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID",
                "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank",
                "HGVS.c", "HGVS.p", "cDNA.pos/cDNA.length", "CDS.pos/CDS.length",
                "AA.pos/AA.length", "Distance", "ERRORS/WARNINGS/INFO"
            )
            
            if (file.exists(sv_somatic_file) && file.exists(sv_all_file)) {
                sv_somatic <- read_tsv(sv_somatic_file, show_col_types = FALSE)
                
                vaf_col <- paste0(sample_name, ".VAF")
                hvaf_col <- paste0(sample_name, ".hVAF")
                dr_col <- paste0(sample_name, ".DR")
                dv_col <- paste0(sample_name, ".DV")
                
                clean_tbl <- sv_somatic %>%
                    dplyr::select(where(~ !all(is.na(.x) | .x == ""))) %>% # remove fully empty columns
                    mutate(across(where(is.factor), as.character)) %>%
                    mutate(across(where(is.list), ~ unlist(.x)))
                
                # split the annotation column, create base severus ID
                df <- clean_tbl %>%
                    separate_rows(ANN, sep = ",") %>%
                    separate(ANN, into = column_names, sep = "\\|", fill="right", extra="drop") %>%
                    dplyr::select(
                        CHROM, POS, ID, QUAL, FILTER, REF, ALT, PRECISE, IMPRECISE, SVTYPE, SVLEN, END,
                        DETAILED_TYPE, MAPQ, Annotation, Annotation_Impact, Gene_Name, Gene_ID,
                        Feature_Type, Transcript_BioType,
                        all_of(vaf_col), all_of(hvaf_col), all_of(dr_col), all_of(dv_col)
                    ) %>%
                    distinct() %>%
                    mutate(base_ID = sub("(_\\d+)$", "", ID)) %>% 
                    mutate(across(where(is.character), ~na_if(., ""))) 
                
                collapse_all <- function(x) paste(x, collapse = ";")
                
                df_collapsed <- df %>% 
                    separate(ALT, into = c("CHROM2", "POS2"), sep = ":", remove = FALSE) %>% 
                    distinct() %>% 
                    mutate(
                        CHROM2 = str_extract(CHROM2, "CHR[\\w]+"),
                        POS2   = as.numeric(str_extract(POS2, "\\d+")),
                        CHROM2 = str_replace(CHROM2, "CHR", "chr"),
                        CHROM2 = coalesce(CHROM2, CHROM),      # fallback for CHROM2
                        POS2   = coalesce(POS2, END, POS),      # fallback chain: POS2 → END → POS
                    ) %>% 
                    # dplyr::select(-END) %>% 
                    group_by(base_ID) %>% 
                    dplyr::select(-END) %>% 
                    summarise(
                        CHROM1 = CHROM[1],
                        POS1 = POS[1],
                        CHROM2 = CHROM2[1],
                        POS2 = POS2[1],
                        SVTYPE = SVTYPE[1],
                        QUAL = QUAL[1],
                        FILTER = FILTER[1],
                        REF = REF[1],
                        ALT = ALT[1],
                        MAPQ = MAPQ[1],
                        SVLEN = SVLEN[1],
                        PRECISE = PRECISE[1],
                        IMPRECISE = IMPRECISE[1],
                        DETAILED_TYPE = DETAILED_TYPE[1],
                        across(c(Annotation, Annotation_Impact, Gene_Name, Gene_ID, 
                                 Feature_Type, Transcript_BioType), collapse_all),
                        !!vaf_col := (!!sym(vaf_col))[1],
                        !!hvaf_col := (!!sym(hvaf_col))[1],
                        !!dr_col := (!!sym(dr_col))[1],
                        !!dv_col := (!!sym(dv_col))[1],
                        .groups = "drop"
                    ) %>% 
                    mutate(across(where(is.character), ~na_if(., ""))) %>% 
                    relocate(base_ID, CHROM1, POS1, CHROM2, POS2, SVTYPE, SVLEN) %>% 
                    filter(CHROM1 %in% paste0("chr", c(1:22, "X", "Y")),
                           CHROM2 %in% paste0("chr", c(1:22, "X", "Y")))
                
                sv_all <- read_tsv(
                    sv_all_file,
                    comment = "##",
                    col_select = c("ID", "tumor", "normal")
                ) %>%
                    mutate(base_ID = sub("(_\\d+)$", "", ID)) %>%
                    dplyr::select(-ID) %>%
                    distinct()
                
                sv_w_normal <- df_collapsed %>%
                    left_join(sv_all, by = "base_ID") %>%
                    separate(
                        normal,
                        into = c(
                            paste0(sample_name, "_N.GT"),
                            paste0(sample_name, "_N.GQ"),
                            paste0(sample_name, "_N.VAF"),
                            paste0(sample_name, "_N.hVAF"),
                            paste0(sample_name, "_N.DR"),
                            paste0(sample_name, "_N.DV")
                        ),
                        sep = ":"
                    ) %>%
                    separate(
                        tumor,
                        into = c(
                            paste0(sample_name, "_T.GT"),
                            paste0(sample_name, "_T.GQ"),
                            paste0(sample_name, "_T.VAF"),
                            paste0(sample_name, "_T.hVAF"),
                            paste0(sample_name, "_T.DR"),
                            paste0(sample_name, "_T.DV")
                        ),
                        sep = ":"
                    ) %>%
                    dplyr::select(-ends_with(".GQ")) %>%
                    mutate(
                        across(where(is.character), ~ if (is_all_numeric(.)) as.numeric(.) else .),
                        across(where(is.numeric), ~ round(., 2))
                    )
                
            } else {
                sv_w_normal <- NULL
                message("SV files not found for sample: ", sample_name, ". Skipping SV plots.")
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
            
            list( # this is data_list
                bed = bed,
                cnv = cnv_df_cov,
                coverage_ranges = coverage_ranges,
                sv = sv_w_normal
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
            
            # Reset SV filters to defaults
            updateNumericInput(session, "T_DV", value = 6)
            updateNumericInput(session, "T_vaf", value = 0.1)
            updateNumericInput(session, "mapq", value = 50)
            updateNumericInput(session, "N_DV", value = 2)
            updateNumericInput(session, "N_DR", value = 5)
            
            # Reset checkboxes to defaults
            updateCheckboxInput(session, "circos_cnv", value = FALSE)
            updateCheckboxInput(session, "circos_sv", value = TRUE)
            updateCheckboxInput(session, "filter_tumor", value = TRUE)
            updateCheckboxInput(session, "filter_normal", value = FALSE)
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
            circos_cnv   = reactive(input$circos_cnv),
            circos_sv    = reactive(input$circos_sv),
            filter_tumor = reactive(input$filter_tumor),
            min_tumor_DV  = reactive(input$T_DV),
            min_tumor_VAF = reactive(input$T_vaf),
            min_MAPQ      = reactive(input$mapq),
            filter_normal = reactive(input$filter_normal),
            max_normal_DV = reactive(input$N_DV),
            min_normal_DR = reactive(input$N_DR),
            sample       = sample_info$sample
        ))
    })
}
