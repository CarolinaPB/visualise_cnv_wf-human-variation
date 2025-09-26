library(shiny)
library(tidyverse)
library(plotly)
library(VariantAnnotation)
library(here)
library(config)
library(shinyFiles)
library(ggplot2)
library(circlize)


i_am("cnv_nanopore/app.R")

source(here("cnv_nanopore/utils.R"))
source(here("cnv_nanopore/modules/select_dir.R"))
source(here("cnv_nanopore/modules/select_sample.R"))
source(here("cnv_nanopore/modules/controls.R"))
source(here("cnv_nanopore/modules/plots.R"))
source(here("cnv_nanopore/modules/circlize.R"))


# sv_command <- paste(
#     'gatk VariantsToTable \\',
#     '-V "$vcf" \\',
#     '-O "$out" \\',
#     '-F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER \\',
#     '-F PRECISE -F IMPRECISE -F SVTYPE -F SVLEN -F CHR2 -F END -F STRANDS \\',
#     '-F DETAILED_TYPE -F INSLEN -F MAPQ -F PHASESETID -F HP -F CLUSTERID \\',
#     '-F INSSEQ -F MATE_ID -F INSIDE_VNTR -F ALINGED_POS \\',
#     '-F ANN -F LOF -F NMD \\',
#     '-F DBVARID -F ALLELEID -F CLNSIG -F CLNVCSO -F SCIDNINCL -F CLNREVSTAT \\',
#     '-F ONCREVSTAT -F RS -F CLNDNINCL -F ONC -F ORIGIN -F ONCINCL -F ONCDNINCL \\',
#     '-F ONCDISDB -F SCIREVSTAT -F ONCDISDBINCL -F MC -F CLNDN -F ONCCONF \\',
#     '-F CLNVC -F SCIDISDB -F CLNVI -F AF_EXAC -F ONCDN -F AF_ESP -F CLNSIGINCL \\',
#     '-F CLNDISDB -F GENEINFO -F CLNDISDBINCL -F AF_TGP -F CLNSIGCONF \\',
#     '-F SCIDISDBINCL -F CLNHGVS -F SCIINCL -F SCIDN -F SCI \\',
#     '-GF GT -GF GQ -GF DR -GF DV -GF VAF -GF hVAF',
#     sep = "\n"
# )


ui <- fluidPage(
    titlePanel("EPI2ME results Visualisation"),
    
    sidebarLayout(
        sidebarPanel(
            dirSelectUI("select_dir"),
            sample_selector_ui("sample"),
            mod_controls_ui("controls"),
            width = 2
        ),
        
        mainPanel(
            tabsetPanel(
                # tabPanel(
                    # "Intro",
                    # p("For more details see the README")
    #                 br(),
    #                 h3("Overview"),
    #                 p("This app allows visualization of CNV and SV results from Nanopore sequencing. 
    # Select a directory containing your data, choose a sample, and adjust the controls to customize the plots."),
    #                 
    #                 h3("Input Files"),
    #                 tags$ul(
    #                     tags$li(strong("wf-human-variation:"),
    #                             tags$ul(
    #                                 tags$li(tags$code("SAMPLE.wf_cnv.vcf.gz")),
    #                                 tags$li(tags$code("SAMPLE.regions.bed.gz"))
    #                             )
    #                     ),
    #                     tags$li(strong("wf-somatic-variation:"),
    #                             tags$ul(
    #                                 tags$li(tags$code("SAMPLE.SV_raw.tsv")),
    #                                 tags$li(tags$code("SAMPLE/sv/severus-output/all_SVs/SAMPLE.severus_all.vcf"))
    #                             )
    #                     )
    #                 ),
    #                 
    #                 h3("Processing SV Files"),
    #                 p("Use the following command to convert the SV VCF file SAMPLE.wf-somatic-sv.vcf.gz to SAMPLE.SV_raw.tsv"),
    #                 code(sv_command),
    #                 
    #                 h3("Circos plots"),
    #                 p("The CNV plotting in the circos plots is off by default since it takes a long time to render.")
                # ),
                tabPanel("CNVs", mod_plots_ui("plots")),
                tabPanel("Circos Plot", 
                         br(),
                         mod_circos_circlize_ui("circosPlot"))
            ),
            width = 9
        )
    )
)



server <- function(input, output, session) {
    root_dir <- dirSelectServer("select_dir")
    sample_info <- sample_selector_server("sample", root_dir)
    control_vals <- mod_controls_server("controls", root_dir, sample_info)
    plots_res <- mod_plots_server("plots", control_vals)
    mod_circos_circlize_server("circosPlot", control_vals, plots_res)
}

shinyApp(ui, server)
