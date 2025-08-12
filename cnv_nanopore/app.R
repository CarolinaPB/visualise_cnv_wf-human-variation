library(shiny)
library(tidyverse)
library(plotly)
library(VariantAnnotation)
library(here)
library(config)
library(shinyFiles)
library(ggplot2)


i_am("cnv_nanopore/app.R")

source(here("cnv_nanopore/utils.R"))
source(here("cnv_nanopore/modules/select_dir.R"))
source(here("cnv_nanopore/modules/select_sample.R"))
source(here("cnv_nanopore/modules/controls.R"))
source(here("cnv_nanopore/modules/plots.R"))



ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            dirSelectUI("select_dir"),
            sample_selector_ui("sample"),
            mod_controls_ui("controls"),width = 2,
            
        ),
        mainPanel(
            mod_plots_ui("plots"), width = 10,
        )
    )
)

server <- function(input, output, session) {
    root_dir <- dirSelectServer("select_dir")
    sample_info <- sample_selector_server("sample", root_dir)
    control_vals <- mod_controls_server("controls", root_dir, sample_info)
    mod_plots_server("plots", control_vals)
}

shinyApp(ui, server)
