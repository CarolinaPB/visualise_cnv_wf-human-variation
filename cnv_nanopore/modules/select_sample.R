sample_selector_ui <- function(id) {
    ns <- NS(id)
    tagList(
        # verbatimTextOutput(ns("dir")),
        uiOutput(ns("sample_ui"))
    )
}

sample_selector_server <- function(id, root_dir) {
    moduleServer(id, function(input, output, session) {
        
        
        # Get available samples whenever directory changes
        samples <- reactive({
            req(root_dir(), dir.exists(root_dir()))
            get_samples(root_dir())
        })
        
        # Render the selectInput dynamically
        output$sample_ui <- renderUI({
            ns <- session$ns
            samp <- samples()
            if (length(samp) == 0) {
                return(tags$em("No samples found in this directory"))
            }
            selectInput(ns("sample"), "Select sample:", choices = samp, selected = samp[1])
        })
        
        # Return selected sample + list of all samples
        return(list(
            sample = reactive(input$sample),
            samples = samples
        ))
    })
}
