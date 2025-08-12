dirSelectUI <- function(id) {
    ns <- NS(id)
    tagList(
        shinyFiles::shinyDirButton(ns("dir"), "Choose directory", "Please select a folder"),
        verbatimTextOutput(ns("path"))
    )
}

dirSelectServer <- function(id, config_section = "default") {
    moduleServer(id, function(input, output, session) {
        volumes <- c(Home = fs::path_home(), Root = "/")
        shinyFiles::shinyDirChoose(input, "dir", roots = volumes, session = session)
        
        print(here)
        config <- config::get(config = config_section, file = here::here("cnv_nanopore/configs/config.yaml"))
        print(config)
        default_dir <- here(config$data_path)
        
        dirPath <- reactive({
            if (is.integer(input$dir)) {
                return(default_dir)
            } else {
                parsed_path <- shinyFiles::parseDirPath(volumes, input$dir)
                return(parsed_path)
            }
        })
        
        output$path <- renderText({
            req(dirPath())
            if (!dir.exists(dirPath())) {
                warning_msg <- paste("Warning: Directory does not exist:", dirPath())
                showNotification(warning_msg, type = "warning")
                return(warning_msg)
            }
            dirPath()
        })
        
        return(dirPath)
    })
}
