require(shiny)
require(shinyjs)
require(D3Scatter)
require(GSVA)
require(Rtsne)
require(RColorBrewer)

ui <- navbarPage(title = "tSNEfier",
    tabPanel(title = "File loading",
        fluidRow(
            includeCSS(system.file('www','tSNEfier',package='tSNEfier')),
            column(6,
                fileInput("eSet_file", 
                          "Load gene expression data (.rds)",
                          accept = c(".rds",".RDS")),
                fileInput("gmt_file", 
                          "Load gene set file (.gmt)",
                          accept = c(".gmt")),
                actionButton("runGSVA", "Run GSVA"),
                actionButton("runtSNE", "Run tSNE")
            ),
            column(6,
                fileInput("state_file", 
                          "Load a saved state (.rds)",
                          accept = c(".rds",".RDS")),
                downloadButton("SaveStateButton", "Save current state")
            )
        )    
    ),
    tabPanel(title = "Display panel",
        fluidRow(
            column(6,
                D3ScatterOutput("scatter", width = "100%", height = "600"),
                fluidRow(
                    column(4,
                        selectizeInput("gene_color", 
                                    label = h5("Select gene"), 
                                    choices = NULL)
                    ),
                    column(4,
                        selectizeInput("pathway_color", 
                                    label = h5("Select pathway"), 
                                    choices = NULL)
                    ),
                    column(4,
                    selectizeInput("covariate_color", 
                                label = h5("Select covariate"), 
                                choices = NULL)
                    )
                )
            ),
            column(6,
                   verbatimTextOutput("currentOutput")
            )
        )
    )
)

                 



                                    
                           
