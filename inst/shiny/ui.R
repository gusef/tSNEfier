require(shiny)
require(shinyjs)
require(d3Toolbox)
require(GSVA)
require(Rtsne)
require(RColorBrewer)
require(limma)
require(SummarizedExperiment)

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
                actionButton("runGSVA", "Run GSVA")
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
                d3ScatterOutput("scatter", width = "100%", height = "600"),
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
                ),
                actionButton("runtSNE", "Run tSNE"),
                selectizeInput("gene_filter", 
                               label = h5("Select gene space for the filter"), 
                               choices = NULL)
            ),
            column(6,
                   actionButton("diffGenes", "Differential Genes"),
                   actionButton("diffPathways", "Differential Pathways"),
                   dataTableOutput('diffTable')

            )
        )
    )
)
