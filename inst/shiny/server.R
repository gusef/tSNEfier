rm(list=ls())
gc()


update_inputs <- function(session, input, values){

    #update gene select
    selection <- rownames(values$eSet)
    if (input$gene_filter != 'All'){

        selection <- selection[selection %in% values$gene_sets[[input$gene_filter]]] 
    }
    genes <- lapply(selection,function(x)x)
    names(genes) <- selection
    updateSelectizeInput(session,
                         'gene_color', 
                         choices = genes)
    
    #update pathway select
    updateSelectizeInput(session,
                         'pathway_color', 
                         choices = values$gene_set_select[-1])
    
    #update covariate select
    cova <- lapply(names(colData(values$eSet)),function(x)x)
    names(cova) <- names(colData(values$eSet))
    updateSelectizeInput(session,
                         'covariate_color',
                         choices = cova)
}

server <- function(input, output, session) {
    options(shiny.maxRequestSize=500*1024^2)
    
    #reactive values
    values <- reactiveValues(diffTable = NULL,
                             eSet = NULL,
                             gene_sets = NULL,
                             gene_set_desc = NULL,
                             gene_set_select = NULL,
                             gsva_scores = NULL,
                             tSNE = NULL,
                             tSNE_color = "black",
                             tSNE_legend = NULL,
                             tSNE_title = "",
                             tSNE_space = "All",
                             verbose = NULL)
     
    #loading ExpressionSet
    observeEvent(input$eSet_file,{ 
        if (is.null(input$eSet_file)){
            return (NULL)
        }
        #read in the gene expression object and figure out whether it is a eSet or Summarized experiment
        dat <- readRDS(input$eSet_file$datapath)
                
        if (class(dat) == "ExpressionSet"){
            dat <- makeSummarizedExperimentFromExpressionSet(dat)
        }else if(class(dat) != "SummarizedExperiment"){
            showModal(modalDialog(
                title = "Error",
                "Gene expression set set to be an 'ExpressionSet' or 'SummarizedExperiment' class"
            ))
            return(NULL)
        } 
        values$gsva <- NULL
        values$eSet <- dat
    })
    
    #loading gmt file
    observeEvent(input$gmt_file,{ 
        if (is.null(values$eSet)){
            return(NULL)
        }
        values$gsva <- NULL
        
        #preprocess sets
        sets <- readLines(input$gmt_file$datapath)
        sets <- strsplit(sets,'\t')
        names(sets) <- sapply(sets,function(x)x[1])
        desc <- sapply(sets,function(x)x[2])
        sets <- sapply(sets,function(x)x[-(1:2)])
        
        #make sure they are actually represented
        sets <- lapply(sets,function(x,nams)x[x %in% nams],rownames(values$eSet))
        filter <- sapply(sets,length) > 0
        if (sum(filter) == 0){
            showModal(modalDialog(
                title = "Error",
                "None of the genes in the gene sets were represented in the expression set."
            ))
        } else{
            values$gene_set_desc <- desc[filter]
            values$gene_sets <- sets[filter]
            
            #make a list so selectize is easier
            values$gene_set_select <- c('All', lapply(names(sets),function(x)x))
            names(values$gene_set_select) <- c('All', names(sets))
            
            #update gene space filter
            updateSelectizeInput(session,
                                 'gene_filter', 
                                 choices = values$gene_set_select,
                                 selected = values$gene_set_select[2])
        }
    })
    
    #if the runGSVA button is pressed
    observeEvent(input$runGSVA, {
        if (is.null(values$gene_sets) || is.null(values$eSet)){
            return(NULL)
        }
        
        set.seed(123) 
        withProgress(message = 'Running GSVA ... ',
                     detail = 'This may take a while...',
                     value = 0.5, {
                         
            values$gsva_scores <- gsva(assay(values$eSet), values$gene_sets)$es.obs
        })
    })
    
    #run tSNE button
    observeEvent(input$runtSNE, {
        if (is.null(values$eSet)){
            return (NULL)
        } 
        if (is.null(values$gene_sets)){
            values$gene_set_select <- list('All' = "All") 
        }
        showModal(modalDialog(
            size='l',
            title = "Running t-SNE",
            selectInput("tSNE_pathway", 
                        label = h3("Select the gene set space"), 
                        choices = values$gene_set_select,
                        selected = values$gene_set_select[2]),
            numericInput(inputId = 'tSNE_iter',
                         value = 1000,
                         max = 5000,
                         min = 100,
                         label = 'tSNE iterations'),
            numericInput(inputId = 'tSNE_theta',
                         value = 0.5,
                         max = 0.1,
                         min = 1,
                         label = 'Theta'),
            numericInput(inputId = 'tSNE_perplexity',
                         value = 25,
                         max = 2,
                         min = 500,
                         label = 'perplexity'),
            easyClose = TRUE,
            footer = tagList(
                modalButton("Cancel"),
                actionButton("tsne_ok", "Run"))                
        ))
    })
    
    #if run tsne button has been pushed
    observeEvent(input$tsne_ok, {
        mat <- t(assay(values$eSet))
        
        #indicate what space was used to derive the tSNE
        values$tSNE_space <- input$tSNE_pathway
        
        if (input$tSNE_pathway != 'All'){
            mat <- mat[,values$gene_sets[[input$tSNE_pathway]]]     
        }

        withProgress(message = 'Running t-SNE ... ',
                     detail = 'This may take a while...',
                     value = 0.5, {
                        set.seed(123) 
                        values$tSNE <- Rtsne(mat,
                                             max_iter = input$tSNE_iter,
                                             perplexity = input$tSNE_perplexity,
                                             theta = input$tSNE_theta,
                                             verbose = TRUE)
                     })
        removeModal()
        update_inputs(session, input, values)
    })
    
    #loading a saved state
    observeEvent(input$state_file,{ 
        vals <- readRDS(input$state_file$datapath)
        for (idx in names(vals)){
            values[[idx]] <- vals[[idx]]
        }     
        if (!is.null(values$tSNE)){
            update_inputs(session, input, values)
        }
        
        #update gene space filter
        if (!is.null(values$gene_set_select)){
            updateSelectizeInput(session,
                                 'gene_filter', 
                                 choices = values$gene_set_select,
                                 selected = values$gene_set_select[2])
        }
    })
    
    observeEvent(input$gene_filter,{
        if (!is.null(values$eSet)){
            update_inputs(session, input, values)
        }
    })
    
    #save the current state
    output$SaveStateButton <- downloadHandler(
        filename = function() {
            paste("State_", Sys.Date(), ".rds", sep="")
        },
        content = function(file) {
            state <- list()
            for (idx in names(values)){
                state[[idx]] <- values[[idx]]
            }
            saveRDS(state,file=file)},
        contentType='rds')
    
    
    observeEvent(input$gene_color,{
        if (is.null(values$eSet) ||
            input$gene_color == ""){
            return(NULL)
        }
        values$tSNE_legend <- NULL
        values$tSNE_color <- assay(values$eSet)[input$gene_color,]
        values$tSNE_title <- input$gene_color
    })
    
    observeEvent(input$pathway_color,{
        if (is.null(values$gsva_scores) ||
            input$pathway_color == ""){
            return (NULL)
        }
        values$tSNE_legend <- NULL
        values$tSNE_color <- values$gsva_scores[input$pathway_color,]
        values$tSNE_title <- input$pathway_color
    })
    
    observeEvent(input$covariate_color,{
        if (is.null(values$eSet) ||
            input$covariate_color == ""){
            return(NULL)
        }
        #set the title
        values$tSNE_title <- input$covariate_color
        
        #extract values
        val <- values$eSet[[input$covariate_color]]
        
        #check the data type
        if (is.numeric(val)){
            values$tSNE_color <- val
        } else {
            val <- as.factor(val)
            recist <- c('#f03b20','#ffeda0','#99d8c9','#2ca25f')
            names(recist) <- c('PD','SD','PR','CR')

            if (all(levels(val) %in%  names(recist))){
                #RECIST covariate
                values$tSNE_color <- recist[match(val,names(recist))]
                values$tSNE_legend <- data.frame(col=recist,
                                                 name=names(recist))
            } else {
                if (length(levels(val))<=12){
                    #if there is <= than 12 classes in the covariate
                    pal <- brewer.pal(length(levels(val)),'Paired')
                    values$tSNE_legend <- data.frame(col=pal,
                                                     name=levels(val))
                } else {
                    #if there is more than 12
                    pal <- rainbow(length(levels(val)))
                    values$tSNE_legend <- NULL
                }
                #set up the coloring
                values$tSNE_color <- pal[as.numeric(val)]

            }
        }
    })
    
    output$scatter <- renderd3Scatter({
        if (!is.null(values$tSNE)){

            dat <- data.frame(x=values$tSNE$Y[,1],
                              y=values$tSNE$Y[,2],
                              names=colnames(values$eSet))
            d3Scatter(dat,
                      col=values$tSNE_color,
                      dotsize = 5,
                      xlab=paste(values$tSNE_space,'- 1st axis'),
                      ylab=paste(values$tSNE_space,'- 2nd axis'),
                      title=values$tSNE_title,
                      tooltip = c('names'),
                      legend = values$tSNE_legend,
                      callback='ScatterSelection')
        }
    })
    
    observeEvent(input$diffGenes, {
        if (is.null(input$ScatterSelection)){
            return(NULL)
        }    
        showModal(modalDialog(
            title = "Differential Expression",
            selectInput("diff_pathway", 
                        label = h3("Select the gene set space"), 
                        choices = values$gene_set_select,
                        selected = values$gene_set_select[2]),
            easyClose = TRUE,
            footer = tagList(
                modalButton("Cancel"),
                actionButton("diff_ok", "Run"))
        ))
    })    

    observeEvent(input$diff_ok, {    
        dat <- assay(values$eSet)
        
        if (input$diff_pathway != 'All'){
            dat <- dat[values$gene_sets[[input$diff_pathway]],]     
        }
        
        removeModal()
        
        #generate label according to selection
        label <- (1:ncol(values$eSet)) %in% as.numeric(input$ScatterSelection)
        label <- c('notSelected','Selected')[label+1]
        
        #differential expression
        design <- model.matrix(~1+label)
        fit <- lmFit(dat, design)
        fit <- eBayes(fit)
        
        topTable <- topTable(fit, 
                             number=500, 
                             adjust = "fdr",
                             coef=colnames(fit$coefficients)[2])
        
        nams <- colnames(topTable)
        topTable <- cbind(rownames(topTable),apply(topTable,2,format,digits=2))
        colnames(topTable) <- c('Gene','logFC','avg. expression',
                                't-stat','p-value','FDR','B')
        topTable <- topTable[,-ncol(topTable)]
        values$diffTable <- topTable
    })
    

    
    observeEvent(input$diffPathways, {    
        if (is.null(input$ScatterSelection)){
            return(NULL)
        }    
        
        dat <- values$gsva_scores
        
        #generate label according to selection
        label <- (1:ncol(values$eSet)) %in% as.numeric(input$ScatterSelection)
        label <- c('notSelected','Selected')[label+1]
        
        #extract p-value
        res <- sapply(1:nrow(dat),
                      function(x, dat, lab)
                          wilcox.test(dat[x, lab == unique(lab)[1]],
                                      dat[x, lab == unique(lab)[2]])$p.value,
                      dat,
                      label)
        
        #W stat
        stat <- sapply(1:nrow(dat),
                        function(x, dat, lab)
                            wilcox.test(dat[x, lab == unique(lab)[1]],
                                         dat[x, lab == unique(lab)[2]])$statistic,
                        dat,
                        label)
        
        #results
        wilcox_res <- data.frame(
            Pathway = gsub('_',' ',rownames(dat)),
            W = stat,
            p.value = res,
            FDR = p.adjust(res,
                           method = 'BH',
                           n = length(res))
        )
        wilcox_res <- wilcox_res[order(wilcox_res$p.value,decreasing = F), ]
        wilcox_res$p.value <- format(wilcox_res$p.value, digits = 4 )
        wilcox_res$FDR <- format(wilcox_res$FDR, digits = 4 )
        
        values$diffTable <- wilcox_res
    })
        
    output$diffTable = renderDataTable({
        values$diffTable
    })
    

}



