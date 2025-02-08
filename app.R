# app.R

library(shiny)
library(edgeR)
library(limma)
library(org.Ss.eg.db)
library(RColorBrewer)
library(gplots)
library(tidyverse)
library(AnnotationDbi)

ui <- fluidPage(
    titlePanel("RNAseq Differential Expression Analysis"),
    sidebarLayout(
        sidebarPanel(
            helpText("Upload your count data and sample information files (CSV format)."),
            fileInput("countsFile", "Choose Count Data CSV File", accept = ".csv"),
            fileInput("sampleFile", "Choose Sample Info CSV File", accept = ".csv"),
            actionButton("run", "Run Analysis"),
            hr(),
            downloadButton("downloadDEG", "Download DEG Table")
        ),
        mainPanel(
            tabsetPanel(
                tabPanel("Library Sizes", plotOutput("libSizesPlot")),
                tabPanel("MD Plot", plotOutput("mdPlot")),
                tabPanel("Boxplots", plotOutput("boxplotPlot")),
                tabPanel("MDS Plot", plotOutput("mdsPlot")),
                tabPanel("Heatmap", plotOutput("heatmapPlot")),
                tabPanel("DEG Table", tableOutput("degTable")),
                tabPanel("Full Model Output", tableOutput("fullFitTable"))
            )
        )
    )
)

server <- function(input, output, session) {
    
    # Reactive expression to read the count data file
    countsData <- reactive({
        req(input$countsFile)
        tryCatch({
            read.csv(input$countsFile$datapath, stringsAsFactors = FALSE, row.names = 1)
        }, error = function(e) {
            showNotification("Error reading count data file", type = "error")
            return(NULL)
        })
    })
    
    # Reactive expression to read the sample information file
    sampleInfo <- reactive({
        req(input$sampleFile)
        tryCatch({
            read.csv(input$sampleFile$datapath, stringsAsFactors = FALSE)
        }, error = function(e) {
            showNotification("Error reading sample info file", type = "error")
            return(NULL)
        })
    })
    
    # Run the analysis when the button is clicked
    analysisResults <- eventReactive(input$run, {
        req(countsData(), sampleInfo())
        
        # Load data
        seqdata <- countsData()
        sampInfo <- sampleInfo()
        
        # Create design matrix from the "Condition" column
        condition <- factor(sampInfo$Condition)
        design <- model.matrix(~ 0 + condition)
        
        # Create DGEList object (group is based on condition)
        dge <- DGEList(counts = seqdata, group = condition, remove.zeros = FALSE)
        
        # Annotation: Map gene SYMBOLs to ENTREZ IDs
        ENTREZID <- mapIds(org.Ss.eg.db, rownames(dge$counts), keytype = "SYMBOL", column = "ENTREZID")
        rownames(dge$counts) <- ENTREZID
        ann <- AnnotationDbi::select(org.Ss.eg.db, keys = rownames(dge$counts),
                                     columns = c("ENTREZID", "SYMBOL", "GENENAME"))
        dge$genes <- ann
        
        # Remove genes with missing annotation
        keepAnnot <- !is.na(dge$genes$ENTREZID)
        dge <- dge[keepAnnot, ]
        dge <- dge[!is.na(rownames(dge$counts)), ]
        
        # Filtering low-expressed genes and normalization
        keep <- filterByExpr(dge, design)
        dge <- dge[keep, , keep.lib.size = FALSE]
        dge <- calcNormFactors(dge)
        
        # Voom transformation with quality weights (without plotting)
        v <- voomWithQualityWeights(dge, design = design, plot = FALSE)
        
        # Define contrasts (adjust these if your conditions differ)
        my.contrasts <- makeContrasts(
            Treated_vs_Control = conditionHF_Treated - conditionControl,
            Placebo_vs_Control = conditionHF_Placebo - conditionControl,
            Treated_vs_Placebo = conditionHF_Treated - conditionHF_Placebo,
            levels = design
        )
        
        # Fit linear model and apply contrasts
        vfit <- lmFit(v, design)
        vfit <- contrasts.fit(vfit, contrasts = my.contrasts)
        vfit <- eBayes(vfit, robust = TRUE)
        
        # Differential expression results for the first contrast (Treated_vs_Control)
        DEG <- topTable(vfit, coef = "Treated_vs_Control", number = Inf, sort.by = "none")
        DEG_sig <- DEG %>% filter(adj.P.Val < 0.05)
        
        # Extract full model output using topTable (using sort.by = "none")
        fullFit <- topTable(vfit, number = Inf, sort.by = "none")
        
        # Prepare diagnostic objects:
        # Library sizes (in millions)
        libSizes <- dge$samples$lib.size / 1e6
        
        # MD Plot: We'll use the first sample as an example
        mdPlotData <- list(vfit = vfit, dt = decideTests(vfit))
        
        # Boxplots: Unnormalized and normalized logCPM
        logCPM_unnorm <- cpm(dge, log = TRUE)
        logCPM_norm   <- cpm(dge, log = TRUE, normalized.lib.sizes = TRUE)
        
        # MDS Plot: Use pseudoCounts (log2(counts + 1))
        pseudoCounts <- log2(dge$counts + 1)
        
        # Heatmap: Top 100 most variable genes based on normalized logCPM
        logCPM_all <- cpm(dge, log = TRUE)
        var_genes <- apply(logCPM_all, 1, var)
        select_var <- names(sort(var_genes, decreasing = TRUE))[1:min(100, length(var_genes))]
        heatmapData <- logCPM_all[select_var, ]
        
        list(
            dge = dge,
            design = design,
            v = v,
            vfit = vfit,
            DEG = DEG_sig,
            fullFit = fullFit,
            libSizes = libSizes,
            mdPlotData = mdPlotData,
            logCPM_unnorm = logCPM_unnorm,
            logCPM_norm = logCPM_norm,
            pseudoCounts = pseudoCounts,
            heatmapData = heatmapData,
            sampInfo = sampInfo
        )
    })
    
    ### Outputs ###
    
    # 1. Library Sizes Barplot
    output$libSizesPlot <- renderPlot({
        req(analysisResults())
        barplot(analysisResults()$libSizes,
                names.arg = colnames(analysisResults()$dge$counts),
                las = 2,
                col = "lightskyblue",
                space = 0.5,
                ylab = "Library size (millions)")
        title("Barplot of Library Sizes")
    })
    
    # 2. MD Plot for the first sample (example)
    output$mdPlot <- renderPlot({
        req(analysisResults())
        plotMD(analysisResults()$vfit, column = 1,
               xlab = "Average log CPM",
               ylab = "log-fold change",
               main = "MD Plot (Sample 1)")
        abline(h = 0, col = "red", lty = 2, lwd = 2)
    })
    
    # 3. Boxplots: Unnormalized vs Normalized logCPM
    output$boxplotPlot <- renderPlot({
        req(analysisResults())
        par(mfrow = c(1, 2))
        boxplot(analysisResults()$logCPM_unnorm,
                las = 2,
                main = "Unnormalized logCPM",
                ylab = "Log2 CPM")
        abline(h = median(as.matrix(analysisResults()$logCPM_unnorm)), col = "blue")
        
        boxplot(analysisResults()$logCPM_norm,
                las = 2,
                main = "Normalized logCPM",
                ylab = "Log2 CPM")
        abline(h = median(as.matrix(analysisResults()$logCPM_norm)), col = "blue")
        par(mfrow = c(1, 1))
    })
    
    # 4. MDS Plot
    output$mdsPlot <- renderPlot({
        req(analysisResults())
        grp <- analysisResults()$dge$samples$group
        cols <- brewer.pal(max(3, length(levels(grp))), "Set1")[as.numeric(grp)]
        plotMDS(analysisResults()$pseudoCounts, col = cols, pch = 16, main = "MDS Plot")
        legend("topright", legend = levels(grp),
               col = brewer.pal(max(3, length(levels(grp))), "Set1")[1:length(levels(grp))],
               pch = 16)
    })
    
    # 5. Heatmap of Top Variable Genes
    output$heatmapPlot <- renderPlot({
        req(analysisResults())
        mypalette <- brewer.pal(11, "RdYlBu")
        morecols <- colorRampPalette(mypalette)
        heatmap.2(analysisResults()$heatmapData,
                  col = rev(morecols(50)),
                  trace = "none",
                  scale = "row",
                  dendrogram = "both",
                  margins = c(12, 8),
                  srtCol = 45,
                  main = "Heatmap of Top Variable Genes")
    })
    
    # 6. DEG Table (Differentially Expressed Genes)
    output$degTable <- renderTable({
        req(analysisResults())
        analysisResults()$DEG
    }, rownames = TRUE)
    
    # 7. Full Model Output (results from topTable on vfit)
    output$fullFitTable <- renderTable({
        req(analysisResults())
        head(analysisResults()$fullFit)
    }, rownames = TRUE)
    
    # 8. Download Handler for DEG table
    output$downloadDEG <- downloadHandler(
        filename = function() {
            paste("DEG_table_", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
            req(analysisResults())
            write.csv(analysisResults()$DEG, file, row.names = TRUE)
        }
    )
    
}

shinyApp(ui = ui, server = server)
