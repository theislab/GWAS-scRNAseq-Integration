# loading libraries

library(shiny)
library(Seurat)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(melt)

# loading data
setwd('C:/Users/user/GWAS-scRNAseq-Integration/data/') # change here

if(!length(list.files(path = './getMarkers/')) > 0)
  unzip(zipfile = './getMarkers.zip')

tissues <- gsub('.csv', '', list.files(path = './getMarkers/'), fixed = TRUE)

if(!file.exists('gwas_catalog_v1.0.2-associations_e92_r2018-05-12.tsv'))
  unzip(zipfile = './gwas_catalog_v1.0.2-associations_e92_r2018-05-12.tsv.zip')

gwas <- read.delim('./gwas_catalog_v1.0.2-associations_e92_r2018-05-12.tsv')

diseases <- as.character(unique(gwas$DISEASE.TRAIT))

geneList <- unique(unlist(lapply(as.character(gwas$REPORTED.GENE.S.), function(x) ifelse(grepl(',',x),strsplit(x, ",", fixed = TRUE),x))))
if('intergenic' %in% geneList){
  geneList <-  geneList[-which(geneList=='intergenic')] 
}
genes <- as.character(geneList)

# Define UI for application that draws a histogram
ui <- fluidPage(theme="simplex.min.css",
          tags$style(type="text/css", 
                     "label {font-size: 12px;}", 
                     ".recalculating {opacity: 1.0;}"
                 ),
   
  
  # Application title
  tags$h2("Integrating GWAS with Mouse Cell Atlas"),
  p("Finding enrichments in mouse cell types from ",
    tags$a(href="http://bis.zju.edu.cn/MCA/", "MCA"),
    "for published rare and common variants associated with disease traits from ",
    tags$a(href="https://www.ebi.ac.uk/gwas/", "GWAS Catalog")),
  hr(),

  # Sidebar
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = 'tissue', label = 'Select a tissue for exploring MCA cell-type signatures:', choices = tissues),
      
      sliderInput("no_of_markers", "Select number of markers to display per cell-type:",
                  min = 0, max = 15,
                  value = 5),
      
      radioButtons(inputId = "do", label = "Show results for",
                   choices = c("tSNE representation" = "tSNE",
                               "Heatmap for top 5 markers per cell-type" = "heatmap"
                               ),
                   selected = "heatmap"),
      helpText('Note: The heatmap does not update if the slider for number of markers is changed'),
      # Horizontal line ----
      tags$hr(),
      p("Select a combination of GWAS disease trait and mouse tissue for calculating enrichments"),
      
      
      selectizeInput(inputId = 'disease', label = 'Select a disease trait for exploring GWAS associations:', choices = diseases, options = list(maxOptions = length(diseases))),
      radioButtons(inputId = "what", label = "Show results for",
                   choices = c("Studies" = "studies",
                               "Associations" = "associations",
                               "Catalog traits" = "traits"),
                   selected = "studies"),
      tags$hr(),
      selectizeInput(inputId = 'gene', label = 'Select a gene for exploring its GWAS associations and expression patterns in MCA cell-types:', choices = genes, options = list(maxOptions = length(genes))),
      radioButtons(inputId = "which", label = "Show results for",
                   choices = c("GWAS hits" = "GWAS_hits",
                               "Expression profiles" = "exp"),
                   selected = "GWAS_hits"),
      tags$hr()
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("MCA table", fluidRow(
                    DT::dataTableOutput("contents")
                  )),
                  tabPanel("MCA Visualization", plotOutput("plot")),
                  tabPanel("GWAS table", fluidRow(
                    DT::dataTableOutput("GWAS"))),
                  tabPanel("GWAS Enrichments", fluidRow(DT::dataTableOutput("enrichments"))),
                  #tabPanel("GWAS Visualization", plotOutput("GWASplot")),
                  tabPanel("GWAS Visualization", fluidRow(
                    column(12, align = 'center', tags$h2("GWAS and MCA scRNA-seq enrichments"))
                  ),fluidRow(
                    column(6, align = 'center', plotOutput("GWASplot_all", height = "600px")))
                  , fluidRow(
                    column(6, align = 'center', plotOutput("GWASplot", height = "600px"))
                  )),
                  tabPanel("Gene-wise Visualization",
                    fluidRow(
                    column(12, plotOutput("GWAShits", height = "300px")))
                  , fluidRow(
                    column(12, plotOutput("GWAShits_all", height = "300px"))
                  ))
                  
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  MCA_preprocessing <- function(markers){
    # removing NA columns
    markers[sapply(markers, function(x) all(is.na(x)))] <- NULL
    # removing NA rows
    markers <- markers[complete.cases(markers),]
    # some files have an additional adjusted p value column, replace and remove
    if('p_val_adj' %in% colnames(markers)){
      markers$p_val <- markers$p_val_adj
      markers <- markers[-which(colnames(markers) == 'p_val_adj')]
    }
    # rename columns
    colnames(markers) <- c('p_val','avg_diff','pct.1','pct.2','cluster','gene','celltype')
    # drop factors, necessary for plotting
    if(is.factor(markers$avg_diff))
      markers$avg_diff <- as.numeric(as.vector(markers$avg_diff))
    return(markers)
  }
  
  getData1 <- reactive({
    
    markers <- read.csv(paste0('./getMarkers/',input$tissue,'.csv'))
    markers <- MCA_preprocessing(markers)
    
    celltypes <- as.character(unique(markers$celltype))
    celltypes <- celltypes[celltypes != ""]
    celltypes <- celltypes[celltypes != " "]
    
    topn_markers <- markers %>% group_by(cluster) %>% top_n(input$no_of_markers, avg_diff)
    topn_markers$celltype = rep(celltypes, table(topn_markers$cluster))
    
    return(topn_markers)
  })
  
  getData2 <- reactive({
    
    markers <- read.csv(paste0('./getMarkers/',input$tissue,'.csv'))
    markers <- MCA_preprocessing(markers)
    
    celltypes <- as.character(unique(markers$celltype))
    celltypes <- celltypes[celltypes != ""]
    celltypes <- celltypes[celltypes != " "]
    
    top5_markers <- markers %>% group_by(cluster) %>% top_n(5, avg_diff)
        
    # preparing heatmap
    mat <- as.data.frame(matrix(data = 0, nrow = length(unique(top5_markers$gene)), ncol = length(unique(markers$cluster))))
    rownames(mat) <- unique(top5_markers$gene)
    colnames(mat) <- unique(markers$cluster)
    
    for(cluster in colnames(mat)){
      for(gene in rownames(mat)){
        if(length(markers$avg_diff[which(markers$gene == gene & markers$cluster == cluster)]))
          mat[gene,cluster] <- markers$avg_diff[which(markers$gene == gene & markers$cluster == cluster)]
      }
    }
    colnames(mat) <- celltypes
    
    return(mat)
    
  })
  
  getData3 <- reactive({
    
    gwas_sub <- gwas[which(gwas$DISEASE.TRAIT == input$disease),]
    return(gwas_sub)
    
  })
  
  getData4 <- reactive({
  
    # lung-GWAS integration
    
    trait <- gwas[which(gwas$DISEASE.TRAIT == input$disease),]
    
    # genes associated with 'trait' (reported)
    geneList <- unique(unlist(lapply(as.character(trait$REPORTED.GENE.S.), function(x) ifelse(grepl(',',x),strsplit(x, ",", fixed = TRUE),x))))
    if('intergenic' %in% geneList){
      geneList <-  geneList[-which(geneList=='intergenic')] 
    }
    
    # reading in MCA markers
    markers <- read.csv(paste0('./getMarkers/',input$tissue,'.csv'))
    markers <- MCA_preprocessing(markers)
    
    celltypes <- as.character(unique(markers$celltype))
    celltypes <- celltypes[celltypes != ""]
    celltypes <- celltypes[celltypes != " "]
    
    df <- as.data.frame(matrix(data = 0, nrow = length(unique(markers$gene)), ncol = length(celltypes)))
    rownames(df) <- unique(markers$gene)
    colnames(df) <- celltypes
    
    for(x in 1:ncol(df)){
      df[unique(as.character(markers$gene[which(markers$cluster == x)])),x] = 1
    }

    rownames(df) <- toupper(rownames(df))
    
    # applying Fisher enrichment test
    matrix_fisher <- data.frame()
    
    v <- sapply(colnames(df), FUN = function(x){
      present_disease <- geneList # trait
      present_anno <- rownames(df)[which(df[,x]==1)]
      intersection <- intersect(present_anno, present_disease)
      # make a contingency table
      cont_tbl <- matrix(c(abs(length(rownames(df)) - length(unique(union(present_anno, present_disease)))),
                           length(present_anno) - length(intersection),
                           length(present_disease) - length(intersection),
                           length(intersection)), ncol = 2)
      rownames(cont_tbl) <- c('not_in_anno','in_anno')
      colnames(cont_tbl) <- c('not_in_disease','in_disease')
      # Performing Fisher's test
      res <- fisher.test(x = cont_tbl)#S, alternative = 'greater')
      # Saving p-value, odds ratio, common genes between disease markers from GWAS catalog 
      # and celltype signature, no. of genes in annotation and
      # no. of markers in cluster
      vec <- c(as.numeric(res$p.value), 
               res$estimate, 
               length(present_anno), 
               length(present_disease), 
               ifelse(test = length(intersection) > 1, 
                      paste(intersection,collapse = '|'), intersection), length(intersection)
      )
      return(vec)
    })
    
    matrix_fisher <- rbind(matrix_fisher,cbind(rownames(t(v)),data.frame(t(v), row.names=NULL)))
    colnames(matrix_fisher) <- c('annotation','p.value', 'estimate.of.odds.ratio', 'present.in.celltype', 
                                 'present.in.disease', 'genes.in.intersection','length.of.intersection')
    matrix_fisher$p.value <- as.numeric(as.vector(matrix_fisher$p.value))
    return(matrix_fisher)
    
  })
  
  # get enrichments for selected trait across all tissue|cell-types, display top 15 hits
  getData5 <- reactive({
    load("./GWAS_enrichments_allcomparisons.RData")
    tmp<- l[[which(names(l) == input$disease)]]
    tmp <- tmp[order(tmp$p.value, decreasing = FALSE)[1:15],] # top 15 hits
    return(tmp)
  })
  # get GWAS hits for selected gene
  getData6 <- reactive({
    geneList <- unlist(lapply(as.character(gwas$REPORTED.GENE.S.), function(x) ifelse(grepl(',',x),strsplit(x, ",", fixed = TRUE),x)))
    if('intergenic' %in% geneList ){
      geneList <-  geneList[-which(geneList=='intergenic')] 
    }
    if('Intergenic' %in% geneList ){
      geneList <-  geneList[-which(geneList=='Intergenic')] 
    }
    a <- table(geneList)
    a[names(a) == input$gene]
    
    x <- lapply(split(1:nrow(gwas), gwas$DISEASE.TRAIT), function(y){
      tmp <- gwas[y,]
      tmpList <- unlist(lapply(as.character(tmp$REPORTED.GENE.S.), function(x) ifelse(grepl(',',x),strsplit(x, ",", fixed = TRUE),x)))
      if('intergenic' %in% tmpList){
        tmpList <-  tmpList[-which(tmpList=='intergenic')] 
      }
      if(input$gene %in% tmpList)
        unique(tmp$DISEASE.TRAIT)
    })
    v <- c(a[names(a) == input$gene],length(x[!unlist(lapply(x, is.null))]))
    names(v) <- c("no. of GWAS hits","no. of GWAS disease traits implicated")
    return(v)
  })
  
  getData7 <- reactive({
    x <- lapply(split(1:nrow(gwas), gwas$DISEASE.TRAIT), function(y){
    tmp <- gwas[y,]
    tmpList <- unlist(lapply(as.character(tmp$REPORTED.GENE.S.), function(x) ifelse(grepl(',',x),strsplit(x, ",", fixed = TRUE),x)))
    if('intergenic' %in% tmpList){
      tmpList <-  tmpList[-which(tmpList=='intergenic')] 
    }
    if(input$gene %in% tmpList)
      c(tmp$P.VALUE[which(tmp$REPORTED.GENE.S == input$gene)][1])
  })
    z <-x[!unlist(lapply(x, is.null))]
    df <- melt(z)
    df <- df[complete.cases(df),]
    df$pvalue <- -log10(df$value)
    df <- df[order(df$pvalue, decreasing = TRUE),]
    
    return(df)
  })
  
  # displaying MCA cell-type markers table for selected tissue
  output$contents <- DT::renderDataTable(DT::datatable({
    
    topn_markers <- getData1()
    colnames(topn_markers) <- c('p value','avg difference','pct.1','pct.2','cluster','gene','celltype')
    topn_markers
  }))
  
  # displaying heatmap (p values of marker genes across cell-types) & tSNE (to do) for selected MCA tissue
  output$plot <- renderPlot({
    
    if(input$do == "heatmap"){
      pheatmap(mat = getData2(), 
               legend = TRUE, 
               cluster_rows = TRUE,
               cluster_cols = TRUE, 
               border_color = 'red', 
               main = "Heatmap of fold changes for top 5 markers across clusters", 
               color = brewer.pal(n = 8, "YlGn"))
    }},height = 2000, width = 1500)
    
  # displaying GWAS Catalog for selected GWAS disease trait
  output$GWAS <- DT::renderDataTable(DT::datatable({
    gwas_sub <- getData3()
    if(input$what == "studies"){
        tmp <- gwas_sub[,c('FIRST.AUTHOR','PUBMEDID','STUDY.ACCESSION','DATE','JOURNAL','STUDY','DISEASE.TRAIT')]
        colnames(tmp) <- c('Author','PMID','Study accession','Publication Date','Journal','Title','Reported trait')
    }
    if(input$what == "associations"){
      tmp <- gwas_sub[,c('STRONGEST.SNP.RISK.ALLELE','RISK.ALLELE.FREQUENCY','P.VALUE','OR.or.BETA','X95..CI..TEXT.','REGION','CHR_ID','CHR_POS','CONTEXT','REPORTED.GENE.S.','MAPPED_GENE','DISEASE.TRAIT','STUDY.ACCESSION','PUBMEDID','FIRST.AUTHOR')]
      colnames(tmp) <- c('SNP','RAF','p-value','OR/Beta','CI','Region','CHR ID','CHR Pos','Functional class','Reported gene(s)','Mapped gene(s)','Reported trait','Study accession','PMID','Author')
    }
    if(input$what == "traits"){
      tmp <- gwas_sub[,c('MAPPED_TRAIT','DISEASE.TRAIT','PUBMEDID')]
      colnames(tmp) <- c('Trait','Reported trait','PMID')
    }
    
    tmp
  }))
  
  # displaying Fisher's Exact test enrichment results between selected disease trait
  # and MCA tissue
  output$enrichments <- DT::renderDataTable(DT::datatable({
    matrix_fisher <- getData4()
    colnames(matrix_fisher) <- c('Cell-type','p value', 'Estimate of odds ratio', 'Present in cell-type', 
                                 'Present in trait', 'Genes in intersection','Length of intersection')
    matrix_fisher
  }))
  
  output$GWASplot <- renderPlot({
    # bold.12.text <- element_text(face = "bold", color = "black", size = 12)
    ggplot(getData4(), aes(x = reorder(annotation, -p.value), y = -log10(p.value))) + 
      geom_bar(stat = "identity", fill="royal blue", colour="black") + 
      labs(x = "cell-types (selected tissue)", y = "-log10(p value)") + 
      coord_flip() + 
      theme(plot.title = element_text(hjust = 2)) + geom_hline(yintercept = -log10(0.05))
  }, height = 500, width = 500)
  
  output$GWASplot_all <- renderPlot({
    ggplot(getData5(), aes(x = reorder(annotation, -p.value), y = -log10(p.value))) + 
      geom_bar(stat = "identity", fill="tomato", colour="black") + 
      labs(x = "cell-types (all tissues)", y = "-log10(p value)") + 
      coord_flip() +
      theme(plot.title = element_text(hjust = 2)) + geom_hline(yintercept = -log10(0.05))
  }, height = 400, width = 400)
  
  output$GWAShits <- renderPlot({
    barplot(getData6(), horiz = FALSE, col = c('royal blue','tomato'), main = paste0('Gene: ',input$gene), ylab = 'frequency')
  })
  
  output$GWAShits_all <- renderPlot({
    ggplot(getData7(), aes(x = reorder(L1, pvalue), y = pvalue)) + 
      geom_bar(stat = "identity", fill="tomato", colour="black") + 
      labs(x = "GWAS disease traits implicated in gene", y = "-log10(p value)") + 
      coord_flip() + 
      theme(plot.title = element_text(hjust = 2)) + geom_hline(yintercept = -log10(0.05))
  }, height = 600, width = 1200)


}

# Run the application 
shinyApp(ui = ui, server = server)