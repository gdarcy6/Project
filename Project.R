
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(htmlwidgets)
library(R.utils)
library(shinycssloaders)
library(dplyr) # Manipulating datasets
library(DT)
library(tidyverse)
library(affy)
library(affyio)
library(GEOquery)
library(ggplot2)
library(plotly) # making graphs interactive
library(stringr)
library(ggthemes)
library(extrafont)
library(limma)
library(umap)
library(xgboost)
library(caret) # for confusion matirx
library(DiagrammeR) # visualising trees
library(rpart.plot) # visualising trees
library(pROC) # for ROC and AUC graphs
library(plotROC)
library(png)
#------------------------------------------------------------------------------


# load in the data
load("~/Project/normData_GSE36059.RData")
df <- normData

# make proper column names to match toptable 
fvarLabels(df) <- make.names(fvarLabels(df))


# group membership for all samples df
gsms1 <- paste0("11111111010000011111111100011001001001110111000011",
                "11110000110111110010100111011111000010110110111101",
                "00001011100110111111111111101000111110111001111100",
                "11111111111111110101011111111111110110110010101111",
                "10100110111101110111010100111101001101011011101110",
                "11011011110001111111111011111100111101001111111100",
                "10011111011110101011111110001110111111111111111111",
                "11111111110011110110110111101101111110111101111001",
                "100111000")
sml1 <- strsplit(gsms1, split="")[[1]]



# log2 transformation df
expr <- exprs(df)
qx <- as.numeric(quantile(expr, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { expr[which(expr <= 0)] <- NaN
exprs(df) <- log(expr) }


# assign samples to groups and set up design matrix 
gs1 <- factor(sml1)
groups <- make.names(c("rejection","non-rejection"))
levels(gs1) <- groups
df$group <- gs1
design1 <- model.matrix(~group + 0, df)
colnames(design1) <- levels(gs1)


# transpose the expression matrix
expr <- t(expr) 
# make expression matrix a data frame
expr <- as.data.frame(expr) 


# get the feature data with the assigned groups and remove the index
feat.dat <- as.data.frame(pData(df))
feat.dat$index = NULL


# merge the expression data and the feature data with assigned groups by Row.names
expr <- expr %>%
  merge(., feat.dat, by =0, all = TRUE)

# rename Row.names to samples
expr <- expr %>%
  rename(samples = Row.names)


# move group to the front of the expression matrix and make samples the row names
expr <- expr %>%
  relocate(., 'group', .after = 'samples') %>%
  column_to_rownames(., var = 'samples')

# createDataPartition() function from the caret package to split the original dataset into a training and testing set and split data into training (80%) and testing set (20%)
parts = createDataPartition(expr$group, p = 0.8, list = F)
train = expr[parts, ]
test = expr[-parts, ]


# seperate the training and test data into matrix and labels
train_labels <- train[,1]
train_data <-  data.matrix(train[,-1])

test_labels <- test[,1]
test_data <- data.matrix(test[,-1])


# create a matrix format that xgboost needs for training and test data
dtrain <- xgb.DMatrix(data = train_data, label = train_labels)
dtest <- xgb.DMatrix(data = test_data, label = test_labels)


# train a model using our training data
model <- xgboost(data = dtrain,                    # the data   
                 max.depth=3,                      # max depth 
                 nrounds=50)                       # max number of boosting iterations



# use model to make predictions on test data
pred_test = predict(model, dtest)


# convert prediction to factor type
pred_test[(pred_test>3)] = 3
pred_y = as.factor((levels(test_labels))[round(pred_test)])


# confusion matrix
conf_mat = confusionMatrix(test_labels, pred_y)

# xgb tree
xgbtree <- xgb.plot.tree(model = model, trees = 0:2, show_node_id = TRUE)

# auc plot
auc_xgb <- roc(test$group, pred_test)
auc_val <- round(auc(test$group, pred_test),4)


ROC_xgb_auc <- auc(auc_xgb)

# fit linear model
fit <- lmFit(df, design1) 

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2],sep ="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design1)
fit2 <- contrasts.fit(fit, cont.matrix)


# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)


# load in series matrix and get the feature data
gse <- getGEO(filename = "~/Project/GSE36059_series_matrix.txt.gz")

feature.data <- gse@featureData@data
# subset
feature.data <- feature.data[,c(1,11,10)]

# join table of top significant genes and feature data 
tT <- tT %>%
  rownames_to_column(var = 'ID') %>%
  inner_join(., feature.data, by = 'ID') 


tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene Symbol","Gene Title"))
tT <- column_to_rownames(tT, var = 'ID')





#------------------------------------------------------------------------------

### End of pre-processing

# Define UI for application that draws a histogram
ui <- dashboardPage(
                    skin = "blue",

    # Application title
    dashboardHeader(title = "Kidney Rejection Prediction"),

    # Sidebar 
    dashboardSidebar(
      sidebarMenu(
        menuItem(text = "Home", tabName = "home", icon = icon("home")),
        
        menuItem(text = "Prediction Model", tabName = "xgboost", icon = icon("th"), startExpanded = FALSE,
                 menuSubItem("AUC", tabName = "auc"),
                 menuSubItem("Tree", tabName = "tree")),
        
        menuItem(text = "Genomic Data", tabName = "genomicdata", icon = icon("dna"), startExpanded = FALSE,
                 menuSubItem("Gene Table", tabName = "genetable")),        
        
        menuItem(text = "About", tabName = "about", icon = icon("paper-plane"))
    )
),  
    # Main body
    dashboardBody(
      tabItems(
        tabItem("home",
                
                tags$h1(strong("Home")),
                
                br(),
                
                fluidPage(
                  splitLayout( cellWidths = 600,
                               box(title ="CEL FILES",
                                   solidHeader = T, status = "info", width = 6,
                                   
                                   fileInput("file", "Load CEL files(.CEL)",
                                             multiple = TRUE,
                                             accept = ".CEL",
                                             placeholder = "Please, insert the datapath of  CEL files")
                                  
                               ),
                 
                             
                  
                  box(
                    tableOutput("contents")))
                  )
                  ),
        
        
        
        

        
        tabItem("auc",
                fluidRow(
                  box(
                    title = "AUC for prediction model",
                    width = 6,
                    status = "info",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    plotOutput("xgb_auc")),
                  
                  box(
                    title = "Prediction Results Using Testing data",
                    width = 4,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    status= "danger",
                    verbatimTextOutput('confusion_matrix')
                           
                      )  
                  
                  )
                ),
        
        tabItem("tree",
                fluidRow(
                  box(
                    title = "Sample Tree of Prediction Model",
                    width = 5,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    status= "danger",
                    imageOutput("tree", height = "auto")),
                  
                  box(
                    title = "Genes of Importance in Prediction Model",
                    width = 5,
                    height = 12,
                    collapsible = TRUE,
                    plotlyOutput("gene_importance")
                  )
                )),
                
                

        
        
tabItem("genetable",
        fluidRow(
          box(
            
            title = "Genes Panel",
            width = 12,
            status = "warning",
            solidHeader = TRUE,
            collapsible = TRUE,
            label = boxLabel(
              text = 4,
              status = "danger",
              style = "circle"
            ),
            
            tabPanel("Genes",
                     status = "success",
                     dataTableOutput("genomic_plot") %>% withSpinner(color = "red"))))),
        
        tabItem("about",
                
                
                tags$h1(strong("Prediction of Kidney Rejection")),
                
                br(),
                
                fluidRow(
                  
                  box(
                  
                
                tags$h2(strong("Overview")),
                "The primary objective of this clinical web application is to give 
                clinicians real time predictions on kidney graft rejection."),
                
                box(
                  title = "About me",
                  status = "danger",
                  solidHeader = TRUE,
                  tags$p(
                    class = "text-center",
                    tags$img(class = "img-responsive img-rounded center-block", 
                             src = "", 
                             style = "max-width: 250px;")
                  ),
                  tags$p(
                    class = "text-center",
                    tags$strong("Hi! I'm Gabriel."),
                    HTML(paste0("(", 
                                tags$a(href = "https://twitter.com/gabrieldarcy", 
                                       "@gabrieldarcy"), ")"))
                  ),
                  tags$p(
                    "I'm a data scientist from Co. Galway, Ireland.",
                    br()
                  ),
                  tags$p(
                    "Get in touch with me on Twitter at",
                    HTML(paste0(
                      "(", tags$a(href = "", 
                                  "@gabrieldarcy", target = "_blank"), "),")),
                    
                    "LinkedIn at",
                    HTML(paste0(
                      "(", tags$a(href = "www.linkedin.com/in/gabriel-darcy", 
                                  "@gabriel-darcy", target = "_blank"), "),")),
                    
                    
                    "or by email at",
                    HTML(paste0(
                      tags$a(href = "gabrieldarcy94@gmail.com", 
                             "gabrieldarcy94@gmail.com"), "."))
                  )
                )
        
        
                )))
    )
  
  )
# Define server logic required to draw a histogram
server <- function(input, output) {
  
  rawData <- reactive({
   # validate(need(input$file != "", "Please Load CEL files"))
    
    read.celfile(input$file$datapath)
    
  })
  
  
  output$contents <- DT::renderDataTable({
    req(rawData())
    # normalise the input data
    eset <- rma(rawData)
    # get the expression of the input data and make into a dataframe
    exprd <- exprs(eset)
    exprd <- as.data.frame(exprd)
    
    datatable(exprd)
   
  })
  

  output$xgb_auc <- renderPlot({
    
    plot(auc_xgb, legacy.axis = TRUE, plot = TRUE, col = "#377eb8", print.auc = TRUE, grid = TRUE, asp = NA)
    
    
  })
  
  output$confusion_matrix <- renderPrint({
    
    print(conf_mat)
 })
  
  
  output$picture <- renderImage({
    return(list(src = "~/Project/WWW//tree.png",contentType = "image/png",alt = "tree"))
    
 },deleteFile = FALSE) #where the src is wherever you have the picture)
  
  output$gene_importance <- renderPlotly({
    importance_matrix = xgb.importance(colnames(train_data), model = model)
    
    feature.data1 <- feature.data %>%
      rename(., "Feature" = "ID")
    
    importance_matrix <- importance_matrix %>%
      inner_join(., feature.data1, by = 'Feature')
    
    importance_matrix$Feature <- paste(importance_matrix$Feature, importance_matrix$`Gene Symbol`)
    
    gene_importance <- xgb.ggplot.importance(importance_matrix[1:25,])
    
    ggplotly(gene_importance)
    
  })
  
 output$genomic_plot <- DT::renderDataTable({
   
   datatable(tT, rownames = FALSE)
 })
} 

# Run the application 
shinyApp(ui = ui, server = server)
