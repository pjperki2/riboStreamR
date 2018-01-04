library(shiny)
library(ggplot2)
library(riboSeqR)
library(Rsamtools)
library(systemPipeR)
library(GenomicFeatures)
library(cowplot)
library(grid)
shinyUI(navbarPage(inverse = T,"riboStreamR",
                tabPanel("Start", #Start page for welcome message, uploading files, and data preprocessing
                         sidebarLayout(
                           sidebarPanel(
                             #actionButton("go", "Go"),
                             #numericInput("n", "n", 50),
                             fileInput( "RNAfile", "Please Choose Your RNA-seq BAM Files", multiple = T),
                             fileInput( "RPFfile", "Please Choose Footprinting BAM Files", multiple = T),
                             fileInput( "ANNOfile", "Please Choose Annotation GFF3 File", multiple = T),
                             fileInput( "FASTAfile", "Please Choose FASTA File", multiple = T),
                             #checkboxInput("prepro", "Preprocess")
                             actionButton("prepro", "Preprocess",icon("paper-plane"), 
                                          style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                             #submitButton("Submit")
                           ),
                           mainPanel(
                             h3("Bringing Your Data to Life", align = "center"),
                             br(),
                             span("Welcome to ", strong("riboStreamR") , ", a platform for ribo-seq data analysis."),
                             p("Get started by uploading your files in the panel to the left and hitting ", span("Preprocess.", style = "color:blue")),
                             p("Once you recieve conifrmation below that your files have been successfully uploaded and processed, you can move onto viewing our collection of QC plots, available in the tabs above."),
                             hr(),
                             textOutput("upload"),
                             textOutput("process"),
                             textOutput("cheat"),
                             #plotOutput("plot"),
                             dataTableOutput(outputId = "starter")
                             
                             
                           )
                         )
                ),
                tabPanel("Read Length", #Read Length Distribution plot
                         sidebarLayout(
                           sidebarPanel(
                             uiOutput("choose_dataset1"), 
                             uiOutput("choose_dataset15"),                           
                             selectInput(inputId = "strand", label = "Choose the Strands", choices = c("+","-"),multiple = T, selected = c("+","-")),
                             selectInput(inputId = "mapping", label = "Choose which Mapped Types", choices = c("Multi","Unique"),multiple = T, selected = c("Multi","Unique")),
                             selectInput(inputId = "features", label = "Choose which Features", choices = c("mRNA", "miRNA", "ncRNA", "tRNA", "pseudogenic_transcript", "snoRNA", "snRNA", "rRNA", "promoter", "intron", "exon", "cds", "fiveUTR", "threeUTR", "intergenic"),multiple = T, selected = c("cds", "tRNA","rRNA")),
                             selectInput(inputId = "sepFeats", label = "Choose which Features to Separate", choices = c("mRNA", "miRNA", "ncRNA", "tRNA", "pseudogenic_transcript", "snoRNA", "snRNA", "rRNA", "promoter", "intron", "exon", "cds", "fiveUTR", "threeUTR", "intergenic"),multiple = T,selected = c("cds", "tRNA","rRNA")),
                             uiOutput("choose_gene"),   
                             selectInput(inputId = "frame", label = "Choose which Frames", choices = c("none", 0,1, 2),multiple = T,selected = c("none", 0,1,2)),
                             selectInput(inputId = "splitLine", label = "Choose the Factor to Separate Lines", choices = c("none","sample","strands", "mapping", "feature", "gene", "frame"),selected = "sample"),
                             selectInput(inputId = "splitPlot", label = "Choose the Factor to Separate Plots", choices = c("none","sample","strands", "mapping", "feature", "gene", "frame"),selected = "feature"),
                             sliderInput(inputId = "bw", label = "Choose Smoothing Bandwidth", min = 0.01, max = 2, value = 1),
                             sliderInput(inputId = "xmin", label = "Select X axis min", min = 0, max = 100, value = 0),
                             sliderInput(inputId = "xmax", label = "Select X axis max", min = 0, max = 100, value = 60),
                             selectInput("cols1", "Choose your color palette", choices = c("Spectral", "Blues", "BuGn", "GnBu", "Greens", "Greys", "Oranges", "Reds", "RdBu", "Set1"), multiple = F, selected = "Set1"),
                             actionButton("submit", "Submit"),
                             downloadButton("DownloadRld", "Download")
                           ),
                           mainPanel(
                             textOutput(outputId = "nulls"),
                             plotOutput(outputId = "dens", height = "600px")
                           )
                         )
                ),
                tabPanel("Components", #Feature component plot
                         sidebarLayout(
                           sidebarPanel(
                             uiOutput("choose_dataset2"),
                             uiOutput("choose_dataset25"),                             
                             selectInput(inputId = "strand2", label = "Choose the Strands", choices = c("+","-"),multiple = T, selected = c("+","-")),
                             selectInput(inputId = "mapping2", label = "Choose which Mapped Types", choices = c("Multi","Unique"),multiple = T, selected = c("Multi","Unique")),
                             selectInput(inputId = "features2", label = "Choose which Features", choices = c("mRNA", "miRNA", "ncRNA", "tRNA", "pseudogenic_transcript", "snoRNA", "snRNA", "rRNA", "promoter", "intron", "exon", "cds", "fiveUTR", "threeUTR", "intergenic"),multiple = T, selected = c("cds", "tRNA","rRNA")),
                             selectInput(inputId = "sepFeats2", label = "Choose which Features to Separate", choices = c("mRNA", "miRNA", "ncRNA", "tRNA", "pseudogenic_transcript", "snoRNA", "snRNA", "rRNA", "promoter", "intron", "exon", "cds", "fiveUTR", "threeUTR", "intergenic"),multiple = T,selected = c("cds", "tRNA","rRNA")),
                             sliderInput(inputId = "length", label = "Choose Read Length Range", value = c(20,40), min = 0, max = 100),
                             uiOutput("choose_gene2"),
                             selectInput(inputId = "frame2", label = "Choose which Frames", choices = c("none", 0,1, 2),multiple = T,selected = c("none", 0,1,2)),
                             selectInput(inputId = "splitBar", label = "Choose the Factor to Separate Bars", choices = c("none","sample","strands", "mapping", "feature", "widths","gene", "frame"),selected = "sample"),
                             selectInput(inputId = "splitCol", label = "Choose the Factor to Separate Colors", choices = c("none","sample","strands", "mapping", "feature", "widths","gene", "frame"),selected = "feature"),
                             selectInput(inputId = "splitGroup", label = "Choose the Factor to Separate Groups", choices = c("none","sample","strands", "mapping", "feature", "widths","gene", "frame"),selected = "none"),
                             selectInput("cols2", "Choose your color palette", choices = c("Spectral", "Blues", "BuGn", "GnBu", "Greens", "Greys", "Oranges", "Reds", "RdBu", "Set1"), multiple = F, selected = "Set1"),
                             
                             #selectInput(inputId = "splitBar", label = "Choose the Factor to Separate Bars", choices = c("none","sample","strands", "mapping", "feature", "widths","gene", "frame"),selected = "sample"),
                             #uiOutput(outputId = "splitCol"),
                             #uiOutput(outputId = "splitGroup"),
                             #uiOutput(outputId = "splitPlot2"),
                             actionButton("submit2", "Submit"),
                             downloadButton("DownloadComp", "Download")
                           ),
                           mainPanel(
                             plotOutput(outputId = "comp", height = "600px")
                           )
                         )
                ),
                tabPanel("GC", #GC% distribution plot
                         sidebarLayout(
                           sidebarPanel(
                             uiOutput("choose_dataset3"),
                             uiOutput("choose_dataset35"),                             
                             selectInput(inputId = "strand3", label = "Choose the strands", choices = c("+","-"),multiple = T, selected = c("+","-")),
                             sliderInput(inputId = "length3", label = "Choose Read Length Range", value = c(20,40), min = 0, max = 100),
                             selectInput(inputId = "mapping3", label = "Choose which Mapped Types", choices = c("Multi","Unique"),multiple = T, selected = c("Multi","Unique")),
                             selectInput(inputId = "features3", label = "Choose which Features", choices = c("mRNA", "miRNA", "ncRNA", "tRNA", "pseudogenic_transcript", "snoRNA", "snRNA", "rRNA", "promoter", "intron", "exon", "cds", "fiveUTR", "threeUTR", "intergenic"),multiple = T, selected = c("mRNA", "miRNA", "ncRNA", "tRNA", "pseudogenic_transcript", "snoRNA", "snRNA", "rRNA", "promoter", "intron", "exon", "cds", "fiveUTR", "threeUTR", "intergenic")),
                             selectInput(inputId = "sepFeats3", label = "Choose which Features to Separate", choices = c("mRNA", "miRNA", "ncRNA", "tRNA", "pseudogenic_transcript", "snoRNA", "snRNA", "rRNA", "promoter", "intron", "exon", "cds", "fiveUTR", "threeUTR", "intergenic"),multiple = T,selected = c("mRNA", "miRNA", "ncRNA", "tRNA", "pseudogenic_transcript", "snoRNA", "snRNA", "rRNA", "promoter", "intron", "exon", "cds", "fiveUTR", "threeUTR", "intergenic")),
                             uiOutput("choose_gene3"),  
                             selectInput(inputId = "frame3", label = "Choose which Frames", choices = c("none", 0,1, 2),multiple = T,selected = c("none", 0,1,2)),
                             selectInput(inputId = "splitLine3", label = "Choose the Factor to Separate Lines", choices = c("none","sample","strands", "mapping", "feature","gene", "frame", "widths"),selected = "sample"),
                             selectInput(inputId = "splitPlot3", label = "Choose the Factor to Separate Plots", choices = c("none","sample","strands", "mapping", "feature","gene", "frame", "widths"),selected = "none"),
                             selectInput("cols3", "Choose your color palette", choices = c("Spectral", "Blues", "BuGn", "GnBu", "Greens", "Greys", "Oranges", "Reds", "RdBu", "Set1"), multiple = F, selected = "Set1"),
                             actionButton("submit3", "Submit"),
                             downloadButton("DownloadGc", "Download")
                             
                           ),
                           mainPanel(
                             plotOutput(outputId = "gc", height = "600px")
                           )
                         )
                ),
                tabPanel("Summary Table", #Summary table for each sample
                         dataTableOutput(outputId = "table"),
                         #h3("Alerts", align = "left", style = "color:red"),
                         tableOutput(outputId = "alert"),
                         hr(),
                         fluidRow(
                           column(4,uiOutput("choose_dataset4"),
                                  uiOutput("choose_dataset45"),                             
                                  selectInput(inputId = "strand4", label = "Choose the Strands", choices = c("+","-"),multiple = T, selected = c("+","-")),
                                  selectInput(inputId = "mapping4", label = "Choose which Mapped Types", choices = c("Multi","Unique"),multiple = T, selected = c("Multi","Unique"))
                           ),
                           column(4, selectInput(inputId = "features4", label = "Choose which Features", choices = c("mRNA", "miRNA", "ncRNA", "tRNA", "pseudogenic_transcript", "snoRNA", "snRNA", "rRNA", "promoter", "intron", "exon", "cds", "fiveUTR", "threeUTR", "intergenic"),multiple = T, selected = c("mRNA", "miRNA", "ncRNA", "tRNA", "pseudogenic_transcript", "snoRNA", "snRNA", "rRNA", "promoter", "intron", "exon", "cds", "fiveUTR", "threeUTR", "intergenic")),
                                  sliderInput(inputId = "length4", label = "Choose Read Length Range", value = c(20,40), min = 0, max = 100),
                                  uiOutput("choose_gene4")
                           ),
                           column(4,selectInput(inputId = "frame4", label = "Choose which Frames", choices = c("none", 0,1, 2),multiple = T,selected = c("none",0,1,2)),
                                  selectInput(inputId = "cols", label = "Choose which Features Columns to Include", choices = c("sample","Aligned Reads","Unique%", "cds%", "rRNA%", "tRNA%", "3UTR%", "5UTR%", "intergenic%", "promoter%", "complexity", "periodicity"),multiple = T, selected = c("sample","Aligned Reads","Unique%", "cds%", "rRNA%", "tRNA%", "3UTR%", "5UTR%", "intergenic%", "intron%", "promoter%", "complexity", "periodicity")),
                                  actionButton("submit4", "Submit")
                           )
                         )
                ),
                tabPanel("Meta-gene", #Metagene read distributions around start and stop codons
                         sidebarLayout(
                           sidebarPanel(
                             uiOutput("choose_dataset5"),
                             uiOutput("choose_dataset55"),
                             selectInput(inputId = "strand5", label = "Choose the Strands", choices = c("+","-"),multiple = T, selected = c("+","-")),
                             selectInput(inputId = "mapping5", label = "Choose which Mapped Types", choices = c("Multi","Unique"),multiple = T, selected = c("Multi","Unique")),
                             sliderInput(inputId = "length5", label = "Choose Read Length Range", value = c(20,40), min = 0, max = 100),
                             uiOutput("choose_gene5"),
                             selectInput(inputId = "frame5", label = "Choose which Frames", choices = c("none", 0,1, 2),multiple = T,selected = c("none", 0,1,2)),
                             sliderInput(inputId = "window5", label = "Choose your window size", value = 100, min = 20, max = 300),
                             selectInput(inputId = "splitLine5", label = "Choose the Factor to Separate Lines", choices = c("none","sample","strands", "mapping", "gene", "frame", "widths"),selected = "sample"),
                             selectInput("cols4", "Choose your color palette", choices = c("Spectral", "Blues", "BuGn", "GnBu", "Greens", "Greys", "Oranges", "Reds", "RdBu", "Set1"), multiple = F, selected = "Set1"),
                             actionButton("submit5", "Submit"),
                             downloadButton("DownloadMeta", "Download")
                           ),
                           mainPanel(
                             plotOutput(outputId = "meta", height = "600px")
                           )
                         )
                ),
                tabPanel("Gene View", #read distributions for individual genes
                         sidebarLayout(
                           sidebarPanel(
                             uiOutput("choose_dataset6"),
                             uiOutput("choose_dataset65"),
                             selectInput(inputId = "strand6", label = "Choose the Strands", choices = c("+","-"),multiple = T, selected = c("+","-")),
                             selectInput(inputId = "mapping6", label = "Choose which Mapped Types", choices = c("Multi","Unique"),multiple = T, selected = c("Multi","Unique")),
                             sliderInput(inputId = "length6", label = "Choose Read Length Range", value = c(20,40), min = 0, max = 100),
                             uiOutput("choose_gene6"),
                             selectInput(inputId = "frame6", label = "Choose which Frames", choices = c("none", 0,1, 2),multiple = T,selected = c("none", 0,1,2)),
                             selectInput(inputId = "splitLine6", label = "Choose the Factor to Separate Lines", choices = c("none","sample","strands", "mapping", "frame", "widths"),selected = "sample"),
                             selectInput(inputId = "splitPlot6", label = "Choose the Factor to Separate Plots", choices = c("none","sample","strands", "mapping", "frame"),selected = "frame"),
                             actionButton("submit6", "Submit"),
                             downloadButton("DownloadGview", "Download")
                             
                           ),
                           mainPanel(
                             plotOutput(outputId = "gview", height = "800px")
                           )
                         )
                ),
                tabPanel("Meta-Periodicity",
                         sidebarLayout(
                           sidebarPanel(
                             uiOutput("choose_dataset7"),
                             uiOutput("choose_dataset75"),
                             selectInput(inputId = "strand7", label = "Choose the Strands", choices = c("+","-"),multiple = T, selected = c("+","-")),
                             selectInput(inputId = "mapping7", label = "Choose which Mapped Types", choices = c("Multi","Unique"),multiple = T, selected = c("Multi","Unique")),
                             sliderInput(inputId = "length7", label = "Choose Read Length Range", value = c(20,40), min = 0, max = 100),
                             uiOutput("choose_gene7"),
                             #selectInput(inputId = "frame7", label = "Choose which Frames", choices = c("none", 0,1, 2),multiple = T,selected = c("none", 0,1,2)),
                             sliderInput(inputId = "window7", label = "Choose your window size", value = 100, min = 20, max = 300),
                             selectInput(inputId = "splitPlot7", label = "Choose the Factor to Separate Plots", choices = c("none","sample","strands", "mapping", "gene", "frame", "widths"),selected = "sample"),
                             selectInput("cols5", "Choose your color palette", choices = c("Spectral", "Blues", "BuGn", "GnBu", "Greens", "Greys", "Oranges", "Reds", "RdBu", "Set1"), multiple = F, selected = "Set1"),
                             actionButton("submit7", "Submit"),
                             downloadButton("DownloadMetap", "Download")
                             
                           ),
                           mainPanel(
                             plotOutput(outputId = "meta_perio", height = "600px"),
                             plotOutput("offsetPlot")
                           )
                         )
                ),
                tabPanel("Length Periodicity",
                         sidebarLayout(
                           sidebarPanel(
                             uiOutput("choose_dataset8"),
                             uiOutput("choose_dataset85"),
                             selectInput(inputId = "strand8", label = "Choose the Strands", choices = c("+","-"),multiple = T, selected = c("+","-")),
                             selectInput(inputId = "mapping8", label = "Choose which Mapped Types", choices = c("Multi","Unique"),multiple = T, selected = "Unique"),
                             selectInput(inputId = "features8", label = "Choose which Features", choices = c("mRNA", "miRNA", "ncRNA", "tRNA", "pseudogenic_transcript", "snoRNA", "snRNA", "rRNA", "promoter", "intron", "exon", "cds", "fiveUTR", "threeUTR", "intergenic"),multiple = T, selected = "cds"),
                             sliderInput(inputId = "length8", label = "Choose Read Length Range", value = c(25,31), min = 0, max = 100),
                             selectInput(inputId = "splitGroup8", label = "Choose how to separate groups",choices = c("none","sample","strands", "mapping", "feature", "widths", "gene", "frame"),selected = "sample"),
                             uiOutput("choose_gene8"),
                             selectInput("cols6", "Choose your color palette", choices = c("Spectral", "Blues", "BuGn", "GnBu", "Greens", "Greys", "Oranges", "Reds", "RdBu", "Set1"), multiple = F, selected = "Set1"),
                             selectInput(inputId = "frame8", label = "Choose which Frames", choices = c("none", 0,1, 2),multiple = T,selected = c(0,1,2)),
                             actionButton("submit8", "Submit"),
                             downloadButton("DownloadPerio", "Download")
                           ),
                           mainPanel(
                             plotOutput(outputId = "perio", height = "600px"),
                             plotOutput("offsetPlot2")
                             
                           )
                         )
                ),
                tabPanel("Differential",
                         sidebarLayout(
                           sidebarPanel(
                             selectInput("select", "Type of differential analysis:", choices = c("Expression","Translation", "Efficiency"), multiple = F, selected = "Translation"),
                             selectInput("expers", label = "How many different Experiment Types",choices = c(1,2), selected = 2),
                             numericInput("treats", label = "How many different Treatments?", value = 2),
                             #actionButton("update", "Update UI"),
                             uiOutput("inputvals1"),
                             uiOutput("inputvals2"),
                             #selectInput(inputId = "RNAtreat1", label = "Choose the RNA seq files for treatment 1", choices = list, multiple = T, selected = list[1:2]),
                             #selectInput(inputId = "RPFtreat1", label = "Choose the Ribo-seq files for treatment 1", choices = list, multiple = T, selected = list[5:6]),
                             #selectInput(inputId = "RNAtreat2", label = "Choose the RNA seq files for treatment 2", choices = list, multiple = T, selected = list[3:4]),
                             #selectInput(inputId = "RPFtreat2", label = "Choose the Ribo-seq files for treatment 2", choices = list, multiple = T, selected = list[7:8]),
                             #submitButton("Submit"),
                             downloadButton('downloadCount', 'Download Counts'),
                             downloadButton('downloadRPKM', 'Download RPKMs'),
                             downloadButton('downloadDiff', 'Download edgeR'),
                             actionButton("submit9", "Submit")
                             
                           ),
                           mainPanel(
                             dataTableOutput(outputId = "countTable"),
                             plotOutput("subseq")
                           )
                         )
                ),
                tabPanel("Report",
                         sidebarLayout(
                           sidebarPanel(
                             checkboxGroupInput("plots", "Choose Plots", choices = c("Summary Table", "Read Length Distribution", "GC Distribution", "Feature Components", "Meta Gene", "Single Gene", "Meta periodicity", "Periodicity")),
                             downloadButton("gen", "Generate",icon("paper-plane"), 
                                            style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                             
                           ),
                           mainPanel(
                             h3("QC Report Generator", align = "center"),
                             br(),
                             span("Choose which plots to include in your report using the panel to the left and press ", span("Generate ", style = "color:blue"), p("to create the report.")),
                             hr()
                           )
                         )
                         
                )
                
))


