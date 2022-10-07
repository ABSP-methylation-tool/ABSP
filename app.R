#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#



#────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

# Here, to simplify, a short list of genomes is displayed but all BSgenome can be used
# To get the list of all genomes of BSgenome package run : 'BSgenome::available.genomes()', more information on genomes at https://genome.ucsc.edu/cgi-bin/hgGateway
# A new genome can be added to the list displayed just below :
list_genomes <- c(
    "BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Mmusculus.UCSC.mm10", "BSgenome.Mmusculus.UCSC.mm39", 
    "BSgenome.Rnorvegicus.UCSC.rn6", "BSgenome.Rnorvegicus.UCSC.rn7", "BSgenome.Cfamiliaris.UCSC.canFam3", "BSgenome.Mmulatta.UCSC.rheMac8",
    "BSgenome.Ggallus.UCSC.galGal6" , "BSgenome.Drerio.UCSC.danRer11", "BSgenome.Celegans.UCSC.ce11", "BSgenome.Dmelanogaster.UCSC.dm6")

#────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────


#────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

# Install packages 
packages <- c("arrangements","BiocManager","compareGroups","DiagrammeR","dplyr","formattable","GenomeInfoDb",
              "ggdendro","ggplot2","ggpubr","htmltools","htmlwidgets","knitr","openxlsx","pdftools","plotly",
              "png","purrr","RColorBrewer","readr","rlist","rmarkdown","Rmisc","rstatix","seqinr","shiny",
              "shinythemes","webshot")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(new.packages)>0) {
  install.packages(new.packages)
}
# Install packages from Bioconductor
Biocpackages <-c("BiocGenerics","Biostrings","BSgenome","GenomicRanges","Gviz","sangeranalyseR","sangerseqR")
new.Biocpackages <- Biocpackages[!(Biocpackages %in% installed.packages()[,"Package"])]
if (length(new.Biocpackages)>0) {
  BiocManager::install(new.Biocpackages)
}

# Add packages to library
suppressWarnings(lapply(packages, library, character.only = TRUE)) #, quietly = TRUE
suppressWarnings(lapply(Biocpackages, library, character.only = TRUE))

#────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────





ui <- fluidPage(
    
    # Theme
    theme = shinytheme("flatly"),
    includeCSS("www/custom_app.css"),
    
    
    titlePanel(windowTitle = "ABSP",
               title = 
                   div(img(src = "logo.svg", style = "height:60px ; padding-right:20px"),
                       span("Analysis ", style="color:#40b7a0"), span("of Bisulfite Sequencing PCR", style="color:#1E2934"))
    ),
    
    
    navbarPage("",
               
               #─────────────────────────────────────────────────────────────────────────────────────────────────────────
               # MAIN PAGE ----------------------------------------------------------------------------------------------
               
               tabPanel("Main",
                        
                        column(width=5, 
                               
                               div(
                                   h3("What is ABSP ?"),
                                   
                                   br(),
                                   
                                   p("ABSP is a R based tool to analyze Bisulfite Sequencing PCR (BSP) experiment results, which is why it is the acronym for 'Analysis of Bisulfite Sequencing PCR'."), 
                                   p("It was developed to help researchers estimate and compare CpG methylation percentages of a DNA region studied using BSP experiments."),
                                   p("It provides a complete automated workflow, from trace file sequencing results to data visualization and statistics."),
                                   
                                   br()
                                   
                               ),
                               
                               
                               div(class="container", style="width: 100%; display:table; padding:0px;",
                                   
                                   div(class="box", style="display:table;",
                                       
                                       div(class="box-row", style="display:table-row;",
                                           
                                           div(class="box-cell box1",
                                               span("Complete workflow from raw data to statistics", style="vertical-align: middle;"),
                                               style="color:#ffffff; background-color:#2c3e50; border:5px solid white; border-radius: 12px;
                                               vertical-align: middle; text-align: center; display:table-cell; width:25%; padding: 10px;"
                                           ),
                                           
                                           div(class="box-cell box2",
                                               span("Fully automated and user-friendly", style="vertical-align: middle;"),
                                               style="color:#ffffff; background-color:#2c3e50; border:5px solid white; border-radius: 12px;
                                               vertical-align: middle; text-align: center; display:table-cell; width:25%; padding: 10px;"
                                           ),
                                           
                                           div(class="box-cell box3",
                                               span("Analysis of Direct-BSP and Cloning-BSP", style="vertical-align: middle;"),
                                               style="color:#ffffff; background-color:#2c3e50; border:5px solid white; border-radius: 12px;
                                               vertical-align: middle; text-align: center; display:table-cell; width:25%; padding: 10px;"
                                           ),
                                           
                                           div(class="box-cell box4",
                                               span("Accessibility and flexibility using the R language", style="vertical-align: middle;"),
                                               style="color:#ffffff; background-color:#2c3e50; border:5px solid white; border-radius: 12px;
                                               vertical-align: middle; text-align: center; display:table-cell; width:25%; padding: 10px;"
                                           ) 
                                       )
                                   )
                               ),
                               
                               
                               br(),
                               
                               div(
                                 
                                 h3("Please, cite"),
                                 
                                 br()
                                 
                                 
                               ),
                               
                               div(
                                 
                                 p(icon(name="newspaper", lib = "font-awesome", style="padding-right:8px; font-size:1.2em;"),"Denoulet ", em("et al."),", 2022"),
                                 
                                 br(),
                                 br()
                                 
                               ),
                               
                               div(
                                 
                                 h3("Ressources"),
                                 
                                 br(),
                                 
                                 p(icon(name="download", lib = "font-awesome", style="padding-right:8px; font-size:1.2em;"),"The ABSP tool is available for download on github at", htmltools::tags$a(href="https://github.com/ABSP-methylation-tool/ABSP", "https://github.com/ABSP-methylation-tool/ABSP",target="_blank"),"."),
                                 p(icon(name="file", lib = "font-awesome", style="padding-right:8px; font-size:1.2em;"),"For detailed instructions, please find the user guide in your ABSP folder."),
                                 
                                 br(),
                                 br()
                                 
                               ),
                               
                               div(
                                 
                                 h3("Contact"),
                                 
                                 br(),
                                 
                                 p(icon(name="envelope", lib = "font-awesome", style="padding-right:8px; font-size:1.2em;"), htmltools::tags$a(href="mailto:absp@univ-lille.fr","absp@univ-lille.fr"))
                               )
                        ),
                        
                        
                        column(width=7,
                               
                               
                               column(width = 4,
                                      div(img(src="ABSP - BSP.svg", style = "width:100%;"), 
                                          style="text-align:center")
                               ),
                               
                               column(width = 4,
                                      div(img(src="ABSP - Analysis.svg", style = "width:100%;"), 
                                          style="text-align:center")
                               ),
                               
                               column(width = 4,
                                      div(img(src="logo.svg", style = "width:50%;"), 
                                          style="text-align:center"),
                                      div(img(src="ABSP - Workflow simple.svg", style = "width:100%;"), 
                                          style="text-align:center")
                               )
                        )
               ),
               
    
               #─────────────────────────────────────────────────────────────────────────────────────────────────────────
               # TAB PANEL INDIVIDUAL ANALYSIS --------------------------------------------------------------------------
               
               tabPanel("Individual analysis",
                        
                        #─────────────────────────────────────────────────────────────────────────────────────────────────────────
                        # Side bar indiv -----------------------------------------------------------------------------------------
                    
                        sidebarPanel(width=5,
                            
                            h4("Sample information", style="color:#374e64 ; font-weight: bold ; font-size:18px"),
                            
                            fluidRow(
                                column(width = 7,
                                       # input$folderexisting
                                       selectInput("folderexisting", label = "Select existing folder",
                                                   choices = c(basename(list.dirs(file.path(getwd(),"results"), recursive = F)),"Create new folder")
                                                   )),
                                column(width = 5, 
                                       # input$foldernew
                                       uiOutput("foldernew"))
                            ),
                            
                            fluidRow(
                                column(width = 7,
                                       # input$seqexisting
                                       selectInput("seqexisting", label = "Select existing sequence folder",
                                                   choices = "Create new sequence folder")),
                                column(width = 5, 
                                       # input$seqnew
                                       uiOutput("seqnew"))
                            ),   
                            
                            # input$collection
                            textInput("collection","Enter collection name","collection1"),
                            
                            # input$group
                            textInput("group","Enter group name","group1"),
                            
                            
                            fluidRow(
                                column(width = 6,
                                       # input$replicate
                                       numericInput("replicate", span("Direct-BSP",br(),"Enter replicate number"), "1")),
                                
                                column(width = 6,
                                       # input$clone
                                       numericInput("clone",span("Cloning-BSP",br(),"Enter clone number"), "0"))
                            ),
                            
                            br(),
                            
                            h4("Reference DNA sequence", style="color:#374e64 ; font-weight: bold ; font-size:18px"),
                            
                            fluidRow(
                                column(width = 8,
                                       # input$genomeI
                                       selectInput("genomeI",label ="Select genome",choices=list_genomes,selected="BSgenome.Hsapiens.UCSC.hg19")
                                ),
                                
                                column(width = 4,
                                       # input$install_genomeI : Pre install genome
                                       actionButton("install_genomeI","Pre-install genome", 
                                                    icon=shiny::icon(name="download", lib = "font-awesome", style="padding-right:8px;"),
                                                    style="margin-top:25px; white-space: normal; padding:5px;"),
                                       uiOutput("installed_genomeI")
                                )
                                
                            ),
                                       
                            # input$DNA_seq
                            fileInput("DNA_seq","Select .fasta file of reference DNA sequence"), # from plus strand with header containing coordinates and strand (plus or minus) used for primer design
                            
                            br(),
                            
                            h4("Sequencing results", style="color:#374e64 ; font-weight: bold ; font-size:18px"),
                            

                            fluidRow(
                                column(width = 6,
                                       # input$date_s1
                                       dateInput("date_s1","Select date of sequencing #1",Sys.Date()),
                                       
                                       
                                       fluidRow(
                                                column(width = 9,
                                                       fileInput("ab1_s1","Select sequencing file #1")
                                                ),
                                                column(width = 3,
                                                       # input$reset_s1
                                                       actionButton("reset_s1", "Reset",  style="padding:10px; margin-top:25px;") 
                                                )
                                                
                                       ),
                                       
                                       fluidRow(
                                         div(verbatimTextOutput("summary_ab1_s1"), style="padding: 0px 15px 15px 15px;")
                                       )
                                ),
                                
                                column(width = 6,
                                       # input$date_s2
                                       dateInput("date_s2","Select date of sequencing #2",Sys.Date()),
                                       
                                       
                                       fluidRow(
                                                column(width = 9,
                                                       fileInput("ab1_s2","Select sequencing file #2")
                                                       #fileInput("ab1_s2","Select sequencing file #2")
                                                ),
                                                column(width = 3,
                                                       # input$reset_s2
                                                       actionButton("reset_s2", "Reset",  style="padding:10px; margin-top:25px;") 
                                                )
                                                
                                       ),
                                       
                                       fluidRow(
                                         div(verbatimTextOutput("summary_ab1_s2"), style="padding: 0px 15px 15px 15px;")
                                       )
                                       
                                )
                            ),
                            
                            
                            
                            #─────────────────────────────────────────────────────────────────────────────────────────────────────────
                            # Action button indiv ------------------------------------------------------------------------------------
                            
                            # input$indiv
                            fluidRow( 
                                align="center",
                                actionButton("indiv","Run individual analysis", 
                                             icon=shiny::icon("play", lib = "font-awesome", style="padding-right:8px ;"),
                                             style="background-color: #40b7a0; border-color:#40b7a0; font-size:17px;")
                            ),
                            
                            uiOutput("indivReport")
                            
                        ),
                        
                        #─────────────────────────────────────────────────────────────────────────────────────────────────────────
                        # Main panel indiv ---------------------------------------------------------------------------------------
                        
                        mainPanel(width=7,
                            
                            
                            
                            tabsetPanel(
                                type = "tabs",
                                
                                #─────────────────────────────────────────────────────────────────────────────────────────────────────────
                                # Tabset panel - experiments infos -----------------------------------------------------------------------
                                
                                
                                tabPanel(
                                    
                                    # Tab title
                                    span("Sample information", style="font-size:18px"),
                                    
                                    
                                    br(),
                                    
                                    # WARNINGS
                                    
                                    h5("Be aware that these entries must not contain any special characters :", 
                                       htmltools::tags$code("/ \\ : * ? ! \" ' ` < > | & % @ # + = { }", style="color:grey ; background-color:transparent"),
                                       style="color:#bf3232"),
                                    
                                    h5("Allowed characters (that must be avoided if not necessary) :", 
                                       htmltools::tags$code(" space - _ ", style="color:grey ; background-color:transparent"), 
                                       style="color:#bf3232"),
                                    
                                    h5("All entries must be named consistently between the different analysis 
                                                    (e.g. the group names must strickly be identical between samples of the same group).", 
                                       style="color:#bf3232"),
                                    
                                    hr(),
                                    
                                    div( # FOLDER
                                        h4("Select existing folder", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("Select an existing folder within the ABSP results folder to locate all of the analysis results.
                                                       Having different folders of results can be use to separate the different analysis by projects, 
                                                       experiments or users."),
                                        p("To create a new folder, select the 'Create new folder' entry and enter the name 
                                                       of the new folder in the text input. Note that the six first letters will appear in the report file name.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # SEQUENCE
                                        
                                        h4("Select existing sequence folder", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("Select an existing folder corresponding to the name of your sequence.
                                           This is the name that will be used on tables and plots to refer to the sequence."),
                                        p("To analyze a new sequence and therefore create a new sequence folder,
                                          select the 'Create new sequence folder' entry and enter the new sequence name in the text input.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # COLLECTION
                                        h4("Enter collection name", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("The collection corresponds to a separation of samples above groups. 
                                          Samples from different collections can not be compared, even if they belong to the same group. 
                                          For example, collections can be different cell lines, organs, or patients, in which the same 
                                          groups are compared but not between the different collections. 
                                          To compare these types of samples, consider them as groups."),
                                        em("If you do not want to specify any collection, leave empty or enter '0'. 
                                           Make sure the collection name is strictly identical for all samples of the same collection.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # GROUP
                                        h4("Enter group name", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("The group corresponds to the condition you want to compare in the grouped analysis. 
                                          For example, groups can be the 'control' and 'treated' conditions."),
                                        em("The group field entry is required. Make sure the group name is strictly identical for all samples of the same group.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # REPLICATE
                                        h4("Direct-BSP - Enter replicate number", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("In the case of direct sequencing of PCR products only (direct-BSP). 
                                          The replicate number corresponds to the number of repetition for one group, for experiment reproducibility 
                                          and statistic significance determination. A minimum of 3 replicates is generally recommended.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # CLONE
                                        h4("Cloning-BSP - Enter clone number", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("In the case of clone sequencing only (cloning-BSP). 
                                          The clone number corresponds to the identification number of each clone for one group. 
                                          In general, the sequencing of 10 clones is recommended for results reliability.")
                                    )
                                ),
                                
                                #─────────────────────────────────────────────────────────────────────────────────────────────────────────
                                # Tabset panel - Ref DNA sequence ------------------------------------------------------------------------
                                
                                tabPanel(
                                    # Tab title
                                    span("Reference DNA sequence", style="font-size:18px"),
                                    
                                    br(),
                                    
                                    div( # GENOME
                                        h4("Select genome", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("Select the reference genome. Only used to display the genomic sequence in the genomic plot."),
                                        p("Make sure to click on 'Pre-install genome' button if the selected genome is used for the first time.", 
                                          style="color:#bf3232"),
                                        em("Only a few of the available genomes are listed in the selection tab, but more genomes can be used. 
                                           If your studied genome does not appear please refer to the user guide.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # FASTA
                                        h4("Select .fasta file of reference DNA sequence", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("Select the file '.fasta' in your folders to provide the reference DNA sequence from the plus strand of genome.
                                                       The fasta file must contain in the header two elements :"),
                                        htmltools::tags$ul(
                                          htmltools::tags$li("The precise genomic coordinates of the sequence contained in the file,in the strict format : ", br(),
                                                    span("chr#:######-###### (e.g. chr16:68771087-68771462)", style="color:#40b7a0")),
                                            
                                          htmltools::tags$li("The strand chosen for primer design : the strand complementary to primers after bisulfite 
                                                                 conversion, in the strict format : ", br(),
                                                    span("primers=plus or primers=minus", style="color:#40b7a0"))
                                        ),
                                        em("Note that any other information in the header, as the sequence name for example, can be added 
                                                        without consequences if they do not interfere with the previously described formats."),
                                        
                                        br(),
                                        
                                        h5("Example of reference DNA fasta file format :"),
                                        div(img(src="ABSP - fasta file.png")) 
                                    )
                                ),
                                
                                #─────────────────────────────────────────────────────────────────────────────────────────────────────────
                                # Tabset panel - sequencing results ----------------------------------------------------------------------
                                
                                tabPanel(
                                    # Tab title
                                    span("Sequencing results", style="font-size:18px"),
                                    
                                    br(),
                                    
                                    div( # DATE
                                        h4("Select dates of sequencing", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("In a matter of traceability, you can select the dates when the sequencing runs were performed."),
                                        em("If one of the sequencing files is not provided, the date entry must be empty.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # SEQ 
                                        h4("Select .ab1 sequencing files", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("Select the files in your folders corresponding to the .ab1 sequencing files. 
                                          The sequencing reads from one direction (#1) and the other direction (#2) should both be provided, although the analysis can be runned with only one sequencing read provided."),
                                        p("The directions (forward or reverse) will be determined during the analysis.")
                                    )
                                )
                            )
                         )
               ),
               
               #─────────────────────────────────────────────────────────────────────────────────────────────────────────
               # TAB PANEL GROUPED ANALYSIS -----------------------------------------------------------------------------
               
               tabPanel("Grouped analysis",
                        
                        #─────────────────────────────────────────────────────────────────────────────────────────────────────────
                        # Side bar grouped ---------------------------------------------------------------------------------------
                        
                        sidebarPanel(width=5,
                            
                            h4("Experiment information", style="color:#374e64 ; font-weight: bold ; font-size:18px"),
                            
                            # input$folderselect
                            selectInput("folderselect",label ="Select folder",
                                        choices=basename(list.dirs(file.path(getwd(),"results"), recursive = F)),
                                        selected=basename(list.dirs(file.path(getwd(),"results"), recursive = F))[1]),
                            
                            # input$seqselect
                            selectInput("seqselect",label ="Select sequence",
                                        choices="Select folder before"),
                            
                            # input$genomeG
                            selectInput("genomeG",label ="Select genome",choices=list_genomes,selected="BSgenome.Hsapiens.UCSC.hg19"),
                            
                            # input$exptype
                            radioButtons("exptype",label ="Select experiment type",choices=c("Direct-BSP","Cloning-BSP"), inline=T),
                            
                            br(),
                            
                            h4("Plot parameters", style="color:#374e64 ; font-weight: bold ; font-size:18px"),
                            
                            # input$pos_labels
                            radioButtons("pos_labels",label ="Select position labels for plots",
                                         choices=c("CpG coordinates","CpG numbers","None"), inline=T),
                            
                            p("Collection separation", style="font-weight: bold ; font-size:16px"), 
                            
                            # input$coll_sep
                            checkboxInput("coll_sep",label="Separate plots by collection", value=F),
                            
                            # # input$group_order
                            # textInput("group_order",label ="Indicate groups to compare in the correct order for display (comma separated)",
                            #           "group1, group2, group3"),
                            
                            # input$group_order
                            selectInput("group_order",label ="Indicate the order of groups for display (all must be selected)",
                                      choices="Select folder and sequence first to display groups (or modify experiment type)", multiple=T),
                            
                            # input$sample_order
                            selectInput("sample_order",label = "Select the types of sample ordering for plots (multiple choice allowed)",
                                        choices = c("As it is", "By groups", "By methylation levels", "By clusters"), multiple=T),
                            
                            
                            #─────────────────────────────────────────────────────────────────────────────────────────────────────────
                            # Action button grouped ----------------------------------------------------------------------------------
                            
                            # input$grouped
                            fluidRow( 
                                align="center",
                                actionButton("grouped","Run grouped analysis", 
                                             icon=shiny::icon("play", lib = "font-awesome", style="padding-right:8px ;"),
                                             style="background-color: #40b7a0; border-color: #40b7a0; font-size:17px;")
                            ),
                            
                                uiOutput("groupedReport")
                        ),
                        
                        #─────────────────────────────────────────────────────────────────────────────────────────────────────────
                        # Main panel grouped -------------------------------------------------------------------------------------
                        
                        mainPanel(width=7,
                            
                            tabsetPanel(
                                type = "tabs",
                                
                                #─────────────────────────────────────────────────────────────────────────────────────────────────────────
                                # Tabset panel - sample infos ----------------------------------------------------------------------------
                                
                                tabPanel(
                                    # Tab title 
                                    span("Experiment information", style="font-size:18px"),
                                    
                                    br(),
                                    
                                    div( # FOLDER
                                        h4("Select folder", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("Select an existing folder within the ABSP results folder to retrieve the previously 
                                                       generated methylation data files after individual analyses, and to locate the newly 
                                                       generated results of the grouped analysis.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # SEQUENCE
                                        h4("Select sequence", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("Select an existing folder corresponding to the name of your sequence, to retrieve the previously 
                                                       generated methylation data files after individual analyses, and to locate the newly 
                                                       generated results of the grouped analysis.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # GENOME
                                        h4("Select genome", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("Select the reference genome. Only used to display the genomic sequence in the genomic plots.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # EXP
                                        h4("Select experiment type", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("Select the experiment type, either 'Direct-BSP' or 'Cloning-BSP', for proper 
                                                       retrieval of methylation data and analysis specification.") 
                                    )
                                    
                                ),
                                
                                #─────────────────────────────────────────────────────────────────────────────────────────────────────────
                                # Tabset panel - plot parameters -------------------------------------------------------------------------
                                
                                tabPanel(
                                    # Tab title
                                    span("Plot parameters", style="font-size:18px"),
                                    
                                    br(),
                                    
                                    div( # POS LABS
                                        h4("Select position labels for plots", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("Select which label type you want to display as the CpG positions in plots, among 'CpG coordinates', 'CpG numbers' and 'None'.
                                               Note that in case of extremely close CpG positions, labels may overlap, therefore the 'None' label type can provide a good alternative.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # COLL SEP
                                        h4("Collection separation", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("Parameter to completely separate collections for display in plots or not, in both cases the groups of different collection are not compared.
                                          If unchecked, plots will represents all samples of all collections, and if checked each collection will be represented in a separated plot.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # GROUP ORDER
                                        h4("Indicate order of groups for display", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("The names of groups are retrieved from files in the selected folder above. 
                                          Select them all in the correct order, chosen to display samples in plots.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # SAMPLE ORDER
                                        h4("Select the types of sample ordering for plots", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("Choose from one to four of the sample ordering available for plots :"),
                                        htmltools::tags$ul(
                                          htmltools::tags$li("'As it is' arranges samples by alphabetic order of collections. 
                                                    If none or one collection are present, this order is equivalent as the 'By groups' one."),
                                          htmltools::tags$li("'By groups' arranges samples by the provided group order above."),
                                          htmltools::tags$li("'By methylation levels' arranges samples depending on their methylation mean."),
                                          htmltools::tags$li("'By clusters' arranges samples depending on the hierarchical clustering calculated and represented by an associated dendrogram.")
                                        )
                                    )
                                )
                            )
                        )
               ),
               
               
               #─────────────────────────────────────────────────────────────────────────────────────────────────────────
               # TAB PANEL MULTIPLE ANALYSIS ----------------------------------------------------------------------------

               tabPanel("Multiple analyses",
                        
                        sidebarPanel(width=5,
                            
                            
                            fluidRow(
                                column(width = 7,
                                       
                                       # input$foldermainexisting
                                       selectInput("foldermainexisting", label = "Select existing folder",
                                                   choices = c(basename(list.dirs(file.path(getwd(),"results"), recursive = F)),"Create new folder")
                                       )),
                                column(width = 5, 
                                       # input$foldermainnew
                                       uiOutput("foldermainnew"))
                            ),
                            
                            fluidRow(
                                
                                column(width=9, 
                                       fileInput("exptable","Select individual analyses table file")
                                       
                                ),
                                column(width=3,       
                                       # input$reset_exptable
                                       actionButton("reset_exptable", "Reset", style="margin-top:25px;")
                                ),
                                
                                
                            ),
                            
                            fluidRow(
                                
                                column(width=9, 
                                       verbatimTextOutput("summary_exptable")
                                )
                                
                            ),
                            
                            br(),
                            
                            fluidRow(
                                
                                column(width=9,
                                       fileInput("groupedparams","Select grouped analyses table file")
                                       
                                ),
                                column(width=3,       
                                       # input$reset_groupedparams
                                       actionButton("reset_groupedparams", "Reset", style="margin-top:25px;")
                                )
                                
                            ),
                            
                            fluidRow(
                                
                                column(width=9, 
                                       verbatimTextOutput("summary_groupedparams")
                                )
                                
                            ),
                            
                            
                            
                            #─────────────────────────────────────────────────────────────────────────────────────────────────────────
                            # Action button multiple ---------------------------------------------------------------------------------
                            

                            # input$multiple
                            fluidRow( 
                                align="center",
                                actionButton("multiple"," Run analyses", 
                                             icon=shiny::icon("play", lib = "font-awesome", style="padding-right:8px ;"),
                                             style="background-color: #40b7a0; border-color: #40b7a0; font-size:17px;")
                            ), 
                            
                            fluidRow(
                                align="center",
                                uiOutput("multipleIReport"),
                                uiOutput("multipleGReport")
                            )
                            
                        ),
                        
                        
                        mainPanel(width=7,
                                  
                                  tabsetPanel(
                                      type = "tabs",
                                      
                                      tabPanel(
                                          # Tab title 
                                          span("Principle", style="font-size:18px"),
                                          
                                          br(),
                                          
                                          h4("Automated launch of multiple analyses"),
                                          p("Several analyses can be launched from this tab, in only one click, using pre-filled tables containing the input entries."),
                                          
                                          br(),
                                          
                                          h4("How to proceed ?"),
                                          p("1. Fill one or the two input tables below with the sample information and the choice of parameters 
                                            (you can open the tables files with the buttons below)."),
                                          
                                          p("2. In the left panel, select an existing folder within the ABSP results folder to locate all of the analyses results."),
                                          p("To create a new folder, select the 'Create new folder' entry and enter the name of the new folder in the text input. Note that the six first letters will appear in the report file name."),
                                          
                                          
                                          p("3. In the left panel, select your filled table as input."), 
                                          p("Both the sample data table for multiple individual analyses and the parameters table for multiple grouped analyses can be provided at the same time to launch the individual analyses followed by the grouped analyses,
                                            or only one of the two tables can be provided and will launch the corresponding analyses, either individual analyses or grouped analyses."),
                                          
                                          p("4. Launch the analyses by clicking on the bottom button 'Run analyses'. "),
                                          p("All the input analyses will be launch one after the other, do not close the app until the end."),
                                        
                                          em("Note: The tables files are located in the 'documents' folder. The input file format must be either '.xlsx' or '.csv'."),
                                          
                                          hr(),
                                          
                                          h4("Sample data for individual analyses"),
                                          p("The table regroups the individual analyses inputs."),
                                          
                                          actionButton("exp_tab_open"," Open multiple individual analysis table file", 
                                                       icon=shiny::icon("folder-open", lib = "font-awesome", style="padding-right:8px ;"),
                                                       style="font-size:17px; margin:10px;"),
                                          
                                          hr(),
                                          
                                          h4("Parameters for the grouped analyses"),
                                          p("The table regroups the grouped analyses inputs."),
                                          
                                          actionButton("group_tab_open"," Open multiple grouped analyses table file", 
                                                       icon=shiny::icon("folder-open", lib = "font-awesome", style="padding-right:8px ;"),
                                                       style="font-size:17px;")
                                          
                                          
                                      ),
                                      
                                      tabPanel(
                                          # Tab title 
                                          span("Diagram of analysis launch", style="font-size:18px"),
                                          
                                          div(img(src="ABSP - Launch analysis.svg", style = "width:80% ; padding:10px"), 
                                              style="text-align:center")
                                      ) 
                                  )     
                        )
               )
    ),
    
    br(),
    br(),
    br(),
    br(),
    
    #─────────────────────────────────────────────────────────────────────────────────────────────────────────
    # Close button -------------------------------------------------------------------------------------------
    
    div(class="navbar navbar-default navbar-fixed-bottom",
        style="min-height: 40px;",
        fluidRow(
            div(style="padding-right: 15px; padding-left: 15px;", 
                
                column(width=6,
                       p("ABSP v1.0.0 - Copyright © 2022 CANTHER laboratory, released under the GPL-3 license", style="font-size:14px ; color:#ffffff ; padding-top:15px")),
                
                column(width=6,
                       htmltools::tags$button(
                           id = 'close',
                           type = "button",
                           class = "btn btn-default action-button btn-primary navbar-btn pull-right shiny-bound-input",
                           style = "background-color: #40b7a0; border-color: #40b7a0 ; float: right!important; margin-right:10px; border-width: 0px; padding: 6px 15px;",
                           onclick = "setTimeout(function(){window.close();},500);",  # close browser
                           "Close"
                       ))
            )
            
        )
        
        
    )
    
    
    
)



##############################################################################################################



server <- function(input, output, session) {
  
  
  #─────────────────────────────────────────────────────────────────────────────────────────────────────────
  # Close button -------------------------------------------------------------------------------------------
  observe({
    if (input$close > 0) stopApp() # stop shiny
  })
  # Close app when browser page is closed
  session$onSessionEnded(stopApp)
  
  
  
  #─────────────────────────────────────────────────────────────────────────────────────────────────────────
  # Multiple analyses tables to open -----------------------------------------------------------------------
  observeEvent(input$exp_tab_open,{
    shell.exec(file.path(getwd(),"documents","multiple_individual_analyses_table.xlsx"))
  })
  
  observeEvent(input$group_tab_open,{
    shell.exec(file.path(getwd(),"documents","multiple_grouped_analyses_table.xlsx"))
  })
  
  
  
  
  
  #─────────────────────────────────────────────────────────────────────────────────────────────────────────
  # INDIVIDUAL ANALYSIS ------------------------------------------------------------------------------------
  #─────────────────────────────────────────────────────────────────────────────────────────────────────────
  
  
  #─────────────────────────────────────────────────────────────────────────────────────────────────────────
  # Files input --------------------------------------------------------------------------------------------
  
  upload_DNA_seq<-reactive({
    inFile <- input$DNA_seq
    if (is.null(inFile))
      return(NULL)
    return(inFile$datapath)
  })
  
  
  # State values
  values <- reactiveValues(
    ab1_s1_state = NULL,
    ab1_s2_state = NULL
  )
  
  # If file uploaded
  observeEvent(input$ab1_s1, {
    values$ab1_s1_state <- 'uploaded'
  })
  observeEvent(input$ab1_s2, {
    values$ab1_s2_state <- 'uploaded'
  })
  
  # If file is reset
  observeEvent(input$reset_s1, {
    values$ab1_s1_state <- 'reset'
  })
  observeEvent(input$reset_s2, {
    values$ab1_s2_state <- 'reset'
  })
  
  # Sample table for individual analyses
  ab1_s1_input <- reactive({
    if (is.null(values$ab1_s1_state)) {
      return(NULL)
    } else if (values$ab1_s1_state == 'uploaded') {
      return(input$ab1_s1)
    } else if (values$ab1_s1_state == 'reset') {
      return(NULL)
    }
  })
  ab1_s1_file <- reactive({
    inFile <- ab1_s1_input()
    return(inFile$datapath)
  })
  
  # Grouped parameters for grouped analyses
  ab1_s2_input <- reactive({
    if (is.null(values$ab1_s2_state)) {
      return(NULL)
    } else if (values$ab1_s2_state == 'uploaded') {
      return(input$ab1_s2)
    } else if (values$ab1_s2_state == 'reset') {
      return(NULL)
    }
  })
  ab1_s2_file <- reactive({
    inFile <- ab1_s2_input()
    return(inFile$datapath)
  })
  
  # Display the uploaded file names
  output$summary_ab1_s1 <- renderText({
    return(paste("Uploaded file:", ab1_s1_input()$name))
  })
  output$summary_ab1_s2 <- renderText({
    return(paste("Uploaded file:", ab1_s2_input()$name))
  })
  
  
  
  #─────────────────────────────────────────────────────────────────────────────────────────────────────────
  # Add input to enter new folder name
  observeEvent(input$folderexisting, {
    
    if(input$folderexisting=="Create new folder") {
      output$foldernew <- renderUI({
        # input$foldernew
        textInput("foldernew", label = "Enter new folder name","folder name") })
    } else {
      output$foldernew <- NULL
    }
  })
  
  #─────────────────────────────────────────────────────────────────────────────────────────────────────────
  # Update seqname folder 'input$folderselect' depending on input$folderselect
  observeEvent(input$folderexisting, {
    
    if(input$folderexisting=="Create new folder") {
      updateSelectInput(session = getDefaultReactiveDomain(), "seqexisting",
                        choices = "Create new sequence folder",
                        selected = "Create new sequence folder")
      
    } else { 
      updateSelectInput(session = getDefaultReactiveDomain(), "seqexisting",
                        choices = c(basename(list.dirs(file.path(getwd(), "results", input$folderexisting), recursive = F)), "Create new sequence folder"))
    }
  })
  
  
  #─────────────────────────────────────────────────────────────────────────────────────────────────────────
  # Add input to enter new folder name
  observeEvent(input$seqexisting, {
    
    if(input$seqexisting=="Create new sequence folder") {
      output$seqnew <- renderUI({
        # input$seqnew
        textInput("seqnew", label = "Enter new sequence folder name","sequence name") })
    } else {
      output$seqnew <- NULL
    }
  })
  
  
  
  #─────────────────────────────────────────────────────────────────────────────────────────────────────────
  # Pre-install genome -------------------------------------------------------------------------------------
  
  observeEvent(input$install_genomeI, {
    
    if (!(input$genomeI %in% installed.packages()[,"Package"])) {
      
      withProgress(message = "Genome pre-installation in progress, please wait !", value=0, {
        shiny::setProgress(0.1)
        BiocManager::install(input$genomeI, ask = FALSE)
        shiny::setProgress(0.8)
        suppressMessages(lapply(input$genomeI, library, character.only = TRUE, quietly = TRUE))  
        on.exit(
          output$installed_genomeI <- renderUI({
            fluidRow(p("Genome successfully pre-installed !", style="font-size:15px ;")) 
          }) 
        )
      })
    } else {
      output$installed_genomeI <- renderUI({
        fluidRow(p("Genome already pre-installed !", style="font-size:15px ;")) 
      })
    }
  })
  
  
  #─────────────────────────────────────────────────────────────────────────────────────────────────────────
  # Render Rmarkdown report of analysis --------------------------------------------------------------------
  
  observeEvent(input$indiv, {
    
    output$indivReport <- renderUI(NULL)
    
    # Get folder name
    if(input$folderexisting == "Create new folder") { foldername <- input$foldernew } else { foldername <- input$folderexisting }
    
    # Get sequence folder name
    if(input$seqexisting == "Create new sequence folder") { seqname <- input$seqnew } else { seqname <- input$seqexisting }
    
    
    # Sample name based on input
    indiv_name <- paste(substr(foldername, 1, 6), # folder name 6 first characters
                        seqname,
                        ifelse(input$collection=="0" | input$collection=="", "", input$collection), # collection name if exists
                        input$group,
                        ifelse(input$replicate=="0" | input$replicate=="", "", input$replicate),
                        ifelse(input$clone=="0" | input$clone=="", "", input$clone),sep="-")
    indiv_name <- gsub("--","-",indiv_name)
    indiv_name <- gsub("-$","",indiv_name)
    
    # Name of report file [ABSP-Indiv]_[folder-seqname-collection-group-replicate-clone]_[YYYYMMDD-HHMM].html
    indiv_file <- paste0(paste("ABSP-Indiv", indiv_name, format(Sys.time(), "%Y%m%d-%H%M"), sep = "_"),".html")
    
    # Path of report file 
    if (!dir.exists(file.path(getwd(),"reports"))) {suppressWarnings(dir.create(file.path(getwd(),"reports")))}
    indiv_path <- file.path(getwd(),"reports",indiv_file)
    
    # Set up parameters to pass to Rmd document
    indiv_params <- list(foldername = foldername,
                         seqname = seqname,
                         collection = input$collection,
                         group = input$group,
                         replicate = input$replicate,
                         clone = input$clone,
                         date_s1 = input$date_s1,
                         date_s2 = input$date_s2,
                         genome = input$genomeI,
                         DNA_seq = upload_DNA_seq(),
                         ab1_s1 = ab1_s1_file(),
                         ab1_s2 = ab1_s2_file())
    
    indiv_RMD <- list.files(path=file.path(getwd(),"scripts"), pattern ="ABSP_individual_analysis.Rmd$", full.names = T)
    
    withProgress(message = "Analysis in progress, please wait!", {
      rmarkdown::render(indiv_RMD, encoding = 'UTF-8', output_file = indiv_path, params = indiv_params, envir = new.env(parent = globalenv()))
      on.exit( {
        output$indivReport <- renderUI({
          fluidRow(
            align="center",
            br(),
            p("Analysis done !", style="font-size:17px ;"),
            p("The report file ",em(basename(indiv_path)), " has been generated in the ",em("reports")," directory", style="font-size:15px ;")
            
          )
        })
      })
    })
  })
  #─────────────────────────────────────────────────────────────────────────────────────────────────────────
  
  
  
  
  
  #─────────────────────────────────────────────────────────────────────────────────────────────────────────
  # GROUPED ANALYSIS ---------------------------------------------------------------------------------------
  #─────────────────────────────────────────────────────────────────────────────────────────────────────────
  
  
  
  # Update seqname folder 'input$folderselect' depending on input$folderselect
  observeEvent(input$folderselect, {
    updateSelectInput(session = getDefaultReactiveDomain(), "seqselect",
                      choices = basename(list.dirs(file.path(getwd(), "results", input$folderselect), recursive = F)),
                      selected = basename(list.dirs(file.path(getwd(), "results", input$folderselect), recursive = F))[1])
    
  })
  
  
  # Update group order selection with groups
  observeEvent( c(input$seqselect,input$exptype), {
    
    output$groupedReport <- renderUI(NULL) # Clear the displayed text of the previous report
    
    data_files <- reactive( {
      if(input$exptype =="Cloning-BSP") {
        data_path <- file.path(getwd(),"results",input$folderselect,input$seqselect,"individual_results_cloning","data")
      }
      if(input$exptype =="Direct-BSP") {
        data_path <- file.path(getwd(),"results",input$folderselect,input$seqselect,"individual_results_direct","data")
      }
      
      data_files <- as.list(list.files(path = data_path, pattern = "*_methdata.csv$",
                                       all.files = T, full.names = T, recursive = T, ignore.case = F, include.dirs = T))
      return(data_files)
    })
    
    
    if(length(data_files())!=0) {
      data_all <- lapply(data_files(), read.table, header=TRUE, check.names = FALSE)
      all <- purrr::reduce(data_all, rbind)
      listgroup <- unique(all$group)  
      updateSelectInput(session = getDefaultReactiveDomain(), "group_order",choices = listgroup)
    } else {
      updateSelectInput(session = getDefaultReactiveDomain(), "group_order",choices = "No files found for the provided sequence and experiment type")
    }
    
  })
  
  #─────────────────────────────────────────────────────────────────────────────────────────────────────────
  # Render Rmarkdown report of analysis --------------------------------------------------------------------
  
  observeEvent(input$grouped, {
    
    # Experiment name based on input
    grouped_name <- paste(substr(input$folderselect, 1, 6), # folder name 6 first characters
                          input$seqselect,
                          ifelse(input$exptype=="Direct-BSP", "directBSP", "cloningBSP"),sep="-")
    
    # Name of report file [ABSP-Indiv]_[folder-seqname-collection-group-replicate-clone]_[YYYYMMDD-HHMM].html
    grouped_file <- paste0(paste("ABSP-Grouped", grouped_name, format(Sys.time(), "%Y%m%d-%H%M"), sep = "_"),".html")
    
    # Path of report file 
    
    if (!dir.exists(file.path(getwd(),"reports"))) {suppressWarnings(dir.create(file.path(getwd(),"reports")))}
    grouped_path <- file.path(getwd(),"reports",grouped_file)
    
    # Set up parameters to pass to Rmd document
    grouped_params <- list(foldername = input$folderselect,
                           seqname = input$seqselect,
                           genome = input$genomeG,
                           cloning = input$exptype,
                           coll_sep = input$coll_sep,
                           pos_labels = input$pos_labels,
                           group_order = paste(input$group_order,collapse = ","),
                           sample_order = input$sample_order)
    
    grouped_RMD <- list.files(path=file.path(getwd(),"scripts"), pattern ="ABSP_grouped_analysis.Rmd$", full.names = T)
    
    
    withProgress(message = "Analysis in progress, please wait!", {
      rmarkdown::render(grouped_RMD, encoding = 'UTF-8', output_file = grouped_path, params = grouped_params, envir = new.env(parent = globalenv()))
      on.exit( {
        output$groupedReport <- renderUI({
          fluidRow(
            align="center",
            br(),
            p("Analysis done !", style="font-size:17px ;"),
            p("The report file ",em(basename(grouped_path)), " has been generated in the ",em("reports")," directory", style="font-size:15px ;")
            
          )
        })
      })
    })
  })
  #─────────────────────────────────────────────────────────────────────────────────────────────────────────
  
  
  
  
  
  #─────────────────────────────────────────────────────────────────────────────────────────────────────────
  # MULTIPLE ANALYSES --------------------------------------------------------------------------------------
  #─────────────────────────────────────────────────────────────────────────────────────────────────────────
  
  
  
  #─────────────────────────────────────────────────────────────────────────────────────────────────────────
  # Files input --------------------------------------------------------------------------------------------
  
  # Add input to enter new folder name
  observeEvent(input$foldermainexisting, {
    
    if(input$foldermainexisting=="Create new folder") {
      output$foldermainnew <- renderUI({
        # input$foldermainnew
        textInput("foldermainnew", label = "Enter new folder name","folder name") })
    } else {
      output$foldermainnew <- NULL
    }
  })
  
  # State values
  values <- reactiveValues(
    exptable_state = NULL,
    groupedparams_state = NULL
  )
  
  # If file uploaded
  observeEvent(input$exptable, {
    values$exptable_state <- 'uploaded'
  })
  observeEvent(input$groupedparams, {
    values$groupedparams_state <- 'uploaded'
  })
  
  # If file is reset
  observeEvent(input$reset_exptable, {
    values$exptable_state <- 'reset'
  })
  observeEvent(input$reset_groupedparams, {
    values$groupedparams_state <- 'reset'
  })
  
  # Sample table for individual analyses
  exptable_input <- reactive({
    if (is.null(values$exptable_state)) {
      return(NULL)
    } else if (values$exptable_state == 'uploaded') {
      return(input$exptable)
    } else if (values$exptable_state == 'reset') {
      return(NULL)
    }
  })
  exptable_file <- reactive({
    inFile <- exptable_input()
    return(inFile$datapath)
  })
  
  # Grouped parameters for grouped analyses
  groupedparams_input <- reactive({
    if (is.null(values$groupedparams_state)) {
      return(NULL)
    } else if (values$groupedparams_state == 'uploaded') {
      return(input$groupedparams)
    } else if (values$groupedparams_state == 'reset') {
      return(NULL)
    }
  })
  groupedparams_file <- reactive({
    inFile <- groupedparams_input()
    return(inFile$datapath)
  })
  
  # Display the uploaded file names
  output$summary_exptable <- renderText({
    return(paste("Uploaded file:", exptable_input()$name))
  })
  output$summary_groupedparams <- renderText({
    return(paste("Uploaded file:", groupedparams_input()$name))
  })
  
  
  
  #─────────────────────────────────────────────────────────────────────────────────────────────────────────
  # Action button obs event --------------------------------------------------------------------------------
  
  observeEvent(input$multiple, {
    
    # Get folder name
    if(input$foldermainexisting == "Create new folder") { foldermain <- input$foldermainnew } else { foldermain <- input$foldermainexisting }
    
    
    #─────────────────────────────────────────────────────────────────────────────────────────────────────────
    # Create directories
    
    output_res <- file.path(getwd(),"results")
    if (!dir.exists(output_res)) {suppressWarnings(dir.create(output_res))}
    
    output_rep <- file.path(getwd(),"reports")
    if (!dir.exists(output_rep)) {suppressWarnings(dir.create(output_rep))}
    
    dir_main <- file.path(getwd(),"results", foldermain)
    if (!dir.exists(dir_main)) {suppressWarnings(dir.create(dir_main))}
    
    
    
    #─────────────────────────────────────────────────────────────────────────────────────────────────────────
    # All individual analysis --------------------------------------------------------------------------------
    #─────────────────────────────────────────────────────────────────────────────────────────────────────────
    
    if (!is.null(exptable_file())) {
      # Read file
      if (grepl(".xlsx$", exptable_file())==T) { exp_table <- read.xlsx(exptable_file(), cols= c(1:11), detectDates = T) }
      if (grepl(".csv$", exptable_file())==T) { exp_table <- read.csv2(exptable_file(), na.strings="") }
      
      
      exp_table <- exp_table[,1:11] # keep only 11 first columns
      exp_table <- exp_table[which(complete.cases(exp_table[,1])),] # remove rows without sequence name
      
      # Remove extra character if needed
      exp_table[,7] <- file.path(gsub("\"", "",exp_table[,7]))
      exp_table[,9] <- file.path(gsub("\"", "",exp_table[,9]))
      exp_table[,11] <- file.path(gsub("\"", "",exp_table[,11]))
      
      # Remove additionnal rows
      exp_table <- exp_table[!exp_table[,1]=="",]
      
      # Date format
      if (class(exp_table[,8])!="Date") {exp_table[,8] <- as.Date(exp_table[,8], tryFormats = c("%Y-%m-%d","%m/%d/%y"))}
      if (class(exp_table[,10])!="Date") {exp_table[,10] <- as.Date(exp_table[,10], tryFormats = c("%Y-%m-%d","%m/%d/%y"))}
      exp_table[,8] <- as.character(exp_table[,8])
      exp_table[,10] <- as.character(exp_table[,10])
      
      
      #─────────────────────────────────────────────────────────────────────────────────────────────────────────
      # Pre-install genomes
      
      for (i in unique(exp_table[,6])) {
        if (!(i %in% installed.packages()[,"Package"])) {
          BiocManager::install(i, ask = FALSE)
        } 
        suppressMessages(lapply(i, library, character.only = TRUE, quietly = TRUE))
      }
      
      
      #─────────────────────────────────────────────────────────────────────────────────────────────────────────
      # Run all individual analysis from the sample data table -------------------------------------------------
      
      list_Ireports <- c()
      
      for (i in 1:nrow(exp_table)) {
        
        
        # Sample name based on idata table
        indiv_name <- paste(substr(foldermain, 1, 6), # folder name 6 first characters
                            exp_table[i,1],
                            ifelse(exp_table[i,2]=="0" | is.na(exp_table[i,2]), "", exp_table[i,2]), # collection name if exists
                            exp_table[i,3],
                            ifelse(exp_table[i,4]=="0" | is.na(exp_table[i,4]), "", exp_table[i,4]),
                            ifelse(exp_table[i,5]=="0" | is.na(exp_table[i,5]), "", exp_table[i,5]), sep="-")
        indiv_name <- gsub("--","-",indiv_name)
        indiv_name <- gsub("-$","",indiv_name)
        
        # Name of report file [ABSP-Indiv]_[folder-seqname-collection-group-replicate-clone]_[YYYYMMDD-HHMM].html
        indiv_file <- paste0(paste("ABSP-Indiv", indiv_name, format(Sys.time(), "%Y%m%d-%H%M"), sep = "_"),".html")
        
        # Path of report file
        indiv_path <- file.path(output_rep,indiv_file)
        
        # Set up parameters to pass to Rmd document
        indiv_params <- list(foldername = foldermain,
                             seqname = exp_table[i,1],
                             collection = ifelse(!is.na(exp_table[i,2]),exp_table[i,2],""),
                             group = exp_table[i,3],
                             replicate = as.integer(ifelse(!is.na(exp_table[i,4]),exp_table[i,4],"")),
                             clone = as.integer(ifelse(!is.na(exp_table[i,5]),exp_table[i,5],"")),
                             genome = exp_table[i,6],
                             DNA_seq = exp_table[i,7],
                             date_s1 =  exp_table[i,8],
                             ab1_s1 = exp_table[i,9],
                             date_s2 = exp_table[i,10],
                             ab1_s2 = exp_table[i,11])
        
        indiv_RMD <- list.files(path=file.path(getwd(),"scripts"), pattern ="ABSP_individual_analysis.Rmd$", full.names = T)
        
        withProgress(message = "Analysis in progress, please wait!", {
          rmarkdown::render(indiv_RMD, encoding = 'UTF-8', output_file = indiv_path, params = indiv_params)
          on.exit( {
            output$multipleIReport <- renderUI({
              fluidRow(
                align="center",
                br(),
                p(basename(indiv_path)," has been generated in the ", em("reports")," directory", style="font-size:15px ;")
              )
            })
          })
        })
      }
    }
    
    #─────────────────────────────────────────────────────────────────────────────────────────────────────────
    # All grouped analysis -----------------------------------------------------------------------------------
    #─────────────────────────────────────────────────────────────────────────────────────────────────────────
    
    if (!is.null(groupedparams_file())) {
      
      # Read file
      if (grepl(".xlsx$", groupedparams_file())==T) { grouped_table <- read.xlsx(groupedparams_file()) }
      if (grepl(".csv$", groupedparams_file())==T) { grouped_table <- read.csv2(groupedparams_file(), na.strings="") }
      grouped_table <- grouped_table[,1:7] # 7 first columns, remove options table
      grouped_table <- grouped_table[which(complete.cases(grouped_table)),] # remove NA rows
      grouped_table <- grouped_table[which(complete.cases(grouped_table[,1])),] # remove rows without sequence name
      
      
      
      #─────────────────────────────────────────────────────────────────────────────────────────────────────────
      # Pre-install genomes
      
      for (i in unique(grouped_table[,2])) {
        if (!(i %in% installed.packages()[,"Package"])) {
          BiocManager::install(i, ask = FALSE)
        } 
        suppressMessages(lapply(i, library, character.only = TRUE, quietly = TRUE))
      }
      
      
      
      #─────────────────────────────────────────────────────────────────────────────────────────────────────────
      # Run all grouped analysis of previously analysed experiments --------------------------------------------
      
      list_Greports <- c()
      
      for (i in 1:nrow(grouped_table)) {
        
        # Experiment name based on input
        grouped_name <- paste(substr(foldermain, 1, 6), # folder name 6 first characters
                              grouped_table[i,1],
                              ifelse(grouped_table[i,3]=="Direct-BSP", "directBSP", "cloningBSP"),sep="-")
        
        # Name of report file [ABSP-Indiv]_[folder-seqname-collection-group-replicate-clone]_[YYYYMMDD-HHMM].html
        grouped_file <- paste0(paste("ABSP-Grouped", grouped_name, format(Sys.time(), "%Y%m%d-%H%M"), sep = "_"),".html")
        
        # Path of report file
        grouped_path <- file.path(getwd(),"reports",grouped_file)
        
        # Split string by commas to get list of sample orders
        sample_order <- strsplit(grouped_table[i,7],",")[[1]] 
        sample_order <- gsub("^ ","", sample_order) # remove first spaces
        sample_order <- gsub(" $","", sample_order) # remove last spaces
        
        # Set up parameters to pass to Rmd document
        grouped_params <- list(foldername = foldermain,
                               seqname = grouped_table[i,1],
                               genome = grouped_table[i,2],
                               cloning = grouped_table[i,3],
                               pos_labels = grouped_table[i,4],
                               coll_sep = grouped_table[i,5],
                               group_order = grouped_table[i,6],
                               sample_order = sample_order)
        
        grouped_RMD <- list.files(path=file.path(getwd(),"scripts"), pattern ="ABSP_grouped_analysis.Rmd$", full.names = T)
        
        withProgress(message = "Analysis in progress, please wait!", {
          rmarkdown::render(grouped_RMD, encoding = 'UTF-8', output_file = grouped_path, params = grouped_params, envir = new.env(parent = globalenv()))
          on.exit( {
            output$multipleIReport <- renderUI(NULL)
            output$multipleGReport <- renderUI({
              fluidRow(
                align="center",
                br(),
                p(basename(grouped_path)," has been generated in the ", em("reports")," directory", style="font-size:15px ;")
              )
            })
          })
        })
        
      }
    }
    
  })
  
  
}





# Run the application 
shinyApp(ui = ui, server = server)


# > sessionInfo()
# R version 4.2.1 (2022-06-23 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=French_France.utf8  LC_CTYPE=French_France.utf8    LC_MONETARY=French_France.utf8
# [4] LC_NUMERIC=C                   LC_TIME=French_France.utf8    
# 
# attached base packages:
#   [1] tools     parallel  grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] sangeranalyseR_1.6.1  logger_0.2.2          BiocStyle_2.24.0      shinyWidgets_0.7.3    shinycssloaders_1.0.0
# [6] excelR_0.4.0          zeallot_0.1.0         DT_0.25               data.table_1.14.2     shinyjs_2.1.0        
# [11] shinydashboard_0.7.2  gridExtra_2.3         sangerseqR_1.32.0     phangorn_2.10.0       reshape2_1.4.4       
# [16] DECIPHER_2.24.0       RSQLite_2.2.17        ape_5.6-2             stringr_1.4.1         Gviz_1.40.1          
# [21] BSgenome_1.64.0       rtracklayer_1.56.1    GenomicRanges_1.48.0  Biostrings_2.64.0     XVector_0.36.0       
# [26] webshot_0.5.3         shinythemes_1.2.0     shiny_1.7.2           seqinr_4.2-16         rstatix_0.7.0        
# [31] Rmisc_1.5.1           plyr_1.8.7            lattice_0.20-45       rmarkdown_2.16        rlist_0.4.6.2        
# [36] readr_2.1.2           RColorBrewer_1.1-3    purrr_0.3.4           png_0.1-7             plotly_4.10.0        
# [41] pdftools_3.3.0        openxlsx_4.2.5        knitr_1.40            htmlwidgets_1.5.4     htmltools_0.5.3      
# [46] ggpubr_0.4.0          ggplot2_3.3.6         ggdendro_0.1.23       GenomeInfoDb_1.32.2   IRanges_2.30.0       
# [51] S4Vectors_0.34.0      BiocGenerics_0.42.0   formattable_0.2.1     dplyr_1.0.10          DiagrammeR_1.0.9     
# [56] compareGroups_4.5.1   BiocManager_1.30.18   arrangements_1.1.9   
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.2.2                  tidyselect_1.1.2            AnnotationDbi_1.58.0       
# [4] BiocParallel_1.30.3         gmp_0.6-6                   munsell_0.5.0              
# [7] codetools_0.2-18            interp_1.1-3                chron_2.3-57               
# [10] withr_2.5.0                 colorspace_2.0-3            Biobase_2.56.0             
# [13] filelock_1.0.2              uuid_1.1-0                  rstudioapi_0.14            
# [16] ggsignif_0.6.3              officer_0.4.4               MatrixGenerics_1.8.1       
# [19] GenomeInfoDbData_1.2.8      bit64_4.0.5                 vctrs_0.4.1                
# [22] generics_0.1.3              xfun_0.33                   biovizBase_1.44.0          
# [25] BiocFileCache_2.4.0         qpdf_1.2.0                  R6_2.5.1                   
# [28] AnnotationFilter_1.20.0     bitops_1.0-7                cachem_1.0.6               
# [31] DelayedArray_0.22.0         assertthat_0.2.1            promises_1.2.0.1           
# [34] BiocIO_1.6.0                scales_1.2.1                nnet_7.3-17                
# [37] gtable_0.3.1                ensembldb_2.20.2            rlang_1.0.5                
# [40] systemfonts_1.0.4           splines_4.2.1               lazyeval_0.2.2             
# [43] dichromat_2.0-0.1           checkmate_2.1.0             broom_1.0.1                
# [46] yaml_2.3.5                  abind_1.4-5                 GenomicFeatures_1.48.3     
# [49] backports_1.4.1             httpuv_1.6.6                HardyWeinberg_1.7.5        
# [52] Hmisc_4.7-1                 ellipsis_0.3.2              kableExtra_1.3.4           
# [55] jquerylib_0.1.4             Rsolnp_1.16                 Rcpp_1.0.9                 
# [58] base64enc_0.1-3             visNetwork_2.1.0            progress_1.2.2             
# [61] zlibbioc_1.42.0             RCurl_1.98-1.8              prettyunits_1.1.1          
# [64] rpart_4.1.16                deldir_1.0-6                SummarizedExperiment_1.26.1
# [67] cluster_2.1.4               magrittr_2.0.3              flextable_0.8.1            
# [70] truncnorm_1.0-8             ProtGenerics_1.28.0         matrixStats_0.62.0         
# [73] hms_1.1.2                   mime_0.12                   evaluate_0.16              
# [76] xtable_1.8-4                XML_3.99-0.10               jpeg_0.1-9                 
# [79] compiler_4.2.1              biomaRt_2.52.0              tibble_3.1.8               
# [82] mice_3.14.0                 writexl_1.4.0               crayon_1.5.1               
# [85] later_1.3.0                 tzdb_0.3.0                  Formula_1.2-4              
# [88] tidyr_1.2.1                 DBI_1.1.3                   dbplyr_2.2.1               
# [91] MASS_7.3-58.1               rappdirs_0.3.3              Matrix_1.5-1               
# [94] ade4_1.7-19                 car_3.1-0                   cli_3.3.0                  
# [97] quadprog_1.5-8              igraph_1.3.4                pkgconfig_2.0.3            
# [100] GenomicAlignments_1.32.0    foreign_0.8-82              xml2_1.3.3                 
# [103] svglite_2.1.0               bslib_0.4.0                 rvest_1.0.3                
# [106] VariantAnnotation_1.42.1    digest_0.6.29               fastmatch_1.1-3            
# [109] htmlTable_2.4.1             gdtools_0.2.4               restfulr_0.0.15            
# [112] curl_4.3.2                  Rsamtools_2.12.0            rjson_0.2.21               
# [115] nlme_3.1-159                lifecycle_1.0.2             jsonlite_1.8.0             
# [118] carData_3.0-5               viridisLite_0.4.1           askpass_1.1                
# [121] fansi_1.0.3                 pillar_1.8.1                KEGGREST_1.36.2            
# [124] fastmap_1.1.0               httr_1.4.4                  survival_3.4-0             
# [127] glue_1.6.2                  zip_2.2.1                   bit_4.0.4                  
# [130] sass_0.4.2                  stringi_1.7.8               blob_1.2.3                 
# [133] latticeExtra_0.6-30         memoise_2.0.1               renv_0.15.5 

