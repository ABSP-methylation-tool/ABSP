#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# ABSP: Analysis of Bisulfite Sequencing PCR
#
# Download ABSP at https://github.com/ABSP-methylation-tool/ABSP
# Find more information in the User Guide



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
options(timeout = 3000)
packages <- c("arrangements","BiocManager","compareGroups","DiagrammeR","dplyr","formattable","GenomeInfoDb",
              "ggdendro","ggplot2","ggpubr","htmltools","htmlwidgets","knitr","openxlsx","pdftools","plotly",
              "png","purrr","RColorBrewer","readr","rlist","rmarkdown","Rmisc","rstatix","seqinr","shiny", "shinybusy",
              "shinythemes","webshot")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(new.packages)>0) {
  install.packages(new.packages)}

# Install packages from Bioconductor
Biocpackages <-c("BiocGenerics","Biostrings","BSgenome","GenomicRanges","Gviz","sangeranalyseR","sangerseqR")
new.Biocpackages <- Biocpackages[!(Biocpackages %in% installed.packages()[,"Package"])]
if (length(new.Biocpackages)>0) {
  BiocManager::install(new.Biocpackages)}

# Add packages to library
suppressWarnings(lapply(packages, library, character.only = TRUE)) #, quietly = TRUE
suppressWarnings(lapply(Biocpackages, library, character.only = TRUE))

source(file = "./scripts/ABSP_functions.R")

#────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────





ui <- fluidPage(
    
    # Theme
    theme = shinytheme("flatly"),
    includeCSS("www/custom_app.css"),
    
    add_busy_spinner(spin = "fading-circle",color="#40b7a0",position="bottom-right", margins = c(60, 10)),
    
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
                                   
                                   p("ABSP is a R-based tool designed to analyze results from Bisulfite Sequencing PCR (BSP) experiments, hence its acronym 'Analysis of Bisulfite Sequencing PCR'."), 
                                   p("It was developed to assist researchers in estimating and comparing CpG methylation percentages of DNA regions studied using BSP experiments."),
                                   p("It offers a comprehensive automated workflow, spanning from trace file sequencing results to data visualization and statistical analysis."),
                                   
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
                                 
                                 p(icon(name="newspaper", lib = "font-awesome", style="padding-right:8px; font-size:1.2em;"),"Marie Denoulet, Mathilde Brulé, François Anquez, Audrey Vincent, Julie Schnipper, Eric Adriaenssens, Robert-Alain Toillon, Xuefen Le Bourhis, Chann Lagadec, ",
                                   span("ABSP: an automated R tool to efficiently analyze region-specific CpG methylation from bisulfite sequencing PCR", style= "font-weight: bold"),", ",
                                   em("Bioinformatics"),", Volume 39, Issue 1, January 2023, btad008, ", htmltools::tags$a(href="https://doi.org/10.1093/bioinformatics/btad008", "https://doi.org/10.1093/bioinformatics/btad008",target="_blank")),
                                 
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
                                    
                                    h5("Be aware that these entries must not contain any special characters such as :", 
                                       htmltools::tags$code("/ \\ : * ? ! \" ' ` < > | & % @ # + = { }", style="color:grey ; background-color:transparent"),
                                       style="color:#bf3232"),
                                    
                                    h5("Allowed characters (that must be avoided if not necessary) :", 
                                       htmltools::tags$code(" space - _ ", style="color:grey ; background-color:transparent"), 
                                       style="color:#bf3232"),
                                    
                                    h5("All entries must be named consistently across the different analyses 
                                                    (e.g. the group names must strickly be identical between samples belonging to the same group).", 
                                       style="color:#bf3232"),
                                    
                                    hr(),
                                    
                                    div( # FOLDER
                                        h4("Select existing folder", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("Select an existing folder within the ABSP results folder to locate all of the analysis results.
                                                       Organizing results into separate folders can help differentiate analyses by projects, 
                                                       experiments or users."),
                                        p("To create a new folder, select the 'Create new folder' option and enter the desired 
                                                        name in the text input. Please note that the six first letters of the folder name will appear in the report file name.")
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
                                        p("The term 'collection' refers to a segragation of samples above the level of 'groups'. 
                                          Samples from different collections cannot be compared, even if they belong to the same group. 
                                          For instance, collections might represent different cell lines, organs, or patients, where comparisons 
                                          between groups are made within each collection but not between different collections.  
                                          To actually compare these types of samples, consider them as groups instead collections."),
                                        em("If you do not want to specify any collection, leave empty or enter '0'. 
                                           Make sure the collection name is strictly identical for all samples of the same collection.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # GROUP
                                        h4("Enter group name", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("The term 'group' represents the condition you intend to compare in the grouped analysis. 
                                          For example, groups can be the 'control' and 'treated' conditions."),
                                        em("The group field entry is required. Ensure that the group name is consistent across all samples belonging to the same group.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # REPLICATE
                                        h4("Direct-BSP - Enter replicate number", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("In the case of direct sequencing of PCR products only (direct-BSP). 
                                          The 'replicate number' refers to the number of repetitions within a single group, for experiment reproducibility 
                                          and statistic significance determination. It is generally recommended to have a minimum of 3 replicates for statistical analysis.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # CLONE
                                        h4("Cloning-BSP - Enter clone number", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("In the case of clone sequencing only (cloning-BSP). 
                                          The 'clone number' represents the identification number assigned to each clone within a group. 
                                          Typically, sequencing around 10 clones is advised to ensure the reliability of results.")
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
                                        p("Select the reference genome. This is solely used for displaying the genomic sequence in the genomic plot."),
                                        p("Ensure to click on the 'Pre-install genome' button if the selected genome is being used for the first time.", 
                                          style="color:#bf3232"),
                                        em("Note that only a limited number of genomes are listed in the selection tab, but more genomes can be used. 
                                           If your studied genome is not listed please refer to the user guide for further instructions.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # FASTA
                                        h4("Select .fasta file of reference DNA sequence", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("Select the '.fasta' file from your folders to provide the reference DNA sequence from the plus strand of genome.
                                                       The fasta file header must contain two elements :"),
                                        htmltools::tags$ul(
                                          htmltools::tags$li("The precise genomic coordinates of the sequence contained in the file, following this strict format : ", br(),
                                                    span("chr#:######-###### (e.g. chr16:68771087-68771462)", style="color:#40b7a0")),
                                            
                                          htmltools::tags$li("The strand chosen for primer design : the strand complementary to primers after bisulfite 
                                                                 conversion, following this strict format : ", br(),
                                                    span("primers=plus or primers=minus", style="color:#40b7a0"))
                                        ),
                                        em("Note that any additional information in the header, such as the sequence name for example,
                                                        can be included without consequences as long as it does not interfere with the formats described above."),
                                        
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
                                        h4("Select sequencing dates", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("For traceability purposes, you can select the dates when the sequencing runs were performed."),
                                        em("If one of the sequencing files is not provided, leave the date entry empty.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # SEQ 
                                        h4("Select .ab1 sequencing files", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("Select the files in your folders corresponding to the .ab1 sequencing files. 
                                          Both the sequencing reads from one direction (#1) and the other direction (#2) should be provided, although the analysis can be run with only one sequencing read provided."),
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
                            
                            # input$group_order
                            selectInput("group_order",label ="Select groups in the desired order for display",
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
                                        p("Select an existing folder within the ABSP results folder to retrieve previously 
                                                       generated methylation data files after individual analyses, and to locate the newly 
                                                       generated results from the grouped analysis.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # SEQUENCE
                                        h4("Select sequence", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("Select an existing folder corresponding to the name of your sequence, to retrieve previously 
                                                       generated methylation data files after individual analyses, and to locate the newly 
                                                       generated results of the grouped analysis.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # GENOME
                                        h4("Select genome", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("Select the reference genome. This is solely used to display the genomic sequence in the genomic plots.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # EXP
                                        h4("Select experiment type", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("Select the experiment type as either 'Direct-BSP' or 'Cloning-BSP', to ensure proper 
                                                       retrieval of methylation data and specification of analysis.") 
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
                                        p("Select the label type to display as CpG positions in plots, options include 'CpG coordinates', 'CpG numbers' and 'None'.
                                               Note that in cases of extremely close CpG positions, labels may overlap, so selecting 'None' can be a suitable alternative.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # COLL SEP
                                        h4("Collection separation", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("This parameter determines whether collections are completely separated for display in plots or not. In both cases, groups of different collection are not compared.
                                          If unchecked, plots will represents all samples of all collections together. If checked each collection will be represented in a separated plot.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # GROUP ORDER
                                        h4("Select groups in the desired order for display", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("The names of groups are retrieved from files in the selected folder above. 
                                          Select the desired groups for analysis in the correct order to display samples in plots.")
                                    ),
                                    
                                    hr(),
                                    
                                    div( # SAMPLE ORDER
                                        h4("Select the types of sample ordering for plots", style="color:#374e64 ; font-weight: bold ; font-size:17px"),
                                        p("Choose from one to four available options for sample ordering in plots :"),
                                        htmltools::tags$ul(
                                          htmltools::tags$li("'As it is' arranges samples by alphabetic order of collections. 
                                                    If none or only one collection is present, this order is equivalent to 'By groups'."),
                                          htmltools::tags$li("'By groups' arranges samples by the provided group order above."),
                                          htmltools::tags$li("'By methylation levels' arranges samples based on their mean methylation levels."),
                                          htmltools::tags$li("'By clusters' arranges samples based on hierarchical clustering, which is depicted in an associated dendrogram.")
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
                                          p("Several analyses can be launched from this tab with just one click, using pre-filled tables containing the input entries."),
                                          
                                          br(),
                                          
                                          h4("How to proceed ?"),
                                          p("1. Fill out one or the input tables below with the sample information and the chosen parameters 
                                            (you can open the tables files with the buttons provided below)."),
                                          
                                          p("2. In the left panel, select an existing folder within the ABSP results folder to store all of the analyses results."),
                                          p("To create a new folder, select the 'Create new folder' option and enter desired name in the text input. Please note that the six first letters of the folder name will appear in the report file name."),
                                          
                                          
                                          p("3. In the left panel, select your filled table as input."), 
                                          p("Both the sample data table for multiple individual analyses and the parameter table for multiple grouped analyses can be provided concomitantly to launch individual analyses followed by grouped analyses.
                                             Alternatively, you can provide only one of the two tables, which will launch the corresponding analyses, either individual or grouped."),
                                          
                                          p("4. Launch the analyses by clicking on the 'Run analyses' button. "),
                                          p("All the input analyses will be launch one after the other, do not close the app until all analyses have completed."),
                                        
                                          em("Note: The tables files are located in the 'documents' folder. The input file format must be either '.xlsx' or '.csv'."),
                                          
                                          hr(),
                                          
                                          h4("Sample data for individual analyses"),
                                          p("The table compiles the inputs for individual analyses."),
                                          
                                          actionButton("exp_tab_open"," Open multiple individual analysis table file", 
                                                       icon=shiny::icon("folder-open", lib = "font-awesome", style="padding-right:8px ;"),
                                                       style="font-size:17px; margin:10px;"),
                                          
                                          hr(),
                                          
                                          h4("Parameters for the grouped analyses"),
                                          p("The table compiles the inputs for grouped analyses."),
                                          
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
    
    div(class="navbar navbar-default navbar-fixed-bottom",
        style="min-height: 40px;",
        fluidRow(
            div(style="padding-right: 15px; padding-left: 15px;", 
                
                column(width=6,
                       p("ABSP v1.2.2 - Copyright © 2023 CANTHER laboratory, released under the GPL-3 license", style="font-size:14px ; color:#ffffff ; padding-top:15px"))
                
            )
            
        )
        
        
    )
    
    
    
)



##############################################################################################################



server <- function(input, output, session) {
  
  
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
# R version 4.3.3 (2024-02-29 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 11 x64 (build 22631)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=French_France.utf8  LC_CTYPE=French_France.utf8    LC_MONETARY=French_France.utf8 LC_NUMERIC=C                  
# [5] LC_TIME=French_France.utf8    
# 
# time zone: Europe/Paris
# tzcode source: internal
# 
# attached base packages:
#   [1] tools     parallel  grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] sangeranalyseR_1.12.0 logger_0.3.0          BiocStyle_2.30.0      shinyWidgets_0.8.3    shinycssloaders_1.0.0
# [6] excelR_0.4.0          zeallot_0.1.0         DT_0.32               data.table_1.15.4     shinyjs_2.1.0        
# [11] shinydashboard_0.7.2  gridExtra_2.3         sangerseqR_1.38.0     phangorn_2.11.1       reshape2_1.4.4       
# [16] DECIPHER_2.30.0       RSQLite_2.3.6         ape_5.7-1             stringr_1.5.1         Gviz_1.46.1          
# [21] BSgenome_1.70.2       rtracklayer_1.62.0    BiocIO_1.12.0         GenomicRanges_1.54.1  Biostrings_2.70.3    
# [26] XVector_0.42.0        webshot_0.5.5         shinythemes_1.2.0     shinybusy_0.3.3       seqinr_4.2-36        
# [31] rstatix_0.7.2         Rmisc_1.5.1           plyr_1.8.9            lattice_0.22-5        rmarkdown_2.26       
# [36] rlist_0.4.6.2         readr_2.1.5           RColorBrewer_1.1-3    purrr_1.0.2           png_0.1-8            
# [41] plotly_4.10.4         pdftools_3.4.0        openxlsx_4.2.5.2      knitr_1.45            htmlwidgets_1.6.4    
# [46] htmltools_0.5.8       ggpubr_0.6.0          ggplot2_3.5.0         ggdendro_0.2.0        GenomeInfoDb_1.38.8  
# [51] IRanges_2.36.0        S4Vectors_0.40.2      BiocGenerics_0.48.1   formattable_0.2.1     dplyr_1.1.4          
# [56] DiagrammeR_1.0.11     compareGroups_4.8.0   BiocManager_1.30.22   arrangements_1.1.9    shiny_1.8.1.1        
# 
# loaded via a namespace (and not attached):
#   [1] ProtGenerics_1.34.0         matrixStats_1.2.0           bitops_1.0-7                fontawesome_0.5.2          
# [5] httr_1.4.7                  backports_1.4.1             utf8_1.2.4                  R6_2.5.1                   
# [9] lazyeval_0.2.2              jomo_2.7-6                  withr_3.0.0                 prettyunits_1.2.0          
# [13] cli_3.6.2                   Biobase_2.62.0              textshaping_0.3.7           officer_0.6.5              
# [17] sass_0.4.9                  askpass_1.2.0               Rsamtools_2.18.0            systemfonts_1.0.6          
# [21] foreign_0.8-86              gfonts_0.2.0                svglite_2.1.3               dichromat_2.0-0.1          
# [25] rstudioapi_0.16.0           httpcode_0.3.0              visNetwork_2.1.2            generics_0.1.3             
# [29] shape_1.4.6.1               car_3.1-2                   zip_2.3.1                   Matrix_1.6-5               
# [33] interp_1.1-6                fansi_1.0.6                 abind_1.4-5                 lifecycle_1.0.4            
# [37] yaml_2.3.8                  carData_3.0-5               SummarizedExperiment_1.32.0 SparseArray_1.2.4          
# [41] BiocFileCache_2.10.2        blob_1.2.4                  promises_1.2.1              crayon_1.5.2               
# [45] mitml_0.4-5                 GenomicFeatures_1.54.4      KEGGREST_1.42.0             pillar_1.9.0               
# [49] rjson_0.2.21                boot_1.3-29                 codetools_0.2-19            fastmatch_1.1-4            
# [53] pan_1.9                     glue_1.7.0                  fontLiberation_0.1.0        qpdf_1.3.3                 
# [57] vctrs_0.6.5                 gtable_0.3.4                cachem_1.0.8                xfun_0.43                  
# [61] S4Arrays_1.2.1              mime_0.12                   survival_3.5-8              iterators_1.0.14           
# [65] gmp_0.7-4                   nlme_3.1-164                bit64_4.0.5                 fontquiver_0.2.1           
# [69] progress_1.2.3              filelock_1.0.3              bslib_0.7.0                 rpart_4.1.23               
# [73] colorspace_2.1-0            DBI_1.2.2                   Hmisc_5.1-2                 nnet_7.3-19                
# [77] ade4_1.7-22                 tidyselect_1.2.1            bit_4.0.5                   compiler_4.3.3             
# [81] curl_5.2.1                  chron_2.3-61                glmnet_4.1-8                htmlTable_2.4.2            
# [85] HardyWeinberg_1.7.7         flextable_0.9.5             mice_3.16.0                 xml2_1.3.6                 
# [89] fontBitstreamVera_0.1.1     DelayedArray_0.28.0         checkmate_2.3.1             scales_1.3.0               
# [93] quadprog_1.5-8              rappdirs_0.3.3              digest_0.6.35               minqa_1.2.6                
# [97] pkgconfig_2.0.3             jpeg_0.1-10                 base64enc_0.1-3             lme4_1.1-35.2              
# [101] MatrixGenerics_1.14.0       dbplyr_2.5.0                fastmap_1.1.1               ensembldb_2.26.0           
# [105] rlang_1.1.3                 jquerylib_0.1.4             jsonlite_1.8.8              BiocParallel_1.36.0        
# [109] VariantAnnotation_1.48.1    RCurl_1.98-1.14             magrittr_2.0.3              kableExtra_1.4.0           
# [113] Formula_1.2-5               GenomeInfoDbData_1.2.11     munsell_0.5.1               Rcpp_1.0.12                
# [117] gdtools_0.3.7               stringi_1.8.3               zlibbioc_1.48.2             MASS_7.3-60.0.1            
# [121] deldir_2.0-4                splines_4.3.3               hms_1.1.3                   igraph_2.0.3               
# [125] uuid_1.2-0                  ggsignif_0.6.4              biomaRt_2.58.2              crul_1.4.0                 
# [129] XML_3.99-0.16.1             evaluate_0.23               latticeExtra_0.6-30         biovizBase_1.50.0          
# [133] nloptr_2.0.3                tzdb_0.4.0                  foreach_1.5.2               httpuv_1.6.15              
# [137] tidyr_1.3.1                 openssl_2.1.1               broom_1.0.5                 xtable_1.8-4               
# [141] restfulr_0.0.15             AnnotationFilter_1.26.0     Rsolnp_1.16                 later_1.3.2                
# [145] viridisLite_0.4.2           ragg_1.3.0                  truncnorm_1.0-9             tibble_3.2.1               
# [149] memoise_2.0.1               AnnotationDbi_1.64.1        GenomicAlignments_1.38.2    writexl_1.5.0              
# [153] cluster_2.1.6 