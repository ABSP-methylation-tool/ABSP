source("renv/activate.R")


options(timeout = 600)


renv::restore(prompt=F)


lib <- c("arrangements","BiocGenerics","BiocManager","Biostrings","BSgenome","compareGroups",
         "DiagrammeR","dplyr","formattable", "GenomeInfoDb","GenomicRanges","ggdendro",
         "ggplot2", "ggpubr","Gviz","htmltools","htmlwidgets",
         "knitr","openxlsx","pdftools","plotly",
         "png","purrr","RColorBrewer", "readr","renv","rlist","rmarkdown",
         "Rmisc","rstatix","sangeranalyseR","sangerseqR","seqinr", "shiny",
         "shinythemes","webshot")
suppressMessages(lapply(lib, library, character.only = TRUE, quietly = TRUE))

source("./scripts/ABSP_functions.R")
