# How to use the examples of data ?

Two data sets are available:
* Data from a real sequence named **"CDH1"** corresponding to real direct-BSP data.  
Low and high methylated DNA were bisulfite converted, the CDH1 sequence amplified by PCR and PCR products were directly sequenced three times, using forward and reverse primers.
* Data from a mock sequence named **"Test sequence"** generated *in silico*.  
The methylation data tables of several samples were generated to reproduce both direct-BSP and cloning-BSP data: from two collections "collection1" and collection2", two groups "group1" and "group2" and 3 replicates for direct-BSP or 10 clones for cloning-BSP.


## "CDH1" data set

### Individual analysis 
Find in the "input" folder: 
* The .ab1 sequencing files of CDH1  
    + 2 groups: "High methylated DNA" and "Low methylated DNA"  
    + 3 replicates: "rep1", "rep2" and "rep3"  
    + 2 sequencing files per sample: "Forward" and "Reverse"  
* The CDH1_sequence.fasta file corresponding to the reference DNA sequence  
    + The reference genome is "BSgenome.Hsapiens.UCSC.hg19"   
  
Use these files and information as input in the "individual analysis" tab and run analysis   

### Grouped analysis 
After running all the individual analyses of samples:  
* The grouped analysis can be runned from the "grouped analysis" tab  
* Select the folder corresponding to the outputs of the individual results  
* Select the experiment information and the plot parameters of your choice  
  
Run analysis

### Multiple analyses
To launch several the analyses in one click:
* Copy the "multiple_individual_analyses_table-CDH1" and "multiple_grouped_analyses_table-CDH1" files before modifying them
* Use the orginal files as an example
* Modify the copied files:  
    + Change the file path of reference sequence .fasta files and sequencing .ab1 files with your own file path of the same files present in the "input" folder
    + Change the parameters of your choice.   
  
Use the modified table files as input in the "multiple analyses" tab and run analyses


## "Test sequence" data set

### Individual analysis  
The methylation data files, resulting from the individual analyses, have been generated *in silico* 
* Move the folders 
"ABSP/examples/results/Example data/Test sequence" to the "ABSP/results" folder with the same structure 
-> "ABSP/results/Example data/Test sequence"  
  
(If the "ABSP/results" folder does not already exist, it can be created)   
(The "ABSP/" folder corresponds to the main ABSP folder, with the version identifier in the folder name)  

### Grouped analysis  
Once the methylation data files from individual analyses are in the "ABSP/results" folder:
* The grouped analysis can be runned from the "grouped analysis" tab
* Select the folder corresponding to the outputs of the individual results 
* Select the experiment information and the plot parameters of your choice  
  
Run analysis

### Multiple analyses  
Once the methylation data files from individual analyses are in the "ABSP/results" folder:
* The multiple grouped analyses can be runned from the "multiple analysis" tab
* The "multiple_grouped_analyses_table-Test sequence" file can be modified or directly used as input  
  
Use it as input in the "multiple analyses" tab and run analyses
