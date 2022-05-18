How to use the examples of data ?

Two data sets are available:
- A real sequence named "CDH1", direct-BSP experiment
- A mock sequence named "Test sequence", both direct-BSP and cloning-BSP experiments


* "CDH1" data set

- Individual analysis:
Find in the "input" folder the .ab1 sequencing files of CDH1
2 groups: "High methylated DNA" and "Low methylated DNA"
3 replicates: "rep1", "rep2" and "rep3"
2 sequencing files per sample: "Forward" and "Reverse"
and the CDH1_sequence.fasta file corresponding to the reference DNA sequence
the reference genome is "BSgenome.Hsapiens.UCSC.hg19"
Use these files and information as input in the "individual analysis" tab.

- Grouped analysis:
After running all the individual analyses of samples, the grouped analysis can be runned from the "grouped analysis" tab.
Select the folder corresponding to the outputs of the individual results.
Select the plot parameters of your choice.

- Multiple analyses:
Copy the "multiple_individual_analyses_table-CDH1" and "multiple_grouped_analyses_table-CDH1" files before modifying them.
Use the orginal files as an example.
Modify the copied files: 
Change the file path of reference sequence .fasta files and sequencing .ab1 files with your own file path of the same files present in the "input" folder. 
Change the parameters of your choice.
Use the modified table files as input in the "multiple analyses" tab.


* "Test sequence data" set

- Individual analysis:
The methylation data files, resulting from the individual analyses, have been generated in silico.
The folders 
"ABSP/examples/results/Example data/Test sequence/individual_results_cloning" and "ABSP/examples/results/Example data/Test sequence/individual_results_direct" must be moved to the "ABSP/results" folder with the same structure 
-> "ABSP/results/Example data/Test sequence/individual_results_cloning" and "ABSP/results/Example data/Test sequence/individual_results_direct"

- Grouped analysis:
Once the methylation data files from individual analyses are in the "ABSP/results" folder, the grouped analysis can be runned from the "grouped analysis" tab.
Select the folder corresponding to the outputs of the individual results.
Select the plot parameters of your choice.

- Multiple analyses:
Once the methylation data files from individual analyses are in the "ABSP/results" folder, the "multiple_grouped_analyses_table-Test sequence" file can be modified or directly used as input in the "multiple analyses" tab.
