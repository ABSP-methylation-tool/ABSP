# ABSP - Analysis of Bisulfite Sequencing PCR

[<img src="./documents/logo.png" width="200" />](./documents/logo.png)

ABSP, standing for *"Analysis of Bisulfite Sequencing PCR"*, is an R based tool to analyze CpG methylation profiles using data from Bisulfite Sequencing PCR (BSP) experiment results. 
It was developed to help researchers to estimate and compare methylation percentages of a DNA region studied using BSP experiments. 
It provides a complete automated workflow, from trace file sequencing results to data visualization and statistics.

For more detailed information, please read the "ABSP User Guide" document.


## Quick start guide

### Installation

* Download all the ABSP files ("Code", "Download Zip")
* Install R and RStudio

### Open the ABSP app

* Open the "ABSP Rproject.Rproject" file with Rstudio (it might take a few minutes to open)
* Open the "app.R" file with Rstudio
* Click on the "Run App" button on the top right corner (select "Run external" to open it in web browser)

### Analysis with ABSP

#### Individual sample analysis: 
* In the "Individual analysis" tab, fill the entries and add input files (reference sequence .fasta file and sequencing .ab1 files)
* Run the analysis to compute the CpG methylation levels of your sample 
* View the results in the .html report generated in your "reports" folder
* Look for output files in your "results" folder

#### Grouped samples analysis:
* Once all of your samples have been individually analyzed, in the "Grouped analysis" tab, fill the entries and choose your parameters for plotting
* Run the analysis to gather samples, generate visualization plots (lollipop plots) and compare methylation data between groups by comparative statistics
* View the results in the .html report generated in your "reports" folder
* Look for output files in your "results" folder

[<img src="./examples/results/Example data/Test sequence/grouped_results_direct/lollipop_plots/plots_replicates/lollipop_replicates_as-is_proportional.png" width="400" />](./examples/results/Example data/Test sequence/grouped_results_direct/lollipop_plots/plots_replicates/lollipop_replicates_as-is_proportional.png)



#### Launch multiple analysis:  
The "Multiple analysis" tab is useful to analyze multiple samples and/or for multiple grouping analysis, launched in one click, using as input tables filled with the required input entries.


## Availability

Currently, ABSP is not functional on MacOS devices due to compatibility issues of the renv package with the 4.2.0 R version on MacOS.

## License

