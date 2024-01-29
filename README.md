### Instruction for 16S rRNA pipeline on CPU server

#####################################################################################
### PREREQUISITES

1. Create a parent directory for your 16S rRNA data:

	`mkdir my_folder`	#where my_folder is a name of your parent directory
	`cd my_folder`
	
2. Create a folder all_data for your *.fastq.gz files (IMPORTANT STEP):

	`mkdir all_data`	
	
3. Copy all your fastq.gz files in all_data folders. IMPORTANT: copy only your R1 fastq.gz files and the endings of your files MUST be fastq.gz


####################
### DRY RUN

1. Activate conda base environment:

	`conda activate base`

2. Activate conda snakemake environment:

	`conda activate snakemake`

3. To check if all prerequisetes were done correct run the follow command:

	`snakemake --dry-run --snakefile 16S_snakefile.smk`
	
4. Use the options from section RUN PIPELINE on steps 2 and 3

5. If the dry run is completed successfully go to the section RUN PIPELINE step 1.


####################
### RUN PIPELINE

1. Go to screen mode:

	screen (press ENTER) # more information about screen: https://linuxize.com/post/how-to-use-linux-screen/

2. Activate conda base environment:

	`conda activate base`

3. Activate conda snakemake environment:

	`conda activate snakemake`
	
4. Run pipeline 16S_snakefile_for_all.smk:
	
	`snakemake --cores 20 --rerun-incomplete --snakefile 16S_snakefile.smk`
	
5. Pipeline asks you an absolute pathway for your parent directory:

	Set a path to directory with samples:
	
	Write the absolute pathway to your parent directory
	EXAMPLE: "Set a path to directory with samples:/home/lam2/Seaweed" (IMPORTANT: NO SLASH AT THE END)
	
6. Pipeline asks you a threshold to remove noise ASVs: 

	Set a threshold for filtration OTU frequency:
	
	Write a number indicating the threshold (Usually, it is 1, but if you have number of sample more than 500 it is better to increase the threshold up to 10)
	EXAMPLE: Set a threshold for filtration OTU frequency: 1
	
7. Close your screen entering CTRL+A and then D
	
8. Wait for the results (4-24 h depending on the amount of files).

9. In your parent folder the tables with OTU frequencies and taxonomy is located in "absolute_path_to_your_parent_folder"/results. These tables are used for further analysis.

10. To open your screen again write the following command and press TAB, and then ENTER:

	`screen -r` 
	
11. To kill your screen, open your screen and press CTRL+A and then press K

### 16S results
The folder contains OTUs and taxonomy tables after running 16S_snakefile.smk

### R scripts
The folder contains R scripts used for the 16S analysis including some test R scripts.

### Tables
The folder contains tables with kit ranking data, Shannon index vs contamination, and filtered OTUs with different decontam parameters.
	
####################################################################################
### Analysis of the results
Script phyloseq.R demonstrates the example of analysis of the result tables.
More tutorials for R analysis of OTU tables:

Phyloseq: https://micca.readthedocs.io/en/latest/phyloseq.html

Alpha Diversity: https://microbiome.github.io/course_2021_radboud/alpha-diversity.html

Beta Diversity: https://microbiome.github.io/course_2021_radboud/beta-diversity.html

Differential abundance analysis: https://microbiome.github.io/course_2021_radboud/differential-abundance-analysis.html
