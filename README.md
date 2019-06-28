# AnnotateVariants
#### Phillip Richmond (Phillip.A.Richmond@gmail.com)

> This is a pipeline for variant annotation in the diagnosis of rare genetic disorders. It relies on open source data and has instructions for software installs.

## Overview
1. [Pipeline Summary & Diagram](#pipeline-summary-and-diagram)

2. [Set-up](#set-Up) 
	- Prepare Datasets and databases
	- Install Necessary Tools

3. [Run Test](#run-test)

4. [Run Sample](#run-sample)


### Pipeline Summary And Diagram
This pipeline was designed by Phillip Richmond in order to analyze & prioritize variants in rare genetic disease cases. Currently, the pipeline uses the following list of software in order to accomplish this task, much to the thanks of tools produced and maintained by the lab of Aaron Quinlan:
+ GEMINI
+ VCFAnno
+ VCF2DB

Furthermore, this pipeline utilizes open source datasets within it's annotation framework, including:
+ CADD 
+ gnomAD
+ OMIM\*
+ ClinVar
+ UCSC RefGene
+ Entrez Gene Summary
+ HPO Term Mapping
+ MeSHOPs
+ pLI
+ RVIS
+ FATHMM-XF
+ Eigen
+ FunSeq2
+ Platinum Genomes ConfidentRegions
+ UCSC Segmental Duplications

\* OMIM requires a license for use of the API/downloadable databases, which must be applied for through their website.

Currently, the pipeline is hard coded for a specific cluster that uses the Torque-Moab scheduler. However, I will expand upon this to include other schedulers such as SLURM. Also, generalizing for software install locations, or developing a single install-script via bioconda will also be performed later in 2018.


![](https://github.com/Phillip-a-richmond/AnnotateVariants/blob/master/Figures/Figure3-NewInterpretationPipeline.png)



### Set-up
*THERE IS A LOT OF SETTING UP TO DO!*  

NOTE - (2018-09-10) Updates to this are coming. We will have a conda install environment, and a unified acquisition script for the databases.

But, once you get set up, then things run nice and smooth.

1. Install Necessary software, details in InstallTools.sh
+ bgzip and tabix
+ vt
+ vcftools/bcftools
+ snpEff
+ vcfanno
+ vcf2db
+ gemini
+ MToolBox
+ In-house Scripts 


2. Prepare Third-party Datasets/databases, there are a few scripts to help do this within [here](https://github.com/Phillip-a-richmond/AnnotateVariants/tree/master/GetThirdPartyDatabases):
+ All the gemini databases
+ Polyphen2  https://github.com/quinlan-lab/pathoscore/blob/master/score-sets/GRCh37/polyphen2/make.sh
+ CADD  http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz 
+ ReMM  http://remm.visze.de/files/ReMM.v0.3.1.tsv.gz
+ gnomAD http://gnomad.broadinstitute.org/downloads


### Run Test
You can run a test analysis by following the instructions within the [Test](https://github.com/Phillip-a-richmond/AnnotateVariants/tree/master/Test) directory.

### Run Sample


### Known Issues
+ There is no inheritance model for de-novo + compound het (GEMINI Limitation). E.g. where a de novo is the second variant of a compound het pair
+ There is no guarantee/unit tests for deletion + snv (old pipeline hemizygous?)
+ Currently hard-coded, and not system agnostic. As the pipeline stands it needs a new version with configuration files before it can be ported to a new system (evident with my hard-coded paths).



### Improvements
+ GEMINI ROH - Needs to wait on Brent Pedersen to fix
+ Better management of duo and singleton analyses. For now, they just can be mined from the General Damaging tab
+ Better utilization of duo-affected pairs

### The Future 
+ TIDEX-Tool (jacques & alice) for better analysis

### Additional Contributors
Bhavi Modi, Robin van der Lee, and Solenne Correard are contributors to this project development within the Wasserman lab at BCCHR/CMMT in Vancouver, BC, Canada.



