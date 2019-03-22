#Install Guide

## Requirements: 
- Miniconda or Anaconda3
- snakemake installed 
- Java 1.8

## Steps 
1) install base packages: 
`conda remove --name AnnotateVariants --all`
`conda env create -f Install/complete_install/environment.yml`

2) change snakemake config to add install path 
3) run snakemake script
    - installs java jar
    - installs gemini
    - installs GATK
    - downloads dbs
