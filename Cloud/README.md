# Cloud
> This directory contains work from using the RONIN cloud at UBC Advanced Research Computing. 

## Overview

The analysis carried out in this directory includes:
1. Identifying the 1kg trios and putting together code that will download CRAMs (mapped to GRCh38) for each trio.
2. Running SV calling on the trios with Manta and Smoove, using conda (for Manta), and docker (for Smoove). 
3. Running Excord on the 1kg data


## Setting up

I've got a couple of set-up deployments that will enable analysis of data from the 1kG project.

Once a new cluster is set up:
1. Clone this repo.
```
git clone https://github.com/Phillip-a-richmond/AnnotateVariants
```

2. Run this script:
```
cd AnnotateVariants/Cloud/
bash DeployBasicGenomics.sh
```
NOTE: It will take a little while. This script will get you:
    - Miniconda
    - Mamba -> Manta, Samtools, Htslib, FGBio, Picard

3. Now you're good to go!


## Running analysis

### DeepVariant vs. DeepTrio

This analysis was a doozie. I gotta say I did not properly calculate the runtime, I probably shouldn't have run this on CPU, and holy hell was it expensive. 

That being said, I think I still got some valuable content out of this. 

Here's what I wanted to have:
Trio CRAMs pulled from S3 to onboard NVMe drives. For each trio, run DeepVariant, and run DeepTrio. Save the VCFs, and GVCFs from each, and then upload to S3. 

When I started running this I was a bit aggressive. I had run it once, didn't track too closely the time required to run it (hard to work w/o seff command, RONIN Ganglia is garbage). 


### Manta + Smoove

This one went a little smoother than the DeepVariant business. I have here calculated the runtime for each of the runs I had. Some were 3hr, some were >1d. The volatility in runtime made the budget hard here.

Need to dig into that more.



### Excord

Running this now (2022-04-11).







