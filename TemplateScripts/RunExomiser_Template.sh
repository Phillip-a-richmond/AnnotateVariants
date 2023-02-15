#!/bin/bash

#SBATCH --partition=defq

## Change to be your email address
#SBATCH --mail-user=email_address
#SBATCH --mail-type=ALL

## CPU Usage
## 16 Gb of RAM for the whole job
#SBATCH --mem=16G

## Using 2 CPUs
#SBATCH --cpus-per-task=2

## Running for a max time of 1 hour
#SBATCH --time=1:00:00

## Using only a single node
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

# Load java
## change to CVMFS
source /cm/shared/BCCHR-apps/env_vars/unset_BCM.sh
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh

## Load StdEnv/2020, then java
module load StdEnv/2020
module load java/17

# Set variables
EXOMISER=exomiser_jar
EXOMISER_DIR=exomiser_dir
WORKING_DIR=working_dir
cd $WORKING_DIR
FAMILY_ID=family_id
CONFIG_YML=$WORKING_DIR/${FAMILY_ID}_Exomiser_Template.yml

java -Xmx4g -jar -Djava.io.tmpdir=$PWD $EXOMISER_DIR/$EXOMISER \
         --analysis $CONFIG_YML \
        --spring.config.location=$EXOMISER_DIR
