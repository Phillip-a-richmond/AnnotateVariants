# Running for this sample, Ideally would automate this from batch script to pull from file
Sample_ID=NA21143

## Set working space
Working_Dir=/scratch/$Sample_ID
mkdir -p $Working_Dir
cd $Working_Dir
CRAM_Dir=$Working_Dir/${Sample_ID}_CRAM
mkdir -p $CRAM_Dir
Sample_CRAM=${Sample_ID}.alt_bwamem_GRCh38DH.20150718.GIH.low_coverage.cram


Seq_Type=WGS
Sample_VCF=${Sample_ID}.vcf.gz
Sample_GVCF=${Sample_ID}.gvcf.gz
Output_Dir=$Working_Dir/Variants/
mkdir -p $Output_Dir
Final_Dir=/shared/DeepVariantOutput/

#### GRCh38 #### 
echo "GRCh38 genome"
Genome=GRCh38
Fasta_Dir=/shared/Genomes/
Fasta_file=GRCh38_full_analysis_set_plus_decoy_hla.fa


# Pull data from s3
aws s3 cp --no-sign-request s3://1000genomes/data/${Sample_ID}/alignment/ $CRAM_Dir

