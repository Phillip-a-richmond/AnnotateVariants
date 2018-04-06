#!/bin/bash

#SBATCH --account=rrg-wyeth

## Mail Options
#SBATCH --mail-user=prichmond@cmmt.ubc.ca
#SBATCH --mail-type=ALL
#SBATCH --mem=100G
#SBATCH --cpus-per-task=32
#SBATCH --time=2-0:00
#SBATCH --nodes=1
## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error


WORKING_DIR=/scratch/richmonp/DATABASES/
SNPEFFCONFIG=/home/richmonp/TOOLS/anaconda/pkgs/snpeff-4.3.1p-1/share/snpeff-4.3.1p-1/snpEff.config
NORMVCFGZ=Clinvar_norm_eff.vcf.gz
NORMVCF=Clinvar_norm_eff.vcf
ANNONORMVCF=Clinvar_norm_eff_anno.vcf
TMPDIR=$WORKING_DIR/tmp
GENOME_FASTA=/scratch/richmonp/GENOME/ucsc.hg19.fasta
rm -rf $TMPDIR
mkdir $TMPDIR

cd $WORKING_DIR

# Get ClinVar Summary
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
gunzip -c ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz


# Get ClinVar
wget -c ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
gunzip -c clinvar.vcf.gz > clinvar.vcf
# Fix ClinVar
python AddCHRtoVCF.py clinvar.vcf clinvar_chromfix.vcf 
grep -v "NW_003" clinvar_chromfix.vcf > clinvar_fixed.vcf 
bgzip clinvar_fixed.vcf
tabix clinvar_fixed.vcf.gz
VCF=clinvar_fixed.vcf

zless $VCF \
	| vt decompose -s - \
	| vt normalize -r $GENOME_FASTA - \
	| snpEff -c $SNPEFFCONFIG GRCh37.75 - \
	| bgzip -c > $NORMVCFGZ
tabix -p vcf $NORMVCFGZ

gunzip -c $NORMVCFGZ > $NORMVCF

table_annovar.pl  $NORMVCF \
	/scratch/richmonp/TOOLS/annovar/humandb -buildver hg19 -out $ANNONORMVCF \
	-nastring . -remove -vcfinput -otherinfo -protocol refGene,cytoBand,revel,gnomad_genome,clinvar_20170130,mitimpact24,eigen,avsnp147,fathmm,dbnsfp33a \
	-operation  g,r,f,f,f,f,f,f,f,f

