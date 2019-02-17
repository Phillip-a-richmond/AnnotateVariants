source /opt/tools/hpcenv.sh
source activate SVcalling

BAM1=/mnt/causes-vnx2/TIDE/PROCESS/GENOME_TIDEX/T308/TestSmoove/12D3695_dedup.bam
BAM2=/mnt/causes-vnx2/TIDE/PROCESS/GENOME_TIDEX/T308/TestSmoove/12D3698_dedup.bam
BAM3=/mnt/causes-vnx2/TIDE/PROCESS/GENOME_TIDEX/T308/TestSmoove/17D8290_dedup.bam
NAME=T308_Trio_SecondTry
OUTDIR=/mnt/causes-vnx2/TIDE/PROCESS/GENOME_TIDEX/T308/${NAME}_Smoove
PROCS=16
GENOME_FASTA=/mnt/causes-vnx1/GENOMES/GSC/GRCh37-lite-smooveTest.fa 

mkdir $OUTDIR

smoove call \
	--name $NAME -p $PROCS \
	--outdir $OUTDIR \
	--genotype -f $GENOME_FASTA \
	$BAM1 $BAM2 $BAM3

