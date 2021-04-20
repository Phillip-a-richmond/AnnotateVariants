
echo "Mitochondrial Variant Analysis Started"
date

# Run MToolBox for mitochondrial variant analysis 
SAMPLE='EPGEN055-01'
MTOOLBOX_PATH=/mnt/common/MToolBox/ 
PATH=$MTOOLBOX_PATH:$MTOOLBOX_PATH:$PATH 
MTOOLBOX_WORKING_DIR=$WORKING_DIR/MToolBox_${SAMPLE}/ 
mkdir -p $MTOOLBOX_WORKING_DIR/ 
 
MTOOLBOX_CONFIG_FILE_ORIGINAL='/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/MToolBox_config_files/MToolBox_RSRS_config_with_markdup_and_indelrealign_RvdL.sh'
MTOOLBOX_CONFIG_FILE_ORIGINAL_BASENAME=$(basename $MTOOLBOX_CONFIG_FILE_ORIGINAL) 
MTOOLBOX_CONFIG_FILE=$MTOOLBOX_WORKING_DIR/$MTOOLBOX_CONFIG_FILE_ORIGINAL_BASENAME 
 
#edit the MToolBox config file template so that is specifies the MToolBox results/working directory for the current analysis 
cp $MTOOLBOX_CONFIG_FILE_ORIGINAL $MTOOLBOX_CONFIG_FILE 
sed -i "s#^output_name\=\.#output_name=$MTOOLBOX_WORKING_DIR#" $MTOOLBOX_CONFIG_FILE 
 
#link the raw fastq files to the MToolBox working directory, and name them as required by MToolBox: \<sample\_name\>.R1.fastq, \<sample\_name\>.R2.fastq 
ln -sf ${FASTQR1} $MTOOLBOX_WORKING_DIR/${SAMPLE}.R1.fastq.gz 
ln -sf ${FASTQR2} $MTOOLBOX_WORKING_DIR/${SAMPLE}.R2.fastq.gz 
 
echo "Changing working directory to $MTOOLBOX_WORKING_DIR and running MToolBox..." 
PWD_CURRENT=`pwd` 
cd $MTOOLBOX_WORKING_DIR 
$MTOOLBOX_PATH/MToolBox/MToolBox.sh -i ${MTOOLBOX_CONFIG_FILE} 
echo "Changing working directory to back to $PWD_CURRENT..." 
cd $PWD_CURRENT 
 
