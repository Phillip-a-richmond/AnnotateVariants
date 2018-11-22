# By default, MToolBox adopts the RSRS (Reconstructed Sapiens Reference Sequence, PMID: 22482806) as mitochondrial reference genome and hg19 as nuclear reference genome
input_type=fastq
ref=RSRS
vcf_name=RSRS_with_markdup_and_indelrealign
output_name=. # output dir within current working dir

UseMarkDuplicates=true
UseIndelRealigner=true
