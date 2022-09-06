import sys,argparse

def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument("-I","--Input",help="Input AnnotSV tsv",required=True)
    parser.add_argument("-O","--Output",help="Filtered AnnotSV file")
    parser.add_argument("-F","--FatherID",help="Father ID", required=True)
    parser.add_argument("-M","--MotherID",help="Mother ID", required=True)
    parser.add_argument("-P","--ProbandID",help="Proband ID", required=True)
    args = parser.parse_args()

    return args.Input,args.Output,args.FatherID,args.MotherID,args.ProbandID

# This currently just filters for unique variants to a single sample (should be mostly de novos)
    # infile looks like:
    # AnnotSV_ID	SV_chrom	SV_start	SV_end	SV_length	SV_type	Samples_ID	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	EPGEN012-01_GRCh38	EPGEN012-02_GRCh38	EPGEN012-03_GRCh38	EPGEN029-01_GRCh38	EPGEN029-02_GRCh38	EPGEN029-03_GRCh38	EPGEN067-01_GRCh38	EPGEN067-02_GRCh38	EPGEN067-03_GRCh38	EPGEN083-01_GRCh38	EPGEN083-02_GRCh38	EPGEN083-03_GRCh38	EPGEN090-01_GRCh38	EPGEN090-02_GRCh38	EPGEN090-03_GRCh38	EPGEN101-01_GRCh38	EPGEN101-02_GRCh38	EPGEN101-03_GRCh38	EPGEN159-01_GRCh38	EPGEN159-02_GRCh38	EPGEN159-03_GRCh38	EPGEN216-01_GRCh38	EPGEN216-02_GRCh38	EPGEN216-03_GRCh38	Annotation_mode	Gene_name	Gene_count	Tx	Tx_start	Tx_end	Overlapped_tx_length	Overlapped_CDS_length	Overlapped_CDS_percent	Frameshift	Exon_count	Location	Location2	Dist_nearest_SS	Nearest_SS_type	Intersect_start	Intersect_end	RE_gene	P_gain_phen	P_gain_hpo	P_gain_source	P_gain_coord	P_loss_phen	P_loss_hpo	P_loss_source	P_loss_coord	P_ins_phen	P_ins_hpo	P_ins_source	P_ins_coord	P_snvindel_nb	P_snvindel_phen	B_gain_source	B_gain_coord	B_loss_source	B_loss_coord	B_ins_source	B_ins_coord	B_inv_source	B_inv_coord	gnomAD_AF-altinsertion	InHouseDB_MEI_Count	gnomAD_NHOMALT	InHouseDB_SV_Freq	gnomAD_NHOMALT-altinsertion	gnomAD_AF	InHouseDB_SV_Count	InHouseDB_MEI_Freq	TAD_coordinate	ENCODE_experiment	gnomAD_AF-altinsertion	InHouseDB_MEI_Count	gnomAD_NHOMALT	InHouseDB_SV_Freq	gnomAD_NHOMALT-altinsertion	gnomAD_AF	InHouseDB_SV_Count	InHouseDB_MEI_Freq	GC_content_left	GC_content_right	Repeat_coord_left	Repeat_type_left	Repeat_coord_right	Repeat_type_right	Gap_left	Gap_right	SegDup_left	SegDup_right	ENCODE_blacklist_left	ENCODE_blacklist_characteristics_left	ENCODE_blacklist_right	ENCODE_blacklist_characteristics_right	ACMG	HI	TS	DDD_HI_percent	DDD_status	DDD_mode	DDD_consequence	DDD_disease	DDD_pmid	ExAC_delZ	ExAC_dupZ	ExAC_cnvZ	ExAC_synZ	ExAC_misZ	OMIM_ID	OMIM_phenotype	OMIM_inheritance	OMIM_morbid	OMIM_morbid_candidate	LOEUF_bin	GnomAD_pLI	ExAC_pLI	AnnotSV_ranking_score	AnnotSV_ranking_criteria	ACMG_class
    # 1_789465_224014609_DUP_1	1	789465	224014609	223225129	DUP	EPGEN029-01_GRCh38,EPGEN029-02_GRCh38,EPGEN029-03_GRCh38,EPGEN067-01_GRCh38,EPGEN067-02_GRCh38,EPGEN067-03_GRCh38,EPGEN083-01_GRCh38,EPGEN083-02_GRCh38,EPGEN083-03_GRCh38,EPGEN090-01_GRCh38,EPGEN090-02_GRCh38,EPGEN090-03_GRCh38,EPGEN101-01_GRCh38,EPGEN101-02_GRCh38,EPGEN101-03_GRCh38,EPGEN159-01_GRCh38,EPGEN159-02_GRCh38,EPGEN159-03_GRCh38,EPGEN216-01_GRCh38,EPGEN216-02_GRCh38,EPGEN216-03_GRCh38	384	N	<DUP>


def FilterAnnotSV(inputfile,outputfile,father_id,mother_id,proband_id):
    infile = open(inputfile,'r')
    outfile = open(outputfile,'w')
    for line in infile.readlines():
        cols = line.strip('\n').split('\t')
        Samples_ID=cols[6]
        # parse Samples_ID to remove any hit in more than one sample
        if ',' in Samples_ID:
            continue
        if father_id in Samples_ID:
            continue
        if mother_id in Samples_ID:
            continue

        # now we should have only single sample hits, removing anything in either parent
        outfile.write(line)




def Main():
    print("hey buddy")
    inputfilename,outputfilename,father_id,mother_id,proband_id = GetOptions()
    FilterAnnotSV(inputfilename,outputfilename,father_id,mother_id,proband_id)

    sys.exit()





if __name__=="__main__":
    Main()
