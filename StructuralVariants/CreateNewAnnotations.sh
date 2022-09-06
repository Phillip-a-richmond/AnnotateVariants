# This script details how I went from the downloaded GRCh37 file, into GRCh38 files.

# The sub-files in this process are for either alternative coordinates with insertions (ins.shifted) or just reformatted
# This script is meant to run from:
cd /mnt/common/DATABASES/REFERENCES/GRCh37/GNOMADSV/V2.1/

# Step 1 - run python script to split out and reformat 2 sub files, shifted and reformatted
python CreateAltPositionGnomAD.py -I gnomad_v2.1_sv.sites.bed -O gnomad_v2.1_sv.sites.bed.altered

# This makes these files:
#-rwxrwx--- 1 prichmond domain users  12M Aug 25 14:03 gnomad_v2.1_sv.sites.bed.altered.ins.shifted.bed
#-rwxrwx--- 1 prichmond domain users  12M Aug 25 14:03 gnomad_v2.1_sv.sites.bed.altered.ins.shifted.chr.bed
#-rwxrwx--- 1 prichmond domain users   96 Aug 25 14:03 gnomad_v2.1_sv.sites.bed.altered.ins.shifted.chr.bed.header
#drwxrwx--- 2 prichmond domain users 1.2K Aug 25 14:03 .
#-rwxrwx--- 1 prichmond domain users  42M Aug 25 14:03 gnomad_v2.1_sv.sites.bed.altered.reformatted.bed
#-rwxrwx--- 1 prichmond domain users  43M Aug 25 14:03 gnomad_v2.1_sv.sites.bed.altered.reformatted.chr.bed
#-rwxrwx--- 1 prichmond domain users  36M Aug 25 14:03 gnomad_v2.1_sv.sites.bed.altered.reformatted.chr.bed.header

# Step 2 - liftOver
# FOr details on this go look at that script
sh /mnt/common/Precision/LiftOver/RunLiftOver.sh

# This made these files:
# In: /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/
# gnomad_v2.1_sv.sites.altered.ins.shifted.chr.hg38.bed
# gnomad_v2.1_sv.sites.altered.reformatted.chr.hg38.bed


# Step 3 - Re-header lifted-over files
cat /mnt/common/DATABASES/REFERENCES/GRCh37/GNOMADSV/V2.1/gnomad_v2.1_sv.sites.bed.altered.reformatted.chr.bed.header \
	/mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/gnomad_v2.1_sv.sites.altered.reformatted.chr.hg38.bed \
	> /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/gnomad_v2.1_sv.sites.bed.altered.reformatted.chr.hg38.headered.bed


cat /mnt/common/DATABASES/REFERENCES/GRCh37/GNOMADSV/V2.1/gnomad_v2.1_sv.sites.bed.altered.ins.shifted.chr.bed.header \
	/mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/gnomad_v2.1_sv.sites.altered.ins.shifted.chr.hg38.bed \
	> /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/gnomad_v2.1_sv.sites.bed.altered.ins.shifted.chr.hg38.headered.bed


# Step 4 trim the 'chr'
sed -e 's/^chr//g' /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/gnomad_v2.1_sv.sites.bed.altered.ins.shifted.chr.hg38.headered.bed \
	> /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/gnomad_v2.1_sv.sites.bed.altered.ins.shifted.hg38.bed

sed -e 's/^chr//g' /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/gnomad_v2.1_sv.sites.bed.altered.reformatted.chr.hg38.headered.bed \
	> /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/gnomad_v2.1_sv.sites.bed.altered.reformatted.hg38.bed

# Step 5 split into sub-files for AnnotSV
# Full Reformatted
# gnomAD_AF
cut -f1,2,3,9 /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/gnomad_v2.1_sv.sites.bed.altered.reformatted.hg38.bed > \
	/mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/UserGNOMADSV-AF.bed
# gnomAD_NHOMALT
cut -f1,2,3,10 /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/gnomad_v2.1_sv.sites.bed.altered.reformatted.hg38.bed > \
	/mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/UserGNOMADSV-NHOMALT.bed


# Shifted Insertions
# gnomAD_AF
cut -f1,2,3,9 /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/gnomad_v2.1_sv.sites.bed.altered.reformatted.hg38.bed > \
        /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/UserGNOMADSV-AF-insertions.bed
# gnomAD_NHOMALT
cut -f1,2,3,10 /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/gnomad_v2.1_sv.sites.bed.altered.reformatted.hg38.bed > \
        /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/UserGNOMADSV-NHOMALT-insertions.bed


# Step 6 copy over to AnnotSV
ANNOTSVDIR=/mnt/common/Precision/AnnotSV/share/AnnotSV/Annotations_Human/Users/GRCh38/
# FtIncludedInSV
# Full reformatted
cp /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/UserGNOMADSV-AF.bed  \
	$ANNOTSVDIR/FtIncludedInSV/UserGNOMADSV-AF.bed
cp /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/UserGNOMADSV-NHOMALT.bed  \
	$ANNOTSVDIR/FtIncludedInSV/UserGNOMADSV-NHOMALT.bed

# Shifted insertions
cp /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/UserGNOMADSV-AF-insertions.bed  \
	$ANNOTSVDIR/FtIncludedInSV/UserGNOMADSV-AF-insertions.bed
cp /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/UserGNOMADSV-NHOMALT-insertions.bed  \
	$ANNOTSVDIR/FtIncludedInSV/UserGNOMADSV-NHOMALT-insertions.bed

#SVincludedInFt
# Full reformatted
cp /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/UserGNOMADSV-AF.bed  \
	$ANNOTSVDIR/SVincludedInFt/UserGNOMADSV-AF.bed
cp /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/UserGNOMADSV-NHOMALT.bed  \
	$ANNOTSVDIR/SVincludedInFt/UserGNOMADSV-NHOMALT.bed

# Shifted insertions
cp /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/UserGNOMADSV-AF-insertions.bed  \
	$ANNOTSVDIR/SVincludedInFt/UserGNOMADSV-AF-insertions.bed
cp /mnt/common/DATABASES/REFERENCES/GRCh38/GNOMADSV/V2.1/UserGNOMADSV-NHOMALT-insertions.bed  \
	$ANNOTSVDIR/SVincludedInFt/UserGNOMADSV-NHOMALT-insertions.bed


