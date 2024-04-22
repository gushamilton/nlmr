

module load apps/plink
module load apps/bgen


UKB_MERGED=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/directly_genotyped/released/2017-07-04/data/derived/merged_chr1-22
UKBIOBANK_DATA=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen
SCRATCH_NJT=/mnt/storage/private/njt_grp_space/FH_UKB

cd /user/work/fh6520/bmi/

cut -f 1 bmi_locke.tsv > rsid_list_genome.tsv

cmd=""
for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22
do
  bgenix -g $UKBIOBANK_DATA/data.chr$i.bgen -incl-rsids rsid_list_genome.tsv > chr_${i}.bgen
  cmd=$cmd"chr_${i}.bgen "
done

# Combine the .bgen files for each chromosome into one
cat-bgen -g  $cmd -og initial_chr.bgen -clobber
# Write index file .bgen.bgi
bgenix -g initial_chr.bgen -index -clobber

# Remove the individual chromosome files
for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22
do
  rm chr_${i}.bgen
done

plink2 --bgen initial_chr.bgen --sample $UKBIOBANK_DATA/data.chr1-22_plink.sample --rm-dup exclude-all --score bmi_locke.tsv 1 2 4


# now generate individual PRS:


# Define the input files
BGEN_FILE="initial_chr.bgen"
SAMPLE_FILE="$UKBIOBANK_DATA/data.chr1-22_plink.sample"
SCORE_FILE="bmi_locke.tsv"

# Extract SNP names from the score file (assuming SNP names are in the first column)
cut -f1 $SCORE_FILE > snp_list.txt

# Loop through each SNP
while read SNP; do
  # Generate a temporary score file excluding the current SNP
  awk -v snp=$SNP '$1 != snp' $SCORE_FILE > temp_score_file.tsv
  
  # Run PLINK to calculate the PRS without the current SNP
  plink2 --bgen $BGEN_FILE --sample $SAMPLE_FILE --rm-dup exclude-all --score temp_score_file.tsv 1 2 4 --out prs_wo_$SNP
  
done < snp_list.txt

# Clean up temporary files
rm temp_score_file.tsv snp_list.txt

# then plink to generate text datasxe

  plink2 --bgen $BGEN_FILE --sample $SAMPLE_FILE --rm-dup exclude-all --recode AD
