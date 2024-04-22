
for chr in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22;

do echo "chr${chr}"

cat > $chr.sh <<'endmsg'
#!/bin/bash
#SBATCH --job-name=${chr}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=4:00:00
#SBATCH --mem=8G
#SBATCH --account=sscm013902



# UKB_MERGED=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/directly_genotyped/released/2017-07-04/data/derived/merged_chr1-22
# UKBIOBANK_DATA=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen
# INPUT_DIR=/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/phenotypes/fh6520/input
# SCRATCH_NJT=/mnt/storage/private/njt_grp_space/FH_UKB



UKB_MERGED=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/directly_genotyped/released/2017-07-04/data/derived/merged_chr1-22
UKBIOBANK_DATA=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/
SCRATCH_NJT=/user/work/fh6520/stratified
cd /user/work/fh6520/stratified

cd /user/work/fh6520/stratified
module load apps/regenie

endmsg


echo >> $chr.sh

echo "regenie \
  --step 2 \
  --bgen \$UKBIOBANK_DATA/data.chr$chr.bgen \
  --sample \$UKBIOBANK_DATA/data.chr1-22.sample  \
  --pred \$SCRATCH_NJT/pre_sub/strata_pred.list\
  --ref-first \
  --minMAC 20 \
--covarFile data.covariates.plink.txt \
   --phenoFile stratified_pre_submission.tsv \
  --covarCol PC{1:10},sex,chip \
  --bsize 400 \
  --force-qt \
  --out \$SCRATCH_NJT/results/$chr \
  --gz" >>  $chr.sh ;

# sbatch $chr.sh
chmod +x $chr.sh
done
