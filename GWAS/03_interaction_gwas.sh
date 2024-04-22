#!/bin/bash
#SBATCH --job-name=qc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --account=sscm013902

# do step one - have to do this across whole genome

# RTUN O N BC%


#UKB_MERGED=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/directly_genotyped/released/2017-07-04/data/derived/merged_chr1-22
#UKBIOBANK_DATA=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen
#SCRATCH_NJT=/mnt/storage/private/njt_grp_space/FH_UKB


UKB_MERGED=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/directly_genotyped/released/2017-07-04/data/derived/merged_chr1-22
UKBIOBANK_DATA=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/
SCRATCH_NJT=/user/work/fh6520/stratified

cd /user/work/fh6520/stratified


# ./plink2 \
#   --bfile $UKB_MERGED/chr1-22_merged \
#   --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \
#   --mind 0.1 \
#   --write-snplist --write-samples --no-id-header \
#   --out qc_pass


conda activate regenie_env





/user/home/fh6520/.conda/envs/regenie_env/bin/regenie \
  --step 1 \
  --bed $UKB_MERGED/chr1-22_merged \
  --extract $SCRATCH_NJT/qc_pass.snplist \
  --keep $SCRATCH_NJT/qc_pass.id \
  --covarFile data.covariates.plink.txt \
  --phenoFile phenotype_interact.tsv \
  --covarCol PC{1:10},sex,chip \
  --bsize 500 \
  --lowmem \
  --lowmem-prefix $SCRATCH_NJT/tmpdir \
  --out $SCRATCH_NJT/results/interaction \
  --force-qt \
  --gz

for chr in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22;

do echo "chr${chr}"

cat > $chr.sh <<'endmsg'
#!/bin/bash
#SBATCH --job-name=${chr}
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=16
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH --account=sscm013902



# UKB_MERGED=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/directly_genotyped/released/2017-07-04/data/derived/merged_chr1-22
# UKBIOBANK_DATA=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen
# INPUT_DIR=/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/phenotypes/fh6520/input
# SCRATCH_NJT=/mnt/storage/private/njt_grp_space/FH_UKB




UKB_MERGED=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/directly_genotyped/released/2017-07-04/data/derived/merged_chr1-22
UKBIOBANK_DATA=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/
SCRATCH_NJT=/user/work/fh6520/stratified
cd /user/work/fh6520/stratified

module load apps/regenie

endmsg


echo >> $chr.sh

echo "regenie \
  --step 2 \
  --bgen \$UKBIOBANK_DATA/data.chr$chr.bgen \
  --sample \$UKBIOBANK_DATA/data.chr1-22.sample  \
  --pred \$SCRATCH_NJT/results/interaction_pred.list\
  --ref-first \
  --minMAC 20 \
    --covarFile covariates_interact.tsv \
  --phenoFile phenotype_interact.tsv \
  --covarCol PC{1:10},sex,chip \
  --bsize 400  --no-split \
  --interaction bmi_prs \
  --no-condtl \
  --out $SCRATCH_NJT/results \
  --gz" >>  $chr.sh ;

sbatch $chr.sh
chmod +x $chr.sh
done





