#!/bin/bash
#########################################################################################################################################
#																	#
#  				  *** Downstream analysis ***										#
#										Script dated: 19 May, 2019			     	#
#										Version: 1.0						#
#  This phase includes the following steps:											     	#
#																	#
#   1. Combine all the gVCF.gz files into chromosome by chromosome. (USE: gatk GenomicsDBImport ) 					#
#   2. Generate GenotypeGVCFs for each chromosome files. (USE: gatk GenotypeGVCFs) 							#
#   3. Extract the SNPs and apply filter to the SNPs. (USE: (a) gatk SelectVariants and (b) gatk VariantFiltration			#
#   4. Extract the INDELs and apply filter to the INDELs. (USE: (a) gatk SelectVariants and (b) gatk VariantFiltration 			#
#																	#
# As a result, this script will generate SNPs and INDELs per chromosome across all the samples.						#
#########################################################################################################################################
#
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=115gb
#SBATCH --array=1-19
#SBATCH --output=logs/Chr-%A.%a.out
#SBATCH --error=logs/Chr-%A.%a.err
#SBATCH --time=1-00:00:00 
#SBATCH --job-name=GenomicDBImport

## Software Modules 
module load java/8u131
module load gatk/4.1.1.0

## Reference files
export REF=/ibex/scratch/kathirn/work/project/for_elodie/ref/CQ41.fa

## Sample/Data Variables 
export PROJECT=/ibex/scratch/reyel/1000quinoa/naga/GenomicDB/import
export INPUT_DIR=/ibex/scratch/reyel/1000quinoa/naga/VCF
export gVCF=${PROJECT}/gVCF;
export Chr=${PROJECT}/Chr;
export SNPs=${PROJECT}/SNPs;
export INDELs=${PROJECT}/INDELs;

mkdir -p $gVCF;
mkdir -p $Chr;
mkdir -p $SNPs;
mkdir -p $INDELs;
mkdir -p logs;
 
## Add all the *.g.vcf files for Genotyping 
set INPUT_TMP
for i in `ls -l ${INPUT_DIR}/*.g.vcf.gz | awk '{print $9}'`
do
 INPUT_TMP+="$i -V "; 
done
INPUT=${INPUT_TMP::-4}
#echo $INPUT

 ##### 1.CombineGVCF 2.GenotypeGVCFs 3.(a). Select SNP 3.(b). Apply Filter to SNPs 4.(a). Slect INDELs 4.(b). Apply filter to INDELs.

if [ ${SLURM_ARRAY_TASK_ID} -eq 1 ]
then
 time -p gatk GenomicsDBImport -V $INPUT --genomicsdb-workspace-path $gVCF/Chr00 --intervals Chr00 --reader-threads 32;
 time -p gatk GenotypeGVCFs --variant gendb://$gVCF/Chr00 -R $REF --output $Chr/Chr00.vcf.gz -L Chr00;
 time -p gatk SelectVariants --variant $Chr/Chr00.vcf.gz --reference $REF -select-type SNP --output $SNPs/raw_snps.Chr00.vcf;
 time -p gatk VariantFiltration --variant $SNPs/raw_snps.Chr00.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name snp_filter --output $SNPs/hard_filtered_snps.Chr00.vcf;
 time -p gatk SelectVariants --variant $Chr/Chr00.vcf.gz --reference $REF -select-type INDEL --output $INDELs/raw_indels.Chr00.vcf;
 time -p gatk VariantFiltration --variant $INDELs/raw_indels.Chr00.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name indel_filter --output $INDELs/hard_filtered_indels.Chr00.vcf;
elif [ ${SLURM_ARRAY_TASK_ID} -eq 2 ]
then
 time -p gatk GenomicsDBImport -V $INPUT --genomicsdb-workspace-path $gVCF/Chr01 --intervals Chr01 --reader-threads 32;
 time -p gatk GenotypeGVCFs --variant gendb://$gVCF/Chr01 -R $REF --output $Chr/Chr01.vcf.gz -L Chr01;
 time -p gatk SelectVariants --variant $Chr/Chr01.vcf.gz --reference $REF -select-type SNP --output $SNPs/raw_snps.Chr01.vcf;
 time -p gatk VariantFiltration --variant $SNPs/raw_snps.Chr01.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name snp_filter --output $SNPs/hard_filtered_snps.Chr01.vcf;
 time -p gatk SelectVariants --variant $Chr/Chr01.vcf.gz --reference $REF -select-type INDEL --output $INDELs/raw_indels.Chr01.vcf;
 time -p gatk VariantFiltration --variant $INDELs/raw_indels.Chr01.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name indel_filter --output $INDELs/hard_filtered_indels.Chr01.vcf;
elif [ ${SLURM_ARRAY_TASK_ID} -eq 3 ]
then
 time -p gatk GenomicsDBImport -V $INPUT --genomicsdb-workspace-path $gVCF/Chr02 --intervals Chr02 --reader-threads 32;
 time -p gatk GenotypeGVCFs --variant gendb://$gVCF/Chr02 -R $REF --output $Chr/Chr02.vcf.gz -L Chr02;
 time -p gatk SelectVariants --variant $Chr/Chr02.vcf.gz --reference $REF -select-type SNP --output $SNPs/raw_snps.Chr02.vcf;
 time -p gatk VariantFiltration --variant $SNPs/raw_snps.Chr02.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name snp_filter --output $SNPs/hard_filtered_snps.Chr02.vcf;
 time -p gatk SelectVariants --variant $Chr/Chr02.vcf.gz --reference $REF -select-type INDEL --output $INDELs/raw_indels.Chr02.vcf;
 time -p gatk VariantFiltration --variant $INDELs/raw_indels.Chr02.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name indel_filter --output $INDELs/hard_filtered_indels.Chr02.vcf;
elif [ ${SLURM_ARRAY_TASK_ID} -eq 4 ]
then
 time -p gatk GenomicsDBImport -V $INPUT --genomicsdb-workspace-path $gVCF/Chr03 --intervals Chr03 --reader-threads 32;
 time -p gatk GenotypeGVCFs --variant gendb://$gVCF/Chr03 -R $REF --output $Chr/Chr03.vcf.gz -L Chr03;
 time -p gatk SelectVariants --variant $Chr/Chr03.vcf.gz --reference $REF -select-type SNP --output $SNPs/raw_snps.Chr03.vcf;
 time -p gatk VariantFiltration --variant $SNPs/raw_snps.Chr03.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name snp_filter --output $SNPs/hard_filtered_snps.Chr03.vcf;
 time -p gatk SelectVariants --variant $Chr/Chr03.vcf.gz --reference $REF -select-type INDEL --output $INDELs/raw_indels.Chr03.vcf;
 time -p gatk VariantFiltration --variant $INDELs/raw_indels.Chr03.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name indel_filter --output $INDELs/hard_filtered_indels.Chr03.vcf;
elif [ ${SLURM_ARRAY_TASK_ID} -eq 5 ]
then
 time -p gatk GenomicsDBImport -V $INPUT --genomicsdb-workspace-path $gVCF/Chr04 --intervals Chr04 --reader-threads 32;
 time -p gatk GenotypeGVCFs --variant gendb://$gVCF/Chr04 -R $REF --output $Chr/Chr04.vcf.gz -L Chr04;
 time -p gatk SelectVariants --variant $Chr/Chr04.vcf.gz --reference $REF -select-type SNP --output $SNPs/raw_snps.Chr04.vcf;
 time -p gatk VariantFiltration --variant $SNPs/raw_snps.Chr04.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name snp_filter --output $SNPs/hard_filtered_snps.Chr04.vcf;
 time -p gatk SelectVariants --variant $Chr/Chr04.vcf.gz --reference $REF -select-type INDEL --output $INDELs/raw_indels.Chr04.vcf;
 time -p gatk VariantFiltration --variant $INDELs/raw_indels.Chr04.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name indel_filter --output $INDELs/hard_filtered_indels.Chr04.vcf;
elif [ ${SLURM_ARRAY_TASK_ID} -eq 6 ]
then
 time -p gatk GenomicsDBImport -V $INPUT --genomicsdb-workspace-path $gVCF/Chr05 --intervals Chr05 --reader-threads 32;
 time -p gatk GenotypeGVCFs --variant gendb://$gVCF/Chr05 -R $REF --output $Chr/Chr05.vcf.gz -L Chr05;
 time -p gatk SelectVariants --variant $Chr/Chr05.vcf.gz --reference $REF -select-type SNP --output $SNPs/raw_snps.Chr05.vcf;
 time -p gatk VariantFiltration --variant $SNPs/raw_snps.Chr05.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name snp_filter --output $SNPs/hard_filtered_snps.Chr05.vcf;
 time -p gatk SelectVariants --variant $Chr/Chr05.vcf.gz --reference $REF -select-type INDEL --output $INDELs/raw_indels.Chr05.vcf;
 time -p gatk VariantFiltration --variant $INDELs/raw_indels.Chr05.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name indel_filter --output $INDELs/hard_filtered_indels.Chr05.vcf;
elif [ ${SLURM_ARRAY_TASK_ID} -eq 7 ]
then
 time -p gatk GenomicsDBImport -V $INPUT --genomicsdb-workspace-path $gVCF/Chr06 --intervals Chr06 --reader-threads 32;
 time -p gatk GenotypeGVCFs --variant gendb://$gVCF/Chr06 -R $REF --output $Chr/Chr06.vcf.gz -L Chr06;
 time -p gatk SelectVariants --variant $Chr/Chr06.vcf.gz --reference $REF -select-type SNP --output $SNPs/raw_snps.Chr06.vcf;
 time -p gatk VariantFiltration --variant $SNPs/raw_snps.Chr06.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name snp_filter --output $SNPs/hard_filtered_snps.Chr06.vcf;
 time -p gatk SelectVariants --variant $Chr/Chr06.vcf.gz --reference $REF -select-type INDEL --output $INDELs/raw_indels.Chr06.vcf;
 time -p gatk VariantFiltration --variant $INDELs/raw_indels.Chr06.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name indel_filter --output $INDELs/hard_filtered_indels.Chr06.vcf;
elif [ ${SLURM_ARRAY_TASK_ID} -eq 8 ]
then
 time -p gatk GenomicsDBImport -V $INPUT --genomicsdb-workspace-path $gVCF/Chr07 --intervals Chr07 --reader-threads 32;
 time -p gatk GenotypeGVCFs --variant gendb://$gVCF/Chr07 -R $REF --output $Chr/Chr07.vcf.gz -L Chr07;
 time -p gatk SelectVariants --variant $Chr/Chr07.vcf.gz --reference $REF -select-type SNP --output $SNPs/raw_snps.Chr07.vcf;
 time -p gatk VariantFiltration --variant $SNPs/raw_snps.Chr07.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name snp_filter --output $SNPs/hard_filtered_snps.Chr07.vcf;
 time -p gatk SelectVariants --variant $Chr/Chr07.vcf.gz --reference $REF -select-type INDEL --output $INDELs/raw_indels.Chr07.vcf;
 time -p gatk VariantFiltration --variant $INDELs/raw_indels.Chr07.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name indel_filter --output $INDELs/hard_filtered_indels.Chr07.vcf;
elif [ ${SLURM_ARRAY_TASK_ID} -eq 9 ]
then
 time -p gatk GenomicsDBImport -V $INPUT --genomicsdb-workspace-path $gVCF/Chr08 --intervals Chr08 --reader-threads 32;
 time -p gatk GenotypeGVCFs --variant gendb://$gVCF/Chr08 -R $REF --output $Chr/Chr08.vcf.gz -L Chr08;
 time -p gatk SelectVariants --variant $Chr/Chr08.vcf.gz --reference $REF -select-type SNP --output $SNPs/raw_snps.Chr08.vcf;
 time -p gatk VariantFiltration --variant $SNPs/raw_snps.Chr08.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name snp_filter --output $SNPs/hard_filtered_snps.Chr08.vcf;
 time -p gatk SelectVariants --variant $Chr/Chr08.vcf.gz --reference $REF -select-type INDEL --output $INDELs/raw_indels.Chr08.vcf;
 time -p gatk VariantFiltration --variant $INDELs/raw_indels.Chr08.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name indel_filter --output $INDELs/hard_filtered_indels.Chr08.vcf;
elif [ ${SLURM_ARRAY_TASK_ID} -eq 10 ]
then
 time -p gatk GenomicsDBImport -V $INPUT --genomicsdb-workspace-path $gVCF/Chr09 --intervals Chr09 --reader-threads 32;
 time -p gatk GenotypeGVCFs --variant gendb://$gVCF/Chr09 -R $REF --output $Chr/Chr09.vcf.gz -L Chr09;
 time -p gatk SelectVariants --variant $Chr/Chr09.vcf.gz --reference $REF -select-type SNP --output $SNPs/raw_snps.Chr09.vcf;
 time -p gatk VariantFiltration --variant $SNPs/raw_snps.Chr09.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name snp_filter --output $SNPs/hard_filtered_snps.Chr09.vcf;
 time -p gatk SelectVariants --variant $Chr/Chr09.vcf.gz --reference $REF -select-type INDEL --output $INDELs/raw_indels.Chr09.vcf;
 time -p gatk VariantFiltration --variant $INDELs/raw_indels.Chr09.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name indel_filter --output $INDELs/hard_filtered_indels.Chr09.vcf;
elif [ ${SLURM_ARRAY_TASK_ID} -eq 11 ]
then
 time -p gatk GenomicsDBImport -V $INPUT --genomicsdb-workspace-path $gVCF/Chr10 --intervals Chr10 --reader-threads 32;
 time -p gatk GenotypeGVCFs --variant gendb://$gVCF/Chr10 -R $REF --output $Chr/Chr10.vcf.gz -L Chr10;
 time -p gatk SelectVariants --variant $Chr/Chr10.vcf.gz --reference $REF -select-type SNP --output $SNPs/raw_snps.Chr10.vcf;
 time -p gatk VariantFiltration --variant $SNPs/raw_snps.Chr10.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name snp_filter --output $SNPs/hard_filtered_snps.Chr10.vcf;
 time -p gatk SelectVariants --variant $Chr/Chr10.vcf.gz --reference $REF -select-type INDEL --output $INDELs/raw_indels.Chr10.vcf;
 time -p gatk VariantFiltration --variant $INDELs/raw_indels.Chr10.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name indel_filter --output $INDELs/hard_filtered_indels.Chr10.vcf;
elif [ ${SLURM_ARRAY_TASK_ID} -eq 12 ]
then
 time -p gatk GenomicsDBImport -V $INPUT --genomicsdb-workspace-path $gVCF/Chr11 --intervals Chr11 --reader-threads 32;
 time -p gatk GenotypeGVCFs --variant gendb://$gVCF/Chr11 -R $REF --output $Chr/Chr11.vcf.gz -L Chr11;
 time -p gatk SelectVariants --variant $Chr/Chr11.vcf.gz --reference $REF -select-type SNP --output $SNPs/raw_snps.Chr11.vcf;
 time -p gatk VariantFiltration --variant $SNPs/raw_snps.Chr11.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name snp_filter --output $SNPs/hard_filtered_snps.Chr11.vcf;
 time -p gatk SelectVariants --variant $Chr/Chr11.vcf.gz --reference $REF -select-type INDEL --output $INDELs/raw_indels.Chr11.vcf;
 time -p gatk VariantFiltration --variant $INDELs/raw_indels.Chr11.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name indel_filter --output $INDELs/hard_filtered_indels.Chr11.vcf;
elif [ ${SLURM_ARRAY_TASK_ID} -eq 13 ]
then
 time -p gatk GenomicsDBImport -V $INPUT --genomicsdb-workspace-path $gVCF/Chr12 --intervals Chr12 --reader-threads 32;
 time -p gatk GenotypeGVCFs --variant gendb://$gVCF/Chr12 -R $REF --output $Chr/Chr12.vcf.gz -L Chr12;
 time -p gatk SelectVariants --variant $Chr/Chr12.vcf.gz --reference $REF -select-type SNP --output $SNPs/raw_snps.Chr12.vcf;
 time -p gatk VariantFiltration --variant $SNPs/raw_snps.Chr12.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name snp_filter --output $SNPs/hard_filtered_snps.Chr12.vcf;
 time -p gatk SelectVariants --variant $Chr/Chr12.vcf.gz --reference $REF -select-type INDEL --output $INDELs/raw_indels.Chr12.vcf;
 time -p gatk VariantFiltration --variant $INDELs/raw_indels.Chr12.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name indel_filter --output $INDELs/hard_filtered_indels.Chr12.vcf;
elif [ ${SLURM_ARRAY_TASK_ID} -eq 14 ]
then
 time -p gatk GenomicsDBImport -V $INPUT --genomicsdb-workspace-path $gVCF/Chr13 --intervals Chr13 --reader-threads 32;
 time -p gatk GenotypeGVCFs --variant gendb://$gVCF/Chr13 -R $REF --output $Chr/Chr13.vcf.gz -L Chr13;
 time -p gatk SelectVariants --variant $Chr/Chr13.vcf.gz --reference $REF -select-type SNP --output $SNPs/raw_snps.Chr13.vcf;
 time -p gatk VariantFiltration --variant $SNPs/raw_snps.Chr13.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name snp_filter --output $SNPs/hard_filtered_snps.Chr13.vcf;
 time -p gatk SelectVariants --variant $Chr/Chr13.vcf.gz --reference $REF -select-type INDEL --output $INDELs/raw_indels.Chr13.vcf;
 time -p gatk VariantFiltration --variant $INDELs/raw_indels.Chr13.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name indel_filter --output $INDELs/hard_filtered_indels.Chr13.vcf;
elif [ ${SLURM_ARRAY_TASK_ID} -eq 15 ]
then
 time -p gatk GenomicsDBImport -V $INPUT --genomicsdb-workspace-path $gVCF/Chr14 --intervals Chr14 --reader-threads 32;
 time -p gatk GenotypeGVCFs --variant gendb://$gVCF/Chr14 -R $REF --output $Chr/Chr14.vcf.gz -L Chr14;
 time -p gatk SelectVariants --variant $Chr/Chr14.vcf.gz --reference $REF -select-type SNP --output $SNPs/raw_snps.Chr14.vcf;
 time -p gatk VariantFiltration --variant $SNPs/raw_snps.Chr14.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name snp_filter --output $SNPs/hard_filtered_snps.Chr14.vcf;
 time -p gatk SelectVariants --variant $Chr/Chr14.vcf.gz --reference $REF -select-type INDEL --output $INDELs/raw_indels.Chr14.vcf;
 time -p gatk VariantFiltration --variant $INDELs/raw_indels.Chr14.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name indel_filter --output $INDELs/hard_filtered_indels.Chr14.vcf;
elif [ ${SLURM_ARRAY_TASK_ID} -eq 16 ]
then
 time -p gatk GenomicsDBImport -V $INPUT --genomicsdb-workspace-path $gVCF/Chr15 --intervals Chr15 --reader-threads 32;
 time -p gatk GenotypeGVCFs --variant gendb://$gVCF/Chr15 -R $REF --output $Chr/Chr15.vcf.gz -L Chr15;
 time -p gatk SelectVariants --variant $Chr/Chr15.vcf.gz --reference $REF -select-type SNP --output $SNPs/raw_snps.Chr15.vcf;
 time -p gatk VariantFiltration --variant $SNPs/raw_snps.Chr15.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name snp_filter --output $SNPs/hard_filtered_snps.Chr15.vcf;
 time -p gatk SelectVariants --variant $Chr/Chr15.vcf.gz --reference $REF -select-type INDEL --output $INDELs/raw_indels.Chr15.vcf;
 time -p gatk VariantFiltration --variant $INDELs/raw_indels.Chr15.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name indel_filter --output $INDELs/hard_filtered_indels.Chr15.vcf;
elif [ ${SLURM_ARRAY_TASK_ID} -eq 17 ]
then
 time -p gatk GenomicsDBImport -V $INPUT --genomicsdb-workspace-path $gVCF/Chr16 --intervals Chr16 --reader-threads 32;
 time -p gatk GenotypeGVCFs --variant gendb://$gVCF/Chr16 -R $REF --output $Chr/Chr16.vcf.gz -L Chr16;
 time -p gatk SelectVariants --variant $Chr/Chr16.vcf.gz --reference $REF -select-type SNP --output $SNPs/raw_snps.Chr16.vcf;
 time -p gatk VariantFiltration --variant $SNPs/raw_snps.Chr16.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name snp_filter --output $SNPs/hard_filtered_snps.Chr16.vcf;
 time -p gatk SelectVariants --variant $Chr/Chr16.vcf.gz --reference $REF -select-type INDEL --output $INDELs/raw_indels.Chr16.vcf;
 time -p gatk VariantFiltration --variant $INDELs/raw_indels.Chr16.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name indel_filter --output $INDELs/hard_filtered_indels.Chr16.vcf;
elif [ ${SLURM_ARRAY_TASK_ID} -eq 18 ]
then
 time -p gatk GenomicsDBImport -V $INPUT --genomicsdb-workspace-path $gVCF/Chr17 --intervals Chr17 --reader-threads 32;
 time -p gatk GenotypeGVCFs --variant gendb://$gVCF/Chr17 -R $REF --output $Chr/Chr17.vcf.gz -L Chr17;
 time -p gatk SelectVariants --variant $Chr/Chr17.vcf.gz --reference $REF -select-type SNP --output $SNPs/raw_snps.Chr17.vcf;
 time -p gatk VariantFiltration --variant $SNPs/raw_snps.Chr17.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name snp_filter --output $SNPs/hard_filtered_snps.Chr17.vcf;
 time -p gatk SelectVariants --variant $Chr/Chr17.vcf.gz --reference $REF -select-type INDEL --output $INDELs/raw_indels.Chr17.vcf;
 time -p gatk VariantFiltration --variant $INDELs/raw_indels.Chr17.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name indel_filter --output $INDELs/hard_filtered_indels.Chr17.vcf;
else
 time -p gatk GenomicsDBImport -V $INPUT --genomicsdb-workspace-path $gVCF/Chr18 --intervals Chr18 --reader-threads 32;
 time -p gatk GenotypeGVCFs --variant gendb://$gVCF/Chr18 -R $REF --output $Chr/Chr18.vcf.gz -L Chr18;
 time -p gatk SelectVariants --variant $Chr/Chr18.vcf.gz --reference $REF -select-type SNP --output $SNPs/raw_snps.Chr18.vcf;
 time -p gatk VariantFiltration --variant $SNPs/raw_snps.Chr18.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filter-name snp_filter --output $SNPs/hard_filtered_snps.Chr18.vcf;
 time -p gatk SelectVariants --variant $Chr/Chr18.vcf.gz --reference $REF -select-type INDEL --output $INDELs/raw_indels.Chr18.vcf;
 time -p gatk VariantFiltration --variant $INDELs/raw_indels.Chr18.vcf --reference $REF --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name indel_filter --output $INDELs/hard_filtered_indels.Chr18.vcf;
fi
