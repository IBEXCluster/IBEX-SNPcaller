#!/bin/bash
#######################################################################################################################################################
#  Trimmomatic + GATK 4 (Picard functionality) + GATK 3.8 HaplotypeCaller  + bgzip ( the *.g.VCF files) & tabix-index based workflow 
#     Version 5.0 dated 22 Oct 2019
#     Modification required for any users:
#      1. INPUT (directory location)
#      2. PROJECT (directory location)
#
########################################################################################################################################################
## Software Modules 

function Software_Modules()
{
module load trimmomatic/0.38 bwa/0.7.17/gnu-6.4.0 samtools/1.8 gatk/4.0.1.1 tabix/0.2.6
export GATK="/ibex/scratch/kathirn/work/project/for_elodie/ALL/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"

}

### Execute the workflow 
function Workflow()
{
printf "\n Provide complete path of your genome reference file (e.g.: /ibex/reference/KSL/human_ref/human_g1k_v37.fasta) \n" ;
read REF;
if [ ! -f $REF ]; then
  printf " Your Genome reference file does not exist \n" ;
  break;
fi;

printf "\n Provide your FASTQ files directory (absulute path, e.g.: /ibex/scratch/kathirn/my_samples/) \n" ;
read INPUT;
if [ ! -d $INPUT ];then
 printf " Your INPUT [ $INPUT ] directory does not exist \n" ;
 break;
else
 MY_FIRST_FILE=`ls -lrta $INPUT/*1.fq.gz | head -1 | awk '{print $9}'` ;
fi

if [ ! -f $MY_FIRST_FILE ]; then
 printf " Their is no FASTQ (*.fq) samples in the INPUT [ $INPUT ] directory \n";
 break;
else
 export TOTAL_SAMPLE=`ls -lrta $INPUT/*1.fq.gz | wc -l` ;
fi 

printf "\n Provide your PROJECT directory, where all the data processing (BAM, VCF, SLURM-LOGS, etc.) will happen (e.g.: /ibex/scratch/kathirn/1000genome/ ) \n" 
; 
read PROJECT;

Software_Modules ;
  printf "\033c"
  printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
  printf "\n Genome Reference file: $REF" ;
  printf "\n Total number of samples: ${TOTAL_SAMPLE}" ;
  printf "\n BWA + GATK workflow will be executed now \n" ;
   printf "\n Workflow starts ..... \n" ;
  printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
#  echo " Workflow start "
#  break ;

 

## Directory Variables 

export BAM=${PROJECT}/BAM;
export VCF=${PROJECT}/VCF;
export LOGS=${PROJECT}/LOGS;
export DEBUG=${PROJECT}/DEBUG;


mkdir -p $BAM;
mkdir -p $VCF;
mkdir -p $LOGS;
mkdir -p $DEBUG; 

rm -Rf $DEBUG/job_status.txt ;
touch $DEBUG/job_status.txt ;



## For Sample count 
set COUNT=0;

## Workflow steps 
for SAMPLE in `ls $INPUT/*1.fq.gz`;
do 
  PREFIX=`basename $SAMPLE _1.fq.gz` ;
  LOCATION=${SAMPLE%/*};
#  echo "$PREFIX and $LOCATION" 

 #### Step 1. trimming of reads
    MEM="32gb"
    CORES=4
    JOB1_NAME="Trimming"
    JOB1_TYPE="sbatch --partition=batch --job-name=${JOB1_NAME}.${PREFIX} --time=2:00:00 --output=$LOGS/${JOB1_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB1_NAME}.
${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB1_CMD="time -p java -XX:+UseParallelGC -XX:ParallelGCThreads=${CORES} -jar $TRIMMOMATIC_JAR PE -phred33 $LOCATION/$PREFIX-P1.fq $LOCATION/$PREFIX-P2.fq $B
AM/$PREFIX.trimmed.P1.fastq $BAM/$PREFIX.up.1.fast $BAM/$PREFIX.trimmed.P2.fastq $BAM/$PREFIX.up.2.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:50" ;
    JOB1_ID=$(${JOB1_TYPE} --parsable --wrap="${JOB1_CMD}");
    echo "$PREFIX sample with the job id=$JOB1_ID and Job Name=$JOB1_NAME submitted"
    printf "%-40d \t %-40s \t %-40s\n" "$JOB1_ID" "$JOB1_NAME" "$PREFIX" >> $DEBUG/job_status.txt
   
 #### Step 2. BWA MEM
    MEM="115gb"
    CORES=16
    JOB2_NAME="bwa-mem"
    JOB2_TYPE="sbatch --partition=batch --job-name=${JOB2_NAME}.${PREFIX} --time=12:00:00 --output=$LOGS/${JOB2_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB2_NAME}
.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB2_CMD="time -p bwa mem -M -k 30 -t $CORES $REF $BAM/$PREFIX.trimmed.P1.fastq $BAM/$PREFIX.trimmed.P2.fastq | samtools view -@ $CORES -b -S -h -q 30 - | sa
mtools sort - > $BAM/$PREFIX.sorted.bam"
    JOB2_ID=$(${JOB2_TYPE} --parsable --dependency=afterok:${JOB1_ID} --wrap="${JOB2_CMD}");
    echo "$PREFIX sample with the job id=$JOB2_ID and Job Name=$JOB2_NAME submitted"   
    printf "%-40d \t %-40s \t %-40s\n" "$JOB2_ID" "$JOB2_NAME" "$PREFIX" >> $DEBUG/job_status.txt

 #### Step 3. MarkDuplicates 
    MEM="64gb"
    CORES=1
    JOB3_NAME="MarkDuplicate"
    JOB3_TYPE="sbatch --partition=batch --job-name=${JOB3_NAME}.${PREFIX} --time=4:00:00 --output=$LOGS/${JOB3_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB3_NAME}.
${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB3_CMD="time -p gatk MarkDuplicates --INPUT=$BAM/$PREFIX.sorted.bam --METRICS_FILE=$BAM/$PREFIX.metrics.txt --OUTPUT=$BAM/$PREFIX.rmdup.bam"
    JOB3_ID=$(${JOB3_TYPE} --parsable --dependency=afterok:${JOB2_ID} --wrap="${JOB3_CMD}");
    echo "$PREFIX sample with the job id=$JOB3_ID and Job Name=$JOB3_NAME submitted"   
    printf "%-40d \t %-40s \t %-40s\n" "$JOB3_ID" "$JOB3_NAME" "$PREFIX" >> $DEBUG/job_status.txt

 ##### 4. AddOrReplace
    ## Note: VALIDATION_STRINGENCY=LENIENT missing in GATK 4.0 
    MEM="32gb"
    CORES=1
    JOB4_NAME="AddOrReplace"
    JOB4_TYPE="sbatch --partition=batch --job-name=${JOB4_NAME}.${PREFIX} --time=3:00:00 --output=$LOGS/${JOB4_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB4_NAME}.
${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB4_CMD="time -p gatk AddOrReplaceReadGroups --INPUT=$BAM/$PREFIX.rmdup.bam --OUTPUT=$BAM/$PREFIX.rgroup.bam --SORT_ORDER=coordinate --RGSM=$PREFIX --RGPU=n
one --RGID=1 --RGLB=lib1 --RGPL=Illumina"
    JOB4_ID=$(${JOB4_TYPE} --parsable --dependency=afterok:${JOB3_ID} --wrap="${JOB4_CMD}");
    echo "$PREFIX sample with the job id=$JOB4_ID and Job Name=$JOB4_NAME submitted"   
    printf "%-40d \t %-40s \t %-40s\n" "$JOB4_ID" "$JOB4_NAME" "$PREFIX" >> $DEBUG/job_status.txt
    
 ##### 5. Samtools Index
    MEM="32gb"
    CORES=1
    JOB5_NAME="Samtool-Index"
    JOB5_TYPE="sbatch --partition=batch --job-name=${JOB5_NAME}.${PREFIX} --time=90:00 --output=$LOGS/${JOB5_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB5_NAME}.${
PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB5_CMD="time -p samtools index $BAM/$PREFIX.rgroup.bam"
    JOB5_ID=$(${JOB5_TYPE} --parsable --dependency=afterok:${JOB4_ID} --wrap="${JOB5_CMD}");
    echo "$PREFIX sample with the job id=$JOB5_ID and Job Name=$JOB5_NAME submitted"   
    printf "%-40d \t %-40s \t %-40s\n" "$JOB5_ID" "$JOB5_NAME" "$PREFIX" >> $DEBUG/job_status.txt
 
 ##### 6. HaplotypeCaller
    MEM="115gb"
    CORES=16
    JOB6_NAME="HaplotypeCaller"
    JOB6_TYPE="sbatch --partition=batch --job-name=${JOB6_NAME}.${PREFIX} --time=30:00:00 --output=$LOGS/${JOB6_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB6_NAME}
.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB6_CMD="time -p java -XX:+UseParallelGC -XX:ParallelGCThreads=${CORES} -jar $GATK -T HaplotypeCaller -nct $CORES -pairHMM VECTOR_LOGLESS_CACHING -R $REF -I
 $BAM/$PREFIX.rgroup.bam --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o $VCF/$PREFIX.snps.indels.g.vcf" ;
    JOB6_ID=$(${JOB6_TYPE} --parsable --dependency=afterok:${JOB5_ID} --wrap="${JOB6_CMD}");
    echo "$PREFIX sample with the job id=$JOB6_ID and Job Name=$JOB6_NAME submitted"   
    printf "%-40d \t %-40s \t %-40s\n" "$JOB6_ID" "$JOB6_NAME" "$PREFIX" >> $DEBUG/job_status.txt

 ##### 7. Compress the g.VCF file using bgzip
    MEM="32gb"
    CORES=1
    JOB7_NAME="bgzip"
    JOB7_TYPE="sbatch --partition=batch --job-name=${JOB7_NAME}.${PREFIX} --time=2:00:00 --output=$LOGS/${JOB7_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB7_NAME}.
${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB7_CMD="time -p bgzip -o $VCF/$PREFIX.snps.indels.g.vcf" ;
    JOB7_ID=$(${JOB7_TYPE} --parsable --dependency=afterok:${JOB6_ID} --wrap="${JOB7_CMD}");
    echo "$PREFIX sample with the job id=$JOB7_ID and Job Name=$JOB7_NAME submitted"   
    printf "%-40d \t %-40s \t %-40s\n" "$JOB7_ID" "$JOB7_NAME" "$PREFIX" >> $DEBUG/job_status.txt

 ##### 8. Create Tabix-Index for the g.VCF.gz file
    MEM="32gb"
    CORES=1
    JOB8_NAME="tabix"
    JOB8_TYPE="sbatch --partition=batch --job-name=${JOB8_NAME}.${PREFIX} --time=30:00 --output=$LOGS/${JOB8_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB8_NAME}.${
PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB8_CMD="time -p tabix -o $VCF/$PREFIX.snps.indels.g.vcf.gz" ;
    JOB8_ID=$(${JOB8_TYPE} --parsable --dependency=afterok:${JOB7_ID} --wrap="${JOB8_CMD}");
    echo "$PREFIX sample with the job id=$JOB8_ID and Job Name=$JOB8_NAME submitted"   
    printf "%-40d \t %-40s \t %-40s\n" "$JOB8_ID" "$JOB8_NAME" "$PREFIX" >> $DEBUG/job_status.txt

 
 COUNT=$((COUNT + 1))
done
 
echo "========================================"
echo "Total Number of samples submitted: ${COUNT} in 8 Steps"
echo "========================================" 

}

function main()
{
 menu ;
}
function menu()
{
PS3='Please enter your choice: '
options=("Verify the Reference and Index files" "Execute the Workflow" "List the status of success or failed jobs" "Quit")
select opt in "${options[@]}"
do
    case $opt in
        "Verify the Reference and Index files")
            echo "Will verify your Reference"
            reference ;
	    break ;
            ;;
        "Execute the Workflow")
            echo "Will prepare to execute BWA + GATK workflow"
 	    Workflow ;
            break ;
            ;;
        "List the status of success or failed jobs")
            echo "Will monitor your job stastics from your DEBUG directory"
	    Job_Monitor ;
            break ;
            ;;
        "Quit")
            break
            ;;
        *) echo "invalid option $REPLY";;
    esac
done

} 

function reference()
{
 MEM="115gb" ;
 CORES=1;
 printf "\n\nPlease enter your full path of reference directory name"
 printf "\n \t Example: /ibex/reference/KSL/human_ref/human_g1k_v37.fasta \n\n"
 read REF;
 DIR=`echo "${REF%/*}"`
 REF_FILE=`basename "$REF"` ;
 REF_NAME=${REF_FILE%.*} ;
 echo ${REF_FILE} "and" $DIR "and" ${REF_NAME}; 

 printf "\n\n Validating the reference: ${REF_NAME} ... in the Directory: $DIR \n"

## Check BWA Indexing required or may not

 
 printf "\n ***** CHECKING BWA INDEX FILES ...... **** \n";

 AMB="$REF.amb" ; 
 ANN="$REF.ann"
 BWT="$REF.bwt"
 PAC="$REF.pac"
 SA="$REF.sa"
 if [ -f $AMB ] && [ -f $ANN ] && [ -f $BWT ] && [ -f $PAC ] && [ -f $SA ]; then
  printf "\n\t [ $AMB ] is the text file, to record appearance of N (or other non-ATGC) in the ref fasta .... EXIST[\xE2\x9C\x94]" ;
  printf "\n\t [ $ANN ] is text file, to record ref sequences, name, length, etc. ... EXIST[\xE2\x9C\x94]" ;
  printf "\n\t [ $BWT ] is binary, the Burrows-Wheeler transformed sequence. ... EXIST[\xE2\x9C\x94]" ;
  printf "\n\t [ $PAC ] is binary, packaged sequence (four base pairs encode one byte).... EXIST[\xE2\x9C\x94]" ;
  printf "\n\t [ $SA  ] is binary, suffix array index..... EXIST[\xE2\x9C\x94]" ;
  printf "\n\t ****** All reference index files are exist and BWA INDEX NOT REQUIRED ******* \n\n" ;
 else
  echo "BWA Indexing required" ;
  module load bwa/0.7.17/gnu-6.4.0 ;

  JOB_NAME="BWA_Index"
  JOB_TYPE="sbatch --partition=batch --job-name=${JOB_NAME}.${REF_NAME} --time=10:00:00 --output=$DIR/${JOB_NAME}.${REF_NAME}.%J.out --error=$DIR/${JOB_NAME}.${R
EF_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
  JOB_CMD="time -p bwa index $REF" 
  echo " \n\t BWA INDEX Job $JOB_CMD will be submitted as a batch job ...... \n" 
  JOB_ID=$(${JOB_TYPE} --parsable --wrap="${JOB_CMD}");
  echo "$REF_NAME reference is indexed using BWA MEM, the job id=$JOB_ID and Job Name=$JOB_NAME submitted"   
 fi

## Check SAMTOOLS Index 

 printf "\n ***** CHECKING SAMTOOLS INDEX FILES ...... **** \n";

 FAI="$REF.fai";
 if [ -f $FAI ]; then
  printf "\n\t [ $FAI ] is the text file and Samtools index enabling random access to FASTA and FASTQ files .... EXIST[\xE2\x9C\x94]" ;
  printf "\n\t ****** All reference index files are exist and SAMTOOLS INDEX NOT REQUIRED ******* \n\n" ;
 else 
 module load samtools/1.8
 JOB_NAME="Samtools_Index"
 JOB_TYPE="sbatch --partition=batch --job-name=${JOB_NAME}.${REF_NAME} --time=4:00:00 --output=$DIR/${JOB_NAME}.${REF_NAME}.%J.out --error=$DIR/${JOB_NAME}.${REF
_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
 JOB_CMD="time -p samtools faidx $REF" 
 echo " \n\t SAMTOOLS INDEX Job $JOB_CMD will be submitted as a batch job ...... \n" 
 JOB_ID=$(${JOB_TYPE} --parsable --wrap="${JOB_CMD}");
 echo "$REF_NAME reference is indexed using SAMTOOLS, the job id=$JOB_ID and Job Name=$JOB_NAME submitted"   
fi

## Check Picard Index 

 printf "\n\n \n ***** CHECKING PICARD (GATK) SEQUENCE DIRECTORY FILES ...... **** \n";
  
 DICT="$REF.dict"
  if [ -f $DICT ]; then
    printf "\n\t [ $FAI ] is the text file and Picard (GATK) Sequence Dictionary file  .... EXIST[\xE2\x9C\x94]" ;
    printf "\n\t ****** All reference index files are exist and PICARD (GATK) SEQUENCE DIRECTORY NOT REQUIRED ******* \n\n" ;
 else
 module load gatk/4.0.1.1
 JOB_NAME="Picard_Index"
 JOB_TYPE="sbatch --partition=batch --job-name=${JOB_NAME}.${REF_NAME} --time=30:00 --output=$DIR/${JOB_NAME}.${REF_NAME}.%J.out --error=$DIR/${JOB_NAME}.${REF_N
AME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
 JOB_CMD="time -p gatk CreateSequenceDictionary --REFERENCE $REF --OUTPUT $REF.dict"
 printf " \n\t PICARD SEQUENCE DIRECTORY Job $JOB_CMD will be submitted as a batch job ...... \n" 
 JOB_ID=$(${JOB_TYPE} --parsable --wrap="${JOB_CMD}");
 echo "$REF_NAME reference is indexed using PICARD (GATK), the job id=$JOB_ID and Job Name=$JOB_NAME submitted"   
 fi


}

function Job_Monitor()
{

printf "\n Please provide your PROJECT directory location \n" ;
read PROJ;

if [ ! -d $PROJ ]; then
  printf "\n PROJECT [ $PROJ ] directory does not exist \n" ;
else

while read line
do
    COLUMN=($line);
    STATUS=`sacct -j ${COLUMN[0]} --format=state%13 | sed -n '3 p' | xargs | sed 's/\s.*$//' `;

   case $STATUS in
      COMPLETED )
        	printf "\n ${COLUMN[2]} successfully completed " 
		;;
      FAILED )
		printf "\n ******** ${COLUMN[2]} Failed ********" 
		;;
      OUT_OF_MEMORY )
		printf "\n ##### ${COLUMN[2]} require more memory ####" 
		;;
      TIMEOUT	)
		printf "\n ^^^^^^^^ ${COLUMN[2]} require more itun time, please increase TIME in your job script file ^^^^^^^^" 
		;;
      CANCELLED	)
		printf "\n -C-C-C-C-C- ${COLUMN[2]} job was cancelled by the user -C-C-C-C-C-"
		;;
      PENDING )
		printf "\n ${COLUMN[2]} job waiting in the queue"
		;;
      RUNNING )
		printf "\n ${COLUMN[2]} job running in Ibex"
		;;
       * )
		printf "\n ${COLUMN[2]} not successful, please check with ibex team"
		;;
   esac
done < $PROJ/DEBUG/job_status.txt ;
printf "\n\n" ;
fi 

}

main 
