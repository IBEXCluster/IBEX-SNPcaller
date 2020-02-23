

# Principal investigators (PI)
**Prof. Mark A. Tester**, Associate Director, Center for Desert Agriculture <br/>
**Prof. Simon g. Krattinger**, Prof. of Plant Science <br/>
Center for Desert Agriculture, <br/>
4700 King Abdullah University of Science and Technology <br/>
Thuwal 23955-6900 <br/>
Kingdom of Saudi Arabia <br/>


# Authors:
Nagarajan Kathiresan <nagarajan.kathiresan@kaust.edu.sa> <br/>
Michael Abrouk <michael.abrouk@kaust.edu.sa> <br/>
Elodie L. Rey <elodie.rey@kaust.edu.sa> <br/>


# About Ibex cluster

Ibex is a heterogeneous group of nodes, a mix of AMD, Intel and Nvidia GPUs with different architectures that gives the users a variety of options to work on. Overall, Ibex is made up of 488+ nodes togeter has a heterogeneous cluster and the workload is managed by the SLURM scheduler. More information is available in https://www.hpc.kaust.edu.sa/ibex <br/>

Operating System on nodes: CentOS 7.6 <br/>
Scheduler : SLURM version 19.05.2 <br/>



# Ibex-SNPCaller pipeline

The objective of this SNP Caller pipeline is to automate and optimize the various job steps across multiple samples. To simplify the pipeline for various project requirements, we separated the pipeline into two parts: (i) Data processing and (2) Downstream analysis. <br/> 

## 1. Data processing 

We followed different steps for genome data processing (as part of best practices pipeline) that includes (a) Read trimming (b) Read mappring (c) Mark Duplicate (d) Add/Replace read groups (e) HaplotypeCalling and (f) Compress & Index the gVCF files. All these steps in the data processing pipeline are automated based on the job dependency conditions from SLURM workload scheduler and the automated scripts will accept all the samples from the given INPUT file directory. Further, the software and/or the job steps can be modified based on the various requirements of the project. We selected the optimal number of cores for each job steps based on our vaious case studies. This automated data processing script called "workflow.sh" is available for your experiments and the pipeline stages are demonstrated in Figure (a) Pipeline steps in Data processing.     

**List of software** <br/>
trimmomatic version 0.38 <br/>
bwa version 0.7.17  <br/>
samtools version 1.8 <br/>
gatk version 4.0.1.1 <br/>
tabix version 0.2.6 <br/>

![](https://www.hpc.kaust.edu.sa/sites/default/files/files/public/workflows/HaplotypeCaller_workflow.png)

<p align="center"> Figure (a) Pipeline steps in Data processing </p>

## 2. Downstream analysis
During this downstream analysis, the SNPs & INDELs are separated from the above data processing pipeline. This downstream analysis become more complex and highly Compute-intensive operation when multiple samples are merged into a single gVCF file called "Combined gVCF". Most of the population genetics project required this Combine gVCF file because, the SNPs and INDELs are separated based on chromosome intervals. The latest GATK version 4.x has the featues for importing the genomic data from multiple samples based on the chromosome interval list. Additionally, the SNPs and INDELs are filtered based on various filtering expressions called "hard filters". This automated pipeline called "downstream_analysis.sh" can be modified based on the project requirements and/or the chromosomes are separated based on the reference files. More details are described in Figure (b) Pipeline steps in downstream analysis.      

**List of software** <br/>
Java 8.131 <br/>
GATK 4.1.1.0 <br/>

![](https://www.hpc.kaust.edu.sa/sites/default/files/files/public/workflows/Downstream.png)
<p align="center"> Figure (b) Pipeline steps in downstream analysis <p/>

## Add-on features
To simplify the job monitoring on multiple samples, we calssified the job statues into 7 categories: (1) Completed (2) Failed (3) Out of Memory (4) Timeout (5) Cancelled (6) Pending and (7) Running. This job classifications are useful (as part of SLURM job management) to handle larger number of samples across multiple job steps. <br/>

