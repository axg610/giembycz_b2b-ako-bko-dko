#!/bin/bash

#SBATCH --job-name=k
#SBATCH --chdir="/work/giembycz_lab/tamkeenko_analysis/out/"
#SBATCH --error="%j_fastqc.err"
#SBATCH --output="%j_fastqc.out"
#SBATCH --time=00-05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40GB
#SBATCH --partition=cpu2023,cpu2022

echo ==============================================
echo Hostname: `hostname` 
echo Job ID: $SLURM_JOBID
echo Node: $SLURM_NODELIST
echo Sample: $1
echo Task: FastQC Kallisto

start=`date +%Y-%m-%d\ %H:%M:%S`

module load openmpi/4.1.1-gnu
module load python/anaconda3-2018.12
source activate rnaseq

cd /work/giembycz_lab/tamkeenko_analysis/fastq/$1

mkdir -p /work/giembycz_lab/tamkeenko_analysis/test_kallisto2/$1
# fastqc -q -o /work/giembycz_lab/tamkeenko_analysis/test_qc/$1 *_R1*.fastq.gz *_R2*.fastq.gz

kallisto quant -i /work/giembycz_lab/ALI_analysis/kallisto_index_cdna_GRCh38_p13 -b 100 -t 10 *L001_R1*.fastq.gz *L001_R2*.fastq.gz *L002_R1*.fastq.gz *L002_R2*.fastq.gz -o /work/giembycz_lab/tamkeenko_analysis/test_kallisto2/$1


echo ==============================================
echo Start: $start
echo End: $end
