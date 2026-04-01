ssh alex.gao1@arc.ucalgary.ca

## ARC SETUP
```bash ARC SETUP
salloc -c 4 --mem 32GB --time 05:00:00

MULTIQC="/home/alex.gao1/tools/multiqc_latest.sif"
apptainer exec --bind /work:/work "$MULTIQC" multiqc --version  #v1.19

module load kallisto
kallisto --version  #0.46.1

module load R   #4.4.1

mkdir -p /work/giembycz_lab/B2B_dKO_analysis
cd /work/giembycz_lab/B2B_dKO_analysis
```

## MOVE FILES
```bash MOVE FILES
#files from CHGI
ls '/work/giembycz_lab/Job4182/fastq'
ls '/work/giembycz_lab/Job4182/fastqc'

#copy run files
for r1 in /work/giembycz_lab/Job4182/fastq/*_R1_001.fastq.gz; do
    sample=$(basename "$r1" | sed -E 's/^Li[0-9]+-([A-Za-z0-9]+)_S[0-9]+_R1_001\.fastq\.gz/\1/')
    mkdir -p "/work/giembycz_lab/B2B_dKO_analysis/fastq/$sample"

    echo "Copying $sample..."

    cp "${r1}" "/work/giembycz_lab/B2B_dKO_analysis/fastq/$sample/${sample}_R1.fastq.gz"
    cp "${r1/_R1_/_R2_}" "/work/giembycz_lab/B2B_dKO_analysis/fastq/$sample/${sample}_R2.fastq.gz"
done

#copy fastqc files and run multiQC
cp -r '/work/giembycz_lab/Job4182/fastqc' fastqc/
apptainer exec --bind /work:/work "$MULTIQC" multiqc fastqc/ -o multiqc_reads
```

## RUN KALLISTO
```bash RUN KALLISTO
#grab kallisto index
cp /work/giembycz_lab/ALI_analysis/kallisto_index_cdna_GRCh38_p13 .

#make kallisto run file
cat <<'EOF' > run_kallisto.sh
#!/bin/bash
#SBATCH --job-name=kallisto
#SBATCH --chdir=/work/giembycz_lab/B2B_dKO_analysis
#SBATCH --output=./logs/kallisto_%j.out
#SBATCH --error=./logs/kallisto_%j.err
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=8G

module load kallisto

SAMPLE="$1"

mkdir -p logs
mkdir -p kallisto/$SAMPLE

kallisto quant -i /work/giembycz_lab/ALI_analysis/kallisto_index_cdna_GRCh38_p13 \
    -b 100 -t 10 \
    -o kallisto/$SAMPLE \
    fastq/$SAMPLE/*_R1.fastq.gz fastq/$SAMPLE/*_R2.fastq.gz
EOF

#submit jobs
for s in fastq/*; do
  sbatch run_kallisto.sh $(basename "$s")
done

#run multiqc
apptainer exec --bind /work:/work "$MULTIQC" multiqc logs/ -o multiqc_kallisto
```

## RUN SLEUTH
```bash BOOT R
R
```

```R standard sleuth analysis
.libPaths("/home/alex.gao1/R")
setwd("/work/giembycz_lab/B2B_dKO_analysis")

library(dplyr)
library(readr)
library(sleuth)

#======================
# PREP INPUTS
#======================

t2g = read_tsv("Homo_sapiens.GRCh38.p13.cdna.all.070121.mart_export.txt") %>%
  select(target_id = 1, Gene = 2) %>%
  na.omit() %>% 
  distinct()
t2g <- t2g[!duplicated(t2g$target_id), ]

s2c = read_tsv("meta.txt") %>%
    mutate(
        treatment = factor(treatment, levels = c("NS", "Form")),
        condition = factor(condition, levels = c("WT", "DKO")))

s2c_wt = filter(s2c, condition == "WT")
s2c_dko = filter(s2c, condition == "DKO")

s2c_ns = filter(s2c, treatment == "NS")

new_filter <- function(row, min_reads = 5, min_prop = 0.2){mean(row >= min_reads) >= min_prop}

#======================
# CREATE SLEUTH OBJECTS
#======================

so = sleuth_prep(
    sample_to_covariates = s2c,
    target_mapping = t2g,
    gene_mode = TRUE,
    aggregation_column = "Gene",
    filter_fun = new_filter,
    num_cores = 1
) %>%
sleuth_fit(~condition*treatment, "full") %>%
sleuth_fit(~1, "reduced") %>%
sleuth_fit(~condition, "condition") %>%
sleuth_fit(~treatment, "treatment") %>%
sleuth_lrt("reduced", "full") %>%
sleuth_lrt("condition", "full") %>%
sleuth_lrt("treatment", "full") %>%
sleuth_wt("conditionDKO:treatmentForm")

saveRDS(so, "sleuth/so.rds")
# so = readRDS("sleuth/so.rds")

so_wt = sleuth_prep(
    sample_to_covariates = s2c_wt,
    target_mapping = t2g,
    gene_mode = TRUE,
    aggregation_column = "Gene",
    filter_fun = new_filter,
    num_cores = 1
) %>%
sleuth_fit(~treatment) %>%
sleuth_wt("treatmentForm")

saveRDS(so_wt, "sleuth/so_wt.rds")
# so_wt = readRDS("sleuth/so_wt.rds")

so_dko = sleuth_prep(
    sample_to_covariates = s2c_dko,
    target_mapping = t2g,
    gene_mode = TRUE,
    aggregation_column = "Gene",
    filter_fun = new_filter,
    num_cores = 1
) %>%
sleuth_fit(~treatment) %>%
sleuth_wt("treatmentForm")

saveRDS(so_dko, "sleuth/so_dko.rds")
# so_dko = readRDS("sleuth/so_dko.rds")

so_ns = sleuth_prep(
    sample_to_covariates = s2c_ns,
    target_mapping = t2g,
    gene_mode = TRUE,
    aggregation_column = "Gene",
    filter_fun = new_filter,
    num_cores = 1
) %>%
sleuth_fit(~condition) %>%
sleuth_wt("conditionDKO")

saveRDS(so_ns, "sleuth/so_ns.rds")
# so_ns = readRDS("sleuth/so_ns.rds")

#======================
# PULL TPMS
#======================
results_tpm <- kallisto_table(so, use_filtered = FALSE) %>%
  select(Gene = target_id, sample, rep, clone, condition, treatment, time, tpm)

write_tsv(results_tpm, "sleuth/B2B_DKO_tpm.txt")

#==========================
# ASSEMBLE MULTIVARIATE DEA
#==========================
a <- sleuth_results(so, "reduced:full", "lrt") %>%
  select(Gene = target_id, FDR_full = qval) %>%
  distinct()
b <- sleuth_results(so, "condition:full", "lrt") %>%
  select(Gene = target_id, FDR_treatment = qval) %>%
  distinct()
c <- sleuth_results(so, "treatment:full", "lrt") %>%
  select(Gene = target_id, FDR_condition = qval) %>%
  distinct()
d <- sleuth_results(so, "conditionDKO:treatmentForm", "wt") %>%
  select(Gene = target_id, diff_DKO_Form = b, FDR_DKO_Form = qval) %>%
  distinct()

results_dea_multi <- left_join(a, b, by = "Gene") %>%
  left_join(c, by = "Gene") %>%
  left_join(d, by = "Gene") %>%
  arrange(Gene) %>%
  mutate(across(contains("FDR"), ~ if_else(is.na(.), 1, .))) %>%
  mutate(across(contains("diff"), ~ if_else(is.na(.), 0, .)))

rm(a,b,c,d)

write_tsv(results_dea_multi, "sleuth/B2B_DKO_multivariateDEA.txt")

#==========================
# ASSEMBLE UNIVARIATE DEA
#==========================
a <- sleuth_results(so_wt, "treatmentForm") %>%
  select(Gene = target_id, log2fold = b, FDR = qval) %>%
  mutate(log2fold = log2fold / log(2)) %>%
  distinct() %>%
  mutate(condition = "WT", treatment = "Form", .before = 2)

b <- sleuth_results(so_dko, "treatmentForm") %>%
  select(Gene = target_id, log2fold = b, FDR = qval) %>%
  mutate(log2fold = log2fold / log(2)) %>%
  distinct() %>%
  mutate(condition = "DKO", treatment = "Form", .before = 2)

results_dea_treatment <- rbind(a,b) %>%
  arrange(Gene) %>%
  mutate(log2fold = if_else(is.na(log2fold), 0, log2fold),
         FDR = if_else(is.na(FDR), 1, FDR)) %>%
  mutate(fold = 2^log2fold, .before = 4) %>%
  mutate(sig = if_else(fold >= 1.5 & FDR < 0.05, "up",
                       if_else(fold <= 0.67 & FDR < 0.05, "dn", "ns"))) %>%
  group_by(Gene) %>%
  filter(!all(fold == 1)) %>%
  ungroup()

rm(a,b)

write_tsv(results_dea_treatment, "sleuth/B2B_DKO_univariateDEA.txt")

#================================
# ASSEMBLE BASAL EXPRESSION CHECK
#================================
a = sleuth_results(so_ns, "conditionDKO") %>%
  select(Gene = target_id, log2fold = b, FDR = qval) %>%
  mutate(log2fold = log2fold / log(2)) %>%
  distinct() %>%
  mutate(treatment = "NS", condition = "DKO", .before = 2)

results_dea_basal = a %>%
  arrange(Gene) %>%
  mutate(log2fold = if_else(is.na(log2fold), 0, log2fold),
         FDR = if_else(is.na(FDR), 1, FDR)) %>%
  mutate(fold = 2^log2fold, .before = 4) %>%
  mutate(sig = if_else(fold >= 1.5 & FDR < 0.05, "up",
                       if_else(fold <= 0.67 & FDR < 0.05, "dn", "ns"))) %>%
  group_by(Gene) %>%
  filter(!all(fold == 1)) %>%
  ungroup()

write_tsv(results_dea_basal, "sleuth/B2B_NS_univariateDEA (effect of dKO on basal expression).txt")
```

```R transcript level analysis
.libPaths("/home/alex.gao1/R")
setwd("/work/giembycz_lab/B2B_dKO_analysis")

library(dplyr)
library(readr)
library(sleuth)

#======================
# PREP INPUTS
#======================

t2g = read_tsv("Homo_sapiens.GRCh38.p13.cdna.all.070121.mart_export.txt") %>%
  select(target_id = 1, Gene = 2) %>%
  na.omit() %>% 
  distinct()
t2g <- t2g[!duplicated(t2g$target_id), ]

s2c = read_tsv("meta.txt") %>%
    mutate(
        treatment = factor(treatment, levels = c("NS", "Form")),
        condition = factor(condition, levels = c("WT", "DKO")))

new_filter <- function(row, min_reads = 5, min_prop = 0.2){mean(row >= min_reads) >= min_prop}

#======================
# CREATE SLEUTH OBJECTS
#======================

so_transcripts = sleuth_prep(
    sample_to_covariates = s2c,
    target_mapping = t2g,
    filter_fun = new_filter,
    num_cores = 4
)

saveRDS(so_transcripts, "sleuth/so_transcripts.rds")

#======================
# PULL TPMS
#======================
results_tpm <- kallisto_table(so_transcripts, use_filtered = FALSE) %>%
                              select(target_id, sample, rep, clone, condition, treatment, time, tpm) %>%
                              left_join(t2g, by = "target_id") %>%
                              filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
                              select(Gene, transcript = target_id, sample, rep, clone, condition, treatment, time, tpm) %>%
                              arrange(Gene, transcript, condition, treatment, time)

write_tsv(results_tpm, "sleuth/B2B_DKO_transcript_tpm.txt")
```

## EXPORTS

cd C:\Users\alexg\"OneDrive - University of Calgary"\"seqProject B2B_AKO_BKO_dKO"\data_dKO\raw
scp alex.gao1@arc.ucalgary.ca:/work/giembycz_lab/B2B_dKO_analysis/sleuth/B2B_DKO_transcript_tpm.txt ./
