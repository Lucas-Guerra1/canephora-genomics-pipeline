# Coffee Genomics Pipeline — SNP Calling & Phylogenetic Analysis

A bioinformatics pipeline for SNP detection, variant filtering, and phylogenetic inference in *Coffea canephora*, developed as part of the M.Sc. dissertation *"Genomic Characterization and Phylogenetic Analysis of Coffea canephora Accessions from the UFLA Germplasm Bank"* (Federal University of Lavras — UFLA, 2026).

The pipeline integrates data from two genotyping platforms — Illumina NGS (FASTQ) and Axiom 26K SNP microarray — enabling cross-platform comparison and phylogenetic placement of germplasm accessions alongside diversity center references.

> **Portfolio notice:** This repository documents the bioinformatics pipeline developed during my M.Sc. dissertation for portfolio and reproducibility demonstration purposes. The underlying genomic data (BAG/UFLA germplasm and Embrapa cultivar genotypes) are confidential and are not included. All rights reserved — see [License](#license).

---

## Table of Contents

- [Overview](#overview)
- [Workflow](#workflow)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input Data](#input-data)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
- [Output Files](#output-files)
- [Reproducibility Validation](#reproducibility-validation)
- [Third-Party Scripts](#third-party-scripts)
- [NCBI Accessions](#ncbi-accessions)
- [Limitations](#limitations)
- [Citation](#citation)
- [License](#license)

---

## Overview

This pipeline was designed to:

1. Process raw Illumina paired-end FASTQ files from *C. canephora* diversity center accessions
2. Align reads to the *C. canephora* CC1.8 reference genome
3. Call SNPs and InDels using GATK (v3.8)
4. Filter variants using BCFtools and VCFtools
5. Extract SNP positions matching the Axiom® 26K chip (EMBR5SNP) to enable cross-platform integration
6. Generate neighbor-joining phylogenetic trees using TASSEL
7. Validate pipeline reproducibility via Cohen's Kappa analysis (R)

The pipeline achieved **~90.6% recovery** of the 25,456 SNPs present in the Axiom® 26K chip from NGS data, with a **Kappa coefficient of 0.92** between NGS-derived and chip-derived genotypes, confirming near-perfect methodological agreement.

---

## Workflow

```
FASTQ (R1/R2)
     │
     ▼
Adapter Trimming (TrimAdaptor.pl — Gong, 2012)
     │
     ▼
Quality Filtering (Filter.pl — Guignon, CIRAD)
     │
     ▼
Pair Adjustment (Compare.pl — Sabot, IRD)
     │
     ▼
BWA Index → BWA-MEM Alignment → SAM→BAM (Samtools)
     │
     ▼
QC: flagstat → Remove unmapped reads → Sort → Add Read Groups (Picard/GATK)
     │
     ▼
Merge BAMs (MergeSamFiles) → Index BAM + Reference (samtools faidx, GATK dict)
     │
     ▼
InDel Realignment (RealignerTargetCreator → IndelRealigner — GATK 3.8)
     │
     ▼
Variant Calling (UnifiedGenotyper — GATK 3.8)
     │
     ▼
Quality Filtering (BCFtools: QD, FS, MQ, MQRankSum, ReadPosRankSum)
     │
     ▼
Depth Filtering + Biallelic selection (VCFtools: --minDP 10)
     │
     ▼
Position Extraction (VCFtools: --positions — Axiom 26K loci)
     │
     ▼
MAF Filtering (VCFtools: --maf 0.4)
     │
     ▼
High-Quality SNP Matrix → TASSEL Neighbor-Joining Tree
```

---

## Requirements

### System

- Linux (Ubuntu 18.04+ recommended)
- Docker (for GATK 3.8)
- Perl ≥ 5.x

### Core tools

| Tool | Version | Purpose |
|---|---|---|
| BWA | ≥ 0.7.17 | Reference indexing and alignment |
| Samtools | ≥ 1.9 | BAM processing, QC, indexing |
| GATK | **3.8-1** (via Docker) | InDel realignment, variant calling |
| BCFtools | ≥ 1.9 | Variant quality filtering |
| VCFtools | ≥ 0.1.16 | Depth filtering, position selection, MAF |
| TASSEL | ≥ 5.0 | Neighbor-joining phylogenetic trees |
| Perl | ≥ 5.x | Pre-processing scripts |
| R | ≥ 4.0 | Kappa validation (epiR package) |

> **Important:** This pipeline uses **GATK 3.8**, not GATK 4.x. The tools `RealignerTargetCreator`, `IndelRealigner`, and `UnifiedGenotyper` were deprecated and removed in GATK 4. Use the Docker image specified below.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/Lucas-Guerra1/coffee-genomics-pipeline.git
cd coffee-genomics-pipeline
```

### 2. Install BWA, Samtools, BCFtools, VCFtools via Conda

```bash
conda env create -f envs/bioinfo.yml
conda activate bioinfo
```

### 3. Pull GATK 3.8 Docker image

```bash
docker pull broadinstitute/gatk:3.8-1
```

Verify:

```bash
sudo docker run -v ~/gatk:/gatk/my_data -it broadinstitute/gatk:3.8-1 \
  java -jar /usr/GenomeAnalysisTK.jar --version
```

### 4. Install TASSEL

Download from [tassel.sourceforge.net](https://tassel.sourceforge.net/) and add to PATH.

### 5. Install R dependencies

```r
install.packages("epiR")
```

### 6. Obtain third-party Perl scripts

The pre-processing scripts (TrimAdaptor.pl, Filter.pl, Compare.pl) were developed by external authors and are not redistributed here. See [Third-Party Scripts](#third-party-scripts) for contact information.

---

## Input Data

### Reference genome

The **CC1.8** genome assembly (Salojarvi et al., 2024) was selected as reference after benchmarking against CC Science (Denoeud et al., 2014). CC1.8 showed superior mapping quality:

| Metric | CC1.8 | CC Science |
|---|---|---|
| Properly paired reads | 84.89% | 79.64% |
| Singletons | 0.58% | 0.81% |
| Secondary alignments | 414,728 | 503,487 |

Place the reference FASTA at:

```
data/CC1.8_v2_pseudomolecule_cat.fasta
```

### FASTQ reads

Paired-end FASTQ files from *C. canephora* diversity center accessions. NCBI SRA accessions used in this study are listed in the [NCBI Accessions](#ncbi-accessions) section.

Place files at:

```
data/reads/SAMPLEID_R1.fastq
data/reads/SAMPLEID_R2.fastq
```

### Axiom 26K SNP positions

The file `data/Affy_Positions_CC1.8_v2_Filter.txt` contains the genomic positions of the 25,456 SNPs from the Axiom® EMBR5SNP chip (Carneiro et al., 2019), used for cross-platform SNP extraction.

---

## Usage

Each pipeline step is documented as a numbered shell script in `pipeline/`. Run them sequentially:

```bash
bash pipeline/01_adapter_trimming.sh
bash pipeline/02_quality_filter.sh
bash pipeline/03_pair_adjustment.sh
bash pipeline/04_bwa_index.sh
bash pipeline/05_bwa_align.sh
bash pipeline/06_samtools_qc.sh
bash pipeline/07_gatk_preprocessing.sh
bash pipeline/08_gatk_variant_calling.sh
bash pipeline/09_bcftools_vcftools_filter.sh
```

For phylogenetic analysis, open TASSEL and import the final VCF, or use the command-line interface:

```bash
bash pipeline/10_tassel_phylogeny.sh
```

For Kappa validation:

```r
source("r_analysis/kappa_validation.R")
```

---

## Pipeline Steps

### Step 1 — Adapter trimming

Removes sequencing adapters using `TrimAdaptor.pl` (Gong, 2012):

```bash
perl scripts/TrimAdaptor.pl \
  -f1 YOURFILE_R1.fastq \
  -f2 YOURFILE_R2.fastq \
  -a1 AGATCGGAAGAGCAC \
  -a2 AGATCGGAAGAGCGT \
  -trim N -m 40 > counts.txt
```

### Step 2 — Quality filtering

Filters reads by mean Phred quality ≥ 30 and minimum length ≥ 35 using `Filter.pl` (Guignon, CIRAD):

```bash
perl scripts/Filter.pl \
  -f Clones/Filtered-YOURFILE_R1.fastq \
  -o Clones/Filtered-R1_001.fastq
```

### Step 3 — Pair adjustment

Re-pairs filtered reads using `Compare.pl` (Sabot, IRD):

```bash
perl scripts/Compare.pl \
  -f Clones/Filtered-R1_001.fastq \
  -r Clones/Filtered-R2_001.fastq \
  -of Clones/Paired_Filtered-R1.fastq \
  -or Clones/Paired_Filtered-R2.fastq \
  -os Clones/Single-Paired_Filtered-R1
```

### Step 4 — Reference indexing

```bash
bwa index CC1.8_v2_pseudomolecule_cat.fasta
```

### Step 5 — Alignment

```bash
bwa mem CC1.8/CC1.8_v2_pseudomolecule_cat.fasta \
  -M -B 4 \
  Clones/Paired_Filtered-YOURFILE_R1.fastq \
  Clones/Paired_Filtered-YOURFILE_R2.fastq \
  > Clones/BWA/Sam-CC1.8_Paired_Filtered_YOURFILE.sam
```

### Step 6 — SAM→BAM, QC, filter unmapped reads

```bash
# Convert and raw QC
samtools view -h -b -S Sam-CC1.8.sam -o Clones/BWA/raw-CC1.8.bam
samtools flagstat -O tsv raw-CC1.8.bam > flagstat-raw-CC1.8.txt

# Remove unmapped reads
samtools view -h -b -S -f 0x2 raw-CC1.8.bam -o Clones/BWA/Clean-Sam-CC1.8.bam

# Sort
samtools sort Clean-Sam-CC1.8.bam -o Clones/BWA/Sorted-Clean-Sam-CC1.8.bam
```

### Step 7 — GATK preprocessing

```bash
# Add Read Groups (Picard via GATK Docker)
sudo docker run -v ~/gatk:/gatk/my_data -it broadinstitute/gatk:3.8-1 \
  gatk AddOrReplaceReadGroups \
  -I /gatk/my_data/Clones/BWA/Sorted-Clean-Sam-CC1.8.bam \
  -O /gatk/my_data/RG-Sorted-Clean-Sam-CC1.8.bam \
  -RGLB lib1 -RGPL Illumina -RGPU unit1 -RGSM YOURFILE \
  --VALIDATION_STRINGENCY SILENT

# Merge all samples
gatk MergeSamFiles \
  -I /gatk/my_data/RG-Sorted-Clean-Sam-CC1.8.bam \
  -O /gatk/my_data/YOURFILE+YOURFILE2_Merged.bam

# Index merged BAM and reference
samtools index ~/gatk/YOURFILE+YOURFILE2_Merged.bam
samtools faidx ~/gatk/CC1.8_v2_pseudomolecule_cat.fasta

# Create sequence dictionary
gatk CreateSequenceDictionary \
  -R /gatk/my_data/CC1.8_v2_pseudomolecule_cat.fasta \
  -O /gatk/my_data/CC1.8_v2_pseudomolecule_cat.dict
```

### Step 8 — Variant calling (GATK 3.8)

```bash
# Identify InDel regions
sudo docker run -v ~/gatk:/gatk/my_data -it broadinstitute/gatk:3.8-1 \
  java -jar /usr/GenomeAnalysisTK.jar \
  -T RealignerTargetCreator \
  -R /gatk/my_data/CC1.8_v2_pseudomolecule_cat.fasta \
  -I /gatk/my_data/YOURFILE+YOURFILE2_Merged.bam \
  -o /gatk/my_data/YOURFILE+YOURFILE2_INDEL.intervals

# Local realignment
java -jar /usr/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  -R /gatk/my_data/CC1.8_v2_pseudomolecule_cat.fasta \
  -I /gatk/my_data/YOURFILE+YOURFILE2_Merged.bam \
  -targetIntervals /gatk/my_data/YOURFILE+YOURFILE2_INDEL.intervals \
  -o /gatk/my_data/Realign_YOURFILE+YOURFILE2_Merged.bam

# SNP + InDel calling
java -jar /usr/GenomeAnalysisTK.jar \
  -T UnifiedGenotyper \
  -I /gatk/my_data/Realign_YOURFILE+YOURFILE2_Merged.bam \
  -R /gatk/my_data/CC1.8_v2_pseudomolecule_cat.fasta \
  -o /gatk/my_data/YOURFILE+YOURFILE2_Calling.vcf \
  -glm BOTH
```

### Step 9 — Variant filtering

```bash
# BCFtools quality filter
bcftools filter \
  -e "INFO/QD < 2.0 | INFO/FS > 60.0 | INFO/MQ < 40.0 | INFO/MQRankSum < -12.5 | INFO/ReadPosRankSum < -8.0" \
  -o YOURFILE+YOURFILE2_Calling_filtGATK.vcf \
  YOURFILE+YOURFILE2_Calling.vcf

# VCFtools: minimum depth filter
vcftools --gzvcf YOURFILE+YOURFILE2_Calling_filtGATK.vcf.gz \
  --minDP 10 \
  --out YOURFILE+YOURFILE2_DP10 \
  --recode-INFO-all --recode

# VCFtools: extract Axiom 26K positions
vcftools --vcf YOURFILE+YOURFILE2_DP10.recode.vcf \
  --positions data/Affy_Positions_CC1.8_v2_Filter.txt \
  --remove-indels \
  --min-alleles 2 --max-alleles 2 \
  --out data/2alleles_noIndels \
  --recode-INFO-all --recode

# VCFtools: MAF filter
vcftools --vcf data/2alleles_noIndels.recode.vcf \
  --maf 0.4 \
  --out data/YOURFILE+YOURFILE2_Affy_test2 \
  --recode-INFO-all --recode
```

---

## Output Files

| File | Description |
|---|---|
| `*.sorted.bam` | Coordinate-sorted BAM file |
| `*.sorted.bam.bai` | BAM index |
| `flagstat-*.txt` | Samtools alignment statistics |
| `*_Calling.vcf` | Raw variant calls (SNPs + InDels) |
| `*_filtGATK.vcf.gz` | Quality-filtered variants |
| `*_DP10.recode.vcf` | Depth-filtered variants |
| `2alleles_noIndels.recode.vcf` | Biallelic SNPs at Axiom 26K positions |
| `*_Affy_test2.recode.vcf` | Final high-quality SNP matrix (MAF ≥ 0.4) |
| `*.nwk` | Neighbor-joining phylogenetic tree (TASSEL) |
| `kappa_results.csv` | Cohen's Kappa cross-platform concordance |

---

## Reproducibility Validation

The pipeline was validated by comparing SNP calls from NGS (this pipeline) against chip-based genotypes from the same sample (Clone 14, Axiom® 26K).

| Comparison | Kappa |
|---|---|
| Chip replicate 1 vs 2 | 0.99 |
| Chip replicate 1 vs 3 | 0.99 |
| Chip replicate 2 vs 3 | 0.99 |
| **NGS pipeline vs chip** | **0.92** |

A Kappa of 0.92 indicates near-perfect agreement (Landis & Koch, 1977), confirming that the pipeline reliably reproduces chip-based genotyping results from NGS data.

Run the Kappa analysis:

```r
source("r_analysis/kappa_validation.R")
```

Requires the `epiR` package: `install.packages("epiR")`

---

## Third-Party Scripts

The pre-processing Perl scripts used in this pipeline were developed by external authors and are **not redistributed** in this repository. Please contact the authors directly to obtain them:

| Script | Author | Institution | Contact |
|---|---|---|---|
| `TrimAdaptor.pl` | George Gong | — | — |
| `Filter.pl` (Filter_Fastq_On_Mean_Quality.pl) | Valentin Guignon | CIRAD | valentin.guignon@cirad.fr |
| `Compare.pl` (Compare_fastq_paired) | François Sabot | IRD | francois.sabot@ird.fr |

Once obtained, place the scripts in the `scripts/` directory.

---

## NCBI Accessions

Raw sequencing data from diversity center accessions are available at NCBI SRA under project **PRJNA803612** (Tournebize et al., 2022):

| NCBI ID | SRA Run | Diversity Group |
|---|---|---|
| 20135 | SRR18336404 | D |
| 20145 | SRR18336405 | D |
| Ramasoni3 | SRR18336420 | B |
| ZO-16 | SRR18336422 | O |
| K4.17.2 | SRR18336425 | R |
| 20731 | SRR18336445 | E |
| 20723 | SRR18336447 | A |
| 20716 | SRR18336448 | C |

Download with SRA Toolkit:

```bash
prefetch SRR18336404
fastq-dump --split-files SRR18336404
```

---

## Limitations

- Uses **GATK 3.8** (legacy); tools `RealignerTargetCreator`, `IndelRealigner`, and `UnifiedGenotyper` are not available in GATK 4.x
- Designed for **diploid** organisms only
- SNP position extraction requires the Axiom® EMBR5SNP chip position file — contact Embrapa Café or the original authors (Carneiro et al., 2019) for access
- BAG/UFLA germplasm genotyping data and Embrapa cultivar data are confidential and are not included in this repository
- Pre-processing Perl scripts must be obtained directly from their authors (see [Third-Party Scripts](#third-party-scripts))

---

## Citation

If you use this pipeline, please cite:

> Guerra, L. R. S. *Genomic characterization and phylogenetic analysis of Coffea canephora accessions from the UFLA germplasm bank*. M.Sc. dissertation, Federal University of Lavras (UFLA), Lavras, Brazil, 2026.

Please also cite the tools used:

- Li & Durbin (2009) — BWA. *Bioinformatics*, 25:1754–1760
- McKenna et al. (2010) — GATK. *Genome Research*, 20:1297–1303
- Danecek et al. (2011) — VCFtools. *Bioinformatics*, 27:2156–2158
- Bradbury et al. (2007) — TASSEL. *Bioinformatics*, 23:2633–2635
- Carneiro et al. (2019) — Axiom 26K chip. UFLA doctoral thesis
- Tournebize et al. (2022) — Diversity center accessions. *Global Change Biology*, 28:4124–4142

---

## Author

**Lucas Ribeiro de Souza Guerra**
Federal University of Lavras (UFLA) — Brazil
Ph.D. candidate | M.Sc. Plant Biotechnology

[GitHub](https://github.com/Lucas-Guerra1)

---

## License

© 2026 Lucas Ribeiro de Souza Guerra. All rights reserved.

This repository is made publicly visible for portfolio purposes only. No part of this code, pipeline, or documentation may be used, copied, modified, merged, published, distributed, sublicensed, or sold without explicit written permission from the author.

For collaboration or licensing inquiries, contact via [GitHub](https://github.com/Lucas-Guerra1).
