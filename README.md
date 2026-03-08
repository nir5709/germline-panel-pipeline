# germline-panel-pipeline

![GitHub repo size](https://img.shields.io/github/repo-size/nir5709/germline-panel-pipeline)
![GitHub last commit](https://img.shields.io/github/last-commit/nir5709/germline-panel-pipeline)
![GitHub license](https://img.shields.io/github/license/nir5709/germline-panel-pipeline)
![Platform](https://img.shields.io/badge/platform-Linux-blue)
![Genome](https://img.shields.io/badge/genome-hg19%20%2F%20GRCh37-green)

A targeted **germline variant calling and QC workflow** for paired-end sequencing data using BWA-MEM, GATK HaplotypeCaller, and integrated panel-enrichment metrics.

## Overview

This workflow processes paired-end FASTQ files and performs:

1. adapter trimming and quality filtering using **fastp**
2. alignment to **hg19 / GRCh37** using **BWA-MEM**
3. BAM sorting and indexing using **SAMtools**
4. duplicate marking using **GATK MarkDuplicatesSpark**
5. base quality score recalibration using **GATK BQSR**
6. per-sample germline variant calling using **GATK HaplotypeCaller**
7. target-enrichment and coverage QC calculation
8. cohort-level QC table generation

## Pipeline Workflow

```text
FASTQ
↓
fastp trimming
↓
BWA-MEM alignment
↓
SAMtools sorting
↓
MarkDuplicatesSpark
↓
BQSR
↓
HaplotypeCaller
↓
Per-sample QC metrics
↓
Cohort QC summary
```

## Intended Use

This pipeline is designed for:

- targeted germline panel sequencing
- per-sample germline SNV / indel calling
- post-alignment QC and enrichment reporting
- research workflows and assay-development pipelines

This is **not** a joint-genotyping workflow. Each sample is processed independently with HaplotypeCaller.

## Key Features

- per-sample germline variant calling
- integrated QC and enrichment reporting
- target capture metrics using Picard HsMetrics
- mosdepth-based uniformity and enrichment calculations
- per-sample and cohort-level summary TSV files
- config-driven setup without hardcoded machine-specific paths

## Repository Structure

```text
germline-panel-pipeline/
├── run_germline_pipeline.sh
├── README.md
├── LICENSE
├── .gitignore
├── config/
│   └── config.example.sh
├── templates/
│   └── sample_naming.example.txt
├── docs/
│   └── metrics_explained.md
└── env/
    └── tool_versions.txt
```

## Input Expectations

The pipeline expects FASTQ files in `INPUT_DIR` with names like:

```text
SAMPLE001_R1_001.fastq.gz
SAMPLE001_R2_001.fastq.gz
```

Sample IDs are derived from the FASTQ filenames.

## Outputs

For each sample, the workflow creates:

### Processed reads and alignment
- trimmed FASTQ files
- sorted BAM
- duplicate-marked BAM
- recalibrated BAM

### Variant calls
- per-sample raw VCF (`*.raw.vcf.gz`)

### QC and metrics
- fastp HTML and JSON reports
- duplicate metrics
- insert size metrics
- HsMetrics output
- per-sample QC TSV

### Cohort summary
- `work/metrics/all_samples.qc_metrics.tsv`

## QC Metrics Reported

The cohort and sample-level summary tables include:

- `sequencing_depth_gb`
- `unique_read_enrichment`
- `pct_duplicate_aligned_reads`
- `unique_base_enrichment`
- `total_data_gb`
- `total_aligned_bases`
- `mean_target_coverage_depth`
- `uniformity_pct_gt_0p2x_mean`
- `pct_target_bases_50x`
- `median_insert_size`

Detailed metric definitions are described in `docs/metrics_explained.md`.

## Methodological Approach

This pipeline combines alignment, duplicate marking, recalibration, and per-sample calling with panel-focused QC evaluation.

Coverage and enrichment assessment includes:

- total aligned bases from `samtools stats`
- duplicate fraction from `MarkDuplicatesSpark`
- insert size distribution from Picard
- target capture metrics from Picard `CollectHsMetrics`
- mosdepth-based target depth and uniformity estimation
- unique read and unique base enrichment on targets

Uniformity is reported as:

> percentage of target bases with depth ≥ 0.2 × mean target coverage

## Software Requirements

This workflow assumes access to:

- fastp
- bwa
- samtools
- gatk
- bcftools
- picard
- mosdepth
- bedtools
- awk
- gzip / zcat

The pipeline is intended to run with **latest stable tool versions**.

## Reference Build

This workflow is written for:

```text
Human genome reference: hg19 / GRCh37
```

If adapting to GRCh38, update:

- reference FASTA
- reference indexes
- dbSNP resource
- target BED / interval resource
- any expected downstream annotation resources

## Installation

```bash
git clone https://github.com/nir5709/germline-panel-pipeline.git
cd germline-panel-pipeline
```

## Configuration

```bash
cp config/config.example.sh config/my_run_config.sh
nano config/my_run_config.sh
```

Update all paths before running.

## Running the Pipeline

```bash
bash run_germline_pipeline.sh config/my_run_config.sh
```

## Notes on Target Metrics

If `TARGETS` is left empty, the workflow will still run variant calling, but target-based enrichment and uniformity metrics will be reported as `NA`.

## Disclaimer

This pipeline is intended for **research use and assay development**.

Clinical deployment requires:

- local validation
- QC thresholds
- regulatory compliance
- reporting SOPs
- appropriate accreditation and oversight

## Author

**Nihar Garg**  
Bioinformatics | Cancer Genomics

## Citation

If you use this pipeline in research, please cite:

Nihar Garg.  
Germline Panel Pipeline.  
GitHub repository.  
https://github.com/nir5709/germline-panel-pipeline

## License

MIT License
