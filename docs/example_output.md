# Example Cohort QC Output

The germline panel pipeline generates a cohort-level QC metrics table summarizing sequencing performance, enrichment efficiency, and coverage characteristics for each processed sample.

This table can be used to identify samples suitable for downstream germline variant interpretation.

## Example Output

| Sample ID | Unique Read Enrichment | Mean Target Coverage (×) | Uniformity (%) | Target Bases ≥50× (%) | QC Result |
| --------- | ---------------------- | ------------------------ | -------------- | --------------------- | --------- |
| RH_01     | 0.484                  | 71.50                    | 98.76          | 82.59                 | FAIL      |
| RH_02     | 0.567                  | 41.02                    | 98.49          | 28.42                 | FAIL      |
| RH_03     | 0.590                  | 79.14                    | 98.86          | 88.03                 | FAIL      |
| RH_06     | 0.619                  | 112.98                   | 98.73          | 95.25                 | PASS      |
| RH_08     | 0.739                  | 111.75                   | 98.79          | 96.04                 | PASS      |
| RH_10     | 0.562                  | 108.17                   | 98.73          | 94.92                 | PASS      |
| RH_15     | 0.730                  | 108.29                   | 98.62          | 94.03                 | PASS      |
| RH_16     | 0.524                  | 146.02                   | 98.52          | 97.39                 | PASS      |
| RH_18     | 0.460                  | 116.82                   | 97.90          | 93.12                 | PASS      |

## Description

The cohort QC table summarizes sequencing and enrichment metrics across all processed samples.

Metrics include:

* sequencing depth
* duplicate fraction
* enrichment efficiency
* coverage uniformity
* target coverage thresholds
* insert size statistics

Samples passing overall QC assessment are retained for downstream germline variant analysis.
