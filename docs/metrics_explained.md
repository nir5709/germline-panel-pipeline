# Metrics explained

## sequencing_depth_gb
Approximate total sequenced bases in gigabases, estimated from `samtools stats` using:
raw total sequences × average read length

## unique_read_enrichment
Fraction of aligned primary reads that are high-quality, non-duplicate, and overlap target regions.

## pct_duplicate_aligned_reads
Percent duplication estimated from `MarkDuplicatesSpark` metrics (`PERCENT_DUPLICATION × 100`).

## unique_base_enrichment
Fraction of all aligned bases represented by high-quality, non-duplicate bases overlapping target regions.

## total_data_gb
Approximate total sequenced bases in gigabases.

## total_aligned_bases
Aligned bases estimated from `samtools stats` (`bases mapped (cigar)`).

## mean_target_coverage_depth
Mean on-target coverage depth from Picard `CollectHsMetrics`.

## uniformity_pct_gt_0p2x_mean
Percent of target bases with depth ≥ 0.2 × mean target coverage.

## pct_target_bases_50x
Percent of target bases covered at ≥50×, reported on a 0–100 scale.

## median_insert_size
Median insert size from Picard `CollectInsertSizeMetrics`.
