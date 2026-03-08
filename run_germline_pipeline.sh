#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Germline Panel Pipeline
# Author: Nihar Garg
#
# Overview:
#   Paired-end FASTQ -> trimming -> alignment -> duplicate marking ->
#   BQSR -> per-sample germline calling -> QC and target enrichment metrics
#
# Usage:
#   bash run_germline_pipeline.sh config/config.example.sh
#
# Notes:
#   - Designed for targeted germline panel sequencing.
#   - Current implementation performs per-sample calling with HaplotypeCaller.
#   - This is not a joint-genotyping workflow.
#   - Written for hg19 / GRCh37-style resources.
###############################################################################

usage() {
  cat <<EOF
Usage:
  $0 <config.sh>

Example:
  $0 config/config.example.sh
EOF
}

[[ $# -eq 1 ]] || { usage; exit 1; }

CONFIG="$1"
[[ -f "$CONFIG" ]] || { echo "[FATAL] Config not found: $CONFIG" >&2; exit 1; }

# shellcheck source=/dev/null
source "$CONFIG"

###############################################################################
# PRECHECKS
###############################################################################
log(){ echo "[$(date '+%F %T')] $*"; }
need(){ [[ -s "$1" ]] || { echo "[FATAL] Missing: $1" >&2; exit 1; }; }

: "${PROJECT_DIR:?Missing PROJECT_DIR in config}"
: "${INPUT_DIR:?Missing INPUT_DIR in config}"
: "${REF_DIR:?Missing REF_DIR in config}"
: "${REF:?Missing REF in config}"
: "${DBSNP:?Missing DBSNP in config}"
: "${WORK:?Missing WORK in config}"

THREADS="${THREADS:-16}"
FASTP_THREADS="${FASTP_THREADS:-4}"
PAIRHMM_THREADS="${PAIRHMM_THREADS:-4}"

TARGETS="${TARGETS:-}"

need "$REF"
need "${REF}.fai"
need "$DBSNP"
if [[ -n "$TARGETS" ]]; then
  need "$TARGETS"
fi

mkdir -p "$WORK"/{qc,trim,bam,metrics,bqsr,vcf,logs,mosdepth,tmp}

COHORT_TSV="$WORK/metrics/all_samples.qc_metrics.tsv"
if [[ ! -s "$COHORT_TSV" ]]; then
  echo -e "sample\tsequencing_depth_gb\tunique_read_enrichment\tpct_duplicate_aligned_reads\tunique_base_enrichment\ttotal_data_gb\ttotal_aligned_bases\tmean_target_coverage_depth\tuniformity_pct_gt_0p2x_mean\tpct_target_bases_50x\tmedian_insert_size" > "$COHORT_TSV"
fi

for R1 in "$INPUT_DIR"/*_R1_001.fastq.gz; do
  [[ -e "$R1" ]] || continue
  R2="${R1/_R1_001/_R2_001}"
  SID="$(basename "$R1" | sed 's/_R1_001\.fastq\.gz//')"
  [[ -s "$R2" ]] || { log "[FATAL] Missing pair FASTQ for $SID"; exit 1; }

  log "=== Processing $SID ==="

  fastp \
    -i "$R1" -I "$R2" \
    -o "$WORK/trim/${SID}_R1.trim.fq.gz" \
    -O "$WORK/trim/${SID}_R2.trim.fq.gz" \
    --detect_adapter_for_pe \
    --thread "$FASTP_THREADS" \
    --html "$WORK/qc/${SID}.fastp.html" \
    --json "$WORK/qc/${SID}.fastp.json"

  bwa mem -t "$THREADS" -R "@RG\tID:${SID}\tSM:${SID}\tPL:ILLUMINA\tLB:lib1" \
      "$REF" "$WORK/trim/${SID}_R1.trim.fq.gz" "$WORK/trim/${SID}_R2.trim.fq.gz" \
    | samtools sort -@ 4 -o "$WORK/bam/${SID}.sorted.bam" -
  samtools index "$WORK/bam/${SID}.sorted.bam"

  gatk --java-options "-Xmx8g" MarkDuplicatesSpark \
      -I "$WORK/bam/${SID}.sorted.bam" \
      -O "$WORK/bam/${SID}.dedup.bam" \
      -M "$WORK/metrics/${SID}.dedup.metrics.txt" \
      --conf 'spark.executor.cores=1'
  samtools index "$WORK/bam/${SID}.dedup.bam"

  gatk --java-options "-Xmx8g" BaseRecalibrator \
      -R "$REF" -I "$WORK/bam/${SID}.dedup.bam" \
      --known-sites "$DBSNP" \
      -O "$WORK/bqsr/${SID}.recal.table"

  gatk --java-options "-Xmx8g" ApplyBQSR \
      -R "$REF" -I "$WORK/bam/${SID}.dedup.bam" \
      --bqsr-recal-file "$WORK/bqsr/${SID}.recal.table" \
      -O "$WORK/bam/${SID}.recal.bam"
  samtools index "$WORK/bam/${SID}.recal.bam"

  EXTRA_ARGS=()
  if [[ -n "$TARGETS" && -s "$TARGETS" ]]; then
    EXTRA_ARGS=(-L "$TARGETS")
  fi

  gatk --java-options "-Xmx8g" HaplotypeCaller \
      -R "$REF" -I "$WORK/bam/${SID}.recal.bam" \
      -O "$WORK/vcf/${SID}.raw.vcf.gz" \
      --native-pair-hmm-threads "$PAIRHMM_THREADS" \
      "${EXTRA_ARGS[@]}"
  bcftools index -t -f "$WORK/vcf/${SID}.raw.vcf.gz"

  BAM="$WORK/bam/${SID}.recal.bam"
  DEDUP_MET="$WORK/metrics/${SID}.dedup.metrics.txt"
  METDIR="$WORK/metrics"
  MOSDIR="$WORK/mosdepth"
  TMPDIR="$WORK/tmp"

  samtools stats -@ 8 "$BAM" > "$TMPDIR/${SID}.stats.txt"
  RAW_READS=$(grep "^raw total sequences:" "$TMPDIR/${SID}.stats.txt" | awk '{print $4}')
  AVG_LEN=$(grep "^average length:" "$TMPDIR/${SID}.stats.txt" | awk '{print $4}')
  TOTAL_BASES=$(awk -v r="$RAW_READS" -v l="$AVG_LEN" 'BEGIN{printf "%.0f", r*l}')
  TOTAL_GB=$(awk -v b="$TOTAL_BASES" 'BEGIN{printf "%.3f", b/1e9}')
  ALN_BASES=$(grep "^bases mapped (cigar):" "$TMPDIR/${SID}.stats.txt" | awk '{gsub(",","",$5); print $5}')

  PCT_DUP=$(awk 'BEGIN{FS="\t"} $1!~/^#/ && NR>1 {print $9}' "$DEDUP_MET" 2>/dev/null || echo "0")
  PCT_DUP100=$(awk -v p="$PCT_DUP" 'BEGIN{printf "%.3f", 100.0*p}')

  picard CollectInsertSizeMetrics \
    I="$BAM" O="$METDIR/${SID}.insert_metrics.txt" H="$METDIR/${SID}.insert_hist.pdf" M=0.5 \
    VALIDATION_STRINGENCY=SILENT >/dev/null 2>&1
  MED_INS=$(awk 'BEGIN{FS="\t"} $1!~/^#/ && NR>1 {print $1}' "$METDIR/${SID}.insert_metrics.txt")

  MEAN_COV="NA"; PCT50X="NA"; ON_TAR_BASES="NA"; PF_BASES_ALIGNED="NA"
  if [[ -n "$TARGETS" && -s "$TARGETS" ]]; then
    picard CollectHsMetrics \
      I="$BAM" O="$METDIR/${SID}.hs_metrics.txt" R="$REF" \
      BAIT_INTERVALS="$TARGETS" TARGET_INTERVALS="$TARGETS" \
      VALIDATION_STRINGENCY=SILENT >/dev/null 2>&1

    read MEAN_COV PCT50X ON_TAR_BASES PF_BASES_ALIGNED <<<"$(
      awk 'BEGIN{FS="\t"}
        $0 ~ /^#/{next}
        NR==1 {hdr=$0; next}
        NR==2{
          split(hdr,h,"\t"); for(i=1;i<=NF;i++){idx[$i]=i}
          printf "%.3f\t%.6f\t%d\t%d\n", $(idx["MEAN_TARGET_COVERAGE"]), $(idx["PCT_TARGET_BASES_50X"]), $(idx["ON_TARGET_BASES"]), $(idx["PF_BASES_ALIGNED"])
        }' "$METDIR/${SID}.hs_metrics.txt"
    )"
  fi
  PCT50X_100=$(awk -v v="$PCT50X" 'BEGIN{ if(v=="NA"||v=="") print "NA"; else printf "%.3f", 100.0*v }')

  UNIFORMITY="NA"; UNIQUE_READ_ENRICH="NA"; UNIQUE_BASE_ENRICH="NA"
  if [[ -n "$TARGETS" && -s "$TARGETS" ]]; then
    mkdir -p "$MOSDIR"

    mosdepth --by "$TARGETS" --thresholds 50 --fast-mode --mapq 0 --threads 8 \
      "$MOSDIR/${SID}.all" "$BAM" >/dev/null 2>&1

    samtools view -@ 8 -b -q 30 -F 1024 "$BAM" > "$TMPDIR/${SID}.uniq.bam"
    samtools index -@ 8 "$TMPDIR/${SID}.uniq.bam"

    mosdepth --by "$TARGETS" --fast-mode --mapq 0 --threads 8 \
      "$MOSDIR/${SID}.uniq" "$TMPDIR/${SID}.uniq.bam" >/dev/null 2>&1

    UQ_ONTARGET_READS=$(samtools view -@ 8 -c -L "$TARGETS" "$TMPDIR/${SID}.uniq.bam")
    TOTAL_ALIGNED_READS=$(samtools view -@ 8 -c -F 0x904 "$BAM")
    UNIQUE_READ_ENRICH=$(awk -v a="$UQ_ONTARGET_READS" -v b="$TOTAL_ALIGNED_READS" 'BEGIN{ if(b>0) printf "%.6f", a/b; else print "NA"}')

    UQ_BASES_ON_TARGET=$(zcat "$MOSDIR/${SID}.uniq.regions.bed.gz" | awk '{len=$3-$2; d=$4; sum+=len*d} END{printf "%.0f", sum}')
    UNIQUE_BASE_ENRICH=$(awk -v u="$UQ_BASES_ON_TARGET" -v a="$ALN_BASES" 'BEGIN{ if(a>0) printf "%.6f", u/a; else print "NA"}')

    if [[ "$MEAN_COV" != "NA" ]]; then
      THR=$(awk -v m="$MEAN_COV" 'BEGIN{printf "%d", (0.2*m)+0.5}')
      mosdepth --by "$TARGETS" --thresholds "$THR" --fast-mode --mapq 0 --threads 8 \
        "$MOSDIR/${SID}.uni" "$BAM" >/dev/null 2>&1

      UNIFORMITY=$(zcat "$MOSDIR/${SID}.uni.thresholds.bed.gz" | \
        awk -v thrcol=5 'BEGIN{tot=0; ge=0} {len=$3-$2; tot+=len; if($thrcol!="NA"){ge += len*$thrcol}} END{ if(tot>0) printf "%.3f", 100.0*ge/tot; else print "NA"}')
    fi
  fi

  SAMPLE_TSV="$METDIR/${SID}.qc_metrics.tsv"
  echo -e "sample\tsequencing_depth_gb\tunique_read_enrichment\tpct_duplicate_aligned_reads\tunique_base_enrichment\ttotal_data_gb\ttotal_aligned_bases\tmean_target_coverage_depth\tuniformity_pct_gt_0p2x_mean\tpct_target_bases_50x\tmedian_insert_size" > "$SAMPLE_TSV"
  echo -e "${SID}\t${TOTAL_GB}\t${UNIQUE_READ_ENRICH}\t${PCT_DUP100}\t${UNIQUE_BASE_ENRICH}\t${TOTAL_GB}\t${ALN_BASES}\t${MEAN_COV}\t${UNIFORMITY}\t${PCT50X_100}\t${MED_INS}" >> "$SAMPLE_TSV"

  tail -n 1 "$SAMPLE_TSV" >> "$COHORT_TSV"

  rm -f "$TMPDIR/${SID}.uniq.bam" "$TMPDIR/${SID}.uniq.bam.bai" 2>/dev/null || true

  log "=== Finished $SID ==="
done

log "All done. Outputs are in $WORK/"
log "Per-sample QC: $WORK/metrics/<SID>.qc_metrics.tsv"
log "Cohort QC:     $WORK/metrics/all_samples.qc_metrics.tsv"
