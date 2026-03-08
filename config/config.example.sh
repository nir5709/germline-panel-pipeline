# Example configuration for run_germline_pipeline.sh
# Copy this file and edit all paths for your environment.

PROJECT_DIR="/path/to/project"
INPUT_DIR="${PROJECT_DIR}/read"
REF_DIR="${PROJECT_DIR}/reference"

REF="${REF_DIR}/hg19.fa"
DBSNP="${REF_DIR}/dbsnp/dbsnp_138.hg19.nochrM.fixhdr.vcf.gz"

# Leave empty if running without target restriction / target metrics
TARGETS="${REF_DIR}/design.hg19.bed"

WORK="${PROJECT_DIR}/work"

THREADS=16
FASTP_THREADS=4
PAIRHMM_THREADS=4
