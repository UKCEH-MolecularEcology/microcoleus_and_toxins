#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

AA_FASTA="Dataset2.AA.microcoleus.and.GenBank.fa"
ASSEMBLY_LINK="assembly"
OUTDIR="abundance_out_from_bams"

TOTAL_THREADS=64
INDEX_JOBS=16
INDEX_THREADS=4
THREADS_PER_SAMPLE=8
JOBS=$((TOTAL_THREADS / THREADS_PER_SAMPLE))

EVALUE="1e-10"
DIAMOND_QUERY_COVER=50

DBDIR="${OUTDIR}/diamond_db"
TMPDIR="${OUTDIR}/tmp"
PERSAMPLE="${OUTDIR}/per_sample"
MATRIXDIR="${OUTDIR}/matrices"
LOGFILE="${OUTDIR}/pipeline.log"

mkdir -p "$DBDIR" "$TMPDIR" "$PERSAMPLE" "$MATRIXDIR"

log(){
  echo "[$(date '+%F %T')] $*" | tee -a "$LOGFILE"
}

need_cmd(){ command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing '$1'"; exit 1; }; }

need_cmd diamond
need_cmd samtools
need_cmd parallel
need_cmd awk
need_cmd readlink

log "========== MICROCOLEUS ABUNDANCE PIPELINE =========="

############################################################
# Resolve assembly symlink
############################################################
ASSEMBLY_DIR="$(readlink -f "$ASSEMBLY_LINK")"
log "Assembly directory resolved to: $ASSEMBLY_DIR"

mapfile -t FASTAS < <(find "$ASSEMBLY_DIR" -maxdepth 1 -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \) | sort)
TOTAL_SAMPLES=${#FASTAS[@]}

if [[ "$TOTAL_SAMPLES" -eq 0 ]]; then
  log "ERROR: No assembly FASTA files found."
  exit 1
fi

log "Detected $TOTAL_SAMPLES assemblies."

############################################################
# Genome IDs
############################################################
GENOME_LIST="${TMPDIR}/genome_ids.txt"
grep '^>' "$AA_FASTA" | sed 's/^>//' | awk '{print $1}' | sort -u > "$GENOME_LIST"
GENOME_COUNT=$(wc -l < "$GENOME_LIST")
log "Genome markers detected: $GENOME_COUNT"

############################################################
# Build DIAMOND DB (if needed)
############################################################
DB_PREFIX="${DBDIR}/Dataset2_AA"
if [[ ! -f "${DB_PREFIX}.dmnd" ]]; then
  log "Building DIAMOND database..."
  diamond makedb --in "$AA_FASTA" -d "$DB_PREFIX"
else
  log "DIAMOND database already exists."
fi

############################################################
# Index BAMs
############################################################
log "Indexing BAMs (parallel: ${INDEX_JOBS} jobs x ${INDEX_THREADS} threads)..."

find "$ASSEMBLY_DIR" -maxdepth 1 -type f -name "*.bam" | sort | \
parallel --bar -j "${INDEX_JOBS}" '
  bam="{}"
  if [[ -f "${bam}.bai" || -f "${bam%.bam}.bai" || -f "${bam}.csi" ]]; then
    exit 0
  fi
  samtools index -@ '"${INDEX_THREADS}"' "$bam"
'

log "BAM indexing complete."

############################################################
# Worker script
############################################################
WORKER="${TMPDIR}/worker.sh"

cat > "$WORKER" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

fasta="$1"
ASSEMBLY_DIR="$2"
DB_PREFIX="$3"
GENOME_LIST="$4"
TMPDIR="$5"
PERSAMPLE="$6"
EVALUE="$7"
DIAMOND_QUERY_COVER="$8"
THREADS_PER_SAMPLE="$9"

sample="$(basename "$fasta")"
sample="${sample%.*}"
bam="${ASSEMBLY_DIR}/${sample}.bam"

echo "START $sample"

if [[ ! -f "$bam" ]]; then
  echo "SKIP $sample (missing BAM)"
  exit 0
fi

# Ensure index
if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" && ! -f "${bam}.csi" ]]; then
  samtools index -@ 2 "$bam"
fi

m8="${TMPDIR}/${sample}.diamond.m8"
c2g="${TMPDIR}/${sample}.contig2genome.tsv"
out="${PERSAMPLE}/${sample}.genome_abundance.tsv"

if [[ ! -s "$c2g" ]]; then
  diamond blastx \
    -d "$DB_PREFIX" \
    -q "$fasta" \
    -o "$m8" \
    --outfmt 6 qseqid sseqid bitscore \
    --max-target-seqs 1 \
    --evalue "$EVALUE" \
    --query-cover "$DIAMOND_QUERY_COVER" \
    --threads "$THREADS_PER_SAMPLE"

  awk 'BEGIN{OFS="\t"}{
      if(!(best[$1]) || $3>best[$1]){
        best[$1]=$3; hit[$1]=$2
      }
    }
    END{for(c in hit) print c, hit[c]}' "$m8" > "$c2g"
fi

if [[ ! -s "$c2g" ]]; then
  echo "DONE $sample (no hits)"
  exit 0
fi

samtools coverage -H "$bam" > "${TMPDIR}/${sample}.cov"

echo "DONE $sample"
EOF

chmod +x "$WORKER"

############################################################
# Run per-sample in parallel with progress bar
############################################################
log "Starting per-sample processing..."
printf "%s\n" "${FASTAS[@]}" | \
parallel --bar -j "${JOBS}" \
"$WORKER" {} "$ASSEMBLY_DIR" "$DB_PREFIX" "$GENOME_LIST" "$TMPDIR" "$PERSAMPLE" \
"$EVALUE" "$DIAMOND_QUERY_COVER" "$THREADS_PER_SAMPLE"

log "Per-sample processing complete."

############################################################
# Completion summary
############################################################
COMPLETED=$(ls -1 "$PERSAMPLE"/*.tsv 2>/dev/null | wc -l || true)

log "========== SUMMARY =========="
log "Total assemblies: $TOTAL_SAMPLES"
log "Completed outputs: $COMPLETED"
log "Output directory: $PERSAMPLE"
log "Pipeline finished successfully."


