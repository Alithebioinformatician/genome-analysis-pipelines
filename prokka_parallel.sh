#!/bin/bash
# =============================================================================
# ESKAPE Prokka Annotation Pipeline — Parallel Edition
# =============================================================================
# Pathogens  : AB | KP | SA | PA | EF | Ent
# Parallelism: GNU parallel — 8 concurrent jobs x 10 CPUs/job
# ETA        : Live on terminal via --eta | per-folder timing in master log
# Fixes vs v1: locus-tag collision, counter thread-safety, --addgenes removed
# =============================================================================
#
# USAGE
#   conda activate prokka_env
#   bash eskape_prokka_parallel.sh
#
# REQUIREMENTS
#   prokka      in PATH
#   GNU parallel in PATH   (parallel --version >= 2016)
#   ~32 GB RAM free        (8 jobs × ~4 GB peak each)
#   ~80 CPU threads        (8 jobs × 10 threads each)
# =============================================================================

# --------------------------------------------------------------------------- #
#  STRICT MODE                                                                 #
#  -u  : unbound variable reference → immediate fatal error                   #
#  -o pipefail : pipeline exit = last non-zero stage exit                     #
#  NO -e : we handle per-job failures ourselves; -e would abort on any        #
#          non-zero inside the loop                                            #
# --------------------------------------------------------------------------- #
set -uo pipefail

# --------------------------------------------------------------------------- #
#  PATHS  —  edit these                                                        #
# --------------------------------------------------------------------------- #
readonly BASE_DIR="/home/bioinfo/Desktop/Data_5TB/BSI"
readonly OUT_BASE="/home/bioinfo/Desktop/Data_5TB/BSIBSI_Prokka_Results"
readonly LOG_DIR="$OUT_BASE/logs"

# --------------------------------------------------------------------------- #
#  PARALLELISM CONFIG  —  edit PARALLEL_JOBS if you want fewer/more slots     #
# --------------------------------------------------------------------------- #
PARALLEL_JOBS=6
THREADS_PER_JOB=12                                   # 6 × 12 = 72 cores
TOTAL_CPUS=$(nproc)

# --------------------------------------------------------------------------- #
#  TIMESTAMP & MASTER LOG                                                      #
# --------------------------------------------------------------------------- #
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
mkdir -p "$LOG_DIR"
readonly MASTER_LOG="$LOG_DIR/pipeline_${TIMESTAMP}.log"

# --------------------------------------------------------------------------- #
#  ATOMIC COUNTER DIRECTORY                                                    #
#  Each worker appends one genome name per event to ok / skip / fail files.   #
#  Appends under 4 096 bytes are atomic on Linux — no flock required.         #
#  Directory is deleted on any exit (normal, error, or signal).               #
# --------------------------------------------------------------------------- #
COUNTER_DIR=$(mktemp -d)
WORKER_SCRIPT=$(mktemp /tmp/prokka_worker_XXXXXX.sh)
trap '[[ -d "${COUNTER_DIR:-}"  ]] && rm -rf "$COUNTER_DIR";
      [[ -f "${WORKER_SCRIPT:-}" ]] && rm -f  "$WORKER_SCRIPT"' EXIT

# --------------------------------------------------------------------------- #
#  SUPPRESS GNU PARALLEL CITATION NAG                                          #
#  Creates the sentinel file all parallel versions look for.                  #
# --------------------------------------------------------------------------- #
mkdir -p "$HOME/.parallel"
touch    "$HOME/.parallel/will-cite"

# --------------------------------------------------------------------------- #
#  LOGGING HELPER  —  MAIN PROCESS ONLY                                       #
#  Workers never call this; they write to their own per-job log files.        #
#  This keeps all master-log writes single-threaded → no interleaving.        #
# --------------------------------------------------------------------------- #
log() {
    local level="$1"; shift
    printf '[%s] [%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$level" "$*" \
        | tee -a "$MASTER_LOG"
}

# --------------------------------------------------------------------------- #
#  DURATION FORMATTER  —  seconds → human-readable string                     #
# --------------------------------------------------------------------------- #
format_duration() {
    local secs=$1
    local h=$(( secs / 3600 ))
    local m=$(( (secs % 3600) / 60 ))
    local s=$(( secs % 60 ))
    if   [[ $h -gt 0 ]]; then printf '%dh %02dm %02ds' "$h" "$m" "$s"
    elif [[ $m -gt 0 ]]; then printf '%dm %02ds'        "$m" "$s"
    else                       printf '%ds'              "$s"
    fi
}

# --------------------------------------------------------------------------- #
#  PER-PATHOGEN CONFIGURATION                                                  #
#  NOTE: Associative arrays cannot be bash-exported to subprocesses.          #
#  They are used only by the MAIN PROCESS to look up values, which are then   #
#  passed as explicit positional arguments to the worker function.             #
# --------------------------------------------------------------------------- #
declare -A GENUS=(
    [AB]="Acinetobacter"   [KP]="Klebsiella"
    [SA]="Staphylococcus"  [PA]="Pseudomonas"
    [EF]="Enterococcus"    [Ent]="Enterobacter"
)
declare -A SPECIES=(
    [AB]="baumannii"   [KP]="pneumoniae"
    [SA]="aureus"      [PA]="aeruginosa"
    [EF]="faecium"     [Ent]=""          # empty = genus-only (12 mixed spp.)
)
declare -A GRAM=(
    [AB]="neg"   [KP]="neg"
    [SA]="pos"   [PA]="neg"
    [EF]="pos"   [Ent]="neg"
)
declare -A LOCUS_PREFIX=(
    [AB]="ABSP"    [KP]="KPSP"
    [SA]="SASP"    [PA]="PASP"
    [EF]="EFSP"    [Ent]="ENTSP"
)

FOLDERS=("AB" "KP" "SA" "PA" "EF" "Ent")

# If a folder key is passed as argument, process only that species
# Usage:  bash eskape_prokka_parallel.sh KP
#         bash eskape_prokka_parallel.sh       ← runs all six
if [[ $# -gt 0 ]]; then
    FOLDERS=("$1")
    log_override=1   # flag used in header below
fi

# --------------------------------------------------------------------------- #
#  OUTPUT FILES TO KEEP AFTER ANNOTATION                                      #
#  (hardcoded inside worker too — cannot rely on exported arrays)             #
#   .gff  → Roary / Panaroo pan-genome                                        #
#   .faa  → AMRFinderPlus / ResFinder protein search                          #
#   .fna  → MLST / reference mapping                                          #
#   .ffn  → core-SNP phylogenomics (nucleotide CDSs)                          #
#   .gbk  → full feature archive, manual inspection                           #
#   .txt  → QC summary (CDS count, rRNA, tRNA, hypothetical %)               #
#   .tsv  → tab-delimited annotation table                                    #
# --------------------------------------------------------------------------- #
KEEP=("gbk" "gff" "fna" "ffn" "faa")

# --------------------------------------------------------------------------- #
#  WORKER FUNCTION                                                             #
#  Self-contained — zero dependency on associative arrays or global vars.     #
#  Everything it needs is passed as positional arguments.                     #
#                                                                              #
#  Args                                                                        #
#    $1  folder        pathogen key               e.g.  AB                    #
#    $2  file          full path to .fna assembly                              #
#    $3  idx           sequential genome index (1..N), used for locus tag     #
#    $4  threads       CPUs assigned to this Prokka instance                  #
#    $5  genus         Prokka --genus value                                    #
#    $6  species       Prokka --species value  (empty string = genus-only)    #
#    $7  gram          pos | neg                                               #
#    $8  locus_prefix  short alphanumeric prefix  e.g. ABSP                   #
#    $9  out_base      root output directory                                   #
#   $10  log_dir       directory for per-job log files                         #
#   $11  counter_dir   temp dir for atomic ok / skip / fail counters          #
# --------------------------------------------------------------------------- #
run_prokka_parallel() {
    local folder="$1"
    local file="$2"
    local idx="$3"
    local threads="$4"
    local genus="$5"
    local species="$6"
    local gram="$7"
    local locus_prefix="$8"
    local out_base="$9"
    local log_dir="${10}"
    local counter_dir="${11}"

    local base
    base=$(basename "$file" .fna)

    local outdir="$out_base/$folder/$base"
    local joblog="$log_dir/${folder}_${base}.log"

    # ------------------------------------------------------------------ #
    #  Per-job log helper                                                  #
    #  Writes ONLY to this genome's private log file.                     #
    #  Never touches the master log — all master-log writes are in main.  #
    #  $joblog is in scope because bash uses dynamic (call-stack) scoping  #
    #  for locals, so nested functions see the caller's locals.            #
    # ------------------------------------------------------------------ #
    jlog() {
        printf '[%s] [%s] %s\n' \
            "$(date '+%Y-%m-%d %H:%M:%S')" "$1" "${*:2}" \
            >> "$joblog"
    }

    # ------------------------------------------------------------------ #
    #  Skip if already annotated  (idempotent — safe to re-run pipeline)  #
    # ------------------------------------------------------------------ #
    if [[ -f "$outdir/$base.gff" ]] && grep -q "^##FASTA" "$outdir/$base.gff"; then
    jlog "SKIP" "$folder/$base — complete GFF exists, skipping"
    echo "$base" >> "$counter_dir/$folder/skip"
    return 0
    elif [[ -f "$outdir/$base.gff" ]]; then
	     jlog "WARN " "$folder/$base — incomplete GFF found (no ##FASTA), re-annotating"
         rm -f "$outdir/$base.gff"   # let --force handle the rest
    fi

    mkdir -p "$outdir"

    # ------------------------------------------------------------------ #
    #  Locus tag:  PREFIX + zero-padded sequential index                  #
    #                                                                      #
    #  Format:  ABSP0001 … ABSP9999   (8 chars max for ABSP prefix)      #
    #           ENTSP0001 … ENTSP9999  (9 chars max for ENTSP prefix)     #
    #                                                                      #
    #  Prokka appends _00001, _00002 … to this prefix for individual      #
    #  gene locus tags.  The prefix itself must be unique per genome.     #
    #                                                                      #
    #  This replaces the old approach of truncating PREFIX_ACCESSION to   #
    #  20 chars, which silently created identical prefixes across genomes  #
    #  when accession strings shared a long common substring, corrupting   #
    #  Panaroo / Roary input downstream.                                  #
    # ------------------------------------------------------------------ #
    local locus_tag
    locus_tag=$(printf '%s%04d' "$locus_prefix" "$idx")

    # ------------------------------------------------------------------ #
    #  Build Prokka command                                                #
    #                                                                      #
    #  --addgenes removed: deprecated in Prokka >= 1.14; behavior is now  #
    #  the default or handled via --compliant. Leaving it in generates a  #
    #  warning that contaminates logs.                                     #
    #                                                                      #
    #  --force: overwrite any partial output from a previous crashed run  #
    #  on this genome. Safe because skip logic above guards completed      #
    #  runs (GFF present) — --force only fires on partial outputs.        #
    # ------------------------------------------------------------------ #
    local cmd=(
        prokka
        --outdir        "$outdir"
        --prefix        "$base"
        --kingdom       Bacteria
        --genus         "$genus"
        --usegenus
        --locustag      "$locus_tag"
        --mincontiglen  200
        --cpus          "$threads"
        --force
        # --gram removed: requires signalp (proprietary, not installed)
        # Signal peptide prediction is irrelevant for pan-genome/AMR/phylogenomics
    )

    # Add --species only when defined (skipped for Ent: mixed-species folder)
    [[ -n "$species" ]] && cmd+=(--species "$species")

    # Input file — always last
    cmd+=("$file")

    # ------------------------------------------------------------------ #
    #  Execute Prokka                                                      #
    # ------------------------------------------------------------------ #
    jlog "START" "$folder/$base | locus_tag=$locus_tag | cpus=$threads"

    if "${cmd[@]}" >> "$joblog" 2>&1; then

        jlog "OK   " "$folder/$base — annotation complete"
        echo "$base" >> "$counter_dir/$folder/ok"

        # Keep only what downstream tools need:
        #   .gff  → Panaroo / Roary pan-genome
        #   .faa  → AMRFinderPlus / ResFinder protein search
        #   .fna  → MLST / reference mapping
        #   .ffn  → core-SNP phylogenomics (nucleotide CDSs)
        #   .gbk  → full feature archive, manual inspection
        local keep_exts=("gbk" "gff" "fna" "ffn" "faa")
        local find_cmd=(find "$outdir" -type f)
        for ext in "${keep_exts[@]}"; do
            find_cmd+=(! -name "$base.$ext")
        done
        find_cmd+=(-delete)
        "${find_cmd[@]}"

    else
        jlog "FAIL " "$folder/$base — Prokka failed (full output above)"
        echo "$base" >> "$counter_dir/$folder/fail"
        return 1    # non-zero exit is recorded in --joblog TSV
    fi
}

# --------------------------------------------------------------------------- #
#  WORKER TEMP SCRIPT                                                         #
#  GNU parallel 20160622 does NOT pick up bash-exported functions via         #
#  `export -f`. It spawns /usr/bin/bash in a way that ignores bash's         #
#  BASH_FUNC_ environment variables entirely — hence the error:              #
#    /usr/bin/bash: line 1: run_prokka_parallel: command not found           #
#                                                                             #
#  Fix: write the function definition to a temp .sh file and call            #
#  `bash "$WORKER_SCRIPT"` directly. No environment variable mechanism       #
#  needed — works on every bash/parallel version combination.               #
# --------------------------------------------------------------------------- #
{
    declare -f run_prokka_parallel   # emit the full function definition
    echo 'run_prokka_parallel "$@"'  # append: call it with all positional args
} > "$WORKER_SCRIPT"

# --------------------------------------------------------------------------- #
#  SANITY CHECKS                                                               #
# --------------------------------------------------------------------------- #
for tool in prokka parallel paste seq; do
    if ! command -v "$tool" &>/dev/null; then
        echo "[ERROR] Required tool not found in PATH: $tool" >&2
        case "$tool" in
            prokka)   echo "  → conda activate prokka_env" >&2 ;;
            parallel) echo "  → conda install -c conda-forge parallel" >&2 ;;
        esac
        exit 1
    fi
done

# --------------------------------------------------------------------------- #
#  PER-FOLDER RESULT TRACKERS  (read from counter files after each run)       #
# --------------------------------------------------------------------------- #
declare -A FOLDER_TOTAL FOLDER_OK FOLDER_SKIP FOLDER_FAIL FOLDER_ELAPSED

# --------------------------------------------------------------------------- #
#  PIPELINE HEADER                                                             #
# --------------------------------------------------------------------------- #
log "INFO " "================================================================"
log "INFO " "  ESKAPE Prokka Pipeline (Parallel Edition) — $(date)"
log "INFO " "================================================================"
log "INFO " "  Parallel jobs   : $PARALLEL_JOBS"
log "INFO " "  Threads / job   : $THREADS_PER_JOB   (total CPUs: $TOTAL_CPUS)"
log "INFO " "  Base dir        : $BASE_DIR"
log "INFO " "  Output dir      : $OUT_BASE"
log "INFO " "  Master log      : $MASTER_LOG"
log "INFO " "  ETA display     : live on terminal (stderr) via --eta"
log "INFO " "================================================================"

# --------------------------------------------------------------------------- #
#  MAIN LOOP                                                                   #
#  Pathogens run sequentially.                                                 #
#  Genomes within each pathogen run in parallel (PARALLEL_JOBS slots).        #
# --------------------------------------------------------------------------- #
for folder in "${FOLDERS[@]}"; do

    input_dir="$BASE_DIR/$folder"

    # ------------------------------------------------------------------ #
    #  Validate input directory                                            #
    # ------------------------------------------------------------------ #
    if [[ ! -d "$input_dir" ]]; then
        log "WARN " "Directory not found: $input_dir — skipping $folder entirely"
        continue
    fi

    # Collect .fna files into a sorted array
    mapfile -t fna_files < <(find "$input_dir" -maxdepth 1 -name "*.fna" | sort)

    if [[ ${#fna_files[@]} -eq 0 ]]; then
        log "WARN " "No .fna files found in $input_dir — skipping $folder"
        continue
    fi

    total=${#fna_files[@]}
    FOLDER_TOTAL[$folder]=$total

    # Create counter subdirectory and output directory for this folder
    mkdir -p "$COUNTER_DIR/$folder" "$OUT_BASE/$folder"

    # Resolve species string for log display
    sp_display="${SPECIES[$folder]:-}"
    [[ -z "$sp_display" ]] && sp_display="(genus-only — mixed species)"

    log "INFO " "----------------------------------------------------------------"
    log "INFO " "  Pathogen   : $folder"
    log "INFO " "  Organism   : ${GENUS[$folder]} $sp_display"
    log "INFO " "  Gram       : ${GRAM[$folder]}"
    log "INFO " "  Locus pfx  : ${LOCUS_PREFIX[$folder]}  →  tags: ${LOCUS_PREFIX[$folder]}0001 … ${LOCUS_PREFIX[$folder]}$(printf '%04d' "$total")"
    log "INFO " "  Genomes    : $total"
    log "INFO " "  Concurrency: $PARALLEL_JOBS jobs × $THREADS_PER_JOB CPUs each"
    log "INFO " "----------------------------------------------------------------"
    log "INFO " "  Launching parallel workers — watch terminal for live ETA ..."

    folder_start=$(date +%s)

    # ------------------------------------------------------------------ #
    #  PARALLEL EXECUTION BLOCK                                            #
    #                                                                      #
    #  Input construction                                                  #
    #    seq 1 N          → one integer per line  (genome index)          #
    #    printf '%s\n'    → one file path per line                        #
    #    paste            → merges into: idx <TAB> filepath               #
    #                                                                      #
    #  GNU parallel reads each line, splits on TAB:                       #
    #    {1} = index                                                       #
    #    {2} = file path                                                   #
    #                                                                      #
    #  Flags                                                               #
    #    --jobs N         : N concurrent worker slots                     #
    #    --eta            : live ETA + progress bar on terminal (stderr)  #
    #    --colsep $'\t'   : split input lines on tab                      #
    #    --joblog FILE    : TSV log with per-job runtime + exit code      #
    #    --halt never     : never abort on individual job failure         #
    #                                                                      #
    #  All per-pathogen config (genus, species, gram, prefix) is passed   #
    #  as fixed arguments — avoids the bash associative-array export      #
    #  limitation (declare -A arrays are NOT exportable via `export`).    #
    # ------------------------------------------------------------------ #
    paste \
        <(seq 1 "$total") \
        <(printf '%s\n' "${fna_files[@]}") \
    | parallel \
        --jobs        "$PARALLEL_JOBS"                       \
        --eta                                                \
        --colsep      $'\t'                                  \
        --joblog      "$LOG_DIR/${folder}_parallel.joblog"   \
        --halt        never                                  \
        bash "$WORKER_SCRIPT"                                \
            "$folder"                                        \
            "{2}"                                            \
            "{1}"                                            \
            "$THREADS_PER_JOB"                               \
            "${GENUS[$folder]}"                              \
            "${SPECIES[$folder]}"                            \
            "${GRAM[$folder]}"                               \
            "${LOCUS_PREFIX[$folder]}"                       \
            "$OUT_BASE"                                      \
            "$LOG_DIR"                                       \
            "$COUNTER_DIR"

    # Capture parallel's exit code before anything else resets PIPESTATUS
    # With --halt never: parallel exits with number-of-failed-jobs (0 = all OK)
    # With no -e in our strict mode: non-zero here does NOT abort the script
    par_rc=${PIPESTATUS[1]}

    folder_end=$(date +%s)
    elapsed=$(( folder_end - folder_start ))
    FOLDER_ELAPSED[$folder]=$elapsed

    # ------------------------------------------------------------------ #
    #  Read atomic counters                                                #
    #  Each file has one genome name per line.                            #
    #  wc -l gives the exact genome count for each outcome.              #
    # ------------------------------------------------------------------ #
    ok=0
    skip=0
    fail=0
    [[ -f "$COUNTER_DIR/$folder/ok"   ]] && ok=$(wc   -l < "$COUNTER_DIR/$folder/ok")
    [[ -f "$COUNTER_DIR/$folder/skip" ]] && skip=$(wc -l < "$COUNTER_DIR/$folder/skip")
    [[ -f "$COUNTER_DIR/$folder/fail" ]] && fail=$(wc -l < "$COUNTER_DIR/$folder/fail")

    FOLDER_OK[$folder]=$ok
    FOLDER_SKIP[$folder]=$skip
    FOLDER_FAIL[$folder]=$fail

    # ------------------------------------------------------------------ #
    #  Throughput calculation                                              #
    # ------------------------------------------------------------------ #
    dur_str=$(format_duration "$elapsed")
    if [[ "$elapsed" -gt 0 ]]; then
        throughput=$(awk "BEGIN { printf \"%.1f\", $total / ($elapsed / 60.0) }")
    else
        throughput="<1"
    fi

    # ------------------------------------------------------------------ #
    #  Folder completion log                                               #
    # ------------------------------------------------------------------ #
    log "INFO " "  ✓ $folder complete : $total genomes | $dur_str | ${throughput} genome/min"
    log "INFO " "  $folder counts     : OK=$ok  SKIP=$skip  FAIL=$fail"

    if [[ "$par_rc" -ne 0 ]]; then
        log "WARN " "  parallel exit code = $par_rc for $folder ($fail failed job(s))"
    fi

    # ------------------------------------------------------------------ #
    #  List failed genome names (read from counter file)                  #
    #  Per-job logs are at: $LOG_DIR/${folder}_${genome}.log              #
    # ------------------------------------------------------------------ #
    if [[ "$fail" -gt 0 ]] && [[ -f "$COUNTER_DIR/$folder/fail" ]]; then
        log "WARN " "  Failed genomes in $folder — check individual logs in $LOG_DIR :"
        while IFS= read -r genome_name; do
            log "WARN " "    ✗  $genome_name"
        done < "$COUNTER_DIR/$folder/fail"
    fi

    log "INFO " "----------------------------------------------------------------"

done

# --------------------------------------------------------------------------- #
#  GRAND SUMMARY TABLE                                                         #
# --------------------------------------------------------------------------- #
log "INFO " "================================================================"
log "INFO " "                     PIPELINE SUMMARY"
log "INFO " "================================================================"

# Print table header
{
    printf '%-8s  %-7s  %-7s  %-7s  %-7s  %-16s\n' \
        "FOLDER" "TOTAL" "OK" "SKIP" "FAIL" "ELAPSED"
    printf '%-8s  %-7s  %-7s  %-7s  %-7s  %-16s\n' \
        "------" "-----" "--" "----" "----" "-------"
} | tee -a "$MASTER_LOG"

grand_total=0
grand_ok=0
grand_skip=0
grand_fail=0

for folder in "${FOLDERS[@]}"; do
    # Skip folders that were never processed (missing directories etc.)
    [[ -z "${FOLDER_TOTAL[$folder]:-}" ]] && continue

    t=${FOLDER_TOTAL[$folder]:-0}
    o=${FOLDER_OK[$folder]:-0}
    s=${FOLDER_SKIP[$folder]:-0}
    f=${FOLDER_FAIL[$folder]:-0}
    dur=$(format_duration "${FOLDER_ELAPSED[$folder]:-0}")

    printf '%-8s  %-7s  %-7s  %-7s  %-7s  %-16s\n' \
        "$folder" "$t" "$o" "$s" "$f" "$dur" | tee -a "$MASTER_LOG"

    # (( expr )) returns 1 when result is 0 — || true prevents false abort
    (( grand_total += t )) || true
    (( grand_ok    += o )) || true
    (( grand_skip  += s )) || true
    (( grand_fail  += f )) || true
done

# Print table footer
{
    printf '%-8s  %-7s  %-7s  %-7s  %-7s\n' \
        "------" "-----" "--" "----" "----"
    printf '%-8s  %-7s  %-7s  %-7s  %-7s\n' \
        "TOTAL" "$grand_total" "$grand_ok" "$grand_skip" "$grand_fail"
} | tee -a "$MASTER_LOG"

log "INFO " "================================================================"

if [[ "$grand_fail" -gt 0 ]]; then
    log "WARN " "$grand_fail genome(s) failed across all pathogens."
    log "WARN " "Per-job logs   : $LOG_DIR/<FOLDER>_<genome>.log"
    log "WARN " "Parallel joblogs: $LOG_DIR/<FOLDER>_parallel.joblog"
fi

log "INFO " "Pipeline complete — $(date)"
log "INFO " "================================================================"
