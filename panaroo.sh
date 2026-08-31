#!/bin/bash
set -euo pipefail

BASE="/mnt/d/ESKAPE_Prokka_Results/Panaroo_Input"
OUT_BASE="/mnt/d/ESKAPE_Prokka_Results/Panaroo_Results"

mkdir -p "$OUT_BASE"

PATHOGENS=(KP AB PA SA EF Ent)

for P in "${PATHOGENS[@]}"; do
    INPUT_DIR="${BASE}/${P}_Panaroo_Input"
    OUTPUT_DIR="${OUT_BASE}/${P}_Panaroo"
    
    echo ""
    echo "========================================"
    echo "  Running Panaroo: ${P}"
    echo "  Input : ${INPUT_DIR}"
    echo "  Output: ${OUTPUT_DIR}"
    echo "  GFF count: $(ls ${INPUT_DIR}/*.gff 2>/dev/null | wc -l)"
    echo "========================================"
    
    panaroo \
        -i "${INPUT_DIR}"/*.gff \
        -o "${OUTPUT_DIR}" \
        -c 0.95 \
        -f 0.80 \
        --clean-mode strict \
        --remove-invalid-genes \
        --refind-mode strict \
        --core_threshold 1.0 \
        --merge_paralogs \
        -a core \
        --aligner mafft \
        -t 56 \
        2>&1 | tee "${OUT_BASE}/${P}_panaroo.log"
    
    echo "[DONE] ${P} — $(date)"
done

echo ""
echo "All 6 ESKAPE Panaroo runs complete."