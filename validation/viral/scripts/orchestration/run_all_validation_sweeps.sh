#!/bin/bash
set -euo pipefail

# Master script to run all validation sweeps for publication
# Run this on the remote host with sufficient GPU memory

cd ~/RNA-Validation/Biosphere_codec
source venv/bin/activate

DATA_DIR="$HOME/RNA-Validation/Biosphere_codec/data/validation_datasets"
SWEEP_SCRIPT="hex_only_sweep_wrapper.py"

# Grid and seeds
KAPPA_GRID="1.20,1.25,1.30,1.35,1.40,1.45,1.50,1.55,1.60"
SEEDS="13,17,23"

echo "========================================="
echo "VALIDATION SWEEP MASTER RUN"
echo "========================================="
echo "Start: $(date)"
echo

# Function to run a sweep
run_sweep() {
    local NAME=$1
    local FASTA=$2
    local LOG_DIR=$3
    local MIN_LEN=${4:-1000}
    local MAX_LEN=${5:-50000}
    
    if [ ! -f "$FASTA" ]; then
        echo "⚠️  $NAME: FASTA not found: $FASTA"
        return 1
    fi
    
    local N_SEQS=$(grep -c '^>' "$FASTA" || echo 0)
    if [ "$N_SEQS" -lt 50 ]; then
        echo "⚠️  $NAME: Too few sequences ($N_SEQS), skipping"
        return 1
    fi
    
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "SWEEP: $NAME ($N_SEQS sequences)"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    python3 "$SWEEP_SCRIPT" \
        --data_path "$FASTA" \
        --log_dir "$LOG_DIR" \
        --kappa_grid "$KAPPA_GRID" \
        --seeds "$SEEDS" \
        --min_length "$MIN_LEN" \
        --max_length "$MAX_LEN" \
        --max_sequences 10000 \
        --hex_weight 0.8 \
        --batch_size 4 \
        --steps_per_kappa 3000 \
        2>&1 | tee "${LOG_DIR}_sweep.log"
    
    echo "✓ $NAME sweep complete"
    echo
}

# Flaviviridae (RNA, expect κ ~ 1.35-1.55)
run_sweep "Zika" "$DATA_DIR/zika_multi.fasta" "logs_sweep_zika" 9000 13000
run_sweep "WestNile" "$DATA_DIR/west_nile_multi.fasta" "logs_sweep_westnile" 9000 13000
run_sweep "YellowFever" "$DATA_DIR/yellow_fever_multi.fasta" "logs_sweep_yellowfever" 9000 13000

# Picornaviridae (RNA, expect κ ~ 1.30-1.45)
run_sweep "Poliovirus" "$DATA_DIR/poliovirus_multi.fasta" "logs_sweep_poliovirus" 7000 8000
run_sweep "Enterovirus" "$DATA_DIR/enterovirus_multi.fasta" "logs_sweep_enterovirus" 7000 8000

# Paramyxoviridae (RNA, expect κ ~ 1.30-1.45)
run_sweep "Measles" "$DATA_DIR/measles_multi.fasta" "logs_sweep_measles" 15000 17000
run_sweep "Mumps" "$DATA_DIR/mumps_multi.fasta" "logs_sweep_mumps" 15000 17000

# Rhabdoviridae (RNA, expect κ ~ 1.30-1.40)
run_sweep "Rabies" "$DATA_DIR/rabies_multi.fasta" "logs_sweep_rabies" 11000 13000

# Filoviridae (RNA, expect κ ~ 1.35-1.50)
run_sweep "Ebola" "$DATA_DIR/ebola_multi.fasta" "logs_sweep_ebola" 18000 20000

# DNA controls (expect κ ~ 1.20-1.28)
run_sweep "HSV1" "$DATA_DIR/hsv1_multi.fasta" "logs_sweep_hsv1_DNA" 100000 200000
run_sweep "CMV" "$DATA_DIR/cmv_multi.fasta" "logs_sweep_cmv_DNA" 200000 250000

echo "========================================="
echo "ALL SWEEPS COMPLETE"
echo "End: $(date)"
echo "========================================="

# Generate summary
python3 summarize_sweeps.py \
    logs_sweep_zika \
    logs_sweep_westnile \
    logs_sweep_yellowfever \
    logs_sweep_poliovirus \
    logs_sweep_enterovirus \
    logs_sweep_measles \
    logs_sweep_mumps \
    logs_sweep_rabies \
    logs_sweep_ebola \
    logs_sweep_hsv1_DNA \
    logs_sweep_cmv_DNA \
    --output publication_outputs/complete_validation_summary.md

echo "✓ Summary generated: publication_outputs/complete_validation_summary.md"
