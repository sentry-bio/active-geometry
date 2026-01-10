#!/usr/bin/env bash
# =============================================================================
# Active Geometry: Run All Verifications
# =============================================================================
# One command to verify everything:
#   1. Lean formal proofs (machine-checked)
#   2. Python validation notebooks (numerical)
#   3. Consistency checks (constants.yaml alignment)
#
# Usage:
#   ./run_all_verifications.sh           # Run all
#   SKIP_LEAN=1 ./run_all_verifications.sh  # Skip Lean (faster)
# =============================================================================

set -euo pipefail

echo "=============================================="
echo "  Active Geometry: Verification Suite"
echo "=============================================="
echo ""

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$REPO_ROOT"

PASS=0
FAIL=0
SKIP=0

# Colors (if terminal supports)
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log_pass() { echo -e "${GREEN}PASS${NC} $1"; PASS=$((PASS+1)); }
log_fail() { echo -e "${RED}FAIL${NC} $1"; FAIL=$((FAIL+1)); }
log_skip() { echo -e "${YELLOW}SKIP${NC} $1"; SKIP=$((SKIP+1)); }

# -----------------------------------------------------------------------------
# 1. Lean Formal Proofs
# -----------------------------------------------------------------------------
echo "----------------------------------------------"
echo "1. Verifying Lean formal proofs..."
echo "----------------------------------------------"

if [[ "${SKIP_LEAN:-0}" == "1" ]]; then
    log_skip "Lean proofs (SKIP_LEAN=1)"
else
    if command -v lake &> /dev/null; then
        if [[ -d "$REPO_ROOT/theory/lean" ]]; then
            cd "$REPO_ROOT/theory/lean"
            if lake build 2>&1; then
                log_pass "Lean proofs compiled successfully"
            else
                log_fail "Lean proofs failed to compile"
            fi
            cd "$REPO_ROOT"
        else
            log_skip "Lean directory not found"
        fi
    else
        log_skip "Lean/lake not installed"
    fi
fi
echo ""

# -----------------------------------------------------------------------------
# 2. Python Validation Notebooks
# -----------------------------------------------------------------------------
echo "----------------------------------------------"
echo "2. Executing validation notebooks..."
echo "----------------------------------------------"

NOTEBOOKS_DIR="$REPO_ROOT/validation/notebooks"
EXEC_DIR="$REPO_ROOT/validation/_executed"
mkdir -p "$EXEC_DIR"

if ! command -v jupyter &> /dev/null; then
    echo "Installing Jupyter..."
    pip install -q jupyter nbconvert nbclient nbformat
fi

run_notebook() {
    local nb_path="$1"
    local nb_name="$(basename "$nb_path")"

    if [[ ! -f "$nb_path" ]]; then
        log_skip "$nb_name (not found)"
        return 2
    fi

    echo "  Running: $nb_name"
    local out_path="$EXEC_DIR/$nb_name"

    if jupyter nbconvert --to notebook --execute "$nb_path" \
        --output "$out_path" \
        --ExecutePreprocessor.timeout=600 \
        --ExecutePreprocessor.kernel_name=python3 \
        2>&1 | tail -5; then
        log_pass "$nb_name"
        return 0
    else
        log_fail "$nb_name"
        return 1
    fi
}

# Execute notebooks in order
if [[ -d "$NOTEBOOKS_DIR" ]]; then
    for nb in "$NOTEBOOKS_DIR"/*.ipynb; do
        [[ -e "$nb" ]] || continue
        run_notebook "$nb" || true
    done
else
    log_skip "No notebooks directory found"
fi
echo ""

# -----------------------------------------------------------------------------
# 3. Constants Consistency Check
# -----------------------------------------------------------------------------
echo "----------------------------------------------"
echo "3. Checking constants consistency..."
echo "----------------------------------------------"

if [[ -f "$REPO_ROOT/constants.yaml" ]]; then
    # Quick sanity check: kappa values present and reasonable
    if python3 - <<'PY'
import yaml
from pathlib import Path

constants = yaml.safe_load(Path("constants.yaml").read_text())
kappa_emp = constants["curvature"]["kappa_empirical"]
kappa_theory = constants["curvature"]["kappa_theory"]

# Check values are in expected range
assert 1.0 < kappa_emp < 1.5, f"kappa_empirical={kappa_emp} out of range"
assert 1.0 < kappa_theory < 1.5, f"kappa_theory={kappa_theory} out of range"

# Check agreement
agreement = abs(kappa_emp - kappa_theory) / kappa_emp * 100
assert agreement < 5, f"Theory-empirical agreement {agreement:.1f}% > 5%"

print(f"  kappa_empirical = {kappa_emp}")
print(f"  kappa_theory = {kappa_theory}")
print(f"  Agreement: {agreement:.1f}%")
PY
    then
        log_pass "constants.yaml consistency"
    else
        log_fail "constants.yaml consistency"
    fi
else
    log_skip "constants.yaml not found"
fi
echo ""

# -----------------------------------------------------------------------------
# 4. Generate Verification Report
# -----------------------------------------------------------------------------
echo "----------------------------------------------"
echo "4. Generating verification report..."
echo "----------------------------------------------"

REPORT_PATH="$REPO_ROOT/verification_report.json"

python3 - <<PY
import json
import hashlib
from pathlib import Path
from datetime import datetime

def sha256(path):
    if not path.exists():
        return None
    return hashlib.sha256(path.read_bytes()).hexdigest()

report = {
    "timestamp": datetime.utcnow().isoformat() + "Z",
    "counts": {
        "pass": $PASS,
        "fail": $FAIL,
        "skip": $SKIP
    },
    "files": {
        "manifest_sha256": sha256(Path("manifest.yaml")),
        "constants_sha256": sha256(Path("constants.yaml")),
    },
    "executed_notebooks": []
}

exec_dir = Path("validation/_executed")
if exec_dir.exists():
    for nb in sorted(exec_dir.glob("*.ipynb")):
        report["executed_notebooks"].append({
            "name": nb.name,
            "sha256": sha256(nb)
        })

Path("$REPORT_PATH").write_text(json.dumps(report, indent=2))
print(f"  Report written to: verification_report.json")
PY

echo ""

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo "=============================================="
echo "  SUMMARY"
echo "=============================================="
echo -e "  ${GREEN}PASS${NC}: $PASS"
echo -e "  ${RED}FAIL${NC}: $FAIL"
echo -e "  ${YELLOW}SKIP${NC}: $SKIP"
echo "=============================================="

if [[ "$FAIL" -gt 0 ]]; then
    echo -e "${RED}Some verifications failed.${NC}"
    exit 1
fi

echo -e "${GREEN}All verifications passed.${NC}"
exit 0
