#!/bin/bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIM_BASE="$SCRIPT_DIR/../simulations/cutoff2.5_standard_LJ"

ROOT_DIRS=(
	"$SIM_BASE/eps0.50_co_phase_separated"
	"$SIM_BASE/eps0.25_distinct_condensates"
	"$SIM_BASE/eps0.10_distinct_condensates"
	"$SIM_BASE/eps0.40_touching_condensates"
)

files=(
)


n_ok=0
n_missing=0
n_failed=0
failed_list=()

process_dir() {
    local target="$1"
    if [ ! -d "$target" ]; then
        echo "  Skipping: $target (not a directory)"
        return
    fi

    echo "Entering: $target"
    pushd "$target" >/dev/null || { echo "  [FAIL] could not cd into $target"; return; }

    for f in "${files[@]}"; do
        if [ -f "$f" ]; then
            echo "  -> Executing $f"
            if jupyter nbconvert \
                --to notebook \
                --execute "$f" \
                --inplace \
                --ExecutePreprocessor.timeout=6000; then
                echo "     [OK] $f"
                n_ok=$((n_ok + 1))
            else
                echo "     [FAIL] $f (see nbconvert output above for the error)"
                n_failed=$((n_failed + 1))
                failed_list+=("$target/$f")
            fi
        else
            echo "  -> $f not found in $target"
            n_missing=$((n_missing + 1))
        fi
    done

    popd >/dev/null
}

for root in "${ROOT_DIRS[@]}"; do
    echo "Processing root: $root"

    if [ ! -d "$root" ]; then
        echo "  [WARN] Root directory not found, skipping: $root"
        continue
    fi

    shopt -s nullglob
    n_dirs=("$root"/N_*/)
    shopt -u nullglob

    if [ ${#n_dirs[@]} -eq 0 ]; then
        echo "  [WARN] No N_* folders found under $root"
        continue
    fi

    for n_dir in "${n_dirs[@]}"; do
        process_dir "$n_dir"
    done
done

echo ""
echo "=== Summary ==="
echo "  Notebooks executed successfully : $n_ok"
echo "  Notebooks missing               : $n_missing"
echo "  Notebooks failed                : $n_failed"

if [ ${#failed_list[@]} -gt 0 ]; then
    echo ""
    echo "Failed notebooks:"
    for f in "${failed_list[@]}"; do
        echo "  - $f"
    done
fi

echo ""
echo "Done! All specified directories processed."
