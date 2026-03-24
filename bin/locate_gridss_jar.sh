#!/bin/bash
set -euo pipefail

# --------------------------------------------------------------------------
# Locate the GRIDSS / Picard JAR.
#
# Strategy (first match wins):
#   1. Look next to the resolved `gridss` wrapper script (works in containers)
#   2. Search $CONDA_PREFIX/share/gridss-* (works in conda environments)
#   3. Search $CONDA_PREFIX/share/picard-* as a last resort
# --------------------------------------------------------------------------

find_jar_in_dir() {
    # Usage: find_jar_in_dir <directory>
    # Prints the first matching JAR (sorted for determinism) or nothing.
    local dir="$1"
    if [ -d "$dir" ]; then
        find "$dir" -maxdepth 1 \( -name "gridss*.jar" -o -name "picard*.jar" \) 2>/dev/null \
            | sort | head -n 1
    fi
}

# --- Strategy 1: directory of the resolved `gridss` wrapper ---
gridssShellPath=$(which gridss 2>/dev/null || true)

if [ -n "$gridssShellPath" ]; then
    SOURCE="$gridssShellPath"
    while [ -h "$SOURCE" ]; do
        DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
        SOURCE="$(readlink "$SOURCE")"
        [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
    done
    DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

    jar=$(find_jar_in_dir "$DIR")
    if [ -n "$jar" ]; then
        echo "$jar"
        exit 0
    fi
fi

# --- Strategy 2: conda share directories (gridss, then picard) ---
if [ -n "${CONDA_PREFIX:-}" ]; then
    for pattern in "$CONDA_PREFIX"/share/gridss-* "$CONDA_PREFIX"/share/picard-*; do
        jar=$(find_jar_in_dir "$pattern")
        if [ -n "$jar" ]; then
            echo "$jar"
            exit 0
        fi
    done
fi

echo "ERROR: could not locate gridss or picard JAR" >&2
exit 1
