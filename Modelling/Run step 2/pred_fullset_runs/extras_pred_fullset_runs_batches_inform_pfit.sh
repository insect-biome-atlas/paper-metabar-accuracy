#!/bin/zsh
# chmod +x extras_pred_fullset_runs_batches_inform_pfit.sh
# ./extras_pred_fullset_runs_batches_inform_pfit.sh

set -euo pipefail

echo "Starting: run model on precomputed per-species batches - inform pfit"

SCRIPT_DIR="$(cd "$(dirname "${(%):-%x}")" && pwd)"
MODELLING_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$MODELLING_ROOT/config.sh"

tpplc "$MODELS_DIR/predict/optmodel15_Pred_match_inform_pfit.tppl" \
  -m smc-apf \
  --particles 100000 \
  --subsample \
  --subsample-size 1 \
  --resample manual \
  --output "$BIN_DIR/optmodel15_Pred_match_inform_pfit_100000"

MODEL_BIN="$BIN_DIR/optmodel15_Pred_match_inform_pfit_100000"
ORIG_JSON="$DATA_DIR/15_hom_priors.json"

SWEEPS=${SWEEPS:-100}
OUT_ROOT="${OUT_ROOT:-$PRED_OUT_DIR/fullset/inform_pfit}"
BATCH_DIR="${BATCH_DIR:-$FULLSET_BATCH_DIR}"
RUN_DIR="${OUT_ROOT}/runs"
MANIFEST="${OUT_ROOT}/manifest.json"

OVERWRITE=${OVERWRITE:-1}         # 0=skip when output exists; 1=rerun
ONLY_NONEMPTY=${ONLY_NONEMPTY:-1}  # 1=skip all-zero species by inspecting the batch file
START_AT=${START_AT:-1}             # optional, used if INDEX_LIST is empty and you want a range
END_AT=${END_AT:-472}                 # optional, used with START_AT
INDEX_LIST="8, 18, 31, 142, 151, 189, 190, 192, 257, 273, 283, 295, 347, 352, 371, 384, 392"        # optional explicit list "1, 2, 7 9"

mkdir -p "$RUN_DIR" "$OUT_ROOT"

command -v jq >/dev/null 2>&1 || { echo "ERROR: jq is required"; exit 1; }
[[ -x "$MODEL_BIN" ]] || { echo "ERROR: model binary not found/executable: $MODEL_BIN"; exit 1; }
[[ -d "$BATCH_DIR" ]] || { echo "ERROR: batch dir not found: $BATCH_DIR"; exit 1; }

# Gather available batch files
typeset -a BATCH_FILES
BATCH_FILES=(${BATCH_DIR}/pred_in_species_*.json(N))   # (N) = nullglob in zsh
(( ${#BATCH_FILES[@]} > 0 )) || { echo "ERROR: no batch files found in $BATCH_DIR"; exit 1; }

# Helper: extract 1-based species index from filename
get_idx_from_file() {
  local f="$1"
  basename "$f" | sed -E 's/^[^0-9]*([0-9]+)\.json$/\1/'
}

# Build map: idx -> file
typeset -A FILE_BY_IDX
typeset -a ALL_IDXS
for f in "${BATCH_FILES[@]}"; do
  idx="$(get_idx_from_file "$f")" || true
  [[ "$idx" == <-> ]] || { echo "WARN: skip unrecognized batch file: $f" >&2; continue; }
  FILE_BY_IDX[$idx]="$f"
  ALL_IDXS+="$idx"
done

(( ${#ALL_IDXS[@]} > 0 )) || { echo "ERROR: no valid pred_in_species_*.json files"; exit 1; }

# Determine TOTAL_SAMPLES from any one batch; also sanity-check width==1
SAMPLE_PROBE_FILE="${FILE_BY_IDX[${ALL_IDXS[1]}]}"
TOTAL_SAMPLES=$(jq '(.dataset | length)' "$SAMPLE_PROBE_FILE")
WIDTH_OK=$(jq -e '[.dataset[] | length] | all(. == 1)' "$SAMPLE_PROBE_FILE" >/dev/null; print $?)
if (( WIDTH_OK != 0 )); then
  echo "ERROR: batch file $SAMPLE_PROBE_FILE is not samples×1 column"; exit 1;
fi
echo "Detected ${#ALL_IDXS[@]} available batch files; sample rows per batch: ${TOTAL_SAMPLES}"

# --- Decide which indices to run ---
# Normalize INDEX_LIST if provided
_index_norm="$(printf '%s' "${INDEX_LIST}" | tr ',' ' ' | tr -s '[:space:]' ' ' | sed 's/^ //; s/ $//')"
typeset -a CANDIDATE_IDXS
if [[ -n "$_index_norm" ]]; then
  echo "Mode: explicit vector (INDEX_LIST)"
  CANDIDATE_IDXS=(${=_index_norm})
elif [[ -n "${START_AT}" && -n "${END_AT}" ]]; then
  echo "Mode: range (${START_AT}..${END_AT}) intersected with available batches"
  integer i
  for (( i=START_AT; i<=END_AT; i++ )); do CANDIDATE_IDXS+="$i"; done
else
  echo "Mode: all indices found in ${BATCH_DIR}"
  CANDIDATE_IDXS=("${ALL_IDXS[@]}")
fi

# Validate candidates: numeric and have a matching file
typeset -a CLEAN_INDICES
for idx in "${CANDIDATE_IDXS[@]}"; do
  [[ "$idx" == <-> ]] || { echo "Skip non-numeric index: '$idx'" >&2; continue; }
  if [[ -z "${FILE_BY_IDX[$idx]-}" ]]; then
    echo "Skip index $idx: no batch file found in $BATCH_DIR" >&2
    continue
  fi
  CLEAN_INDICES+="$idx"
done

(( ${#CLEAN_INDICES[@]} > 0 )) || { echo "Nothing to process after validation."; exit 0; }
echo "Will process species (1-based): ${CLEAN_INDICES[*]}"

# --- Manifest seed ---
mkdir -p "$OUT_ROOT"
echo '[]' > "$MANIFEST.tmp"

# Helper: test if a batch has any non-zero across samples
is_batch_nonempty() {
  local file="$1"
  # dataset is samples×1; pull the single column and check any != 0
  jq -e '
    ( [ .dataset[] | .[0] ] | any(. != 0) )
  ' "$file" >/dev/null
}

# --- Main loop (reusing existing batches) ---
for col1 in "${CLEAN_INDICES[@]}"; do
  in_json="${FILE_BY_IDX[$col1]}"
  out_json="${RUN_DIR}/pred_out_species_${col1}.json"

  # Optional emptiness filter
  NONEMPTY=1
  if (( ONLY_NONEMPTY == 1 )); then
    if ! is_batch_nonempty "$in_json"; then
      NONEMPTY=0
    fi
  fi

  if (( NONEMPTY == 0 )); then
    echo "→ Skip species ${col1}: all-zero across ${TOTAL_SAMPLES} samples."
    jq --argjson idx "$col1" '. + [{species: $idx, status:"skipped_empty"}]' "$MANIFEST.tmp" > "$MANIFEST.tmp2" && mv "$MANIFEST.tmp2" "$MANIFEST.tmp"
    continue
  fi

  # Validate: samples × width=1
  jq -e --argjson ns "$TOTAL_SAMPLES" '
    (.dataset | length) == $ns and
    ([.dataset[] | length] | all(. == 1))
  ' "$in_json" >/dev/null || { echo "ERROR: validation failed for species ${col1} (file: $in_json)"; exit 1; }

  # Skip or run
  if (( OVERWRITE == 0 )) && [[ -f "$out_json" ]]; then
    echo "→ Skip species ${col1}: output exists (use OVERWRITE=1 to rerun)."
    jq --argjson idx "$col1" '. + [{species: $idx, status:"skipped_exists"}]' "$MANIFEST.tmp" > "$MANIFEST.tmp2" && mv "$MANIFEST.tmp2" "$MANIFEST.tmp"
    continue
  fi

  echo "→ Running species ${col1} …"
  "$MODEL_BIN" "$in_json" --sweeps "$SWEEPS" > "$out_json"

  jq --argjson idx "$col1" '. + [{species: $idx, status:"done"}]' "$MANIFEST.tmp" > "$MANIFEST.tmp2" && mv "$MANIFEST.tmp2" "$MANIFEST.tmp"
done

# --- Finalize manifest ---
jq --argjson total_samples "$TOTAL_SAMPLES" '
  { meta: { total_samples: $total_samples },
    runs: .
  }
  | .kept    = [ .runs[] | select(.status=="done") | .species ]
  | .skipped = [ .runs[] | select(.status!="done") | {species, status} ]
' "$MANIFEST.tmp" > "$MANIFEST"
rm -f "$MANIFEST.tmp"

echo "Done. Manifest: $MANIFEST"
echo "Inputs (reused): $BATCH_DIR"
echo "Outputs:         $RUN_DIR"
echo "Tip: set INDEX_LIST=\"...\" or START_AT/END_AT, OVERWRITE=1 to rerun."
