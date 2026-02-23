#!/bin/zsh

#chmod +x pred_spikeins_runs_batches_6k.sh
#./pred_spikeins_runs_batches_6k.sh

set -euo pipefail

echo "Starting script...pred spikein set (batched by species)"

SCRIPT_DIR="$(cd "$(dirname "${(%):-%x}")" && pwd)"
MODELLING_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
source "$MODELLING_ROOT/config.sh"

tpplc "$MODELS_DIR/predict/optmodel15_Pred_match_inform_pfit.tppl" \
  -m smc-apf \
  --particles 2000 \
  --subsample \
  --subsample-size 1 \
  --resample manual \
  --output "$BIN_DIR/optmodel15_Pred_match_inform_pfit_2000"

MODEL_BIN="$BIN_DIR/optmodel15_Pred_match_inform_pfit_2000"
ORIG_JSON="$DATA_DIR/modified_pred_in_8k_mcmc.json" #spikein data replaces full dataset, 8k model

SWEEPS=${SWEEPS:-100}
OUT_ROOT="${OUT_ROOT:-$PRED_OUT_DIR/spikeins/inform_pfit/8k}"
BATCH_DIR="${BATCH_DIR:-${OUT_ROOT}/batches}"
RUN_DIR="${OUT_ROOT}/runs"
MANIFEST="${OUT_ROOT}/manifest.json"

# --- Resume/limits/flags (env override) ---
# 1-based, inclusive range to run; default = all species
START_AT=${START_AT:-1}
END_AT=${END_AT:-""}       # empty means auto-detect total species
OVERWRITE=${OVERWRITE:-1}  # 0=skip species with existing output; 1=rerun and overwrite
ONLY_NONEMPTY=${ONLY_NONEMPTY:-1} # 1=skip all-zero species; 0=run all

mkdir -p "$RUN_DIR"

command -v jq >/dev/null 2>&1 || { echo "ERROR: jq is required"; exit 1; }
[[ -f "$ORIG_JSON" ]] || { echo "ERROR: input JSON not found: $ORIG_JSON"; exit 1; }
[[ -x "$MODEL_BIN" ]] || { echo "ERROR: model binary not found/executable: $MODEL_BIN"; exit 1; }

# --- Detect shape ---
TOTAL_SPECIES=$(jq '(.dataset[0] | length)' "$ORIG_JSON")
TOTAL_SAMPLES=$(jq '(.dataset | length)' "$ORIG_JSON")
[[ "$TOTAL_SPECIES" =~ ^[0-9]+$ ]] || { echo "ERROR: cannot read species count"; exit 1; }
[[ "$TOTAL_SAMPLES" =~ ^[0-9]+$ ]] || { echo "ERROR: cannot read sample count"; exit 1; }
echo "Detected ${TOTAL_SPECIES} species × ${TOTAL_SAMPLES} samples."

# Compute END_AT if empty
if [[ -z "$END_AT" ]]; then END_AT=$TOTAL_SPECIES; fi
if (( START_AT < 1 || END_AT < START_AT || END_AT > TOTAL_SPECIES )); then
  echo "ERROR: invalid START_AT/END_AT ($START_AT..$END_AT)"; exit 1
fi
echo "Will process species (1-based): ${START_AT}..${END_AT}"

# --- Precompute which species are non-empty (any nonzero across samples) ---
# Produces arrays of 1-based indices for kept/skipped species; also per-species sums.
NONEMPTY_JSON=$(jq '
  (.dataset | transpose) as $cols
  | [ range(0; ($cols|length)) | {
      idx1: (.+1),
      sum: ($cols[.] | add),
      any_nonzero: ($cols[.] | any(. != 0))
    } ]
' "$ORIG_JSON")

# Helper to query nonempty flag fast using jq
is_nonempty() {
  local idx1="$1"
  echo "$NONEMPTY_JSON" | jq -e --argjson i "$idx1" '
    .[] | select(.idx1 == $i) | .any_nonzero
  ' >/dev/null
}

# --- Slice one species (keeps model input shape: samples × 1 column) ---
slice_one_species() {
  local col0="$1" in="$2" out="$3"
  local e=$(( col0 + 1 ))
  jq --argjson s "$col0" --argjson e "$e" '
    def slice_cols(a): (a | transpose | .[$s:$e] | transpose);
    .dataset = slice_cols(.dataset)
    | .species_start = ($s + 1)
    | .species_end   = ($s + 1)
  ' "$in" > "$out"
}

# --- Manifest seed ---
# We’ll append records as we go, then compact/format at end.
echo '[]' > "$MANIFEST.tmp"

# --- Main loop ---
for (( col1=START_AT; col1<=END_AT; col1++ )); do
  col0=$(( col1 - 1 ))
  in_json="${BATCH_DIR}/pred_in_species_${col1}.json"
  out_json="${RUN_DIR}/pred_out_species_${col1}.json"

  # Check emptiness
  NONEMPTY=1
  if (( ONLY_NONEMPTY == 1 )); then
    if ! is_nonempty "$col1"; then
      NONEMPTY=0
    fi
  fi

  # Skip if empty and ONLY_NONEMPTY
  if (( NONEMPTY == 0 )); then
    echo "→ Skip species ${col1}: all-zero across ${TOTAL_SAMPLES} samples."
    # Append manifest record
    jq --argjson idx "$col1" '. + [{species: $idx, status:"skipped_empty"}]' "$MANIFEST.tmp" > "$MANIFEST.tmp2" && mv "$MANIFEST.tmp2" "$MANIFEST.tmp"
    continue
  fi

  # Skip if output exists and not overwriting
  if (( OVERWRITE == 0 )) && [[ -f "$out_json" ]]; then
    echo "→ Skip species ${col1}: output exists (use OVERWRITE=1 to rerun)."
    jq --argjson idx "$col1" '. + [{species: $idx, status:"skipped_exists"}]' "$MANIFEST.tmp" > "$MANIFEST.tmp2" && mv "$MANIFEST.tmp2" "$MANIFEST.tmp"
    continue
  fi

  echo "→ Running species ${col1}/${TOTAL_SPECIES} …"
  slice_one_species "$col0" "$ORIG_JSON" "$in_json"

  # Validate: samples × width=1
  jq -e --argjson ns "$TOTAL_SAMPLES" '
    (.dataset | length) == $ns and
    ([.dataset[] | length] | all(. == 1))
  ' "$in_json" >/dev/null || { echo "ERROR: validation failed for species ${col1}"; exit 1; }

  # Run model
  "$MODEL_BIN" "$in_json" --sweeps "$SWEEPS" > "$out_json"

  # Append manifest record
  jq --argjson idx "$col1" '. + [{species: $idx, status:"done"}]' "$MANIFEST.tmp" > "$MANIFEST.tmp2" && mv "$MANIFEST.tmp2" "$MANIFEST.tmp"
done

# --- Finalize manifest: add metadata and split lists ---
jq --argjson total_species "$TOTAL_SPECIES" --argjson total_samples "$TOTAL_SAMPLES" '
  { meta: { total_species: $total_species, total_samples: $total_samples },
    runs: .
  }
  | .kept    = [ .runs[] | select(.status=="done") | .species ]
  | .skipped = [ .runs[] | select(.status!="done") | {species, status} ]
' "$MANIFEST.tmp" > "$MANIFEST"
rm -f "$MANIFEST.tmp"

echo "Done. Manifest: $MANIFEST"
echo "Inputs:  $BATCH_DIR"
echo "Outputs: $RUN_DIR"
echo "Tip: set START_AT=… END_AT=… to resume a subset; set OVERWRITE=1 to rerun existing."
