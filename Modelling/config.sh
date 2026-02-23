#!/usr/bin/env zsh
# =============================================================================
# config.sh — central path configuration for all run scripts
#
# Usage: source this file from any run script AFTER setting MODELLING_ROOT.
#
# Each script computes MODELLING_ROOT relative to its own location, then:
#   source "$MODELLING_ROOT/config.sh"
#
# Prerequisites:
#   - tpplc (TreePPL compiler) must be in PATH.
#     Build from: https://github.com/miking-lang/miking-treeppl
#     Follow the build instructions there; no other external repo is needed.
# =============================================================================

# --- Core directories ---
MODELS_DIR="$MODELLING_ROOT/Models"
DATA_DIR="$MODELLING_ROOT/Data"
RESULTS_DIR="$MODELLING_ROOT/Results"
BIN_DIR="$MODELLING_ROOT/bin"          # compiled TreePPL binaries land here

# --- Step 1: model selection ---
SMC_LOGZ_DIR="$RESULTS_DIR/step1_model_selection/smc_logZ"
MCMC_CHAINS_DIR="$RESULTS_DIR/step1_model_selection/mcmc/chains"

# --- Step 2: predictions ---
PRED_OUT_DIR="$RESULTS_DIR/step2_predictions/prediction_outputs"
STITCHED_DIR="$RESULTS_DIR/step2_predictions/stitched"
SPIKEIN_INPUT_DIR="$RESULTS_DIR/step2_predictions/prediction_inputs/spikeins"
FULLSET_BATCH_DIR="$MODELLING_ROOT/Run step 2/pred_fullset_runs/batches"

# --- Ensure bin/ exists ---
mkdir -p "$BIN_DIR"
