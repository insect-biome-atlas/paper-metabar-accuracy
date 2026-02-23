#!/bin/zsh
# Run MCMC chains for all fixed-theta models (3 chains each)
#chmod +x mcmc_runs.sh
#./mcmc_runs.sh

SCRIPT_DIR="$(cd "$(dirname "${(%):-%x}")" && pwd)"
MODELLING_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$MODELLING_ROOT/config.sh"

mkdir -p "$MCMC_CHAINS_DIR"

# =============================================================================
# Compile models
# =============================================================================
echo "Starting MCMC script... compiling"
tpplc "$MODELS_DIR/mcmc/optmodel15_combined_fixTheta.tppl"    -m 'mcmc-lightweight' --kernel --mcmc-lw-gprob 0 --sampling-period 100 --particles 2000000 --output "$BIN_DIR/optmodel_mcmc_comb_fixTheta"
tpplc "$MODELS_DIR/mcmc/optmodel15_biospikeins_fixTheta.tppl" -m 'mcmc-lightweight' --kernel --mcmc-lw-gprob 0 --sampling-period 100 --particles 2000000 --output "$BIN_DIR/optmodel15_mcmc_biospikeins_fixTheta"
tpplc "$MODELS_DIR/mcmc/optmodel15_6k_fixTheta.tppl"          -m 'mcmc-lightweight' --kernel --mcmc-lw-gprob 0 --sampling-period 100 --particles 5000000 --output "$BIN_DIR/optmodel_mcmc_6k_fixTheta"
tpplc "$MODELS_DIR/mcmc/optmodel15_8k_fixTheta.tppl"          -m 'mcmc-lightweight' --kernel --mcmc-lw-gprob 0 --sampling-period 100 --particles 5000000 --output "$BIN_DIR/optmodel_mcmc_8k_fixTheta"
tpplc "$MODELS_DIR/mcmc/optmodel15_art_fixTheta.tppl"         -m 'mcmc-lightweight' --kernel --mcmc-lw-gprob 0 --sampling-period 100 --particles 2000000 --output "$BIN_DIR/optmodel_mcmc_art_fixTheta"

# =============================================================================
# Run chains (3 per model)
# =============================================================================
echo "Starting MCMC script... running chains"
"$BIN_DIR/optmodel_mcmc_comb_fixTheta"          "$DATA_DIR/15_hom_priors.json" > "$MCMC_CHAINS_DIR/optmodel_mcmc_comb_fixTheta_1.json"
"$BIN_DIR/optmodel_mcmc_comb_fixTheta"          "$DATA_DIR/15_hom_priors.json" > "$MCMC_CHAINS_DIR/optmodel_mcmc_comb_fixTheta_2.json"
"$BIN_DIR/optmodel_mcmc_comb_fixTheta"          "$DATA_DIR/15_hom_priors.json" > "$MCMC_CHAINS_DIR/optmodel_mcmc_comb_fixTheta_3.json"
"$BIN_DIR/optmodel15_mcmc_biospikeins_fixTheta" "$DATA_DIR/15_hom_priors.json" > "$MCMC_CHAINS_DIR/optmodel_mcmc_biospikeins_fixTheta_1.json"
"$BIN_DIR/optmodel15_mcmc_biospikeins_fixTheta" "$DATA_DIR/15_hom_priors.json" > "$MCMC_CHAINS_DIR/optmodel_mcmc_biospikeins_fixTheta_2.json"
"$BIN_DIR/optmodel15_mcmc_biospikeins_fixTheta" "$DATA_DIR/15_hom_priors.json" > "$MCMC_CHAINS_DIR/optmodel_mcmc_biospikeins_fixTheta_3.json"
"$BIN_DIR/optmodel_mcmc_6k_fixTheta"            "$DATA_DIR/15_hom_priors.json" > "$MCMC_CHAINS_DIR/optmodel_mcmc_6k_fixTheta_1.json"
"$BIN_DIR/optmodel_mcmc_6k_fixTheta"            "$DATA_DIR/15_hom_priors.json" > "$MCMC_CHAINS_DIR/optmodel_mcmc_6k_fixTheta_2.json"
"$BIN_DIR/optmodel_mcmc_6k_fixTheta"            "$DATA_DIR/15_hom_priors.json" > "$MCMC_CHAINS_DIR/optmodel_mcmc_6k_fixTheta_3.json"
"$BIN_DIR/optmodel_mcmc_8k_fixTheta"            "$DATA_DIR/15_hom_priors.json" > "$MCMC_CHAINS_DIR/optmodel_mcmc_8k_fixTheta_1.json"
"$BIN_DIR/optmodel_mcmc_8k_fixTheta"            "$DATA_DIR/15_hom_priors.json" > "$MCMC_CHAINS_DIR/optmodel_mcmc_8k_fixTheta_2.json"
"$BIN_DIR/optmodel_mcmc_8k_fixTheta"            "$DATA_DIR/15_hom_priors.json" > "$MCMC_CHAINS_DIR/optmodel_mcmc_8k_fixTheta_3.json"
"$BIN_DIR/optmodel_mcmc_art_fixTheta"           "$DATA_DIR/15_hom_priors.json" > "$MCMC_CHAINS_DIR/optmodel_mcmc_art_fixTheta_1.json"
"$BIN_DIR/optmodel_mcmc_art_fixTheta"           "$DATA_DIR/15_hom_priors.json" > "$MCMC_CHAINS_DIR/optmodel_mcmc_art_fixTheta_2.json"
"$BIN_DIR/optmodel_mcmc_art_fixTheta"           "$DATA_DIR/15_hom_priors.json" > "$MCMC_CHAINS_DIR/optmodel_mcmc_art_fixTheta_3.json"

echo "Done!"
