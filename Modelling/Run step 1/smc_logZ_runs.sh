#!/bin/zsh
# Run SMC log normalising constant for all models
#chmod +x smc_logZ_runs.sh
#./smc_logZ_runs.sh

SCRIPT_DIR="$(cd "$(dirname "${(%):-%x}")" && pwd)"
MODELLING_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$MODELLING_ROOT/config.sh"

# =============================================================================
# Art models
# =============================================================================
echo "Starting script... art models - log Z"
tpplc "$MODELS_DIR/art/optmodel_art.tppl"        -m 'smc-apf' --resample manual --particles 500000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel_art"
tpplc "$MODELS_DIR/art/optmodel15_1c.tppl"       -m 'smc-apf' --resample manual --particles 500000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel15_1c"
tpplc "$MODELS_DIR/art/optmodel15_2k.tppl"       -m 'smc-apf' --resample manual --particles 500000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel15_2k"
tpplc "$MODELS_DIR/art/optmodel15_2theta.tppl"   -m 'smc-apf' --resample manual --particles 500000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel15_2theta"
tpplc "$MODELS_DIR/art/optmodel15_2k2theta.tppl" -m 'smc-apf' --resample manual --particles 500000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel15_2k2theta"
OUTPUT_DIR="$SMC_LOGZ_DIR/art_logZ"
mkdir -p "$OUTPUT_DIR"
"$BIN_DIR/optmodel_art"        "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel_art.json"
"$BIN_DIR/optmodel15_1c"       "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel15_1c.json"
"$BIN_DIR/optmodel15_2k"       "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel15_2k.json"
"$BIN_DIR/optmodel15_2theta"   "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel15_2theta.json"
"$BIN_DIR/optmodel15_2k2theta" "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel15_2k2theta.json"
echo "Script done! art models - log Z"

# =============================================================================
# Bio models
# =============================================================================
echo "Starting script... bio models - log Z"
tpplc "$MODELS_DIR/bio/optmodel15_1c_biospikeins.tppl" -m 'smc-apf' --resample manual --particles 1000000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel15_1c_biospikeins"
tpplc "$MODELS_DIR/bio/optmodel15_6k.tppl"             -m 'smc-apf' --resample manual --particles 1000000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel15_6k"
tpplc "$MODELS_DIR/bio/optmodel15_6theta.tppl"         -m 'smc-apf' --resample manual --particles 1000000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel15_6theta"
tpplc "$MODELS_DIR/bio/optmodel15_6k6theta.tppl"       -m 'smc-apf' --resample manual --particles 1000000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel15_6k6theta"
tpplc "$MODELS_DIR/bio/optmodel15_biospikeins.tppl"    -m 'smc-apf' --resample manual --particles 1000000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel15_biospikeins"
OUTPUT_DIR="$SMC_LOGZ_DIR/bio_logZ"
mkdir -p "$OUTPUT_DIR"
"$BIN_DIR/optmodel15_1c_biospikeins" "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel15_1c_biospikeins.json"
"$BIN_DIR/optmodel15_6k"             "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel15_6k.json"
"$BIN_DIR/optmodel15_6theta"         "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel15_6theta.json"
"$BIN_DIR/optmodel15_6k6theta"       "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel15_6k6theta.json"
"$BIN_DIR/optmodel15_biospikeins"    "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel15_biospikeins.json"
echo "Script done! bio models - log Z"

# =============================================================================
# Combined models
# =============================================================================
echo "Starting script... combined models - log Z"
tpplc "$MODELS_DIR/combined/optmodel15_1c_combined.tppl"    -m 'smc-apf' --resample manual --particles 1000000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel15_1c_combined"
tpplc "$MODELS_DIR/combined/optmodel15_combined.tppl"       -m 'smc-apf' --resample manual --particles 1000000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel15_combined"
tpplc "$MODELS_DIR/combined/optmodel15_8k.tppl"             -m 'smc-apf' --resample manual --particles 1000000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel15_8k"
tpplc "$MODELS_DIR/combined/optmodel15_8k7theta.tppl"       -m 'smc-apf' --resample manual --particles 1000000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel15_8k7theta"
tpplc "$MODELS_DIR/combined/optmodel15_8k8theta.tppl"       -m 'smc-apf' --resample manual --particles 1000000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel15_8k8theta"
tpplc "$MODELS_DIR/combined/optmodel15_combined_2theta.tppl" -m 'smc-apf' --resample manual --particles 1000000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel15_combined_2theta"
OUTPUT_DIR="$SMC_LOGZ_DIR/combined_logZ"
mkdir -p "$OUTPUT_DIR"
"$BIN_DIR/optmodel15_1c_combined"     "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel15_1c_combined.json"
"$BIN_DIR/optmodel15_combined"        "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel15_combined.json"
"$BIN_DIR/optmodel15_8k"              "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel15_8k.json"
"$BIN_DIR/optmodel15_8k7theta"        "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel15_8k7theta.json"
"$BIN_DIR/optmodel15_8k8theta"        "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel15_8k8theta.json"
"$BIN_DIR/optmodel15_combined_2theta" "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel15_combined_2theta.json"
echo "Script done! combined models - log Z"

# =============================================================================
# Fixed theta models
# =============================================================================
echo "Starting script... fixed theta models - log Z"
tpplc "$MODELS_DIR/fixTheta/optmodel15_combined_fixTheta.tppl"    -m 'smc-apf' --resample manual --particles 1000000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel_smc_comb_fixTheta"
tpplc "$MODELS_DIR/fixTheta/optmodel15_biospikeins_fixTheta.tppl" -m 'smc-apf' --resample manual --particles 1000000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel_smc_biospikeins_fixTheta"
tpplc "$MODELS_DIR/fixTheta/optmodel15_art_fixTheta.tppl"         -m 'smc-apf' --resample manual --particles 1000000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel_smc_art_fixTheta"
tpplc "$MODELS_DIR/fixTheta/optmodel15_6k_fixTheta.tppl"          -m 'smc-apf' --resample manual --particles 1000000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel_smc_6k_fixTheta"
tpplc "$MODELS_DIR/fixTheta/optmodel15_8k_fixTheta.tppl"          -m 'smc-apf' --resample manual --particles 1000000 --subsample --subsample-size 1 --output "$BIN_DIR/optmodel_smc_8k_fixTheta"
OUTPUT_DIR="$SMC_LOGZ_DIR/fix_theta_logZ"
mkdir -p "$OUTPUT_DIR"
"$BIN_DIR/optmodel_smc_comb_fixTheta"        "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel_smc_comb_fixTheta.json"
"$BIN_DIR/optmodel_smc_biospikeins_fixTheta" "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel_smc_biospikeins_fixTheta.json"
"$BIN_DIR/optmodel_smc_art_fixTheta"         "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel_smc_art_fixTheta.json"
"$BIN_DIR/optmodel_smc_6k_fixTheta"          "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel_smc_6k_fixTheta.json"
"$BIN_DIR/optmodel_smc_8k_fixTheta"          "$DATA_DIR/15_hom_priors.json" --sweeps 100 > "${OUTPUT_DIR}/optmodel_smc_8k_fixTheta.json"
echo "Script done! fixed theta models - log Z"

# =============================================================================
# Amend JSON output (wrap NDJSON into JSON array)
# =============================================================================
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/art_logZ/optmodel_art.json"
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/art_logZ/optmodel15_1c.json"
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/art_logZ/optmodel15_2k.json"
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/art_logZ/optmodel15_2theta.json"
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/art_logZ/optmodel15_2k2theta.json"
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/bio_logZ/optmodel15_1c_biospikeins.json"
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/bio_logZ/optmodel15_6k.json"
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/bio_logZ/optmodel15_6theta.json"
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/bio_logZ/optmodel15_6k6theta.json"
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/bio_logZ/optmodel15_biospikeins.json"
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/combined_logZ/optmodel15_1c_combined.json"
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/combined_logZ/optmodel15_combined.json"
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/combined_logZ/optmodel15_8k.json"
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/combined_logZ/optmodel15_8k7theta.json"
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/combined_logZ/optmodel15_8k8theta.json"
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/combined_logZ/optmodel15_combined_2theta.json"
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/fix_theta_logZ/optmodel_smc_comb_fixTheta.json"
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/fix_theta_logZ/optmodel_smc_biospikeins_fixTheta.json"
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/fix_theta_logZ/optmodel_smc_art_fixTheta.json"
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/fix_theta_logZ/optmodel_smc_6k_fixTheta.json"
"$MODELLING_ROOT/json_amend.zsh" "$SMC_LOGZ_DIR/fix_theta_logZ/optmodel_smc_8k_fixTheta.json"

echo "All done!"
