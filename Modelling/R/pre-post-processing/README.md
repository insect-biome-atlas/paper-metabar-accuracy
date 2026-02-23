# TreePPL Abundance Estimates - supporting functions in R

**Project**: Metabarcoding abundance estimation using TreePPL probabilistic programming  
**Author**: Emma  
**Last Updated**: February 2026

---

## 📋 OVERVIEW

This repository contains a complete R pipeline for:
1. **Training** Bayesian models via MCMC in TreePPL
2. **Validating** convergence and model quality
3. **Predicting** species abundance from metabarcoding data
4. **Processing** prediction outputs for downstream analysis

**Total Files**: 11 production-ready R scripts organized in 2 main pipelines

---

## 🗂️ REPOSITORY STRUCTURE

```
.
├── README.md                                    # This file - Master overview
├── README_MCMC_TRAINING.md                      # Pipeline 1: Training & postprocessing
├── README_STITCHING.md                          # Pipeline 2: Prediction processing
│
├── MCMC_TRAINING_PIPELINE/
│   ├── mcmc_convergence.R                       # Core: Convergence diagnostics
│   ├── run_mcmc_convergence.R                   # Runner: Batch convergence checks
│   ├── postprocessing_after_mcmc.R              # Core: MCMC → prediction transform
│   ├── run_postprocessing_after_mcmc.R          # Runner: Batch postprocessing
│   └── create_species_mapping.R                 # Utility: Species ID mapping
│
├── PREDICTION_STITCHING_PIPELINE/
│   ├── stitch_simple_n_models_spikeins.R        # Core: Stitch simple models
│   ├── stitch_match_models_spikeins.R           # Core: Stitch match models
│   ├── stitch_predictions_fullset.R             # Core: Stitch with QC (advanced)
│   ├── run_stitcher_spikeins.r                  # Runner: Batch spike-in stitching
│   └── run_stitcher_fullset.r                   # Runner: Full dataset stitching
│
└── data/
    ├── species_map.csv                          # Species ID → name mapping
    ├── cleaned_nochimera_MATCHED_cluster_counts_*.csv
    └── species_map_diagnostics/                 # Validation outputs
```

---

## 🔄 COMPLETE WORKFLOW

### The Big Picture

```
┌─────────────────┐
│   TreePPL MCMC  │  Run Bayesian MCMC training
│   (3+ chains)   │
└────────┬────────┘
         │
         ↓
┌─────────────────────────────────────────────────┐
│  PIPELINE 1: MCMC TRAINING & POSTPROCESSING     │
├─────────────────────────────────────────────────┤
│  1. Check convergence (R-hat, ESS, Geweke)      │
│  2. Validate all chains converged               │
│  3. Transform samples for prediction            │
│  4. Generate prediction input JSONs             │
└────────┬────────────────────────────────────────┘
         │
         ↓
┌─────────────────┐
│ TreePPL Predict │  Run predictions on test data
│ (per species)   │
└────────┬────────┘
         │
         ↓
┌─────────────────────────────────────────────────┐
│  PIPELINE 2: PREDICTION STITCHING              │
├─────────────────────────────────────────────────┤
│  1. Create species mapping                      │
│  2. Stitch NDJSON outputs                       │
│  3. Validate data quality                       │
│  4. Generate RDS for analysis                   │
└────────┬────────────────────────────────────────┘
         │
         ↓
┌─────────────────┐
│ Analysis/Plots  │  Downstream analysis & visualization
│  (Separate)     │  (Plotting functions in separate folder)
└─────────────────┘
```

---

## 🚀 QUICK START GUIDE

### Prerequisites

```r
# Install required packages
install.packages(c(
  "jsonlite", "dplyr", "tidyr", "stringr", "purrr", 
  "tibble", "coda", "ggplot2", "readr"
))
```

### Step-by-Step Workflow

#### **Stage 1: MCMC Training** (After TreePPL MCMC runs)

```r
# 1. Check convergence for all models
source("MCMC_TRAINING_PIPELINE/run_mcmc_convergence.R")
# → Creates mcmc_check_*/ directories with diagnostics

# 2. Review convergence
# Check mcmc_check_biospikeins/convergence_summary.csv (and others)
# Ensure all parameters show "converged" verdict

# 3. Postprocess MCMC outputs
source("MCMC_TRAINING_PIPELINE/run_postprocessing_after_mcmc.R")
# → Creates pred_in_*.json files for TreePPL prediction step
```

**Outputs**: 5 JSON files ready for prediction (art, bio, comb, 6k, 8k)

---

#### **Stage 2: Prediction Processing** (After TreePPL predictions)

**Option A: Spike-in Validation**

```r
# 1. Stitch all spike-in validation models
source("PREDICTION_STITCHING_PIPELINE/run_stitcher_spikeins.r")
# → Creates stitched_pred_spikeins_*.rds files

# 2. Ready for plotting/analysis (use plotting pipeline)
```

**Option B: Full Dataset with QC**

```r
# 1. Create species mapping
source("PREDICTION_STITCHING_PIPELINE/create_species_mapping.R")
# → Creates data/species_map.csv + diagnostics

# 2. Review mapping diagnostics
# Check data/species_map_diagnostics/*.csv for mismatches

# 3. Stitch predictions with quality control
source("PREDICTION_STITCHING_PIPELINE/stitch_predictions_fullset.R")
# → Creates stitched_pred_fullset.rds + QC diagnostics

# 4. Review stitching diagnostics
# Check diag_sweeps.csv, diag_species_lengths.csv, etc.

# 5. Ready for analysis
```

---

## 📚 DETAILED DOCUMENTATION

### Pipeline 1: MCMC Training & Postprocessing

**See**: `README_MCMC_TRAINING.md`

**Key Files**:
- `mcmc_convergence.R` - Computes R-hat, ESS, Geweke diagnostics
- `postprocessing_after_mcmc.R` - Transforms MCMC samples for prediction
- `create_species_mapping.R` - Builds TreePPL-compatible species mapping

**What it does**:
- Validates MCMC convergence across multiple chains
- Checks R-hat ≤ 1.05, ESS ≥ 400, |Geweke| < 2
- Generates trace plots, density plots, diagnostic CSVs
- Thins samples to 500 per model
- Formats JSON for TreePPL prediction input
- Creates species ID mapping matching TreePPL indexing

**Inputs**: TreePPL MCMC output JSONs (3+ chains per model)  
**Outputs**: Prediction input JSONs + convergence diagnostics + species mapping

---

### Pipeline 2: Prediction Stitching

**See**: `README_STITCHING.md` (to be created based on README_STITCHING_PLOTTING.md)

**Key Files**:
- `stitch_simple_n_models_spikeins.R` - For bio/art/comb models
- `stitch_match_models_spikeins.R` - For 6k/8k models (with m parameter)
- `stitch_predictions_fullset.R` - Advanced stitcher with QC

**What it does**:
- Reads NDJSON prediction outputs (pred_out_species_*.json)
- Aggregates predictions across all species
- Validates against count data (fullset version)
- Enforces expected vector lengths
- Generates extensive diagnostics
- Outputs compact RDS for analysis

**Inputs**: TreePPL prediction outputs (NDJSON per species)  
**Outputs**: Stitched RDS files + diagnostic CSVs

---

## 🎯 COMMON WORKFLOWS

### Workflow 1: Standard 5-Model Training & Validation

**Goal**: Train, validate, and prepare 5 models for prediction

```bash
# After TreePPL MCMC training completes:
```

```r
# Step 1: Check all models converged
source("MCMC_TRAINING_PIPELINE/run_mcmc_convergence.R")

# Step 2: Review convergence (all should show "converged")
read.csv("mcmc_check_biospikeins/convergence_summary.csv")
read.csv("mcmc_check_comb/convergence_summary.csv")
read.csv("mcmc_check_6k/convergence_summary.csv")
read.csv("mcmc_check_8k/convergence_summary.csv")
read.csv("mcmc_check_art/convergence_summary.csv")

# Step 3: If all converged, generate prediction inputs
source("MCMC_TRAINING_PIPELINE/run_postprocessing_after_mcmc.R")

# Ready for TreePPL prediction step!
```

---

### Workflow 2: Spike-in Validation Analysis

**Goal**: Compare model performance on 6 known spike-ins

```bash
# After TreePPL predictions complete:
```

```r
# Stitch all 5 models
source("PREDICTION_STITCHING_PIPELINE/run_stitcher_spikeins.r")

# Files created:
# - stitched_pred_spikeins_bio.rds
# - stitched_pred_spikeins_art.rds
# - stitched_pred_spikeins_comb.rds
# - stitched_pred_spikeins_6k.rds
# - stitched_pred_spikeins_8k.rds

# Use plotting pipeline (separate folder) to visualize results
```

---

### Workflow 3: Full Dataset Production Analysis

**Goal**: Process all species with quality control

```bash
# After TreePPL predictions complete:
```

```r
# Step 1: Create species mapping
source("PREDICTION_STITCHING_PIPELINE/create_species_mapping.R")

# Step 2: Review mapping diagnostics
list.files("data/species_map_diagnostics/")
# Check for:
# - runs_ids_not_in_map.csv (unexpected IDs)
# - map_ids_missing_files.csv (missing outputs)

# Step 3: Stitch with quality control
source("PREDICTION_STITCHING_PIPELINE/stitch_predictions_fullset.R")

# Step 4: Review stitching diagnostics
read.csv("output_dir/diag_sweeps.csv")          # Per-sweep decisions
read.csv("output_dir/diag_species_lengths.csv") # Summary by species
read.csv("output_dir/diag_bad_lengths.csv")     # Dropped sweeps

# Step 5: Use stitched_pred_fullset.rds for analysis
```

---

## 🔧 CONFIGURATION

### Quick Configuration Checklist

**Before running MCMC convergence**:
- [ ] Update `chain_paths` in `run_mcmc_convergence.R`
- [ ] Set appropriate `burn_in` (default: 0.5)
- [ ] Set `rhat_thresh` (default: 1.05)
- [ ] Set `ess_thresh` (default: 300)

**Before postprocessing**:
- [ ] Verify convergence passed for all models
- [ ] Update `dir_in` in `run_postprocessing_after_mcmc.R`
- [ ] Update `dir_out` for prediction inputs
- [ ] Update `prior_data` path

**Before stitching spike-ins**:
- [ ] Update `runs_dir` paths in `run_stitcher_spikeins.r`
- [ ] Update `manifest_path` for each model
- [ ] Update `out_rds` paths
- [ ] Verify `expected_len_per_sweep = 15` is correct

**Before stitching fullset**:
- [ ] Update `runs_dir` in `stitch_predictions_fullset.R`
- [ ] Update `counts_path` to your count matrix
- [ ] Set `drop_first_n_rows` to match TreePPL preprocessing
- [ ] Set `drop_spikeins = TRUE/FALSE` to match TreePPL
- [ ] Update `biospikeins` list if needed
- [ ] Choose `length_mode = "strict"` or `"truncate"`

---

## 📊 KEY OUTPUTS EXPLAINED

### MCMC Convergence Diagnostics

**Location**: `mcmc_check_*/convergence_summary.csv`

**What to look for**:
- ✅ All `verdict` = "converged"
- ✅ `rhat` ≤ 1.05 for all parameters
- ✅ `ess` ≥ 300 (preferably >400)
- ⚠️ If any "needs_attention": increase burn-in or run longer

**Example**:
```csv
parameter,rhat,ess,geweke_z,verdict
theta,1.002,1245,-0.34,converged
k_1,1.089,234,1.87,needs_attention  ← FIX THIS
c_1,1.01,678,0.12,converged
```

---

### Prediction Input JSON

**Location**: Output from postprocessing (e.g., `pred_in_bio_mcmc.json`)

**Structure**:
```json
{
  "k_bio_dist": [2.3, 1.8, 2.1, ...],     // 500 samples
  "c_dist": [[0.1, ...], [0.11, ...], ...], // 500 × 15
  "theta_list": [0.45, 0.48, ...],        // 500 samples
  "weights": [0.002, 0.002, ...],         // All equal
  // ... prior parameters
}
```

**Use**: Input to TreePPL prediction step

---

### Species Mapping

**Location**: `data/species_map.csv`

**Structure**:
```csv
species_id,species,canonical_index
1,Acarina_cluster1,1
2,Aleyrodidae_cluster1,2
5,Apidae_cluster1,5
```

**Note**: Gaps in species_id (e.g., 3, 4 missing) are normal if zero-count species were dropped

**Diagnostics**: Check `data/species_map_diagnostics/` for validation

---

### Stitched Predictions

**Location**: `stitched_pred_spikeins_*.rds` or `stitched_pred_fullset.rds`

**Structure**:
```r
list(
  n               = numeric(),  # All predictions concatenated
  offsets         = integer(),  # Start index per sweep
  lengths         = integer(),  # Length per sweep
  species_indices = integer(),  # Species ID per sweep
  m               = integer(),  # Match parameter (6k/8k only)
  normconst       = numeric()   # Normalization constant
)
```

**Use**: Input for analysis and plotting

---

## 🐛 COMMON ISSUES & SOLUTIONS

### Issue: MCMC hasn't converged

**Symptoms**:
- R-hat > 1.05
- ESS < 300
- Verdict = "needs_attention"

**Solutions**:
1. Increase burn-in (try 0.6 or 0.7)
2. Run MCMC for more iterations
3. Check trace plots for mixing issues
4. Verify initial values are reasonable

---

### Issue: Species mapping mismatches

**Symptoms**:
- Diagnostic CSVs show unexpected IDs
- `runs_ids_not_in_map.csv` has many entries
- `map_ids_missing_files.csv` has many entries

**Solutions**:
1. Verify `drop_first_n_rows` matches TreePPL exactly
2. Verify `drop_spikeins` matches TreePPL exactly
3. Check that `biospikeins` list is complete
4. Ensure counts file is the one TreePPL used

---

### Issue: Stitching drops many sweeps

**Symptoms**:
- `diag_bad_lengths.csv` has many rows
- Many sweeps show "drop" decision

**Solutions**:
1. Check `expected_len_per_sweep` is correct
2. Use `length_mode = "truncate"` instead of `"strict"`
3. Review prediction outputs - are they malformed?
4. Check TreePPL prediction logs for errors

---

### Issue: Wrong abundances in results

**Symptoms**:
- Predictions don't match spike-in known values
- Results seem shifted or misaligned

**Solutions**:
1. Regenerate species mapping
2. Verify TreePPL preprocessing steps
3. Check manifest.json has correct species IDs
4. Ensure prediction outputs are from correct model

---

## 💡 BEST PRACTICES

### MCMC Training
✅ **DO**:
- Always run 3+ chains for R-hat calculation
- Check convergence before proceeding
- Save all diagnostic plots
- Document settings used

❌ **DON'T**:
- Skip convergence checks
- Use only 1-2 chains
- Delete raw MCMC outputs
- Ignore "needs_attention" warnings

### Prediction Processing
✅ **DO**:
- Create species mapping before stitching
- Validate mapping diagnostics thoroughly
- Review stitching diagnostics
- Keep preprocessing consistent with TreePPL

❌ **DON'T**:
- Skip diagnostic review
- Change preprocessing mid-project
- Ignore length mismatches
- Mix data from different runs

### General
✅ **DO**:
- Version control your scripts
- Document directory paths used
- Keep README files updated
- Use descriptive output directory names

❌ **DON'T**:
- Hardcode temporary paths
- Skip documentation
- Mix different TreePPL runs
- Overwrite important outputs

---

## 📦 FILE DEPENDENCIES

### MCMC Training Pipeline

```
run_mcmc_convergence.R
  └── requires: mcmc_convergence.R

run_postprocessing_after_mcmc.R
  └── requires: postprocessing_after_mcmc.R
  └── requires: jsonlite package

create_species_mapping.R
  └── standalone (no dependencies)
```

### Prediction Stitching Pipeline

```
run_stitcher_spikeins.r
  ├── requires: stitch_match_models_spikeins.R (for 6k, 8k)
  └── requires: stitch_simple_n_models_spikeins.R (for bio, art, comb)

stitch_predictions_fullset.R
  ├── requires: data/species_map.csv
  ├── requires: counts CSV file
  └── standalone otherwise
```

---

## 🔢 PIPELINE STATISTICS

### Processing Scale
- **Models**: 5 configurations (bio, art, comb, 6k, 8k)
- **MCMC chains**: 3 per model (15 total)
- **Parameters tracked**: ~20-30 per model
- **Samples per chain**: Varies (typically 5,000-10,000)
- **Thinned samples**: 500 per model (for prediction)
- **Species processed**: Configurable (6 for spike-ins, 100s for fullset)
- **Sweeps per species**: Varies (prediction-dependent)

### Output Sizes (typical)
- Convergence diagnostics: ~50 KB CSV + ~2 MB plots per model
- Prediction input JSON: ~1-5 MB per model
- Stitched RDS: ~5-50 MB (depending on species count)
- Diagnostic CSVs: ~10-100 KB total

---

## 🆘 GETTING HELP

### Documentation Hierarchy

1. **This README** - Overall workflow and quick start
2. **README_MCMC_TRAINING.md** - Detailed MCMC pipeline documentation
3. **README_STITCHING.md** - Detailed stitching pipeline documentation
4. **Script comments** - Function-level documentation in each file

### Debugging Strategy

1. **Check file paths** - Most errors are path-related
2. **Review diagnostics** - CSVs contain detailed error info
3. **Inspect intermediate outputs** - Don't skip validation steps
4. **Check TreePPL logs** - Upstream errors affect downstream
5. **Verify preprocessing** - Consistency is critical

### Support Resources

- TreePPL documentation: For MCMC and prediction setup
- R CODA package docs: For convergence diagnostic interpretation
- Script header comments: For function-specific guidance

---

## 📅 TYPICAL TIMELINE

**For 5 models with 3 chains each**:

| Stage | Time | Notes |
|-------|------|-------|
| TreePPL MCMC training | Hours-days | Depends on model complexity |
| Convergence checking | 5-10 min | Automated |
| Review diagnostics | 10-20 min | Manual inspection |
| Postprocessing | 1-2 min | Automated |
| TreePPL predictions | Minutes-hours | Depends on data size |
| Species mapping | 1-2 min | One-time setup |
| Stitching | 2-5 min | Depends on species count |
| Review stitching diagnostics | 5-10 min | Manual inspection |

**Total hands-on time**: ~30-60 minutes (plus TreePPL compute time)

---

## 🎓 LEARNING PATH

### For New Users

**Week 1**: MCMC Training Pipeline
1. Read `README_MCMC_TRAINING.md`
2. Run convergence checks on example data
3. Practice interpreting R-hat, ESS values
4. Review trace plots

**Week 2**: Prediction Processing
1. Read `README_STITCHING.md`
2. Create species mapping
3. Stitch spike-in validation data
4. Review diagnostics

**Week 3**: Full Workflow
1. Run complete pipeline end-to-end
2. Troubleshoot any issues
3. Customize for your data

**Week 4**: Advanced Topics
1. Tune convergence parameters
2. Use fullset stitcher with QC
3. Interpret complex diagnostic outputs
4. Optimize for your specific use case

---

## 🔄 VERSION HISTORY

**Current Version**: v1.0 (February 2026)
- Complete MCMC training pipeline
- Spike-in validation stitching
- Full dataset stitching with QC
- Species mapping utility
- Comprehensive documentation

**Recent Updates**:
- ✅ Fixed plotting script (merged complete functions)
- ✅ Added clear file naming (`_spikeins` suffix)
- ✅ Created species mapping validation
- ✅ Added extensive diagnostics
- ✅ Separated plotting to its own folder

---

## 📝 FILE CHECKLIST

Before running analysis, ensure you have:

**MCMC Training Pipeline** (5 files):
- [ ] `mcmc_convergence.R`
- [ ] `run_mcmc_convergence.R`
- [ ] `postprocessing_after_mcmc.R`
- [ ] `run_postprocessing_after_mcmc.R`
- [ ] `create_species_mapping.R`

**Prediction Stitching Pipeline** (3-4 files):
- [ ] `stitch_simple_n_models_spikeins.R`
- [ ] `stitch_match_models_spikeins.R`
- [ ] `run_stitcher_spikeins.r`
- [ ] `stitch_predictions_fullset.R` (optional, for fullset)

**Data Files**:
- [ ] Count matrix CSV
- [ ] Prior parameters JSON
- [ ] Manifest JSON files (per model)

**TreePPL Outputs**:
- [ ] MCMC chain JSONs (3+ per model)
- [ ] Prediction output JSONs (per species)

---

## 🎯 SUCCESS CRITERIA

Your pipeline is working correctly if:

✅ **After convergence checking**:
- All parameters show "converged" verdict
- R-hat values ≤ 1.05
- ESS values ≥ 400
- Trace plots show good mixing

✅ **After postprocessing**:
- 5 prediction input JSON files created
- Files are valid JSON (test with `fromJSON()`)
- Each contains 500 samples

✅ **After species mapping**:
- `species_map.csv` created
- Diagnostic CSVs show minimal mismatches
- IDs align with TreePPL outputs

✅ **After stitching**:
- RDS files contain expected species
- Diagnostic CSVs show acceptable drop rates
- Data structure is correct (test with `str()`)

---

## 🚦 FINAL CHECKLIST

**Before running your analysis**:
- [ ] Read this README completely
- [ ] Read pipeline-specific READMEs
- [ ] Update all file paths in scripts
- [ ] Install all required R packages
- [ ] Verify TreePPL outputs are complete
- [ ] Create backup of original data
- [ ] Set up output directories
- [ ] Document your configuration

