# README

## Versions used

| Component | Branch | Commit |
|---|---|---|
| **miking** | `develop` | `658b91a15ce620b1659546084a967ea48eb13ac0` |
| **miking-dppl** | `master` | `8fa394182bd41b94d53e3c7bc2b696ea4fa00ceb` |
| **treeppl** | `main` | `e480014432ee946ae1b126f0af72e03a50242af5` |

---

## Install TreePPL

### Option 1: Install locally

Please follow the instructions on **treeppl.org** (for developers) to install TreePPL and all dependencies: 
https://treeppl.org/docs/for-developers 

When you reach the part for installing **miking**, **miking-dppl** and **treeppl**, check out the specific branches/commits listed above (fresh clones):

```bash
# miking
git clone https://github.com/miking-lang/miking.git
cd miking
git checkout develop
git checkout 658b91a15ce620b1659546084a967ea48eb13ac0
cd ..

# miking-dppl
git clone https://github.com/miking-lang/miking-dppl.git
cd miking-dppl
git checkout master
git checkout 8fa394182bd41b94d53e3c7bc2b696ea4fa00ceb
cd ..

# treeppl
git clone https://github.com/treeppl/treeppl.git
cd treeppl
git checkout main
git checkout e480014432ee946ae1b126f0af72e03a50242af5
cd ..
```

---

### Option 2: Use Docker containers [TO BE COMPLETED]

Use the available build Docker containers to run the models.

**The images are:**
- `[image-1]`
- `[image-2]`
- `[image-3]`

**They can be run by adjusting the compile and run command from:**
- `XXX`

**to:**
- `XXX`

---

## Results

The `Results/` directory is empty after cloning — it is populated by running the scripts.

A snapshot of pre-computed results is included as **`earlier_results.zip`** (36 MB) in the repo root. To use these instead of rerunning:

```bash
unzip earlier_results.zip
```

This extracts into `Results/` and allows you to skip directly to post-processing and plotting steps.

---

## Structure of directories

| Directory | Contents |
|---|---|
| `Models/` | TreePPL models |
| `Data/` | Input data for abundance modelling |
| `Run step 1/` | Scripts for model selection (SMC logZ + MCMC) |
| `Run step 2/` | Scripts for predictions (fullset + spikeins) |
| `Results/` | Generated output (gitignored; use scripts or unzip `earlier_results.zip`) |
| `R/` | Pre/post-processing and plotting scripts |

---

## Running the analysis

All scripts are self-contained — no path adjustments needed after cloning. The only prerequisite is `tpplc` in your PATH (see Install TreePPL above).

**Step 1 — Model selection:**
1. `Run step 1/smc_logZ_runs.sh`
2. `Run step 1/mcmc_runs.sh`

**Between steps (R):**

3. `R/pre-post-processing/run_mcmc_convergence.R`
4. `R/pre-post-processing/run_postprocessing_after_mcmc.R`

**Step 2 — Predictions:**

5. `Run step 2/pred_spikein_runs/prior_*/pred_spikeins_runs_batches_*.sh`
6. `Run step 2/pred_fullset_runs/pred_fullset_runs_batches_*.sh`

**Post-processing and plotting (R):**

7. `R/pre-post-processing/run_stitcher_spikeins.r`
8. `R/pre-post-processing/stitch_predictions_fullset.R`
9. `R/plotting/`

