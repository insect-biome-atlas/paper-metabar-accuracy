library(here)
repo_root   <- here::here()
pred_out    <- file.path(repo_root, "Results/step2_predictions/prediction_outputs/spikeins")
stitched    <- file.path(repo_root, "Results/step2_predictions/stitched/spikeins")

source(file.path(repo_root, "R/pre-post-processing/stitch_one_species_ndjson_all_sweeps.R"))

# 6k
runs_dir      = file.path(pred_out, "uninform/6k/runs")
manifest_path = file.path(pred_out, "uninform/6k/manifest.json")
out_rds       = file.path(stitched, "uninform/stitched_pred_spikeins_6k.rds")
res <- stitch_species_autodetect(runs_dir, manifest_path, expected_len_per_sweep = 15)
saveRDS(res, out_rds)

# 8k
runs_dir      = file.path(pred_out, "uninform/8k/runs")
manifest_path = file.path(pred_out, "uninform/8k/manifest.json")
out_rds       = file.path(stitched, "uninform/stitched_pred_spikeins_8k.rds")
res <- stitch_species_autodetect(runs_dir, manifest_path, expected_len_per_sweep = 15)
saveRDS(res, out_rds)

source(file.path(repo_root, "R/pre-post-processing/stitch_simple_n_models.R"))

# art
runs_dir      = file.path(pred_out, "uninform/art/runs")
manifest_path = file.path(pred_out, "uninform/art/manifest.json")
out_rds       = file.path(stitched, "uninform/stitched_pred_spikeins_art.rds")
res <- stitch_simple_n_models(runs_dir, manifest_path, expected_len_per_sweep = 15)
saveRDS(res, out_rds)

# bio
runs_dir      = file.path(pred_out, "uninform/bio/runs")
manifest_path = file.path(pred_out, "uninform/bio/manifest.json")
out_rds       = file.path(stitched, "uninform/stitched_pred_spikeins_bio.rds")
res <- stitch_simple_n_models(runs_dir, manifest_path, expected_len_per_sweep = 15)
saveRDS(res, out_rds)

# comb
runs_dir      = file.path(pred_out, "uninform/comb/runs")
manifest_path = file.path(pred_out, "uninform/comb/manifest.json")
out_rds       = file.path(stitched, "uninform/stitched_pred_spikeins_comb.rds")
res <- stitch_simple_n_models(runs_dir, manifest_path, expected_len_per_sweep = 15)
saveRDS(res, out_rds)
