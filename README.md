# Data_analysis

MATLAB scripts for preprocessing sorted electrophysiology data, extracting trial-wise stimulus metadata, preparing DLAG-ready inputs, training DLAG models, and analyzing inferred latent structure.

This repository is organized around a practical end-to-end workflow:

  1. align spikes/events to trials
  2. extract stimulus information from Expo XML
  3. compute run-level and condition-level unit metrics
  4. bin and normalize neural activity
  5. build `model_data_allruns.mat`
  6. train DLAG models
  7. analyze inferred latents, condition effects, and DLAG subspace geometry

* * *

## Overview

This repo is a personal working codebase for multi-probe neural data analysis centered on DLAG (Delayed Latents Across Groups).

It combines:

  * session-level preprocessing from SpikeGLX / CatGT / Kilosort outputs
  * trial-wise stimulus parsing from Expo XML
  * unit-level quality and response summaries
  * construction of DLAG-ready MATLAB data structures
  * DLAG training in pooled or condition-specific modes
  * post hoc analysis of latent reproducibility, latent categories, shared variance, condition dependence, RF/size effects, and subspace similarity

The code is currently organized as script-based workflows, not as a packaged toolbox. In most cases, you should open a script, edit the user parameters near the top, and run it for your session.

* * *

## Repository structure

    data_seg/
      sorting_check/
        get_spikt_evt_mrkidx_seg.m
        unit_statistic_by_condition.m
      generate_lfp_seg_trial.m
      processing_to_count_and_fr.m

    readstim_info/
      stiminfo_pertrial.m

    dlag_mainheng/
      prepar_data/
        model_data_prepar.m
        model_data_prepar_with_trialshuffle.m
        plot_size_rawcount_scatter.m
      train/
        train_dlag.m
        train_dlag_by_condition.m
      analysis/
        plot_dlag_results.m
        Latents_compare.m
        Anova_latents_for_all_conds_used_dlag.m
        RF_analysis.m
        plot_size_effect.m
        subspace_similarity_dlag.m

* * *

## Main workflow

### 1. Preprocess sorted data

Use scripts in `data_seg/` to align sorted spike data and event markers to stimulus runs and trials.

Main scripts:

  * `get_spikt_evt_mrkidx_seg.m` — load sorted spike/event data and create trial-aligned spike structures
  * `unit_statistic_by_condition.m` — compute unit response statistics by stimulus condition
  * `generate_lfp_seg_trial.m` — generate trial-aligned LFP segments
  * `processing_to_count_and_fr.m` — convert trial spike data into count, firing-rate, z-scored, or demeaned formats

### 2. Extract stimulus metadata

Use `readstim_info/stiminfo_pertrial.m` to parse Expo XML files and save trial-wise stimulus information into `stiminfo.mat`.

The exact fields in `stiminfo` depend on the stimulus protocol used in each run.

### 3. Prepare DLAG model input

Use scripts in `dlag_mainheng/prepar_data/` to build model-ready MATLAB data files.

Main scripts:

  * `model_data_prepar.m` — prepare `model_data_allruns.mat` for DLAG training
  * `model_data_prepar_with_trialshuffle.m` — generate trial-shuffled control datasets
  * `plot_size_rawcount_scatter.m` — inspect size-related raw-count effects before modeling

### 4. Train DLAG models

Use scripts in `dlag_mainheng/train/`:

  * `train_dlag.m` — train DLAG models using pooled all-condition data
  * `train_dlag_by_condition.m` — train separate DLAG models for individual conditions

Training outputs are saved under folders such as:

    FA_Dlag_<data_content>/mat_results/runXXX/
    FA_Dlag_<data_content>_condition<condition>/mat_results/runXXX/

### 5. Analyze DLAG results

Use scripts in `dlag_mainheng/analysis/` after model training.

Main scripts:

  * `plot_dlag_results.m` — summarize fitted DLAG models and save best-model outputs
  * `Latents_compare.m` — compare latent categories, DSL, shared variance, and condition-related latent structure
  * `Anova_latents_for_all_conds_used_dlag.m` — test condition effects in inferred latents
  * `RF_analysis.m` — analyze receptive-field related structure
  * `plot_size_effect.m` — summarize size effects
  * `subspace_similarity_dlag.m` — compute within-group DLAG subspace principal angles and directional subspace overlap, including across-vs-within and feedforward-vs-feedback comparisons

* * *

## Notes

This repository is mainly for personal research analysis. File paths, stimulus conventions, and parameter choices are expected to be edited inside each script before running.
