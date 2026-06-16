# Data_analysis

MATLAB scripts for preprocessing sorted electrophysiology data, extracting trial-wise stimulus metadata, preparing DLAG-ready inputs, training DLAG models, and analyzing inferred latent structure.

This repository is a personal working codebase for multi-probe neural population analysis centered on DLAG (Delayed Latents Across Groups). It is organized as a script-based workflow rather than a packaged MATLAB toolbox. In most cases, open the relevant script, edit the user parameters near the top, and run it for the dataset being analyzed.

## Overview

The code supports an end-to-end analysis pipeline:

1. Align sorted spikes and event markers to trials.
2. Extract trial-wise stimulus metadata from Expo XML files.
3. Compute run-level and condition-level unit metrics.
4. Bin spike trains and generate count, firing-rate, z-scored, and demeaned data formats.
5. Merge probe-specific data into DLAG-ready `model_data_allruns.mat`.
6. Train DLAG models in pooled all-condition mode or condition-specific mode.
7. Analyze fitted DLAG models, latent structure, reconstruction quality, size effects, and subspace geometry.

The current workflow is designed for SpikeGLX / CatGT / Kilosort / Phy style outputs, Expo stimulus files, and MATLAB-based DLAG model fitting.

## Repository structure

```text
Data_analysis/
├── data_seg/
│   ├── sorting_check/
│   │   ├── get_spikt_evt_mrkidx_seg.m
│   │   └── unit_statistic_by_condition.m
│   ├── generate_lfp_seg_trial.m
│   └── processing_to_count_and_fr.m
│
├── readstim_info/
│   └── stiminfo_pertrial.m
│
└── dlag_mainheng/
    ├── prepar_data/
    │   ├── model_data_prepar.m
    │   ├── model_data_prepar_with_trialshuffle.m
    │   └── plot_size_rawcount_scatter.m
    │
    ├── train/
    │   ├── train_dlag.m
    │   └── train_dlag_by_condition.m
    │
    └── analysis/
        ├── plot_dlag_results.m
        ├── Latents_compare.m
        ├── Anova_latents_for_all_conds_used_dlag.m
        ├── RF_analysis.m
        ├── plot_size_effect.m
        ├── size_effect_comparation.m
        ├── pick_latents_by_svexp.m
        ├── subspace_similarity_dlag.m
        ├── data_reconstruction.m
        ├── calculate_pooled_conditions_R2.m
        ├── plot_reconstruction_R2.m
        ├── reconstruction_visualization.m
        └── plotCSVEvsDim.m