# Data_analysis

MATLAB workflow for multi-probe electrophysiology preprocessing, stimulus metadata extraction, DLAG model training, and post hoc latent-variable analysis.

This repository contains script-based analysis code for converting sorted neural recordings into trial-aligned population data, preparing model-ready inputs, training DLAG models, and analyzing the inferred latent structure across stimulus conditions, brain areas, and model variants.

The code is currently organized as a research workflow rather than a packaged MATLAB toolbox. Most scripts are intended to be run after editing the user-parameter section at the top of each file.

---

## Overview

This repository supports a full analysis pipeline:

1. Extract event markers and align SpikeGLX / Kilosort outputs to stimulus runs.
2. Parse Expo XML files into trial-wise stimulus metadata.
3. Segment spikes and optionally LFPs into trial-aligned data.
4. Compute unit-level run metrics and condition-level response statistics.
5. Bin spike trains into count, firing-rate, z-scored, and demeaned response formats.
6. Build DLAG-ready `model_data_allruns.mat` datasets.
7. Train DLAG models either across all conditions or separately by condition.
8. Analyze latent variables, delay directionality, shared variance, RF structure, condition effects, and size effects.
9. Generate trial-shuffled control datasets for comparison against true simultaneously recorded trial structure.

The main target use case is multi-area neural population data analyzed with DLAG-style latent variable models.

---

## Repository structure

```text
Data_analysis/
├── data_seg/
│   ├── generate_lfp_seg_trial.m
│   ├── processing_to_count_and_fr.m
│   └── sorting_check/
│       ├── get_spikt_evt_mrkidx_seg.m
│       └── unit_statistic_by_condition.m
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
        └── plot_size_effect.m