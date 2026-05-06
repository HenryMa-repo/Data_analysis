%% =========================================================================
% model_data_prepar with trial shuffle
%
% Purpose:
%   For selected stimTag runs, prepare model-ready trial data from multiple
%   probe-specific kilosort folders.
%
%   This script now also creates a within-condition trial-shuffled control
%   dataset immediately after the original model_data_allruns dataset is
%   generated and saved.
%
% For each selected stim_tag:
%   1) Select units from each probe using one of two user-defined methods.
%   2) Read eight data types from bined_data_allruns.
%   3) Keep only selected units in each probe.
%   4) Merge probes along the unit dimension for each trial.
%   5) Remove NaN-containing trials or neurons using user-defined strategies.
%   6) Convert each data type from unit x trial x T into a 1 x Ntrial
%      struct array with fields:
%          - trialId
%          - T
%          - y   (unit x T)
%   7) Group each processed data type again by condition and store
%      condition-wise trial struct arrays.
%   8) Save the original model_data_allruns.mat.
%   9) Create a within-condition, group-specific trial-shuffled version.
%  10) Save the shuffled control dataset as a separate .mat file.
%
% Inputs required in each kilosort folder:
%   - unit_run_metrics.mat
%   - unit_condition_metrics.mat
%   - bined_data_allruns.mat
%
% Original output saved in the common CatGT folder:
%   - model_data_allruns.mat
%
% Original output variables:
%   model_data_allruns
%
% Shuffled output saved in the common CatGT folder:
%   - model_data_allruns_trialshuffled_withincondition.mat
%
% Shuffled output variables:
%   model_data_allruns
%   trial_shuffle_withincondition_info
%
% Output structure:
%   model_data_allruns{i}
%       .stim_tag
%       .analysis_window
%       .bin_size
%       .fr_threshold
%       .ff_threshold
%       .unit_selection_method
%       .nan_trial_strategy
%       .groupd
%       .probe0_usedunit_ids
%       .probe1_usedunit_ids
%       ...
%       .condition_fields
%       .condition_index_per_trial_full
%       .conditions_full
%       .n_trials_full
%       .raw_count
%       .raw_fr
%       .z_within_trial
%       .z_within_condition
%       .z_across_conditions
%       .demean_count_within_trial
%       .demean_fr_within_trial
%       .demean_pooledsd_within_condition
%
%   For each of the eight main data fields above, an additional field is
%   also stored:
%       .<data_field>_by_condition
%
%   Each .<data_field>_by_condition is a 1 x Ncond struct array.
%   Each condition struct contains:
%       .condition_index
%       .trial_indices_full
%       .trial_ids_present
%       .n_trials_full
%       .n_trials_present
%       .trials
%       plus all condition identity fields from conditions_full except
%       trial_indices
%
%   If nan_trial_strategy == 4, the following additional fields are stored:
%       .raw_count_nanmask
%       .raw_fr_nanmask
%       .z_within_trial_nanmask
%       .z_within_condition_nanmask
%       .z_across_conditions_nanmask
%       .demean_count_within_trial_nanmask
%       .demean_fr_within_trial_nanmask
%       .demean_pooledsd_within_condition_nanmask
%
%   If nan_trial_strategy == 4, each mask field also has:
%       .<mask_field>_by_condition
%
%   If nan_trial_strategy == 6, the following additional field patterns are
%   stored for each main data field:
%       .<data_field>_groupd
%       .<data_field>_probe0_usedunit_ids
%       .<data_field>_probe1_usedunit_ids
%       ...
%       .<data_field>_kept_neuron_global
%
% Unit selection methods:
%   method 1:
%       keep units with
%           unit_run_metrics.fr_stim >= fr_threshold
%           unit_run_metrics.fano_factor < ff_threshold
%
%   method 2:
%       keep units with
%           unit_condition_metrics.min_fr_stim >= fr_threshold
%           unit_condition_metrics.max_fano_factor < ff_threshold
%
% NaN trial strategies:
%   strategy 1:
%       remove a trial from ALL data types if that trial contains a
%       NaN in ANY one of the merged data types
%
%   strategy 2:
%       remove NaN-containing trials independently for each data type,
%       without forcing the outputs to stay trial-aligned
%
%   strategy 3:
%       do not remove any trials, even if NaN is present
%
%   strategy 4:
%       do not remove any trials; replace NaN with 0 and store a mask
%
%   strategy 5:
%       do not remove any trials; remove any neuron that contains NaN in
%       any one of the merged data types, then update probe-specific
%       used unit IDs and groupd accordingly
%
%   strategy 6:
%       do not remove any trials; for each data type independently,
%       remove neurons that contain NaN in that data type only
%
% Trial-shuffle control:
%   The shuffled dataset is generated from the original model_data_allruns
%   already in memory, after model_data_allruns.mat is saved.
%
%   Within each run and each condition, each neural group receives its own
%   trial permutation. For a target trial t:
%
%       group 1 response can come from source trial perm_group1(t)
%       group 2 response can come from source trial perm_group2(t)
%
%   Therefore, each group keeps its own complete population x time response
%   from real trials, but the trial-by-trial pairing between groups is
%   broken within the same condition.
%
%   The shuffle is applied to all main data versions:
%       raw_count
%       raw_fr
%       z_within_trial
%       z_within_condition
%       z_across_conditions
%       demean_count_within_trial
%       demean_fr_within_trial
%       demean_pooledsd_within_condition
%
%   If nan_trial_strategy == 4, each *_nanmask field follows the same
%   shuffle as its corresponding base data field.
%
%   For nan_trial_strategy == 2 only, different data versions may contain
%   different trialId lists. In that case, the script checks the main data
%   fields:
%
%       - if all main data fields have identical trialId lists, all data
%         versions share the same condition-wise group-specific shuffle;
%
%       - if any main data field has a different trialId list, each data
%         version gets its own condition-wise group-specific shuffle based
%         on its own retained trialId list.
%
%   For all other nan_trial_strategy values, the trial lists are produced
%   by this same model_data_prepar workflow as aligned trial lists, so the
%   script uses a shared shuffle across data versions.
%
%   nTrials <= 1 within a condition is left unchanged silently.
%
%   For nTrials == 2 and two groups, it is impossible for both groups to be
%   non-identity and also different. The script uses swap for group 1 and
%   identity for group 2, so cross-group pairing is still broken.
%
%   By default, trial_shuffle_random_seed = [], which uses rng('shuffle').
%   Set trial_shuffle_random_seed to a numeric scalar for reproducible
%   shuffled output.
%
% Notes:
%   1. All probe folders must belong to the same CatGT folder.
%   2. For a given stim_tag, all probes must have the same analysis_window,
%      bin_size, trial count, bin count, and trial identities.
%   3. For within-trial z-score, NaN from zero variance is common;
%      strategy 4 and 6 are often more suitable if this data type is used.
%   4. The shuffled .mat file intentionally saves the shuffled data using
%      the same variable name, model_data_allruns, so downstream DLAG
%      training scripts can load it without changing variable names.
%   5. trial_shuffle_withincondition_info records the shuffle mode, random
%      seed mode, RNG state, and condition-wise group-specific source trial
%      mappings used to create the shuffled control dataset.
% =========================================================================

clc;
clear;

addpath(genpath(fullfile('.', 'expo_tools')));
addpath(genpath(fullfile('.', 'utils')));

%% ----------------------- User parameters -----------------------

probe_ksDirs = { ...
    'I:\np_data\RafiL001p0120_g1\catgt_RafiL001p0120_g1\RafiL001p0120_g1_imec0\kilosort4_10_dedup_phy', ...
    'I:\np_data\RafiL001p0120_g1\catgt_RafiL001p0120_g1\RafiL001p0120_g1_imec1\kilosort4_2_dedup_phy' ...
};

stimTag = { ...
    '[RFG_coarse2dg_99_4_150isi]', ...
    '[dir12_gpl_2_200isi_fixedphase]', ...
    '_2[Gpl2_2c_2sz_400_2_200isi]'};

fr_threshold = 0.5;
ff_threshold = 5;

% Unit selection method:
%   1 = use unit_run_metrics
%   2 = use unit_condition_metrics
unit_selection_method = 2;

% NaN trial strategies:
%   1 = global aligned deletion across all data types
%   2 = each data type removes its own NaN trials independently
%   3 = do not remove any trials
%   4 = do not remove any trials; replace NaN with 0 and store a mask
%   5 = do not remove any trials; remove neurons globally across all data types
%   6 = do not remove any trials; remove neurons independently for each data type
nan_trial_strategy = 4;

% Trial-shuffle control output:
%   [] or omitted behavior = rng('shuffle') each time this script runs.
%   Numeric scalar, for example 123, gives reproducible shuffle output.
trial_shuffle_random_seed = [];

%% ----------------------- Validate user parameters -----------------------

if ~iscell(probe_ksDirs) || isempty(probe_ksDirs)
    error('probe_ksDirs must be a non-empty cell array of kilosort folder paths.');
end

for p = 1:numel(probe_ksDirs)
    if ~ischar(probe_ksDirs{p}) && ~isstring(probe_ksDirs{p})
        error('probe_ksDirs{%d} must be a char or string path.', p);
    end
    probe_ksDirs{p} = char(probe_ksDirs{p});
end

if ~isscalar(fr_threshold) || ~isnumeric(fr_threshold) || ~isfinite(fr_threshold)
    error('fr_threshold must be a finite numeric scalar.');
end

if ~isscalar(ff_threshold) || ~isnumeric(ff_threshold) || ~isfinite(ff_threshold)
    error('ff_threshold must be a finite numeric scalar.');
end

if ~(isequal(unit_selection_method, 1) || isequal(unit_selection_method, 2))
    error('unit_selection_method must be 1 or 2.');
end

if ~(isequal(nan_trial_strategy, 1) || isequal(nan_trial_strategy, 2) || ...
        isequal(nan_trial_strategy, 3) || isequal(nan_trial_strategy, 4) || ...
        isequal(nan_trial_strategy, 5) || isequal(nan_trial_strategy, 6))
    error('nan_trial_strategy must be 1, 2, 3, 4, 5, or 6.');
end

if ~(isempty(trial_shuffle_random_seed) || ...
        (isscalar(trial_shuffle_random_seed) && isnumeric(trial_shuffle_random_seed) && isfinite(trial_shuffle_random_seed)))
    error('trial_shuffle_random_seed must be [] or a finite numeric scalar.');
end

%% ----------------------- Determine common CatGT folder -----------------------

catgt_folder = get_common_catgt_folder(probe_ksDirs);
fprintf('Common catgt_folder: %s\n', catgt_folder);

%% ----------------------- Process all selected stim tags -----------------------

model_data_allruns = cell(numel(stimTag), 1);

for s = 1:numel(stimTag)
    this_stim_tag = stimTag{s};

    fprintf('\n============================================================\n');
    fprintf('Processing stim_tag: %s\n', this_stim_tag);
    fprintf('============================================================\n');

    probe_data = cell(numel(probe_ksDirs), 1);
    groupd = zeros(1, numel(probe_ksDirs));

    for p = 1:numel(probe_ksDirs)
        ksDir = probe_ksDirs{p};
        fprintf(' Probe %d ksDir: %s\n', p-1, ksDir);

        probe_data{p} = process_one_probe_one_run( ...
            ksDir, this_stim_tag, fr_threshold, ff_threshold, unit_selection_method);

        groupd(p) = numel(probe_data{p}.used_unit_ids);
    end

    ref = probe_data{1};
    for p = 2:numel(probe_ksDirs)
        validate_probe_alignment(ref, probe_data{p}, this_stim_tag, p-1);
    end

    merged = merge_probes_for_one_run(probe_data);

    model_data_allruns{s} = build_model_output_for_one_run( ...
        merged, probe_data, this_stim_tag, fr_threshold, ff_threshold, ...
        unit_selection_method, nan_trial_strategy, groupd);
end

%% ----------------------- Save original output -----------------------

original_output_file = fullfile(catgt_folder, 'model_data_allruns.mat');
save(original_output_file, 'model_data_allruns', '-v7.3');

fprintf('\nSaved original model data:\n');
fprintf(' %s\n', original_output_file);

%% ----------------------- Create and save trial-shuffled output -----------------------

[model_data_allruns_shuffled, trial_shuffle_withincondition_info] = ...
    make_trialshuffled_model_data_allruns( ...
        model_data_allruns, original_output_file, ...
        fullfile(catgt_folder, 'model_data_allruns_trialshuffled_withincondition.mat'), ...
        trial_shuffle_random_seed);

model_data_allruns = model_data_allruns_shuffled;
shuffle_output_file = trial_shuffle_withincondition_info.output_file;

save(shuffle_output_file, 'model_data_allruns', 'trial_shuffle_withincondition_info', '-v7.3');

fprintf('\nSaved trial-shuffled model data:\n');
fprintf(' %s\n', shuffle_output_file);


fprintf('\nDone.\n');

%% ======================= Local functions =======================

function probe_out = process_one_probe_one_run(ksDir, stim_tag, fr_threshold, ff_threshold, unit_selection_method)
if ~isfolder(ksDir)
    error('kilosort folder does not exist: %s', ksDir);
end

run_file = fullfile(ksDir, 'unit_run_metrics.mat');
cond_file = fullfile(ksDir, 'unit_condition_metrics.mat');
bined_file = fullfile(ksDir, 'bined_data_allruns.mat');

if ~isfile(run_file)
    error('Missing file: %s', run_file);
end
if ~isfile(cond_file)
    error('Missing file: %s', cond_file);
end
if ~isfile(bined_file)
    error('Missing file: %s', bined_file);
end

Srun = load(run_file, 'unit_run_metrics');
Scond = load(cond_file, 'unit_condition_metrics');
Sbined = load(bined_file, 'bined_data_allruns');

if ~isfield(Srun, 'unit_run_metrics')
    error('unit_run_metrics not found in %s', run_file);
end
if ~isfield(Scond, 'unit_condition_metrics')
    error('unit_condition_metrics not found in %s', cond_file);
end
if ~isfield(Sbined, 'bined_data_allruns')
    error('bined_data_allruns not found in %s', bined_file);
end

unit_run_metrics = Srun.unit_run_metrics;
unit_condition_metrics = Scond.unit_condition_metrics;
bined_data_allruns = Sbined.bined_data_allruns;

run_idx_in_run_metrics = find_run_index_by_stim_tag(unit_run_metrics, stim_tag);
run_idx_in_cond_metrics = find_run_index_by_stim_tag(unit_condition_metrics, stim_tag);
run_idx_in_bined = find_run_index_by_stim_tag(bined_data_allruns, stim_tag);

run_metrics = unit_run_metrics{run_idx_in_run_metrics};
cond_metrics = unit_condition_metrics{run_idx_in_cond_metrics};
bined_entry = bined_data_allruns{run_idx_in_bined};

[used_unit_ids, ~] = select_units( ...
    run_metrics, cond_metrics, fr_threshold, ff_threshold, unit_selection_method);

if ~isfield(bined_entry, 'unit_ids')
    error('unit_ids missing in bined_data_allruns entry for stim_tag %s', stim_tag);
end

bined_unit_ids = bined_entry.unit_ids(:);
[tf, idx_in_bined] = ismember(used_unit_ids, bined_unit_ids);

if ~all(tf)
    missing_ids = used_unit_ids(~tf);
    error('Selected unit_ids not found in bined_data_allruns for stim_tag %s: %s', ...
        stim_tag, mat2str(missing_ids(:)'));
end

data_fields = get_data_field_list();

probe_out = struct();
probe_out.stim_tag = stim_tag;
probe_out.analysis_window = bined_entry.analysis_window;
probe_out.bin_size = bined_entry.bin_size;
probe_out.bin_edges = bined_entry.bin_edges;
probe_out.bin_centers = bined_entry.bin_centers;
probe_out.unit_ids_all = bined_unit_ids;
probe_out.used_unit_ids = used_unit_ids(:);
probe_out.condition_fields = bined_entry.condition_fields;
probe_out.condition_index_per_trial = bined_entry.condition_index_per_trial(:);
probe_out.conditions = bined_entry.conditions;

for k = 1:numel(data_fields)
    f = data_fields{k};
    if ~isfield(bined_entry, f)
        error('Field %s missing in bined_data_allruns entry for stim_tag %s', f, stim_tag);
    end
    X = bined_entry.(f);
    if size(X, 1) ~= numel(bined_unit_ids)
        error('Field %s has unit dimension mismatch in ksDir %s', f, ksDir);
    end
    probe_out.(f) = X(idx_in_bined, :, :);
end
end

function [used_unit_ids, keep_idx] = select_units(run_metrics, cond_metrics, fr_threshold, ff_threshold, unit_selection_method)
switch unit_selection_method
    case 1
        if ~isfield(run_metrics, 'unit_ids') || ~isfield(run_metrics, 'fr_stim') || ~isfield(run_metrics, 'fano_factor')
            error('unit_run_metrics entry is missing required fields for method 1.');
        end
        unit_ids = run_metrics.unit_ids(:);
        fr_metric = run_metrics.fr_stim(:);
        ff_metric = run_metrics.fano_factor(:);
        keep_idx = isfinite(fr_metric) & isfinite(ff_metric) & ...
            (fr_metric >= fr_threshold) & (ff_metric < ff_threshold);
        used_unit_ids = unit_ids(keep_idx);

    case 2
        if ~isfield(cond_metrics, 'unit_ids') || ~isfield(cond_metrics, 'min_fr_stim') || ~isfield(cond_metrics, 'max_fano_factor')
            error('unit_condition_metrics entry is missing required fields for method 2.');
        end
        unit_ids = cond_metrics.unit_ids(:);
        fr_metric = cond_metrics.min_fr_stim(:);
        ff_metric = cond_metrics.max_fano_factor(:);
        keep_idx = isfinite(fr_metric) & isfinite(ff_metric) & ...
            (fr_metric >= fr_threshold) & (ff_metric < ff_threshold);
        used_unit_ids = unit_ids(keep_idx);

    otherwise
        error('Unknown unit_selection_method.');
end
end

function validate_probe_alignment(ref, cur, stim_tag, probe_index)
if ~isequal(ref.analysis_window, cur.analysis_window)
    error('analysis_window mismatch across probes for stim_tag %s (probe %d).', stim_tag, probe_index);
end
if ~isequal(ref.bin_size, cur.bin_size)
    error('bin_size mismatch across probes for stim_tag %s (probe %d).', stim_tag, probe_index);
end
if ~isequal(ref.bin_edges, cur.bin_edges)
    error('bin_edges mismatch across probes for stim_tag %s (probe %d).', stim_tag, probe_index);
end
if ~isequal(ref.bin_centers, cur.bin_centers)
    error('bin_centers mismatch across probes for stim_tag %s (probe %d).', stim_tag, probe_index);
end
if ~isequal(ref.condition_index_per_trial, cur.condition_index_per_trial)
    error('condition_index_per_trial mismatch across probes for stim_tag %s (probe %d).', stim_tag, probe_index);
end
if numel(ref.conditions) ~= numel(cur.conditions)
    error('Condition count mismatch across probes for stim_tag %s (probe %d).', stim_tag, probe_index);
end

data_fields = get_data_field_list();
for k = 1:numel(data_fields)
    f = data_fields{k};
    Xref = ref.(f);
    Xcur = cur.(f);
    if size(Xref, 2) ~= size(Xcur, 2)
        error('Trial count mismatch in %s across probes for stim_tag %s (probe %d).', f, stim_tag, probe_index);
    end
    if size(Xref, 3) ~= size(Xcur, 3)
        error('Bin count mismatch in %s across probes for stim_tag %s (probe %d).', f, stim_tag, probe_index);
    end
end
end

function merged = merge_probes_for_one_run(probe_data)
data_fields = get_data_field_list();

merged = struct();
merged.stim_tag = probe_data{1}.stim_tag;
merged.analysis_window = probe_data{1}.analysis_window;
merged.bin_size = probe_data{1}.bin_size;
merged.bin_edges = probe_data{1}.bin_edges;
merged.bin_centers = probe_data{1}.bin_centers;
merged.condition_fields = probe_data{1}.condition_fields;
merged.condition_index_per_trial_full = probe_data{1}.condition_index_per_trial;
merged.conditions_full = probe_data{1}.conditions;

for k = 1:numel(data_fields)
    f = data_fields{k};
    Xcat = probe_data{1}.(f);
    for p = 2:numel(probe_data)
        Xcat = cat(1, Xcat, probe_data{p}.(f));
    end
    merged.(f) = Xcat;
end
end

function out = build_model_output_for_one_run(merged, probe_data, stim_tag, fr_threshold, ff_threshold, unit_selection_method, nan_trial_strategy, groupd)
data_fields = get_data_field_list();

out = struct();
out.stim_tag = stim_tag;
out.analysis_window = merged.analysis_window;
out.bin_size = merged.bin_size;
out.bin_edges = merged.bin_edges;
out.bin_centers = merged.bin_centers;
out.fr_threshold = fr_threshold;
out.ff_threshold = ff_threshold;
out.unit_selection_method = unit_selection_method;
out.nan_trial_strategy = nan_trial_strategy;

if nan_trial_strategy ~= 6
    out.groupd = groupd(:)';
    for p = 1:numel(probe_data)
        field_name = sprintf('probe%d_usedunit_ids', p-1);
        out.(field_name) = probe_data{p}.used_unit_ids(:);
    end
end

out.condition_fields = merged.condition_fields;
out.condition_index_per_trial_full = merged.condition_index_per_trial_full;
out.conditions_full = merged.conditions_full;
out.n_trials_full = numel(merged.condition_index_per_trial_full);

switch nan_trial_strategy
    case 1
        bad_trial = false(1, out.n_trials_full);
        for k = 1:numel(data_fields)
            f = data_fields{k};
            bad_trial = bad_trial | get_bad_trial_mask(merged.(f));
        end
        keep_trial_ids = find(~bad_trial);
        out.kept_trial_ids_global = keep_trial_ids(:)';
        for k = 1:numel(data_fields)
            f = data_fields{k};
            out.(f) = build_trial_struct_array(merged.(f), keep_trial_ids);
        end

    case 2
        for k = 1:numel(data_fields)
            f = data_fields{k};
            bad_trial = get_bad_trial_mask(merged.(f));
            keep_trial_ids = find(~bad_trial);
            keep_field_name = sprintf('%s_kept_trial_ids', f);
            out.(keep_field_name) = keep_trial_ids(:)';
            out.(f) = build_trial_struct_array(merged.(f), keep_trial_ids);
        end

    case 3
        keep_trial_ids = 1:out.n_trials_full;
        out.kept_trial_ids_global = keep_trial_ids;
        for k = 1:numel(data_fields)
            f = data_fields{k};
            out.(f) = build_trial_struct_array(merged.(f), keep_trial_ids);
        end

    case 4
        keep_trial_ids = 1:out.n_trials_full;
        out.kept_trial_ids_global = keep_trial_ids;
        for k = 1:numel(data_fields)
            f = data_fields{k};
            X = merged.(f);
            nanmask = isnan(X);
            X_filled = X;
            X_filled(nanmask) = 0;
            out.(f) = build_trial_struct_array(X_filled, keep_trial_ids);
            mask_field_name = sprintf('%s_nanmask', f);
            out.(mask_field_name) = build_trial_struct_array(nanmask, keep_trial_ids);
        end

    case 5
        bad_neuron = false(sum(groupd), 1);
        for k = 1:numel(data_fields)
            f = data_fields{k};
            bad_neuron = bad_neuron | get_bad_neuron_mask(merged.(f));
        end
        keep_neuron = ~bad_neuron;
        if ~any(keep_neuron)
            error('After nan_trial_strategy = 5, no neurons remain for stim_tag %s.', stim_tag);
        end
        [new_probe_unit_ids, new_groupd] = update_probe_unit_ids_after_neuron_removal( ...
            probe_data, groupd, keep_neuron);
        out.groupd = new_groupd(:)';
        for p = 1:numel(new_probe_unit_ids)
            field_name = sprintf('probe%d_usedunit_ids', p-1);
            out.(field_name) = new_probe_unit_ids{p}(:);
        end
        keep_trial_ids = 1:out.n_trials_full;
        out.kept_trial_ids_global = keep_trial_ids;
        out.kept_neuron_global = find(keep_neuron(:))';
        for k = 1:numel(data_fields)
            f = data_fields{k};
            X = merged.(f);
            X = X(keep_neuron, :, :);
            out.(f) = build_trial_struct_array(X, keep_trial_ids);
        end

    case 6
        keep_trial_ids = 1:out.n_trials_full;
        out.kept_trial_ids_global = keep_trial_ids;
        for k = 1:numel(data_fields)
            f = data_fields{k};
            X = merged.(f);
            bad_neuron = get_bad_neuron_mask(X);
            keep_neuron = ~bad_neuron;
            [new_probe_unit_ids, new_groupd] = update_probe_unit_ids_after_neuron_removal( ...
                probe_data, groupd, keep_neuron);
            groupd_field = sprintf('%s_groupd', f);
            out.(groupd_field) = new_groupd(:)';
            for p = 1:numel(new_probe_unit_ids)
                unit_field_name = sprintf('%s_probe%d_usedunit_ids', f, p-1);
                out.(unit_field_name) = new_probe_unit_ids{p}(:);
            end
            kept_neuron_field = sprintf('%s_kept_neuron_global', f);
            out.(kept_neuron_field) = find(keep_neuron(:))';
            X = X(keep_neuron, :, :);
            out.(f) = build_trial_struct_array(X, keep_trial_ids);
        end

    otherwise
        error('Unknown nan_trial_strategy.');
end

out = add_condition_groupings_to_output(out, data_fields);

if nan_trial_strategy == 4
    mask_fields = get_mask_field_list(data_fields);
    out = add_condition_groupings_to_output(out, mask_fields);
end
end

function out = add_condition_groupings_to_output(out, fields_to_group)
for k = 1:numel(fields_to_group)
    f = fields_to_group{k};
    if ~isfield(out, f)
        error('Field %s is missing when trying to add condition groupings.', f);
    end
    by_field_name = sprintf('%s_by_condition', f);
    out.(by_field_name) = build_by_condition_struct( ...
        out.(f), out.conditions_full, out.condition_index_per_trial_full);
end
end

function cond_struct = build_by_condition_struct(trial_struct_array, conditions_full, condition_index_per_trial_full)
nCond = numel(conditions_full);
cond_struct = repmat(struct(), 1, nCond);

present_trial_ids = zeros(1, numel(trial_struct_array));
for i = 1:numel(trial_struct_array)
    present_trial_ids(i) = trial_struct_array(i).trialId;
end

for c = 1:nCond
    if isfield(conditions_full(c), 'trial_indices')
        trial_indices_full = conditions_full(c).trial_indices(:)';
    else
        trial_indices_full = find(condition_index_per_trial_full == c);
    end

    keep_mask = false(1, numel(trial_struct_array));
    for i = 1:numel(trial_struct_array)
        tid = trial_struct_array(i).trialId;
        if tid < 1 || tid > numel(condition_index_per_trial_full)
            error('trialId %d is out of range for condition_index_per_trial_full.', tid);
        end
        keep_mask(i) = (condition_index_per_trial_full(tid) == c);
    end

    trial_ids_present = present_trial_ids(keep_mask);
    trials_present = trial_struct_array(keep_mask);

    cond_struct(c).condition_index = c;
    cond_struct(c).trial_indices_full = trial_indices_full(:)';
    cond_struct(c).trial_ids_present = trial_ids_present(:)';
    cond_struct(c).n_trials_full = numel(trial_indices_full);
    cond_struct(c).n_trials_present = numel(trial_ids_present);
    cond_struct(c).trials = trials_present;

    fn = fieldnames(conditions_full(c));
    for j = 1:numel(fn)
        this_field = fn{j};
        if strcmp(this_field, 'trial_indices')
            continue;
        end
        cond_struct(c).(this_field) = conditions_full(c).(this_field);
    end
end
end

function mask_fields = get_mask_field_list(data_fields)
mask_fields = cell(size(data_fields));
for k = 1:numel(data_fields)
    mask_fields{k} = sprintf('%s_nanmask', data_fields{k});
end
end

function bad_trial = get_bad_trial_mask(X)
bad_trial = squeeze(any(any(isnan(X), 1), 3));
bad_trial = reshape(bad_trial, 1, []);
end

function S = build_trial_struct_array(X, keep_trial_ids)
nKeep = numel(keep_trial_ids);
nUnit = size(X, 1);
nBin = size(X, 3);
S = repmat(struct('trialId', [], 'T', [], 'y', []), 1, nKeep);

for i = 1:nKeep
    tr = keep_trial_ids(i);
    S(i).trialId = tr;
    S(i).T = nBin;
    S(i).y = reshape(X(:, tr, :), nUnit, nBin);
end
end

function idx = find_run_index_by_stim_tag(cell_of_structs, stim_tag)
all_tags = cell(numel(cell_of_structs), 1);
for i = 1:numel(cell_of_structs)
    if ~isfield(cell_of_structs{i}, 'stim_tag')
        error('Entry %d is missing stim_tag.', i);
    end
    all_tags{i} = cell_of_structs{i}.stim_tag;
end

idx = find(strcmp(all_tags, stim_tag));
if isempty(idx)
    error('Requested stim_tag not found: %s', stim_tag);
end
if numel(idx) > 1
    error('Duplicate stim_tag found: %s', stim_tag);
end
end

function fields = get_data_field_list()
fields = { ...
    'raw_count', ...
    'raw_fr', ...
    'z_within_trial', ...
    'z_within_condition', ...
    'z_across_conditions', ...
    'demean_count_within_trial', ...
    'demean_fr_within_trial', ...
    'demean_pooledsd_within_condition'};
end

function catgt_folder = get_common_catgt_folder(probe_ksDirs)
catgt_list = cell(numel(probe_ksDirs), 1);
for i = 1:numel(probe_ksDirs)
    ksDir = probe_ksDirs{i};
    probe_folder = fileparts(ksDir);
    catgt_folder_i = fileparts(probe_folder);
    catgt_list{i} = catgt_folder_i;
end

catgt_folder = catgt_list{1};
for i = 2:numel(catgt_list)
    if ~strcmp(catgt_folder, catgt_list{i})
        error(['All probe kilosort folders must belong to the same CatGT folder.\n' ...
               'Got:\n%s\n%s'], catgt_folder, catgt_list{i});
    end
end
end

function bad_neuron = get_bad_neuron_mask(X)
bad_neuron = squeeze(any(any(isnan(X), 2), 3));
bad_neuron = bad_neuron(:);
end

function [new_probe_unit_ids, new_groupd] = update_probe_unit_ids_after_neuron_removal(probe_data, old_groupd, keep_neuron)
new_probe_unit_ids = cell(numel(probe_data), 1);
new_groupd = zeros(1, numel(probe_data));
row_start = 1;

for p = 1:numel(probe_data)
    row_end = row_start + old_groupd(p) - 1;
    this_keep = keep_neuron(row_start:row_end);
    this_ids = probe_data{p}.used_unit_ids(:);
    new_probe_unit_ids{p} = this_ids(this_keep);
    new_groupd(p) = numel(new_probe_unit_ids{p});
    row_start = row_end + 1;
end
end

%% ======================= Trial-shuffle local functions =======================

function [model_data_allruns_shuffled, info] = make_trialshuffled_model_data_allruns(model_data_allruns, source_file, output_file, random_seed)
if ~iscell(model_data_allruns)
    error('model_data_allruns must be a cell array.');
end

if nargin < 2
    source_file = '';
end
if nargin < 3 || isempty(output_file)
    output_file = fullfile(pwd, 'model_data_allruns_trialshuffled_withincondition.mat');
end
if nargin < 4
    random_seed = [];
end

if isempty(random_seed)
    rng('shuffle');
    random_seed_mode = 'shuffle';
else
    rng(random_seed);
    random_seed_mode = 'fixed';
end
rng_state_after_seed_setup = rng;

data_fields = get_data_field_list();
model_data_allruns_shuffled = model_data_allruns;

info = struct();
info.source_file = source_file;
info.output_file = output_file;
info.random_seed = random_seed;
info.random_seed_mode = random_seed_mode;
info.rng_state_after_seed_setup = rng_state_after_seed_setup;
info.created_by = mfilename;
info.created_on = datestr(now);
info.rule = [ ...
    'Within each run and condition, each group receives its own trial permutation. ', ...
    'The group-specific permutations are shared across data versions unless ', ...
    'nan_trial_strategy == 2 produced non-identical trialId lists across main data fields.'];
info.runs = cell(size(model_data_allruns));

fprintf('\n============================================================\n');
fprintf('Creating within-condition group-specific trial-shuffled control\n');
fprintf('============================================================\n');

for s = 1:numel(model_data_allruns_shuffled)
    this_run = model_data_allruns_shuffled{s};

    if isfield(this_run, 'stim_tag')
        fprintf(' Shuffling run %d/%d: %s\n', s, numel(model_data_allruns_shuffled), this_run.stim_tag);
    else
        fprintf(' Shuffling run %d/%d\n', s, numel(model_data_allruns_shuffled));
    end

    validate_run_has_required_condition_fields_for_shuffle(this_run, s);

    main_fields = data_fields(cellfun(@(f) isfield(this_run, f), data_fields));
    if isempty(main_fields)
        error('No main data fields were found in model_data_allruns{%d}.', s);
    end

    mask_fields = get_existing_mask_fields_for_shuffle(this_run, data_fields);
    fields_to_group_after_shuffle = [main_fields, mask_fields];

    run_info = struct();
    if isfield(this_run, 'stim_tag')
        run_info.stim_tag = this_run.stim_tag;
    else
        run_info.stim_tag = '';
    end
    if isfield(this_run, 'nan_trial_strategy')
        run_info.nan_trial_strategy = this_run.nan_trial_strategy;
    else
        run_info.nan_trial_strategy = [];
    end

    use_shared_shuffle = true;
    trial_id_check = struct();
    trial_id_check.checked = false;
    trial_id_check.all_main_fields_identical = true;
    trial_id_check.reference_field = main_fields{1};
    trial_id_check.different_fields = {};

    if isfield(this_run, 'nan_trial_strategy') && isequal(this_run.nan_trial_strategy, 2)
        trial_id_check.checked = true;
        [all_same, different_fields] = check_main_field_trial_ids_identical(this_run, main_fields);
        trial_id_check.all_main_fields_identical = all_same;
        trial_id_check.different_fields = different_fields;
        if ~all_same
            use_shared_shuffle = false;
        end
    end

    run_info.trial_id_check = trial_id_check;

    if use_shared_shuffle
        run_info.shuffle_mode = 'shared_across_data_versions';
        ref_field = main_fields{1};
        ref_trial_ids = get_trial_ids_for_shuffle(this_run.(ref_field));
        ref_groupd = get_groupd_for_shuffle(this_run, ref_field);
        n_groups = numel(ref_groupd);
        shared_info = build_run_shuffle_info_for_field(this_run, ref_trial_ids, n_groups, ref_field);
        run_info.shared_shuffle = shared_info;

        for k = 1:numel(main_fields)
            f = main_fields{k};
            groupd = get_groupd_for_shuffle(this_run, f);
            this_run.(f) = apply_shuffle_to_trial_struct_for_field( ...
                this_run.(f), groupd, shared_info, f, s);
        end

        for k = 1:numel(mask_fields)
            f = mask_fields{k};
            groupd = get_groupd_for_shuffle(this_run, f);
            this_run.(f) = apply_shuffle_to_trial_struct_for_field( ...
                this_run.(f), groupd, shared_info, f, s);
        end
    else
        run_info.shuffle_mode = 'per_data_version';
        run_info.field_shuffle = struct();

        for k = 1:numel(main_fields)
            f = main_fields{k};
            trial_ids = get_trial_ids_for_shuffle(this_run.(f));
            groupd = get_groupd_for_shuffle(this_run, f);
            field_info = build_run_shuffle_info_for_field(this_run, trial_ids, numel(groupd), f);
            run_info.field_shuffle.(f) = field_info;

            this_run.(f) = apply_shuffle_to_trial_struct_for_field( ...
                this_run.(f), groupd, field_info, f, s);

            mask_field = sprintf('%s_nanmask', f);
            if isfield(this_run, mask_field)
                mask_groupd = get_groupd_for_shuffle(this_run, mask_field);
                this_run.(mask_field) = apply_shuffle_to_trial_struct_for_field( ...
                    this_run.(mask_field), mask_groupd, field_info, mask_field, s);
            end
        end
    end

    this_run = rebuild_by_condition_fields_after_shuffle(this_run, fields_to_group_after_shuffle);

    model_data_allruns_shuffled{s} = this_run;
    info.runs{s} = run_info;
end
end

function fields = get_existing_mask_fields_for_shuffle(run_entry, data_fields)
fields = {};
for k = 1:numel(data_fields)
    f = sprintf('%s_nanmask', data_fields{k});
    if isfield(run_entry, f)
        fields{end+1} = f; %#ok<AGROW>
    end
end
end

function validate_run_has_required_condition_fields_for_shuffle(run_entry, run_index)
if ~isfield(run_entry, 'condition_index_per_trial_full')
    error('condition_index_per_trial_full is missing in model_data_allruns{%d}.', run_index);
end
if ~isfield(run_entry, 'conditions_full')
    error('conditions_full is missing in model_data_allruns{%d}.', run_index);
end
end

function [all_same, different_fields] = check_main_field_trial_ids_identical(run_entry, main_fields)
ref_field = main_fields{1};
ref_ids = get_trial_ids_for_shuffle(run_entry.(ref_field));
all_same = true;
different_fields = {};

for k = 2:numel(main_fields)
    f = main_fields{k};
    ids = get_trial_ids_for_shuffle(run_entry.(f));
    if ~isequal(ids, ref_ids)
        all_same = false;
        different_fields{end+1} = f; %#ok<AGROW>
    end
end
end

function trial_ids = get_trial_ids_for_shuffle(trial_struct_array)
if ~isstruct(trial_struct_array) || ~isfield(trial_struct_array, 'trialId')
    error('Expected a trial struct array with field trialId.');
end

trial_ids = zeros(1, numel(trial_struct_array));
for i = 1:numel(trial_struct_array)
    trial_ids(i) = trial_struct_array(i).trialId;
end
end

function groupd = get_groupd_for_shuffle(run_entry, field_name)
base_field = strip_suffix_for_shuffle(field_name, '_nanmask');
specific_groupd_field = sprintf('%s_groupd', base_field);

if isfield(run_entry, specific_groupd_field)
    groupd = run_entry.(specific_groupd_field);
elseif isfield(run_entry, 'groupd')
    groupd = run_entry.groupd;
else
    error('Cannot find groupd or %s.', specific_groupd_field);
end

groupd = double(groupd(:)');
if any(groupd < 0) || any(~isfinite(groupd)) || any(groupd ~= round(groupd))
    error('Invalid groupd for field %s.', field_name);
end
end

function out = strip_suffix_for_shuffle(in, suffix)
out = in;
if numel(in) >= numel(suffix) && strcmp(in(end-numel(suffix)+1:end), suffix)
    out = in(1:end-numel(suffix));
end
end

function run_shuffle = build_run_shuffle_info_for_field(run_entry, trial_ids, n_groups, field_name)
cond_idx = run_entry.condition_index_per_trial_full(:)';
n_cond = numel(run_entry.conditions_full);

run_shuffle = struct();
run_shuffle.field_name = field_name;
run_shuffle.n_groups = n_groups;
run_shuffle.conditions = repmat(struct( ...
    'condition_index', [], ...
    'trial_ids_present', [], ...
    'n_trials_present', [], ...
    'local_permutation_by_group', [], ...
    'source_trial_ids_by_group', []), 1, n_cond);

for c = 1:n_cond
    trial_ids_this_condition = trial_ids(cond_idx(trial_ids) == c);
    n_trials = numel(trial_ids_this_condition);
    local_perms = make_group_permutations_for_shuffle(n_trials, n_groups);
    source_ids_by_group = cell(1, n_groups);
    for g = 1:n_groups
        source_ids_by_group{g} = trial_ids_this_condition(local_perms{g});
    end

    run_shuffle.conditions(c).condition_index = c;
    run_shuffle.conditions(c).trial_ids_present = trial_ids_this_condition(:)';
    run_shuffle.conditions(c).n_trials_present = n_trials;
    run_shuffle.conditions(c).local_permutation_by_group = local_perms;
    run_shuffle.conditions(c).source_trial_ids_by_group = source_ids_by_group;
end
end

function perms = make_group_permutations_for_shuffle(n_trials, n_groups)
base = 1:n_trials;
perms = cell(1, n_groups);

if n_trials <= 1
    for g = 1:n_groups
        perms{g} = base;
    end
    return;
end

if n_trials == 2
    if n_groups == 1
        perms{1} = [2 1];
    else
        perms{1} = [2 1];
        perms{2} = [1 2];
        for g = 3:n_groups
            if mod(g, 2) == 1
                perms{g} = [2 1];
            else
                perms{g} = [1 2];
            end
        end
    end
    return;
end

perms{1} = random_derangement_for_shuffle(n_trials);
for g = 2:n_groups
    perms{g} = random_derangement_different_from_previous_for_shuffle(n_trials, perms(1:g-1));
end
end

function p = random_derangement_for_shuffle(n_trials)
base = 1:n_trials;
max_attempts = 10000;

for attempt = 1:max_attempts %#ok<NASGU>
    p = randperm(n_trials);
    if all(p ~= base)
        return;
    end
end

p = [2:n_trials 1];
end

function p = random_derangement_different_from_previous_for_shuffle(n_trials, previous_perms)
base = 1:n_trials;
max_attempts = 10000;

for attempt = 1:max_attempts %#ok<NASGU>
    p = randperm(n_trials);
    if any(p == base)
        continue;
    end

    ok = true;
    for j = 1:numel(previous_perms)
        if any(p == previous_perms{j})
            ok = false;
            break;
        end
    end

    if ok
        return;
    end
end

p = n_trials:-1:1;
if any(p == base)
    p = [2:n_trials 1];
end

for j = 1:numel(previous_perms)
    if any(p == previous_perms{j})
        error(['Could not construct a group-specific permutation that is both ', ...
               'non-identity and elementwise different from previous groups for n_trials = %d.'], n_trials);
    end
end
end

function trial_struct_array = apply_shuffle_to_trial_struct_for_field(trial_struct_array, groupd, run_shuffle, field_name, run_index)
if isempty(trial_struct_array)
    return;
end

if numel(groupd) ~= run_shuffle.n_groups
    error(['Group count mismatch while shuffling field %s in model_data_allruns{%d}: ', ...
           'numel(groupd) = %d, expected %d.'], ...
        field_name, run_index, numel(groupd), run_shuffle.n_groups);
end

n_rows_expected = sum(groupd);
n_rows_actual = size(trial_struct_array(1).y, 1);
if n_rows_actual ~= n_rows_expected
    error(['Row count mismatch in model_data_allruns{%d}.%s: ', ...
           'size(y,1) = %d, sum(groupd) = %d.'], ...
        run_index, field_name, n_rows_actual, n_rows_expected);
end

row_ranges = make_group_row_ranges_for_shuffle(groupd);
trial_id_to_position = make_trial_id_to_position_map_for_shuffle(trial_struct_array);
original_trial_struct_array = trial_struct_array;

for c = 1:numel(run_shuffle.conditions)
    target_trial_ids = run_shuffle.conditions(c).trial_ids_present;
    n_trials = numel(target_trial_ids);

    if n_trials <= 1
        continue;
    end

    for i = 1:n_trials
        target_id = target_trial_ids(i);
        target_pos = trial_id_to_position(target_id);
        y_new = original_trial_struct_array(target_pos).y;

        for g = 1:numel(groupd)
            rows = row_ranges{g};
            if isempty(rows)
                continue;
            end

            source_ids = run_shuffle.conditions(c).source_trial_ids_by_group{g};
            source_id = source_ids(i);
            source_pos = trial_id_to_position(source_id);
            y_source = original_trial_struct_array(source_pos).y;

            if size(y_source, 2) ~= size(y_new, 2)
                error(['Time-bin mismatch in model_data_allruns{%d}.%s between ', ...
                       'target trialId %d and source trialId %d.'], ...
                    run_index, field_name, target_id, source_id);
            end

            y_new(rows, :) = y_source(rows, :);
        end

        trial_struct_array(target_pos).y = y_new;
    end
end
end

function row_ranges = make_group_row_ranges_for_shuffle(groupd)
row_ranges = cell(1, numel(groupd));
row_start = 1;

for g = 1:numel(groupd)
    row_end = row_start + groupd(g) - 1;
    if groupd(g) == 0
        row_ranges{g} = [];
    else
        row_ranges{g} = row_start:row_end;
    end
    row_start = row_end + 1;
end
end

function M = make_trial_id_to_position_map_for_shuffle(trial_struct_array)
M = containers.Map('KeyType', 'double', 'ValueType', 'double');

for i = 1:numel(trial_struct_array)
    tid = trial_struct_array(i).trialId;
    if isKey(M, tid)
        error('Duplicate trialId found: %d.', tid);
    end
    M(tid) = i;
end
end

function run_entry = rebuild_by_condition_fields_after_shuffle(run_entry, fields_to_group)
for k = 1:numel(fields_to_group)
    f = fields_to_group{k};
    by_field_name = sprintf('%s_by_condition', f);
    run_entry.(by_field_name) = build_by_condition_struct( ...
        run_entry.(f), run_entry.conditions_full, run_entry.condition_index_per_trial_full);
end
end
