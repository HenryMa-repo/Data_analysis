%% compute_psd_timescale.m
%
% Compute PSD-based timescale summaries from saved DLAG posterior latent
% trajectories in bestmodel*.mat.
%
% The script reads seqEst(n).xsm from each bestmodel*.mat, computes one
% PSD-based summary for each across-area latent and each within-area latent,
% then saves psd_timescale_stats.mat in the corresponding run folder.
%
% Saved local order for each group:
%   PSD_timescale.local(g).period_ms =
%       [across latents, within-area latents for group g]
%
% This local order matches the latent order used by plot_timescale_compare.

clear;
clc;

%% ======================= USER SETTINGS =================================

% Session folder. If this script is run from the catgt_* session folder,
% leave this as pwd.
session_dir = pwd;

data_content = 'raw_count';
% Common options:
% 'raw_count'
% 'raw_fr'
% 'z_within_trial'
% 'z_within_condition'
% 'z_across_conditions'
% 'demean_count_within_trial'
% 'demean_fr_within_trial'
% 'demean_pooledsd_within_condition'

data_condition = [];
% []   : all-condition model, folder FA_Dlag_<data_content>
% 1:16 : condition-specific models, folders FA_Dlag_<data_content>_conditionN

runIdx = 1;

% Stimulus tag used to find the matching run in model_data_allruns.
stim_tag = '_2[Gpl2_2c_2sz_400_2_200isi]';

% All-condition model PSD averaging mode.
% This option is used only when data_condition = [].
all_condition_psd_average_mode = 'condition_balanced_raw';
% options:
% 'trial_pooled'
% 'condition_balanced_raw'
% 'condition_balanced_norm'

% Welch PSD window length in time bins.
pwelch_window_len_bins = 20;

% Prominent peak threshold for PSD unimodality check.
peak_prominence_frac = 0.01;

save_mat = true;
overwrite_existing = true;

print_full_error_report = false;

%% ======================= VALIDATE SETTINGS ==============================

if ~isfolder(session_dir)
    error('session_dir does not exist:\n%s', session_dir);
end

if isempty(data_content) || ~(ischar(data_content) || isstring(data_content))
    error('data_content must be a character vector or string.');
end
data_content = char(data_content);

if ~(ischar(stim_tag) || isstring(stim_tag))
    error('stim_tag must be a character vector or string.');
end
stim_tag = char(stim_tag);

if ~isempty(data_condition)
    if ~isnumeric(data_condition) || any(~isfinite(data_condition(:)))
        error('data_condition must be [] or a numeric vector.');
    end
    data_condition = reshape(data_condition, 1, []);
end

if ~isscalar(runIdx) || ~isnumeric(runIdx) || ~isfinite(runIdx) || runIdx < 1
    error('runIdx must be a positive scalar.');
end
runIdx = round(runIdx);

valid_average_modes = { ...
    'trial_pooled', ...
    'condition_balanced_raw', ...
    'condition_balanced_norm'};

all_condition_psd_average_mode = char(all_condition_psd_average_mode);

if ~ismember(all_condition_psd_average_mode, valid_average_modes)
    error('Unknown all_condition_psd_average_mode: %s', all_condition_psd_average_mode);
end

if ~isscalar(pwelch_window_len_bins) || ...
        ~isnumeric(pwelch_window_len_bins) || ...
        ~isfinite(pwelch_window_len_bins) || ...
        pwelch_window_len_bins < 2
    error('pwelch_window_len_bins must be a scalar integer >= 2.');
end
pwelch_window_len_bins = round(pwelch_window_len_bins);

if ~isscalar(peak_prominence_frac) || ...
        ~isnumeric(peak_prominence_frac) || ...
        ~isfinite(peak_prominence_frac) || ...
        peak_prominence_frac < 0
    error('peak_prominence_frac must be a non-negative scalar.');
end

if exist('pwelch', 'file') ~= 2
    error('pwelch was not found. This script requires MATLAB Signal Processing Toolbox.');
end

if exist('findpeaks', 'file') ~= 2
    error('findpeaks was not found. This script requires MATLAB Signal Processing Toolbox.');
end

%% ======================= BUILD TARGET LIST ==============================

targets = buildTargetListLocal(session_dir, data_content, data_condition, runIdx);

fprintf('\ncompute_psd_timescale\n');
fprintf('Session dir: %s\n', session_dir);
fprintf('Data content: %s\n', data_content);
fprintf('Run index: %03d\n', runIdx);
fprintf('Stim tag: %s\n', stim_tag);
fprintf('All-condition PSD average mode: %s\n', all_condition_psd_average_mode);
fprintf('Welch window length: %d bins\n', pwelch_window_len_bins);
fprintf('Peak prominence fraction: %.6g\n', peak_prominence_frac);
fprintf('Number of target model folders: %d\n', numel(targets));

%% ======================= PROCESS TARGETS ================================

outputs = struct( ...
    'target_label', {}, ...
    'run_dir', {}, ...
    'bestmodel_file', {}, ...
    'output_file', {}, ...
    'effective_psd_average_mode', {}, ...
    'n_latents_total', {}, ...
    'n_multi_peak', {}, ...
    'frac_multi_peak', {});

skipped = struct( ...
    'target_label', {}, ...
    'run_dir', {}, ...
    'reason', {});

for ti = 1:numel(targets)
    target = targets(ti);

    fprintf('\n[%d/%d] Processing %s\n', ti, numel(targets), target.target_label);
    fprintf('Run folder:\n  %s\n', target.run_dir);

    try
        PSD_timescale = computePsdTimescaleForRunLocal( ...
            target.run_dir, ...
            session_dir, ...
            data_content, ...
            target.model_mode, ...
            target.condition_index, ...
            runIdx, ...
            stim_tag, ...
            all_condition_psd_average_mode, ...
            pwelch_window_len_bins, ...
            peak_prominence_frac);

        output_file = fullfile(target.run_dir, 'psd_timescale_stats.mat');

        if exist(output_file, 'file') && ~overwrite_existing
            error('Output file already exists and overwrite_existing=false:\n%s', output_file);
        end

        if save_mat
            save(output_file, 'PSD_timescale', '-v7.3');
            fprintf('Saved:\n  %s\n', output_file);
        end

        U = PSD_timescale.unimodality;

        fprintf('PSD multi-peak check: %d / %d latents = %.2f%% have >1 prominent peak.\n', ...
            U.n_multi_peak, U.n_latents_total, 100 * U.frac_multi_peak);

        fprintf('Effective PSD average mode: %s\n', ...
            PSD_timescale.meta.effective_psd_average_mode);

        if ~isempty(PSD_timescale.meta.condition_ids_used)
            fprintf('Condition ids used: [%s]\n', ...
                num2str(PSD_timescale.meta.condition_ids_used));
            fprintf('Trials per condition: [%s]\n', ...
                num2str(PSD_timescale.meta.n_trials_per_condition));
        end

        fprintf('Frequency info: Fs = %.6g Hz, Nyquist = %.6g Hz, grid spacing = %.6g Hz, rough resolution = %.6g Hz.\n', ...
            PSD_timescale.meta.Fs_hz, ...
            PSD_timescale.meta.nyquist_hz, ...
            PSD_timescale.meta.freq_grid_spacing_hz, ...
            PSD_timescale.meta.rough_freq_resolution_hz);

        fprintf('Welch segments per trial: %d\n', ...
            PSD_timescale.meta.pwelch_n_segments_per_trial);

        outputs(end+1).target_label = target.target_label; %#ok<SAGROW>
        outputs(end).run_dir = target.run_dir;
        outputs(end).bestmodel_file = PSD_timescale.meta.bestmodel_file;
        outputs(end).output_file = output_file;
        outputs(end).effective_psd_average_mode = PSD_timescale.meta.effective_psd_average_mode;
        outputs(end).n_latents_total = U.n_latents_total;
        outputs(end).n_multi_peak = U.n_multi_peak;
        outputs(end).frac_multi_peak = U.frac_multi_peak;

    catch ME
        skipped(end+1).target_label = target.target_label; %#ok<SAGROW>
        skipped(end).run_dir = target.run_dir;
        skipped(end).reason = ME.message;

        warning('Skipped %s\nReason: %s', target.target_label, ME.message);

        if print_full_error_report
            fprintf('\nFull error report:\n');
            fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
        end
    end
end

%% ======================= SUMMARY PRINT ==================================

fprintf('\nDone.\n');
fprintf('Completed targets: %d\n', numel(outputs));
fprintf('Skipped targets: %d\n', numel(skipped));

if ~isempty(outputs)
    fprintf('\nCompleted output files:\n');
    for i = 1:numel(outputs)
        fprintf('  %s  [%s]\n    %s\n', ...
            outputs(i).target_label, ...
            outputs(i).effective_psd_average_mode, ...
            outputs(i).output_file);
    end
end

if ~isempty(skipped)
    fprintf('\nSkipped targets:\n');
    for i = 1:numel(skipped)
        fprintf('  %s\n    %s\n', skipped(i).target_label, skipped(i).reason);
    end
end

%% ========================================================================
% Local functions
% ========================================================================

function targets = buildTargetListLocal(session_dir, data_content, data_condition, runIdx)

targets = struct( ...
    'target_label', {}, ...
    'model_mode', {}, ...
    'condition_index', {}, ...
    'model_dir', {}, ...
    'run_dir', {});

if isempty(data_condition)
    model_dir = fullfile(session_dir, sprintf('FA_Dlag_%s', data_content));
    run_dir = fullfile(model_dir, 'mat_results', sprintf('run%03d', runIdx));

    targets(1).target_label = 'all_condition_model';
    targets(1).model_mode = 'all_condition_model';
    targets(1).condition_index = [];
    targets(1).model_dir = model_dir;
    targets(1).run_dir = run_dir;
else
    for i = 1:numel(data_condition)
        c = round(data_condition(i));

        model_dir = fullfile(session_dir, ...
            sprintf('FA_Dlag_%s_condition%d', data_content, c));
        run_dir = fullfile(model_dir, 'mat_results', sprintf('run%03d', runIdx));

        targets(i).target_label = sprintf('condition%d', c);
        targets(i).model_mode = 'condition_specific_models';
        targets(i).condition_index = c;
        targets(i).model_dir = model_dir;
        targets(i).run_dir = run_dir;
    end
end
end

function PSD_timescale = computePsdTimescaleForRunLocal( ...
    run_dir, session_dir, data_content, model_mode, condition_index, runIdx, ...
    stim_tag, all_condition_psd_average_mode, pwelch_window_len_bins, ...
    peak_prominence_frac)

if ~isfolder(run_dir)
    error('Run folder does not exist:\n%s', run_dir);
end

bestmodel_file = findBestmodelFileLocal(run_dir);

S = load(bestmodel_file);

if ~isfield(S, 'bestModel')
    error('bestModel not found in:\n%s', bestmodel_file);
end

if ~isfield(S, 'seqEst')
    error('seqEst not found in:\n%s\nRun plot_dlag_results first or save seqEst in bestmodel*.mat.', ...
        bestmodel_file);
end

bestModel = S.bestModel;
seqEst = S.seqEst;

if isempty(seqEst) || ~isstruct(seqEst) || ~isfield(seqEst, 'xsm')
    error('seqEst is empty or does not contain seqEst(n).xsm in:\n%s', bestmodel_file);
end

xDim_across = getXDimAcrossLocal(bestModel);
xDim_within = getXDimWithinLocal(bestModel);

numGroups = numel(xDim_within);

bin_width_ms = getBinWidthMsLocal(S);
Fs_hz = 1000 / bin_width_ms;

[latent_layout, X_first] = detectLatentLayoutLocal(seqEst, xDim_across, xDim_within);
n_time_bins = size(X_first, 2);
n_trials = numel(seqEst);

if pwelch_window_len_bins > n_time_bins
    error('pwelch_window_len_bins = %d exceeds number of time bins = %d.', ...
        pwelch_window_len_bins, n_time_bins);
end

[effective_psd_average_mode, condition_info] = prepareAveragingInfoLocal( ...
    session_dir, model_mode, seqEst, stim_tag, all_condition_psd_average_mode);

pwelch_overlap_bins = floor(pwelch_window_len_bins / 2);
nfft = max(256, 2^nextpow2(pwelch_window_len_bins));
nyquist_hz = Fs_hz / 2;
freq_grid_spacing_hz = Fs_hz / nfft;
rough_freq_resolution_hz = Fs_hz / pwelch_window_len_bins;

pwelch_n_segments = computeWelchSegmentCountLocal( ...
    n_time_bins, pwelch_window_len_bins, pwelch_overlap_bins);

fprintf('Bestmodel:\n  %s\n', bestmodel_file);
fprintf('Latent layout: %s\n', latent_layout.layout_type);
fprintf('Latent dims: xDim_across = %d, xDim_within = [%s]\n', ...
    xDim_across, num2str(xDim_within));
fprintf('seqEst: %d trials, %d time bins, bin width %.6g ms\n', ...
    n_trials, n_time_bins, bin_width_ms);
fprintf('PSD average mode used for this target: %s\n', effective_psd_average_mode);

win = makeHanningWindowLocal(pwelch_window_len_bins);

across_rows_group1 = latent_layout.across_rows{1};

across_stats = computeRowsPsdLocal( ...
    seqEst, across_rows_group1, latent_layout.total_latent_dim, ...
    Fs_hz, win, pwelch_overlap_bins, nfft, peak_prominence_frac, ...
    effective_psd_average_mode, condition_info.condition_index_per_trial);

freq_hz = across_stats.freq_hz;

within_stats = cell(1, numGroups);

for g = 1:numGroups
    within_stats{g} = computeRowsPsdLocal( ...
        seqEst, latent_layout.within_rows{g}, latent_layout.total_latent_dim, ...
        Fs_hz, win, pwelch_overlap_bins, nfft, peak_prominence_frac, ...
        effective_psd_average_mode, condition_info.condition_index_per_trial);

    if isempty(freq_hz) && ~isempty(within_stats{g}.freq_hz)
        freq_hz = within_stats{g}.freq_hz;
    end
end

if isempty(freq_hz)
    error('No latent trajectories were available for PSD calculation.');
end

local_stats = struct( ...
    'center_frequency_hz', cell(1, numGroups), ...
    'period_ms', cell(1, numGroups), ...
    'n_prominent_peaks', cell(1, numGroups), ...
    'is_multi_peak', cell(1, numGroups));

for g = 1:numGroups
    local_stats(g).center_frequency_hz = [ ...
        across_stats.center_frequency_hz, ...
        within_stats{g}.center_frequency_hz];

    local_stats(g).period_ms = [ ...
        across_stats.period_ms, ...
        within_stats{g}.period_ms];

    local_stats(g).n_prominent_peaks = [ ...
        across_stats.n_prominent_peaks, ...
        within_stats{g}.n_prominent_peaks];

    local_stats(g).is_multi_peak = [ ...
        across_stats.is_multi_peak, ...
        within_stats{g}.is_multi_peak];
end

unimodality = summarizeUnimodalityLocal(across_stats, within_stats, peak_prominence_frac);

PSD_timescale = struct();

PSD_timescale.meta = struct();
PSD_timescale.meta.created_by = 'compute_psd_timescale.m';
PSD_timescale.meta.created_on = datestr(now);
PSD_timescale.meta.session_dir = session_dir;
PSD_timescale.meta.data_content = char(data_content);
PSD_timescale.meta.model_mode = char(model_mode);
PSD_timescale.meta.condition_index = condition_index;
PSD_timescale.meta.runIdx = runIdx;
PSD_timescale.meta.stim_tag = stim_tag;
PSD_timescale.meta.bestmodel_file = bestmodel_file;
PSD_timescale.meta.run_dir = run_dir;

PSD_timescale.meta.timescale_source = 'psd-timescale';
PSD_timescale.meta.latent_field = 'xsm';
PSD_timescale.meta.across_latent_psd_source = 'group1 across rows in seqEst(n).xsm';
PSD_timescale.meta.local_latent_order = '[across latents, within-area latents for this group]';

PSD_timescale.meta.requested_all_condition_psd_average_mode = all_condition_psd_average_mode;
PSD_timescale.meta.effective_psd_average_mode = effective_psd_average_mode;
PSD_timescale.meta.condition_index_source = condition_info.source;
PSD_timescale.meta.condition_ids_used = condition_info.condition_ids;
PSD_timescale.meta.n_trials_per_condition = condition_info.n_trials_per_condition;
PSD_timescale.meta.condition_full = condition_info.condition_full;

PSD_timescale.meta.bin_width_ms = bin_width_ms;
PSD_timescale.meta.Fs_hz = Fs_hz;
PSD_timescale.meta.nyquist_hz = nyquist_hz;
PSD_timescale.meta.n_trials = n_trials;
PSD_timescale.meta.n_time_bins = n_time_bins;
PSD_timescale.meta.total_window_duration_sec = n_time_bins * bin_width_ms / 1000;

PSD_timescale.meta.pwelch_window_len_bins = pwelch_window_len_bins;
PSD_timescale.meta.pwelch_window_len_sec = pwelch_window_len_bins * bin_width_ms / 1000;
PSD_timescale.meta.pwelch_overlap_bins = pwelch_overlap_bins;
PSD_timescale.meta.pwelch_overlap_sec = pwelch_overlap_bins * bin_width_ms / 1000;
PSD_timescale.meta.pwelch_n_segments_per_trial = pwelch_n_segments;
PSD_timescale.meta.nfft = nfft;
PSD_timescale.meta.freq_grid_spacing_hz = freq_grid_spacing_hz;
PSD_timescale.meta.rough_freq_resolution_hz = rough_freq_resolution_hz;

PSD_timescale.meta.psd_method = ...
    'pwelch is computed for each trial; trial PSDs are averaged according to effective_psd_average_mode; aggregate PSD is normalized to area 1.';
PSD_timescale.meta.center_frequency_definition = ...
    'center_frequency_hz = trapz(freq_hz, freq_hz .* psd_norm) / trapz(freq_hz, psd_norm)';
PSD_timescale.meta.period_ms_definition = 'period_ms = 1000 / center_frequency_hz';

PSD_timescale.meta.peak_prominence_frac = peak_prominence_frac;
PSD_timescale.meta.multi_peak_definition = ...
    'is_multi_peak is true when more than one PSD peak has prominence >= peak_prominence_frac * max(prominence).';

PSD_timescale.dimension_info = struct();
PSD_timescale.dimension_info.xDim_across = xDim_across;
PSD_timescale.dimension_info.xDim_within = xDim_within;
PSD_timescale.dimension_info.numGroups = numGroups;
PSD_timescale.dimension_info.total_latent_dim = latent_layout.total_latent_dim;
PSD_timescale.dimension_info.layout_type = latent_layout.layout_type;
PSD_timescale.dimension_info.across_rows = latent_layout.across_rows;
PSD_timescale.dimension_info.within_rows = latent_layout.within_rows;
PSD_timescale.dimension_info.local_rows = latent_layout.local_rows;

PSD_timescale.freq_hz = freq_hz;

PSD_timescale.across = across_stats;
PSD_timescale.within = within_stats;
PSD_timescale.local = local_stats;
PSD_timescale.unimodality = unimodality;

end

function [effective_mode, condition_info] = prepareAveragingInfoLocal( ...
    session_dir, model_mode, seqEst, stim_tag, all_condition_psd_average_mode)

condition_info = emptyConditionInfoLocal();

switch char(model_mode)
    case 'all_condition_model'
        effective_mode = char(all_condition_psd_average_mode);

        if strcmp(effective_mode, 'trial_pooled')
            return;
        end

        condition_info = loadConditionInfoFromModelDataLocal(session_dir, stim_tag, seqEst);

    case 'condition_specific_models'
        effective_mode = 'trial_pooled';

    otherwise
        error('Unknown model_mode: %s', model_mode);
end

end

function condition_info = emptyConditionInfoLocal()

condition_info = struct();
condition_info.condition_index_per_trial = [];
condition_info.condition_ids = [];
condition_info.n_trials_per_condition = [];
condition_info.condition_full = [];
condition_info.source = '';

end

function condition_info = loadConditionInfoFromModelDataLocal(session_dir, stim_tag, seqEst)

[model_data_allruns, model_data_file] = loadModelDataAllrunsLocal(session_dir);

all_run_tags = getAllRunTagsLocal(model_data_allruns);

if isempty(stim_tag)
    if numel(model_data_allruns) == 1
        run_idx = 1;
    else
        error('stim_tag is empty and model_data_allruns contains %d runs.', ...
            numel(model_data_allruns));
    end
else
    run_idx = find(strcmp(all_run_tags, stim_tag));
end

if isempty(run_idx)
    fprintf('\nAvailable stim_tag values in model_data_allruns:\n');
    for i = 1:numel(all_run_tags)
        fprintf('  %d: %s\n', i, all_run_tags{i});
    end
    error('Requested stim_tag not found: %s', stim_tag);
end

if numel(run_idx) > 1
    error('Duplicate stim_tag found: %s', stim_tag);
end

run_data = model_data_allruns{run_idx};

if ~isfield(run_data, 'conditions_full')
    error('model_data_allruns{%d} does not contain conditions_full.', run_idx);
end

condition_full = run_data.conditions_full;

trial_condition_ids = extractTrialConditionIdsFromConditionFullLocal(condition_full, seqEst);

condition_ids = unique(trial_condition_ids(isfinite(trial_condition_ids)));
condition_ids = reshape(condition_ids, 1, []);

n_trials_per_condition = nan(1, numel(condition_ids));

for ci = 1:numel(condition_ids)
    n_trials_per_condition(ci) = sum(trial_condition_ids == condition_ids(ci));
end

condition_info = struct();
condition_info.condition_index_per_trial = reshape(trial_condition_ids, [], 1);
condition_info.condition_ids = condition_ids;
condition_info.n_trials_per_condition = n_trials_per_condition;
condition_info.condition_full = condition_full;
condition_info.source = sprintf('%s:model_data_allruns{%d}.conditions_full.trial_indices', ...
    model_data_file, run_idx);

fprintf('Loaded trial condition ids from:\n  %s\n', condition_info.source);

end

function [model_data_allruns, model_data_file] = loadModelDataAllrunsLocal(session_dir)

candidate_files = { ...
    fullfile(session_dir, 'model_data_allruns.mat'), ...
    fullfile(session_dir, 'model_data_allruns')};

model_data_file = '';

for i = 1:numel(candidate_files)
    if exist(candidate_files{i}, 'file')
        model_data_file = candidate_files{i};
        break;
    end
end

if isempty(model_data_file)
    d = dir(fullfile(session_dir, '*model_data_allruns*.mat'));
    if ~isempty(d)
        model_data_file = fullfile(d(1).folder, d(1).name);
    end
end

if isempty(model_data_file)
    error('Could not find model_data_allruns or model_data_allruns.mat in:\n%s', session_dir);
end

S = load(model_data_file, 'model_data_allruns');

if ~isfield(S, 'model_data_allruns')
    error('model_data_allruns variable not found in:\n%s', model_data_file);
end

model_data_allruns = S.model_data_allruns;

if ~iscell(model_data_allruns) || isempty(model_data_allruns)
    error('model_data_allruns must be a non-empty cell array.');
end

end

function all_tags = getAllRunTagsLocal(model_data_allruns)

all_tags = cell(numel(model_data_allruns), 1);

for j = 1:numel(model_data_allruns)
    if ~isfield(model_data_allruns{j}, 'stim_tag')
        error('stim_tag missing in model_data_allruns{%d}.', j);
    end
    all_tags{j} = model_data_allruns{j}.stim_tag;
end

end

function trial_condition_ids = extractTrialConditionIdsFromConditionFullLocal(condition_full, seqEst)

if isempty(condition_full)
    error('condition_full is empty.');
end

maxTrialIndex = 0;

for k = 1:numel(condition_full)
    if ~isfield(condition_full(k), 'trial_indices')
        error('condition_full(%d) missing field trial_indices.', k);
    end

    idx = condition_full(k).trial_indices(:);

    if ~isempty(idx)
        maxTrialIndex = max(maxTrialIndex, max(idx));
    end
end

if maxTrialIndex < 1
    error('No valid trial indices found in condition_full.trial_indices.');
end

trial_to_condition = nan(1, maxTrialIndex);

for condID = 1:numel(condition_full)
    idx = condition_full(condID).trial_indices(:)';

    if isempty(idx)
        continue;
    end

    if any(idx < 1) || any(mod(idx, 1) ~= 0)
        error('condition_full(%d).trial_indices contains invalid entries.', condID);
    end

    alreadyAssigned = ~isnan(trial_to_condition(idx));

    if any(alreadyAssigned)
        dupIdx = idx(find(alreadyAssigned, 1));
        error('Trial index %d appears in multiple conditions.', dupIdx);
    end

    trial_to_condition(idx) = condID;
end

if isfield(seqEst, 'trialId')
    seq_trial_ids = [seqEst.trialId];
else
    seq_trial_ids = 1:numel(seqEst);
end

seq_trial_ids = double(seq_trial_ids(:)');

if any(seq_trial_ids < 1) || any(mod(seq_trial_ids, 1) ~= 0)
    error('seqEst trial IDs are invalid.');
end

if max(seq_trial_ids) > numel(trial_to_condition)
    error('seqEst contains trial IDs larger than max trial index in condition_full.trial_indices.');
end

trial_condition_ids = trial_to_condition(seq_trial_ids);

if any(isnan(trial_condition_ids))
    missingTrial = seq_trial_ids(find(isnan(trial_condition_ids), 1));
    error('Could not map seqEst trial %d to any condition using condition_full.trial_indices.', ...
        missingTrial);
end

trial_condition_ids = reshape(trial_condition_ids, 1, []);

end

function bestmodel_file = findBestmodelFileLocal(run_dir)

files = dir(fullfile(run_dir, 'bestmodel*.mat'));

if isempty(files)
    files = dir(fullfile(run_dir, '*bestmodel*.mat'));
end

if isempty(files)
    error('No bestmodel*.mat file found in:\n%s', run_dir);
end

[~, idx] = max([files.datenum]);
bestmodel_file = fullfile(files(idx).folder, files(idx).name);

end

function xDim_across = getXDimAcrossLocal(bestModel)

candidate_fields = {'xDim_across', 'xDimAcross'};

for i = 1:numel(candidate_fields)
    f = candidate_fields{i};

    if isfield(bestModel, f)
        xDim_across = double(bestModel.(f));
        xDim_across = xDim_across(:)';

        if isempty(xDim_across)
            xDim_across = 0;
        end

        if numel(xDim_across) ~= 1
            error('bestModel.%s must be scalar.', f);
        end

        xDim_across = round(xDim_across);
        return;
    end
end

if isfield(bestModel, 'gp_params') && isfield(bestModel.gp_params, 'tau_across')
    xDim_across = numel(bestModel.gp_params.tau_across);
    return;
end

error('Could not determine xDim_across from bestModel.');

end

function xDim_within = getXDimWithinLocal(bestModel)

candidate_fields = {'xDim_within', 'xDimWithin'};

for i = 1:numel(candidate_fields)
    f = candidate_fields{i};

    if isfield(bestModel, f)
        xDim_within = bestModel.(f);
        xDim_within = normalizeXDimWithinLocal(xDim_within);
        return;
    end
end

if isfield(bestModel, 'gp_params') && isfield(bestModel.gp_params, 'tau_within')
    tau_within = bestModel.gp_params.tau_within;

    if iscell(tau_within)
        xDim_within = cellfun(@numel, tau_within);
    else
        error('bestModel.gp_params.tau_within exists but is not a cell array.');
    end

    xDim_within = reshape(round(double(xDim_within)), 1, []);
    return;
end

error('Could not determine xDim_within from bestModel.');

end

function xDim_within = normalizeXDimWithinLocal(xDim_within)

if iscell(xDim_within)
    xDim_within = cellfun(@double, xDim_within);
end

xDim_within = double(xDim_within);
xDim_within = reshape(xDim_within, 1, []);
xDim_within = round(xDim_within);

if any(~isfinite(xDim_within)) || any(xDim_within < 0)
    error('xDim_within must contain non-negative finite values.');
end

end

function bin_width_ms = getBinWidthMsLocal(S)

candidate_values = {};

if isfield(S, 'res') && isstruct(S.res)
    candidate_values = addCandidateFieldLocal(candidate_values, S.res, 'binWidth');
    candidate_values = addCandidateFieldLocal(candidate_values, S.res, 'bin_width');
    candidate_values = addCandidateFieldLocal(candidate_values, S.res, 'binWidth_ms');
    candidate_values = addCandidateFieldLocal(candidate_values, S.res, 'bin_width_ms');
end

if isfield(S, 'bestModel') && isstruct(S.bestModel)
    candidate_values = addCandidateFieldLocal(candidate_values, S.bestModel, 'binWidth');
    candidate_values = addCandidateFieldLocal(candidate_values, S.bestModel, 'bin_width');
    candidate_values = addCandidateFieldLocal(candidate_values, S.bestModel, 'binWidth_ms');
    candidate_values = addCandidateFieldLocal(candidate_values, S.bestModel, 'bin_width_ms');
end

candidate_values = addCandidateFieldLocal(candidate_values, S, 'binWidth');
candidate_values = addCandidateFieldLocal(candidate_values, S, 'bin_width');
candidate_values = addCandidateFieldLocal(candidate_values, S, 'binWidth_ms');
candidate_values = addCandidateFieldLocal(candidate_values, S, 'bin_width_ms');

bin_width_ms = [];

for i = 1:numel(candidate_values)
    v = candidate_values{i};

    if isempty(v) || ~isnumeric(v) || ~isscalar(v) || ~isfinite(v) || v <= 0
        continue;
    end

    v = double(v);

    if v < 1
        bin_width_ms = v * 1000;
    else
        bin_width_ms = v;
    end

    break;
end

if isempty(bin_width_ms)
    error(['Could not determine bin width from bestmodel file. ', ...
        'Expected a field like res.binWidth or bestModel.binWidth.']);
end

end

function candidate_values = addCandidateFieldLocal(candidate_values, S, field_name)

if isstruct(S) && isfield(S, field_name)
    candidate_values{end+1} = S.(field_name); %#ok<AGROW>
end

end

function [latent_layout, X_first] = detectLatentLayoutLocal(seqEst, xDim_across, xDim_within)

numGroups = numel(xDim_within);

group_block_dim = numGroups * xDim_across + sum(xDim_within);
unique_across_dim = xDim_across + sum(xDim_within);

X_raw = seqEst(1).xsm;

if ~isnumeric(X_raw) || ndims(X_raw) ~= 2
    error('seqEst(1).xsm must be a numeric 2D matrix.');
end

row_dim = size(X_raw, 1);
col_dim = size(X_raw, 2);

layout_type = '';
total_latent_dim = NaN;
transpose_needed = false;

if row_dim == group_block_dim
    layout_type = 'group_blocks';
    total_latent_dim = group_block_dim;
elseif row_dim == unique_across_dim
    layout_type = 'unique_across';
    total_latent_dim = unique_across_dim;
elseif col_dim == group_block_dim
    layout_type = 'group_blocks';
    total_latent_dim = group_block_dim;
    transpose_needed = true;
elseif col_dim == unique_across_dim
    layout_type = 'unique_across';
    total_latent_dim = unique_across_dim;
    transpose_needed = true;
else
    error(['Could not match seqEst(1).xsm dimensions to expected DLAG latent dimensions.\n', ...
        'size(xsm) = [%d %d]\n', ...
        'Expected group-block total = %d or unique-across total = %d.'], ...
        row_dim, col_dim, group_block_dim, unique_across_dim);
end

if transpose_needed
    X_first = X_raw';
else
    X_first = X_raw;
end

latent_layout = makeLatentRowInfoLocal(layout_type, total_latent_dim, xDim_across, xDim_within);

end

function latent_layout = makeLatentRowInfoLocal(layout_type, total_latent_dim, xDim_across, xDim_within)

numGroups = numel(xDim_within);

across_rows = cell(1, numGroups);
within_rows = cell(1, numGroups);
local_rows = cell(1, numGroups);

switch char(layout_type)
    case 'group_blocks'
        pos = 1;

        for g = 1:numGroups
            across_rows{g} = pos:(pos + xDim_across - 1);
            pos = pos + xDim_across;

            within_rows{g} = pos:(pos + xDim_within(g) - 1);
            pos = pos + xDim_within(g);

            local_rows{g} = [across_rows{g}, within_rows{g}];
        end

    case 'unique_across'
        shared_across_rows = 1:xDim_across;
        pos = xDim_across + 1;

        for g = 1:numGroups
            across_rows{g} = shared_across_rows;

            within_rows{g} = pos:(pos + xDim_within(g) - 1);
            pos = pos + xDim_within(g);

            local_rows{g} = [across_rows{g}, within_rows{g}];
        end

    otherwise
        error('Unknown latent layout type: %s', layout_type);
end

latent_layout = struct();
latent_layout.layout_type = layout_type;
latent_layout.total_latent_dim = total_latent_dim;
latent_layout.across_rows = across_rows;
latent_layout.within_rows = within_rows;
latent_layout.local_rows = local_rows;

end

function stats = computeRowsPsdLocal( ...
    seqEst, row_list, total_latent_dim, Fs_hz, win, overlap_bins, nfft, ...
    peak_prominence_frac, psd_average_mode, condition_index_per_trial)

row_list = reshape(row_list, 1, []);
row_list = row_list(isfinite(row_list));

nRows = numel(row_list);

stats = emptyPsdStatsLocal(nRows);

if nRows == 0
    return;
end

for j = 1:nRows
    row_idx = row_list(j);

    latent_trials = collectLatentTrialsLocal(seqEst, row_idx, total_latent_dim);

    one = computeOneLatentPsdLocal( ...
        latent_trials, Fs_hz, win, overlap_bins, nfft, ...
        peak_prominence_frac, psd_average_mode, condition_index_per_trial);

    if j == 1
        nFreq = numel(one.freq_hz);

        stats.freq_hz = one.freq_hz;
        stats.mean_psd = nan(nFreq, nRows);
        stats.psd_norm = nan(nFreq, nRows);
        stats.center_frequency_hz = nan(1, nRows);
        stats.period_ms = nan(1, nRows);
        stats.total_power_before_normalization = nan(1, nRows);
        stats.n_prominent_peaks = nan(1, nRows);
        stats.is_multi_peak = false(1, nRows);
        stats.peak_freq_hz = cell(1, nRows);
        stats.peak_prominence = cell(1, nRows);
        stats.row_index = row_list;
        stats.psd_average_mode = psd_average_mode;

        if ~isempty(one.condition_ids)
            nCond = numel(one.condition_ids);

            stats.condition_ids = one.condition_ids;
            stats.condition_mean_psd_raw = nan(nFreq, nCond, nRows);
            stats.condition_psd_norm = nan(nFreq, nCond, nRows);
            stats.condition_center_frequency_hz = nan(nCond, nRows);
            stats.condition_period_ms = nan(nCond, nRows);
            stats.condition_total_power = nan(nCond, nRows);
            stats.condition_n_trials = nan(nCond, nRows);
            stats.n_conditions_used_for_average = nan(1, nRows);
        end
    end

    stats.mean_psd(:, j) = one.mean_psd;
    stats.psd_norm(:, j) = one.psd_norm;
    stats.center_frequency_hz(j) = one.center_frequency_hz;
    stats.period_ms(j) = one.period_ms;
    stats.total_power_before_normalization(j) = one.total_power_before_normalization;
    stats.n_prominent_peaks(j) = one.n_prominent_peaks;
    stats.is_multi_peak(j) = one.is_multi_peak;
    stats.peak_freq_hz{j} = one.peak_freq_hz;
    stats.peak_prominence{j} = one.peak_prominence;

    if ~isempty(one.condition_ids)
        if isempty(stats.condition_ids) || ...
                numel(stats.condition_ids) ~= numel(one.condition_ids) || ...
                any(stats.condition_ids ~= one.condition_ids)
            error('Condition ids differ across latent rows.');
        end

        stats.condition_mean_psd_raw(:, :, j) = one.condition_mean_psd_raw;
        stats.condition_psd_norm(:, :, j) = one.condition_psd_norm;
        stats.condition_center_frequency_hz(:, j) = one.condition_center_frequency_hz(:);
        stats.condition_period_ms(:, j) = one.condition_period_ms(:);
        stats.condition_total_power(:, j) = one.condition_total_power(:);
        stats.condition_n_trials(:, j) = one.condition_n_trials(:);
        stats.n_conditions_used_for_average(j) = one.n_conditions_used_for_average;
    end
end

end

function stats = emptyPsdStatsLocal(nRows)

stats = struct();
stats.freq_hz = [];
stats.row_index = nan(1, nRows);
stats.center_frequency_hz = nan(1, nRows);
stats.period_ms = nan(1, nRows);
stats.total_power_before_normalization = nan(1, nRows);
stats.mean_psd = [];
stats.psd_norm = [];
stats.n_prominent_peaks = nan(1, nRows);
stats.is_multi_peak = false(1, nRows);
stats.peak_freq_hz = cell(1, nRows);
stats.peak_prominence = cell(1, nRows);
stats.psd_average_mode = '';

stats.condition_ids = [];
stats.condition_mean_psd_raw = [];
stats.condition_psd_norm = [];
stats.condition_center_frequency_hz = [];
stats.condition_period_ms = [];
stats.condition_total_power = [];
stats.condition_n_trials = [];
stats.n_conditions_used_for_average = [];

end

function latent_trials = collectLatentTrialsLocal(seqEst, row_idx, total_latent_dim)

nTrials = numel(seqEst);

X0 = standardizeLatentMatrixLocal(seqEst(1).xsm, total_latent_dim);
nTime = size(X0, 2);

latent_trials = nan(nTrials, nTime);

for n = 1:nTrials
    if ~isfield(seqEst(n), 'xsm')
        error('seqEst(%d) does not contain xsm.', n);
    end

    X = standardizeLatentMatrixLocal(seqEst(n).xsm, total_latent_dim);

    if size(X, 2) ~= nTime
        error('seqEst(%d).xsm has %d time bins, expected %d.', ...
            n, size(X, 2), nTime);
    end

    if row_idx < 1 || row_idx > size(X, 1)
        error('Requested latent row %d exceeds xsm row count %d.', ...
            row_idx, size(X, 1));
    end

    latent_trials(n, :) = double(X(row_idx, :));
end

end

function X = standardizeLatentMatrixLocal(X_raw, total_latent_dim)

if size(X_raw, 1) == total_latent_dim
    X = double(X_raw);
elseif size(X_raw, 2) == total_latent_dim
    X = double(X_raw');
else
    error('xsm size [%d %d] does not match total latent dim %d.', ...
        size(X_raw, 1), size(X_raw, 2), total_latent_dim);
end

end

function one = computeOneLatentPsdLocal( ...
    latent_trials, Fs_hz, win, overlap_bins, nfft, ...
    peak_prominence_frac, psd_average_mode, condition_index_per_trial)

nTrials = size(latent_trials, 1);

pxx_trials = [];
freq_hz = [];

valid_count = 0;
valid_trial_indices = nan(nTrials, 1);

for tr = 1:nTrials
    x = latent_trials(tr, :);
    x = double(x(:));

    if any(~isfinite(x))
        continue;
    end

    [pxx, f] = pwelch(x, win, overlap_bins, nfft, Fs_hz);

    pxx = double(pxx(:));
    f = double(f(:));

    if isempty(pxx_trials)
        pxx_trials = nan(numel(pxx), nTrials);
        freq_hz = f;
    else
        if numel(pxx) ~= size(pxx_trials, 1) || any(abs(f - freq_hz) > 1e-12)
            error('pwelch returned inconsistent frequency grids across trials.');
        end
    end

    valid_count = valid_count + 1;
    pxx_trials(:, valid_count) = pxx;
    valid_trial_indices(valid_count) = tr;
end

if valid_count == 0
    error('No valid trials available for PSD calculation.');
end

pxx_trials = pxx_trials(:, 1:valid_count);
valid_trial_indices = valid_trial_indices(1:valid_count);

[mean_psd, condition_summary] = averagePsdAcrossTrialsLocal( ...
    pxx_trials, freq_hz, valid_trial_indices, psd_average_mode, condition_index_per_trial);

[psd_norm, total_power, center_frequency_hz, period_ms] = ...
    normalizeAndSummarizePsdLocal(freq_hz, mean_psd);

[n_prominent_peaks, is_multi_peak, peak_freq_hz, peak_prominence] = ...
    countProminentPeaksLocal(freq_hz, psd_norm, peak_prominence_frac);

one = struct();
one.freq_hz = freq_hz;
one.mean_psd = mean_psd;
one.psd_norm = psd_norm;
one.center_frequency_hz = center_frequency_hz;
one.period_ms = period_ms;
one.total_power_before_normalization = total_power;
one.n_prominent_peaks = n_prominent_peaks;
one.is_multi_peak = is_multi_peak;
one.peak_freq_hz = peak_freq_hz;
one.peak_prominence = peak_prominence;

one.condition_ids = condition_summary.condition_ids;
one.condition_mean_psd_raw = condition_summary.condition_mean_psd_raw;
one.condition_psd_norm = condition_summary.condition_psd_norm;
one.condition_center_frequency_hz = condition_summary.condition_center_frequency_hz;
one.condition_period_ms = condition_summary.condition_period_ms;
one.condition_total_power = condition_summary.condition_total_power;
one.condition_n_trials = condition_summary.condition_n_trials;
one.n_conditions_used_for_average = condition_summary.n_conditions_used_for_average;

end

function [mean_psd, condition_summary] = averagePsdAcrossTrialsLocal( ...
    pxx_trials, freq_hz, valid_trial_indices, psd_average_mode, condition_index_per_trial)

condition_summary = emptyConditionSummaryLocal();

switch char(psd_average_mode)
    case 'trial_pooled'
        mean_psd = mean(pxx_trials, 2, 'omitnan');

    case {'condition_balanced_raw', 'condition_balanced_norm'}
        if isempty(condition_index_per_trial)
            error('condition_index_per_trial is required for %s.', psd_average_mode);
        end

        condition_index_per_trial = double(condition_index_per_trial(:));

        if numel(condition_index_per_trial) < max(valid_trial_indices)
            error('condition_index_per_trial is shorter than the number of seqEst trials.');
        end

        labels = condition_index_per_trial(valid_trial_indices);
        labels = labels(:);

        if any(~isfinite(labels))
            error('condition_index_per_trial contains non-finite labels for valid trials.');
        end

        condition_ids = unique(labels(isfinite(labels)));
        condition_ids = reshape(condition_ids, 1, []);

        nFreq = size(pxx_trials, 1);
        nCond = numel(condition_ids);

        condition_mean_psd_raw = nan(nFreq, nCond);
        condition_psd_norm = nan(nFreq, nCond);
        condition_center_frequency_hz = nan(nCond, 1);
        condition_period_ms = nan(nCond, 1);
        condition_total_power = nan(nCond, 1);
        condition_n_trials = nan(nCond, 1);

        psd_for_condition_average = nan(nFreq, nCond);

        for ci = 1:nCond
            c = condition_ids(ci);
            idx = labels == c;

            condition_n_trials(ci) = sum(idx);

            raw_mean = mean(pxx_trials(:, idx), 2, 'omitnan');
            condition_mean_psd_raw(:, ci) = raw_mean;

            [cond_psd_norm, cond_power, cond_cf, cond_period] = ...
                normalizeAndSummarizePsdLocal(freq_hz, raw_mean);

            condition_psd_norm(:, ci) = cond_psd_norm;
            condition_total_power(ci) = cond_power;
            condition_center_frequency_hz(ci) = cond_cf;
            condition_period_ms(ci) = cond_period;

            switch char(psd_average_mode)
                case 'condition_balanced_raw'
                    psd_for_condition_average(:, ci) = raw_mean;

                case 'condition_balanced_norm'
                    psd_for_condition_average(:, ci) = cond_psd_norm;
            end
        end

        mean_psd = mean(psd_for_condition_average, 2, 'omitnan');

        condition_summary.condition_ids = condition_ids;
        condition_summary.condition_mean_psd_raw = condition_mean_psd_raw;
        condition_summary.condition_psd_norm = condition_psd_norm;
        condition_summary.condition_center_frequency_hz = condition_center_frequency_hz;
        condition_summary.condition_period_ms = condition_period_ms;
        condition_summary.condition_total_power = condition_total_power;
        condition_summary.condition_n_trials = condition_n_trials;
        condition_summary.n_conditions_used_for_average = ...
            sum(any(isfinite(psd_for_condition_average), 1));

    otherwise
        error('Unknown psd_average_mode: %s', psd_average_mode);
end

end

function condition_summary = emptyConditionSummaryLocal()

condition_summary = struct();
condition_summary.condition_ids = [];
condition_summary.condition_mean_psd_raw = [];
condition_summary.condition_psd_norm = [];
condition_summary.condition_center_frequency_hz = [];
condition_summary.condition_period_ms = [];
condition_summary.condition_total_power = [];
condition_summary.condition_n_trials = [];
condition_summary.n_conditions_used_for_average = [];

end

function [psd_norm, total_power, center_frequency_hz, period_ms] = ...
    normalizeAndSummarizePsdLocal(freq_hz, psd_in)

psd_in = double(psd_in(:));
freq_hz = double(freq_hz(:));

total_power = trapz(freq_hz, psd_in);

if ~isfinite(total_power) || total_power <= 0
    psd_norm = nan(size(psd_in));
    center_frequency_hz = NaN;
    period_ms = NaN;
    return;
end

psd_norm = psd_in ./ total_power;

denom = trapz(freq_hz, psd_norm);

if ~isfinite(denom) || denom <= 0
    center_frequency_hz = NaN;
    period_ms = NaN;
    return;
end

center_frequency_hz = trapz(freq_hz, freq_hz .* psd_norm) ./ denom;

if isfinite(center_frequency_hz) && center_frequency_hz > 0
    period_ms = 1000 ./ center_frequency_hz;
elseif center_frequency_hz == 0
    period_ms = Inf;
else
    period_ms = NaN;
end

end

function [n_prominent_peaks, is_multi_peak, peak_freq_hz, peak_prominence] = ...
    countProminentPeaksLocal(freq_hz, psd_norm, peak_prominence_frac)

n_prominent_peaks = 0;
is_multi_peak = false;
peak_freq_hz = [];
peak_prominence = [];

if isempty(freq_hz) || isempty(psd_norm)
    return;
end

if any(~isfinite(freq_hz)) || any(~isfinite(psd_norm))
    return;
end

try
    [~, locs, ~, prominences] = findpeaks(psd_norm, freq_hz);
catch
    return;
end

if isempty(prominences)
    return;
end

max_prominence = max(prominences);

if ~isfinite(max_prominence) || max_prominence <= 0
    return;
end

keep = prominences >= peak_prominence_frac * max_prominence;

peak_freq_hz = reshape(locs(keep), 1, []);
peak_prominence = reshape(prominences(keep), 1, []);

n_prominent_peaks = numel(peak_freq_hz);
is_multi_peak = n_prominent_peaks > 1;

end

function unimodality = summarizeUnimodalityLocal(across_stats, within_stats, peak_prominence_frac)

all_n_peaks = [];

if ~isempty(across_stats.n_prominent_peaks)
    all_n_peaks = [all_n_peaks, across_stats.n_prominent_peaks]; %#ok<AGROW>
end

for g = 1:numel(within_stats)
    if ~isempty(within_stats{g}.n_prominent_peaks)
        all_n_peaks = [all_n_peaks, within_stats{g}.n_prominent_peaks]; %#ok<AGROW>
    end
end

all_n_peaks = all_n_peaks(isfinite(all_n_peaks));

n_latents_total = numel(all_n_peaks);
n_multi_peak = sum(all_n_peaks > 1);

if n_latents_total > 0
    frac_multi_peak = n_multi_peak / n_latents_total;
else
    frac_multi_peak = NaN;
end

unimodality = struct();
unimodality.n_latents_total = n_latents_total;
unimodality.n_multi_peak = n_multi_peak;
unimodality.frac_multi_peak = frac_multi_peak;
unimodality.peak_prominence_frac = peak_prominence_frac;
unimodality.n_prominent_peaks_all_latents = all_n_peaks;
unimodality.is_multi_peak_all_latents = all_n_peaks > 1;

end

function win = makeHanningWindowLocal(n)

if exist('hann', 'file') == 2
    try
        win = hann(n);
        return;
    catch
    end
end

if exist('hanning', 'file') == 2
    win = hanning(n);
    return;
end

idx = (0:(n-1))';
win = 0.5 - 0.5 .* cos(2 .* pi .* idx ./ max(1, n - 1));

end

function nSegments = computeWelchSegmentCountLocal(nTime, winLen, overlap)

if nTime < winLen
    nSegments = 0;
    return;
end

step = winLen - overlap;

if step <= 0
    nSegments = 1;
else
    nSegments = 1 + floor((nTime - winLen) / step);
end

end