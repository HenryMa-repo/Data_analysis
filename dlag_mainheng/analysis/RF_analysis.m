%% =========================================================================
% RF_analysis
%
% Purpose:
% First-step RF analysis based on bined_data_allruns.mat.
% Optional second-step target overlay analysis based on model_data_allruns.
%
% Step 1:
% 1. Finds the RF run in bined_data_allruns by rf_stim_tag.
% 2. Sums raw_count across bins to get one response per unit per trial.
% 3. Reconstructs trial-wise xPos, yPos, and stimsize from conditions.
% 4. Computes RF maps: unit x y-position x x-position.
% 5. Fits a 2D Gaussian to each unit RF map.
% 6. Plots RF maps sorted by descending depth.
% 7. Saves unit_rf_results.mat and unit_rf_map_depth_desc.png.
%
% Step 2, only if target_stim_tag is not empty:
% 1. Loads model_data_allruns.
% 2. Finds the target run by target_stim_tag.
% 3. Extracts target stimulus x/y/size combinations.
% 4. Extracts target-selected unit ids for each probe.
% 5. Keeps units that are both selected in target run and pass RF R2 threshold.
% 6. Plots per-probe RF fit centers/sizes over target stimulus layout.
%    Two versions are saved:
%       readable: target layout + RF centers define the view window.
%       full    : target layout + full RF extents define the view window.
% 7. Plots combined probe-average RF over target stimulus layout.
%    Two versions are saved: readable and full.
% 8. Saves target summary mat in cat folder.
%
% Required input in each kilosort folder:
% - bined_data_allruns.mat
%
% Optional input for Step 2:
% - model_data_allruns.mat
%
% Main output saved in each kilosort folder:
% - unit_rf_results.mat
% - unit_rf_map_depth_desc.png
% - unit_rf_target_overlay_readable_<target_safe_name>.png, if target_stim_tag is not empty
% - unit_rf_target_overlay_full_<target_safe_name>.png, if target_stim_tag is not empty
%
% Target output saved in cat folder:
% - unit_rf_target_summary_<target_safe_name>.mat
% - unit_rf_target_combined_readable_<target_safe_name>.png
% - unit_rf_target_combined_full_<target_safe_name>.png
%
% Output variable in unit_rf_results.mat:
% unit_rf.rf_stim_tag
% unit_rf.unit_ids
% unit_rf.unit_depth_um
% unit_rf.unit_channel
% unit_rf.response_count_trial
% unit_rf.xPos_trial
% unit_rf.yPos_trial
% unit_rf.stimsize_trial
% unit_rf.rfs.map
% unit_rf.rfs.x
% unit_rf.rfs.y
% unit_rf.rfstimsize
% unit_rf.fit.center
% unit_rf.fit.size
% unit_rf.fit.rsquare
% unit_rf.fit.params
% unit_rf.plot.depth_desc_order
%
% Notes:
% 1. raw_count can be unit x trial x bin, or unit x trial if only one bin.
% 2. RF response is defined as sum(raw_count, 3), or raw_count itself for
%    single-bin data.
% 3. RF map collapses across any non-position condition dimensions.
%    For example, if ori is also present, trials with the same xPos/yPos
%    are averaged together regardless of ori.
% 4. Depth sorting for plotting:
%    larger depth first, then original unit index.
%    This matches the convention that depth 0 is at the probe tip.
% =========================================================================

clc;
clear;

addpath(genpath(fullfile('.', 'expo_tools')));
addpath(genpath(fullfile('.', 'utils')));

%% ----------------------- User parameters -----------------------

root_folder = 'I:\np_data';
runName = 'RafiL001p0120';
runind = 1;              % run index after -g
probes = [0,1];          % probe indices after -prb

% RF run used for computing unit RF.
rf_stim_tag = '[RFG_coarse2dg_99_4_150isi]';

% Target run for Step 2.
% If [], Step 2 is skipped completely.
target_stim_tag = '_2[Gpl2_2c_2sz_400_2_200isi]';

% model_data_allruns file for Step 2.
% Only used when target_stim_tag is not empty.
dat_file = 'I:\np_data\RafiL001p0120_g1\catgt_RafiL001p0120_g1\model_data_allruns.mat';

% RF R2 threshold for Step 2.
RF_R2_threshold = 0.5;

% If target conditions do not contain a size/stimsize field, this default
% value will be used. Keep [] to force an error when target size is missing.
target_default_stimsize = [];

% Output names for Step 1.
rf_mat_name = 'unit_rf_results.mat';
rf_png_name = 'unit_rf_map_depth_desc.png';

% Probe colors for Step 2.
% Row 1 is probe 0, row 2 is probe 1.
probe_colors = [
    0.0000, 0.4470, 0.7410;
    0.8500, 0.3250, 0.0980
];

%% ----------------------- Build shared session paths -----------------------

run_g = sprintf('%s_g%d', runName, runind);
destDir = fullfile(root_folder, run_g);
cat_folder = fullfile(destDir, ['catgt_' run_g]);

fprintf('destDir    : %s\n', destDir);
fprintf('cat_folder : %s\n', cat_folder);
fprintf('RF stim tag: %s\n', rf_stim_tag);

%% ----------------------- Optional Step 2 setup -----------------------

target_enabled = ~isempty(target_stim_tag);

if target_enabled
    fprintf('\nTarget analysis enabled.\n');
    fprintf('Target stim tag: %s\n', target_stim_tag);
    fprintf('Reading from %s\n', dat_file);

    M = load(dat_file, 'model_data_allruns');
    if ~isfield(M, 'model_data_allruns')
        error('model_data_allruns not found in %s', dat_file);
    end
    model_data_allruns = M.model_data_allruns;

    all_run_tags = get_all_run_tags(model_data_allruns);
    target_run_idx = find(strcmp(all_run_tags, target_stim_tag));

    if isempty(target_run_idx)
        error('Requested target_stim_tag not found: %s', target_stim_tag);
    end
    if numel(target_run_idx) > 1
        error('Duplicate target_stim_tag found: %s', target_stim_tag);
    end

    target_run = model_data_allruns{target_run_idx};
    target_stimulus_xy_size = extract_target_stimulus_xy_size( ...
        target_run, target_default_stimsize);
    target_safe_name = make_filename_safe(target_stim_tag);

    target_rf_summary = struct();
    target_rf_summary.rf_stim_tag = rf_stim_tag;
    target_rf_summary.target_stim_tag = target_stim_tag;
    target_rf_summary.target_safe_name = target_safe_name;
    target_rf_summary.RF_R2_threshold = RF_R2_threshold;
    target_rf_summary.model_data_file = dat_file;
    target_rf_summary.target_run_idx = target_run_idx;
    target_rf_summary.target_stimulus_xy_size = target_stimulus_xy_size;
    target_rf_summary.probe = [];
else
    fprintf('\nTarget analysis disabled because target_stim_tag is empty.\n');
    target_run = [];
    target_stimulus_xy_size = [];
    target_safe_name = '';
    target_rf_summary = [];
end

%% ----------------------- Process each probe folder -----------------------

for ip = 1:numel(probes)
    thisProbe = probes(ip);
    imecStr = sprintf('imec%d', thisProbe);
    probe_folder = fullfile(cat_folder, [run_g '_' imecStr]);

    fprintf('\n============================================================\n');
    fprintf('Processing probe %d\n', thisProbe);
    fprintf('probe_folder: %s\n', probe_folder);
    fprintf('============================================================\n');

    if ~isfolder(probe_folder)
        warning('probe_folder does not exist, skipping probe %d: %s', ...
            thisProbe, probe_folder);
        continue;
    end

    %% ----------------------- Find kilosort folders -----------------------

    d = dir(fullfile(probe_folder, 'kilosort*'));
    d = d([d.isdir]);
    if isempty(d)
        warning('No kilosort* folders found under probe %d: %s', ...
            thisProbe, probe_folder);
        continue;
    end

    [~, idx] = sort(lower({d.name}));
    d = d(idx);

    fprintf('Found %d kilosort folder(s) under probe %d.\n', numel(d), thisProbe);

    %% ----------------------- Process each kilosort folder -----------------------

    for i = 1:numel(d)
        ksDir = fullfile(d(i).folder, d(i).name);
        fprintf('\nProcessing probe %d, ksDir: %s\n', thisProbe, ksDir);

        try
            bined_file = fullfile(ksDir, 'bined_data_allruns.mat');
            if ~isfile(bined_file)
                error('Missing file: %s', bined_file);
            end

            S = load(bined_file, 'bined_data_allruns');
            if ~isfield(S, 'bined_data_allruns')
                error('bined_data_allruns not found in %s', bined_file);
            end
            bined_data_allruns = S.bined_data_allruns;

            rf_idx = find_bined_run_by_stim_tag(bined_data_allruns, rf_stim_tag);
            if isempty(rf_idx)
                error('RF stim_tag not found in bined_data_allruns: %s', rf_stim_tag);
            end
            if numel(rf_idx) > 1
                error('Duplicate RF stim_tag found in bined_data_allruns: %s', rf_stim_tag);
            end

            rf_data = bined_data_allruns{rf_idx};

            %% ----------------------- Step 1: compute RF -----------------------

            unit_rf = compute_unit_rf_from_bined_run(rf_data, rf_stim_tag);

            mat_file = fullfile(ksDir, rf_mat_name);
            png_file = fullfile(ksDir, rf_png_name);

            save(mat_file, 'unit_rf');
            plot_rf_map_depth_desc(unit_rf, png_file);

            fprintf('Saved Step 1:\n');
            fprintf('  %s\n', mat_file);
            fprintf('  %s\n', png_file);

            %% ----------------------- Step 2: target overlay -----------------------

            if target_enabled
                target_result = compute_target_rf_result( ...
                    unit_rf, target_run, thisProbe, ksDir, ...
                    target_stimulus_xy_size, RF_R2_threshold);

                this_color = get_probe_color(thisProbe, probe_colors);

                target_png_readable_name = sprintf('unit_rf_target_overlay_readable_%s.png', target_safe_name);
                target_png_readable_file = fullfile(ksDir, target_png_readable_name);
                plot_target_rf_overlay_per_probe( ...
                    target_stimulus_xy_size, ...
                    target_result.good_center, ...
                    target_result.good_size, ...
                    this_color, ...
                    target_png_readable_file, ...
                    sprintf('probe %d, target %s, readable', thisProbe, target_stim_tag), ...
                    'readable');

                target_png_full_name = sprintf('unit_rf_target_overlay_full_%s.png', target_safe_name);
                target_png_full_file = fullfile(ksDir, target_png_full_name);
                plot_target_rf_overlay_per_probe( ...
                    target_stimulus_xy_size, ...
                    target_result.good_center, ...
                    target_result.good_size, ...
                    this_color, ...
                    target_png_full_file, ...
                    sprintf('probe %d, target %s, full', thisProbe, target_stim_tag), ...
                    'full');

                % Backward-compatible field points to the readable figure.
                target_result.target_overlay_png = target_png_readable_file;
                target_result.target_overlay_readable_png = target_png_readable_file;
                target_result.target_overlay_full_png = target_png_full_file;

                if isempty(target_rf_summary.probe)
                    target_rf_summary.probe = target_result;
                else
                    target_rf_summary.probe(end+1) = target_result; %#ok<SAGROW>
                end

                fprintf('Saved Step 2 per-probe overlays:\n');
                fprintf('  %s\n', target_png_readable_file);
                fprintf('  %s\n', target_png_full_file);
            end

        catch ME
            fprintf(2, 'Error in probe %d, ksDir %s\n', thisProbe, ksDir);
            fprintf(2, '%s\n', ME.message);
        end
    end
end

%% ----------------------- Save Step 2 cat-folder outputs -----------------------

if target_enabled
    target_mat_name = sprintf('unit_rf_target_summary_%s.mat', target_safe_name);
    target_mat_file = fullfile(cat_folder, target_mat_name);

    target_combined_readable_png_name = sprintf('unit_rf_target_combined_readable_%s.png', target_safe_name);
    target_combined_readable_png_file = fullfile(cat_folder, target_combined_readable_png_name);

    target_combined_full_png_name = sprintf('unit_rf_target_combined_full_%s.png', target_safe_name);
    target_combined_full_png_file = fullfile(cat_folder, target_combined_full_png_name);

    save(target_mat_file, 'target_rf_summary');

    plot_target_rf_combined( ...
        target_rf_summary, probe_colors, target_combined_readable_png_file, 'readable');

    plot_target_rf_combined( ...
        target_rf_summary, probe_colors, target_combined_full_png_file, 'full');

    fprintf('\nSaved Step 2 cat-folder outputs:\n');
    fprintf('  %s\n', target_mat_file);
    fprintf('  %s\n', target_combined_readable_png_file);
    fprintf('  %s\n', target_combined_full_png_file);
end

fprintf('\nDone.\n');

%% ======================= Local functions =======================

function rf_idx = find_bined_run_by_stim_tag(bined_data_allruns, rf_stim_tag)
%% =========================================================================
% Find the entry in bined_data_allruns whose .stim_tag matches rf_stim_tag.
% =========================================================================

rf_idx = [];

if ~iscell(bined_data_allruns)
    error('bined_data_allruns must be a cell array.');
end

for r = 1:numel(bined_data_allruns)
    if isempty(bined_data_allruns{r})
        continue;
    end
    if ~isstruct(bined_data_allruns{r})
        continue;
    end
    if ~isfield(bined_data_allruns{r}, 'stim_tag')
        continue;
    end

    if strcmp(bined_data_allruns{r}.stim_tag, rf_stim_tag)
        rf_idx(end+1) = r; %#ok<AGROW>
    end
end
end

function unit_rf = compute_unit_rf_from_bined_run(rf_data, rf_stim_tag)
%% =========================================================================
% Compute unit RF maps and Gaussian fits from one RF run in bined_data_allruns.
%
% Supports:
% raw_count = unit x trial x bin
% raw_count = unit x trial when there is only one bin
% =========================================================================

validate_rf_bined_data(rf_data);

unit_ids = rf_data.unit_ids(:);
raw_count = double(rf_data.raw_count);
nUnit = size(raw_count, 1);
nTrial = size(raw_count, 2);

if numel(unit_ids) ~= nUnit
    error('numel(unit_ids) does not match size(raw_count,1).');
end

if ismatrix(raw_count)
    response_count_trial = raw_count;
else
    response_count_trial = sum(raw_count, 3);
    response_count_trial = reshape(response_count_trial, [nUnit, nTrial]);
end

unit_depth_um = get_optional_unit_vector(rf_data, 'unit_depth_um', nUnit);
unit_channel = get_optional_unit_vector(rf_data, 'unit_channel', nUnit);

[xPos_trial, yPos_trial, stimsize_trial] = ...
    get_trial_position_from_conditions(rf_data.conditions, ...
    rf_data.condition_index_per_trial, ...
    nTrial);

rfstimsize = get_single_stimsize(stimsize_trial);

[RFmap, unique_x, unique_y] = compute_rf_map( ...
    response_count_trial, xPos_trial, yPos_trial);

fit_results = fitGaussianHeatmaps(RFmap, unique_x, unique_y, rfstimsize);

unit_rf = struct();
unit_rf.rf_stim_tag = rf_stim_tag;
unit_rf.unit_ids = unit_ids;
unit_rf.unit_depth_um = unit_depth_um;
unit_rf.unit_channel = unit_channel;
unit_rf.response_count_trial = response_count_trial;
unit_rf.xPos_trial = xPos_trial;
unit_rf.yPos_trial = yPos_trial;
unit_rf.stimsize_trial = stimsize_trial;
unit_rf.rfs = struct();
unit_rf.rfs.map = RFmap;
unit_rf.rfs.x = unique_x;
unit_rf.rfs.y = unique_y;
unit_rf.rfstimsize = rfstimsize;
unit_rf.fit = fit_results;
unit_rf.plot.depth_desc_order = get_depth_desc_order(unit_depth_um, nUnit);
end

function validate_rf_bined_data(rf_data)
%% =========================================================================
% Validate required fields for RF analysis.
% =========================================================================

required_fields = { ...
    'stim_tag', ...
    'unit_ids', ...
    'raw_count', ...
    'condition_index_per_trial', ...
    'conditions'};

for k = 1:numel(required_fields)
    f = required_fields{k};
    if ~isfield(rf_data, f)
        error('Required field missing from RF bined data: %s', f);
    end
end

if ~(ismatrix(rf_data.raw_count) || ndims(rf_data.raw_count) == 3)
    error('rf_data.raw_count must be unit x trial or unit x trial x bin.');
end

if ~isnumeric(rf_data.condition_index_per_trial)
    error('rf_data.condition_index_per_trial must be numeric.');
end

if ~isstruct(rf_data.conditions)
    error('rf_data.conditions must be a struct array.');
end
end

function v = get_optional_unit_vector(rf_data, field_name, nUnit)
%% =========================================================================
% Read an optional unit-level vector from rf_data.
% If the field is missing, return NaN.
% =========================================================================

if isfield(rf_data, field_name)
    v = rf_data.(field_name);
    v = double(v(:));

    if numel(v) ~= nUnit
        warning('%s exists but its length does not match nUnit. Filling with NaN.', field_name);
        v = nan(nUnit, 1);
    end
else
    v = nan(nUnit, 1);
end
end

function [xPos_trial, yPos_trial, stimsize_trial] = get_trial_position_from_conditions(conditions, condition_index_per_trial, nTrial)
%% =========================================================================
% Reconstruct trial-wise xPos, yPos, and stimsize from condition definitions.
% =========================================================================

condition_index_per_trial = condition_index_per_trial(:);
if numel(condition_index_per_trial) ~= nTrial
    error('condition_index_per_trial length does not match number of trials.');
end

x_field = find_field_case_insensitive(conditions, ...
    {'xPos', 'xpos', 'x_pos', 'x', 'centerX', 'center_x'});
y_field = find_field_case_insensitive(conditions, ...
    {'yPos', 'ypos', 'y_pos', 'y', 'centerY', 'center_y'});
size_field = find_field_case_insensitive(conditions, ...
    {'stimsize', 'stim_size', 'stimSize', 'size'});

if isempty(x_field)
    error('No x position field found in rf_data.conditions.');
end
if isempty(y_field)
    error('No y position field found in rf_data.conditions.');
end
if isempty(size_field)
    error('No stimulus size field found in rf_data.conditions.');
end

xPos_trial = nan(nTrial, 1);
yPos_trial = nan(nTrial, 1);
stimsize_trial = nan(nTrial, 1);

nCond = numel(conditions);

for t = 1:nTrial
    c = condition_index_per_trial(t);

    if ~isfinite(c) || c ~= round(c) || c < 1 || c > nCond
        error('Invalid condition index at trial %d: %g', t, c);
    end

    xPos_trial(t) = get_numeric_scalar_from_struct(conditions(c), x_field);
    yPos_trial(t) = get_numeric_scalar_from_struct(conditions(c), y_field);
    stimsize_trial(t) = get_numeric_scalar_from_struct(conditions(c), size_field);
end
end

function field_name = find_field_case_insensitive(S, candidates)
%% =========================================================================
% Case-insensitive field-name lookup.
% =========================================================================

field_name = '';

if isempty(S) || ~isstruct(S)
    return;
end

fn = fieldnames(S);
fn_lower = lower(fn);

for k = 1:numel(candidates)
    cand = lower(candidates{k});
    idx = find(strcmp(fn_lower, cand), 1);
    if ~isempty(idx)
        field_name = fn{idx};
        return;
    end
end
end

function x = get_numeric_scalar_from_struct(S, field_name)
%% =========================================================================
% Extract one numeric scalar from a struct field.
% =========================================================================

val = S.(field_name);

if isnumeric(val) || islogical(val)
    if isempty(val) || numel(val) ~= 1
        error('Field %s is not a scalar.', field_name);
    end
    x = double(val);
    return;
end

if isstring(val)
    if isempty(val) || numel(val) ~= 1
        error('Field %s is not a scalar string.', field_name);
    end
    x = str2double(val);
    return;
end

if ischar(val)
    x = str2double(val);
    return;
end

error('Field %s has unsupported type.', field_name);
end

function rfstimsize = get_single_stimsize(stimsize_trial)
%% =========================================================================
% Require one unique finite stimulus size for the RF run.
% =========================================================================

finite_sizes = stimsize_trial(isfinite(stimsize_trial));

if isempty(finite_sizes)
    error('No finite stimsize values found for RF run.');
end

unique_sizes = unique(finite_sizes);

if numel(unique_sizes) ~= 1
    error(['RF run has multiple stimulus sizes. Current first-step fit ' ...
        'expects one scalar stimsize. Found %d unique sizes.'], ...
        numel(unique_sizes));
end

rfstimsize = unique_sizes(1);
end

function [RFmap, unique_x, unique_y] = compute_rf_map(response_count_trial, xPos_trial, yPos_trial)
%% =========================================================================
% Compute one-eye RF map:
%
% RFmap(unit, y, x) = mean response_count_trial(unit, trials at x/y)
% =========================================================================

[nUnit, nTrial] = size(response_count_trial);

xPos_trial = xPos_trial(:);
yPos_trial = yPos_trial(:);

if numel(xPos_trial) ~= nTrial || numel(yPos_trial) ~= nTrial
    error('xPos_trial/yPos_trial length does not match number of trials.');
end

good_trial = isfinite(xPos_trial) & isfinite(yPos_trial);
if ~any(good_trial)
    error('No finite x/y positions found.');
end

unique_x = unique(xPos_trial(good_trial));
unique_y = unique(yPos_trial(good_trial));

nXpos = numel(unique_x);
nYpos = numel(unique_y);

RFmap = nan(nUnit, nYpos, nXpos);

for x = 1:nXpos
    take_x = xPos_trial == unique_x(x);
    for y = 1:nYpos
        take_xy = take_x & yPos_trial == unique_y(y);
        if any(take_xy)
            RFmap(:, y, x) = mean(response_count_trial(:, take_xy), 2, 'omitnan');
        end
    end
end
end

function results = fitGaussianHeatmaps(data, xc, yc, stimsize)
%% =========================================================================
% Fit 2D Gaussians to RF heatmaps.
%
% Output:
% results.center  : N x 2 [x0, y0]
% results.size    : N x 2 [r95_x, r95_y]
% results.rsquare : N x 1 R^2
% results.params  : N x 6 [amp, x0, y0, sx, sy, offset]
%
% Notes:
% Apply a heuristic lower bound on corrected RF sigma.
% This lower bound is set to stimsize/4, so the minimum reported
% 95% RF width is approximately one stimulus size.
% =========================================================================

if nargin < 4
    error('fitGaussianHeatmaps requires data, xc, yc, and stimsize.');
end

[N, Y, X] = size(data);

xc = xc(:)';
yc = yc(:)';

if numel(xc) ~= X
    error('numel(xc) does not match size(data,3).');
end
if numel(yc) ~= Y
    error('numel(yc) does not match size(data,2).');
end

[XG, YG] = meshgrid(xc, yc);
xy_all = [XG(:), YG(:)];

results = struct();
results.center = nan(N, 2);
results.size = nan(N, 2);
results.rsquare = nan(N, 1);
results.params = nan(N, 6);

gauss2d = @(p,xy) p(1).*exp(-(((xy(:,1)-p(2)).^2)/(2*p(4)^2) + ...
    ((xy(:,2)-p(3)).^2)/(2*p(5)^2))) + p(6);

eps_sigma = 1e-3;
has_lsqcurvefit = exist('lsqcurvefit', 'file') == 2;

if has_lsqcurvefit
    opts_lsq = optimoptions('lsqcurvefit', ...
        'Display', 'off', ...
        'TolFun', 1e-6);
else
    warning('lsqcurvefit not found. Using fminsearch fallback for Gaussian fitting.');
    opts_fmin = optimset('Display', 'off', ...
        'TolFun', 1e-6, ...
        'MaxIter', 2000, ...
        'MaxFunEvals', 10000);
end

for k = 1:N
    Z = squeeze(data(k, :, :));
    zv_all = Z(:);

    finite_mask = isfinite(zv_all) & isfinite(xy_all(:,1)) & isfinite(xy_all(:,2));
    if ~any(finite_mask)
        continue;
    end

    xy = xy_all(finite_mask, :);
    zv = zv_all(finite_mask);

    total = sum(zv);
    if total <= 0
        continue;
    end

    amp0 = max(zv) - min(zv);
    off0 = min(zv);
    if amp0 <= 0
        amp0 = eps_sigma;
    end

    x0 = sum(xy(:,1) .* zv) / total;
    y0 = sum(xy(:,2) .* zv) / total;
    sx0 = sqrt(sum((xy(:,1) - x0).^2 .* zv) / total);
    sy0 = sqrt(sum((xy(:,2) - y0).^2 .* zv) / total);
    sx0 = max(sx0, eps_sigma);
    sy0 = max(sy0, eps_sigma);

    p0 = [amp0, x0, y0, sx0, sy0, off0];

    lb = [0, min(xc), min(yc), eps_sigma, eps_sigma, -Inf];
    ub = [Inf, max(xc), max(yc), Inf, Inf, Inf];

    if has_lsqcurvefit
        try
            [p, ~, ~, flag] = lsqcurvefit(gauss2d, p0, xy, zv, lb, ub, opts_lsq);
        catch
            p = p0;
            flag = -1;
        end

        if flag <= 0
            warning('Gaussian fit for unit index %d did not converge.', k);
        end
    else
        obj = @(p) gaussian_sse_with_penalty(p, xy, zv, gauss2d, lb, ub);
        try
            [p, ~, flag] = fminsearch(obj, p0, opts_fmin);
        catch
            p = p0;
            flag = -1;
        end
        p = enforce_gaussian_bounds(p, lb, ub, eps_sigma);

        if flag <= 0
            warning('Gaussian fminsearch fit for unit index %d did not converge.', k);
        end
    end

    zfit = gauss2d(p, xy);
    SSres = sum((zv - zfit).^2);
    SStot = sum((zv - mean(zv)).^2);

    if SStot > 0
        r2 = 1 - SSres / SStot;
    else
        r2 = NaN;
    end

    varK = stimsize^2 / 12;
    sxr = sqrt(max(p(4)^2 - varK, 0));
    syr = sqrt(max(p(5)^2 - varK, 0));

    sx = max(sxr, stimsize / 4);
    sy = max(syr, stimsize / 4);

    r95 = 2 * 1.96 * [sx, sy];

    results.center(k, :) = [p(2), p(3)];
    results.size(k, :) = r95;
    results.rsquare(k) = r2;
    results.params(k, :) = p;
end
end

function sse = gaussian_sse_with_penalty(p, xy, zv, gauss2d, lb, ub)
%% =========================================================================
% Objective for fminsearch fallback.
% =========================================================================

penalty = 0;

for j = 1:numel(p)
    if isfinite(lb(j)) && p(j) < lb(j)
        penalty = penalty + 1e12 * (lb(j) - p(j))^2;
    end
    if isfinite(ub(j)) && p(j) > ub(j)
        penalty = penalty + 1e12 * (p(j) - ub(j))^2;
    end
end

if p(4) <= 0 || p(5) <= 0 || p(1) < 0
    penalty = penalty + 1e12;
end

try
    zfit = gauss2d(p, xy);
    residual = zv - zfit;
    sse = sum(residual.^2) + penalty;
catch
    sse = Inf;
end
end

function p = enforce_gaussian_bounds(p, lb, ub, eps_sigma)
%% =========================================================================
% Enforce Gaussian parameter bounds after fminsearch fallback.
% =========================================================================

for j = 1:numel(p)
    if isfinite(lb(j))
        p(j) = max(p(j), lb(j));
    end
    if isfinite(ub(j))
        p(j) = min(p(j), ub(j));
    end
end

p(1) = max(p(1), 0);
p(4) = max(p(4), eps_sigma);
p(5) = max(p(5), eps_sigma);
end

function plot_rf_map_depth_desc(unit_rf, png_file)
%% =========================================================================
% Plot RF maps sorted by descending unit_depth_um.
%
% This version builds one RGB mosaic image directly and saves it with imwrite.
% =========================================================================

RFmap = unit_rf.rfs.map;
unit_depth_um = unit_rf.unit_depth_um(:);
nUnit = numel(unit_rf.unit_ids);

if size(RFmap, 1) ~= nUnit
    error('size(unit_rf.rfs.map,1) does not match numel(unit_rf.unit_ids).');
end

order = get_depth_desc_order(unit_depth_um, nUnit);

tile_scale = 5;
gap_pix = 1;
margin_pix = 6;

nCols = ceil(sqrt(nUnit * 1.25));
nRows = ceil(nUnit / nCols);

[~, nY, nX] = size(RFmap);
tile_h = nY * tile_scale;
tile_w = nX * tile_scale;

img_h = 2 * margin_pix + nRows * tile_h + (nRows - 1) * gap_pix;
img_w = 2 * margin_pix + nCols * tile_w + (nCols - 1) * gap_pix;

mosaic_rgb = ones(img_h, img_w, 3);
cmap = parula(256);

for ii = 1:nUnit
    u = order(ii);
    row = floor((ii - 1) / nCols) + 1;
    col = mod(ii - 1, nCols) + 1;

    y0 = margin_pix + (row - 1) * (tile_h + gap_pix) + 1;
    x0 = margin_pix + (col - 1) * (tile_w + gap_pix) + 1;

    yidx = y0:(y0 + tile_h - 1);
    xidx = x0:(x0 + tile_w - 1);

    Z = squeeze(RFmap(u, :, :));
    Z = flipud(Z);

    tile_rgb = rf_tile_to_rgb(Z, cmap, tile_scale);
    mosaic_rgb(yidx, xidx, :) = tile_rgb;
end

imwrite(mosaic_rgb, png_file);
end

function tile_rgb = rf_tile_to_rgb(Z, cmap, tile_scale)
%% =========================================================================
% Convert one RF map to an RGB tile.
% Each unit is scaled independently, similar to using imagesc separately.
% =========================================================================

finite_mask = isfinite(Z);
tile_h = size(Z, 1) * tile_scale;
tile_w = size(Z, 2) * tile_scale;
tile_rgb = ones(tile_h, tile_w, 3);

if ~any(finite_mask(:))
    return;
end

zmin = min(Z(finite_mask));
zmax = max(Z(finite_mask));

if zmax > zmin
    Zn = (Z - zmin) ./ (zmax - zmin);
else
    if zmax == 0
        Zn = zeros(size(Z));
    else
        Zn = 0.5 * ones(size(Z));
    end
end

Zn(~finite_mask) = NaN;

idx = round(Zn * 255) + 1;
idx(idx < 1) = 1;
idx(idx > 256) = 256;

rgb_small = ones(size(Z,1), size(Z,2), 3);
for yy = 1:size(Z,1)
    for xx = 1:size(Z,2)
        if isfinite(idx(yy,xx))
            rgb_small(yy,xx,:) = reshape(cmap(idx(yy,xx), :), [1 1 3]);
        else
            rgb_small(yy,xx,:) = [1 1 1];
        end
    end
end

tile_rgb = repelem(rgb_small, tile_scale, tile_scale, 1);
end

function order = get_depth_desc_order(unit_depth_um, nUnit)
%% =========================================================================
% Sort by descending depth.
% Same depth uses original unit index order.
% NaN depth is placed last.
% =========================================================================

unit_depth_um = unit_depth_um(:);
if numel(unit_depth_um) ~= nUnit
    unit_depth_um = nan(nUnit, 1);
end

unit_index = (1:nUnit)';
finite_depth = isfinite(unit_depth_um);

finite_idx = unit_index(finite_depth);
finite_depth_val = unit_depth_um(finite_depth);

if isempty(finite_idx)
    order = unit_index;
    return;
end

sort_table = [-finite_depth_val, finite_idx];
[~, sort_pos] = sortrows(sort_table, [1 2]);
finite_order = finite_idx(sort_pos);
nan_order = unit_index(~finite_depth);

order = [finite_order; nan_order];
end

%% ======================= Step 2 functions =======================

function target_result = compute_target_rf_result(unit_rf, target_run, probe_id, ksDir, target_stimulus_xy_size, RF_R2_threshold)
%% =========================================================================
% Compute per-probe target overlay result.
%
% Keeps units that:
% 1. are used by target model run for this probe
% 2. have RF fit R2 >= threshold
% 3. have finite RF center and size
% =========================================================================

selected_unit_ids = get_selected_unit_ids_for_probe(target_run, probe_id);

r2 = unit_rf.fit.rsquare(:);
center = unit_rf.fit.center;
rf_size = unit_rf.fit.size;
unit_ids = unit_rf.unit_ids(:);

in_target = ismember(unit_ids, selected_unit_ids(:));
finite_fit = isfinite(r2) ...
    & all(isfinite(center), 2) ...
    & all(isfinite(rf_size), 2);

good = in_target & finite_fit & (r2 >= RF_R2_threshold);

good_idx = find(good);
good_center = center(good_idx, :);
good_size = rf_size(good_idx, :);
good_rsquared = r2(good_idx);
good_unit_ids = unit_ids(good_idx);

if isempty(good_idx)
    mean_center = [NaN NaN];
    mean_size = [NaN NaN];
else
    mean_center = mean(good_center, 1, 'omitnan');
    mean_size = mean(good_size, 1, 'omitnan');
end

target_result = struct();
target_result.probe_id = probe_id;
target_result.ksDir = ksDir;
target_result.selected_unit_ids = selected_unit_ids(:);
target_result.RF_R2_threshold = RF_R2_threshold;
target_result.good_unit_idx = good_idx;
target_result.good_unit_ids = good_unit_ids;
target_result.good_center = good_center;
target_result.good_size = good_size;
target_result.good_rsquared = good_rsquared;
target_result.mean_center = mean_center;
target_result.mean_size = mean_size;
target_result.target_stimulus_xy_size = target_stimulus_xy_size;
target_result.target_overlay_png = '';
target_result.target_overlay_readable_png = '';
target_result.target_overlay_full_png = '';
end

function selected_unit_ids = get_selected_unit_ids_for_probe(target_run, probe_id)
%% =========================================================================
% Dynamically extract selected unit ids for the current probe.
%
% Expected field examples:
% probe0_usedunit_ids
% probe1_usedunit_ids
% =========================================================================

field_name = sprintf('probe%d_usedunit_ids', probe_id);

if ~isfield(target_run, field_name)
    error('Target run does not contain field: %s', field_name);
end

selected_unit_ids = target_run.(field_name);
selected_unit_ids = selected_unit_ids(:);
end

function target_stimulus_xy_size = extract_target_stimulus_xy_size(target_run, target_default_stimsize)
%% =========================================================================
% Extract unique target stimulus [x, y, size] combinations from target run.
%
% Field lookup:
% x    : centerX, xPos, x
% y    : centerY, yPos, y
% size : size, stimsize
%
% If size field is missing, target_default_stimsize is used if non-empty.
% =========================================================================

conditions = get_target_conditions(target_run);

x_field = find_field_case_insensitive(conditions, ...
    {'centerX', 'center_x', 'xPos', 'xpos', 'x_pos', 'x'});
y_field = find_field_case_insensitive(conditions, ...
    {'centerY', 'center_y', 'yPos', 'ypos', 'y_pos', 'y'});
size_field = find_field_case_insensitive(conditions, ...
    {'size', 'stimsize', 'stim_size', 'stimSize'});

if isempty(x_field)
    error('No x-position field found in target run conditions.');
end
if isempty(y_field)
    error('No y-position field found in target run conditions.');
end
if isempty(size_field) && isempty(target_default_stimsize)
    error(['No size/stimsize field found in target run conditions, and ' ...
        'target_default_stimsize is empty.']);
end

nCond = numel(conditions);
xyzs = nan(nCond, 3);

for c = 1:nCond
    xyzs(c, 1) = get_numeric_scalar_from_struct(conditions(c), x_field);
    xyzs(c, 2) = get_numeric_scalar_from_struct(conditions(c), y_field);

    if isempty(size_field)
        xyzs(c, 3) = target_default_stimsize;
    else
        xyzs(c, 3) = get_numeric_scalar_from_struct(conditions(c), size_field);
    end
end

good = all(isfinite(xyzs), 2);
if ~any(good)
    error('No finite target stimulus [x, y, size] combinations found.');
end

target_stimulus_xy_size = unique(xyzs(good, :), 'rows');
end

function conditions = get_target_conditions(target_run)
%% =========================================================================
% Get condition structure from model_data_allruns target run.
% Prefer full conditions when available.
% =========================================================================

if isfield(target_run, 'conditions_full')
    conditions = target_run.conditions_full;
elseif isfield(target_run, 'conditions')
    conditions = target_run.conditions;
else
    error('Target run does not contain conditions_full or conditions.');
end

if ~isstruct(conditions)
    error('Target run conditions must be a struct array.');
end
end

function plot_target_rf_overlay_per_probe(target_stimulus_xy_size, centers, sizes, probe_color, png_file, plot_title, view_mode)
%% =========================================================================
% Plot target stimuli and individual unit RF fit errorbar-like crosses.
%
% view_mode:
%   'readable' : axis limits use target layout + RF centers only.
%                Huge RF size does not expand the window.
%   'full'     : axis limits include full RF extents plus zero meridians.
%                This is useful for QC while keeping HM/VM visible.
% =========================================================================

if nargin < 7 || isempty(view_mode)
    view_mode = 'full';
end

fig = figure('Visible', 'off', ...
    'Color', 'w', ...
    'Position', [100 100 850 800]);
ax = axes('Parent', fig);
hold(ax, 'on');

draw_target_stimulus_circles(ax, target_stimulus_xy_size);
draw_rf_errorbars(ax, centers, sizes, probe_color, 1.1);

axis(ax, 'equal');
box(ax, 'off');
grid(ax, 'off');
xlabel(ax, '');
ylabel(ax, '');
title(ax, plot_title, 'Interpreter', 'none');

set_overlay_axis_limits(ax, target_stimulus_xy_size, centers, sizes, view_mode);
format_rf_overlay_axis(ax);

save_figure_png(fig, png_file, 200);
close(fig);
end

function plot_target_rf_combined(target_rf_summary, probe_colors, png_file, view_mode)
%% =========================================================================
% Plot target stimuli and probe-average RF ellipses.
%
% view_mode:
%   'readable' : axis limits include target layout + full mean RF ellipse extents,
%                so probe-average ellipses are not clipped.
%   'full'     : axis limits include mean RF ellipse extents plus zero meridians.
% =========================================================================

if nargin < 4 || isempty(view_mode)
    view_mode = 'full';
end

fig = figure('Visible', 'off', ...
    'Color', 'w', ...
    'Position', [100 100 850 800]);
ax = axes('Parent', fig);
hold(ax, 'on');

target_stimulus_xy_size = target_rf_summary.target_stimulus_xy_size;
draw_target_stimulus_circles(ax, target_stimulus_xy_size);

all_centers = [];
all_sizes = [];

for i = 1:numel(target_rf_summary.probe)
    probe_id = target_rf_summary.probe(i).probe_id;
    probe_color = get_probe_color(probe_id, probe_colors);

    mean_center = target_rf_summary.probe(i).mean_center;
    mean_size = target_rf_summary.probe(i).mean_size;

    if all(isfinite(mean_center)) && all(isfinite(mean_size))
        draw_rf_ellipse(ax, mean_center, mean_size, probe_color, 2.2);
        plot(ax, mean_center(1), mean_center(2), 'o', ...
            'MarkerFaceColor', probe_color, ...
            'MarkerEdgeColor', probe_color, ...
            'MarkerSize', 5);

        all_centers = [all_centers; mean_center]; %#ok<AGROW>
        all_sizes = [all_sizes; mean_size]; %#ok<AGROW>
    end
end

axis(ax, 'equal');
box(ax, 'off');
grid(ax, 'off');
xlabel(ax, '');
ylabel(ax, '');
title(ax, sprintf('Probe-average RF over target: %s, %s', ...
    target_rf_summary.target_stim_tag, view_mode), ...
    'Interpreter', 'none');

set_overlay_axis_limits(ax, target_stimulus_xy_size, all_centers, all_sizes, view_mode, true);
format_rf_overlay_axis(ax);

save_figure_png(fig, png_file, 200);
close(fig);
end

function draw_target_stimulus_circles(ax, target_stimulus_xy_size)
%% =========================================================================
% Draw target stimuli as filled circles.
% Smaller stimulus is darker; larger stimulus is lighter.
% Larger circles are drawn first so smaller circles remain visible.
% =========================================================================

if isempty(target_stimulus_xy_size)
    return;
end

x = target_stimulus_xy_size(:,1);
y = target_stimulus_xy_size(:,2);
s = target_stimulus_xy_size(:,3);

[~, order] = sort(s, 'descend');
s_min = min(s);
s_max = max(s);

theta = linspace(0, 2*pi, 80);

for ii = 1:numel(order)
    k = order(ii);
    cx = x(k);
    cy = y(k);
    r = s(k) / 2;

    if s_max > s_min
        gray = 0.15 + 0.60 * (s(k) - s_min) / (s_max - s_min);
    else
        gray = 0.45;
    end

    xx = cx + r * cos(theta);
    yy = cy + r * sin(theta);

    patch(ax, xx, yy, gray * [1 1 1], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.75);
end
end

function draw_rf_errorbars(ax, centers, sizes, color_value, line_width)
%% =========================================================================
% Draw RF fit centers and x/y full-width sizes as errorbar-like crosses.
% size(:,1) is full x width, size(:,2) is full y width.
% =========================================================================

if isempty(centers)
    return;
end

for i = 1:size(centers, 1)
    c = centers(i, :);
    sz = sizes(i, :);

    if ~all(isfinite(c)) || ~all(isfinite(sz))
        continue;
    end

    hx = sz(1) / 2;
    hy = sz(2) / 2;

    line(ax, [c(1) - hx, c(1) + hx], [c(2), c(2)], ...
        'Color', color_value, ...
        'LineWidth', line_width);
    line(ax, [c(1), c(1)], [c(2) - hy, c(2) + hy], ...
        'Color', color_value, ...
        'LineWidth', line_width);
end
end

function draw_rf_ellipse(ax, center, rf_size, color_value, line_width)
%% =========================================================================
% Draw one RF ellipse using full x/y width.
% =========================================================================

theta = linspace(0, 2*pi, 160);
rx = rf_size(1) / 2;
ry = rf_size(2) / 2;

x = center(1) + rx * cos(theta);
y = center(2) + ry * sin(theta);

line(ax, x, y, ...
    'Color', color_value, ...
    'LineWidth', line_width);
end

function set_overlay_axis_limits(ax, target_stimulus_xy_size, centers, sizes, view_mode, include_size_in_readable)
%% =========================================================================
% Set axis limits for RF-target overlay.
%
% view_mode = 'full':
%   Include target circles and full RF extents, plus zero meridians.
%   This is the QC-style view while still guaranteeing HM/VM visibility.
%
% view_mode = 'readable':
%   By default, include target circles and RF centers, but do not let huge RF
%   sizes determine the display window. The RF crosses are still drawn and are
%   simply clipped by the axes. For combined probe-average ellipse plots, set
%   include_size_in_readable = true so readable views include the full mean
%   ellipses and do not cut them.
% =========================================================================

if nargin < 5 || isempty(view_mode)
    view_mode = 'full';
end
if nargin < 6 || isempty(include_size_in_readable)
    include_size_in_readable = false;
end

xs = [];
ys = [];

% Always include target circles.
if ~isempty(target_stimulus_xy_size)
    x = target_stimulus_xy_size(:,1);
    y = target_stimulus_xy_size(:,2);
    s = target_stimulus_xy_size(:,3);

    xs = [xs; x - s/2; x + s/2]; %#ok<AGROW>
    ys = [ys; y - s/2; y + s/2]; %#ok<AGROW>
end

% Include RF centers in both modes. Include full RF extents in full mode,
% and optionally in readable mode for combined mean-ellipse plots.
if ~isempty(centers)
    good_center = all(isfinite(centers), 2);

    if strcmpi(view_mode, 'full') || (strcmpi(view_mode, 'readable') && include_size_in_readable)
        good = good_center & all(isfinite(sizes), 2);
        if any(good)
            c = centers(good, :);
            sz = sizes(good, :);
            xs = [xs; c(:,1) - sz(:,1)/2; c(:,1) + sz(:,1)/2]; %#ok<AGROW>
            ys = [ys; c(:,2) - sz(:,2)/2; c(:,2) + sz(:,2)/2]; %#ok<AGROW>
        end
    else
        if any(good_center)
            c = centers(good_center, :);
            xs = [xs; c(:,1)]; %#ok<AGROW>
            ys = [ys; c(:,2)]; %#ok<AGROW>
        end
    end
end

% Always keep HM/VM available as visual references.
% This guarantees that y = 0 horizontal meridian and x = 0 vertical
% meridian are inside the plotting window in both readable and full modes.
xs = [xs; 0]; %#ok<AGROW>
ys = [ys; 0]; %#ok<AGROW>

xs = xs(isfinite(xs));
ys = ys(isfinite(ys));

if isempty(xs) || isempty(ys)
    return;
end

xmin = min(xs);
xmax = max(xs);
ymin = min(ys);
ymax = max(ys);

xrange = xmax - xmin;
yrange = ymax - ymin;

if xrange == 0
    xrange = 1;
end
if yrange == 0
    yrange = 1;
end

if strcmpi(view_mode, 'readable')
    pad = 0.25 * max(xrange, yrange);
    min_half_width = 8; % deg
else
    pad = 0.08 * max(xrange, yrange);
    min_half_width = 0;
end

xlim_candidate = [xmin - pad, xmax + pad];
ylim_candidate = [ymin - pad, ymax + pad];

if min_half_width > 0
    cx = mean([xmin xmax]);
    cy = mean([ymin ymax]);

    if diff(xlim_candidate) < 2 * min_half_width
        xlim_candidate = [cx - min_half_width, cx + min_half_width];
    end
    if diff(ylim_candidate) < 2 * min_half_width
        ylim_candidate = [cy - min_half_width, cy + min_half_width];
    end
end

xlim(ax, xlim_candidate);
ylim(ax, ylim_candidate);
end

function save_figure_png(fig, png_file, dpi)
%% =========================================================================
% Save figure to PNG.
% =========================================================================

drawnow;

try
    if exist('exportgraphics', 'file') == 2
        exportgraphics(fig, png_file, ...
            'Resolution', dpi, ...
            'BackgroundColor', 'white');
    else
        set(fig, 'PaperPositionMode', 'auto');
        print(fig, png_file, '-dpng', sprintf('-r%d', dpi));
    end
catch
    saveas(fig, png_file);
end
end

function color_value = get_probe_color(probe_id, probe_colors)
%% =========================================================================
% Get plotting color for probe id.
% probe_id 0 uses row 1, probe_id 1 uses row 2, etc.
% =========================================================================

row = probe_id + 1;

if row >= 1 && row <= size(probe_colors, 1)
    color_value = probe_colors(row, :);
else
    color_value = [0 0 0];
end
end

function safe_name = make_filename_safe(s)
%% =========================================================================
% Convert stim tag to safe filename component.
% =========================================================================

s = char(s);
safe_name = regexprep(s, '[^A-Za-z0-9_-]', '_');
safe_name = regexprep(safe_name, '_+', '_');
safe_name = regexprep(safe_name, '^_+', '');
safe_name = regexprep(safe_name, '_+$', '');

if isempty(safe_name)
    safe_name = 'target_run';
end

max_len = 80;
if numel(safe_name) > max_len
    safe_name = safe_name(1:max_len);
end
end

function all_tags = get_all_run_tags(model_data_allruns)
%% =========================================================================
% Extract all stim_tag values from model_data_allruns.
% =========================================================================

all_tags = cell(numel(model_data_allruns), 1);

for j = 1:numel(model_data_allruns)
    if ~isfield(model_data_allruns{j}, 'stim_tag')
        error('stim_tag missing in model_data_allruns{%d}.', j);
    end
    all_tags{j} = model_data_allruns{j}.stim_tag;
end
end

function format_rf_overlay_axis(ax)
%% =========================================================================
% Format overlay plot like a visual-field/RF plot.
%
% This version:
% 1. Removes the outer box and default edge axes.
% 2. Uses y=0 as HM, horizontal meridian.
% 3. Uses x=0 as VM, vertical meridian.
% 4. Adds ruler-like tick marks and numeric labels.
% 5. Adds a 2 deg scale bar.
% =========================================================================

hold(ax, 'on');

xl = xlim(ax);
yl = ylim(ax);
xrange = xl(2) - xl(1);
yrange = yl(2) - yl(1);

box(ax, 'off');
grid(ax, 'off');
set(ax, ...
    'XColor', 'none', ...
    'YColor', 'none', ...
    'FontSize', 11, ...
    'LineWidth', 1);

xlabel(ax, '');
ylabel(ax, '');

axis_color = [0 0 0];
axis_lw = 2.4;       % main axis width
tick_lw = 1.4;       % ruler tick width
tick_len_frac = 0.018;
label_font = 9;

has_HM = yl(1) <= 0 && yl(2) >= 0; % y = 0
has_VM = xl(1) <= 0 && xl(2) >= 0; % x = 0

% ---------------------------------------------------------------------
% Draw HM and VM.
% ---------------------------------------------------------------------
if has_HM
    line(ax, xl, [0 0], ...
        'Color', axis_color, ...
        'LineStyle', '-', ...
        'LineWidth', axis_lw);
end

if has_VM
    line(ax, [0 0], yl, ...
        'Color', axis_color, ...
        'LineStyle', '-', ...
        'LineWidth', axis_lw);
end

xlim(ax, xl);
ylim(ax, yl);

% ---------------------------------------------------------------------
% Ruler ticks.
% Use an adaptive tick step so large RF ranges do not create overcrowded labels.
% ---------------------------------------------------------------------
target_n_tick = 8;
tick_step = choose_nice_tick_step(max(xrange, yrange) / target_n_tick);

if has_HM
    xticks_internal = ceil(xl(1) / tick_step) * tick_step : tick_step : floor(xl(2) / tick_step) * tick_step;
    xticks_internal = xticks_internal(abs(xticks_internal) > 1e-9);

    tick_len_y = tick_len_frac * yrange;

    for i = 1:numel(xticks_internal)
        x = xticks_internal(i);
        line(ax, [x x], [-tick_len_y tick_len_y], ...
            'Color', axis_color, ...
            'LineWidth', tick_lw);
        text(ax, x, -2.3 * tick_len_y, ...
            sprintf('%g', x), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', ...
            'FontSize', label_font, ...
            'Color', axis_color);
    end
end

if has_VM
    yticks_internal = ceil(yl(1) / tick_step) * tick_step : tick_step : floor(yl(2) / tick_step) * tick_step;
    yticks_internal = yticks_internal(abs(yticks_internal) > 1e-9);

    tick_len_x = tick_len_frac * xrange;

    for i = 1:numel(yticks_internal)
        y = yticks_internal(i);
        line(ax, [-tick_len_x tick_len_x], [y y], ...
            'Color', axis_color, ...
            'LineWidth', tick_lw);
        text(ax, -2.0 * tick_len_x, y, ...
            sprintf('%g', y), ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', label_font, ...
            'Color', axis_color);
    end
end

% ---------------------------------------------------------------------
% Meridian labels.
% ---------------------------------------------------------------------
if has_HM
    text(ax, xl(2) - 0.02 * xrange, 0 + 0.04 * yrange, ...
        'HM', ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 13, ...
        'FontWeight', 'bold', ...
        'Color', axis_color);
end

if has_VM
    text(ax, 0 + 0.035 * xrange, yl(2) - 0.02 * yrange, ...
        'VM', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 13, ...
        'FontWeight', 'bold', ...
        'Color', axis_color);
end

if has_HM && has_VM
    text(ax, 0 + 0.02 * xrange, 0 - 0.035 * yrange, ...
        '0', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'FontSize', label_font, ...
        'Color', axis_color);
end

draw_scale_bar_deg(ax, 2);
end

function step = choose_nice_tick_step(raw_step)
%% =========================================================================
% Choose a readable tick step: 1, 2, 5, 10, 20, 50, ...
% =========================================================================

if ~isfinite(raw_step) || raw_step <= 0
    step = 2;
    return;
end

base = 10 ^ floor(log10(raw_step));
candidates = base * [1 2 5 10];
idx = find(candidates >= raw_step, 1, 'first');

if isempty(idx)
    step = candidates(end);
else
    step = candidates(idx);
end
end

function draw_scale_bar_deg(ax, scale_len)
%% =========================================================================
% Draw a horizontal scale bar in degree units.
%
% This version places the scale bar farther left to avoid overlap with
% VM tick marks.
% =========================================================================

xl = xlim(ax);
yl = ylim(ax);
xrange = xl(2) - xl(1);
yrange = yl(2) - yl(1);

if xrange <= 0 || yrange <= 0
    return;
end

% Put scale bar lower-left.
y0 = yl(1) + 0.075 * yrange;

% Prefer placing it clearly on the left side of VM if x=0 is visible.
if xl(1) < 0 && xl(2) > 0
    x0 = xl(1) + 0.025 * xrange;

    % Make sure the right end is not too close to VM.
    right_margin_from_VM = 0.08 * xrange;
    if x0 + scale_len > -right_margin_from_VM
        x0 = -right_margin_from_VM - scale_len;
    end

    % If this pushes it outside the plot, fall back to near left edge.
    if x0 < xl(1) + 0.01 * xrange
        x0 = xl(1) + 0.025 * xrange;
    end
else
    x0 = xl(1) + 0.025 * xrange;
end

% If the scale bar exceeds the right edge, pull it back.
if x0 + scale_len > xl(2) - 0.04 * xrange
    x0 = xl(2) - scale_len - 0.08 * xrange;
end

if x0 < xl(1) || x0 + scale_len > xl(2)
    return;
end

line(ax, [x0, x0 + scale_len], [y0, y0], ...
    'Color', 'k', ...
    'LineWidth', 2.2);

text(ax, x0 + scale_len / 2, y0 - 0.035 * yrange, ...
    sprintf('%g deg', scale_len), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 11, ...
    'Color', 'k');
end
