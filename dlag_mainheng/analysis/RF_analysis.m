%% =========================================================================
% RF_analysis
%
% Purpose:
%   First-step RF analysis based on bined_data_allruns.mat.
%
%   This script:
%     1. Finds the RF run in bined_data_allruns by rf_stim_tag.
%     2. Sums raw_count across bins to get one response per unit per trial.
%     3. Reconstructs trial-wise xPos, yPos, and stimsize from conditions.
%     4. Computes RF maps: unit x y-position x x-position.
%     5. Fits a 2D Gaussian to each unit RF map.
%     6. Plots RF maps sorted by descending depth.
%     7. Saves unit_rf_results.mat and unit_rf_map_depth_desc.png.
%
% Required input in each kilosort folder:
%   - bined_data_allruns.mat
%
% Expected RF run fields in bined_data_allruns{r}:
%   - stim_tag
%   - unit_ids
%   - raw_count
%   - condition_index_per_trial
%   - conditions
%   - optionally unit_depth_um
%   - optionally unit_channel
%
% Main output saved in each kilosort folder:
%   - unit_rf_results.mat
%   - unit_rf_map_depth_desc.png
%
% Output variable:
%   unit_rf.rf_stim_tag
%   unit_rf.unit_ids
%   unit_rf.unit_depth_um
%   unit_rf.unit_channel
%   unit_rf.response_count_trial
%   unit_rf.xPos_trial
%   unit_rf.yPos_trial
%   unit_rf.stimsize_trial
%   unit_rf.rfs.map
%   unit_rf.rfs.x
%   unit_rf.rfs.y
%   unit_rf.rfstimsize
%   unit_rf.fit.center
%   unit_rf.fit.size
%   unit_rf.fit.rsquare
%   unit_rf.fit.params
%
% Notes:
%   1. raw_count is assumed to be unit x trial x bin.
%   2. RF response is defined as sum(raw_count, 3).
%   3. RF map collapses across any non-position condition dimensions.
%      For example, if ori is also present, trials with the same xPos/yPos
%      are averaged together regardless of ori.
%   4. Depth sorting for plotting:
%        larger depth first, then original unit index.
%      This matches the convention that depth 0 is at the probe tip.
% =========================================================================

clc;
clear;

addpath(genpath(fullfile('.', 'expo_tools')));
addpath(genpath(fullfile('.', 'utils')));

%% ----------------------- User parameters -----------------------

root_folder = 'I:\np_data';
runName     = 'RafiL001p0120';
runind      = 1;          % run index after -g
probes      = [0];      % probe indices after -prb

% RF run used for computing unit RF.
rf_stim_tag = '[RFG_coarse2dg_99_4_150isi]';

% Target run reserved for later step.
% This script does not use target_stim_tag yet.
target_stim_tag = '_2[Gpl2_2c_2sz_400_2_200isi]'; %#ok<NASGU>

% Output names.
rf_mat_name = 'unit_rf_results.mat';
rf_png_name = 'unit_rf_map_depth_desc.png';

%% ----------------------- Build shared session paths -----------------------

run_g   = sprintf('%s_g%d', runName, runind);
destDir = fullfile(root_folder, run_g);

fprintf('destDir : %s\n', destDir);
fprintf('RF stim tag: %s\n', rf_stim_tag);

%% ----------------------- Process each probe folder -----------------------

for ip = 1:numel(probes)

    thisProbe = probes(ip);
    imecStr   = sprintf('imec%d', thisProbe);

    probe_folder = fullfile(destDir, ['catgt_' run_g], [run_g '_' imecStr]);

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

            unit_rf = compute_unit_rf_from_bined_run(rf_data, rf_stim_tag);

            mat_file = fullfile(ksDir, rf_mat_name);
            png_file = fullfile(ksDir, rf_png_name);

            save(mat_file, 'unit_rf');

            plot_rf_map_depth_desc(unit_rf, png_file);

            fprintf('Saved:\n');
            fprintf('  %s\n', mat_file);
            fprintf('  %s\n', png_file);

        catch ME
            fprintf(2, 'Error in probe %d, ksDir %s\n', thisProbe, ksDir);
            fprintf(2, '%s\n', ME.message);
        end
    end
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
%   raw_count = unit x trial x bin
%   raw_count = unit x trial        when there is only one bin
% =========================================================================

    validate_rf_bined_data(rf_data);

    unit_ids = rf_data.unit_ids(:);
    raw_count = double(rf_data.raw_count);

    nUnit = size(raw_count, 1);
    nTrial = size(raw_count, 2);

    if numel(unit_ids) ~= nUnit
        error('numel(unit_ids) does not match size(raw_count,1).');
    end

    % ---------------------------------------------------------------------
    % Convert raw_count to one response per unit per trial.
    %
    % If raw_count is unit x trial x bin, sum over bins.
    % If raw_count is unit x trial, it already represents the single-bin
    % response count per trial.
    % ---------------------------------------------------------------------
    if ismatrix(raw_count)
        response_count_trial = raw_count;
    else
        response_count_trial = sum(raw_count, 3);
        response_count_trial = reshape(response_count_trial, [nUnit, nTrial]);
    end

    unit_depth_um = get_optional_unit_vector(rf_data, 'unit_depth_um', nUnit);
    unit_channel  = get_optional_unit_vector(rf_data, 'unit_channel',  nUnit);

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
%
% This function intentionally uses only x position, y position, and stimulus
% size. Other condition fields, for example ori, are ignored for RF map
% grouping later.
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
%   RFmap(unit, y, x) = mean response_count_trial(unit, trials at x/y)
%
% This follows the old get_RFmap1eye logic, but uses unit responses from
% bined_data_allruns.
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
                RFmap(:, y, x) = mean_omitnan_dim2(response_count_trial(:, take_xy));
            end
        end
    end
end

function m = mean_omitnan_dim2(A)
%% =========================================================================
% Mean over columns, ignoring NaN, without requiring nanmean.
% =========================================================================

    nRow = size(A, 1);
    m = nan(nRow, 1);

    for i = 1:nRow
        v = A(i, :);
        v = v(isfinite(v));

        if ~isempty(v)
            m(i) = mean(v);
        end
    end
end

function results = fitGaussianHeatmaps(data, xc, yc, stimsize)
%% =========================================================================
% Fit 2D Gaussians to RF heatmaps.
%
% Inputs:
%   data     : N x Y x X RF maps
%   xc       : X x 1 or 1 x X vector of x coordinates
%   yc       : Y x 1 or 1 x Y vector of y coordinates
%   stimsize : scalar stimulus size
%
% Output:
%   results.center  : N x 2 [x0, y0]
%   results.size    : N x 2 [r95_x, r95_y]
%   results.rsquare : N x 1 R^2
%   results.params  : N x 6 [amp, x0, y0, sx, sy, offset]
%
% Notes:
%   This is the no-plot version of the old fitGaussianHeatmaps logic.
%   If lsqcurvefit is available, bounded least-squares is used.
%   If lsqcurvefit is unavailable, fminsearch is used as a fallback.
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

gauss2d = @(p,xy) p(1).*exp(-(((xy(:,1)-p(2)).^2)/(2*p(4)^2) + ((xy(:,2)-p(3)).^2)/(2*p(5)^2))) + p(6);
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

        lb = [0,   min(xc), min(yc), eps_sigma, eps_sigma, -Inf];
        ub = [Inf, max(xc), max(yc), Inf,       Inf,        Inf];

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

        sx = max(sxr, stimsize / 2);
        sy = max(syr, stimsize / 2);

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
% Same depth is sorted by original unit index.
% Missing depth is placed at the end in original unit order.
% =========================================================================

    RFmap = unit_rf.rfs.map;
    unit_ids = unit_rf.unit_ids(:);
    unit_depth_um = unit_rf.unit_depth_um(:);
    unit_channel = unit_rf.unit_channel(:);

    nUnit = numel(unit_ids);

    if size(RFmap, 1) ~= nUnit
        error('size(unit_rf.rfs.map,1) does not match numel(unit_rf.unit_ids).');
    end

    order = get_depth_desc_order(unit_depth_um, nUnit);

    nCols = ceil(sqrt(nUnit));
    nRows = ceil(nUnit / nCols);

    fig = figure('Visible', 'off', ...
                 'Color', 'w', ...
                 'Position', [100 100 1800 1400]);

    colormap(parula);

    for ii = 1:nUnit

        u = order(ii);

        subplot(nRows, nCols, ii);

        Z = squeeze(RFmap(u, :, :));

        imagesc(flipud(Z));

        axis equal tight;
        axis off;

        title_str = sprintf('idx=%d id=%g d=%.1f ch=%.0f', ...
            u, unit_ids(u), unit_depth_um(u), unit_channel(u));

        title(title_str, ...
            'Interpreter', 'none', ...
            'FontSize', 5);
    end

    if exist('sgtitle', 'file') == 2
        sgtitle(sprintf('RF maps sorted by descending depth: %s', unit_rf.rf_stim_tag), ...
            'Interpreter', 'none');
    end

    set(fig, 'PaperPositionMode', 'auto');

    try
        print(fig, png_file, '-dpng', '-r200');
    catch
        saveas(fig, png_file);
    end

    close(fig);
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