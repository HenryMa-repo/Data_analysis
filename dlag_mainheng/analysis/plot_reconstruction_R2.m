%% plot_reconstruction_R2.m
% Plot neuron-wise reconstruction R2 comparisons.
%
% Run this script after data_reconstruction.m has saved reconstruction_R2.mat.
%
% For each condition, or for pooled all-condition mode, this script plots
% one 2 x 4 figure:
%
% Row 1:
%   1) y = use across          , x = use all
%   2) y = across excl within  , x = use all
%   3) y = across excl within  , x = use across
%   4) across delta R2 violin + individual neuron points:
%        all-use, all-excl, use-excl
%
% Row 2:
%   5) y = use within          , x = use all
%   6) y = within excl across  , x = use all
%   7) y = within excl across  , x = use within
%   8) within delta R2 violin + individual neuron points:
%        all-use, all-excl, use-excl
%
% Each scatter point is one neuron.
% Group 1 and Group 2 are plotted using different colors.
%
% Output files, if enabled:
%   reconstruction_R2_comparison.fig
%   reconstruction_R2_comparison.png
%
% Since each condition has its own output folder, output file names do not
% include condition IDs or condition labels.

clc;
clear;

%% ------------------------------------------------------------------------
% User parameters
% -------------------------------------------------------------------------

data_content = 'raw_count';
% options usually include:
% raw_count, raw_fr, z_within_trial, z_within_condition,
% z_across_conditions, demean_count_within_trial, demean_fr_within_trial,
% demean_pooledsd_within_condition

data_condition = [];
% [] for pooled all-condition mode, or e.g. 1:16 for condition mode.

runIdx = 1;

% Used only in condition mode to get G-S-L / G-L-H / etc.
dat_file = fullfile('.', 'model_data_allruns');
stim_tag = '_2[Gpl2_2c_2sz_400_2_200isi]';

% Save switches
save_fig = true;
save_png = true;
close_after_save = false;

% Figure options
figure_visible = 'on';
png_dpi = 300;
marker_size = 30;
use_marker_alpha = true;

% Violin options
violin_width = 0.28;
violin_face_alpha = 0.20;
show_violin_median = true;

% Individual neuron points overlaid on violin plots
show_violin_points = true;
violin_point_size = 3;
violin_point_jitter_width = 0.08;
violin_point_alpha = 0.35;

% Same colors as RF_analysis.m.
% Row 1 is probe 0 / group 1, row 2 is probe 1 / group 2.
probe_colors = [
    0.0000, 0.4470, 0.7410;
    0.8500, 0.3250, 0.0980
];

%% ------------------------------------------------------------------------
% Main setup
% -------------------------------------------------------------------------

if isempty(data_condition)
    use_condition_mode = false;
    condition_list = [];
    numConditions = 1;
else
    use_condition_mode = true;
    condition_list = data_condition(:)';
    numConditions = numel(condition_list);
end

stim_abbrev = cell(1, numConditions);

if use_condition_mode
    fprintf('Reading stimulus metadata from %s\n', dat_file);

    Sdata = loadMatFileFlexibleLocal(dat_file, 'model_data_allruns');

    if ~isfield(Sdata, 'model_data_allruns')
        error('%s does not contain model_data_allruns.', dat_file);
    end

    model_data_allruns = Sdata.model_data_allruns;
    all_run_tags = get_all_run_tagsLocal(model_data_allruns);

    run_idx = find(strcmp(all_run_tags, stim_tag));

    if isempty(run_idx)
        error('Requested stim_tag not found: %s', stim_tag);
    end

    if numel(run_idx) > 1
        error('Duplicate stim_tag found: %s', stim_tag);
    end

    if ~isfield(model_data_allruns{run_idx}, 'conditions_full')
        error('model_data_allruns{%d} is missing conditions_full.', run_idx);
    end

    condition_full = model_data_allruns{run_idx}.conditions_full;
    conditionMap = buildConditionSummaryMapLocal(condition_full, condition_list);
    stim_abbrev = {conditionMap.entries.panelCondShortLabel};
end

%% ------------------------------------------------------------------------
% Main loop
% -------------------------------------------------------------------------

for cond_i = 1:numConditions

    if use_condition_mode
        this_condition = condition_list(cond_i);
        baseDir = ['./FA_Dlag_', data_content, '_condition', num2str(this_condition)];
        titleLabel = stim_abbrev{cond_i};
    else
        this_condition = [];
        baseDir = ['./FA_Dlag_', data_content];
        titleLabel = 'all conditions';
    end

    tempfname = sprintf('%s/mat_results/run%03d', baseDir, runIdx);
    r2File = fullfile(tempfname, 'reconstruction_R2.mat');

    fprintf('\n============================================================\n');

    if isempty(this_condition)
        fprintf('Plotting reconstruction R2: pooled all-condition mode\n');
    else
        fprintf('Plotting reconstruction R2: %s\n', titleLabel);
    end

    fprintf('Reading R2 file: %s\n', r2File);

    if ~exist(r2File, 'file')
        error('R2 file not found: %s. Run data_reconstruction.m first.', r2File);
    end

    S = load(r2File, 'recon_R2');

    if ~isfield(S, 'recon_R2')
        error('%s does not contain recon_R2.', r2File);
    end

    recon_R2 = S.recon_R2;

    checkReconR2FieldsLocal(recon_R2);

    fig = plotOneReconstructionR2FigureLocal( ...
        recon_R2, ...
        titleLabel, ...
        probe_colors, ...
        marker_size, ...
        use_marker_alpha, ...
        figure_visible, ...
        violin_width, ...
        violin_face_alpha, ...
        show_violin_median, ...
        show_violin_points, ...
        violin_point_size, ...
        violin_point_jitter_width, ...
        violin_point_alpha);

    if save_fig
        figFile = fullfile(tempfname, 'reconstruction_R2_comparison.fig');
        savefig(fig, figFile);
        fprintf('Saved FIG: %s\n', figFile);
    end

    if save_png
        pngFile = fullfile(tempfname, 'reconstruction_R2_comparison.png');
        savePngLocal(fig, pngFile, png_dpi);
        fprintf('Saved PNG: %s\n', pngFile);
    end

    if close_after_save
        close(fig);
    end
end

fprintf('\nDone.\n');

%% ========================================================================
% Local functions
% ========================================================================

function fig = plotOneReconstructionR2FigureLocal( ...
    recon_R2, titleLabel, probe_colors, marker_size, use_marker_alpha, ...
    figure_visible, violin_width, violin_face_alpha, show_violin_median, ...
    show_violin_points, violin_point_size, violin_point_jitter_width, ...
    violin_point_alpha)

scatterSpecs = {
    % tile index   x field              x label          y field                 y label
    1,             'use_all',           'use all',        'use_across',           'use across';
    2,             'use_all',           'use all',        'across_excl_within',   'across excl within';
    3,             'use_across',        'use across',     'across_excl_within',   'across excl within';
    5,             'use_all',           'use all',        'use_within',           'use within';
    6,             'use_all',           'use all',        'within_excl_across',   'within excl across';
    7,             'use_within',        'use within',     'within_excl_across',   'within excl across'
};

numScatterPlots = size(scatterSpecs, 1);

x0_by_group = recon_R2.use_all.neuron_by_group;
numGroups = numel(x0_by_group);

if size(probe_colors, 1) < numGroups
    warning(['probe_colors has fewer rows than number of groups. ', ...
        'Extra groups will be plotted in black.']);
end

% Collect all finite R2 values for scatter axes, so all scatter subplots
% use the same x/y limits.
allR2Vals = [];

for p = 1:numScatterPlots
    xField = scatterSpecs{p, 2};
    yField = scatterSpecs{p, 4};

    x_by_group = recon_R2.(xField).neuron_by_group;
    y_by_group = recon_R2.(yField).neuron_by_group;

    if numel(x_by_group) ~= numGroups
        error('R2 field %s has %d groups, but use_all has %d groups.', ...
            xField, numel(x_by_group), numGroups);
    end

    if numel(y_by_group) ~= numGroups
        error('R2 field %s has %d groups, but use_all has %d groups.', ...
            yField, numel(y_by_group), numGroups);
    end

    for g = 1:numGroups
        allR2Vals = [allR2Vals; x_by_group{g}(:); y_by_group{g}(:)]; %#ok<AGROW>
    end
end

allR2Vals = allR2Vals(isfinite(allR2Vals));
r2LimVals = paddedLimitsLocal(allR2Vals, [-1, 1]);

% Build delta R2 values for violin plots.
deltaAcross = computeDeltaR2ByGroupLocal( ...
    recon_R2, 'use_across', 'across_excl_within');

deltaWithin = computeDeltaR2ByGroupLocal( ...
    recon_R2, 'use_within', 'within_excl_across');

allDeltaVals = [];
for c = 1:3
    for g = 1:numGroups
        allDeltaVals = [allDeltaVals; deltaAcross{c, g}(:); deltaWithin{c, g}(:)]; %#ok<AGROW>
    end
end

allDeltaVals = allDeltaVals(isfinite(allDeltaVals));
deltaLimVals = paddedLimitsLocal(allDeltaVals, [-1, 1]);

fig = figure( ...
    'Name', 'reconstruction R2 comparison', ...
    'Color', 'w', ...
    'Visible', figure_visible, ...
    'Position', [100, 100, 1850, 850]);

t = tiledlayout(fig, 2, 4, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

% -------------------------------------------------------------------------
% Scatter plots
% -------------------------------------------------------------------------

for p = 1:numScatterPlots

    tileIdx = scatterSpecs{p, 1};
    xField = scatterSpecs{p, 2};
    xLabelText = scatterSpecs{p, 3};
    yField = scatterSpecs{p, 4};
    yLabelText = scatterSpecs{p, 5};

    x_by_group = recon_R2.(xField).neuron_by_group;
    y_by_group = recon_R2.(yField).neuron_by_group;

    ax = nexttile(t, tileIdx);
    hold(ax, 'on');

    % Diagonal reference line. Do not show this line in legend.
    plot(ax, r2LimVals, r2LimVals, 'k--', ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off');

    % Zero reference lines, only shown when the axis range contains zero.
    % Do not show these lines in legend.
    if r2LimVals(1) < 0 && r2LimVals(2) > 0
        plot(ax, [0 0], r2LimVals, 'k--', ...
            'LineWidth', 1, ...
            'HandleVisibility', 'off');

        plot(ax, r2LimVals, [0 0], 'k--', ...
            'LineWidth', 1, ...
            'HandleVisibility', 'off');
    end

    for g = 1:numGroups

        x = x_by_group{g}(:);
        y = y_by_group{g}(:);

        if numel(x) ~= numel(y)
            error(['Group %d has mismatched neuron counts between %s ', ...
                'and %s: %d vs %d.'], ...
                g, xField, yField, numel(x), numel(y));
        end

        valid = isfinite(x) & isfinite(y);
        x = x(valid);
        y = y(valid);

        if g <= size(probe_colors, 1)
            thisColor = probe_colors(g, :);
        else
            thisColor = [0, 0, 0];
        end

        if isempty(x)
            continue;
        end

        if use_marker_alpha
            scatter(ax, x, y, marker_size, ...
                'MarkerFaceColor', thisColor, ...
                'MarkerEdgeColor', thisColor, ...
                'MarkerFaceAlpha', 0.65, ...
                'MarkerEdgeAlpha', 0.65, ...
                'DisplayName', sprintf('Group %d neurons', g));
        else
            scatter(ax, x, y, marker_size, ...
                'MarkerFaceColor', thisColor, ...
                'MarkerEdgeColor', thisColor, ...
                'DisplayName', sprintf('Group %d neurons', g));
        end
    end

    xlim(ax, r2LimVals);
    ylim(ax, r2LimVals);
    axis(ax, 'square');

    xlabel(ax, sprintf('%s R^2', xLabelText), 'Interpreter', 'none');
    ylabel(ax, sprintf('%s R^2', yLabelText), 'Interpreter', 'none');
    title(ax, sprintf('%s vs %s', yLabelText, xLabelText), 'Interpreter', 'none');

    cleanAxisLocal(ax);

    if p == 1
        legend(ax, 'Location', 'best', 'Interpreter', 'none', 'Box', 'off');
    end
end

% -------------------------------------------------------------------------
% Violin plots with individual points
% -------------------------------------------------------------------------

axAcross = nexttile(t, 4);
plotDeltaViolinLocal( ...
    axAcross, deltaAcross, probe_colors, deltaLimVals, ...
    violin_width, violin_face_alpha, show_violin_median, ...
    show_violin_points, violin_point_size, violin_point_jitter_width, ...
    violin_point_alpha, ...
    'across delta R^2');

axWithin = nexttile(t, 8);
plotDeltaViolinLocal( ...
    axWithin, deltaWithin, probe_colors, deltaLimVals, ...
    violin_width, violin_face_alpha, show_violin_median, ...
    show_violin_points, violin_point_size, violin_point_jitter_width, ...
    violin_point_alpha, ...
    'within delta R^2');

sgtitle(t, sprintf('%s | reconstruction R2 comparison', titleLabel), ...
    'Interpreter', 'none');
end

function deltaVals = computeDeltaR2ByGroupLocal(recon_R2, useField, exclField)
% Compute three delta R2 types by group:
%   1) all-use  = use_all - use
%   2) all-excl = use_all - excl
%   3) use-excl = use     - excl
%
% Output:
%   deltaVals{deltaType, groupIdx}

all_by_group = recon_R2.use_all.neuron_by_group;
use_by_group = recon_R2.(useField).neuron_by_group;
excl_by_group = recon_R2.(exclField).neuron_by_group;

numGroups = numel(all_by_group);
deltaVals = cell(3, numGroups);

for g = 1:numGroups

    r2_all = all_by_group{g}(:);
    r2_use = use_by_group{g}(:);
    r2_excl = excl_by_group{g}(:);

    if numel(r2_all) ~= numel(r2_use) || numel(r2_all) ~= numel(r2_excl)
        error('Group %d has mismatched neuron counts for delta R2 computation.', g);
    end

    deltaVals{1, g} = r2_all - r2_use;
    deltaVals{2, g} = r2_all - r2_excl;
    deltaVals{3, g} = r2_use - r2_excl;
end
end

function plotDeltaViolinLocal( ...
    ax, deltaVals, probe_colors, deltaLimVals, violin_width, ...
    violin_face_alpha, show_violin_median, show_violin_points, ...
    violin_point_size, violin_point_jitter_width, violin_point_alpha, ...
    plotTitle)

hold(ax, 'on');

numDeltaTypes = size(deltaVals, 1);
numGroups = size(deltaVals, 2);

deltaLabels = {'all-use', 'all-excl', 'use-excl'};

if numDeltaTypes ~= numel(deltaLabels)
    error('Expected three delta R2 types.');
end

if numGroups == 1
    groupOffsets = 0;
else
    groupOffsets = linspace(-0.16, 0.16, numGroups);
end

% Horizontal zero line if needed.
if deltaLimVals(1) < 0 && deltaLimVals(2) > 0
    plot(ax, [0.5, numDeltaTypes + 0.5], [0, 0], 'k--', ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off');
end

for c = 1:numDeltaTypes
    for g = 1:numGroups

        vals = deltaVals{c, g};
        vals = vals(:);
        vals = vals(isfinite(vals));

        if isempty(vals)
            continue;
        end

        if g <= size(probe_colors, 1)
            thisColor = probe_colors(g, :);
        else
            thisColor = [0, 0, 0];
        end

        xPos = c + groupOffsets(g);

        drawOneViolinLocal( ...
            ax, vals, xPos, violin_width, thisColor, ...
            violin_face_alpha, show_violin_median, ...
            show_violin_points, violin_point_size, ...
            violin_point_jitter_width, violin_point_alpha);
    end
end

xlim(ax, [0.5, numDeltaTypes + 0.5]);
ylim(ax, deltaLimVals);

set(ax, ...
    'XTick', 1:numDeltaTypes, ...
    'XTickLabel', deltaLabels);

ylabel(ax, '\Delta R^2', 'Interpreter', 'tex');
title(ax, plotTitle, 'Interpreter', 'none');

cleanAxisLocal(ax);
end

function drawOneViolinLocal( ...
    ax, vals, xPos, violin_width, thisColor, violin_face_alpha, ...
    showMedian, showPoints, pointSize, pointJitterWidth, pointAlpha)

vals = vals(:);
vals = vals(isfinite(vals));

if isempty(vals)
    return;
end

% Individual neuron points with density-aware horizontal jitter.
% Points near high-density y values spread wider, so the point cloud
% follows the violin shape.
if showPoints
    xJitter = densityAwareJitterLocal(vals, xPos, pointJitterWidth);

    try
        scatter(ax, xJitter, vals, pointSize, ...
            'MarkerFaceColor', thisColor, ...
            'MarkerEdgeColor', thisColor, ...
            'MarkerFaceAlpha', pointAlpha, ...
            'MarkerEdgeAlpha', pointAlpha, ...
            'HandleVisibility', 'off');
    catch
        scatter(ax, xJitter, vals, pointSize, ...
            'MarkerFaceColor', thisColor, ...
            'MarkerEdgeColor', thisColor, ...
            'HandleVisibility', 'off');
    end
end

if numel(vals) < 2 || max(vals) == min(vals)
    y0 = vals(1);
    plot(ax, [xPos - violin_width * 0.45, xPos + violin_width * 0.45], ...
        [y0, y0], '-', ...
        'Color', thisColor, ...
        'LineWidth', 2, ...
        'HandleVisibility', 'off');

    return;
end

[f, xi] = estimateDensityLocal(vals);

if isempty(f) || isempty(xi) || max(f) <= 0
    medVal = median(vals);
    plot(ax, [xPos - violin_width * 0.45, xPos + violin_width * 0.45], ...
        [medVal, medVal], '-', ...
        'Color', thisColor, ...
        'LineWidth', 2, ...
        'HandleVisibility', 'off');
    return;
end

f = f(:)';
xi = xi(:)';

f = f ./ max(f) .* violin_width;

xPatch = [xPos - f, fliplr(xPos + f)];
yPatch = [xi, fliplr(xi)];

patch(ax, xPatch, yPatch, thisColor, ...
    'FaceAlpha', violin_face_alpha, ...
    'EdgeColor', thisColor, ...
    'LineWidth', 1, ...
    'HandleVisibility', 'off');

if showMedian
    medVal = median(vals);
    plot(ax, [xPos - violin_width * 0.55, xPos + violin_width * 0.55], ...
        [medVal, medVal], '-', ...
        'Color', thisColor, ...
        'LineWidth', 2, ...
        'HandleVisibility', 'off');
end
end

function [f, xi] = estimateDensityLocal(vals)
% Estimate 1-D density for violin plots.
% Prefer ksdensity if available. Fall back to histogram density otherwise.

vals = vals(:);
vals = vals(isfinite(vals));

if isempty(vals)
    f = [];
    xi = [];
    return;
end

if numel(vals) < 2 || max(vals) == min(vals)
    f = [];
    xi = [];
    return;
end

if exist('ksdensity', 'file') == 2
    try
        [f, xi] = ksdensity(vals, 'NumPoints', 100);
        return;
    catch
        % Fall through to histogram fallback.
    end
end

numBins = min(20, max(5, round(sqrt(numel(vals)))));

try
    [counts, edges] = histcounts(vals, numBins, 'Normalization', 'pdf');
catch
    [counts, edges] = histcounts(vals, numBins);
    binWidth = mean(diff(edges));
    counts = counts ./ (sum(counts) * binWidth);
end

xi = edges(1:end-1) + diff(edges) / 2;
f = counts;

valid = isfinite(f) & isfinite(xi);
f = f(valid);
xi = xi(valid);
end

function limVals = paddedLimitsLocal(vals, defaultLim)
% Compute padded limits from finite values.

vals = vals(:);
vals = vals(isfinite(vals));

if isempty(vals)
    limVals = defaultLim;
    return;
end

minVal = min(vals);
maxVal = max(vals);

if minVal == maxVal
    padVal = max(0.1, abs(minVal) * 0.1);
else
    padVal = 0.05 * (maxVal - minVal);
end

limVals = [minVal - padVal, maxVal + padVal];

if limVals(1) == limVals(2)
    limVals = limVals + [-0.1, 0.1];
end
end

function cleanAxisLocal(ax)

grid(ax, 'off');
box(ax, 'off');

set(ax, ...
    'TickDir', 'out', ...
    'LineWidth', 1, ...
    'FontSize', 11);
end

function checkReconR2FieldsLocal(recon_R2)

requiredTopFields = {
    'use_all'
    'use_across'
    'across_excl_within'
    'use_within'
    'within_excl_across'
};

for i = 1:numel(requiredTopFields)
    f = requiredTopFields{i};

    if ~isfield(recon_R2, f)
        error('recon_R2 is missing field %s.', f);
    end

    if ~isfield(recon_R2.(f), 'neuron_by_group')
        error('recon_R2.%s is missing neuron_by_group.', f);
    end

    if ~iscell(recon_R2.(f).neuron_by_group)
        error('recon_R2.%s.neuron_by_group must be a cell array.', f);
    end
end
end

function savePngLocal(fig, pngFile, dpi)

drawnow;

try
    if exist('exportgraphics', 'file') == 2
        exportgraphics(fig, pngFile, ...
            'Resolution', dpi, ...
            'BackgroundColor', 'white');
    else
        set(fig, 'PaperPositionMode', 'auto');
        print(fig, pngFile, '-dpng', sprintf('-r%d', dpi));
    end
catch
    warning('exportgraphics/print failed. Falling back to saveas.');
    saveas(fig, pngFile);
end
end

function S = loadMatFileFlexibleLocal(fileBase, varargin)
% Load a .mat file while allowing fileBase with or without .mat extension.

if exist(fileBase, 'file') == 2
    S = load(fileBase, varargin{:});
    return;
end

if ~endsWith(fileBase, '.mat')
    fileMat = [fileBase, '.mat'];

    if exist(fileMat, 'file') == 2
        S = load(fileMat, varargin{:});
        return;
    end
end

error('File not found: %s or %s.mat', fileBase, fileBase);
end

function all_tags = get_all_run_tagsLocal(model_data_allruns)

all_tags = cell(numel(model_data_allruns), 1);

for j = 1:numel(model_data_allruns)
    if ~isfield(model_data_allruns{j}, 'stim_tag')
        error('stim_tag missing in model_data_allruns{%d}.', j);
    end

    all_tags{j} = model_data_allruns{j}.stim_tag;
end
end

function conditionMap = buildConditionSummaryMapLocal(condition_full, condition_list)
% Build short labels such as G-S-L, G-S-H, G-L-L, etc.
%
% This follows the same condition-label convention used in
% subspace_similarity_dlag.m.

if isempty(condition_full)
    error('condition_full is empty.');
end

nAll = numel(condition_full);

stimNameAll = strings(nAll, 1);
sizeAll = nan(nAll, 1);
contrastAll = nan(nAll, 1);
effDirAll = nan(nAll, 1);

for k = 1:nAll

    if ~isfield(condition_full(k), 'stim_name')
        error('condition_full(%d) missing field stim_name.', k);
    end

    if ~isfield(condition_full(k), 'size')
        error('condition_full(%d) missing field size.', k);
    end

    if ~isfield(condition_full(k), 'contrast')
        error('condition_full(%d) missing field contrast.', k);
    end

    currStim = lower(string(condition_full(k).stim_name));

    stimNameAll(k) = currStim;
    sizeAll(k) = condition_full(k).size;
    contrastAll(k) = condition_full(k).contrast;

    if currStim == "plaid"
        if ~isfield(condition_full(k), 'plaid_dir')
            error('condition_full(%d) missing field plaid_dir.', k);
        end

        effDirAll(k) = condition_full(k).plaid_dir;

    elseif currStim == "grating"
        if ~isfield(condition_full(k), 'grating_dir')
            error('condition_full(%d) missing field grating_dir.', k);
        end

        effDirAll(k) = condition_full(k).grating_dir;

    else
        error('Unsupported stim_name in condition_full(%d): %s', ...
            k, char(currStim));
    end
end

allStim = unique(stimNameAll, 'stable');
allStim = lower(allStim);

if all(ismember(["grating", "plaid"], allStim))
    stimLabels = ["grating", "plaid"];
else
    if numel(allStim) ~= 2
        error('Expected exactly 2 stim levels in condition_full.');
    end

    stimLabels = allStim(:)';
end

sizeVals = unique(sizeAll);
sizeVals = sort(sizeVals(:)');

if numel(sizeVals) ~= 2
    error('Expected exactly 2 size levels in condition_full.');
end

contrastValuesByStim = struct();

for s = 1:2
    idx = (stimNameAll == stimLabels(s));
    cvals = unique(contrastAll(idx));
    cvals = sort(cvals(:)');

    if numel(cvals) ~= 2
        error('Stim %s does not have exactly 2 contrast levels.', ...
            char(stimLabels(s)));
    end

    contrastValuesByStim.(char(stimLabels(s))) = cvals;
end

dirVals = unique(effDirAll);
dirVals = sort(dirVals(:)');

if numel(dirVals) ~= 2
    error('Expected exactly 2 effective direction values in condition_full.');
end

stimDirLabels = {'stim_dir1', 'stim_dir2'};

condLabels = {
    'grating-small-low',  'grating-small-high', ...
    'grating-large-low',  'grating-large-high', ...
    'plaid-small-low',    'plaid-small-high', ...
    'plaid-large-low',    'plaid-large-high'
};

condShortLabels = {
    'G-S-L', 'G-S-H', 'G-L-L', 'G-L-H', ...
    'P-S-L', 'P-S-H', 'P-L-L', 'P-L-H'
};

entries = struct([]);

for ii = 1:numel(condition_list)

    condID = condition_list(ii);

    if condID < 1 || condID > nAll
        error('Condition ID %d is outside condition_full range.', condID);
    end

    currStim = lower(string(condition_full(condID).stim_name));
    currSize = condition_full(condID).size;
    currContrast = condition_full(condID).contrast;

    if currStim == "plaid"
        currDir = condition_full(condID).plaid_dir;
    else
        currDir = condition_full(condID).grating_dir;
    end

    stimCode = find(strcmp(cellstr(stimLabels), char(currStim)), 1);
    sizeCode = find(sizeVals == currSize, 1);

    currContrastLevels = contrastValuesByStim.(char(currStim));
    contrastCode = find(currContrastLevels == currContrast, 1);

    stimDirCode = find(dirVals == currDir, 1);

    if isempty(stimCode) || isempty(sizeCode) || isempty(contrastCode) || isempty(stimDirCode)
        error('Could not map condition ID %d to summary label.', condID);
    end

    panelCondIndex = (stimCode - 1) * 4 + (sizeCode - 1) * 2 + contrastCode;

    entries(ii).conditionId = condID;
    entries(ii).stimName = char(currStim);
    entries(ii).stimCode = stimCode;
    entries(ii).sizeValue = currSize;
    entries(ii).sizeCode = sizeCode;
    entries(ii).sizeLabel = ternary_labelLocal(sizeCode, 'small', 'large');
    entries(ii).contrastValue = currContrast;
    entries(ii).contrastCode = contrastCode;
    entries(ii).contrastLabel = ternary_labelLocal(contrastCode, 'low', 'high');
    entries(ii).stimDirValue = currDir;
    entries(ii).stimDirCode = stimDirCode;
    entries(ii).stimDirLabel = stimDirLabels{stimDirCode};
    entries(ii).panelCondIndex = panelCondIndex;
    entries(ii).panelCondLabel = condLabels{panelCondIndex};
    entries(ii).panelCondShortLabel = condShortLabels{panelCondIndex};
end

conditionMap = struct();
conditionMap.entries = entries;
conditionMap.meta.stimLabels = cellstr(stimLabels);
conditionMap.meta.sizeValues = sizeVals;
conditionMap.meta.contrastValuesByStim = contrastValuesByStim;
conditionMap.meta.stimDirLabels = stimDirLabels;
conditionMap.meta.stimDirValues = dirVals;
conditionMap.meta.panelCondLabels = condLabels;
conditionMap.meta.panelCondShortLabels = condShortLabels;
end

function out = ternary_labelLocal(code, label1, label2)

if isempty(code) || ~isfinite(code)
    out = '';
elseif code == 1
    out = label1;
else
    out = label2;
end
end


function xJitter = densityAwareJitterLocal(vals, xPos, maxJitterWidth)
% Density-aware jitter for violin-overlaid points.
%
% vals near high-density y regions get larger horizontal jitter.
% vals near low-density y regions stay close to xPos.

vals = vals(:);

if isempty(vals)
    xJitter = [];
    return;
end

if numel(vals) < 2 || max(vals) == min(vals)
    xJitter = xPos + zeros(size(vals));
    return;
end

[f, xi] = estimateDensityLocal(vals);

if isempty(f) || isempty(xi) || max(f) <= 0
    xJitter = xPos + (rand(size(vals)) - 0.5) * 2 * maxJitterWidth;
    return;
end

f = f(:);
xi = xi(:);

% Interpolate density at each observed point.
densityAtVals = interp1(xi, f, vals, 'linear', 'extrap');
densityAtVals(~isfinite(densityAtVals)) = 0;
densityAtVals(densityAtVals < 0) = 0;

if max(densityAtVals) > 0
    densityAtVals = densityAtVals ./ max(densityAtVals);
else
    densityAtVals = ones(size(vals));
end

% Optional small floor, so low-density points are not exactly on center.
densityFloor = 0.08;
densityAtVals = densityFloor + (1 - densityFloor) .* densityAtVals;

localWidth = maxJitterWidth .* densityAtVals;

xJitter = xPos + (rand(size(vals)) - 0.5) .* 2 .* localWidth;
end