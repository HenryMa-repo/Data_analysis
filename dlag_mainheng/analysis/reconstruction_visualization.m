%% reconstruction_visualization.m
% Visualize original/reconstructed neural responses after data_reconstruction.m.
%
% This script plots trial-by-time heatmaps for randomly selected example
% neurons. It is intended to be run after data_reconstruction.m has updated
% bestmodel*.mat with seqEst reconstruction fields.
%
% Key behavior:
%   1) analysis_fields can contain one or multiple seqEst fields.
%   2) For each group, the same randomly selected neurons are used for all
%      analysis_fields.
%   3) For each group, the heatmap color limit is shared across all selected
%      analysis_fields, but is not shared across groups.
%   4) Each output figure is one group and one analysis field.
%   5) Each selected neuron is one subplot in the figure.
%   6) x-axis is time, y-axis is condition-sorted trials.
%   7) Condition labels are shown only on the first neuron subplot.
%
% Expected saved data:
%   - model_data_prepar.m creates trial structs with fields trialId, T, y.
%   - model_data_prepar.m saves condition_index_per_trial_full and
%     conditions_full in model_data_allruns.
%   - data_reconstruction.m adds yRecon_* fields to seqEst and overwrites
%     the corresponding bestmodel*.mat.

clc;
clear;

%% ------------------------------------------------------------------------
% User parameters
% -------------------------------------------------------------------------

data_content = 'raw_count';
% Options usually include:
%   raw_count
%   raw_fr
%   z_within_trial
%   z_within_condition
%   z_across_conditions
%   demean_count_within_trial
%   demean_fr_within_trial
%   demean_pooledsd_within_condition

% [] means pooled all-condition model.
% Example: 1:16 means condition-specific models, one model per condition.
data_condition = [];

runIdx = 1;

% Used to map trialId back to condition and to get condition labels such as
% G-S-L, G-L-H, P-S-L, etc.
dat_file = fullfile('.', 'model_data_allruns');
stim_tag = '_2[Gpl2_2c_2sz_400_2_200isi]';

% Fields to plot from seqEst. You can choose one or multiple fields.
% Examples:
%   analysis_fields = {'y'};
%   analysis_fields = {'y', 'yRecon_use_all'};
%   analysis_fields = {'y', 'yRecon_use_all', 'yRecon_use_across'};
%
% Available fields after data_reconstruction.m:
%
%   Original data:
%     y
%
%   Base noiseless reconstructions:
%     yRecon_use_across
%     yRecon_use_within
%     yRecon_use_all
%     yRecon_across_excl_within
%     yRecon_within_excl_across
%
%   If add_R_noise_reconstruction = true in data_reconstruction.m:
%     yRecon_use_across_with_R
%     yRecon_use_within_with_R
%     yRecon_use_all_with_R
%     yRecon_across_excl_within_with_R
%     yRecon_within_excl_across_with_R
%
%   If add_keep_resid_reconstruction = true 
%     yRecon_use_across_keep_resid
%     yRecon_use_within_keep_resid
%     yRecon_use_all_keep_resid
%     yRecon_across_excl_within_keep_resid
%     yRecon_within_excl_across_keep_resid

% if add_directional_reconstruction = true in data_reconstruction.m:
% -------------------------------------------------------------------------
%   yRecon_use_feedback
%   yRecon_feedback_excl_within_ff_ambiguous
%   yRecon_feedback_excl_within
%   yRecon_feedback_excl_ff_ambiguous
%
%   yRecon_use_feedforward
%   yRecon_feedforward_excl_within_fb_ambiguous
%   yRecon_feedforward_excl_within
%   yRecon_feedforward_excl_fb_ambiguous
%
analysis_fields = {'yRecon_use_feedback'};

% Number of randomly selected example neurons per group.
n_example_neurons = 10;

% Random selection control. The same selected neurons are used for all
% analysis_fields within the same group.
use_fixed_random_seed = true;
random_seed = 1;

% Heatmap color scale. The color limit is shared across analysis_fields
% within each group, but not shared across groups.
% Percentile limits are more robust than min/max against outliers.
color_percentiles = [1 99];

% Visualization options.
condition_gap_rows = 1;          % Empty rows inserted between condition blocks.
draw_condition_boxes = true;     % Draw a thin box around each condition block.
draw_condition_separators = true;
colormap_name = 'parula';
figure_visible = 'on';
max_tiles_per_row = 5;

% Save switches.
save_fig = true;
save_png = true;
close_after_save = true;
png_dpi = 300;

%% ------------------------------------------------------------------------
% Main setup
% -------------------------------------------------------------------------

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

analysis_fields = normalizeFieldListLocal(analysis_fields);

if isempty(data_condition)
    use_condition_mode = false;
    condition_list = [];
    modeTag = 'all-condition-model';
else
    use_condition_mode = true;
    condition_list = reshape(data_condition, 1, []);
    modeTag = 'condition-specific-models';
end

if use_fixed_random_seed
    rng(random_seed, 'twister');
else
    rng('shuffle');
end

fprintf('Reading stimulus metadata from %s\n', dat_file);
Sdata = loadMatFileFlexibleLocal(dat_file, 'model_data_allruns');
if ~isfield(Sdata, 'model_data_allruns')
    error('%s does not contain model_data_allruns.', dat_file);
end

model_data_allruns = Sdata.model_data_allruns;
all_run_tags = getAllRunTagsLocal(model_data_allruns);
run_idx = find(strcmp(all_run_tags, stim_tag));

if isempty(run_idx)
    error('Requested stim_tag not found: %s', stim_tag);
end
if numel(run_idx) > 1
    error('Duplicate stim_tag found: %s', stim_tag);
end

this_run = model_data_allruns{run_idx};

if ~isfield(this_run, 'conditions_full')
    error('model_data_allruns{%d} is missing conditions_full.', run_idx);
end
if ~isfield(this_run, 'condition_index_per_trial_full')
    error('model_data_allruns{%d} is missing condition_index_per_trial_full.', run_idx);
end

condition_full = this_run.conditions_full;
condition_index_per_trial_full = this_run.condition_index_per_trial_full(:);

%% ------------------------------------------------------------------------
% Load model data and organize trials into condition blocks
% -------------------------------------------------------------------------

if use_condition_mode
    dataBlocks = repmat(makeEmptyBlockLocal(), 1, numel(condition_list));
    yDims_ref = [];

    condition_labels = getConditionLabelsLocal(condition_full, condition_list);

    for cond_i = 1:numel(condition_list)
        this_condition = condition_list(cond_i);
        baseDir = ['./FA_Dlag_', data_content, '_condition', num2str(this_condition)];
        tempfname = sprintf('%s/mat_results/run%03d', baseDir, runIdx);

        fprintf('\n============================================================\n');
        fprintf('Loading condition-specific model: condition %d\n', this_condition);
        fprintf('Reading from %s\n', tempfname);

        modelData = loadBestModelDataLocal(tempfname);

        if isempty(yDims_ref)
            yDims_ref = modelData.yDims;
        elseif ~isequal(yDims_ref(:)', modelData.yDims(:)')
            error('yDims mismatch between condition-specific models.');
        end

        seqThis = sortSeqByTrialIdLocal(modelData.seqEst);
        warnIfConditionMismatchLocal(seqThis, this_condition, condition_index_per_trial_full);

        dataBlocks(cond_i).conditionId = this_condition;
        dataBlocks(cond_i).label = condition_labels{cond_i};
        dataBlocks(cond_i).seqEst = seqThis;
        dataBlocks(cond_i).sourceFolder = tempfname;
    end

    yDims = yDims_ref;
    saveDir = scriptDir;
else
    baseDir = ['./FA_Dlag_', data_content];
    tempfname = sprintf('%s/mat_results/run%03d', baseDir, runIdx);

    fprintf('\n============================================================\n');
    fprintf('Loading pooled all-condition model\n');
    fprintf('Reading from %s\n', tempfname);

    modelData = loadBestModelDataLocal(tempfname);
    yDims = modelData.yDims;

    condition_order = getConditionOrderFromSeqLocal( ...
        modelData.seqEst, condition_index_per_trial_full);

    condition_labels = getConditionLabelsLocal(condition_full, condition_order);

    dataBlocks = splitSeqByConditionLocal( ...
        modelData.seqEst, condition_order, condition_labels, ...
        condition_index_per_trial_full);

    saveDir = tempfname;
end

if isempty(yDims) || any(yDims <= 0)
    error('Invalid yDims.');
end

numGroups = numel(yDims);

fprintf('\nAnalysis fields to plot:\n');
for f = 1:numel(analysis_fields)
    fprintf('  %s\n', analysis_fields{f});
end
fprintf('Output folder: %s\n', saveDir);

%% ------------------------------------------------------------------------
% Plot one figure per group per analysis field
% -------------------------------------------------------------------------

for groupIdx = 1:numGroups
    groupRows = getGroupRowsLocal(yDims, groupIdx);
    nAvailable = numel(groupRows);

    nSelect = min(n_example_neurons, nAvailable);
    if nSelect < n_example_neurons
        warning('Group %d only has %d neurons. Plotting all available neurons.', ...
            groupIdx, nAvailable);
    end

    selectedLocalNeuronIds = randperm(nAvailable, nSelect);
    selectedGlobalRows = groupRows(selectedLocalNeuronIds);

    fprintf('\n============================================================\n');
    fprintf('Group %d: selected local neurons %s\n', ...
        groupIdx, mat2str(selectedLocalNeuronIds));

    allVals = collectValuesForColorLimitLocal( ...
        dataBlocks, analysis_fields, selectedGlobalRows);

    climVals = robustColorLimitsLocal(allVals, color_percentiles);

    fprintf('Group %d shared color limit across fields: [%g, %g]\n', ...
        groupIdx, climVals(1), climVals(2));

    for f = 1:numel(analysis_fields)
        fieldName = analysis_fields{f};

        figTitle = sprintf('%s_%s_%s_group%d', ...
            data_content, fieldName, modeTag, groupIdx);

        fig = plotOneGroupOneFieldLocal( ...
            dataBlocks, fieldName, selectedGlobalRows, selectedLocalNeuronIds, ...
            figTitle, this_run, climVals, condition_gap_rows, ...
            draw_condition_boxes, draw_condition_separators, colormap_name, ...
            figure_visible, max_tiles_per_row);

        fileBase = sprintf('%s_%dneuron-example', figTitle, nSelect);
        fileBase = sanitizeFileNameLocal(fileBase);

        if save_fig
            figFile = fullfile(saveDir, [fileBase, '.fig']);
            saveFigLocal(fig, figFile);
            fprintf('Saved FIG: %s\n', figFile);
        end

        if save_png
            pngFile = fullfile(saveDir, [fileBase, '.png']);
            savePngLocal(fig, pngFile, png_dpi);
            fprintf('Saved PNG: %s\n', pngFile);
        end

        if close_after_save
            close(fig);
        end
    end
end

fprintf('\nDone.\n');

%% ========================================================================
% Local functions
% ========================================================================

function fields = normalizeFieldListLocal(fields)
    if ischar(fields)
        fields = {fields};
    elseif isstring(fields)
        fields = cellstr(fields(:));
    elseif iscell(fields)
        fields = fields(:);
        for i = 1:numel(fields)
            if isstring(fields{i})
                fields{i} = char(fields{i});
            end
            if ~ischar(fields{i})
                error('analysis_fields{%d} must be a char or string.', i);
            end
        end
    else
        error('analysis_fields must be a char, string array, or cell array.');
    end

    fields = reshape(fields, 1, []);
    if isempty(fields)
        error('analysis_fields is empty.');
    end

    if numel(unique(fields, 'stable')) ~= numel(fields)
        warning('analysis_fields contains duplicates. Duplicates will be removed.');
        fields = unique(fields, 'stable');
    end
end

function block = makeEmptyBlockLocal()
    block = struct();
    block.conditionId = [];
    block.label = '';
    block.seqEst = [];
    block.sourceFolder = '';
end

function modelData = loadBestModelDataLocal(tempfname)
    bestFile = findOneFileLocal(tempfname, 'bestmodel*', true);
    fprintf('Loading best model: %s\n', bestFile);

    Sbest = load(bestFile);

    requiredVars = {'bestModel', 'res', 'seqEst'};
    for i = 1:numel(requiredVars)
        if ~isfield(Sbest, requiredVars{i})
            error('%s is missing variable %s.', bestFile, requiredVars{i});
        end
    end

    seqEst = Sbest.seqEst;
    if isempty(seqEst)
        error('seqEst is empty in %s.', bestFile);
    end

    yDims = getYDimsFromModelLocal(Sbest.bestModel, Sbest.res, seqEst, bestFile);

    modelData = struct();
    modelData.bestFile = bestFile;
    modelData.seqEst = seqEst;
    modelData.yDims = yDims;
end

function yDims = getYDimsFromModelLocal(bestModel, res, seqEst, sourceName)
    yDims = [];

    if isfield(res, 'estParams') && ...
            isfield(res.estParams, 'yDims') && ...
            ~isempty(res.estParams.yDims)
        yDims = res.estParams.yDims;
    elseif isfield(bestModel, 'yDims') && ~isempty(bestModel.yDims)
        yDims = bestModel.yDims;
    end

    if isempty(yDims)
        error('Could not determine yDims from %s.', sourceName);
    end

    yDims = reshape(double(yDims), 1, []);

    if ~isfield(seqEst, 'y')
        error('seqEst.y is missing in %s.', sourceName);
    end

    yDimObserved = size(seqEst(1).y, 1);
    if sum(yDims) ~= yDimObserved
        error('sum(yDims) = %d but size(seqEst(1).y,1) = %d in %s.', ...
            sum(yDims), yDimObserved, sourceName);
    end
end

function seqOut = sortSeqByTrialIdLocal(seqIn)
    if isempty(seqIn)
        seqOut = seqIn;
        return;
    end

    if isfield(seqIn, 'trialId')
        trialIds = arrayfun(@(s) s.trialId, seqIn);
        [~, ord] = sort(trialIds(:));
        seqOut = seqIn(ord);
    else
        seqOut = seqIn;
    end
end

function warnIfConditionMismatchLocal(seqEst, expectedCondition, condition_index_per_trial_full)
    if isempty(seqEst) || ~isfield(seqEst, 'trialId')
        warning('seqEst has no trialId. Cannot verify condition-specific trial membership.');
        return;
    end

    trialIds = arrayfun(@(s) s.trialId, seqEst);
    valid = trialIds >= 1 & trialIds <= numel(condition_index_per_trial_full);

    if ~all(valid)
        warning('Some trialId values are outside condition_index_per_trial_full.');
        return;
    end

    condIds = condition_index_per_trial_full(trialIds);
    if any(condIds ~= expectedCondition)
        warning(['Some trials in condition-specific model condition %d map to ', ...
            'different condition IDs according to trialId metadata.'], expectedCondition);
    end
end

function condition_order = getConditionOrderFromSeqLocal(seqEst, condition_index_per_trial_full)
    if ~isfield(seqEst, 'trialId')
        error(['seqEst is missing trialId. Cannot sort all-condition trials by condition. ', ...
            'model_data_prepar.m should have created trialId for each trial.']);
    end

    trialIds = arrayfun(@(s) s.trialId, seqEst);

    if any(trialIds < 1) || any(trialIds > numel(condition_index_per_trial_full))
        error('Some seqEst trialId values are outside condition_index_per_trial_full.');
    end

    condIds = condition_index_per_trial_full(trialIds);
    condition_order = unique(condIds(:)', 'stable');
    condition_order = sort(condition_order);
end

function dataBlocks = splitSeqByConditionLocal( ...
    seqEst, condition_order, condition_labels, condition_index_per_trial_full)

    dataBlocks = repmat(makeEmptyBlockLocal(), 1, numel(condition_order));

    trialIds = arrayfun(@(s) s.trialId, seqEst);
    condIds = condition_index_per_trial_full(trialIds);

    for i = 1:numel(condition_order)
        c = condition_order(i);
        keep = find(condIds == c);

        if isempty(keep)
            seqThis = seqEst([]);
        else
            [~, ord] = sort(trialIds(keep));
            seqThis = seqEst(keep(ord));
        end

        dataBlocks(i).conditionId = c;
        dataBlocks(i).label = condition_labels{i};
        dataBlocks(i).seqEst = seqThis;
        dataBlocks(i).sourceFolder = '';
    end
end

function vals = collectValuesForColorLimitLocal(dataBlocks, analysis_fields, selectedRows)
    vals = [];

    for f = 1:numel(analysis_fields)
        fieldName = analysis_fields{f};

        for b = 1:numel(dataBlocks)
            seqEst = dataBlocks(b).seqEst;

            if isempty(seqEst)
                continue;
            end

            checkSeqFieldLocal(seqEst, fieldName, selectedRows);

            for tr = 1:numel(seqEst)
                Y = double(seqEst(tr).(fieldName));
                vals = [vals; reshape(Y(selectedRows, :), [], 1)]; %#ok<AGROW>
            end
        end
    end

    vals = vals(:);
    vals = vals(isfinite(vals));

    if isempty(vals)
        error('No finite values found for color limit computation.');
    end
end

function checkSeqFieldLocal(seqEst, fieldName, selectedRows)
    if ~isfield(seqEst, fieldName)
        error('Field %s is missing from seqEst.', fieldName);
    end

    for tr = 1:numel(seqEst)
        Y = seqEst(tr).(fieldName);

        if ~isnumeric(Y) && ~islogical(Y)
            error('seqEst(%d).%s must be numeric or logical.', tr, fieldName);
        end

        if ndims(Y) ~= 2
            error('seqEst(%d).%s must be a 2-D matrix.', tr, fieldName);
        end

        if max(selectedRows) > size(Y, 1)
            error('seqEst(%d).%s has only %d rows, but selected row %d is requested.', ...
                tr, fieldName, size(Y, 1), max(selectedRows));
        end
    end
end

function climVals = robustColorLimitsLocal(vals, pct)
    vals = vals(:);
    vals = vals(isfinite(vals));

    if isempty(vals)
        error('Cannot compute color limits from empty values.');
    end

    pct = double(pct(:)');
    if numel(pct) ~= 2 || pct(1) < 0 || pct(2) > 100 || pct(1) >= pct(2)
        error('color_percentiles must be [low high], with 0 <= low < high <= 100.');
    end

    lo = percentileLocal(vals, pct(1));
    hi = percentileLocal(vals, pct(2));

    if ~isfinite(lo) || ~isfinite(hi)
        lo = min(vals);
        hi = max(vals);
    end

    if lo == hi
        padVal = max(1e-6, abs(lo) * 0.05);
        lo = lo - padVal;
        hi = hi + padVal;
    end

    climVals = [lo hi];
end

function q = percentileLocal(x, p)
    x = sort(x(:));
    x = x(isfinite(x));

    if isempty(x)
        q = NaN;
        return;
    end

    if p <= 0
        q = x(1);
        return;
    end

    if p >= 100
        q = x(end);
        return;
    end

    pos = 1 + (numel(x) - 1) * p / 100;
    lo = floor(pos);
    hi = ceil(pos);

    if lo == hi
        q = x(lo);
    else
        w = pos - lo;
        q = (1 - w) * x(lo) + w * x(hi);
    end
end

function fig = plotOneGroupOneFieldLocal( ...
    dataBlocks, fieldName, selectedGlobalRows, selectedLocalNeuronIds, ...
    figTitle, runMeta, climVals, condition_gap_rows, ...
    draw_condition_boxes, draw_condition_separators, colormap_name, ...
    figure_visible, max_tiles_per_row)

    nNeurons = numel(selectedGlobalRows);
    nCols = min(max_tiles_per_row, nNeurons);
    nRows = ceil(nNeurons / nCols);

    figWidth = max(900, 340 * nCols);
    figHeight = max(560, 340 * nRows);

    fig = figure( ...
        'Name', figTitle, ...
        'Color', 'w', ...
        'Visible', figure_visible, ...
        'Position', [100, 100, figWidth, figHeight]);

    t = tiledlayout(fig, nRows, nCols, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    lastAx = [];

    for n = 1:nNeurons
        ax = nexttile(t);
        lastAx = ax;

        [M, blockInfo] = buildNeuronHeatmapMatrixLocal( ...
            dataBlocks, fieldName, selectedGlobalRows(n), condition_gap_rows);

        [xValues, xLabelText] = getTimeAxisLocal(runMeta, size(M, 2));
        xEdges = estimateXEdgesLocal(xValues);

        h = imagesc(ax, xValues, 1:size(M, 1), M);
        set(ax, 'CLim', climVals);
        set(h, 'AlphaData', isfinite(M));

        set(ax, 'Color', [1 1 1]);
        set(ax, 'YDir', 'reverse');

        applyColormapLocal(ax, colormap_name);
        hold(ax, 'on');

        drawConditionGuidesLocal( ...
            ax, blockInfo, xEdges, draw_condition_boxes, draw_condition_separators);

        xlim(ax, [xEdges(1), xEdges(end)]);
        ylim(ax, [0.5, size(M, 1) + 0.5]);

        title(ax, sprintf('neuron %d', selectedLocalNeuronIds(n)), ...
            'Interpreter', 'none');

        xlabel(ax, xLabelText, 'Interpreter', 'none');

        if n == 1
            set(ax, ...
                'YTick', [blockInfo.centerRow], ...
                'YTickLabel', {blockInfo.label}, ...
                'TickLabelInterpreter', 'none');
            ylabel(ax, 'Condition / trial', 'Interpreter', 'none');
        else
            set(ax, 'YTick', [], 'YTickLabel', {});
            ylabel(ax, '');
        end

        cleanAxisLocal(ax);
    end

    title(t, figTitle, 'Interpreter', 'none');

    if ~isempty(lastAx) && isgraphics(lastAx)
        cb = colorbar(lastAx);
        cb.Label.String = sprintf('%s response', fieldName);
        cb.Label.Interpreter = 'none';
    end
end

function [M, blockInfo] = buildNeuronHeatmapMatrixLocal( ...
    dataBlocks, fieldName, neuronRow, gapRows)

    M = [];
    blockInfo = struct('conditionId', {}, 'label', {}, 'startRow', {}, ...
        'endRow', {}, 'centerRow', {}, 'nTrials', {});

    Tref = [];

    nonEmptyBlocks = find(arrayfun(@(b) ~isempty(b.seqEst), dataBlocks));
    if isempty(nonEmptyBlocks)
        error('No non-empty condition blocks found.');
    end

    for ii = 1:numel(nonEmptyBlocks)
        b = nonEmptyBlocks(ii);
        seqEst = dataBlocks(b).seqEst;

        checkSeqFieldLocal(seqEst, fieldName, neuronRow);

        nTrials = numel(seqEst);
        T = size(seqEst(1).(fieldName), 2);

        if isempty(Tref)
            Tref = T;
        elseif T ~= Tref
            error('Time length mismatch while building heatmap for field %s.', fieldName);
        end

        B = nan(nTrials, Tref);

        for tr = 1:nTrials
            Y = double(seqEst(tr).(fieldName));

            if size(Y, 2) ~= Tref
                error('Time length mismatch in seqEst(%d).%s.', tr, fieldName);
            end

            B(tr, :) = Y(neuronRow, :);
        end

        startRow = size(M, 1) + 1;
        M = [M; B]; %#ok<AGROW>
        endRow = size(M, 1);

        blockInfo(end+1).conditionId = dataBlocks(b).conditionId; %#ok<AGROW>
        blockInfo(end).label = dataBlocks(b).label;
        blockInfo(end).startRow = startRow;
        blockInfo(end).endRow = endRow;
        blockInfo(end).centerRow = (startRow + endRow) / 2;
        blockInfo(end).nTrials = nTrials;

        if ii < numel(nonEmptyBlocks) && gapRows > 0
            M = [M; nan(gapRows, Tref)]; %#ok<AGROW>
        end
    end

    if isempty(M)
        error('No trials found for field %s.', fieldName);
    end
end

function drawConditionGuidesLocal(ax, blockInfo, xEdges, drawBoxes, drawSeparators)
    xLeft = xEdges(1);
    xRight = xEdges(end);
    xWidth = xRight - xLeft;

    for b = 1:numel(blockInfo)
        yTop = blockInfo(b).startRow - 0.5;
        yBottom = blockInfo(b).endRow + 0.5;

        if drawBoxes
            rectangle(ax, ...
                'Position', [xLeft, yTop, xWidth, yBottom - yTop], ...
                'EdgeColor', [0.25 0.25 0.25], ...
                'LineWidth', 0.5, ...
                'HandleVisibility', 'off');
        end

        if drawSeparators && b < numel(blockInfo)
            plot(ax, [xLeft, xRight], [yBottom, yBottom], '-', ...
                'Color', [0.25 0.25 0.25], ...
                'LineWidth', 0.5, ...
                'HandleVisibility', 'off');
        end
    end
end

function [xValues, xLabelText] = getTimeAxisLocal(runMeta, T)
    xValues = 1:T;
    xLabelText = 'Time bin';

    if isfield(runMeta, 'bin_centers') && ...
            isnumeric(runMeta.bin_centers) && ...
            numel(runMeta.bin_centers) == T

        xValues = double(runMeta.bin_centers(:)');
        xLabelText = 'Time';
    end
end

function xEdges = estimateXEdgesLocal(xValues)
    xValues = double(xValues(:)');

    if numel(xValues) == 1
        xEdges = [xValues(1) - 0.5, xValues(1) + 0.5];
        return;
    end

    dx = diff(xValues);
    dx = dx(isfinite(dx) & dx ~= 0);

    if isempty(dx)
        d = 1;
    else
        d = median(abs(dx));
    end

    xEdges = [xValues(1) - d/2, xValues(end) + d/2];
end



function requireFieldLocal(S, fieldName, sourceName)
    if ~isfield(S, fieldName)
        error('%s missing field %s.', sourceName, fieldName);
    end
end

function condition_labels = getConditionLabelsLocal(condition_full, condition_list)
    % Build short condition labels such as G-S-L, G-L-H, P-S-L, P-L-H.
    % No duplicate-label handling.

    if isempty(condition_full)
        error('condition_full is empty.');
    end

    nAll = numel(condition_full);

    stimNameAll = cell(nAll, 1);
    sizeAll = nan(nAll, 1);
    contrastAll = nan(nAll, 1);
    effDirAll = nan(nAll, 1);

    for k = 1:nAll
        requireFieldLocal(condition_full(k), 'stim_name', sprintf('condition_full(%d)', k));
        requireFieldLocal(condition_full(k), 'size', sprintf('condition_full(%d)', k));
        requireFieldLocal(condition_full(k), 'contrast', sprintf('condition_full(%d)', k));

        stimNameAll{k} = lower(char(condition_full(k).stim_name));
        sizeAll(k) = condition_full(k).size;
        contrastAll(k) = condition_full(k).contrast;

        if strcmpi(stimNameAll{k}, 'plaid')
            requireFieldLocal(condition_full(k), 'plaid_dir', sprintf('condition_full(%d)', k));
            effDirAll(k) = condition_full(k).plaid_dir;
        elseif strcmpi(stimNameAll{k}, 'grating')
            requireFieldLocal(condition_full(k), 'grating_dir', sprintf('condition_full(%d)', k));
            effDirAll(k) = condition_full(k).grating_dir;
        else
            error('Unsupported stim_name in condition_full(%d): %s', ...
                k, stimNameAll{k});
        end
    end

    if any(strcmpi(stimNameAll, 'grating')) && any(strcmpi(stimNameAll, 'plaid'))
        stimLabels = {'grating', 'plaid'};
    else
        stimLabels = unique(stimNameAll, 'stable');
        if numel(stimLabels) ~= 2
            error('Expected exactly 2 stim levels in condition_full.');
        end
    end

    sizeVals = unique(sizeAll(isfinite(sizeAll)));
    sizeVals = sort(sizeVals);
    sizeVals = sizeVals(:)';

    if numel(sizeVals) ~= 2
        error('Expected exactly 2 size levels in condition_full.');
    end

    contrastValuesByStim = struct();

    for s = 1:numel(stimLabels)
        idx = strcmpi(stimNameAll, stimLabels{s});

        cvals = unique(contrastAll(idx));
        cvals = cvals(isfinite(cvals));
        cvals = sort(cvals);
        cvals = cvals(:)';

        if numel(cvals) ~= 2
            error('Stim %s does not have exactly 2 contrast levels.', stimLabels{s});
        end

        contrastValuesByStim.(stimLabels{s}) = cvals;
    end

    dirVals = unique(effDirAll(isfinite(effDirAll)));
    dirVals = sort(dirVals);
    dirVals = dirVals(:)';

    if numel(dirVals) ~= 2
        warning('Expected 2 effective direction values, but found %d.', numel(dirVals));
    end

    condShortLabels = { ...
        'G-S-L', 'G-S-H', 'G-L-L', 'G-L-H', ...
        'P-S-L', 'P-S-H', 'P-L-L', 'P-L-H'};

    condition_labels = cell(1, numel(condition_list));

    for ii = 1:numel(condition_list)
        condID = condition_list(ii);

        if condID < 1 || condID > nAll
            error('Condition ID %d is outside condition_full range.', condID);
        end

        currStim = lower(char(condition_full(condID).stim_name));
        currSize = condition_full(condID).size;
        currContrast = condition_full(condID).contrast;

        stimCode = find(strcmpi(stimLabels, currStim), 1);
        sizeCode = find(sizeVals == currSize, 1);

        currContrastLevels = contrastValuesByStim.(currStim);
        contrastCode = find(currContrastLevels == currContrast, 1);

        if isempty(stimCode) || isempty(sizeCode) || isempty(contrastCode)
            error('Could not map condition ID %d to short label.', condID);
        end

        panelCondIndex = (stimCode - 1) * 4 + (sizeCode - 1) * 2 + contrastCode;
        condition_labels{ii} = condShortLabels{panelCondIndex};
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

function rows = getGroupRowsLocal(yDims, groupIdx)
    yDims = reshape(yDims, 1, []);

    if groupIdx < 1 || groupIdx > numel(yDims)
        error('Invalid groupIdx %d.', groupIdx);
    end

    starts = cumsum([1, yDims(1:end-1)]);
    ends = cumsum(yDims);

    rows = starts(groupIdx):ends(groupIdx);
end

function S = loadMatFileFlexibleLocal(fileBase, varargin)
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

function filePath = findOneFileLocal(folderPath, pattern, required)
    if nargin < 3
        required = true;
    end

    if ~isfolder(folderPath)
        if required
            error('Folder not found: %s', folderPath);
        else
            filePath = '';
            return;
        end
    end

    d = dir(fullfile(folderPath, pattern));
    d = d(~[d.isdir]);

    if isempty(d)
        if required
            error('No file matching %s found in %s.', pattern, folderPath);
        else
            filePath = '';
            return;
        end
    end

    if numel(d) > 1
        names = {d.name};
        error('Multiple files matching %s found in %s: %s', ...
            pattern, folderPath, strjoin(names, ', '));
    end

    filePath = fullfile(folderPath, d(1).name);
end

function fileBase = sanitizeFileNameLocal(fileBase)
    badChars = {'/', '\', ':', '*', '?', '"', '<', '>', '|', ' '};

    for i = 1:numel(badChars)
        fileBase = strrep(fileBase, badChars{i}, '_');
    end
end

function saveFigLocal(fig, figFile)
    if exist('savefig', 'file') == 2
        savefig(fig, figFile);
    else
        saveas(fig, figFile);
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

function applyColormapLocal(ax, colormap_name)
    try
        cmap = feval(colormap_name, 256);
        colormap(ax, cmap);
    catch
        warning('Could not apply colormap %s. Using jet.', colormap_name);
        colormap(ax, jet(256));
    end
end

function cleanAxisLocal(ax)
    grid(ax, 'off');
    box(ax, 'off');

    set(ax, ...
        'TickDir', 'out', ...
        'LineWidth', 1, ...
        'FontSize', 10);
end