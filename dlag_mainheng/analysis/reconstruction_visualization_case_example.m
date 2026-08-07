%% reconstruction_visualization_case_example.m
% Plot a fixed-neuron reconstruction case example for one session.
%
% The loading, trial sorting, condition grouping, condition labels, and
% local/global neuron indexing follow reconstruction_visualization.m.
%
% One figure is created for each neural group. Each figure contains:
%   rows    = three user-specified local neuron IDs
%   columns = Original, Within-area, Across-area, Feedforward, Feedback
%
% The 15 panels within one group share one color scale. Color scales are not
% shared across groups. Condition labels and the "Condition / Trial" axis
% label appear only in the upper-left panel. The time-axis scale appears
% only in the third-row, third-column panel.

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

% Used to map trialId back to condition and to obtain short labels such as
% G-S-L, G-L-H, P-S-L, etc.
dat_file = fullfile('.', 'model_data_allruns');
stim_tag = '_2[Gpl2_2c_2sz_400_2_200isi]';

% Three LOCAL neuron IDs for each group. Replace these example values with
% the neuron IDs to be shown. The cell order is group 1, group 2, ... .
case_neuron_ids = { ...
    [1 65 155], ...  % Group 1
    [137 160 218]  ...  % Group 2
    };

% Fixed column order for the case-example figure.
analysis_fields = { ...
    'y', ...
    'yRecon_use_within', ...
    'yRecon_use_across', ...
    'yRecon_use_feedforward', ...
    'yRecon_use_feedback'};

column_labels = { ...
    'Original', ...
    'Within-area', ...
    'Across-area', ...
    'Feedforward', ...
    'Feedback'};

% Heatmap color scale. All 15 panels in a group share the percentile-based
% limits below. Each group is scaled independently.
color_percentiles = [1 99];

% Heatmap and condition-guide options.
condition_gap_rows = 1;
draw_condition_boxes = true;
draw_condition_separators = true;
colormap_name = 'parula';

% Figure style.
figure_visible = 'on';
figure_width = 1800;
figure_height = 760;
font_name = 'Arial';
axis_font_size = 9;
column_font_size = 11;
row_font_size = 10;

% Show the time-axis label, ticks, and outward tick marks in one panel only.
% With the required 3-by-5 layout, [3 3] is the bottom Across-area panel.
time_axis_panel = [3 3];  % [row, column]

% Layout controls. Neuron labels are figure-level text boxes so that all
% three labels share exactly the same right edge, regardless of the
% condition tick labels shown in the upper-left panel.
layout_outer_position = [0.075 0.045 0.915 0.945];
neuron_label_right = 0.062;
neuron_label_width = 0.058;

% Save switches. FIG and SVG are enabled for later vector-graphics editing.
save_fig = true;
save_svg = true;
save_png = false;
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
column_labels = normalizeLabelListLocal(column_labels);

if numel(analysis_fields) ~= numel(column_labels)
    error('analysis_fields and column_labels must have the same length.');
end

if isempty(data_condition)
    use_condition_mode = false;
    condition_list = [];
    modeTag = 'all-condition-model';
else
    use_condition_mode = true;
    condition_list = reshape(data_condition, 1, []);
    modeTag = 'condition-specific-models';
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
        warnIfConditionMismatchLocal( ...
            seqThis, this_condition, condition_index_per_trial_full);

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
case_neuron_ids = validateCaseNeuronIdsLocal( ...
    case_neuron_ids, yDims, numGroups);

fprintf('\nFields and columns:\n');
for f = 1:numel(analysis_fields)
    fprintf('  %-28s -> %s\n', analysis_fields{f}, column_labels{f});
end
fprintf('Output folder: %s\n', saveDir);

%% ------------------------------------------------------------------------
% Plot one 3-by-5 case-example figure per group
% -------------------------------------------------------------------------

for groupIdx = 1:numGroups
    selectedLocalNeuronIds = case_neuron_ids{groupIdx};
    groupRows = getGroupRowsLocal(yDims, groupIdx);
    selectedGlobalRows = groupRows(selectedLocalNeuronIds);
    nNeurons = numel(selectedLocalNeuronIds);

    fprintf('\n============================================================\n');
    fprintf('Group %d: local neuron IDs %s\n', ...
        groupIdx, mat2str(selectedLocalNeuronIds));

    allVals = collectValuesForColorLimitLocal( ...
        dataBlocks, analysis_fields, selectedGlobalRows);
    climVals = robustColorLimitsLocal(allVals, color_percentiles);

    fprintf('Group %d shared color limit across all panels: [%g, %g]\n', ...
        groupIdx, climVals(1), climVals(2));

    figureName = sprintf('%s_case_example_%s_group%d_%dneuron', ...
        data_content, modeTag, groupIdx, nNeurons);

    fig = plotCaseExampleGroupLocal( ...
        dataBlocks, analysis_fields, column_labels, ...
        selectedGlobalRows, selectedLocalNeuronIds, figureName, this_run, ...
        climVals, condition_gap_rows, draw_condition_boxes, ...
        draw_condition_separators, colormap_name, figure_visible, ...
        figure_width, figure_height, font_name, axis_font_size, ...
        column_font_size, row_font_size, time_axis_panel, ...
        layout_outer_position, neuron_label_right, neuron_label_width);

    fileBase = sanitizeFileNameLocal(figureName);

    if save_fig
        figFile = fullfile(saveDir, [fileBase, '.fig']);
        saveFigLocal(fig, figFile);
        fprintf('Saved FIG: %s\n', figFile);
    end

    if save_svg
        svgFile = fullfile(saveDir, [fileBase, '.svg']);
        saveSvgLocal(fig, svgFile);
        fprintf('Saved SVG: %s\n', svgFile);
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
        error('analysis_fields must not contain duplicates.');
    end
end

function labels = normalizeLabelListLocal(labels)
    if ischar(labels)
        labels = {labels};
    elseif isstring(labels)
        labels = cellstr(labels(:));
    elseif iscell(labels)
        labels = labels(:);
        for i = 1:numel(labels)
            if isstring(labels{i})
                labels{i} = char(labels{i});
            end
            if ~ischar(labels{i})
                error('column_labels{%d} must be a char or string.', i);
            end
        end
    else
        error('column_labels must be a char, string array, or cell array.');
    end

    labels = reshape(labels, 1, []);
end

function ids = validateCaseNeuronIdsLocal(ids, yDims, numGroups)
    if ~iscell(ids)
        error('case_neuron_ids must be a cell array with one entry per group.');
    end
    if numel(ids) ~= numGroups
        error(['case_neuron_ids contains %d group entries, but the model has ', ...
            '%d groups.'], numel(ids), numGroups);
    end

    ids = reshape(ids, 1, []);
    for g = 1:numGroups
        thisIds = ids{g};
        if ~isnumeric(thisIds) || ~isvector(thisIds)
            error('case_neuron_ids{%d} must be a numeric vector.', g);
        end

        thisIds = double(reshape(thisIds, 1, []));
        if numel(thisIds) ~= 3
            error(['case_neuron_ids{%d} must contain exactly three local ', ...
                'neuron IDs.'], g);
        end
        if any(~isfinite(thisIds)) || any(thisIds ~= round(thisIds)) || ...
                any(thisIds < 1)
            error('case_neuron_ids{%d} must contain positive integer IDs.', g);
        end
        if numel(unique(thisIds)) ~= numel(thisIds)
            error('case_neuron_ids{%d} contains duplicate neuron IDs.', g);
        end
        if any(thisIds > yDims(g))
            error(['case_neuron_ids{%d} requests neuron %d, but group %d ', ...
                'contains only %d neurons.'], ...
                g, max(thisIds), g, yDims(g));
        end

        ids{g} = thisIds;
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

    yDims = getYDimsFromModelLocal( ...
        Sbest.bestModel, Sbest.res, seqEst, bestFile);

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

function warnIfConditionMismatchLocal( ...
    seqEst, expectedCondition, condition_index_per_trial_full)

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
            'different condition IDs according to trialId metadata.'], ...
            expectedCondition);
    end
end

function condition_order = getConditionOrderFromSeqLocal( ...
    seqEst, condition_index_per_trial_full)

    if ~isfield(seqEst, 'trialId')
        error(['seqEst is missing trialId. Cannot sort all-condition trials ', ...
            'by condition. model_data_prepar.m should have created trialId ', ...
            'for each trial.']);
    end

    trialIds = arrayfun(@(s) s.trialId, seqEst);
    if any(trialIds < 1) || ...
            any(trialIds > numel(condition_index_per_trial_full))
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

function vals = collectValuesForColorLimitLocal( ...
    dataBlocks, analysis_fields, selectedRows)

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
            error(['seqEst(%d).%s has only %d rows, but selected row %d ', ...
                'is requested.'], ...
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
        error(['color_percentiles must be [low high], with ', ...
            '0 <= low < high <= 100.']);
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

function fig = plotCaseExampleGroupLocal( ...
    dataBlocks, analysisFields, columnLabels, ...
    selectedGlobalRows, selectedLocalNeuronIds, figureName, runMeta, ...
    climVals, conditionGapRows, drawConditionBoxes, ...
    drawConditionSeparators, colormapName, figureVisible, ...
    figureWidth, figureHeight, fontName, axisFontSize, ...
    columnFontSize, rowFontSize, timeAxisPanel, ...
    layoutOuterPosition, neuronLabelRight, neuronLabelWidth)

    nNeurons = numel(selectedGlobalRows);
    nFields = numel(analysisFields);

    validateattributes(timeAxisPanel, {'numeric'}, ...
        {'vector', 'numel', 2, 'integer', 'positive', 'finite'}, ...
        mfilename, 'time_axis_panel');
    timeAxisPanel = reshape(timeAxisPanel, 1, []);
    if timeAxisPanel(1) > nNeurons || timeAxisPanel(2) > nFields
        error(['time_axis_panel = [%d %d] is outside the %d-by-%d ', ...
            'panel layout.'], timeAxisPanel(1), timeAxisPanel(2), ...
            nNeurons, nFields);
    end

    validateattributes(layoutOuterPosition, {'numeric'}, ...
        {'vector', 'numel', 4, 'finite'}, ...
        mfilename, 'layout_outer_position');
    validateattributes(neuronLabelRight, {'numeric'}, ...
        {'scalar', 'finite', '>', 0, '<', 1}, ...
        mfilename, 'neuron_label_right');
    validateattributes(neuronLabelWidth, {'numeric'}, ...
        {'scalar', 'finite', 'positive', '<', neuronLabelRight}, ...
        mfilename, 'neuron_label_width');

    fig = figure( ...
        'Name', figureName, ...
        'NumberTitle', 'off', ...
        'Color', 'w', ...
        'Visible', figureVisible, ...
        'Position', [100, 100, figureWidth, figureHeight]);

    set(fig, ...
        'DefaultAxesFontName', fontName, ...
        'DefaultTextFontName', fontName);

    t = tiledlayout(fig, nNeurons, nFields, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');
    try
        t.TileSpacing = 'tight';
    catch
        % Older MATLAB releases use the compact spacing set above.
    end
    try
        t.OuterPosition = layoutOuterPosition;
    catch
        warning(['Could not apply layout_outer_position; using the ', ...
            'default tiled-layout position.']);
    end

    lastAx = [];
    axesGrid = gobjects(nNeurons, nFields);

    for n = 1:nNeurons
        for f = 1:nFields
            fieldName = analysisFields{f};
            ax = nexttile(t, (n - 1) * nFields + f);
            axesGrid(n, f) = ax;
            lastAx = ax;

            [M, blockInfo] = buildNeuronHeatmapMatrixLocal( ...
                dataBlocks, fieldName, selectedGlobalRows(n), ...
                conditionGapRows);

            [xValues, xLabelText] = getTimeAxisLocal(runMeta, size(M, 2));
            xEdges = estimateXEdgesLocal(xValues);

            h = imagesc(ax, xValues, 1:size(M, 1), M);
            set(h, 'AlphaData', isfinite(M));

            set(ax, ...
                'CLim', climVals, ...
                'Color', [1 1 1], ...
                'YDir', 'reverse');

            applyColormapLocal(ax, colormapName);
            hold(ax, 'on');

            drawConditionGuidesLocal( ...
                ax, blockInfo, xEdges, drawConditionBoxes, ...
                drawConditionSeparators);

            xlim(ax, [xEdges(1), xEdges(end)]);
            ylim(ax, [0.5, size(M, 1) + 0.5]);

            if n == 1
                title(ax, columnLabels{f}, ...
                    'Interpreter', 'none', ...
                    'FontName', fontName, ...
                    'FontSize', columnFontSize, ...
                    'FontWeight', 'normal');
            end

            showTimeAxis = ...
                n == timeAxisPanel(1) && f == timeAxisPanel(2);

            if showTimeAxis
                xlabel(ax, xLabelText, ...
                    'Interpreter', 'none', ...
                    'FontName', fontName, ...
                    'FontSize', axisFontSize);
                set(ax, 'XColor', [0 0 0]);
            else
                xlabel(ax, '');
                set(ax, ...
                    'XTick', [], ...
                    'XTickLabel', {}, ...
                    'XColor', 'none');
            end

            % Show condition labels and their axis meaning in the
            % upper-left panel only.
            if n == 1 && f == 1
                set(ax, ...
                    'YTick', [blockInfo.centerRow], ...
                    'YTickLabel', {blockInfo.label}, ...
                    'YColor', [0 0 0], ...
                    'TickLabelInterpreter', 'none');
                ylabel(ax, 'Condition / Trial', ...
                    'Interpreter', 'none', ...
                    'FontName', fontName, ...
                    'FontSize', axisFontSize, ...
                    'FontWeight', 'normal');
            else
                ylabel(ax, '');
                set(ax, ...
                    'YTick', [], ...
                    'YTickLabel', {}, ...
                    'YColor', 'none');
            end

            cleanAxisLocal(ax, fontName, axisFontSize);
        end
    end

    % No tiled-layout title: the requested figure has column headings only.
    if ~isempty(lastAx) && isgraphics(lastAx)
        cb = colorbar(lastAx);
        try
            cb.Layout.Tile = 'east';
        catch
            % Older MATLAB releases keep the colorbar beside lastAx.
        end
        cb.Label.String = 'Response';
        cb.Label.Interpreter = 'none';
        cb.FontName = fontName;
        cb.FontSize = axisFontSize;
        cb.TickDirection = 'out';
    end

    % Add neuron IDs only after the tiled layout and colorbar have reached
    % their final positions. Figure-level annotations keep all IDs aligned.
    drawnow;
    addAlignedNeuronLabelsLocal( ...
        fig, axesGrid(:, 1), selectedLocalNeuronIds, ...
        neuronLabelRight, neuronLabelWidth, fontName, rowFontSize);
end

function addAlignedNeuronLabelsLocal( ...
    fig, firstColumnAxes, neuronIds, labelRight, labelWidth, ...
    fontName, fontSize)

    labelLeft = labelRight - labelWidth;

    for n = 1:numel(firstColumnAxes)
        ax = firstColumnAxes(n);
        oldUnits = ax.Units;
        ax.Units = 'normalized';
        axPos = ax.Position;
        ax.Units = oldUnits;

        labelHeight = min(0.05, 0.6 * axPos(4));
        labelBottom = axPos(2) + 0.5 * axPos(4) - 0.5 * labelHeight;

        annotation(fig, 'textbox', ...
            [labelLeft, labelBottom, labelWidth, labelHeight], ...
            'String', sprintf('Neuron %d', neuronIds(n)), ...
            'Interpreter', 'none', ...
            'FontName', fontName, ...
            'FontSize', fontSize, ...
            'FontWeight', 'normal', ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', ...
            'EdgeColor', 'none', ...
            'Margin', 0, ...
            'FitBoxToText', 'off');
    end
end

function [M, blockInfo] = buildNeuronHeatmapMatrixLocal( ...
    dataBlocks, fieldName, neuronRow, gapRows)

    M = [];
    blockInfo = struct( ...
        'conditionId', {}, ...
        'label', {}, ...
        'startRow', {}, ...
        'endRow', {}, ...
        'centerRow', {}, ...
        'nTrials', {});

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
            error('Time length mismatch while building heatmap for field %s.', ...
                fieldName);
        end

        B = nan(nTrials, Tref);
        for tr = 1:nTrials
            Y = double(seqEst(tr).(fieldName));
            if size(Y, 2) ~= Tref
                error('Time length mismatch in seqEst(%d).%s.', ...
                    tr, fieldName);
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

function drawConditionGuidesLocal( ...
    ax, blockInfo, xEdges, drawBoxes, drawSeparators)

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

function condition_labels = getConditionLabelsLocal( ...
    condition_full, condition_list)
    % Build short labels such as G-S-L, G-L-H, P-S-L, and P-L-H.
    % This follows reconstruction_visualization.m and intentionally performs
    % no duplicate-label handling.

    if isempty(condition_full)
        error('condition_full is empty.');
    end

    nAll = numel(condition_full);
    stimNameAll = cell(nAll, 1);
    sizeAll = nan(nAll, 1);
    contrastAll = nan(nAll, 1);
    effDirAll = nan(nAll, 1);

    for k = 1:nAll
        requireFieldLocal( ...
            condition_full(k), 'stim_name', sprintf('condition_full(%d)', k));
        requireFieldLocal( ...
            condition_full(k), 'size', sprintf('condition_full(%d)', k));
        requireFieldLocal( ...
            condition_full(k), 'contrast', sprintf('condition_full(%d)', k));

        stimNameAll{k} = lower(char(condition_full(k).stim_name));
        sizeAll(k) = condition_full(k).size;
        contrastAll(k) = condition_full(k).contrast;

        if strcmpi(stimNameAll{k}, 'plaid')
            requireFieldLocal( ...
                condition_full(k), 'plaid_dir', ...
                sprintf('condition_full(%d)', k));
            effDirAll(k) = condition_full(k).plaid_dir;
        elseif strcmpi(stimNameAll{k}, 'grating')
            requireFieldLocal( ...
                condition_full(k), 'grating_dir', ...
                sprintf('condition_full(%d)', k));
            effDirAll(k) = condition_full(k).grating_dir;
        else
            error('Unsupported stim_name in condition_full(%d): %s', ...
                k, stimNameAll{k});
        end
    end

    if any(strcmpi(stimNameAll, 'grating')) && ...
            any(strcmpi(stimNameAll, 'plaid'))
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
            error('Stim %s does not have exactly 2 contrast levels.', ...
                stimLabels{s});
        end
        contrastValuesByStim.(stimLabels{s}) = cvals;
    end

    dirVals = unique(effDirAll(isfinite(effDirAll)));
    dirVals = sort(dirVals);
    dirVals = dirVals(:)';

    if numel(dirVals) ~= 2
        warning('Expected 2 effective direction values, but found %d.', ...
            numel(dirVals));
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

        panelCondIndex = ...
            (stimCode - 1) * 4 + (sizeCode - 1) * 2 + contrastCode;
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

function saveSvgLocal(fig, svgFile)
    drawnow;

    try
        if exist('exportgraphics', 'file') == 2
            exportgraphics(fig, svgFile, ...
                'ContentType', 'vector', ...
                'BackgroundColor', 'white');
        else
            print(fig, svgFile, '-dsvg');
        end
    catch
        warning('SVG export failed. Falling back to saveas.');
        saveas(fig, svgFile);
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
        warning('PNG export failed. Falling back to saveas.');
        saveas(fig, pngFile);
    end
end

function applyColormapLocal(ax, colormapName)
    try
        cmap = feval(colormapName, 256);
        colormap(ax, cmap);
    catch
        warning('Could not apply colormap %s. Using jet.', colormapName);
        colormap(ax, jet(256));
    end
end

function cleanAxisLocal(ax, fontName, fontSize)
    grid(ax, 'off');
    box(ax, 'off');

    set(ax, ...
        'TickDir', 'out', ...
        'LineWidth', 1, ...
        'FontName', fontName, ...
        'FontSize', fontSize);
end
