%% plot_timescale_compare.m
% Plot DLAG latent timescale distributions by latent category.

clc;
clear;
close all;

%% ---------------- User parameters --------------------------------------

data_content = 'demean_count_within_trial';
% options:
% raw_count, raw_fr, z_within_trial, z_within_condition,
% z_across_conditions, demean_count_within_trial, demean_fr_within_trial,
% demean_pooledsd_within_condition

data_condition = [];
% []   : all-condition model
% 1:16 : condition-specific models

runIdx = 1;

stim_tag = '_2[Gpl2_2c_2sz_400_2_200isi]';

timescale_source = 'model-timescale';
% options:
% model-timescale
% posterior-timescale

use_dsl_filter = false;

dsl_field = 'logical';
% condition-specific models:
% rawlogical, logical
%
% all-condition model:
% rawlogical, logical, logical_bystimdir, logical_bystimnamedir,
% logical_bycondition

make_plots = true;
save_figures = true;
save_mat = true;

use_log_y = false;
jitter_width = 0.16;
marker_size = 28;
jitter_seed = 1;

save_each_condition_extraction_mat = false;

rng(jitter_seed);

%% ---------------- Setup -------------------------------------------------

if isempty(data_condition)
    use_condition_mode = false;
    condition_list = [];
    numConditions = 1;
    model_mode = 'all_condition_model';
else
    use_condition_mode = true;
    condition_list = data_condition(:)';
    numConditions = numel(condition_list);
    model_mode = 'condition_specific_models';
end

validateDslFieldForModeLocal(use_dsl_filter, dsl_field, use_condition_mode);

source_tag = makeFileTagLocal(timescale_source);
selection_tag = makeSelectionTagLocal(use_dsl_filter, dsl_field);

if use_condition_mode
    dat_file = '.\model_data_allruns';
    fprintf('Reading from %s\n', dat_file);
    load(dat_file, 'model_data_allruns');

    all_run_tags = get_all_run_tags(model_data_allruns);
    run_idx = find(strcmp(all_run_tags, stim_tag));

    if isempty(run_idx)
        error('Requested stim_tag not found: %s', stim_tag);
    end

    if numel(run_idx) > 1
        error('Duplicate stim_tag found: %s', stim_tag);
    end

    condition_full = model_data_allruns{run_idx}.conditions_full;
else
    condition_full = [];
end

AllConditionResults = struct([]);

%% ---------------- Main loop --------------------------------------------

for cond_i = 1:numConditions

    if use_condition_mode
        this_condition = condition_list(cond_i);
        baseDir = ['./FA_Dlag_', data_content, '_condition', num2str(this_condition)];
    else
        this_condition = [];
        baseDir = ['./FA_Dlag_', data_content];
    end

    tempfname = sprintf('%s/mat_results/run%03d', baseDir, runIdx);
    fprintf('\nReading model from %s\n', tempfname);

    bestFile = findFirstFileLocal(tempfname, 'bestmodel*');
    Sbest = load(bestFile, "bestModel", "gp_params");

    if ~isfield(Sbest, 'bestModel')
        error('bestmodel file is missing variable bestModel: %s', bestFile);
    end

    if ~isfield(Sbest, 'gp_params')
        error('bestmodel file is missing variable gp_params: %s', bestFile);
    end

    bootFile = findFirstFileLocal(tempfname, 'bootstrapResults*');
    Sboot = load(bootFile, "ambiguousIdxs");

    if ~isfield(Sboot, 'ambiguousIdxs')
        error('bootstrapResults file is missing variable ambiguousIdxs: %s', bootFile);
    end

    DSL = [];
    dslFile = '';

    if use_dsl_filter
        dslFile = fullfile(tempfname, 'DSL_and_latent_category_stats.mat');

        if ~exist(dslFile, 'file')
            error('Saved DSL file was not found: %s', dslFile);
        end

        Sdsl = load(dslFile, 'DSL');

        if ~isfield(Sdsl, 'DSL')
            error('Saved DSL file does not contain variable DSL: %s', dslFile);
        end

        DSL = Sdsl.DSL;
    end

    TimescaleStats = collectTimescaleStatsLocal( ...
        Sbest.bestModel, ...
        Sbest.gp_params, ...
        Sboot.ambiguousIdxs, ...
        DSL, ...
        use_dsl_filter, ...
        dsl_field, ...
        timescale_source);

    TimescaleStats.meta.data_content = data_content;
    TimescaleStats.meta.data_condition = this_condition;
    TimescaleStats.meta.use_condition_mode = use_condition_mode;
    TimescaleStats.meta.model_mode = model_mode;
    TimescaleStats.meta.condition_index_in_list = cond_i;
    TimescaleStats.meta.condition_list = condition_list;
    TimescaleStats.meta.runIdx = runIdx;
    TimescaleStats.meta.baseDir = baseDir;
    TimescaleStats.meta.tempfname = tempfname;
    TimescaleStats.meta.bestmodel_file = bestFile;
    TimescaleStats.meta.bootstrap_file = bootFile;
    TimescaleStats.meta.dsl_file = dslFile;
    TimescaleStats.meta.selection_tag = selection_tag;
    TimescaleStats.meta.source_tag = source_tag;

    if ~use_condition_mode

        if make_plots
            figHandles = plotTimescaleStatsLocal( ...
                TimescaleStats, ...
                data_content, ...
                model_mode, ...
                timescale_source, ...
                selection_tag, ...
                use_log_y, ...
                jitter_width, ...
                marker_size);

            if save_figures
                for g = 1:numel(figHandles)
                    figName = sprintf('timescale_distribution_%s_%s_group%d', ...
                        source_tag, selection_tag, g);
                    saveAndCloseFigureLocal(figHandles(g), tempfname, figName);
                end
            end
        end

        if save_mat
            matName = sprintf('timescale_and_latent_category_stats_%s_%s.mat', ...
                source_tag, selection_tag);

            save(fullfile(tempfname, matName), ...
                'TimescaleStats', ...
                'data_content', ...
                'data_condition', ...
                'runIdx', ...
                'timescale_source', ...
                'use_dsl_filter', ...
                'dsl_field', ...
                'selection_tag', ...
                'source_tag', ...
                'model_mode');
        end

    else

        AllConditionResults(cond_i).condition = this_condition;
        AllConditionResults(cond_i).baseDir = baseDir;
        AllConditionResults(cond_i).tempfname = tempfname;
        AllConditionResults(cond_i).TimescaleStats = TimescaleStats;

        if save_mat && save_each_condition_extraction_mat
            matName = sprintf('timescale_and_latent_category_stats_%s_%s.mat', ...
                source_tag, selection_tag);

            save(fullfile(tempfname, matName), ...
                'TimescaleStats', ...
                'data_content', ...
                'data_condition', ...
                'this_condition', ...
                'runIdx', ...
                'timescale_source', ...
                'use_dsl_filter', ...
                'dsl_field', ...
                'selection_tag', ...
                'source_tag', ...
                'model_mode');
        end
    end
end

%% ---------------- Sum all conditions -----------------------------------

if use_condition_mode

    [SummaryTimescale, SummaryFigs] = summarizeAllConditionsTimescaleLocal( ...
        AllConditionResults, ...
        condition_list, ...
        condition_full, ...
        data_content, ...
        model_mode, ...
        timescale_source, ...
        use_dsl_filter, ...
        dsl_field, ...
        selection_tag, ...
        source_tag, ...
        use_log_y, ...
        jitter_width, ...
        marker_size);

    if save_figures
        for g = 1:numel(SummaryFigs)
            figName = sprintf('%s_sum_all_conditions_timescale_distribution_%s_%s_group%d', ...
                data_content, source_tag, selection_tag, g);
            saveAndCloseFigureLocal(SummaryFigs(g), '.', figName);
        end
    end

    if save_mat
        matName = sprintf('%s_sum_all_conditions_timescale_and_latent_category_stats_%s_%s.mat', ...
            data_content, source_tag, selection_tag);

        save(matName, ...
            'AllConditionResults', ...
            'SummaryTimescale', ...
            'condition_list', ...
            'condition_full', ...
            'data_content', ...
            'runIdx', ...
            'stim_tag', ...
            'timescale_source', ...
            'use_dsl_filter', ...
            'dsl_field', ...
            'selection_tag', ...
            'source_tag', ...
            'model_mode');
    end
end

if save_figures
    close all;
end


%% ========================================================================
%% Local functions
%% ========================================================================

function TimescaleStats = collectTimescaleStatsLocal( ...
    bestModel, gp_params, ambiguousIdxs, DSL, use_dsl_filter, ...
    dsl_field, timescale_source)

    if ~isfield(bestModel, 'xDim_across')
        error('bestModel is missing xDim_across.');
    end

    if ~isfield(bestModel, 'xDim_within')
        error('bestModel is missing xDim_within.');
    end

    xDim_across = bestModel.xDim_across;
    xDim_within = bestModel.xDim_within(:)';
    numGroups = numel(xDim_within);
    localDims = xDim_across + xDim_within;

    [tau_across, tau_within] = getTimescaleFromBestmodelLocal( ...
        gp_params, xDim_across, xDim_within, timescale_source);

    latentClass = classifyDlagLatentsLocal(xDim_across, gp_params, ambiguousIdxs);
    labels = latentClass.categoryLabels;
    numCategories = numel(labels);

    TimescaleStats = struct();
    TimescaleStats.labels = labels;

    TimescaleStats.meta.timescale_source = timescale_source;
    TimescaleStats.meta.use_dsl_filter = use_dsl_filter;
    TimescaleStats.meta.dsl_field = dsl_field;
    TimescaleStats.meta.numGroups = numGroups;
    TimescaleStats.meta.xDim_across = xDim_across;
    TimescaleStats.meta.xDim_within = xDim_within;
    TimescaleStats.meta.localDims = localDims;
    TimescaleStats.meta.value_unit = 'ms';

    TimescaleStats.classification = latentClass;

    for g = 1:numGroups

        tau_within_g = reshape(tau_within{g}, 1, []);
        tau_local = [reshape(tau_across, 1, []), tau_within_g];

        if numel(tau_local) ~= localDims(g)
            error('Group %d: tau_local length %d does not match localDim %d.', ...
                g, numel(tau_local), localDims(g));
        end

        if use_dsl_filter
            keepMask = getSavedDslMaskLocal(DSL, dsl_field, g, localDims(g));
        else
            keepMask = true(1, localDims(g));
        end

        categoryMasks = makeDlagCategoryMasksLocal( ...
            localDims(g), ...
            xDim_across, ...
            latentClass.feedforwardIdx, ...
            latentClass.feedbackIdx, ...
            latentClass.ambiguousIdx);

        TimescaleStats.group(g).name = sprintf('Group %d', g);
        TimescaleStats.group(g).tau_local = tau_local;
        TimescaleStats.group(g).tau_across = tau_across;
        TimescaleStats.group(g).tau_within = tau_within_g;
        TimescaleStats.group(g).keepMask = keepMask;
        TimescaleStats.group(g).categoryMasks = categoryMasks;

        TimescaleStats.group(g).values = cell(1, numCategories);
        TimescaleStats.group(g).localIdx = cell(1, numCategories);

        for c = 1:numCategories
            thisMask = keepMask & categoryMasks{c} & isfinite(tau_local);
            idx = find(thisMask);

            TimescaleStats.group(g).values{c} = tau_local(idx);
            TimescaleStats.group(g).localIdx{c} = idx;
        end
    end
end

function [tau_across, tau_within] = getTimescaleFromBestmodelLocal( ...
    gp_params, xDim_across, xDim_within, timescale_source)

    switch lower(timescale_source)

        case 'model-timescale'

            if ~isfield(gp_params, 'tau_across')
                error('gp_params is missing tau_across.');
            end

            if ~isfield(gp_params, 'tau_within')
                error('gp_params is missing tau_within.');
            end

            tau_across = reshape(gp_params.tau_across, 1, []);

            if numel(tau_across) < xDim_across
                error('gp_params.tau_across has fewer entries than xDim_across.');
            end

            tau_across = tau_across(1:xDim_across);

            if ~iscell(gp_params.tau_within)
                error('Expected gp_params.tau_within to be a cell array.');
            end

            tau_within = gp_params.tau_within;

            if numel(tau_within) ~= numel(xDim_within)
                error('Length of gp_params.tau_within does not match xDim_within.');
            end

            for g = 1:numel(xDim_within)
                tau_within{g} = reshape(tau_within{g}, 1, []);

                if numel(tau_within{g}) < xDim_within(g)
                    error('gp_params.tau_within{%d} has fewer entries than xDim_within(%d).', g, g);
                end

                tau_within{g} = tau_within{g}(1:xDim_within(g));
            end

        case 'posterior-timescale'

            error('posterior-timescale is not implemented in this version.');

        otherwise

            error('Unknown timescale_source: %s', timescale_source);
    end
end

function keepMask = getSavedDslMaskLocal(DSL, dsl_field, groupIdx, localDim)

    if isempty(DSL) || ~isstruct(DSL)
        error('DSL must be a non-empty struct when use_dsl_filter = true.');
    end

    if ~isfield(DSL, dsl_field)
        error('DSL does not contain field "%s".', dsl_field);
    end

    dslData = DSL.(dsl_field);

    if ~iscell(dslData)
        error('DSL.%s must be a cell array.', dsl_field);
    end

    if isvector(dslData)
        if numel(dslData) < groupIdx
            error('DSL.%s does not contain group %d.', dsl_field, groupIdx);
        end

        keepMask = reshape(dslData{groupIdx}, 1, []) ~= 0;
    else
        error('DSL.%s is not a group-level mask field.', dsl_field);
    end

    if numel(keepMask) ~= localDim
        error('DSL.%s group %d has wrong length. Expected %d, got %d.', ...
            dsl_field, groupIdx, localDim, numel(keepMask));
    end
end

function latentClass = classifyDlagLatentsLocal(xDim_across, gp_params, ambiguousIdxs)

    if ~isfield(gp_params, 'delays')
        error('gp_params must contain field delays.');
    end

    acrossDelay = reshape(gp_params.delays, 1, []);

    if numel(acrossDelay) < xDim_across
        error('gp_params.delays has fewer entries than xDim_across.');
    end

    acrossDelay = acrossDelay(1:xDim_across);

    ambiguousIdxs = unique(ambiguousIdxs(:)');
    ambiguousIdxs = ambiguousIdxs(ambiguousIdxs >= 1 & ambiguousIdxs <= xDim_across);

    acrossIdx = 1:xDim_across;
    zeroOrNaNIdx = acrossIdx((acrossDelay == 0) | isnan(acrossDelay));
    ambiguousAll = unique([ambiguousIdxs, zeroOrNaNIdx]);

    ffIdx = find(acrossDelay > 0);
    fbIdx = find(acrossDelay < 0);

    ffIdx = setdiff(ffIdx, ambiguousAll);
    fbIdx = setdiff(fbIdx, ambiguousAll);

    coveredAcross = unique([ffIdx, fbIdx, ambiguousAll]);
    missingAcross = setdiff(acrossIdx, coveredAcross);

    if ~isempty(missingAcross)
        ambiguousAll = unique([ambiguousAll, missingAcross]);
    end

    latentClass = struct();
    latentClass.categoryLabels = {'Across', 'Within', 'Feedforward', 'Feedback', 'Ambiguous'};
    latentClass.acrossDelay = acrossDelay;
    latentClass.acrossIdx = acrossIdx;
    latentClass.feedforwardIdx = ffIdx;
    latentClass.feedbackIdx = fbIdx;
    latentClass.ambiguousIdx = ambiguousAll;
end

function categoryMasks = makeDlagCategoryMasksLocal( ...
    localDim, xDim_across, ffIdx, fbIdx, ambiguousIdx)

    acrossMask = false(1, localDim);
    acrossMask(1:xDim_across) = true;

    withinMask = false(1, localDim);
    withinMask(xDim_across+1:end) = true;

    ffMask = false(1, localDim);
    ffIdx = ffIdx(ffIdx >= 1 & ffIdx <= localDim);
    ffMask(ffIdx) = true;

    fbMask = false(1, localDim);
    fbIdx = fbIdx(fbIdx >= 1 & fbIdx <= localDim);
    fbMask(fbIdx) = true;

    ambMask = false(1, localDim);
    ambiguousIdx = ambiguousIdx(ambiguousIdx >= 1 & ambiguousIdx <= localDim);
    ambMask(ambiguousIdx) = true;

    categoryMasks = {acrossMask, withinMask, ffMask, fbMask, ambMask};
end

function figHandles = plotTimescaleStatsLocal( ...
    TimescaleStats, data_content, model_mode, timescale_source, ...
    selection_tag, use_log_y, jitter_width, marker_size)

    numGroups = numel(TimescaleStats.group);
    figHandles = gobjects(1, numGroups);

    for g = 1:numGroups
        figTitle = sprintf('%s | %s | %s | %s | Group %d', ...
            data_content, model_mode, timescale_source, selection_tag, g);

        figHandles(g) = plotOneGroupTimescaleDistributionLocal( ...
            TimescaleStats.group(g).values, ...
            TimescaleStats.labels, ...
            figTitle, ...
            use_log_y, ...
            jitter_width, ...
            marker_size);
    end
end


function [SummaryTimescale, figHandles] = summarizeAllConditionsTimescaleLocal( ...
    AllConditionResults, condition_list, condition_full, data_content, ...
    model_mode, timescale_source, use_dsl_filter, dsl_field, selection_tag, ...
    source_tag, use_log_y, jitter_width, marker_size)

    if isempty(AllConditionResults)
        error('AllConditionResults is empty.');
    end

    if isempty(condition_full)
        error('condition_full is required for sum-all-conditions labels.');
    end

    numGroups = numel(AllConditionResults(1).TimescaleStats.group);
    labels = AllConditionResults(1).TimescaleStats.labels;
    numCategories = numel(labels);

    conditionMap = buildConditionSummaryMapLocal(condition_full, condition_list);

    SummaryTimescale = struct();
    SummaryTimescale.labels = labels;

    SummaryTimescale.meta.data_content = data_content;
    SummaryTimescale.meta.model_mode = model_mode;
    SummaryTimescale.meta.timescale_source = timescale_source;
    SummaryTimescale.meta.source_tag = source_tag;
    SummaryTimescale.meta.use_dsl_filter = use_dsl_filter;
    SummaryTimescale.meta.dsl_field = dsl_field;
    SummaryTimescale.meta.selection_tag = selection_tag;
    SummaryTimescale.meta.conditions = condition_list;
    SummaryTimescale.meta.conditionMap = conditionMap;
    SummaryTimescale.meta.panelConditionLabels = conditionMap.meta.panelCondLabels;
    SummaryTimescale.meta.panelConditionShortLabels = conditionMap.meta.panelCondShortLabels;
    SummaryTimescale.meta.stimDirLabels = conditionMap.meta.stimDirLabels;
    SummaryTimescale.meta.stimDirValues = conditionMap.meta.stimDirValues;
    SummaryTimescale.meta.stimDirPanelTitles = { ...
        sprintf('%s = %s', conditionMap.meta.stimDirLabels{1}, ...
            formatSummaryValueLocal(conditionMap.meta.stimDirValues(1))), ...
        sprintf('%s = %s', conditionMap.meta.stimDirLabels{2}, ...
            formatSummaryValueLocal(conditionMap.meta.stimDirValues(2)))};

    figHandles = gobjects(1, numGroups);

    for g = 1:numGroups

        dirValues = cell(1, 2);
        dirConditionIds = cell(1, 2);

        for d = 1:2
            dirValues{d} = cell(8, numCategories);
            dirConditionIds{d} = nan(8, 1);

            for p = 1:8
                for c = 1:numCategories
                    dirValues{d}{p, c} = [];
                end
            end
        end

        for ci = 1:numel(AllConditionResults)
            condID = AllConditionResults(ci).condition;

            mapIdx = find([conditionMap.entries.conditionId] == condID, 1);

            if isempty(mapIdx)
                warning('Condition %d not found in condition_full mapping. Skipping.', condID);
                continue;
            end

            stimDirCode = conditionMap.entries(mapIdx).stimDirCode;
            panelCondIdx = conditionMap.entries(mapIdx).panelCondIndex;

            dirConditionIds{stimDirCode}(panelCondIdx) = condID;

            for c = 1:numCategories
                vals = AllConditionResults(ci).TimescaleStats.group(g).values{c};
                dirValues{stimDirCode}{panelCondIdx, c} = vals;
            end
        end

        for d = 1:2
            SummaryTimescale.group(g).stim_dir(d).label = ...
                conditionMap.meta.stimDirLabels{d};
            SummaryTimescale.group(g).stim_dir(d).value = ...
                conditionMap.meta.stimDirValues(d);
            SummaryTimescale.group(g).stim_dir(d).conditionIds = ...
                dirConditionIds{d};
            SummaryTimescale.group(g).stim_dir(d).values = ...
                dirValues{d};
        end

        figTitle = sprintf('%s | %s | %s | %s | Group %d', ...
            data_content, model_mode, timescale_source, selection_tag, g);

        figHandles(g) = plotSplitDirGroupedTimescalePointsLocal( ...
            dirValues{1}, ...
            dirValues{2}, ...
            labels, ...
            SummaryTimescale.meta.panelConditionShortLabels, ...
            SummaryTimescale.meta.stimDirPanelTitles, ...
            figTitle, ...
            use_log_y, ...
            jitter_width, ...
            marker_size);
    end
end


function placeCategoryLabelsBelowAxisLocal(ax, xCenters, categoryLabels, use_log_y)

    yl = ylim(ax);

    if use_log_y
        logRange = log10(yl(2)) - log10(yl(1));
        yText = 10 .^ (log10(yl(1)) - 0.075 * logRange);
    else
        yText = yl(1) - 0.075 * (yl(2) - yl(1));
    end

    for c = 1:numel(categoryLabels)
        text(ax, xCenters(c), yText, categoryLabels{c}, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', ...
            'FontWeight', 'bold', ...
            'Interpreter', 'none', ...
            'Clipping', 'off');
    end
end

function vals = flattenValueCellsLocal(valueCells)

    vals = [];

    for k = 1:numel(valueCells)
        thisVals = valueCells{k};
        thisVals = thisVals(isfinite(thisVals));
        vals = [vals, reshape(thisVals, 1, [])];
    end
end

function applyTimescaleAxisLocal(ax, values, use_log_y)

    values = values(isfinite(values));

    if use_log_y
        values = values(values > 0);

        if isempty(values)
            set(ax, 'YScale', 'linear');
            ylim(ax, [0, 1]);
            return;
        end

        set(ax, 'YScale', 'log');

        minVal = min(values);
        maxVal = max(values);

        if minVal == maxVal
            minVal = minVal / 2;
            maxVal = maxVal * 2;
        else
            minVal = minVal / 1.4;
            maxVal = maxVal * 1.4;
        end

        minVal = max(minVal, realmin);
        ylim(ax, [minVal, maxVal]);

    else
        if isempty(values)
            ylim(ax, [0, 1]);
            return;
        end

        minVal = min(values);
        maxVal = max(values);

        if minVal == maxVal
            padding = max(abs(maxVal) * 0.2, 1);
            ylim(ax, [minVal - padding, maxVal + padding]);
        else
            padding = 0.08 * (maxVal - minVal);
            ylim(ax, [minVal - padding, maxVal + padding]);
        end
    end
end

function colors = getConditionColorsLocal(nCond)

    if nargin < 1 || isempty(nCond)
        nCond = 8;
    end

    colors = lines(nCond);
end

function conditionMap = buildConditionSummaryMapLocal(condition_full, condition_list)

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
        effDirAll(k) = getConditionEffectiveDirCanonicalLocal(condition_full(k), k);
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

    dirVals = unique(effDirAll(isfinite(effDirAll)));
    dirVals = sort(dirVals(:)');

    if numel(dirVals) ~= 2
        error('Expected exactly 2 effective canonical direction values in condition_full, found %d: %s.', ...
            numel(dirVals), mat2str(dirVals));
    end

    stimDirLabels = {'stim_dir1', 'stim_dir2'};

    condLabels = { ...
        'grating-small-low', 'grating-small-high', ...
        'grating-large-low', 'grating-large-high', ...
        'plaid-small-low', 'plaid-small-high', ...
        'plaid-large-low', 'plaid-large-high'};

    condShortLabels = { ...
        'G-S-L', 'G-S-H', 'G-L-L', 'G-L-H', ...
        'P-S-L', 'P-S-H', 'P-L-L', 'P-L-H'};

    entries = struct([]);

    for ii = 1:numel(condition_list)
        condID = condition_list(ii);

        if condID < 1 || condID > nAll
            error('Condition ID %d is outside condition_full range.', condID);
        end

        currStim = lower(string(condition_full(condID).stim_name));
        currSize = condition_full(condID).size;
        currContrast = condition_full(condID).contrast;
        currDir = getConditionEffectiveDirCanonicalLocal(condition_full(condID), condID);

        stimCode = find(strcmp(cellstr(stimLabels), char(currStim)), 1);
        sizeCode = find(abs(sizeVals - currSize) < 1e-10, 1);

        currContrastLevels = contrastValuesByStim.(char(currStim));
        contrastCode = find(abs(currContrastLevels - currContrast) < 1e-10, 1);

        stimDirCode = find(abs(dirVals - currDir) < 1e-10, 1);

        if isempty(stimCode)
            error('Could not map stim_name %s to stim code.', char(currStim));
        end

        if isempty(sizeCode)
            error('Could not map size value %s to size code.', ...
                formatSummaryValueLocal(currSize));
        end

        if isempty(contrastCode)
            error('Could not map contrast value %s to contrast code for stim %s.', ...
                formatSummaryValueLocal(currContrast), char(currStim));
        end

        if isempty(stimDirCode)
            error('Could not map canonical direction value %s to stim_dir code.', ...
                formatSummaryValueLocal(currDir));
        end

        panelCondIndex = (stimCode - 1) * 4 + (sizeCode - 1) * 2 + contrastCode;

        entries(ii).conditionId = condID;
        entries(ii).stimName = char(currStim);
        entries(ii).stimCode = stimCode;
        entries(ii).sizeValue = currSize;
        entries(ii).sizeCode = sizeCode;
        entries(ii).sizeLabel = ternaryLabelLocal(sizeCode, 'small', 'large');
        entries(ii).contrastValue = currContrast;
        entries(ii).contrastCode = contrastCode;
        entries(ii).contrastLabel = ternaryLabelLocal(contrastCode, 'low', 'high');
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

function d = getConditionEffectiveDirCanonicalLocal(cond, condID)

    if nargin < 2
        condID = NaN;
    end

    if ~isfield(cond, 'stim_name')
        error('condition_full(%d) missing field stim_name.', condID);
    end

    currStim = lower(string(cond.stim_name));

    if currStim == "grating"
        if ~isfield(cond, 'grating_dir')
            error('condition_full(%d) is grating but missing field grating_dir.', condID);
        end

        d = cond.grating_dir;

    elseif currStim == "plaid"
        if ~isfield(cond, 'plaid_dir')
            error('condition_full(%d) is plaid but missing field plaid_dir.', condID);
        end

        d = cond.plaid_dir;

    else
        error('Unsupported stim_name in condition_full(%d): %s', condID, char(currStim));
    end

    d = canonicalAngle360Local(d);
end

function a = canonicalAngle360Local(a)

    a = double(a);
    finiteMask = isfinite(a);

    a(finiteMask) = mod(a(finiteMask), 360);

    tol = 1e-10;
    nearInteger = finiteMask & abs(a - round(a)) < tol;
    a(nearInteger) = round(a(nearInteger));

    a(finiteMask & abs(a) < tol) = 0;
    a(finiteMask & abs(a - 360) < tol) = 0;
end

function filename = findFirstFileLocal(folderPath, pattern)

    files = dir(fullfile(folderPath, pattern));

    if isempty(files)
        error('No %s file found in %s', pattern, folderPath);
    end

    filename = fullfile(folderPath, files(1).name);
end

function saveAndCloseFigureLocal(figHandle, folderPath, baseName)

    if isempty(figHandle) || ~isgraphics(figHandle)
        return;
    end

    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end

    savefig(figHandle, fullfile(folderPath, [baseName, '.fig']));

    try
        exportgraphics(figHandle, fullfile(folderPath, [baseName, '.png']), ...
            'Resolution', 400);
    catch
        saveas(figHandle, fullfile(folderPath, [baseName, '.png']));
    end

    close(figHandle);
end

function tag = makeFileTagLocal(label)

    label = char(string(label));
    tag = lower(label);
    tag = regexprep(tag, '[^a-z0-9]+', '_');
    tag = regexprep(tag, '^_+|_+$', '');
end

function selection_tag = makeSelectionTagLocal(use_dsl_filter, dsl_field)

    if ~use_dsl_filter
        selection_tag = 'all';
    else
        selection_tag = makeFileTagLocal(dsl_field);
    end
end

function validateDslFieldForModeLocal(use_dsl_filter, dsl_field, use_condition_mode)

    if ~use_dsl_filter
        return;
    end

    if use_condition_mode
        allowedFields = { ...
            'rawlogical', ...
            'logical'};
    else
        allowedFields = { ...
            'rawlogical', ...
            'logical', ...
            'logical_bystimdir', ...
            'logical_bystimnamedir', ...
            'logical_bycondition'};
    end

    if ~ismember(dsl_field, allowedFields)
        if use_condition_mode
            modeText = 'condition-specific models';
        else
            modeText = 'all-condition model';
        end

        error('Invalid dsl_field "%s" for %s.', dsl_field, modeText);
    end
end

function s = formatSummaryValueLocal(v)

    if ~isfinite(v)
        s = 'NaN';
    elseif abs(v - round(v)) < 1e-10
        s = sprintf('%d', round(v));
    else
        s = sprintf('%.4g', v);
    end
end

function out = ternaryLabelLocal(code, label1, label2)

    if isempty(code) || ~isfinite(code)
        out = '';
    elseif code == 1
        out = label1;
    else
        out = label2;
    end
end

function all_tags = get_all_run_tags(model_data_allruns)

    all_tags = cell(numel(model_data_allruns), 1);

    for j = 1:numel(model_data_allruns)
        if ~isfield(model_data_allruns{j}, 'stim_tag')
            error('stim_tag missing in model_data_allruns{%d}.', j);
        end

        all_tags{j} = model_data_allruns{j}.stim_tag;
    end
end



function figHandle = plotSplitDirGroupedTimescalePointsLocal( ...
    valuesDir1, valuesDir2, categoryLabels, conditionShortLabels, ...
    panelTitles, figTitle, use_log_y, jitter_width, marker_size)

    figHandle = figure('Color', 'w', 'Position', [50 80 1750 820]);

    tl = tiledlayout(figHandle, 1, 2, ...
        'Padding', 'compact', ...
        'TileSpacing', 'compact');

    allVals = [flattenValueCellsLocal(valuesDir1), flattenValueCellsLocal(valuesDir2)];
    condColors = getConditionColorsLocal(numel(conditionShortLabels));

    ax1 = nexttile(tl, 1);
    legendHandles = plotOneSplitDirPanelLocal(ax1, valuesDir1, categoryLabels, ...
        conditionShortLabels, condColors, panelTitles{1}, use_log_y, ...
        jitter_width, marker_size, allVals);

    ax2 = nexttile(tl, 2);
    plotOneSplitDirPanelLocal(ax2, valuesDir2, categoryLabels, ...
        conditionShortLabels, condColors, panelTitles{2}, use_log_y, ...
        jitter_width, marker_size, allVals);

    sgtitle(tl, figTitle, 'Interpreter', 'none', 'FontWeight', 'bold');

    lgd = legend(ax1, legendHandles, conditionShortLabels, ...
        'Location', 'southoutside', ...
        'Orientation', 'horizontal', ...
        'Interpreter', 'none');
    lgd.Box = 'off';
end

function legendHandles = plotOneSplitDirPanelLocal(ax, valuesCell, categoryLabels, ...
    conditionShortLabels, condColors, panelTitle, use_log_y, jitter_width, ...
    marker_size, allValsForYLim)

    axes(ax);
    hold(ax, 'on');

    numCond = numel(conditionShortLabels);
    numCategories = numel(categoryLabels);

    condSpacing = 0.78;
    categoryGap = 1.45;

    xPositions = nan(numCond, numCategories);

    for c = 1:numCategories
        baseX = (c - 1) * (numCond * condSpacing + categoryGap);
        xPositions(:, c) = baseX + (1:numCond) * condSpacing;
    end

    categoryCenters = mean(xPositions, 1);

    legendHandles = gobjects(1, numCond);
    for p = 1:numCond
        legendHandles(p) = scatter(ax, nan, nan, marker_size, condColors(p, :), ...
            'filled', 'MarkerEdgeColor', 'none');
    end

    for c = 1:numCategories
        for p = 1:numCond
            vals = valuesCell{p, c};
            vals = vals(isfinite(vals));

            if use_log_y
                vals = vals(vals > 0);
            end

            if isempty(vals)
                continue;
            end

            x0 = xPositions(p, c);
            x = x0 + jitter_width * (rand(size(vals)) - 0.5);

            sc = scatter(ax, x, vals, marker_size, condColors(p, :), ...
                'filled', 'MarkerEdgeColor', 'none');

            try
                sc.MarkerFaceAlpha = 0.70;
            catch
            end
        end
    end

    applyTimescaleAxisLocal(ax, allValsForYLim, use_log_y);

    yl = ylim(ax);

    for c = 1:numCategories
        for p = 1:numCond
            xline(ax, xPositions(p, c), ':', ...
                'Color', [0.82 0.82 0.82], ...
                'LineWidth', 0.7);
        end
    end

    for c = 1:(numCategories - 1)
        sepX = (xPositions(end, c) + xPositions(1, c + 1)) / 2;
        xline(ax, sepX, '-', ...
            'Color', [0.55 0.55 0.55], ...
            'LineWidth', 0.9);
    end

    ylim(ax, yl);

    set(ax, ...
        'XTick', categoryCenters, ...
        'XTickLabel', categoryLabels, ...
        'TickDir', 'out');

    xtickangle(ax, 30);

    xlim(ax, [min(xPositions(:)) - 0.8, max(xPositions(:)) + 0.8]);

    ax.XGrid = 'off';
    ax.YGrid = 'off';

    ylabel(ax, 'Timescale (ms)');
    xlabel(ax, 'Latent category');
    title(ax, panelTitle, 'Interpreter', 'none', 'FontWeight', 'bold');

    box(ax, 'off');
    hold(ax, 'off');
end

function figHandle = plotOneGroupTimescaleDistributionLocal( ...
    valueCells, categoryLabels, figTitle, use_log_y, jitter_width, marker_size)

    figHandle = figure('Color', 'w', 'Position', [100 100 900 650]);
    ax = axes(figHandle);
    hold(ax, 'on');

    numCategories = numel(categoryLabels);
    allVals = [];

    categoryX = 1:numCategories;

    for c = 1:numCategories
        vals = valueCells{c};
        vals = vals(isfinite(vals));

        if use_log_y
            vals = vals(vals > 0);
        end

        allVals = [allVals, vals];

        if isempty(vals)
            continue;
        end

        x = categoryX(c) + jitter_width * (rand(size(vals)) - 0.5);

        sc = scatter(ax, x, vals, marker_size, 'filled', ...
            'MarkerEdgeColor', 'none');

        try
            sc.MarkerFaceAlpha = 0.70;
        catch
        end
    end

    applyTimescaleAxisLocal(ax, allVals, use_log_y);

    yl = ylim(ax);

    for c = 1:numCategories
        xline(ax, categoryX(c), ':', ...
            'Color', [0.82 0.82 0.82], ...
            'LineWidth', 0.7);
    end

    ylim(ax, yl);

    set(ax, ...
        'XTick', categoryX, ...
        'XTickLabel', categoryLabels, ...
        'TickDir', 'out');

    xtickangle(ax, 30);

    xlim(ax, [0.4, numCategories + 0.6]);

    ax.XGrid = 'off';
    ax.YGrid = 'off';

    ylabel(ax, 'Timescale (ms)');
    xlabel(ax, 'Latent category');
    title(ax, figTitle, 'Interpreter', 'none');

    box(ax, 'off');
    hold(ax, 'off');
end