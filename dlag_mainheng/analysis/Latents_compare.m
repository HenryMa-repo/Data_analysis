%% compare across, within, feedforward and feedback percentage latents, and shared var explained
%% two options, remove or not remove stim related latents
%% DSL, analyze diff kinds of latents distribution, this code comparable for all conditions trained Dlag and each condition trained Dlag
%% sum-all-conditions category summary is displayed in two panels:
%% stim_dir1 and stim_dir2, with readable stimulus-condition legends
%% ------------------------------------------------------------------------
%% NEW ADDITION:
%% In non-condition mode, additionally compute split DSL under multiple
%% trial grouping schemes:
%%   1) by stim_dir
%%   2) by stim_name x stim_dir (4 groups)
%%   3) by condition id (current paradigm: 16 groups)
%% These are stored in fields such as:
%%   DSL.indiv_bystimdir / DSL.rawlogical_bystimdir / DSL.logical_bystimdir
%%   DSL.indiv_bystimnamedir / DSL.rawlogical_bystimnamedir / DSL.logical_bystimnamedir
%%   DSL.indiv_bycondition / DSL.rawlogical_bycondition / DSL.logical_bycondition
%% ------------------------------------------------------------------------
%% NEW ADDITION 2:
%% In non-condition mode, compute condition-specific posterior-based shared
%% variance explained from the pooled all-condition DLAG model using
%% posterior second moment:
%%   rawSharedVariance_j(c) = ||c_j||^2 * mean(xsm_j.^2 + posteriorVar_j)
%% This is designed to be the most reasonable posterior-based extension of
%% computeVarExp_dlag(total=false) for pooled models.
%% ------------------------------------------------------------------------

clc; clear;

data_content = 'raw_count';
% options:
% raw_count, raw_fr, z_within_trial, z_within_condition,
% z_across_conditions, demean_count_within_trial, demean_fr_within_trial, demean_pooledsd_within_condition
data_condtion = [1:16];
runIdx = 1;
DSL_threshold = 0.3;

% -------------------------------------------------------------------------
% Load all-run stimulus metadata so condition IDs can be decoded later
% -------------------------------------------------------------------------
dat_file = '.\model_data_allruns';
fprintf('Reading from %s\n', dat_file);
load(dat_file, 'model_data_allruns');

stim_tag = '_2[Gpl2_2c_2sz_400_2_200isi]';

all_run_tags = get_all_run_tags(model_data_allruns);
run_idx = find(strcmp(all_run_tags, stim_tag));

if isempty(run_idx)
    error('Requested stim_tag not found: %s', stim_tag);
end
if numel(run_idx) > 1
    error('Duplicate stim_tag found: %s', stim_tag);
end

condition_full = model_data_allruns{run_idx}.conditions_full;

% -------------------------------------------------------------------------
% Condition mode
% -------------------------------------------------------------------------
if isempty(data_condtion)
    use_condition_mode = false;
    condition_list = [];
    numConditions = 1;
else
    use_condition_mode = true;
    condition_list = data_condtion(:)';
    numConditions = numel(condition_list);
end

AllConditionResults = struct([]);

for cond_i = 1:numConditions

    if use_condition_mode
        this_condition = condition_list(cond_i);
        baseDir = ['./FA_Dlag_', data_content, '_condition', num2str(this_condition)];
    else
        this_condition = [];
        baseDir = ['./FA_Dlag_', data_content];
    end

    tempfname = sprintf('%s/mat_results/run%03d', baseDir, runIdx);

    files = dir(fullfile(tempfname, 'bestmodel*'));
    if isempty(files)
        error('No bestmodel* file found in %s', tempfname);
    end
    filename = fullfile(tempfname, files(1).name);
    load(filename, "bestModel", "res", "seqEst", "varexp", "gp_params");

    files = dir(fullfile(tempfname, 'bootstrapResults*'));
    if isempty(files)
        error('No bootstrapResults* file found in %s', tempfname);
    end
    filename = fullfile(tempfname, files(1).name);
    load(filename, "ambiguousIdxs");

    % -------------------------------------------------------------
    % Original DSL
    % -------------------------------------------------------------
    [DSL, DSL_hist_stats, DSL_figs] = computeDSL_dlag( ...
        bestModel.xDim_across, ...
        bestModel.xDim_within, ...
        seqEst, ...
        DSL_threshold);

    % -------------------------------------------------------------
    % NEW: non-condition mode only, split trials by multiple schemes
    % -------------------------------------------------------------
    DSL_hist_stats_bystimdir = struct([]);
    DSL_figs_bystimdir = gobjects(0);
    Stats_bystimdir = [];
    Stats_figs_bystimdir = struct([]);

    DSL_hist_stats_bystimnamedir = struct([]);
    DSL_figs_bystimnamedir = gobjects(0);
    Stats_bystimnamedir = [];
    Stats_figs_bystimnamedir = struct([]);

    DSL_hist_stats_bycondition = struct([]);
    DSL_figs_bycondition = gobjects(0);
    Stats_bycondition = [];
    Stats_figs_bycondition = struct([]);

    trial_condition_ids = [];
    trial_stim_dir_values = [];
    trial_stimnamedir_group_ids = [];
    trial_stimnamedir_group_labels = {};
    trial_condition_group_labels = {};

    % -------------------------------------------------------------
    % NEW: condition-specific posterior shared variance analysis
    % -------------------------------------------------------------
    CondPosteriorVarExp = [];
    CondPosteriorVarExp_figs = struct([]);

    if use_condition_mode
        main_filtered_display_name = 'condition-specific DSL filtered';
    else
        main_filtered_display_name = 'all-trials DSL filtered';
    end
    main_filtered_file_tag = makeFilterModeFileTag(main_filtered_display_name);

    if ~use_condition_mode
        trial_condition_ids = extractTrialConditionIdsFromConditionFull(condition_full, seqEst);
        trial_stim_dir_values = mapTrialConditionToStimDir(trial_condition_ids, condition_full);

        % ---------------------------------------------------------
        % 1) split by stim_dir
        % ---------------------------------------------------------
        [DSL, DSL_hist_stats_bystimdir, DSL_figs_bystimdir] = ...
            computeDSL_dlag_bytrialgroups( ...
            DSL, ...
            bestModel.xDim_across, ...
            bestModel.xDim_within, ...
            seqEst, ...
            DSL_threshold, ...
            trial_stim_dir_values, ...
            {}, ...
            'bystimdir', ...
            'stim dir split');

        DSL_for_bystimdir_stats = DSL;
        DSL_for_bystimdir_stats.logical = DSL.logical_bystimdir;
        [Stats_bystimdir, Stats_figs_bystimdir] = summarizeLatentCategories( ...
            varexp.indiv, ...
            bestModel.xDim_across, ...
            bestModel.xDim_within, ...
            gp_params, ...
            ambiguousIdxs, ...
            DSL_for_bystimdir_stats, ...
            true, ...
            'stim_dir DSL filtered');

        % ---------------------------------------------------------
        % 2) split by stim_name x stim_dir (4 groups)
        % ---------------------------------------------------------
        [trial_stimnamedir_group_ids, trial_stimnamedir_group_labels] = ...
            mapTrialConditionToStimNameDirGroup(trial_condition_ids, condition_full);

        [DSL, DSL_hist_stats_bystimnamedir, DSL_figs_bystimnamedir] = ...
            computeDSL_dlag_bytrialgroups( ...
            DSL, ...
            bestModel.xDim_across, ...
            bestModel.xDim_within, ...
            seqEst, ...
            DSL_threshold, ...
            trial_stimnamedir_group_ids, ...
            trial_stimnamedir_group_labels, ...
            'bystimnamedir', ...
            'stim name x stim dir split');

        DSL_for_bystimnamedir_stats = DSL;
        DSL_for_bystimnamedir_stats.logical = DSL.logical_bystimnamedir;
        [Stats_bystimnamedir, Stats_figs_bystimnamedir] = summarizeLatentCategories( ...
            varexp.indiv, ...
            bestModel.xDim_across, ...
            bestModel.xDim_within, ...
            gp_params, ...
            ambiguousIdxs, ...
            DSL_for_bystimnamedir_stats, ...
            true, ...
            'stim_name x stim_dir DSL filtered');

        % ---------------------------------------------------------
        % 3) split by condition id (current paradigm: 16 groups)
        % ---------------------------------------------------------
        trial_condition_group_labels = buildConditionGroupLabels(condition_full);

        [DSL, DSL_hist_stats_bycondition, DSL_figs_bycondition] = ...
            computeDSL_dlag_bytrialgroups( ...
            DSL, ...
            bestModel.xDim_across, ...
            bestModel.xDim_within, ...
            seqEst, ...
            DSL_threshold, ...
            trial_condition_ids, ...
            trial_condition_group_labels, ...
            'bycondition', ...
            'condition split');

        DSL_for_bycondition_stats = DSL;
        DSL_for_bycondition_stats.logical = DSL.logical_bycondition;
        [Stats_bycondition, Stats_figs_bycondition] = summarizeLatentCategories( ...
            varexp.indiv, ...
            bestModel.xDim_across, ...
            bestModel.xDim_within, ...
            gp_params, ...
            ambiguousIdxs, ...
            DSL_for_bycondition_stats, ...
            true, ...
            'condition DSL filtered');

        [CondPosteriorVarExp, CondPosteriorVarExp_figs] = ...
            analyzeConditionSpecificPosteriorVarExpFromPooled( ...
                res, ...
                varexp.indiv, ...
                seqEst, ...
                condition_full, ...
                trial_condition_ids, ...
                gp_params, ...
                ambiguousIdxs, ...
                DSL, ...
                true);
    end

    % -------------------------------------------------------------
    % Original latent category summary
    % -------------------------------------------------------------
    [Stats, Stats_figs] = summarizeLatentCategories( ...
        varexp.indiv, ...
        bestModel.xDim_across, ...
        bestModel.xDim_within, ...
        gp_params, ...
        ambiguousIdxs, ...
        DSL, ...
        true, ...
        main_filtered_display_name);

    % -------------------------------------------------------------
    % Save figures
    % -------------------------------------------------------------
    saveFigureVector(DSL_figs, tempfname, 'DSL_distribution_group%d');

    saveCategoryFigureSet(Stats_figs, tempfname, struct( ...
        'allLatentPct', sprintf('percentage_latents_all_group%%d'), ...
        'allSharedVarPct', sprintf('percentage_shared_variance_all_group%%d'), ...
        'filteredLatentPct', sprintf('percentage_latents_%s_group%%d', main_filtered_file_tag), ...
        'filteredSharedVarPct', sprintf('percentage_shared_variance_%s_group%%d', main_filtered_file_tag)));

    if ~use_condition_mode
        saveFigureVector(DSL_figs_bystimdir, tempfname, 'stim_dir_splited_DSL_distribution_group%d');
        saveFigureVector(DSL_figs_bystimnamedir, tempfname, 'stim_name_stim_dir_splited_DSL_distribution_group%d');
        saveFigureVector(DSL_figs_bycondition, tempfname, 'condition_splited_DSL_distribution_group%d');

        saveCategoryFigureSet(Stats_figs_bystimdir, tempfname, struct( ...
            'filteredLatentPct', 'stim_dir_splited_percentage_latents_stim_dir_DSL_filtered_group%d', ...
            'filteredSharedVarPct', 'stim_dir_splited_percentage_shared_variance_stim_dir_DSL_filtered_group%d'));

        saveCategoryFigureSet(Stats_figs_bystimnamedir, tempfname, struct( ...
            'filteredLatentPct', 'stim_name_stim_dir_splited_percentage_latents_stim_name_stim_dir_DSL_filtered_group%d', ...
            'filteredSharedVarPct', 'stim_name_stim_dir_splited_percentage_shared_variance_stim_name_stim_dir_DSL_filtered_group%d'));

        saveCategoryFigureSet(Stats_figs_bycondition, tempfname, struct( ...
            'filteredLatentPct', 'condition_splited_percentage_latents_condition_DSL_filtered_group%d', ...
            'filteredSharedVarPct', 'condition_splited_percentage_shared_variance_condition_DSL_filtered_group%d'));

        saveCategoryFigureSet(CondPosteriorVarExp_figs, tempfname, struct( ...
            'allSharedVarPct', 'condition_specific_posterior_shared_varexp_all_group%d', ...
            'allTrialsDSLfilteredSharedVarPct', 'condition_specific_posterior_shared_varexp_all_trials_DSL_filtered_group%d', ...
            'stimDirDSLfilteredSharedVarPct', 'condition_specific_posterior_shared_varexp_stim_dir_DSL_filtered_group%d', ...
            'stimNameDirDSLfilteredSharedVarPct', 'condition_specific_posterior_shared_varexp_stim_name_stim_dir_DSL_filtered_group%d', ...
            'conditionDSLfilteredSharedVarPct', 'condition_specific_posterior_shared_varexp_condition_DSL_filtered_group%d'));
    end

    % -------------------------------------------------------------
    % Save mat
    % add new ones in non-condition mode
    % -------------------------------------------------------------
    if ~use_condition_mode
        save(fullfile(tempfname, 'DSL_and_latent_category_stats.mat'), ...
             'DSL', 'DSL_hist_stats', 'Stats', ...
             'DSL_hist_stats_bystimdir', 'Stats_bystimdir', ...
             'DSL_hist_stats_bystimnamedir', 'Stats_bystimnamedir', ...
             'DSL_hist_stats_bycondition', 'Stats_bycondition', ...
             'CondPosteriorVarExp', ...
             'trial_condition_ids', 'trial_stim_dir_values', ...
             'trial_stimnamedir_group_ids', 'trial_stimnamedir_group_labels', ...
             'trial_condition_group_labels', ...
             'DSL_threshold', 'this_condition', 'data_content', 'runIdx', ...
             'bestModel', 'gp_params', 'ambiguousIdxs');
    else
        save(fullfile(tempfname, 'DSL_and_latent_category_stats.mat'), ...
             'DSL', 'DSL_hist_stats', 'Stats', ...
             'DSL_threshold', 'this_condition', 'data_content', 'runIdx', ...
             'bestModel', 'gp_params', 'ambiguousIdxs');
    end

    % -------------------------------------------------------------
    % Store for sum-all-conditions summary
    % -------------------------------------------------------------
    if use_condition_mode
        AllConditionResults(cond_i).condition = this_condition;
        AllConditionResults(cond_i).baseDir = baseDir;
        AllConditionResults(cond_i).tempfname = tempfname;
        AllConditionResults(cond_i).DSL = DSL;
        AllConditionResults(cond_i).DSL_hist_stats = DSL_hist_stats;
        AllConditionResults(cond_i).Stats = Stats;
        AllConditionResults(cond_i).xDim_across = bestModel.xDim_across;
        AllConditionResults(cond_i).xDim_within = bestModel.xDim_within;
    end
end

if use_condition_mode

    % =============================================================
    % Sum all conditions: DSL distributions
    % =============================================================
    [SummaryDSL, SummaryDSL_figs] = summarizeAllConditionsDSL(AllConditionResults, DSL_threshold);

    saveFigureVector(SummaryDSL_figs, '.', sprintf('%s_sum_all_conditions_DSL_distribution_group%%d', data_content));

    % =============================================================
    % Sum all conditions: percentage latents / shared variance
    % Split into stim_dir1 and stim_dir2 panels using condition_full
    % =============================================================
    [SummaryCategory, SummaryCategory_figs] = summarizeAllConditionsCategories(AllConditionResults, condition_list, condition_full);

    saveCategoryFigureSet(SummaryCategory_figs, '.', struct( ...
        'allLatentPct', sprintf('%s_sum_all_conditions_percentage_latents_all_group%%d', data_content), ...
        'allSharedVarPct', sprintf('%s_sum_all_conditions_percentage_shared_variance_all_group%%d', data_content), ...
        'filteredLatentPct', sprintf('%s_sum_all_conditions_percentage_latents_%s_group%%d', data_content, SummaryCategory.meta.filteredFileTag), ...
        'filteredSharedVarPct', sprintf('%s_sum_all_conditions_percentage_shared_variance_%s_group%%d', data_content, SummaryCategory.meta.filteredFileTag)));

    save(sprintf('%s_sum_all_conditions_DSL_and_latent_category_stats.mat', data_content), ...
         'AllConditionResults', 'SummaryDSL', 'SummaryCategory', ...
         'condition_list', 'DSL_threshold', 'data_content', 'runIdx', 'stim_tag');
end


close all


function [DSL, histStats, figHandles] = computeDSL_dlag(xDim_across, xDim_within, seqEst, DSL_threshold)
% computeDSL_dlag
%
% Compute DSL (latent reproducibility score) for each latent in each group
% from inferred DLAG latent trajectories stored in seqEst(n).xsm.

    if nargin < 4
        error('Usage: DSL = computeDSL_dlag(xDim_across, xDim_within, seqEst, DSL_threshold)');
    end

    if ~isscalar(xDim_across) || xDim_across < 0 || mod(xDim_across,1) ~= 0
        error('xDim_across must be a nonnegative integer scalar.');
    end

    if ~isvector(xDim_within) || any(xDim_within < 0) || any(mod(xDim_within,1) ~= 0)
        error('xDim_within must be a vector of nonnegative integers.');
    end

    numGroups = numel(xDim_within);
    localDims = xDim_across + xDim_within(:)';
    totalDims = sum(localDims);

    if isempty(seqEst)
        error('seqEst is empty.');
    end

    for n = 1:numel(seqEst)
        if ~isfield(seqEst(n), 'xsm')
            error('seqEst(%d) is missing field xsm.', n);
        end
        [nRows, Tn] = size(seqEst(n).xsm);
        if nRows ~= totalDims
            error('seqEst(%d).xsm has %d rows, expected %d.', n, nRows, totalDims);
        end
        if Tn < 2
            error('Each trial must contain at least 2 time points.');
        end
    end

    if isfield(seqEst, 'trialId')
        try
            [~, ord] = sort([seqEst.trialId]);
            seqEst = seqEst(ord);
        catch
        end
    end

    Ntr = numel(seqEst);

    blockStart = cumsum([1, localDims(1:end-1)]);
    blockEnd   = cumsum(localDims);

    allT = arrayfun(@(s) size(s.xsm,2), seqEst);
    Tmin = min(allT);
    Tmax = max(allT);

    lags = -(Tmin-1):(Tmin-1);

    DSL.indiv      = cell(1, numGroups);
    DSL.rawlogical = cell(1, numGroups);
    DSL.logical    = cell(1, numGroups);

    histStatsTemplate = struct( ...
        'values', [], ...
        'edges', [], ...
        'centers', [], ...
        'counts', [], ...
        'percentages', [], ...
        'threshold', []);

    histStats = repmat(histStatsTemplate, 1, numGroups);
    figHandles = gobjects(1, numGroups);

    for g = 1:numGroups
        DSL.indiv{g}      = nan(1, localDims(g));
        DSL.rawlogical{g} = zeros(1, localDims(g));
        DSL.logical{g}    = zeros(1, localDims(g));

        rows_g = blockStart(g):blockEnd(g);

        for j = 1:localDims(g)
            rowIdx = rows_g(j);

            Xpad = zeros(Ntr, Tmax);
            Mpad = false(Ntr, Tmax);

            for tr = 1:Ntr
                xt = seqEst(tr).xsm(rowIdx, :);
                Ttr = numel(xt);
                Xpad(tr, 1:Ttr) = xt;
                Mpad(tr, 1:Ttr) = true;
            end

            ACG_mean = mean_acg_matrix(Xpad, Mpad, lags);

            if Ntr >= 2
                CCG_mean = mean_ccg_adjacent_pairs_matrix(Xpad, Mpad, lags);
            else
                CCG_mean = nan(1, numel(lags));
            end

            if all(isnan(ACG_mean)) || all(isnan(CCG_mean))
                dsl_val = NaN;
            else
                numerator   = trapz(lags, abs(ACG_mean - CCG_mean));
                denominator = trapz(lags, abs(ACG_mean));

                if ~isfinite(denominator) || denominator <= eps
                    dsl_val = NaN;
                else
                    dsl_val = 1 - numerator / denominator;
                end
            end

            DSL.indiv{g}(j) = dsl_val;

            if isnan(dsl_val) || dsl_val < DSL_threshold
                DSL.rawlogical{g}(j) = 1;
            else
                DSL.rawlogical{g}(j) = 0;
            end
        end
    end

    for g = 1:numGroups
        DSL.logical{g} = DSL.rawlogical{g};
    end

    for j = 1:xDim_across
        all_groups_flagged = true;
        for g = 1:numGroups
            all_groups_flagged = all_groups_flagged && (DSL.rawlogical{g}(j) == 1);
        end
        for g = 1:numGroups
            DSL.logical{g}(j) = all_groups_flagged;
        end
    end

    for g = 1:numGroups
        vals = DSL.indiv{g};
        [tmpFig, tmpHist] = plotDSLHistogram(vals, DSL_threshold, ...
            sprintf('Group %d DSL distribution', g));

        figHandles(g) = tmpFig;
        histStats(g) = tmpHist;
    end
end


function [DSL, histStatsBySplit, figHandlesBySplit] = computeDSL_dlag_bytrialgroups( ...
    DSL, xDim_across, xDim_within, seqEst, DSL_threshold, ...
    trial_group_ids, group_labels, field_suffix, split_title)

    if numel(trial_group_ids) ~= numel(seqEst)
        error('Length of trial_group_ids must match numel(seqEst).');
    end

    numGroups = numel(xDim_within);
    localDims = xDim_across + xDim_within(:)';

    trial_group_ids = reshape(trial_group_ids, 1, []);
    validMask = ~isnan(trial_group_ids);
    rawGroupValues = unique(trial_group_ids(validMask));

    if isempty(rawGroupValues)
        error('No valid split groups found for field suffix %s.', field_suffix);
    end

    if isnumeric(rawGroupValues)
        groupValues = sort(rawGroupValues(:)');
    else
        groupValues = unique(trial_group_ids(validMask), 'stable');
    end
    numSplits = numel(groupValues);

    if isempty(group_labels)
        group_labels_use = arrayfun(@(k) ...
            sprintf('group%d = %s', k, formatSummaryValue(groupValues(k))), ...
            1:numSplits, 'UniformOutput', false);
    else
        group_labels = cellstr(string(group_labels));
        if isnumeric(groupValues) && all(isfinite(groupValues)) && ...
                all(abs(groupValues - round(groupValues)) < 1e-10) && max(groupValues) <= numel(group_labels)
            group_labels_use = group_labels(round(groupValues));
        elseif numel(group_labels) == numSplits
            group_labels_use = group_labels(:)';
        else
            error('group_labels length must match the number of unique trial groups or cover the indexed group IDs.');
        end
    end

    indivField   = ['indiv_' field_suffix];
    rawField     = ['rawlogical_' field_suffix];
    logicalField = ['logical_' field_suffix];
    metaField    = ['meta_' field_suffix];

    DSL.(metaField).groupValues = groupValues(:)';
    DSL.(metaField).groupLabels = group_labels_use(:)';
    DSL.(metaField).threshold   = DSL_threshold;

    DSL.(indivField)   = cell(numSplits, numGroups);
    DSL.(rawField)     = cell(numSplits, numGroups);
    DSL.(logicalField) = cell(1, numGroups);

    histTemplate = struct( ...
        'values', [], ...
        'edges', [], ...
        'centers', [], ...
        'counts', [], ...
        'percentages', [], ...
        'threshold', []);
    histStatsBySplit = repmat(histTemplate, numSplits, numGroups);
    figHandlesBySplit = gobjects(1, numGroups);

    % ---------------------------------------------------------
    % Compute DSL separately within each split
    % ---------------------------------------------------------
    for s = 1:numSplits
        thisValue = groupValues(s);
        keepTrials = (trial_group_ids == thisValue);

        if ~any(keepTrials)
            error('No trials found for split value %s.', formatSummaryValue(thisValue));
        end

        seqEst_thisSplit = seqEst(keepTrials);

        [DSL_thisSplit, ~, ~] = computeDSL_dlag( ...
            xDim_across, ...
            xDim_within, ...
            seqEst_thisSplit, ...
            DSL_threshold);

        for g = 1:numGroups
            DSL.(indivField){s, g} = DSL_thisSplit.indiv{g};
            DSL.(rawField){s, g} = DSL_thisSplit.rawlogical{g};
            histStatsBySplit(s, g) = buildDSLHistogramStats(DSL_thisSplit.indiv{g}, DSL_threshold);
        end
    end

    % ---------------------------------------------------------
    % First AND across split groups within each neural group
    % ---------------------------------------------------------
    for g = 1:numGroups
        tmpLogical = true(1, localDims(g));
        for s = 1:numSplits
            tmpLogical = tmpLogical & (DSL.(rawField){s, g} ~= 0);
        end
        DSL.(logicalField){g} = double(tmpLogical);
    end

    % ---------------------------------------------------------
    % Then for across latents, also AND across neural groups
    % ---------------------------------------------------------
    for j = 1:xDim_across
        all_groups_flagged = true;
        for g = 1:numGroups
            all_groups_flagged = all_groups_flagged && (DSL.(logicalField){g}(j) ~= 0);
        end
        for g = 1:numGroups
            DSL.(logicalField){g}(j) = double(all_groups_flagged);
        end
    end

    % ---------------------------------------------------------
    % Raw DSL distribution figures for each split
    % ---------------------------------------------------------
    for g = 1:numGroups
        rawValueCells = cell(1, numSplits);
        for s = 1:numSplits
            rawValueCells{s} = reshape(DSL.(indivField){s, g}, 1, []);
        end

        figHandlesBySplit(g) = plotDSLHistogramPanels( ...
            rawValueCells, ...
            group_labels_use, ...
            DSL_threshold, ...
            sprintf('Group %d DSL distribution (%s)', g, split_title));
    end
end


function trial_condition_ids = extractTrialConditionIdsFromConditionFull(condition_full, seqEst)
% Build trial -> condition mapping using condition_full(k).trial_indices.
%
% If seqEst has trialId, use it to align.
% Otherwise assume seqEst order corresponds to trial indices 1:numel(seqEst).

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

        if any(idx < 1) || any(mod(idx,1) ~= 0)
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

    if any(seq_trial_ids < 1) || any(mod(seq_trial_ids,1) ~= 0)
        error('seqEst trial IDs are invalid.');
    end

    if max(seq_trial_ids) > numel(trial_to_condition)
        error('seqEst contains trial IDs larger than max trial index in condition_full.trial_indices.');
    end

    trial_condition_ids = trial_to_condition(seq_trial_ids);

    if any(isnan(trial_condition_ids))
        missingTrial = seq_trial_ids(find(isnan(trial_condition_ids), 1));
        error('Could not map seqEst trial %d to any condition using condition_full.trial_indices.', missingTrial);
    end
end

function trial_stim_dir_values = mapTrialConditionToStimDir(trial_condition_ids, condition_full)

trial_condition_ids = trial_condition_ids(:)';
trial_stim_dir_values = nan(size(trial_condition_ids));

for tr = 1:numel(trial_condition_ids)
    condID = trial_condition_ids(tr);

    if condID < 1 || condID > numel(condition_full)
        error('Trial %d has invalid condition ID %d.', tr, condID);
    end

    trial_stim_dir_values(tr) = getConditionEffectiveDirCanonicalLocal(condition_full(condID), condID);
end

end

function [trial_group_ids, group_labels] = mapTrialConditionToStimNameDirGroup(trial_condition_ids, condition_full)

dirVals = extractStimDirLevelsFromConditionFull(condition_full);

if numel(dirVals) ~= 2
    error('stim_name x stim_dir split currently expects exactly 2 unique canonical dir values, found %d: %s.', ...
        numel(dirVals), mat2str(dirVals));
end

group_labels = { ...
    sprintf('grating_dir = %s', formatSummaryValue(dirVals(1))), ...
    sprintf('grating_dir = %s', formatSummaryValue(dirVals(2))), ...
    sprintf('plaid_dir = %s', formatSummaryValue(dirVals(1))), ...
    sprintf('plaid_dir = %s', formatSummaryValue(dirVals(2)))};

trial_condition_ids = trial_condition_ids(:)';
trial_group_ids = nan(size(trial_condition_ids));

for tr = 1:numel(trial_condition_ids)
    condID = trial_condition_ids(tr);

    if condID < 1 || condID > numel(condition_full)
        error('Trial %d has invalid condition ID %d.', tr, condID);
    end

    currStim = lower(string(condition_full(condID).stim_name));
    dirVal = getConditionEffectiveDirCanonicalLocal(condition_full(condID), condID);

    if currStim == "grating"
        stimCode = 1;
    elseif currStim == "plaid"
        stimCode = 2;
    else
        error('Unsupported stim_name in condition_full(%d): %s', condID, char(currStim));
    end

    dirCode = find(abs(dirVals - dirVal) < 1e-10, 1);

    if isempty(dirCode)
        error('Could not map canonical dir value %s to stim_name x stim_dir group.', ...
            formatSummaryValue(dirVal));
    end

    trial_group_ids(tr) = (stimCode - 1) * 2 + dirCode;
end

end


function group_labels = buildConditionGroupLabels(condition_full)

nCond = numel(condition_full);
group_labels = cell(1, nCond);

for condID = 1:nCond
    currStim = lower(string(condition_full(condID).stim_name));

    if currStim == "grating"
        stimShort = 'G';
        dirVal = getConditionEffectiveDirCanonicalLocal(condition_full(condID), condID);
    elseif currStim == "plaid"
        stimShort = 'P';
        dirVal = getConditionEffectiveDirCanonicalLocal(condition_full(condID), condID);
    else
        stimShort = char(currStim);
        dirVal = NaN;
    end

    if isfield(condition_full(condID), 'size')
        sizeStr = formatSummaryValue(condition_full(condID).size);
    else
        sizeStr = 'NA';
    end

    if isfield(condition_full(condID), 'contrast')
        contrastStr = formatSummaryValue(condition_full(condID).contrast);
    else
        contrastStr = 'NA';
    end

    group_labels{condID} = sprintf( ...
        'cond%02d-%s-dir%s-sz%s-ct%s', ...
        condID, stimShort, formatSummaryValue(dirVal), sizeStr, contrastStr);
end

end


function dirVals = extractStimDirLevelsFromConditionFull(condition_full)

dirVals = nan(1, numel(condition_full));

for condID = 1:numel(condition_full)
    dirVals(condID) = getConditionEffectiveDirCanonicalLocal(condition_full(condID), condID);
end

dirVals = dirVals(isfinite(dirVals));
dirVals = unique(dirVals);
dirVals = sort(dirVals(:)');

end


function ACG_mean = mean_acg_matrix(X, M, lags)
    Ntr = size(X, 1);
    ACG_mean = nan(1, numel(lags));

    den = sum((X.^2) .* M, 2);
    validDen = den > eps;

    for k = 1:numel(lags)
        ell = lags(k);

        if ell >= 0
            A  = X(:, 1+ell:end);
            B  = X(:, 1:end-ell);
            Vm = M(:, 1+ell:end) & M(:, 1:end-ell);
        else
            e  = -ell;
            A  = X(:, 1:end-e);
            B  = X(:, 1+e:end);
            Vm = M(:, 1:end-e) & M(:, 1+e:end);
        end

        num = sum((A .* B) .* Vm, 2);

        acg_tr = nan(Ntr, 1);
        acg_tr(validDen) = num(validDen) ./ den(validDen);

        ACG_mean(k) = mean(acg_tr, 'omitnan');
    end
end


function CCG_mean = mean_ccg_adjacent_pairs_matrix(X, M, lags)
    X1 = X(1:end-1, :);
    X2 = X(2:end,   :);
    M1 = M(1:end-1, :);
    M2 = M(2:end,   :);

    Npairs = size(X1, 1);
    CCG_mean = nan(1, numel(lags));

    den1 = sum((X1.^2) .* M1, 2);
    den2 = sum((X2.^2) .* M2, 2);
    den  = sqrt(den1 .* den2);
    validDen = den > eps;

    for k = 1:numel(lags)
        ell = lags(k);

        if ell >= 0
            A  = X1(:, 1+ell:end);
            B  = X2(:, 1:end-ell);
            Vm = M1(:, 1+ell:end) & M2(:, 1:end-ell);
        else
            e  = -ell;
            A  = X1(:, 1:end-e);
            B  = X2(:, 1+e:end);
            Vm = M1(:, 1:end-e) & M2(:, 1+e:end);
        end

        num = sum((A .* B) .* Vm, 2);

        ccg_tr = nan(Npairs, 1);
        ccg_tr(validDen) = num(validDen) ./ den(validDen);

        CCG_mean(k) = mean(ccg_tr, 'omitnan');
    end
end


function [Stats, figHandles] = summarizeLatentCategories(varexp_indiv, xDim_across, xDim_within, gp_params, ambiguousIdxs, DSL, makePlots, filteredDisplayName)

    if nargin < 7
        makePlots = true;
    end
    if nargin < 8 || isempty(filteredDisplayName)
        filteredDisplayName = 'DSL filtered';
    end

    if ~iscell(varexp_indiv)
        error('varexp_indiv must be a cell array.');
    end

    if ~isscalar(xDim_across) || xDim_across < 0 || mod(xDim_across,1) ~= 0
        error('xDim_across must be a nonnegative integer scalar.');
    end

    if ~isvector(xDim_within) || any(xDim_within < 0) || any(mod(xDim_within,1) ~= 0)
        error('xDim_within must be a vector of nonnegative integers.');
    end

    numGroups = numel(xDim_within);

    if numel(varexp_indiv) ~= numGroups
        error('Length of varexp_indiv must match length of xDim_within.');
    end

    if ~isstruct(DSL) || ~isfield(DSL, 'logical')
        error('DSL must be a struct containing field DSL.logical.');
    end

    if ~iscell(DSL.logical) || numel(DSL.logical) ~= numGroups
        error('DSL.logical must be a cell array with one entry per group.');
    end

    localDims = xDim_across + xDim_within(:)';

    for g = 1:numGroups
        if numel(varexp_indiv{g}) ~= localDims(g)
            error('Group %d: varexp_indiv{%d} has wrong length.', g, g);
        end
        if numel(DSL.logical{g}) ~= localDims(g)
            error('Group %d: DSL.logical{%d} has wrong length.', g, g);
        end
    end


    latentClass = classifyDlagLatents(xDim_across, gp_params, ambiguousIdxs);
    labels = latentClass.categoryLabels;
    acrossDelay = latentClass.acrossDelay;
    acrossIdx = latentClass.acrossIdx;
    ffIdx = latentClass.feedforwardIdx;
    fbIdx = latentClass.feedbackIdx;
    ambiguousAll = latentClass.ambiguousIdx;

    Stats = struct();
    Stats.labels = labels;
    Stats.meta.numGroups = numGroups;
    Stats.meta.xDim_across = xDim_across;
    Stats.meta.xDim_within = xDim_within(:)';
    Stats.meta.localDims = localDims;
    Stats.meta.filteredRule = 'Keep latents with DSL.logical == 1; remove DSL.logical == 0.';
    Stats.meta.filteredDisplayName = filteredDisplayName;
    Stats.meta.filteredFileTag = makeFilterModeFileTag(filteredDisplayName);
    Stats.classification.acrossDelay = acrossDelay;
    Stats.classification.acrossIdx = acrossIdx;
    Stats.classification.feedforwardIdx = ffIdx;
    Stats.classification.feedbackIdx = fbIdx;
    Stats.classification.ambiguousIdx = ambiguousAll;

    figHandles = struct([]);

    for g = 1:numGroups
        ve = reshape(varexp_indiv{g}, 1, []);
        keepAll = true(1, localDims(g));
        keepFiltered = reshape(DSL.logical{g}, 1, []) ~= 0;

        categoryMasks = makeDlagCategoryMasks(localDims(g), xDim_across, ffIdx, fbIdx, ambiguousAll);

        Stats.group(g).name = sprintf('Group %d', g);
        Stats.group(g).all = computeOneMethod(ve, keepAll, categoryMasks, labels);
        Stats.group(g).filtered = computeOneMethod(ve, keepFiltered, categoryMasks, labels);
        Stats.group(g).keepMask.all = keepAll;
        Stats.group(g).keepMask.filtered = keepFiltered;
    end

    if makePlots
        for g = 1:numGroups
            figHandles(g).allLatentPct = plotOneBarFigure( ...
                Stats.group(g).all.percentLatents, labels, ...
                sprintf('Group %d - Percentage of latents (all latents)', g), ...
                'Latent category', ...
                'Percentage of latents (%)');

            figHandles(g).allSharedVarPct = plotOneBarFigure( ...
                Stats.group(g).all.percentSharedVariance, labels, ...
                sprintf('Group %d - Percentage of shared variance explained (all latents)', g), ...
                'Latent category', ...
                'Percentage of shared variance explained (%)');

            figHandles(g).filteredLatentPct = plotOneBarFigure( ...
                Stats.group(g).filtered.percentLatents, labels, ...
                sprintf('Group %d - Percentage of latents (%s)', g, filteredDisplayName), ...
                'Latent category', ...
                'Percentage of latents (%)');

            figHandles(g).filteredSharedVarPct = plotOneBarFigure( ...
                Stats.group(g).filtered.percentSharedVariance, labels, ...
                sprintf('Group %d - Percentage of shared variance explained (%s)', g, filteredDisplayName), ...
                'Latent category', ...
                'Percentage of shared variance explained (%)');
        end
    end
end



function latentClass = classifyDlagLatents(xDim_across, gp_params, ambiguousIdxs)

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


function categoryMasks = makeDlagCategoryMasks(localDim, xDim_across, ffIdx, fbIdx, ambiguousIdx)

    acrossMask = false(1, localDim);
    acrossMask(1:xDim_across) = true;

    withinMask = false(1, localDim);
    withinMask(xDim_across+1:end) = true;

    ffMask = false(1, localDim);
    ffMask(ffIdx) = true;

    fbMask = false(1, localDim);
    fbMask(fbIdx) = true;

    ambMask = false(1, localDim);
    ambMask(ambiguousIdx) = true;

    categoryMasks = {acrossMask, withinMask, ffMask, fbMask, ambMask};
end


function keepMask = getDslKeepMask(DSL, fieldName, groupIdx, localDim, defaultKeepValue)

    if nargin < 5
        defaultKeepValue = true;
    end

    if isfield(DSL, fieldName) && ~isempty(DSL.(fieldName)) && numel(DSL.(fieldName)) >= groupIdx && ...
            ~isempty(DSL.(fieldName){groupIdx})
        keepMask = reshape(DSL.(fieldName){groupIdx}, 1, []) ~= 0;
    else
        keepMask = repmat(logical(defaultKeepValue), 1, localDim);
    end

    if numel(keepMask) ~= localDim
        error('DSL.%s{%d} has wrong length. Expected %d, got %d.', ...
            fieldName, groupIdx, localDim, numel(keepMask));
    end
end


function S = computeCategorySummaryFromRaw(rawLatByCondition, keepMask, categoryMasks, labels)

    numConditions = size(rawLatByCondition, 1);
    numCategories = numel(labels);

    rawByCondition = nan(numConditions, numCategories);
    fractionByCondition = nan(numConditions, numCategories);
    percentByCondition = nan(numConditions, numCategories);
    totalByCondition = nan(numConditions, 1);

    keepMask = reshape(keepMask, 1, []) ~= 0;

    for condID = 1:numConditions
        rawLat = rawLatByCondition(condID, :);
        denom = sum(rawLat(keepMask), 'omitnan');
        totalByCondition(condID) = denom;

        if ~(denom > 0)
            continue;
        end

        for k = 1:numCategories
            mask = keepMask & categoryMasks{k};
            rawByCondition(condID, k) = sum(rawLat(mask), 'omitnan');
            fractionByCondition(condID, k) = rawByCondition(condID, k) / denom;
            percentByCondition(condID, k) = 100 * rawByCondition(condID, k) / denom;
        end
    end

    S = struct();
    S.labels = labels;
    S.rawByConditionId = rawByCondition;
    S.fractionByConditionId = fractionByCondition;
    S.percentByConditionId = percentByCondition;
    S.totalByConditionId = totalByCondition;
end


function S = computeOneMethod(ve, keepMask, categoryMasks, labels)
    numCategories = numel(labels);

    counts = zeros(1, numCategories);
    veSums = zeros(1, numCategories);

    for k = 1:numCategories
        mask = keepMask & categoryMasks{k};
        counts(k) = sum(mask);
        veSums(k) = sum(ve(mask));
    end

    denomCount = sum(keepMask);
    denomVE = sum(ve(keepMask));

    if denomCount > 0
        percentLatents = 100 * counts / denomCount;
    else
        percentLatents = zeros(1, numCategories);
    end

    if denomVE > 0
        percentSharedVariance = 100 * veSums / denomVE;
    else
        percentSharedVariance = zeros(1, numCategories);
    end

    S = struct();
    S.labels = labels;
    S.denominator.numLatents = denomCount;
    S.denominator.sharedVariance = denomVE;
    S.counts = counts;
    S.sharedVariance = veSums;
    S.percentLatents = percentLatents;
    S.percentSharedVariance = percentSharedVariance;
end


function figHandle = plotOneBarFigure(values, labels, figTitle, xlab, ylab)
    figHandle = figure;
    bar(values);

    set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels);
    xtickangle(30);

    xlabel(xlab);
    ylabel(ylab);
    title(figTitle);

    ylim([0, max(100, max(values) * 1.1 + eps)]);
end



function [figHandle, histStat] = plotDSLHistogram(values, DSL_threshold, figTitle)
    histStat = buildDSLHistogramStats(values, DSL_threshold);

    figHandle = figure;
    bar(histStat.centers, histStat.percentages, 1);
    hold on;
    xline(DSL_threshold, '--', 'LineWidth', 1.5);
    hold off;

    xlabel('DSL');
    ylabel('Percentage of latents (%)');
    title(figTitle);
    xlim([histStat.edges(1), histStat.edges(end)]);
end


function histStat = buildDSLHistogramStats(values, DSL_threshold, edges)
    values = values(isfinite(values));

    if nargin < 3 || isempty(edges)
        if isempty(values)
            allDSL = [0 1];
        else
            allDSL = values;
        end

        binWidth = 0.05;
        xmin = min(floor(min(allDSL)/binWidth)*binWidth, floor(DSL_threshold/binWidth)*binWidth) - binWidth;
        xmax = max(ceil(max(allDSL)/binWidth)*binWidth,  ceil(DSL_threshold/binWidth)*binWidth) + binWidth;
        edges = xmin:binWidth:xmax;
        if numel(edges) < 2
            edges = [xmin, xmin + binWidth];
        end
    end

    centers = edges(1:end-1) + diff(edges)/2;
    counts = histcounts(values, edges);

    if isempty(values)
        pct = zeros(size(counts));
    else
        pct = 100 * counts / numel(values);
    end

    histStat = struct();
    histStat.values = values;
    histStat.edges = edges;
    histStat.centers = centers;
    histStat.counts = counts;
    histStat.percentages = pct;
    histStat.threshold = DSL_threshold;
end


function figHandle = plotDSLHistogramPanels(valueCells, panelTitles, DSL_threshold, figTitle)

    numPanels = numel(valueCells);

    pooledVals = [];
    for k = 1:numPanels
        vals = valueCells{k};
        vals = vals(isfinite(vals));
        pooledVals = [pooledVals, vals];
    end
    pooledHist = buildDSLHistogramStats(pooledVals, DSL_threshold);
    commonEdges = pooledHist.edges;

    panelHist = cell(1, numPanels);
    ymax = 0;
    for k = 1:numPanels
        panelHist{k} = buildDSLHistogramStats(valueCells{k}, DSL_threshold, commonEdges);
        if ~isempty(panelHist{k}.percentages)
            ymax = max(ymax, max(panelHist{k}.percentages));
        end
    end
    ymax = max(100, ymax * 1.1 + eps);

    [nRows, nCols] = chooseHistogramGrid(numPanels);
    [figWidth, figHeight] = chooseHistogramFigureSize(nRows, nCols, numPanels);
    [titleFontSize, axisFontSize, labelFontSize] = chooseHistogramFontSizes(numPanels);

    figHandle = figure('Color', 'w', 'Units', 'pixels', ...
        'Position', [60, 60, figWidth, figHeight]);
    tl = tiledlayout(figHandle, nRows, nCols, 'Padding', 'compact', 'TileSpacing', 'compact');

    for k = 1:numPanels
        ax = nexttile(tl, k);
        bar(ax, panelHist{k}.centers, panelHist{k}.percentages, 1);
        hold(ax, 'on');
        xline(ax, DSL_threshold, '--', 'LineWidth', 1.4);
        hold(ax, 'off');

        xlim(ax, [commonEdges(1), commonEdges(end)]);
        ylim(ax, [0, ymax]);
        grid(ax, 'on');
        box(ax, 'off');
        ax.FontSize = axisFontSize;
        ax.TitleFontWeight = 'bold';
        ax.TitleFontSizeMultiplier = 1;

        title(ax, formatHistogramPanelTitle(panelTitles{k}, numPanels), ...
            'FontSize', titleFontSize, 'Interpreter', 'none');

        if isBottomRowTile(k, nRows, nCols, numPanels)
            xlabel(ax, 'DSL', 'FontSize', labelFontSize);
        else
            xlabel(ax, '');
        end

        if isLeftColumnTile(k, nCols)
            ylabel(ax, 'Percentage of latents (%)', 'FontSize', labelFontSize);
        else
            ylabel(ax, '');
        end
    end

    sgtitle(tl, figTitle, 'FontSize', max(titleFontSize + 3, 12), 'FontWeight', 'bold');
end




function [nRows, nCols] = chooseHistogramGrid(numPanels)

    if numPanels <= 3
        nRows = 1;
        nCols = numPanels;
    elseif numPanels <= 4
        nRows = 2;
        nCols = 2;
    elseif numPanels <= 6
        nRows = 2;
        nCols = 3;
    elseif numPanels <= 8
        nRows = 2;
        nCols = 4;
    elseif numPanels <= 12
        nRows = 3;
        nCols = 4;
    elseif numPanels <= 16
        nRows = 4;
        nCols = 4;
    else
        nCols = ceil(sqrt(numPanels));
        nRows = ceil(numPanels / nCols);
    end
end


function [figWidth, figHeight] = chooseHistogramFigureSize(nRows, nCols, numPanels)

    if numPanels <= 4
        tileWidth = 420;
        tileHeight = 320;
    elseif numPanels <= 8
        tileWidth = 360;
        tileHeight = 290;
    else
        tileWidth = 340;
        tileHeight = 260;
    end

    figWidth = max(1100, nCols * tileWidth);
    figHeight = max(700, nRows * tileHeight + 80);
end


function [titleFontSize, axisFontSize, labelFontSize] = chooseHistogramFontSizes(numPanels)

    if numPanels <= 4
        titleFontSize = 12;
        axisFontSize = 11;
        labelFontSize = 11;
    elseif numPanels <= 8
        titleFontSize = 11;
        axisFontSize = 10;
        labelFontSize = 10;
    else
        titleFontSize = 10;
        axisFontSize = 9;
        labelFontSize = 9;
    end
end


function tf = isBottomRowTile(tileIdx, nRows, nCols, numPanels)

    rowIdx = ceil(tileIdx / nCols);
    lastUsedRow = ceil(numPanels / nCols);
    tf = (rowIdx == lastUsedRow);
end


function tf = isLeftColumnTile(tileIdx, nCols)

    tf = mod(tileIdx - 1, nCols) == 0;
end


function titleStr = formatHistogramPanelTitle(titleStr, numPanels)

    titleStr = char(string(titleStr));

    if numPanels >= 12
        titleStr = strrep(titleStr, '-dir', sprintf('\ndir'));
        titleStr = strrep(titleStr, '-sz', ' | sz');
        titleStr = strrep(titleStr, '-ct', ' ct');
        titleStr = strrep(titleStr, '_', ' ');
    elseif numPanels >= 8
        titleStr = strrep(titleStr, '-dir', ' | dir');
        titleStr = strrep(titleStr, '-sz', ' | sz');
        titleStr = strrep(titleStr, '-ct', ' ct');
        titleStr = strrep(titleStr, '_', ' ');
    end
end

function saveOneFigure(figHandle, folderPath, baseName)
    savefig(figHandle, fullfile(folderPath, [baseName, '.fig']));
    exportgraphics(figHandle, fullfile(folderPath, [baseName, '.png']), 'Resolution', 400);
end


function saveAndCloseFigure(figHandle, folderPath, baseName)
    if isempty(figHandle) || ~isgraphics(figHandle)
        return;
    end

    saveOneFigure(figHandle, folderPath, baseName);
    close(figHandle);
end


function saveFigureVector(figHandles, folderPath, filenamePattern)
    if isempty(figHandles)
        return;
    end

    for g = 1:numel(figHandles)
        saveAndCloseFigure(figHandles(g), folderPath, sprintf(filenamePattern, g));
    end
end


function saveCategoryFigureSet(figStruct, folderPath, fieldPatterns)
    if isempty(figStruct)
        return;
    end

    fieldsToSave = fieldnames(fieldPatterns);
    allFigureFields = {'allLatentPct', 'allSharedVarPct', 'filteredLatentPct', 'filteredSharedVarPct', ...
        'allTrialsDSLfilteredSharedVarPct', 'stimDirDSLfilteredSharedVarPct', ...
        'stimNameDirDSLfilteredSharedVarPct', 'conditionDSLfilteredSharedVarPct', ...
        'DSLfilteredSharedVarPct'};

    for g = 1:numel(figStruct)
        for f = 1:numel(fieldsToSave)
            fieldName = fieldsToSave{f};
            if isfield(figStruct(g), fieldName)
                saveAndCloseFigure(figStruct(g).(fieldName), folderPath, sprintf(fieldPatterns.(fieldName), g));
            end
        end

        % Close any extra figures generated by the same struct but not saved under this config.
        for f = 1:numel(allFigureFields)
            fieldName = allFigureFields{f};
            if isfield(figStruct(g), fieldName)
                figHandle = figStruct(g).(fieldName);
                if ~isempty(figHandle) && isgraphics(figHandle)
                    close(figHandle);
                end
            end
        end
    end
end


function [SummaryDSL, figHandles] = summarizeAllConditionsDSL(AllConditionResults, DSL_threshold)

    if isempty(AllConditionResults)
        error('AllConditionResults is empty.');
    end

    numGroups = numel(AllConditionResults(1).DSL.indiv);

    SummaryDSL = struct();
    SummaryDSL.meta.numConditions = numel(AllConditionResults);
    SummaryDSL.meta.conditions = [AllConditionResults.condition];
    SummaryDSL.meta.threshold = DSL_threshold;

    figHandles = gobjects(1, numGroups);

    for g = 1:numGroups
        mergedValues = [];
        conditionMembership = [];

        for ci = 1:numel(AllConditionResults)
            vals = AllConditionResults(ci).DSL.indiv{g};
            vals = vals(isfinite(vals));

            mergedValues = [mergedValues, vals];
            conditionMembership = [conditionMembership, repmat(AllConditionResults(ci).condition, 1, numel(vals))];
        end

        [figHandles(g), histStat] = plotDSLHistogram(mergedValues, DSL_threshold, ...
            sprintf('Group %d DSL distribution (sum all conditions)', g));

        SummaryDSL.group(g).values = mergedValues;
        SummaryDSL.group(g).conditionMembership = conditionMembership;
        SummaryDSL.group(g).histogram = histStat;
    end
end


function [SummaryCategory, figHandles] = summarizeAllConditionsCategories(AllConditionResults, condition_list, condition_full)

    if isempty(AllConditionResults)
        error('AllConditionResults is empty.');
    end

    if nargin < 3 || isempty(condition_full)
        error('condition_full is required to decode condition IDs into stimulus labels.');
    end

    numConditions = numel(AllConditionResults);
    numGroups = numel(AllConditionResults(1).Stats.group);
    labels = AllConditionResults(1).Stats.labels;
    numCategories = numel(labels);

    conditionMap = buildConditionSummaryMap(condition_full, condition_list);

    SummaryCategory = struct();
    SummaryCategory.labels = labels;
    SummaryCategory.meta.numConditions = numConditions;
    SummaryCategory.meta.conditions = condition_list;
    SummaryCategory.meta.conditionMap = conditionMap;
    SummaryCategory.meta.panelConditionLabels = { ...
        'grating-small-low'
        'grating-small-high'
        'grating-large-low'
        'grating-large-high'
        'plaid-small-low'
        'plaid-small-high'
        'plaid-large-low'
        'plaid-large-high'};
    SummaryCategory.meta.panelConditionShortLabels = { ...
        'G-S-L'
        'G-S-H'
        'G-L-L'
        'G-L-H'
        'P-S-L'
        'P-S-H'
        'P-L-L'
        'P-L-H'};
    SummaryCategory.meta.stimDirLabels = conditionMap.meta.stimDirLabels;
    SummaryCategory.meta.stimDirValues = conditionMap.meta.stimDirValues;
    SummaryCategory.meta.stimDirPanelTitles = { ...
        sprintf('%s = %s', conditionMap.meta.stimDirLabels{1}, formatSummaryValue(conditionMap.meta.stimDirValues(1))), ...
        sprintf('%s = %s', conditionMap.meta.stimDirLabels{2}, formatSummaryValue(conditionMap.meta.stimDirValues(2)))};
    if isfield(AllConditionResults(1).Stats, 'meta') && isfield(AllConditionResults(1).Stats.meta, 'filteredDisplayName')
        SummaryCategory.meta.filteredDisplayName = AllConditionResults(1).Stats.meta.filteredDisplayName;
    else
        SummaryCategory.meta.filteredDisplayName = 'DSL filtered';
    end
    SummaryCategory.meta.filteredFileTag = makeFilterModeFileTag(SummaryCategory.meta.filteredDisplayName);

    figHandles = struct([]);

    for g = 1:numGroups

        all_dir1_percentLatents     = nan(8, numCategories);
        all_dir1_percentSharedVar   = nan(8, numCategories);
        all_dir1_counts             = nan(8, numCategories);
        all_dir1_sharedVariance     = nan(8, numCategories);
        all_dir1_conditionIds       = nan(8, 1);

        all_dir2_percentLatents     = nan(8, numCategories);
        all_dir2_percentSharedVar   = nan(8, numCategories);
        all_dir2_counts             = nan(8, numCategories);
        all_dir2_sharedVariance     = nan(8, numCategories);
        all_dir2_conditionIds       = nan(8, 1);

        filt_dir1_percentLatents    = nan(8, numCategories);
        filt_dir1_percentSharedVar  = nan(8, numCategories);
        filt_dir1_counts            = nan(8, numCategories);
        filt_dir1_sharedVariance    = nan(8, numCategories);
        filt_dir1_conditionIds      = nan(8, 1);

        filt_dir2_percentLatents    = nan(8, numCategories);
        filt_dir2_percentSharedVar  = nan(8, numCategories);
        filt_dir2_counts            = nan(8, numCategories);
        filt_dir2_sharedVariance    = nan(8, numCategories);
        filt_dir2_conditionIds      = nan(8, 1);

        for ci = 1:numel(AllConditionResults)
            condID = AllConditionResults(ci).condition;

            mapIdx = find([conditionMap.entries.conditionId] == condID, 1);
            if isempty(mapIdx)
                warning('Condition %d not found in condition_full mapping. Skipping.', condID);
                continue;
            end

            stimDirCode   = conditionMap.entries(mapIdx).stimDirCode;
            panelCondIdx  = conditionMap.entries(mapIdx).panelCondIndex;

            S_all      = AllConditionResults(ci).Stats.group(g).all;
            S_filtered = AllConditionResults(ci).Stats.group(g).filtered;

            if stimDirCode == 1
                all_dir1_percentLatents(panelCondIdx, :)    = S_all.percentLatents;
                all_dir1_percentSharedVar(panelCondIdx, :)  = S_all.percentSharedVariance;
                all_dir1_counts(panelCondIdx, :)            = S_all.counts;
                all_dir1_sharedVariance(panelCondIdx, :)    = S_all.sharedVariance;
                all_dir1_conditionIds(panelCondIdx)         = condID;

                filt_dir1_percentLatents(panelCondIdx, :)   = S_filtered.percentLatents;
                filt_dir1_percentSharedVar(panelCondIdx, :) = S_filtered.percentSharedVariance;
                filt_dir1_counts(panelCondIdx, :)           = S_filtered.counts;
                filt_dir1_sharedVariance(panelCondIdx, :)   = S_filtered.sharedVariance;
                filt_dir1_conditionIds(panelCondIdx)        = condID;
            else
                all_dir2_percentLatents(panelCondIdx, :)    = S_all.percentLatents;
                all_dir2_percentSharedVar(panelCondIdx, :)  = S_all.percentSharedVariance;
                all_dir2_counts(panelCondIdx, :)            = S_all.counts;
                all_dir2_sharedVariance(panelCondIdx, :)    = S_all.sharedVariance;
                all_dir2_conditionIds(panelCondIdx)         = condID;

                filt_dir2_percentLatents(panelCondIdx, :)   = S_filtered.percentLatents;
                filt_dir2_percentSharedVar(panelCondIdx, :) = S_filtered.percentSharedVariance;
                filt_dir2_counts(panelCondIdx, :)           = S_filtered.counts;
                filt_dir2_sharedVariance(panelCondIdx, :)   = S_filtered.sharedVariance;
                filt_dir2_conditionIds(panelCondIdx)        = condID;
            end
        end

        SummaryCategory.group(g).all.stim_dir(1).label = conditionMap.meta.stimDirLabels{1};
        SummaryCategory.group(g).all.stim_dir(1).value = conditionMap.meta.stimDirValues(1);
        SummaryCategory.group(g).all.stim_dir(1).conditionIds = all_dir1_conditionIds;
        SummaryCategory.group(g).all.stim_dir(1).percentLatents = all_dir1_percentLatents;
        SummaryCategory.group(g).all.stim_dir(1).percentSharedVariance = all_dir1_percentSharedVar;
        SummaryCategory.group(g).all.stim_dir(1).counts = all_dir1_counts;
        SummaryCategory.group(g).all.stim_dir(1).sharedVariance = all_dir1_sharedVariance;

        SummaryCategory.group(g).all.stim_dir(2).label = conditionMap.meta.stimDirLabels{2};
        SummaryCategory.group(g).all.stim_dir(2).value = conditionMap.meta.stimDirValues(2);
        SummaryCategory.group(g).all.stim_dir(2).conditionIds = all_dir2_conditionIds;
        SummaryCategory.group(g).all.stim_dir(2).percentLatents = all_dir2_percentLatents;
        SummaryCategory.group(g).all.stim_dir(2).percentSharedVariance = all_dir2_percentSharedVar;
        SummaryCategory.group(g).all.stim_dir(2).counts = all_dir2_counts;
        SummaryCategory.group(g).all.stim_dir(2).sharedVariance = all_dir2_sharedVariance;

        SummaryCategory.group(g).filtered.stim_dir(1).label = conditionMap.meta.stimDirLabels{1};
        SummaryCategory.group(g).filtered.stim_dir(1).value = conditionMap.meta.stimDirValues(1);
        SummaryCategory.group(g).filtered.stim_dir(1).conditionIds = filt_dir1_conditionIds;
        SummaryCategory.group(g).filtered.stim_dir(1).percentLatents = filt_dir1_percentLatents;
        SummaryCategory.group(g).filtered.stim_dir(1).percentSharedVariance = filt_dir1_percentSharedVar;
        SummaryCategory.group(g).filtered.stim_dir(1).counts = filt_dir1_counts;
        SummaryCategory.group(g).filtered.stim_dir(1).sharedVariance = filt_dir1_sharedVariance;

        SummaryCategory.group(g).filtered.stim_dir(2).label = conditionMap.meta.stimDirLabels{2};
        SummaryCategory.group(g).filtered.stim_dir(2).value = conditionMap.meta.stimDirValues(2);
        SummaryCategory.group(g).filtered.stim_dir(2).conditionIds = filt_dir2_conditionIds;
        SummaryCategory.group(g).filtered.stim_dir(2).percentLatents = filt_dir2_percentLatents;
        SummaryCategory.group(g).filtered.stim_dir(2).percentSharedVariance = filt_dir2_percentSharedVar;
        SummaryCategory.group(g).filtered.stim_dir(2).counts = filt_dir2_counts;
        SummaryCategory.group(g).filtered.stim_dir(2).sharedVariance = filt_dir2_sharedVariance;

        figHandles(g).allLatentPct = plotSplitDirGroupedBar( ...
            all_dir1_percentLatents, all_dir2_percentLatents, ...
            labels, SummaryCategory.meta.panelConditionShortLabels, ...
            SummaryCategory.meta.stimDirPanelTitles, ...
            sprintf('Group %d - Percentage of latents (sum all conditions, all latents)', g), ...
            'Latent category', ...
            'Percentage of latents (%)');

        figHandles(g).allSharedVarPct = plotSplitDirGroupedBar( ...
            all_dir1_percentSharedVar, all_dir2_percentSharedVar, ...
            labels, SummaryCategory.meta.panelConditionShortLabels, ...
            SummaryCategory.meta.stimDirPanelTitles, ...
            sprintf('Group %d - Percentage of shared variance explained (sum all conditions, all latents)', g), ...
            'Latent category', ...
            'Percentage of shared variance explained (%)');

        figHandles(g).filteredLatentPct = plotSplitDirGroupedBar( ...
            filt_dir1_percentLatents, filt_dir2_percentLatents, ...
            labels, SummaryCategory.meta.panelConditionShortLabels, ...
            SummaryCategory.meta.stimDirPanelTitles, ...
            sprintf('Group %d - Percentage of latents (sum all conditions, %s)', g, SummaryCategory.meta.filteredDisplayName), ...
            'Latent category', ...
            'Percentage of latents (%)');

        figHandles(g).filteredSharedVarPct = plotSplitDirGroupedBar( ...
            filt_dir1_percentSharedVar, filt_dir2_percentSharedVar, ...
            labels, SummaryCategory.meta.panelConditionShortLabels, ...
            SummaryCategory.meta.stimDirPanelTitles, ...
            sprintf('Group %d - Percentage of shared variance explained (sum all conditions, %s)', g, SummaryCategory.meta.filteredDisplayName), ...
            'Latent category', ...
            'Percentage of shared variance explained (%)');
    end
end

function conditionMap = buildConditionSummaryMap(condition_full, condition_list)

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

    % After-the-fact correction for old, non-canonical condition_full.
    % This makes equivalent directions such as 360/0 and -10/350 identical
    % for downstream grouping and plotting.
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
        error('Stim %s does not have exactly 2 contrast levels.', char(stimLabels(s)));
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
    'grating-small-low' 'grating-small-high' ...
    'grating-large-low' 'grating-large-high' ...
    'plaid-small-low' 'plaid-small-high' ...
    'plaid-large-low' 'plaid-large-high'};

condShortLabels = { ...
    'G-S-L' 'G-S-H' 'G-L-L' 'G-L-H' ...
    'P-S-L' 'P-S-H' 'P-L-L' 'P-L-H'};

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
        error('Could not map size value %s to size code.', formatSummaryValue(currSize));
    end
    if isempty(contrastCode)
        error('Could not map contrast value %s to contrast code for stim %s.', ...
            formatSummaryValue(currContrast), char(currStim));
    end
    if isempty(stimDirCode)
        error('Could not map canonical direction value %s to stim_dir code.', ...
            formatSummaryValue(currDir));
    end

    panelCondIndex = (stimCode - 1) * 4 + (sizeCode - 1) * 2 + contrastCode;

    entries(ii).conditionId = condID;
    entries(ii).stimName = char(currStim);
    entries(ii).stimCode = stimCode;

    entries(ii).sizeValue = currSize;
    entries(ii).sizeCode = sizeCode;
    entries(ii).sizeLabel = ternary_label(sizeCode, 'small', 'large');

    entries(ii).contrastValue = currContrast;
    entries(ii).contrastCode = contrastCode;
    entries(ii).contrastLabel = ternary_label(contrastCode, 'low', 'high');

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

function s = formatSummaryValue(v)
    if ~isfinite(v)
        s = 'NaN';
    elseif abs(v - round(v)) < 1e-10
        s = sprintf('%d', round(v));
    else
        s = sprintf('%.4g', v);
    end
end


function tag = makeFilterModeFileTag(label)
    label = char(string(label));
    tag = lower(label);
    tag = regexprep(tag, '[^a-z0-9]+', '_');
    tag = regexprep(tag, '^_+|_+$', '');
end


function out = ternary_label(code, label1, label2)
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


function [CondSV, figHandles] = analyzeConditionSpecificPosteriorVarExpFromPooled( ...
    modelInput, model_varexp_indiv, seqEst, condition_full, trial_condition_ids, ...
    gp_params, ambiguousIdxs, DSL, makePlots)
% analyzeConditionSpecificPosteriorVarExpFromPooled
%
% Posterior-based condition-specific shared variance explained for a pooled
% all-condition DLAG model.
%
% This function accepts either:
%   1) res, where the true DLAG params are in res.estParams
%   2) a params struct directly (for example res.estParams)
%
% Core definition for group m, latent j, condition c:
%
%   rawSV_j(c) = ||c_j||^2 * q_j(c)
%
% where
%
%   q_j(c) = mean_{trials,time in condition c} [ xsm_j^2 + Var(x_j|y) ]
%
% Percent values stored in this function are true percentages in [0,100].
%
% Note:
%   Single-latent figures are intentionally not generated. The latent-level
%   data are still stored in CondSV for downstream analysis.

    if nargin < 9
        makePlots = true;
    end

    params = resolvePosteriorSVParams(modelInput);

    if ~isstruct(params) || ~isfield(params, 'C') || ...
            ~isfield(params, 'xDim_across') || ~isfield(params, 'xDim_within') || ...
            ~isfield(params, 'yDims')
        error('DLAG params must contain C, xDim_across, xDim_within, and yDims.');
    end

    if isempty(seqEst)
        error('seqEst is empty.');
    end

    if numel(seqEst) ~= numel(trial_condition_ids)
        error('numel(seqEst) must match numel(trial_condition_ids).');
    end

    xDim_across = params.xDim_across;
    xDim_within = params.xDim_within(:)';
    yDims = params.yDims(:)';

    localDims = xDim_across + xDim_within;
    totalLatents = sum(localDims);
    numGroups = numel(xDim_within);

    trialLengths = arrayfun(@(s) size(s.xsm, 2), seqEst);
    if any(trialLengths ~= trialLengths(1))
        error('All trials must have the same T for this posterior-based pooled analysis.');
    end
    T = trialLengths(1);

    if ~isfield(seqEst, 'Vsm')
        error('seqEst must contain field Vsm.');
    end

    % -------------------------------------------------------------
    % Trial/condition mapping
    % -------------------------------------------------------------
    condition_ids_all = 1:numel(condition_full);
    conditionMap = buildConditionSummaryMap(condition_full, condition_ids_all);

    panelConditionLabels = conditionMap.meta.panelCondLabels;
    panelConditionShortLabels = conditionMap.meta.panelCondShortLabels;
    stimDirLabels = conditionMap.meta.stimDirLabels;
    stimDirValues = conditionMap.meta.stimDirValues;

    condBarX = [1 2 4 5 8 9 11 12];
    condColors = condsv_getConditionColors(numel(panelConditionShortLabels));

    % -------------------------------------------------------------
    % Latent category classification
    % -------------------------------------------------------------
    latentClass = classifyDlagLatents(xDim_across, gp_params, ambiguousIdxs);
    acrossDelay = latentClass.acrossDelay;
    ffIdx = latentClass.feedforwardIdx;
    fbIdx = latentClass.feedbackIdx;
    ambiguousAll = latentClass.ambiguousIdx;
    categoryLabels = latentClass.categoryLabels;

    % -------------------------------------------------------------
    % Block indices for group rows/cols
    % -------------------------------------------------------------
    obsBlocks = condsv_makeBlockIndices(yDims);
    latBlocks = condsv_makeBlockIndices(localDims);

    % -------------------------------------------------------------
    % Posterior variance diagonal mean across time
    % Since all T are equal, Vsm is the same across same-length trials.
    % Use the first trial.
    % -------------------------------------------------------------
    Vsm0 = seqEst(1).Vsm;
    if ~isequal(size(Vsm0), [totalLatents, totalLatents, T])
        error('seqEst(1).Vsm has unexpected size.');
    end

    posteriorVarMean = nan(1, totalLatents);
    for r = 1:totalLatents
        posteriorVarMean(r) = mean(squeeze(Vsm0(r, r, :)), 'omitnan');
    end

    % -------------------------------------------------------------
    % Initialize output
    % -------------------------------------------------------------
    CondSV = struct();
    CondSV.meta.definition = ['rawSharedVariance_j(c) = ||c_j||^2 * mean(xsm_j^2 + posteriorVar_j) ' ...
                              'within condition c'];
    CondSV.meta.usesPosteriorSecondMoment = true;
    CondSV.meta.numGroups = numGroups;
    CondSV.meta.xDim_across = xDim_across;
    CondSV.meta.xDim_within = xDim_within;
    CondSV.meta.yDims = yDims;
    CondSV.meta.localDims = localDims;
    CondSV.meta.conditionMap = conditionMap;
    CondSV.meta.conditionLabels = panelConditionLabels;
    CondSV.meta.conditionShortLabels = panelConditionShortLabels;
    CondSV.meta.conditionBarX = condBarX;
    CondSV.meta.conditionColors = condColors;
    CondSV.meta.stimDirLabels = stimDirLabels;
    CondSV.meta.stimDirValues = stimDirValues;
    CondSV.meta.categoryLabels = categoryLabels;
    CondSV.meta.filterModeNames = {'all', 'filtered_DSL_alltrials', 'filtered_DSL_bystimdir', 'filtered_DSL_bystimnamedir', 'filtered_DSL_bycondition'};
    CondSV.meta.filterModeDisplayNames = {'all latents', 'all-trials DSL filtered', 'stim_dir DSL filtered', 'stim_name x stim_dir DSL filtered', 'condition DSL filtered'};
    CondSV.meta.acrossDelay = acrossDelay;
    CondSV.meta.feedforwardIdx = ffIdx;
    CondSV.meta.feedbackIdx = fbIdx;
    CondSV.meta.ambiguousIdx = ambiguousAll;

    figHandles = struct([]);

    % -------------------------------------------------------------
    % Main computation, group by group
    % -------------------------------------------------------------
    for g = 1:numGroups

        obsIdx = obsBlocks{g}(1):obsBlocks{g}(2);
        latIdx = latBlocks{g}(1):latBlocks{g}(2);

        localDim = localDims(g);

        Cg = params.C(obsIdx, latIdx);
        cNormSq = sum(Cg.^2, 1);

        % Collect all posterior means for this group
        Xg = nan(numel(seqEst), localDim, T);
        for n = 1:numel(seqEst)
            Xg(n, :, :) = seqEst(n).xsm(latIdx, :);
        end

        % Pooled posterior benchmark
        pooledQ = nan(1, localDim);
        for j = 1:localDim
            xj = squeeze(Xg(:, j, :));
            pooledQ(j) = mean(xj(:).^2, 'omitnan') + posteriorVarMean(latIdx(j));
        end

        pooledRaw = cNormSq .* pooledQ;
        pooledFrac = pooledRaw ./ sum(pooledRaw, 'omitnan');
        pooledPct = 100 * pooledFrac;

        CondSV.group(g).name = sprintf('Group %d', g);
        CondSV.group(g).pooledPosterior.loadingNormSq = cNormSq;
        CondSV.group(g).pooledPosterior.posteriorSecondMoment = pooledQ;
        CondSV.group(g).pooledPosterior.rawSharedVariance = pooledRaw;
        CondSV.group(g).pooledPosterior.fractionSharedVariance = pooledFrac;
        CondSV.group(g).pooledPosterior.percentSharedVariance = pooledPct;

        CondSV.group(g).pooledPosterior.modelVarExpFraction = model_varexp_indiv{g};
        CondSV.group(g).pooledPosterior.modelVarExpPercent = 100 * model_varexp_indiv{g};

        CondSV.group(g).pooledPosterior.modelDifferenceFraction = ...
            pooledFrac - model_varexp_indiv{g};
        CondSV.group(g).pooledPosterior.modelDifferencePercent = ...
            pooledPct - 100 * model_varexp_indiv{g};

        % Category masks and DSL keep masks
        categoryMasks = makeDlagCategoryMasks(localDim, xDim_across, ffIdx, fbIdx, ambiguousAll);

        keepAll = true(1, localDim);
        keepDSL = getDslKeepMask(DSL, 'logical', g, localDim, true);
        keepDSLByStimDir = getDslKeepMask(DSL, 'logical_bystimdir', g, localDim, true);
        keepDSLByStimNameDir = getDslKeepMask(DSL, 'logical_bystimnamedir', g, localDim, true);
        keepDSLByCondition = getDslKeepMask(DSL, 'logical_bycondition', g, localDim, true);

        % ---------------------------------------------------------
        % Per-condition latent raw / fraction / percent
        % ---------------------------------------------------------
        numAllConditions = numel(condition_full);
        rawLatByCondition = nan(numAllConditions, localDim);
        fracLatByCondition = nan(numAllConditions, localDim);
        pctLatByCondition = nan(numAllConditions, localDim);
        totalRawByCondition = nan(numAllConditions, 1);

        for condID = 1:numAllConditions
            trIdx = find(trial_condition_ids == condID);

            if isempty(trIdx)
                continue;
            end

            rawLat = nan(1, localDim);

            for j = 1:localDim
                xj = squeeze(Xg(trIdx, j, :));
                qj = mean(xj(:).^2, 'omitnan') + posteriorVarMean(latIdx(j));
                rawLat(j) = cNormSq(j) * qj;
            end

            rawLatByCondition(condID, :) = rawLat;

            denomAll = sum(rawLat(keepAll), 'omitnan');
            totalRawByCondition(condID) = denomAll;

            if denomAll > 0
                fracLatByCondition(condID, :) = rawLat / denomAll;
                pctLatByCondition(condID, :) = 100 * rawLat / denomAll;
            end
        end

        % ---------------------------------------------------------
        % Store latent-level information
        % Single-latent figures are not generated, but all latent data
        % remain available in CondSV.group(g).latent(l)
        % ---------------------------------------------------------
        CondSV.group(g).latent = struct([]);
        for l = 1:localDim
            latentInfo = condsv_makeLatentInfo( ...
                g, l, xDim_across, ffIdx, fbIdx, ambiguousAll, ...
                DSL.logical{g}(l), keepDSLByStimDir(l), keepDSLByStimNameDir(l), keepDSLByCondition(l));

            CondSV.group(g).latent(l).groupIndex = g;
            CondSV.group(g).latent(l).localLatentIndex = l;
            CondSV.group(g).latent(l).rowIndexInXsm = latIdx(l);

            CondSV.group(g).latent(l).latentType = latentInfo.latentType;
            CondSV.group(g).latent(l).acrossIndex = latentInfo.acrossIndex;
            CondSV.group(g).latent(l).withinIndex = latentInfo.withinIndex;
            CondSV.group(g).latent(l).acrossCategory = latentInfo.acrossCategory;
            CondSV.group(g).latent(l).titleLines = latentInfo.titleLines;

            CondSV.group(g).latent(l).DSL.logical = DSL.logical{g}(l);
            CondSV.group(g).latent(l).DSL.label = latentInfo.dslLabel;

            CondSV.group(g).latent(l).DSL_bystimdir.logical = keepDSLByStimDir(l);
            CondSV.group(g).latent(l).DSL_bystimdir.label = latentInfo.dslByStimDirLabel;

            CondSV.group(g).latent(l).DSL_bystimnamedir.logical = keepDSLByStimNameDir(l);
            CondSV.group(g).latent(l).DSL_bystimnamedir.label = latentInfo.dslByStimNameDirLabel;

            CondSV.group(g).latent(l).DSL_bycondition.logical = keepDSLByCondition(l);
            CondSV.group(g).latent(l).DSL_bycondition.label = latentInfo.dslByConditionLabel;

            CondSV.group(g).latent(l).allByConditionId.rawSharedVariance = rawLatByCondition(:, l);
            CondSV.group(g).latent(l).allByConditionId.fractionSharedVariance = fracLatByCondition(:, l);
            CondSV.group(g).latent(l).allByConditionId.percentSharedVariance = pctLatByCondition(:, l);

            for d = 1:2
                CondSV.group(g).latent(l).stim_dir(d).label = stimDirLabels{d};
                CondSV.group(g).latent(l).stim_dir(d).value = stimDirValues(d);
                CondSV.group(g).latent(l).stim_dir(d).conditionIds = nan(8, 1);
                CondSV.group(g).latent(l).stim_dir(d).rawSharedVariance = nan(8, 1);
                CondSV.group(g).latent(l).stim_dir(d).fractionSharedVariance = nan(8, 1);
                CondSV.group(g).latent(l).stim_dir(d).percentSharedVariance = nan(8, 1);
            end
        end

        % Map latent values into stim_dir panels
        for ii = 1:numel(conditionMap.entries)
            condID = conditionMap.entries(ii).conditionId;
            d = conditionMap.entries(ii).stimDirCode;
            p = conditionMap.entries(ii).panelCondIndex;

            for l = 1:localDim
                CondSV.group(g).latent(l).stim_dir(d).conditionIds(p) = condID;
                CondSV.group(g).latent(l).stim_dir(d).rawSharedVariance(p) = rawLatByCondition(condID, l);
                CondSV.group(g).latent(l).stim_dir(d).fractionSharedVariance(p) = fracLatByCondition(condID, l);
                CondSV.group(g).latent(l).stim_dir(d).percentSharedVariance(p) = pctLatByCondition(condID, l);
            end
        end

        % ---------------------------------------------------------
        % Category-level summaries under 5 filter modes
        % ---------------------------------------------------------
        modeConfigs = struct( ...
            'name', {'all', 'filtered_DSL_alltrials', 'filtered_DSL_bystimdir', 'filtered_DSL_bystimnamedir', 'filtered_DSL_bycondition'}, ...
            'keepMask', {keepAll, keepDSL, keepDSLByStimDir, keepDSLByStimNameDir, keepDSLByCondition}, ...
            'displayName', {'all latents', 'all-trials DSL filtered', 'stim_dir DSL filtered', 'stim_name x stim_dir DSL filtered', 'condition DSL filtered'}, ...
            'figField', {'allSharedVarPct', 'allTrialsDSLfilteredSharedVarPct', 'stimDirDSLfilteredSharedVarPct', 'stimNameDirDSLfilteredSharedVarPct', 'conditionDSLfilteredSharedVarPct'});

        for mi = 1:numel(modeConfigs)
            modeName = modeConfigs(mi).name;
            keepMask = modeConfigs(mi).keepMask;

            catSummary = computeCategorySummaryFromRaw(rawLatByCondition, keepMask, categoryMasks, categoryLabels);

            catRawByCondition = catSummary.rawByConditionId;
            catFracByCondition = catSummary.fractionByConditionId;
            catPctByCondition = catSummary.percentByConditionId;
            catTotalByCondition = catSummary.totalByConditionId;

            CondSV.group(g).category.(modeName).displayName = modeConfigs(mi).displayName;
            CondSV.group(g).category.(modeName).figureFieldName = modeConfigs(mi).figField;
            CondSV.group(g).category.(modeName).rawByConditionId = catRawByCondition;
            CondSV.group(g).category.(modeName).fractionByConditionId = catFracByCondition;
            CondSV.group(g).category.(modeName).percentByConditionId = catPctByCondition;
            CondSV.group(g).category.(modeName).totalByConditionId = catTotalByCondition;

            for d = 1:2
                CondSV.group(g).category.(modeName).stim_dir(d).label = stimDirLabels{d};
                CondSV.group(g).category.(modeName).stim_dir(d).value = stimDirValues(d);
                CondSV.group(g).category.(modeName).stim_dir(d).conditionIds = nan(8, 1);
                CondSV.group(g).category.(modeName).stim_dir(d).sharedVariance = nan(8, numel(categoryLabels));
                CondSV.group(g).category.(modeName).stim_dir(d).fractionSharedVariance = nan(8, numel(categoryLabels));
                CondSV.group(g).category.(modeName).stim_dir(d).percentSharedVariance = nan(8, numel(categoryLabels));
                CondSV.group(g).category.(modeName).stim_dir(d).totalSharedVariance = nan(8, 1);
            end

            for ii = 1:numel(conditionMap.entries)
                condID = conditionMap.entries(ii).conditionId;
                d = conditionMap.entries(ii).stimDirCode;
                p = conditionMap.entries(ii).panelCondIndex;

                CondSV.group(g).category.(modeName).stim_dir(d).conditionIds(p) = condID;
                CondSV.group(g).category.(modeName).stim_dir(d).sharedVariance(p, :) = catRawByCondition(condID, :);
                CondSV.group(g).category.(modeName).stim_dir(d).fractionSharedVariance(p, :) = catFracByCondition(condID, :);
                CondSV.group(g).category.(modeName).stim_dir(d).percentSharedVariance(p, :) = catPctByCondition(condID, :);
                CondSV.group(g).category.(modeName).stim_dir(d).totalSharedVariance(p) = catTotalByCondition(condID);
            end
        end
    end

    % -------------------------------------------------------------
    % Plot category-level figures only
    % -------------------------------------------------------------
    if makePlots
        plotModeConfigs = struct( ...
            'fieldName', {'allSharedVarPct', 'allTrialsDSLfilteredSharedVarPct', 'stimDirDSLfilteredSharedVarPct', 'stimNameDirDSLfilteredSharedVarPct', 'conditionDSLfilteredSharedVarPct'}, ...
            'modeName', {'all', 'filtered_DSL_alltrials', 'filtered_DSL_bystimdir', 'filtered_DSL_bystimnamedir', 'filtered_DSL_bycondition'}, ...
            'displayName', {'all latents', 'all-trials DSL filtered', 'stim_dir DSL filtered', 'stim_name x stim_dir DSL filtered', 'condition DSL filtered'});

        for g = 1:numGroups
            for mi = 1:numel(plotModeConfigs)
                modeName = plotModeConfigs(mi).modeName;
                fieldName = plotModeConfigs(mi).fieldName;
                displayName = plotModeConfigs(mi).displayName;

                figHandles(g).(fieldName) = plotSplitDirGroupedBar( ...
                    CondSV.group(g).category.(modeName).stim_dir(1).percentSharedVariance, ...
                    CondSV.group(g).category.(modeName).stim_dir(2).percentSharedVariance, ...
                    categoryLabels, panelConditionShortLabels, ...
                    {sprintf('%s = %s', stimDirLabels{1}, formatSummaryValue(stimDirValues(1))), ...
                     sprintf('%s = %s', stimDirLabels{2}, formatSummaryValue(stimDirValues(2)))}, ...
                    sprintf('Group %d - Condition-specific shared variance explained (%s)', g, displayName), ...
                    'Latent category', ...
                    'Percentage of shared variance explained (%)', ...
                    condColors);
            end

            % backward-compatible alias
            figHandles(g).DSLfilteredSharedVarPct = figHandles(g).allTrialsDSLfilteredSharedVarPct;
        end
    end
end


function figHandle = plotSplitDirGroupedBar(valuesDir1, valuesDir2, labels, legendLabels, panelTitles, figTitle, xlab, ylab, condColors)

    if nargin < 9 || isempty(condColors)
        condColors = condsv_getConditionColors(numel(legendLabels));
    end

    figHandle = figure('Position', [100 100 1500 820]);
    tl = tiledlayout(figHandle, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax1 = nexttile(tl, 1);
    b1 = bar(ax1, valuesDir1');
    condsv_setGroupedBarColors(b1, condColors);
    set(ax1, 'XTick', 1:numel(labels), 'XTickLabel', labels);
    xtickangle(ax1, 30);
    xlabel(ax1, xlab);
    ylabel(ax1, ylab);
    title(ax1, panelTitles{1});

    ax2 = nexttile(tl, 2);
    b2 = bar(ax2, valuesDir2');
    condsv_setGroupedBarColors(b2, condColors);
    set(ax2, 'XTick', 1:numel(labels), 'XTickLabel', labels);
    xtickangle(ax2, 30);
    xlabel(ax2, xlab);
    ylabel(ax2, ylab);
    title(ax2, panelTitles{2});

    ymax = condsv_computeFlexibleYMax([valuesDir1(:); valuesDir2(:)]);
    ylim(ax1, [0 ymax]);
    ylim(ax2, [0 ymax]);

    lgd = legend(ax1, b1, legendLabels, ...
        'Location', 'southoutside', ...
        'Orientation', 'horizontal');
    lgd.NumColumns = 4;
    lgd.Layout.Tile = 'south';

    sgtitle(tl, figTitle);
end


function params = resolvePosteriorSVParams(modelInput)
% resolvePosteriorSVParams
%
% Accept either:
%   - res struct returned by getModel_dlag, where params are in res.estParams
%   - params struct directly

    if ~isstruct(modelInput)
        error('Input must be a struct: either res or params.');
    end

    if isfield(modelInput, 'estParams')
        params = modelInput.estParams;
        return;
    end

    if isfield(modelInput, 'C') && ...
       isfield(modelInput, 'xDim_across') && ...
       isfield(modelInput, 'xDim_within') && ...
       isfield(modelInput, 'yDims')
        params = modelInput;
        return;
    end

    error(['Could not resolve DLAG params. Pass either res or res.estParams. ', ...
           'In your saved bestmodel* file, the correct object is res.estParams.']);
end


function colors = condsv_getConditionColors(nCond)
% Use one stable color set for all condition-wise plots.
    if nargin < 1 || isempty(nCond)
        nCond = 8;
    end
    colors = lines(nCond);
end


function condsv_setGroupedBarColors(barHandles, colors)
% Apply one consistent color per condition series.
    nBar = min(numel(barHandles), size(colors, 1));
    for k = 1:nBar
        barHandles(k).FaceColor = colors(k, :);
        barHandles(k).EdgeColor = 'none';
    end
end


function ymax = condsv_computeFlexibleYMax(values)
% Adaptive y-axis upper bound for percentage plots.
% This avoids forcing everything to 100 when the actual values are small.

    values = values(isfinite(values));

    if isempty(values)
        ymax = 1;
        return;
    end

    maxVal = max(values);

    if maxVal <= 0
        ymax = 1;
    elseif maxVal <= 0.5
        ymax = 0.6;
    elseif maxVal <= 1
        ymax = 1.2;
    elseif maxVal <= 2
        ymax = 2.4;
    elseif maxVal <= 5
        ymax = 6;
    elseif maxVal <= 10
        ymax = 12;
    elseif maxVal <= 20
        ymax = 24;
    elseif maxVal <= 40
        ymax = 45;
    elseif maxVal <= 60
        ymax = 66;
    elseif maxVal <= 80
        ymax = 88;
    else
        ymax = maxVal * 1.08;
    end
end


function latentInfo = condsv_makeLatentInfo(groupIdx, localIdx, xDim_across, ffIdx, fbIdx, ambiguousIdxs, dslLogical, dslByStimDirLogical, dslByStimNameDirLogical, dslByConditionLogical)

    latentInfo = struct();

    if localIdx <= xDim_across
        latentInfo.latentType = 'across';
        latentInfo.acrossIndex = localIdx;
        latentInfo.withinIndex = [];

        if ismember(localIdx, ambiguousIdxs)
            latentInfo.acrossCategory = 'ambiguous';
        elseif ismember(localIdx, ffIdx)
            latentInfo.acrossCategory = 'feedforward';
        elseif ismember(localIdx, fbIdx)
            latentInfo.acrossCategory = 'feedback';
        else
            latentInfo.acrossCategory = 'ambiguous';
        end

        latentLine = sprintf('Across latent %d (%s)', localIdx, latentInfo.acrossCategory);
    else
        latentInfo.latentType = 'within';
        latentInfo.acrossIndex = [];
        latentInfo.withinIndex = localIdx - xDim_across;
        latentInfo.acrossCategory = '';
        latentLine = sprintf('Within latent %d', latentInfo.withinIndex);
    end

    if dslLogical ~= 0
        dslLabel = 'DSL keep';
    else
        dslLabel = 'DSL remove';
    end

    if dslByStimDirLogical ~= 0
        dslByStimDirLabel = 'DSL(by stim_dir) keep';
    else
        dslByStimDirLabel = 'DSL(by stim_dir) remove';
    end

    if dslByStimNameDirLogical ~= 0
        dslByStimNameDirLabel = 'DSL(by stim_name x stim_dir) keep';
    else
        dslByStimNameDirLabel = 'DSL(by stim_name x stim_dir) remove';
    end

    if dslByConditionLogical ~= 0
        dslByConditionLabel = 'DSL(by condition) keep';
    else
        dslByConditionLabel = 'DSL(by condition) remove';
    end

    latentInfo.dslLabel = dslLabel;
    latentInfo.dslByStimDirLabel = dslByStimDirLabel;
    latentInfo.dslByStimNameDirLabel = dslByStimNameDirLabel;
    latentInfo.dslByConditionLabel = dslByConditionLabel;
    latentInfo.latentLine = latentLine;
    latentInfo.titleLines = { ...
        sprintf('Group %d', groupIdx), ...
        latentLine, ...
        dslLabel, ...
        dslByStimDirLabel, ...
        dslByStimNameDirLabel, ...
        dslByConditionLabel};
end


function blocks = condsv_makeBlockIndices(dims)
    dims = dims(:)';
    blocks = cell(1, numel(dims));

    startIdx = 1;
    for k = 1:numel(dims)
        blocks{k} = [startIdx, startIdx + dims(k) - 1];
        startIdx = startIdx + dims(k);
    end
end

function d = getConditionEffectiveDirCanonicalLocal(cond, condID)
% Return canonical effective stimulus direction for one condition.
%
% This is an after-the-fact correction for model_data generated before
% angle normalization was added upstream. It makes equivalent angles such as
% 360/0 and -10/350 identical during downstream analysis.

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
% Convert degree angles to canonical [0, 360).
%
% Examples:
%   360  -> 0
%   -10  -> 350
%   720  -> 0
%
% NaN remains NaN.

a = double(a);
finiteMask = isfinite(a);

a(finiteMask) = mod(a(finiteMask), 360);

tol = 1e-10;

nearInteger = finiteMask & abs(a - round(a)) < tol;
a(nearInteger) = round(a(nearInteger));

a(finiteMask & abs(a) < tol) = 0;
a(finiteMask & abs(a - 360) < tol) = 0;

end

