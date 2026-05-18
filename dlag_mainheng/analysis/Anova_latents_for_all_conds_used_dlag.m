%% for Dlag all cond modes, latents para are same across conditions,
%% so we can compare whether latent response diff across diff conditions
%% at this moment, dir main effect for latent seems unclear and not important,
%% although can align this with populational dir preference, but how processing for plaid dir?
%% Suppose we find some ways, but this main effect is not our main focus point.
%% plan to use 3 way (contrast size and stim name)anova to analysis each latent
%% use label ambiguous and DSL.logical / DSL.logical_bystimdir / DSL.logical_bystimnamedir / DSL.logical_bycondition to label latents
%% use model_data_allruns{}.condition_full to label trials condition identity

% ==================
clc; clear

% Synthetic data generated from a DLAG model
% dat_file = 'I:\np_data\RafiL001p0120_g1\catgt_RafiL001p0120_g1/model_data_allruns';
dat_file = '.\model_data_allruns';
fprintf('Reading from %s \n', dat_file);
load(dat_file);

stim_tag = '_2[Gpl2_2c_2sz_400_2_200isi]';
data_content = 'z_across_conditions';
% options:
% raw_count, raw_fr, z_within_trial, z_within_condition,
% z_across_conditions, demean_count_within_trial, demean_fr_within_trial, demean_pooledsd_within_condition

runIdx = 1;
baseDir = ['./FA_Dlag_', data_content];       % Base directory where results will be saved

% Extract all stim tags
all_run_tags = get_all_run_tags(model_data_allruns);
wanted_tag = stim_tag;

% Find the requested run by exact stim_tag match.
run_idx = find(strcmp(all_run_tags, wanted_tag));

if isempty(run_idx)
    error('Requested stim_tag not found %s', wanted_tag);
end

if numel(run_idx) > 1
    error('Duplicate stim_tag found %s', wanted_tag);
end

% -------------------------------------------------------------------------
% Get the selected run entry stim index
% -------------------------------------------------------------------------
condition_full = model_data_allruns{run_idx}.conditions_full;

% get DSL, bestmodel latent info, gp_params delay info, ambiguous info
tempfname = sprintf('%s/mat_results/run%03d', baseDir, runIdx);

load(fullfile(tempfname, 'DSL_and_latent_category_stats.mat'), ...
    'DSL', 'bestModel', 'gp_params', 'ambiguousIdxs');

% get latent response, xsm
files = dir(fullfile(tempfname, 'bestmodel*'));
filename = fullfile(tempfname, files(1).name);
load(filename, "seqEst")

Results = analyze_dlag_latents_by_condition( ...
    condition_full, seqEst, bestModel.xDim_across, bestModel.xDim_within, ...
    gp_params, ambiguousIdxs, DSL, tempfname);

Results_split_by_dir = analyze_dlag_latents_by_condition_split_by_dir( ...
    condition_full, seqEst, bestModel.xDim_across, bestModel.xDim_within, ...
    gp_params, ambiguousIdxs, DSL, tempfname);

function Results = analyze_dlag_latents_by_condition(condition_full, seqEst, xDim_across, xDim_within, gp_params, ambiguousIdxs, DSL, tempfname)
% analyze_dlag_latents_by_condition
%
% Analyze condition dependence of DLAG latent trajectories using:
%   1) Condition-averaged time courses (2 x 2 x 2 conditions)
%   2) Three-way ANOVA on trial-averaged latent responses
%
% Factors:
%   - stim_name : grating vs plaid
%   - size      : small vs large
%   - contrast  : low vs high
%
% Trials are grouped by stim_name, size, and contrast.
% Direction variables are ignored here.
%
% OUTPUT
%   Results : structure containing all condition grouping information,
%             condition-averaged time courses, ANOVA results, simple effects,
%             and figure paths.
%
% SAVED OUTPUTS
%   Under tempfname/latent_condition_analysis/
%       - timecourse figures (.fig and .png)
%       - ANOVA summary figures (.fig and .png)
%       - latent_condition_analysis_results.mat
%
% ASSUMPTIONS
%   - All trials in seqEst for this run have the same number of time bins.
%   - seqEst(n).trialId matches condition_full(k).trial_indices.
%   - bestModel contains xDim_across and xDim_within.
%   - latent keep/remove labels use:
%       * DSL.logical
%       * DSL.logical_bystimdir      (final AND result, not raw per-dir)
%       * DSL.logical_bystimnamedir  (final AND result across 4 stim_name x stim_dir groups)
%       * DSL.logical_bycondition    (final AND result across all conditions)

    alpha = 0.05;

    % ---------------------------------------------------------------------
    % Validate basic inputs
    % ---------------------------------------------------------------------
    numGroups = numel(xDim_within);
    localDims = xDim_across + xDim_within;

    if ~isstruct(DSL) || ~isfield(DSL, 'logical')
        error('DSL must contain field DSL.logical.');
    end
    if ~isfield(DSL, 'logical_bystimdir')
        error('DSL must contain field DSL.logical_bystimdir.');
    end
    if ~isfield(DSL, 'logical_bystimnamedir')
        error('DSL must contain field DSL.logical_bystimnamedir.');
    end
    if ~isfield(DSL, 'logical_bycondition')
        error('DSL must contain field DSL.logical_bycondition.');
    end

    if numel(DSL.logical) ~= numGroups
        error('DSL.logical must have one cell per group.');
    end
    if numel(DSL.logical_bystimdir) ~= numGroups
        error('DSL.logical_bystimdir must have one cell per group.');
    end
    if numel(DSL.logical_bystimnamedir) ~= numGroups
        error('DSL.logical_bystimnamedir must have one cell per group.');
    end
    if numel(DSL.logical_bycondition) ~= numGroups
        error('DSL.logical_bycondition must have one cell per group.');
    end

    if isempty(seqEst)
        error('seqEst is empty.');
    end

    trialLengths = arrayfun(@(s) size(s.xsm, 2), seqEst);
    if any(trialLengths ~= trialLengths(1))
        error('All seqEst trials must have the same number of time bins for this analysis.');
    end
    T = trialLengths(1);
    tAxis = 1:T;

    for g = 1:numGroups
        if numel(DSL.logical{g}) ~= localDims(g)
            error('DSL.logical{%d} length does not match xDim_across + xDim_within(%d).', g, g);
        end
        if numel(DSL.logical_bystimdir{g}) ~= localDims(g)
            error('DSL.logical_bystimdir{%d} length does not match xDim_across + xDim_within(%d).', g, g);
        end
        if numel(DSL.logical_bystimnamedir{g}) ~= localDims(g)
            error('DSL.logical_bystimnamedir{%d} length does not match xDim_across + xDim_within(%d).', g, g);
        end
        if numel(DSL.logical_bycondition{g}) ~= localDims(g)
            error('DSL.logical_bycondition{%d} length does not match xDim_across + xDim_within(%d).', g, g);
        end
    end

    % ---------------------------------------------------------------------
    % Output folders
    % ---------------------------------------------------------------------
    outDir   = fullfile(tempfname, 'latent_condition_analysis');
    timeDir  = fullfile(outDir, 'timecourses');
    anovaDir = fullfile(outDir, 'anova');
    ensure_dir(outDir);
    ensure_dir(timeDir);
    ensure_dir(anovaDir);

    % ---------------------------------------------------------------------
    % Trial metadata from condition_full
    % ---------------------------------------------------------------------
    trialMeta = build_trial_metadata(condition_full, seqEst);

    % ---------------------------------------------------------------------
    % Across-latent category labels
    % ---------------------------------------------------------------------
    acrossDelay = gp_params.delays;
    acrossCategory = classify_across_latents(acrossDelay, ambiguousIdxs, xDim_across);

    % ---------------------------------------------------------------------
    % Global latent row indexing in seqEst(n).xsm
    % ---------------------------------------------------------------------
    blockStart = cumsum([1, localDims(1:end-1)]);

    % ---------------------------------------------------------------------
    % Initialize results structure
    % ---------------------------------------------------------------------
    Results = struct();
    Results.meta.alpha = alpha;
    Results.meta.timeAxis = tAxis;
    Results.meta.numGroups = numGroups;
    Results.meta.xDim_across = xDim_across;
    Results.meta.xDim_within = xDim_within;
    Results.meta.localDims = localDims;
    Results.meta.outputDir = outDir;
    Results.meta.timeDir = timeDir;
    Results.meta.anovaDir = anovaDir;

    Results.meta.factorOrder = {'stim_name', 'size', 'contrast'};
    Results.meta.stimLevels = trialMeta.stimLabels;
    Results.meta.sizeLevels = {'small', 'large'};
    Results.meta.sizeValues = trialMeta.sizeValues;
    Results.meta.contrastLevels = {'low', 'high'};
    Results.meta.contrastValuesByStim = trialMeta.contrastValuesByStim;

    Results.meta.conditionLabels = trialMeta.condLabels;
    Results.meta.conditionShortLabels = trialMeta.condShortLabels;
    Results.meta.conditionBarX = trialMeta.barX;
    Results.meta.acrossDelay = acrossDelay;
    Results.meta.acrossCategory = acrossCategory;
    Results.meta.ambiguousIdxs = ambiguousIdxs(:)';

    Results.meta.DSLRules.allTrials = 'Use DSL.logical.';
    Results.meta.DSLRules.byStimDir = 'Use DSL.logical_bystimdir (final AND result; do not use raw per-dir values as latent labels).';
    Results.meta.DSLRules.byStimNameDir = 'Use DSL.logical_bystimnamedir (final AND result across stim_name x stim_dir groups).';
    Results.meta.DSLRules.byCondition = 'Use DSL.logical_bycondition (final AND result across all conditions).';

    Results.trialMeta = trialMeta;

    % ---------------------------------------------------------------------
    % Compute all per-group, per-latent condition summaries and ANOVAs
    % ---------------------------------------------------------------------
    for g = 1:numGroups
        Results.group(g).name = sprintf('Group %d', g);
        Results.group(g).latent = struct([]);

        for l = 1:localDims(g)
            rowIdx = blockStart(g) + l - 1;

            X = zeros(numel(seqEst), T);
            for n = 1:numel(seqEst)
                X(n, :) = seqEst(n).xsm(rowIdx, :);
            end

            trialMean = mean(X, 2);

            latentInfo = make_latent_info( ...
                g, l, xDim_across, acrossCategory, ...
                DSL.logical{g}(l), DSL.logical_bystimdir{g}(l), ...
                DSL.logical_bystimnamedir{g}(l), DSL.logical_bycondition{g}(l));

            [meanTime, semTime, respMean, respSEM, respN, condTrialSeqIdx, condTrialIds] = ...
                compute_condition_averages(X, trialMean, trialMeta.condIdx, trialMeta.trialId, 8);

            anovaRes = run_threeway_anova(trialMean, trialMeta.stimCode, trialMeta.sizeCode, trialMeta.contrastCode, alpha);
            mainEffects = compute_main_effect_summaries(trialMean, trialMeta, anovaRes.p);
            simpleEffects = compute_simple_effects(trialMean, trialMeta, alpha);

            Results.group(g).latent(l).groupIndex = g;
            Results.group(g).latent(l).localLatentIndex = l;
            Results.group(g).latent(l).rowIndexInXsm = rowIdx;

            Results.group(g).latent(l).latentType = latentInfo.latentType;
            Results.group(g).latent(l).acrossIndex = latentInfo.acrossIndex;
            Results.group(g).latent(l).withinIndex = latentInfo.withinIndex;
            Results.group(g).latent(l).acrossCategory = latentInfo.acrossCategory;

            Results.group(g).latent(l).DSL.logical = DSL.logical{g}(l);
            Results.group(g).latent(l).DSL.label = latentInfo.dslLabel;
            Results.group(g).latent(l).DSL.logical_bystimdir = DSL.logical_bystimdir{g}(l);
            Results.group(g).latent(l).DSL.label_bystimdir = latentInfo.dslByStimDirLabel;
            Results.group(g).latent(l).DSL.logical_bystimnamedir = DSL.logical_bystimnamedir{g}(l);
            Results.group(g).latent(l).DSL.label_bystimnamedir = latentInfo.dslByStimNameDirLabel;
            Results.group(g).latent(l).DSL.logical_bycondition = DSL.logical_bycondition{g}(l);
            Results.group(g).latent(l).DSL.label_bycondition = latentInfo.dslByConditionLabel;

            Results.group(g).latent(l).titleLines = latentInfo.titleLines;

            Results.group(g).latent(l).condition.labels = trialMeta.condLabels;
            Results.group(g).latent(l).condition.shortLabels = trialMeta.condShortLabels;
            Results.group(g).latent(l).condition.meanTime = meanTime;
            Results.group(g).latent(l).condition.semTime = semTime;
            Results.group(g).latent(l).condition.meanResponse = respMean;
            Results.group(g).latent(l).condition.semResponse = respSEM;
            Results.group(g).latent(l).condition.nTrials = respN;
            Results.group(g).latent(l).condition.trialSeqIdx = condTrialSeqIdx;
            Results.group(g).latent(l).condition.trialIds = condTrialIds;

            Results.group(g).latent(l).anova.responsePerTrial = trialMean;
            Results.group(g).latent(l).anova.p = anovaRes.p;
            Results.group(g).latent(l).anova.table = anovaRes.tbl;
            Results.group(g).latent(l).anova.stats = anovaRes.stats;
            Results.group(g).latent(l).anova.terms = anovaRes.terms;
            Results.group(g).latent(l).anova.termNames = anovaRes.termNames;
            Results.group(g).latent(l).anova.mainEffects = mainEffects;
            Results.group(g).latent(l).anova.simpleEffects = simpleEffects;
        end
    end

    % ---------------------------------------------------------------------
    % Time-course figures
    %   - Across latents: one figure per across latent, subplots by group
    %   - Within latents: one figure per group/latent
    % ---------------------------------------------------------------------
    for a = 1:xDim_across
        [figFile, pngFile] = plot_across_timecourse_figure(Results, a, timeDir);
        for g = 1:numGroups
            Results.group(g).latent(a).timecourse.figFile = figFile;
            Results.group(g).latent(a).timecourse.pngFile = pngFile;
        end
    end

    for g = 1:numGroups
        for w = 1:xDim_within(g)
            l = xDim_across + w;
            [figFile, pngFile] = plot_within_timecourse_figure(Results.group(g).latent(l), Results.meta.timeAxis, timeDir);
            Results.group(g).latent(l).timecourse.figFile = figFile;
            Results.group(g).latent(l).timecourse.pngFile = pngFile;
        end
    end

    % ---------------------------------------------------------------------
    % ANOVA summary figures
    % ---------------------------------------------------------------------
    for g = 1:numGroups
        for l = 1:localDims(g)
            [figFile, pngFile] = plot_anova_summary_figure(Results.group(g).latent(l), Results.meta, anovaDir, alpha);
            Results.group(g).latent(l).anova.figFile = figFile;
            Results.group(g).latent(l).anova.pngFile = pngFile;
        end
    end

    % ---------------------------------------------------------------------
    % Save results
    % ---------------------------------------------------------------------
    save(fullfile(outDir, 'latent_condition_analysis_results.mat'), 'Results', '-v7.3');
end

% =========================================================================
% Build per-trial factor labels from condition_full and seqEst.trialId
% =========================================================================
function trialMeta = build_trial_metadata(condition_full, seqEst)

    Ntr = numel(seqEst);

    if isfield(seqEst, 'trialId')
        trialId = [seqEst.trialId]';
    else
        trialId = (1:Ntr)';
    end

    stimName    = strings(Ntr, 1);
    sizeValue   = nan(Ntr, 1);
    contrastRaw = nan(Ntr, 1);

    for k = 1:numel(condition_full)
        if ~isfield(condition_full(k), 'trial_indices')
            error('condition_full(%d) is missing field trial_indices.', k);
        end
        theseIds = condition_full(k).trial_indices(:);
        tf = ismember(trialId, theseIds);
        stimName(tf) = lower(string(condition_full(k).stim_name));
        sizeValue(tf) = condition_full(k).size;
        contrastRaw(tf) = condition_full(k).contrast;
    end

    assigned = (stimName ~= "") & ~isnan(sizeValue) & ~isnan(contrastRaw);
    if any(~assigned)
        warning('Some seqEst trials could not be assigned to any condition_full entry. They will be excluded.');
    end

    % Stim order: grating first, plaid second when both exist
    allStim = unique(stimName(assigned), 'stable');
    allStim = lower(allStim);

    if all(ismember(["grating","plaid"], allStim))
        stimLabels = ["grating","plaid"];
    else
        if numel(allStim) ~= 2
            error('Expected exactly 2 stim levels after assignment.');
        end
        stimLabels = allStim(:)';
    end

    stimCode = nan(Ntr, 1);
    for s = 1:2
        stimCode(stimName == stimLabels(s)) = s;
    end

    % Size mapping: low numeric value = small, high numeric value = large
    sizeVals = unique(sizeValue(assigned));
    sizeVals = sort(sizeVals(:)');
    if numel(sizeVals) ~= 2
        error('Expected exactly 2 size values.');
    end
    sizeCode = nan(Ntr, 1);
    sizeCode(sizeValue == sizeVals(1)) = 1;
    sizeCode(sizeValue == sizeVals(2)) = 2;

    % Contrast mapping: low/high within each stim_name separately
    contrastCode = nan(Ntr, 1);
    contrastValuesByStim = struct();
    for s = 1:2
        idx = assigned & (stimCode == s);
        cvals = unique(contrastRaw(idx));
        cvals = sort(cvals(:)');
        if numel(cvals) ~= 2
            error('Stim %s does not have exactly 2 contrast levels.', stimLabels(s));
        end
        contrastCode(idx & contrastRaw == cvals(1)) = 1;
        contrastCode(idx & contrastRaw == cvals(2)) = 2;
        contrastValuesByStim.(char(stimLabels(s))) = cvals;
    end

    valid = assigned & isfinite(stimCode) & isfinite(sizeCode) & isfinite(contrastCode);

    condIdx = nan(Ntr, 1);
    condIdx(valid) = (stimCode(valid)-1) * 4 + (sizeCode(valid)-1) * 2 + contrastCode(valid);

    condLabels = { ...
        'grating-small-low'
        'grating-small-high'
        'grating-large-low'
        'grating-large-high'
        'plaid-small-low'
        'plaid-small-high'
        'plaid-large-low'
        'plaid-large-high'};

    condShortLabels = { ...
        'G-S-L'
        'G-S-H'
        'G-L-L'
        'G-L-H'
        'P-S-L'
        'P-S-H'
        'P-L-L'
        'P-L-H'};

    barX = [1 2 4 5 8 9 11 12];

    trialMeta = struct();
    trialMeta.trialId = trialId;
    trialMeta.stimName = cellstr(stimName);
    trialMeta.stimCode = stimCode;
    trialMeta.stimLabels = cellstr(stimLabels);
    trialMeta.sizeValue = sizeValue;
    trialMeta.sizeCode = sizeCode;
    trialMeta.sizeValues = sizeVals;
    trialMeta.contrastRaw = contrastRaw;
    trialMeta.contrastCode = contrastCode;
    trialMeta.contrastValuesByStim = contrastValuesByStim;
    trialMeta.valid = valid;
    trialMeta.condIdx = condIdx;
    trialMeta.condLabels = condLabels;
    trialMeta.condShortLabels = condShortLabels;
    trialMeta.barX = barX;
end

% =========================================================================
% Classify across latents into feedforward / feedback / ambiguous
% =========================================================================
function acrossCategory = classify_across_latents(acrossDelay, ambiguousIdxs, xDim_across)
    acrossCategory = repmat({''}, 1, xDim_across);

    ambiguousIdxs = unique(ambiguousIdxs(:)');
    ambiguousIdxs = ambiguousIdxs(ambiguousIdxs >= 1 & ambiguousIdxs <= xDim_across);

    for a = 1:xDim_across
        if ismember(a, ambiguousIdxs)
            acrossCategory{a} = 'ambiguous';
        else
            if isnan(acrossDelay(a)) || acrossDelay(a) == 0
                acrossCategory{a} = 'ambiguous';
            elseif acrossDelay(a) > 0
                acrossCategory{a} = 'feedforward';
            else
                acrossCategory{a} = 'feedback';
            end
        end
    end
end

% =========================================================================
% Build one latent info struct for titles and labels
% =========================================================================
function latentInfo = make_latent_info(groupIdx, localIdx, xDim_across, acrossCategory, dslLogical, dslLogicalByStimDir, dslLogicalByStimNameDir, dslLogicalByCondition)
    latentInfo = struct();

    if localIdx <= xDim_across
        latentInfo.latentType = 'across';
        latentInfo.acrossIndex = localIdx;
        latentInfo.withinIndex = [];
        latentInfo.acrossCategory = acrossCategory{localIdx};
        latentLine = sprintf('Across latent %d (%s)', localIdx, acrossCategory{localIdx});
    else
        latentInfo.latentType = 'within';
        latentInfo.acrossIndex = [];
        latentInfo.withinIndex = localIdx - xDim_across;
        latentInfo.acrossCategory = '';
        latentLine = sprintf('Within latent %d', latentInfo.withinIndex);
    end

    dslLabel = dsl_keep_remove_label(dslLogical, 'alltrials');
    dslByStimDirLabel = dsl_keep_remove_label(dslLogicalByStimDir, 'bystimdir');
    dslByStimNameDirLabel = dsl_keep_remove_label(dslLogicalByStimNameDir, 'bystimnamedir');
    dslByConditionLabel = dsl_keep_remove_label(dslLogicalByCondition, 'bycondition');

    latentInfo.dslLabel = dslLabel;
    latentInfo.dslByStimDirLabel = dslByStimDirLabel;
    latentInfo.dslByStimNameDirLabel = dslByStimNameDirLabel;
    latentInfo.dslByConditionLabel = dslByConditionLabel;
    latentInfo.titleLines = { ...
        sprintf('Group %d', groupIdx), ...
        latentLine, ...
        dslLabel, ...
        dslByStimDirLabel, ...
        dslByStimNameDirLabel, ...
        dslByConditionLabel};
end

% =========================================================================
% Compute condition-wise time-course means/SEMs and trial-response means/SEMs
% =========================================================================
function [meanTime, semTime, respMean, respSEM, respN, condTrialSeqIdx, condTrialIds] = ...
    compute_condition_averages(X, trialMean, condIdx, trialIds, nCond)

    T = size(X, 2);

    meanTime = nan(nCond, T);
    semTime  = nan(nCond, T);
    respMean = nan(nCond, 1);
    respSEM  = nan(nCond, 1);
    respN    = zeros(nCond, 1);

    condTrialSeqIdx = cell(nCond, 1);
    condTrialIds    = cell(nCond, 1);

    for c = 1:nCond
        idx = find(condIdx == c);
        condTrialSeqIdx{c} = idx;
        condTrialIds{c} = trialIds(idx);

        if isempty(idx)
            continue;
        end

        Xc = X(idx, :);
        yc = trialMean(idx);

        meanTime(c, :) = mean(Xc, 1, 'omitnan');
        if size(Xc, 1) > 1
            semTime(c, :) = std(Xc, 0, 1, 'omitnan') ./ sqrt(size(Xc, 1));
        else
            semTime(c, :) = zeros(1, T);
        end

        respMean(c) = mean(yc, 'omitnan');
        if numel(yc) > 1
            respSEM(c) = std(yc, 0, 'omitnan') ./ sqrt(numel(yc));
        else
            respSEM(c) = 0;
        end
        respN(c) = numel(yc);
    end
end

% =========================================================================
% Run full 3-way ANOVA
% =========================================================================
function anovaRes = run_threeway_anova(y, stimCode, sizeCode, contrastCode, alpha)

    valid = isfinite(y) & isfinite(stimCode) & isfinite(sizeCode) & isfinite(contrastCode);

    anovaRes = struct();
    anovaRes.p = nan(1, 7);
    anovaRes.tbl = [];
    anovaRes.stats = [];
    anovaRes.terms = [];
    anovaRes.termNames = {'Stim', 'Size', 'Contrast', ...
                          'Stim*Size', 'Stim*Contrast', 'Size*Contrast', ...
                          'Stim*Size*Contrast'};
    anovaRes.alpha = alpha;

    if sum(valid) < 8
        return;
    end

    try
        [p, tbl, stats, terms] = anovan( ...
            y(valid), ...
            {stimCode(valid), sizeCode(valid), contrastCode(valid)}, ...
            'model', 'full', ...
            'varnames', {'Stim', 'Size', 'Contrast'}, ...
            'display', 'off');

        anovaRes.p = p(:)';
        anovaRes.tbl = tbl;
        anovaRes.stats = stats;
        anovaRes.terms = terms;
    catch
        % Leave as NaN / empty if ANOVA fails
    end
end

% =========================================================================
% Compute marginal means/SEMs for the 3 main effects
% =========================================================================
function mainEffects = compute_main_effect_summaries(y, trialMeta, pVec)

    valid = trialMeta.valid & isfinite(y);

    mainEffects = struct([]);

    % Stim
    mainEffects(1).factor = 'Stim';
    mainEffects(1).levelLabels = {'grating', 'plaid'};
    mainEffects(1).means = [ ...
        mean(y(valid & trialMeta.stimCode == 1), 'omitnan'), ...
        mean(y(valid & trialMeta.stimCode == 2), 'omitnan')];
    mainEffects(1).sems = [ ...
        compute_sem(y(valid & trialMeta.stimCode == 1)), ...
        compute_sem(y(valid & trialMeta.stimCode == 2))];
    mainEffects(1).p = pick_p(pVec, 1);

    % Size
    mainEffects(2).factor = 'Size';
    mainEffects(2).levelLabels = {'small', 'large'};
    mainEffects(2).means = [ ...
        mean(y(valid & trialMeta.sizeCode == 1), 'omitnan'), ...
        mean(y(valid & trialMeta.sizeCode == 2), 'omitnan')];
    mainEffects(2).sems = [ ...
        compute_sem(y(valid & trialMeta.sizeCode == 1)), ...
        compute_sem(y(valid & trialMeta.sizeCode == 2))];
    mainEffects(2).p = pick_p(pVec, 2);

    % Contrast
    mainEffects(3).factor = 'Contrast';
    mainEffects(3).levelLabels = {'low', 'high'};
    mainEffects(3).means = [ ...
        mean(y(valid & trialMeta.contrastCode == 1), 'omitnan'), ...
        mean(y(valid & trialMeta.contrastCode == 2), 'omitnan')];
    mainEffects(3).sems = [ ...
        compute_sem(y(valid & trialMeta.contrastCode == 1)), ...
        compute_sem(y(valid & trialMeta.contrastCode == 2))];
    mainEffects(3).p = pick_p(pVec, 3);
end

% =========================================================================
% Compute all simple effects for the 2x2x2 design
% Each test compares the 2 levels of one factor while holding the other two fixed
% =========================================================================
function simpleEffects = compute_simple_effects(y, trialMeta, alpha)

    valid = trialMeta.valid & isfinite(y);

    template = struct( ...
        'factor', '', ...
        'context', '', ...
        'condPair', [NaN NaN], ...
        'xPair', [NaN NaN], ...
        'p', NaN, ...
        'n', [NaN NaN], ...
        'significant', false);

    % There are 12 simple-effect comparisons total in a 2x2x2 design:
    %   4 contrast simple effects
    %   4 size simple effects
    %   4 stim simple effects
    simpleEffects = repmat(template, 12, 1);
    ctr = 0;

    % Contrast simple effects within each Stim x Size
    for s = 1:2
        for z = 1:2
            idx1 = valid & trialMeta.stimCode == s & trialMeta.sizeCode == z & trialMeta.contrastCode == 1;
            idx2 = valid & trialMeta.stimCode == s & trialMeta.sizeCode == z & trialMeta.contrastCode == 2;

            ctr = ctr + 1;
            simpleEffects(ctr) = build_simple_effect_entry( ...
                template, ...
                'Contrast', ...
                sprintf('Stim=%d, Size=%d', s, z), ...
                cond_index(s, z, 1), cond_index(s, z, 2), ...
                trialMeta.barX(cond_index(s, z, 1)), ...
                trialMeta.barX(cond_index(s, z, 2)), ...
                y(idx1), y(idx2), alpha);
        end
    end

    % Size simple effects within each Stim x Contrast
    for s = 1:2
        for c = 1:2
            idx1 = valid & trialMeta.stimCode == s & trialMeta.sizeCode == 1 & trialMeta.contrastCode == c;
            idx2 = valid & trialMeta.stimCode == s & trialMeta.sizeCode == 2 & trialMeta.contrastCode == c;

            ctr = ctr + 1;
            simpleEffects(ctr) = build_simple_effect_entry( ...
                template, ...
                'Size', ...
                sprintf('Stim=%d, Contrast=%d', s, c), ...
                cond_index(s, 1, c), cond_index(s, 2, c), ...
                trialMeta.barX(cond_index(s, 1, c)), ...
                trialMeta.barX(cond_index(s, 2, c)), ...
                y(idx1), y(idx2), alpha);
        end
    end

    % Stim simple effects within each Size x Contrast
    for z = 1:2
        for c = 1:2
            idx1 = valid & trialMeta.stimCode == 1 & trialMeta.sizeCode == z & trialMeta.contrastCode == c;
            idx2 = valid & trialMeta.stimCode == 2 & trialMeta.sizeCode == z & trialMeta.contrastCode == c;

            ctr = ctr + 1;
            simpleEffects(ctr) = build_simple_effect_entry( ...
                template, ...
                'Stim', ...
                sprintf('Size=%d, Contrast=%d', z, c), ...
                cond_index(1, z, c), cond_index(2, z, c), ...
                trialMeta.barX(cond_index(1, z, c)), ...
                trialMeta.barX(cond_index(2, z, c)), ...
                y(idx1), y(idx2), alpha);
        end
    end

    simpleEffects = simpleEffects(1:ctr);
end

% =========================================================================
% Build one simple-effect entry
% =========================================================================
function S = build_simple_effect_entry(template, factorName, contextLabel, cond1, cond2, x1, x2, y1, y2, alpha)
    [pVal, n1, n2] = robust_ttest2(y1, y2);

    S = template;
    S.factor = factorName;
    S.context = contextLabel;
    S.condPair = [cond1, cond2];
    S.xPair = [x1, x2];
    S.p = pVal;
    S.n = [n1, n2];
    S.significant = isfinite(pVal) && (pVal < alpha);
end

% =========================================================================
% Robust two-sample t-test
% =========================================================================
function [pVal, n1, n2] = robust_ttest2(y1, y2)
    y1 = y1(isfinite(y1));
    y2 = y2(isfinite(y2));
    n1 = numel(y1);
    n2 = numel(y2);

    if n1 < 2 || n2 < 2
        pVal = NaN;
        return;
    end

    try
        [~, pVal] = ttest2(y1, y2, 'Vartype', 'unequal');
    catch
        pVal = NaN;
    end
end

% =========================================================================
% Plot one multi-panel time-course figure for one across latent
% =========================================================================
function [figFile, pngFile] = plot_across_timecourse_figure(Results, acrossIdx, timeDir)

    numGroups = Results.meta.numGroups;
    condLabels = Results.meta.conditionShortLabels;
    tAxis = Results.meta.timeAxis;

    category = Results.meta.acrossCategory{acrossIdx};

    % Short, stable file name. DSL keep/remove states remain in figure title.
    baseName = sanitize_filename(sprintf( ...
        'A%03d_%s_tc', ...
        acrossIdx, category));

    figFile = fullfile(timeDir, [baseName, '.fig']);
    pngFile = fullfile(timeDir, [baseName, '.png']);

    figW = max(1400, 420 * numGroups);
    fig = figure('Position', [100 100 figW 650]);
    tl = tiledlayout(fig, 1, numGroups, 'Padding', 'compact', 'TileSpacing', 'compact');

    axList = gobjects(numGroups, 1);
    legendHandles = gobjects(8, 1);

    for g = 1:numGroups
        axList(g) = nexttile(tl, g);
        hold(axList(g), 'on');

        meanTime = Results.group(g).latent(acrossIdx).condition.meanTime;
        semTime  = Results.group(g).latent(acrossIdx).condition.semTime;

        [legendHandlesLocal, ~] = plot_mean_sem_curves(axList(g), tAxis, meanTime, semTime);

        if g == 1
            legendHandles = legendHandlesLocal;
        end

        title(axList(g), Results.group(g).latent(acrossIdx).titleLines, 'Interpreter', 'none');
        xlabel(axList(g), 'Time bin');
        ylabel(axList(g), 'Latent response');
        box(axList(g), 'off');
    end

    validLegend = isgraphics(legendHandles);
    if any(validLegend)
        lgd = legend(axList(1), legendHandles(validLegend), condLabels(validLegend), ...
            'Interpreter', 'none', ...
            'Location', 'southoutside', ...
            'Orientation', 'horizontal');
        lgd.NumColumns = 4;
        lgd.Layout.Tile = 'south';
    end

    sgtitle(tl, sprintf('Condition-averaged time courses for across latent %d', acrossIdx), ...
        'Interpreter', 'none');

    savefig(fig, figFile);
    saveas(fig, pngFile);
    close(fig);
end

% =========================================================================
% Plot one single-panel time-course figure for one within latent
% =========================================================================
function [figFile, pngFile] = plot_within_timecourse_figure(latentEntry, tAxis, timeDir)

    g = latentEntry.groupIndex;
    w = latentEntry.withinIndex;

    % Short, stable file name. DSL keep/remove states remain in figure title.
    baseName = sanitize_filename(sprintf( ...
        'G%02d_W%03d_tc', ...
        g, w));

    figFile = fullfile(timeDir, [baseName, '.fig']);
    pngFile = fullfile(timeDir, [baseName, '.png']);

    fig = figure('Position', [100 100 1200 650]);
    tl = tiledlayout(fig, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    ax = nexttile(tl, 1);
    hold(ax, 'on');

    [legendHandles, ~] = plot_mean_sem_curves(ax, tAxis, latentEntry.condition.meanTime, latentEntry.condition.semTime);

    title(ax, latentEntry.titleLines, 'Interpreter', 'none');
    xlabel(ax, 'Time bin');
    ylabel(ax, 'Latent response');
    box(ax, 'off');

    validLegend = isgraphics(legendHandles);
    if any(validLegend)
        lgd = legend(ax, legendHandles(validLegend), latentEntry.condition.shortLabels(validLegend), ...
            'Interpreter', 'none', ...
            'Location', 'southoutside', ...
            'Orientation', 'horizontal');
        lgd.NumColumns = 4;
        lgd.Layout.Tile = 'south';
    end

    savefig(fig, figFile);
    saveas(fig, pngFile);
    close(fig);
end

% =========================================================================
% Plot 8 condition curves with SEM as shaded region
% =========================================================================
function [lineHandles, patchHandles] = plot_mean_sem_curves(ax, tAxis, meanTime, semTime)

    colors = lines(size(meanTime, 1));
    lineHandles = gobjects(size(meanTime, 1), 1);
    patchHandles = gobjects(size(meanTime, 1), 1);

    for c = 1:size(meanTime, 1)
        m = meanTime(c, :);
        s = semTime(c, :);

        if all(~isfinite(m))
            continue;
        end

        upper = m + s;
        lower = m - s;

        patchHandles(c) = fill(ax, [tAxis, fliplr(tAxis)], [upper, fliplr(lower)], colors(c, :), ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none');

        lineHandles(c) = plot(ax, tAxis, m, 'LineWidth', 1.5, 'Color', colors(c, :));
    end
end

% =========================================================================
% Plot the ANOVA summary figure for one group/latent
% =========================================================================
function [figFile, pngFile] = plot_anova_summary_figure(latentEntry, meta, anovaDir, alpha)

    g = latentEntry.groupIndex;

    % Short, stable file name. DSL keep/remove states remain in figure title.
    if strcmp(latentEntry.latentType, 'across')
        baseName = sanitize_filename(sprintf( ...
            'G%02d_A%03d_%s_anova', ...
            g, latentEntry.acrossIndex, latentEntry.acrossCategory));
    else
        baseName = sanitize_filename(sprintf( ...
            'G%02d_W%03d_anova', ...
            g, latentEntry.withinIndex));
    end

    figFile = fullfile(anovaDir, [baseName, '.fig']);
    pngFile = fullfile(anovaDir, [baseName, '.png']);

    fig = figure('Position', [100 100 1300 900]);
    tl = tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    % Left: full 8-condition plot spanning all rows
    axLeft = nexttile(tl, 1, [3 1]);
    plot_condition_bar_panel(axLeft, latentEntry, meta, alpha);

    % Right: main effects
    ax1 = nexttile(tl, 2);
    plot_main_effect_panel(ax1, latentEntry.anova.mainEffects(1), alpha);

    ax2 = nexttile(tl, 4);
    plot_main_effect_panel(ax2, latentEntry.anova.mainEffects(2), alpha);

    ax3 = nexttile(tl, 6);
    plot_main_effect_panel(ax3, latentEntry.anova.mainEffects(3), alpha);

    sgtitle(tl, latentEntry.titleLines, 'Interpreter', 'none');

    savefig(fig, figFile);
    saveas(fig, pngFile);
    close(fig);
end

% =========================================================================
% Plot the left condition bar panel with all significant simple-effect pairs
% =========================================================================
function plot_condition_bar_panel(ax, latentEntry, meta, alpha)

    xPos = meta.conditionBarX;
    means = latentEntry.condition.meanResponse(:)';
    sems  = latentEntry.condition.semResponse(:)';

    bar(ax, xPos, means, 0.8);
    hold(ax, 'on');
    errorbar(ax, xPos, means, sems, 'k.', 'LineWidth', 1);

    xline(ax, 6.5, '--');
    xline(ax, 3.0, ':');
    xline(ax, 10.0, ':');

    set(ax, 'XTick', xPos, 'XTickLabel', meta.conditionShortLabels);
    xtickangle(ax, 35);

    xlabel(ax, 'Condition');
    ylabel(ax, 'Trial-mean response');
    title(ax, 'Condition means \pm SEM');

    [ymin, ymax, yrng] = safe_y_range(means, sems);
    bracketH = 0.03 * yrng;
    yBase = ymax + 0.04 * yrng;
    yStep = 0.07 * yrng;

    sigSE = latentEntry.anova.simpleEffects;
    if ~isempty(sigSE)
        keep = arrayfun(@(s) s.significant, sigSE);
        sigSE = sigSE(keep);

        if ~isempty(sigSE)
            spans = arrayfun(@(s) abs(diff(s.xPair)), sigSE);
            x1all = arrayfun(@(s) s.xPair(1), sigSE);
            [~, ord] = sortrows([spans(:), x1all(:)], [1 2]);
            sigSE = sigSE(ord);

            for k = 1:numel(sigSE)
                starText = p_to_stars(sigSE(k).p);
                yThis = yBase + (k - 1) * yStep;
                add_sig_bracket(ax, sigSE(k).xPair(1), sigSE(k).xPair(2), yThis, starText, bracketH);
            end

            yTop = yBase + (numel(sigSE) - 1) * yStep + 2.2 * bracketH;
            ylim(ax, [ymin - 0.08 * yrng, yTop + 0.05 * yrng]);
        else
            ylim(ax, [ymin - 0.08 * yrng, ymax + 0.12 * yrng]);
        end
    else
        ylim(ax, [ymin - 0.08 * yrng, ymax + 0.12 * yrng]);
    end

    hold(ax, 'off');
end

% =========================================================================
% Plot one main-effect panel
% =========================================================================
function plot_main_effect_panel(ax, mainEffect, alpha)

    means = mainEffect.means(:)';
    sems  = mainEffect.sems(:)';

    bar(ax, 1:2, means, 0.8);
    hold(ax, 'on');
    errorbar(ax, 1:2, means, sems, 'k.', 'LineWidth', 1);

    set(ax, 'XTick', 1:2, 'XTickLabel', mainEffect.levelLabels);
    ylabel(ax, 'Response');
    title(ax, sprintf('%s main effect', mainEffect.factor), 'Interpreter', 'none');

    [ymin, ymax, yrng] = safe_y_range(means, sems);
    bracketH = 0.035 * yrng;
    yBase = ymax + 0.05 * yrng;

    if isfinite(mainEffect.p) && mainEffect.p < alpha
        add_sig_bracket(ax, 1, 2, yBase, p_to_stars(mainEffect.p), bracketH);
        ylim(ax, [ymin - 0.08 * yrng, yBase + 2.2 * bracketH + 0.05 * yrng]);
    else
        ylim(ax, [ymin - 0.08 * yrng, ymax + 0.12 * yrng]);
    end

    hold(ax, 'off');
end

% =========================================================================
% Draw one significance bracket
% =========================================================================
function add_sig_bracket(ax, x1, x2, y, txt, h)
    if nargin < 6 || ~isfinite(h) || h <= 0
        h = 0.05;
    end
    plot(ax, [x1 x1 x2 x2], [y y+h y+h y], 'k-', 'LineWidth', 1);
    text(ax, mean([x1 x2]), y + 1.15*h, txt, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
end

% =========================================================================
% Convert p-value to star label
% =========================================================================
function txt = p_to_stars(p)
    if ~isfinite(p)
        txt = 'n.s.';
    elseif p < 0.001
        txt = '***';
    elseif p < 0.01
        txt = '**';
    elseif p < 0.05
        txt = '*';
    else
        txt = 'n.s.';
    end
end

% =========================================================================
% Condition index helper
% Order:
%   1 grating-small-low
%   2 grating-small-high
%   3 grating-large-low
%   4 grating-large-high
%   5 plaid-small-low
%   6 plaid-small-high
%   7 plaid-large-low
%   8 plaid-large-high
% =========================================================================
function idx = cond_index(stimCode, sizeCode, contrastCode)
    idx = (stimCode - 1) * 4 + (sizeCode - 1) * 2 + contrastCode;
end

% =========================================================================
% SEM helper
% =========================================================================
function s = compute_sem(x)
    x = x(isfinite(x));
    if isempty(x)
        s = NaN;
    elseif numel(x) == 1
        s = 0;
    else
        s = std(x, 0) / sqrt(numel(x));
    end
end

% =========================================================================
% Safe p-value lookup
% =========================================================================
function p = pick_p(pVec, idx)
    if isempty(pVec) || numel(pVec) < idx
        p = NaN;
    else
        p = pVec(idx);
    end
end

% =========================================================================
% Make sure folder exists
% =========================================================================
function ensure_dir(folderPath)
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end
end

% =========================================================================
% Safe y-range helper
% =========================================================================
function [ymin, ymax, yrng] = safe_y_range(means, sems)
    vals = [means(:) + sems(:); means(:) - sems(:)];
    vals = vals(isfinite(vals));

    if isempty(vals)
        ymin = -1;
        ymax = 1;
        yrng = 2;
        return;
    end

    ymin = min(vals);
    ymax = max(vals);
    yrng = ymax - ymin;

    if ~isfinite(yrng) || yrng <= 0
        base = max(max(abs(vals)), 1);
        ymin = min(vals) - 0.5 * base;
        ymax = max(vals) + 0.5 * base;
        yrng = ymax - ymin;
    end
end

% =========================================================================
% Sanitize file names
% =========================================================================
function out = sanitize_filename(in)
    out = regexprep(in, '[^a-zA-Z0-9_\-]', '_');
    out = regexprep(out, '_+', '_');
end

% =========================================================================
% Convert label to token for file names
% =========================================================================
function out = label_to_token(in)
    out = lower(in);
    out = strrep(out, ' ', '_');
    out = sanitize_filename(out);
end

function out = dsl_keep_remove_label(flag, modeName)
    stem = dsl_mode_label_stem(modeName);
    if flag ~= 0
        out = [stem, ' keep'];
    else
        out = [stem, ' remove'];
    end
end

function out = dsl_mode_label_stem(modeName)
    switch lower(modeName)
        case 'alltrials'
            out = 'DSL';
        case 'bystimdir'
            out = 'DSL(by stim_dir)';
        case 'bystimnamedir'
            out = 'DSL(by stim_name x stim_dir)';
        case 'bycondition'
            out = 'DSL(by condition)';
        otherwise
            error('Unknown DSL modeName: %s', modeName);
    end
end

function all_tags = get_all_run_tags(model_data_allruns)
%% =========================================================================
% get_all_run_tags
%
% Purpose:
%   Extract all stim_tag values from model_data_allruns.
%
% Input:
%   model_data_allruns : cell array
%
% Output:
%   all_tags : cell array of char
% =========================================================================

all_tags = cell(numel(model_data_allruns), 1);

for j = 1:numel(model_data_allruns)
    if ~isfield(model_data_allruns{j}, 'stim_tag')
        error('stim_tag missing in model_data_allruns{%d}.', j);
    end
    all_tags{j} = model_data_allruns{j}.stim_tag;
end

end

function Results = analyze_dlag_latents_by_condition_split_by_dir(condition_full, seqEst, xDim_across, xDim_within, gp_params, ambiguousIdxs, DSL, tempfname)
% analyze_dlag_latents_by_condition_split_by_dir
%
% Analyze DLAG latent responses after splitting trials into two stimulus
% direction subsets:
%   - stim_dir1 : smaller effective direction value
%   - stim_dir2 : larger effective direction value
%
% Effective direction is defined per trial as:
%   - plaid   -> use plaid_dir
%   - grating -> use grating_dir
%
% After splitting by effective direction, each subset is analyzed
% independently using:
%   1) Condition-averaged time courses (2 x 2 x 2)
%   2) Three-way ANOVA on trial-averaged latent responses
%
% Three ANOVA factors:
%   - stim_name : grating vs plaid
%   - size      : small vs large
%   - contrast  : low vs high
%
% OUTPUT
%   Results : structure containing:
%       - trial metadata
%       - per-group / per-latent / per-stim_dir condition averages
%       - per-group / per-latent / per-stim_dir ANOVA results
%       - figure paths
%
% SAVED OUTPUTS
%   Under tempfname/latent_condition_analysis_split_by_dir/
%       - timecourse figures (.fig and .png)
%       - ANOVA figures (.fig and .png)
%       - latent_condition_analysis_split_by_dir_results.mat
%
% NOTES
%   - DSL.logical is the all-trials DSL keep/remove label.
%   - DSL.logical_bystimdir is the final AND result for the stim_dir split.
%   - DSL.logical_bystimnamedir is the final AND result for the 4-way
%     stim_name x stim_dir split.
%   - DSL.logical_bycondition is the final AND result across all conditions.

    alpha = 0.05;

    % ---------------------------------------------------------------------
    % Validate inputs
    % ---------------------------------------------------------------------
    numGroups = numel(xDim_within);
    localDims = xDim_across + xDim_within;

    if ~isstruct(DSL) || ~isfield(DSL, 'logical')
        error('DSL must contain field DSL.logical.');
    end
    if ~isfield(DSL, 'logical_bystimdir')
        error('DSL must contain field DSL.logical_bystimdir.');
    end
    if ~isfield(DSL, 'logical_bystimnamedir')
        error('DSL must contain field DSL.logical_bystimnamedir.');
    end
    if ~isfield(DSL, 'logical_bycondition')
        error('DSL must contain field DSL.logical_bycondition.');
    end

    if numel(DSL.logical) ~= numGroups
        error('DSL.logical must have one cell per group.');
    end
    if numel(DSL.logical_bystimdir) ~= numGroups
        error('DSL.logical_bystimdir must have one cell per group.');
    end
    if numel(DSL.logical_bystimnamedir) ~= numGroups
        error('DSL.logical_bystimnamedir must have one cell per group.');
    end
    if numel(DSL.logical_bycondition) ~= numGroups
        error('DSL.logical_bycondition must have one cell per group.');
    end

    if isempty(seqEst)
        error('seqEst is empty.');
    end

    trialLengths = arrayfun(@(s) size(s.xsm, 2), seqEst);
    if any(trialLengths ~= trialLengths(1))
        error('All seqEst trials must have the same number of time bins for this analysis.');
    end

    T = trialLengths(1);
    tAxis = 1:T;

    for g = 1:numGroups
        if numel(DSL.logical{g}) ~= localDims(g)
            error('DSL.logical{%d} length does not match xDim_across + xDim_within(%d).', g, g);
        end
        if numel(DSL.logical_bystimdir{g}) ~= localDims(g)
            error('DSL.logical_bystimdir{%d} length does not match xDim_across + xDim_within(%d).', g, g);
        end
        if numel(DSL.logical_bystimnamedir{g}) ~= localDims(g)
            error('DSL.logical_bystimnamedir{%d} length does not match xDim_across + xDim_within(%d).', g, g);
        end
        if numel(DSL.logical_bycondition{g}) ~= localDims(g)
            error('DSL.logical_bycondition{%d} length does not match xDim_across + xDim_within(%d).', g, g);
        end
    end

    % ---------------------------------------------------------------------
    % Output folders
    % ---------------------------------------------------------------------
    outDir   = fullfile(tempfname, 'latent_condition_analysis_split_by_dir');
    timeDir  = fullfile(outDir, 'timecourses');
    anovaDir = fullfile(outDir, 'anova');

    dirsplit_ensure_dir(outDir);
    dirsplit_ensure_dir(timeDir);
    dirsplit_ensure_dir(anovaDir);

    % ---------------------------------------------------------------------
    % Trial metadata
    % ---------------------------------------------------------------------
    trialMeta = dirsplit_build_trial_metadata(condition_full, seqEst);

    % ---------------------------------------------------------------------
    % Across-latent category labels
    % ---------------------------------------------------------------------
    acrossDelay    = gp_params.delays;
    acrossCategory = dirsplit_classify_across_latents(acrossDelay, ambiguousIdxs, xDim_across);

    % ---------------------------------------------------------------------
    % Global latent row indexing in seqEst(n).xsm
    % ---------------------------------------------------------------------
    blockStart = cumsum([1, localDims(1:end-1)]);

    % ---------------------------------------------------------------------
    % Initialize results structure
    % ---------------------------------------------------------------------
    Results = struct();
    Results.meta.alpha = alpha;
    Results.meta.timeAxis = tAxis;
    Results.meta.numGroups = numGroups;
    Results.meta.xDim_across = xDim_across;
    Results.meta.xDim_within = xDim_within;
    Results.meta.localDims = localDims;

    Results.meta.outputDir = outDir;
    Results.meta.timeDir = timeDir;
    Results.meta.anovaDir = anovaDir;

    Results.meta.factorOrder = {'stim_name', 'size', 'contrast'};
    Results.meta.splitVariable = 'stim_dir';
    Results.meta.stimLevels = trialMeta.stimLabels;
    Results.meta.sizeLevels = {'small', 'large'};
    Results.meta.sizeValues = trialMeta.sizeValues;
    Results.meta.contrastLevels = {'low', 'high'};
    Results.meta.contrastValuesByStim = trialMeta.contrastValuesByStim;

    Results.meta.stimDirLabels = trialMeta.stimDirLabels;
    Results.meta.stimDirValues = trialMeta.stimDirValues;

    Results.meta.conditionLabels = trialMeta.condLabels;
    Results.meta.conditionShortLabels = trialMeta.condShortLabels;
    Results.meta.conditionBarX = trialMeta.barX;

    Results.meta.acrossDelay = acrossDelay;
    Results.meta.acrossCategory = acrossCategory;
    Results.meta.ambiguousIdxs = ambiguousIdxs(:)';

    Results.meta.DSLRules.allTrials = 'Use DSL.logical.';
    Results.meta.DSLRules.byStimDir = 'Use DSL.logical_bystimdir (final AND result; do not use raw per-dir values as latent labels).';
    Results.meta.DSLRules.byStimNameDir = 'Use DSL.logical_bystimnamedir (final AND result across stim_name x stim_dir groups).';
    Results.meta.DSLRules.byCondition = 'Use DSL.logical_bycondition (final AND result across all conditions).';

    Results.trialMeta = trialMeta;

    % ---------------------------------------------------------------------
    % Compute all per-group / per-latent / per-stim_dir summaries
    % ---------------------------------------------------------------------
    for g = 1:numGroups
        Results.group(g).name = sprintf('Group %d', g);
        Results.group(g).latent = struct([]);

        for l = 1:localDims(g)
            rowIdx = blockStart(g) + l - 1;

            X = zeros(numel(seqEst), T);
            for n = 1:numel(seqEst)
                X(n, :) = seqEst(n).xsm(rowIdx, :);
            end

            trialMean = mean(X, 2);

            latentInfo = dirsplit_make_latent_info( ...
                g, l, xDim_across, acrossCategory, ...
                DSL.logical{g}(l), DSL.logical_bystimdir{g}(l), ...
                DSL.logical_bystimnamedir{g}(l), DSL.logical_bycondition{g}(l));

            % Store latent-level information shared across stim_dir subsets
            Results.group(g).latent(l).groupIndex = g;
            Results.group(g).latent(l).localLatentIndex = l;
            Results.group(g).latent(l).rowIndexInXsm = rowIdx;
            Results.group(g).latent(l).latentType = latentInfo.latentType;
            Results.group(g).latent(l).acrossIndex = latentInfo.acrossIndex;
            Results.group(g).latent(l).withinIndex = latentInfo.withinIndex;
            Results.group(g).latent(l).acrossCategory = latentInfo.acrossCategory;
            Results.group(g).latent(l).latentLine = latentInfo.latentLine;
            Results.group(g).latent(l).DSL.logical = DSL.logical{g}(l);
            Results.group(g).latent(l).DSL.label = latentInfo.dslLabel;
            Results.group(g).latent(l).DSL.logical_bystimdir = DSL.logical_bystimdir{g}(l);
            Results.group(g).latent(l).DSL.label_bystimdir = latentInfo.dslByStimDirLabel;
            Results.group(g).latent(l).DSL.logical_bystimnamedir = DSL.logical_bystimnamedir{g}(l);
            Results.group(g).latent(l).DSL.label_bystimnamedir = latentInfo.dslByStimNameDirLabel;
            Results.group(g).latent(l).DSL.logical_bycondition = DSL.logical_bycondition{g}(l);
            Results.group(g).latent(l).DSL.label_bycondition = latentInfo.dslByConditionLabel;

            % Analyze stim_dir1 and stim_dir2 independently
            for d = 1:2
                subsetMask = trialMeta.valid & (trialMeta.stimDirCode == d);

                [meanTime, semTime, respMean, respSEM, respN, condTrialSeqIdx, condTrialIds] = ...
                    dirsplit_compute_condition_averages( ...
                        X, trialMean, trialMeta.condIdx, trialMeta.trialId, subsetMask, 8);

                anovaRes = dirsplit_run_threeway_anova( ...
                    trialMean, trialMeta.stimCode, trialMeta.sizeCode, trialMeta.contrastCode, subsetMask, alpha);

                mainEffects = dirsplit_compute_main_effect_summaries( ...
                    trialMean, trialMeta, anovaRes.p, subsetMask);

                simpleEffects = dirsplit_compute_simple_effects( ...
                    trialMean, trialMeta, subsetMask, alpha);

                dirTitle = sprintf('%s = %s', ...
                    trialMeta.stimDirLabels{d}, ...
                    dirsplit_format_value(trialMeta.stimDirValues(d)));

                titleLines = { ...
                    sprintf('Group %d | %s', g, dirTitle), ...
                    latentInfo.latentLine, ...
                    latentInfo.dslLabel, ...
                    latentInfo.dslByStimDirLabel, ...
                    latentInfo.dslByStimNameDirLabel, ...
                    latentInfo.dslByConditionLabel};

                Results.group(g).latent(l).stim_dir(d).groupIndex = g;
                Results.group(g).latent(l).stim_dir(d).localLatentIndex = l;
                Results.group(g).latent(l).stim_dir(d).rowIndexInXsm = rowIdx;

                Results.group(g).latent(l).stim_dir(d).latentType = latentInfo.latentType;
                Results.group(g).latent(l).stim_dir(d).acrossIndex = latentInfo.acrossIndex;
                Results.group(g).latent(l).stim_dir(d).withinIndex = latentInfo.withinIndex;
                Results.group(g).latent(l).stim_dir(d).acrossCategory = latentInfo.acrossCategory;

                Results.group(g).latent(l).stim_dir(d).DSL.logical = DSL.logical{g}(l);
                Results.group(g).latent(l).stim_dir(d).DSL.label = latentInfo.dslLabel;
                Results.group(g).latent(l).stim_dir(d).DSL.logical_bystimdir = DSL.logical_bystimdir{g}(l);
                Results.group(g).latent(l).stim_dir(d).DSL.label_bystimdir = latentInfo.dslByStimDirLabel;
                Results.group(g).latent(l).stim_dir(d).DSL.logical_bystimnamedir = DSL.logical_bystimnamedir{g}(l);
                Results.group(g).latent(l).stim_dir(d).DSL.label_bystimnamedir = latentInfo.dslByStimNameDirLabel;
                Results.group(g).latent(l).stim_dir(d).DSL.logical_bycondition = DSL.logical_bycondition{g}(l);
                Results.group(g).latent(l).stim_dir(d).DSL.label_bycondition = latentInfo.dslByConditionLabel;

                Results.group(g).latent(l).stim_dir(d).stimDirCode = d;
                Results.group(g).latent(l).stim_dir(d).stimDirLabel = trialMeta.stimDirLabels{d};
                Results.group(g).latent(l).stim_dir(d).stimDirValue = trialMeta.stimDirValues(d);

                Results.group(g).latent(l).stim_dir(d).titleLines = titleLines;

                Results.group(g).latent(l).stim_dir(d).condition.labels = trialMeta.condLabels;
                Results.group(g).latent(l).stim_dir(d).condition.shortLabels = trialMeta.condShortLabels;
                Results.group(g).latent(l).stim_dir(d).condition.meanTime = meanTime;
                Results.group(g).latent(l).stim_dir(d).condition.semTime = semTime;
                Results.group(g).latent(l).stim_dir(d).condition.meanResponse = respMean;
                Results.group(g).latent(l).stim_dir(d).condition.semResponse = respSEM;
                Results.group(g).latent(l).stim_dir(d).condition.nTrials = respN;
                Results.group(g).latent(l).stim_dir(d).condition.trialSeqIdx = condTrialSeqIdx;
                Results.group(g).latent(l).stim_dir(d).condition.trialIds = condTrialIds;

                Results.group(g).latent(l).stim_dir(d).anova.responsePerTrial = trialMean(subsetMask);
                Results.group(g).latent(l).stim_dir(d).anova.p = anovaRes.p;
                Results.group(g).latent(l).stim_dir(d).anova.table = anovaRes.tbl;
                Results.group(g).latent(l).stim_dir(d).anova.stats = anovaRes.stats;
                Results.group(g).latent(l).stim_dir(d).anova.terms = anovaRes.terms;
                Results.group(g).latent(l).stim_dir(d).anova.termNames = anovaRes.termNames;
                Results.group(g).latent(l).stim_dir(d).anova.mainEffects = mainEffects;
                Results.group(g).latent(l).stim_dir(d).anova.simpleEffects = simpleEffects;
            end
        end
    end

    % ---------------------------------------------------------------------
    % Time-course figures
    %   - Across latents: one figure per across latent, rows = stim_dir,
    %     columns = groups
    %   - Within latents: one figure per group/within latent, panels = stim_dir
    % ---------------------------------------------------------------------
    Results.figures = struct();
    Results.figures.across = struct([]);

    for a = 1:xDim_across
        [figFile, pngFile] = dirsplit_plot_across_timecourse_figure(Results, a, timeDir);

        Results.figures.across(a).acrossIndex = a;
        Results.figures.across(a).figFile = figFile;
        Results.figures.across(a).pngFile = pngFile;

        for g = 1:numGroups
            Results.group(g).latent(a).timecourseSplitByDir.figFile = figFile;
            Results.group(g).latent(a).timecourseSplitByDir.pngFile = pngFile;
        end
    end

    for g = 1:numGroups
        for w = 1:xDim_within(g)
            l = xDim_across + w;
            [figFile, pngFile] = dirsplit_plot_within_timecourse_figure(Results.group(g).latent(l), Results.meta, timeDir);
            Results.group(g).latent(l).timecourseSplitByDir.figFile = figFile;
            Results.group(g).latent(l).timecourseSplitByDir.pngFile = pngFile;
        end
    end

    % ---------------------------------------------------------------------
    % ANOVA summary figures
    % ---------------------------------------------------------------------
    for g = 1:numGroups
        for l = 1:localDims(g)
            for d = 1:2
                [figFile, pngFile] = dirsplit_plot_anova_summary_figure( ...
                    Results.group(g).latent(l).stim_dir(d), Results.meta, anovaDir, alpha);

                Results.group(g).latent(l).stim_dir(d).anova.figFile = figFile;
                Results.group(g).latent(l).stim_dir(d).anova.pngFile = pngFile;
            end
        end
    end

    % ---------------------------------------------------------------------
    % Save results
    % ---------------------------------------------------------------------
    save(fullfile(outDir, 'latent_condition_analysis_split_by_dir_results.mat'), 'Results', '-v7.3');
end

% =========================================================================
% Build per-trial metadata including effective stimulus direction
% =========================================================================
function trialMeta = dirsplit_build_trial_metadata(condition_full, seqEst)

    Ntr = numel(seqEst);

    if isfield(seqEst, 'trialId')
        trialId = [seqEst.trialId]';
    else
        trialId = (1:Ntr)';
    end

    stimName       = strings(Ntr, 1);
    sizeValue      = nan(Ntr, 1);
    contrastRaw    = nan(Ntr, 1);
    effectiveDir   = nan(Ntr, 1);

    for k = 1:numel(condition_full)
        if ~isfield(condition_full(k), 'trial_indices')
            error('condition_full(%d) is missing field trial_indices.', k);
        end

        theseIds = condition_full(k).trial_indices(:);
        tf = ismember(trialId, theseIds);

        currStim = lower(string(condition_full(k).stim_name));

        stimName(tf)    = currStim;
        sizeValue(tf)   = condition_full(k).size;
        contrastRaw(tf) = condition_full(k).contrast;

        if strcmp(currStim, "plaid")
            if ~isfield(condition_full(k), 'plaid_dir')
                error('condition_full(%d) is plaid but missing field plaid_dir.', k);
            end
            effectiveDir(tf) = condition_full(k).plaid_dir;
        elseif strcmp(currStim, "grating")
            if ~isfield(condition_full(k), 'grating_dir')
                error('condition_full(%d) is grating but missing field grating_dir.', k);
            end
            effectiveDir(tf) = condition_full(k).grating_dir;
        else
            error('Unsupported stim_name found: %s', char(currStim));
        end
    end

    assigned = (stimName ~= "") & ~isnan(sizeValue) & ~isnan(contrastRaw) & ~isnan(effectiveDir);
    if any(~assigned)
        warning('Some seqEst trials could not be assigned completely. They will be excluded.');
    end

    % Stim labels
    allStim = unique(stimName(assigned), 'stable');
    allStim = lower(allStim);

    if all(ismember(["grating","plaid"], allStim))
        stimLabels = ["grating","plaid"];
    else
        if numel(allStim) ~= 2
            error('Expected exactly 2 stim levels after assignment.');
        end
        stimLabels = allStim(:)';
    end

    stimCode = nan(Ntr, 1);
    for s = 1:2
        stimCode(stimName == stimLabels(s)) = s;
    end

    % Size mapping
    sizeVals = unique(sizeValue(assigned));
    sizeVals = sort(sizeVals(:)');
    if numel(sizeVals) ~= 2
        error('Expected exactly 2 size values.');
    end

    sizeCode = nan(Ntr, 1);
    sizeCode(sizeValue == sizeVals(1)) = 1;
    sizeCode(sizeValue == sizeVals(2)) = 2;

    % Contrast mapping within each stim_name
    contrastCode = nan(Ntr, 1);
    contrastValuesByStim = struct();

    for s = 1:2
        idx = assigned & (stimCode == s);
        cvals = unique(contrastRaw(idx));
        cvals = sort(cvals(:)');
        if numel(cvals) ~= 2
            error('Stim %s does not have exactly 2 contrast levels.', char(stimLabels(s)));
        end

        contrastCode(idx & contrastRaw == cvals(1)) = 1;
        contrastCode(idx & contrastRaw == cvals(2)) = 2;
        contrastValuesByStim.(char(stimLabels(s))) = cvals;
    end

    % Effective stimulus direction mapping
    dirVals = unique(effectiveDir(assigned));
    dirVals = sort(dirVals(:)');

    if numel(dirVals) ~= 2
        error('Expected exactly 2 effective stimulus direction values.');
    end

    stimDirLabels = {'stim_dir1', 'stim_dir2'};
    stimDirCode = nan(Ntr, 1);

    tol = max(1e-10, 1e-8 * max(1, max(abs(dirVals))));
    stimDirCode(abs(effectiveDir - dirVals(1)) < tol) = 1;
    stimDirCode(abs(effectiveDir - dirVals(2)) < tol) = 2;

    valid = assigned & isfinite(stimCode) & isfinite(sizeCode) & isfinite(contrastCode) & isfinite(stimDirCode);

    condIdx = nan(Ntr, 1);
    condIdx(valid) = (stimCode(valid)-1) * 4 + (sizeCode(valid)-1) * 2 + contrastCode(valid);

    condLabels = { ...
        'grating-small-low'
        'grating-small-high'
        'grating-large-low'
        'grating-large-high'
        'plaid-small-low'
        'plaid-small-high'
        'plaid-large-low'
        'plaid-large-high'};

    condShortLabels = { ...
        'G-S-L'
        'G-S-H'
        'G-L-L'
        'G-L-H'
        'P-S-L'
        'P-S-H'
        'P-L-L'
        'P-L-H'};

    barX = [1 2 4 5 8 9 11 12];

    trialMeta = struct();
    trialMeta.trialId = trialId;
    trialMeta.stimName = cellstr(stimName);
    trialMeta.stimCode = stimCode;
    trialMeta.stimLabels = cellstr(stimLabels);

    trialMeta.sizeValue = sizeValue;
    trialMeta.sizeCode = sizeCode;
    trialMeta.sizeValues = sizeVals;

    trialMeta.contrastRaw = contrastRaw;
    trialMeta.contrastCode = contrastCode;
    trialMeta.contrastValuesByStim = contrastValuesByStim;

    trialMeta.effectiveDirValue = effectiveDir;
    trialMeta.stimDirCode = stimDirCode;
    trialMeta.stimDirLabels = stimDirLabels;
    trialMeta.stimDirValues = dirVals;

    trialMeta.valid = valid;
    trialMeta.condIdx = condIdx;

    trialMeta.condLabels = condLabels;
    trialMeta.condShortLabels = condShortLabels;
    trialMeta.barX = barX;
end

% =========================================================================
% Classify across latents into feedforward / feedback / ambiguous
% =========================================================================
function acrossCategory = dirsplit_classify_across_latents(acrossDelay, ambiguousIdxs, xDim_across)
    acrossCategory = repmat({''}, 1, xDim_across);

    ambiguousIdxs = unique(ambiguousIdxs(:)');
    ambiguousIdxs = ambiguousIdxs(ambiguousIdxs >= 1 & ambiguousIdxs <= xDim_across);

    for a = 1:xDim_across
        if ismember(a, ambiguousIdxs)
            acrossCategory{a} = 'ambiguous';
        else
            if isnan(acrossDelay(a)) || acrossDelay(a) == 0
                acrossCategory{a} = 'ambiguous';
            elseif acrossDelay(a) > 0
                acrossCategory{a} = 'feedforward';
            else
                acrossCategory{a} = 'feedback';
            end
        end
    end
end

% =========================================================================
% Build one latent info struct
% =========================================================================
function latentInfo = dirsplit_make_latent_info(groupIdx, localIdx, xDim_across, acrossCategory, dslLogical, dslLogicalByStimDir, dslLogicalByStimNameDir, dslLogicalByCondition)

    latentInfo = struct();
    latentInfo.groupIndex = groupIdx;

    if localIdx <= xDim_across
        latentInfo.latentType = 'across';
        latentInfo.acrossIndex = localIdx;
        latentInfo.withinIndex = [];
        latentInfo.acrossCategory = acrossCategory{localIdx};
        latentInfo.latentLine = sprintf('Across latent %d (%s)', localIdx, acrossCategory{localIdx});
    else
        latentInfo.latentType = 'within';
        latentInfo.acrossIndex = [];
        latentInfo.withinIndex = localIdx - xDim_across;
        latentInfo.acrossCategory = '';
        latentInfo.latentLine = sprintf('Within latent %d', localIdx - xDim_across);
    end

    latentInfo.dslLabel = dsl_keep_remove_label(dslLogical, 'alltrials');
    latentInfo.dslByStimDirLabel = dsl_keep_remove_label(dslLogicalByStimDir, 'bystimdir');
    latentInfo.dslByStimNameDirLabel = dsl_keep_remove_label(dslLogicalByStimNameDir, 'bystimnamedir');
    latentInfo.dslByConditionLabel = dsl_keep_remove_label(dslLogicalByCondition, 'bycondition');
end

% =========================================================================
% Compute condition-wise means/SEMs within one stim_dir subset
% =========================================================================
function [meanTime, semTime, respMean, respSEM, respN, condTrialSeqIdx, condTrialIds] = ...
    dirsplit_compute_condition_averages(X, trialMean, condIdx, trialIds, subsetMask, nCond)

    T = size(X, 2);

    meanTime = nan(nCond, T);
    semTime  = nan(nCond, T);
    respMean = nan(nCond, 1);
    respSEM  = nan(nCond, 1);
    respN    = zeros(nCond, 1);

    condTrialSeqIdx = cell(nCond, 1);
    condTrialIds    = cell(nCond, 1);

    for c = 1:nCond
        idx = find(subsetMask & (condIdx == c));
        condTrialSeqIdx{c} = idx;
        condTrialIds{c} = trialIds(idx);

        if isempty(idx)
            continue;
        end

        Xc = X(idx, :);
        yc = trialMean(idx);

        meanTime(c, :) = mean(Xc, 1, 'omitnan');

        if size(Xc, 1) > 1
            semTime(c, :) = std(Xc, 0, 1, 'omitnan') ./ sqrt(size(Xc, 1));
        else
            semTime(c, :) = zeros(1, T);
        end

        respMean(c) = mean(yc, 'omitnan');
        if numel(yc) > 1
            respSEM(c) = std(yc, 0, 'omitnan') ./ sqrt(numel(yc));
        else
            respSEM(c) = 0;
        end
        respN(c) = numel(yc);
    end
end

% =========================================================================
% Run 3-way ANOVA within one stim_dir subset
% =========================================================================
function anovaRes = dirsplit_run_threeway_anova(y, stimCode, sizeCode, contrastCode, subsetMask, alpha)

    valid = subsetMask & isfinite(y) & isfinite(stimCode) & isfinite(sizeCode) & isfinite(contrastCode);

    anovaRes = struct();
    anovaRes.p = nan(1, 7);
    anovaRes.tbl = [];
    anovaRes.stats = [];
    anovaRes.terms = [];
    anovaRes.termNames = {'Stim', 'Size', 'Contrast', ...
                          'Stim*Size', 'Stim*Contrast', 'Size*Contrast', ...
                          'Stim*Size*Contrast'};
    anovaRes.alpha = alpha;

    if sum(valid) < 8
        return;
    end

    try
        [p, tbl, stats, terms] = anovan( ...
            y(valid), ...
            {stimCode(valid), sizeCode(valid), contrastCode(valid)}, ...
            'model', 'full', ...
            'varnames', {'Stim', 'Size', 'Contrast'}, ...
            'display', 'off');

        anovaRes.p = p(:)';
        anovaRes.tbl = tbl;
        anovaRes.stats = stats;
        anovaRes.terms = terms;
    catch
        % Leave outputs as NaN / empty if ANOVA fails
    end
end

% =========================================================================
% Compute marginal means/SEMs for the three main effects within one subset
% =========================================================================
function mainEffects = dirsplit_compute_main_effect_summaries(y, trialMeta, pVec, subsetMask)

    valid = subsetMask & trialMeta.valid & isfinite(y);

    mainEffects = struct([]);

    % Stim
    mainEffects(1).factor = 'Stim';
    mainEffects(1).levelLabels = {'grating', 'plaid'};
    mainEffects(1).means = [ ...
        mean(y(valid & trialMeta.stimCode == 1), 'omitnan'), ...
        mean(y(valid & trialMeta.stimCode == 2), 'omitnan')];
    mainEffects(1).sems = [ ...
        dirsplit_compute_sem(y(valid & trialMeta.stimCode == 1)), ...
        dirsplit_compute_sem(y(valid & trialMeta.stimCode == 2))];
    mainEffects(1).p = dirsplit_pick_p(pVec, 1);

    % Size
    mainEffects(2).factor = 'Size';
    mainEffects(2).levelLabels = {'small', 'large'};
    mainEffects(2).means = [ ...
        mean(y(valid & trialMeta.sizeCode == 1), 'omitnan'), ...
        mean(y(valid & trialMeta.sizeCode == 2), 'omitnan')];
    mainEffects(2).sems = [ ...
        dirsplit_compute_sem(y(valid & trialMeta.sizeCode == 1)), ...
        dirsplit_compute_sem(y(valid & trialMeta.sizeCode == 2))];
    mainEffects(2).p = dirsplit_pick_p(pVec, 2);

    % Contrast
    mainEffects(3).factor = 'Contrast';
    mainEffects(3).levelLabels = {'low', 'high'};
    mainEffects(3).means = [ ...
        mean(y(valid & trialMeta.contrastCode == 1), 'omitnan'), ...
        mean(y(valid & trialMeta.contrastCode == 2), 'omitnan')];
    mainEffects(3).sems = [ ...
        dirsplit_compute_sem(y(valid & trialMeta.contrastCode == 1)), ...
        dirsplit_compute_sem(y(valid & trialMeta.contrastCode == 2))];
    mainEffects(3).p = dirsplit_pick_p(pVec, 3);
end

% =========================================================================
% Compute all simple effects within one stim_dir subset
% =========================================================================
function simpleEffects = dirsplit_compute_simple_effects(y, trialMeta, subsetMask, alpha)

    valid = subsetMask & trialMeta.valid & isfinite(y);

    template = struct( ...
        'factor', '', ...
        'context', '', ...
        'condPair', [NaN NaN], ...
        'xPair', [NaN NaN], ...
        'p', NaN, ...
        'n', [NaN NaN], ...
        'significant', false);

    simpleEffects = repmat(template, 12, 1);
    ctr = 0;

    % Contrast simple effects within each Stim x Size
    for s = 1:2
        for z = 1:2
            idx1 = valid & trialMeta.stimCode == s & trialMeta.sizeCode == z & trialMeta.contrastCode == 1;
            idx2 = valid & trialMeta.stimCode == s & trialMeta.sizeCode == z & trialMeta.contrastCode == 2;

            ctr = ctr + 1;
            simpleEffects(ctr) = dirsplit_build_simple_effect_entry( ...
                template, ...
                'Contrast', ...
                sprintf('Stim=%d, Size=%d', s, z), ...
                dirsplit_cond_index(s, z, 1), dirsplit_cond_index(s, z, 2), ...
                trialMeta.barX(dirsplit_cond_index(s, z, 1)), ...
                trialMeta.barX(dirsplit_cond_index(s, z, 2)), ...
                y(idx1), y(idx2), alpha);
        end
    end

    % Size simple effects within each Stim x Contrast
    for s = 1:2
        for c = 1:2
            idx1 = valid & trialMeta.stimCode == s & trialMeta.sizeCode == 1 & trialMeta.contrastCode == c;
            idx2 = valid & trialMeta.stimCode == s & trialMeta.sizeCode == 2 & trialMeta.contrastCode == c;

            ctr = ctr + 1;
            simpleEffects(ctr) = dirsplit_build_simple_effect_entry( ...
                template, ...
                'Size', ...
                sprintf('Stim=%d, Contrast=%d', s, c), ...
                dirsplit_cond_index(s, 1, c), dirsplit_cond_index(s, 2, c), ...
                trialMeta.barX(dirsplit_cond_index(s, 1, c)), ...
                trialMeta.barX(dirsplit_cond_index(s, 2, c)), ...
                y(idx1), y(idx2), alpha);
        end
    end

    % Stim simple effects within each Size x Contrast
    for z = 1:2
        for c = 1:2
            idx1 = valid & trialMeta.stimCode == 1 & trialMeta.sizeCode == z & trialMeta.contrastCode == c;
            idx2 = valid & trialMeta.stimCode == 2 & trialMeta.sizeCode == z & trialMeta.contrastCode == c;

            ctr = ctr + 1;
            simpleEffects(ctr) = dirsplit_build_simple_effect_entry( ...
                template, ...
                'Stim', ...
                sprintf('Size=%d, Contrast=%d', z, c), ...
                dirsplit_cond_index(1, z, c), dirsplit_cond_index(2, z, c), ...
                trialMeta.barX(dirsplit_cond_index(1, z, c)), ...
                trialMeta.barX(dirsplit_cond_index(2, z, c)), ...
                y(idx1), y(idx2), alpha);
        end
    end

    simpleEffects = simpleEffects(1:ctr);
end

% =========================================================================
% Build one simple-effect entry
% =========================================================================
function S = dirsplit_build_simple_effect_entry(template, factorName, contextLabel, cond1, cond2, x1, x2, y1, y2, alpha)
    [pVal, n1, n2] = dirsplit_robust_ttest2(y1, y2);

    S = template;
    S.factor = factorName;
    S.context = contextLabel;
    S.condPair = [cond1, cond2];
    S.xPair = [x1, x2];
    S.p = pVal;
    S.n = [n1, n2];
    S.significant = isfinite(pVal) && (pVal < alpha);
end

% =========================================================================
% Robust two-sample t-test
% =========================================================================
function [pVal, n1, n2] = dirsplit_robust_ttest2(y1, y2)
    y1 = y1(isfinite(y1));
    y2 = y2(isfinite(y2));

    n1 = numel(y1);
    n2 = numel(y2);

    if n1 < 2 || n2 < 2
        pVal = NaN;
        return;
    end

    try
        [~, pVal] = ttest2(y1, y2, 'Vartype', 'unequal');
    catch
        pVal = NaN;
    end
end

% =========================================================================
% Plot across-latent time-course figure
% Rows = stim_dir1 / stim_dir2
% Columns = groups
% =========================================================================
function [figFile, pngFile] = dirsplit_plot_across_timecourse_figure(Results, acrossIdx, timeDir)

    numGroups = Results.meta.numGroups;
    condLabels = Results.meta.conditionShortLabels;
    tAxis = Results.meta.timeAxis;
    category = Results.meta.acrossCategory{acrossIdx};

    % Short, stable file name. DSL keep/remove states remain in figure title.
    baseName = dirsplit_sanitize_filename(sprintf( ...
        'A%03d_%s_splitdir_tc', ...
        acrossIdx, category));

    figFile = fullfile(timeDir, [baseName, '.fig']);
    pngFile = fullfile(timeDir, [baseName, '.png']);

    figW = max(1500, 420 * numGroups);
    fig = figure('Position', [100 100 figW 950]);
    tl = tiledlayout(fig, 2, numGroups, 'Padding', 'compact', 'TileSpacing', 'compact');

    axList = gobjects(2 * numGroups, 1);
    legendHandles = gobjects(8, 1);

    for d = 1:2
        for g = 1:numGroups
            tileIdx = (d-1) * numGroups + g;
            ax = nexttile(tl, tileIdx);
            axList(tileIdx) = ax;
            hold(ax, 'on');

            meanTime = Results.group(g).latent(acrossIdx).stim_dir(d).condition.meanTime;
            semTime  = Results.group(g).latent(acrossIdx).stim_dir(d).condition.semTime;

            [legendHandlesLocal, ~] = dirsplit_plot_mean_sem_curves(ax, tAxis, meanTime, semTime);
            if d == 1 && g == 1
                legendHandles = legendHandlesLocal;
            end

            title(ax, Results.group(g).latent(acrossIdx).stim_dir(d).titleLines, 'Interpreter', 'none');
            xlabel(ax, 'Time bin');
            ylabel(ax, 'Latent response');
            box(ax, 'off');
        end
    end

    validLegend = isgraphics(legendHandles);
    if any(validLegend)
        lgd = legend(axList(1), legendHandles(validLegend), condLabels(validLegend), ...
            'Interpreter', 'none', ...
            'Location', 'southoutside', ...
            'Orientation', 'horizontal');
        lgd.NumColumns = 4;
        lgd.Layout.Tile = 'south';
    end

    sgtitle(tl, sprintf('Across latent %d split by stimulus direction', acrossIdx), 'Interpreter', 'none');

    savefig(fig, figFile);
    saveas(fig, pngFile);
    close(fig);
end

% =========================================================================
% Plot within-latent time-course figure
% Panels = stim_dir1 / stim_dir2
% =========================================================================
function [figFile, pngFile] = dirsplit_plot_within_timecourse_figure(latentEntry, meta, timeDir)

    g = latentEntry.groupIndex;
    w = latentEntry.withinIndex;

    % Short, stable file name. DSL keep/remove states remain in figure title.
    baseName = dirsplit_sanitize_filename(sprintf( ...
        'G%02d_W%03d_splitdir_tc', ...
        g, w));

    figFile = fullfile(timeDir, [baseName, '.fig']);
    pngFile = fullfile(timeDir, [baseName, '.png']);

    fig = figure('Position', [100 100 1450 800]);
    tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    axList = gobjects(2,1);
    legendHandles = gobjects(8, 1);

    for d = 1:2
        ax = nexttile(tl, d);
        axList(d) = ax;
        hold(ax, 'on');

        meanTime = latentEntry.stim_dir(d).condition.meanTime;
        semTime  = latentEntry.stim_dir(d).condition.semTime;

        [legendHandlesLocal, ~] = dirsplit_plot_mean_sem_curves(ax, meta.timeAxis, meanTime, semTime);
        if d == 1
            legendHandles = legendHandlesLocal;
        end

        title(ax, latentEntry.stim_dir(d).titleLines, 'Interpreter', 'none');
        xlabel(ax, 'Time bin');
        ylabel(ax, 'Latent response');
        box(ax, 'off');
    end

    validLegend = isgraphics(legendHandles);
    if any(validLegend)
        lgd = legend(axList(1), legendHandles(validLegend), meta.conditionShortLabels(validLegend), ...
            'Interpreter', 'none', ...
            'Location', 'southoutside', ...
            'Orientation', 'horizontal');
        lgd.NumColumns = 4;
        lgd.Layout.Tile = 'south';
    end

    sgtitle(tl, sprintf('Within latent %d split by stimulus direction', w), 'Interpreter', 'none');

    savefig(fig, figFile);
    saveas(fig, pngFile);
    close(fig);
end

% =========================================================================
% Plot condition curves with SEM as shaded region
% =========================================================================
function [lineHandles, patchHandles] = dirsplit_plot_mean_sem_curves(ax, tAxis, meanTime, semTime)

    colors = lines(size(meanTime, 1));
    lineHandles = gobjects(size(meanTime, 1), 1);
    patchHandles = gobjects(size(meanTime, 1), 1);

    for c = 1:size(meanTime, 1)
        m = meanTime(c, :);
        s = semTime(c, :);

        if all(~isfinite(m))
            continue;
        end

        upper = m + s;
        lower = m - s;

        patchHandles(c) = fill(ax, [tAxis, fliplr(tAxis)], [upper, fliplr(lower)], colors(c, :), ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none');

        lineHandles(c) = plot(ax, tAxis, m, 'LineWidth', 1.5, 'Color', colors(c, :));
    end
end

% =========================================================================
% Plot ANOVA summary figure for one group/latent/stim_dir subset
% =========================================================================
function [figFile, pngFile] = dirsplit_plot_anova_summary_figure(latentEntry, meta, anovaDir, alpha)

    g = latentEntry.groupIndex;
    dirLabel = latentEntry.stimDirLabel;
    dirValueStr = dirsplit_format_value(latentEntry.stimDirValue);

    % Short, stable file name. DSL keep/remove states remain in figure title.
    if strcmp(latentEntry.latentType, 'across')
        baseName = dirsplit_sanitize_filename(sprintf( ...
            'G%02d_A%03d_%s_%s%s_anova', ...
            g, latentEntry.acrossIndex, latentEntry.acrossCategory, ...
            dirLabel, dirValueStr));
    else
        baseName = dirsplit_sanitize_filename(sprintf( ...
            'G%02d_W%03d_%s%s_anova', ...
            g, latentEntry.withinIndex, ...
            dirLabel, dirValueStr));
    end

    figFile = fullfile(anovaDir, [baseName, '.fig']);
    pngFile = fullfile(anovaDir, [baseName, '.png']);

    fig = figure('Position', [100 100 1300 900]);
    tl = tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    axLeft = nexttile(tl, 1, [3 1]);
    dirsplit_plot_condition_bar_panel(axLeft, latentEntry, meta, alpha);

    ax1 = nexttile(tl, 2);
    dirsplit_plot_main_effect_panel(ax1, latentEntry.anova.mainEffects(1), alpha);

    ax2 = nexttile(tl, 4);
    dirsplit_plot_main_effect_panel(ax2, latentEntry.anova.mainEffects(2), alpha);

    ax3 = nexttile(tl, 6);
    dirsplit_plot_main_effect_panel(ax3, latentEntry.anova.mainEffects(3), alpha);

    sgtitle(tl, latentEntry.titleLines, 'Interpreter', 'none');

    savefig(fig, figFile);
    saveas(fig, pngFile);
    close(fig);
end

% =========================================================================
% Plot the 8-condition bar panel with all significant simple-effect pairs
% =========================================================================
function dirsplit_plot_condition_bar_panel(ax, latentEntry, meta, alpha)

    xPos = meta.conditionBarX;
    means = latentEntry.condition.meanResponse(:)';
    sems  = latentEntry.condition.semResponse(:)';

    bar(ax, xPos, means, 0.8);
    hold(ax, 'on');
    errorbar(ax, xPos, means, sems, 'k.', 'LineWidth', 1);

    xline(ax, 6.5, '--');
    xline(ax, 3.0, ':');
    xline(ax, 10.0, ':');

    set(ax, 'XTick', xPos, 'XTickLabel', meta.conditionShortLabels);
    xtickangle(ax, 35);

    xlabel(ax, 'Condition');
    ylabel(ax, 'Trial-mean response');
    title(ax, 'Condition means ± SEM');

    [ymin, ymax, yrng] = dirsplit_safe_y_range(means, sems);
    bracketH = 0.03 * yrng;
    yBase = ymax + 0.04 * yrng;
    yStep = 0.07 * yrng;

    sigSE = latentEntry.anova.simpleEffects;
    if ~isempty(sigSE)
        keep = arrayfun(@(s) s.significant, sigSE);
        sigSE = sigSE(keep);

        if ~isempty(sigSE)
            spans = arrayfun(@(s) abs(diff(s.xPair)), sigSE);
            x1all = arrayfun(@(s) s.xPair(1), sigSE);
            [~, ord] = sortrows([spans(:), x1all(:)], [1 2]);
            sigSE = sigSE(ord);

            for k = 1:numel(sigSE)
                starText = dirsplit_p_to_stars(sigSE(k).p);
                yThis = yBase + (k-1) * yStep;
                dirsplit_add_sig_bracket(ax, sigSE(k).xPair(1), sigSE(k).xPair(2), yThis, starText, bracketH);
            end

            yTop = yBase + (numel(sigSE) - 1) * yStep + 2.2 * bracketH;
            ylim(ax, [ymin - 0.08*yrng, yTop + 0.05*yrng]);
        else
            ylim(ax, [ymin - 0.08*yrng, ymax + 0.12*yrng]);
        end
    else
        ylim(ax, [ymin - 0.08*yrng, ymax + 0.12*yrng]);
    end

    hold(ax, 'off');
end

% =========================================================================
% Plot one main-effect panel
% =========================================================================
function dirsplit_plot_main_effect_panel(ax, mainEffect, alpha)

    means = mainEffect.means(:)';
    sems  = mainEffect.sems(:)';

    bar(ax, 1:2, means, 0.8);
    hold(ax, 'on');
    errorbar(ax, 1:2, means, sems, 'k.', 'LineWidth', 1);

    set(ax, 'XTick', 1:2, 'XTickLabel', mainEffect.levelLabels);
    ylabel(ax, 'Response');
    title(ax, sprintf('%s main effect', mainEffect.factor), 'Interpreter', 'none');

    [ymin, ymax, yrng] = dirsplit_safe_y_range(means, sems);
    bracketH = 0.035 * yrng;
    yBase = ymax + 0.05 * yrng;

    if isfinite(mainEffect.p) && mainEffect.p < alpha
        dirsplit_add_sig_bracket(ax, 1, 2, yBase, dirsplit_p_to_stars(mainEffect.p), bracketH);
        ylim(ax, [ymin - 0.08*yrng, yBase + 2.2*bracketH + 0.05*yrng]);
    else
        ylim(ax, [ymin - 0.08*yrng, ymax + 0.12*yrng]);
    end

    hold(ax, 'off');
end

% =========================================================================
% Safe y-range helper
% =========================================================================
function [ymin, ymax, yrng] = dirsplit_safe_y_range(means, sems)
    vals = [means(:) + sems(:); means(:) - sems(:)];
    vals = vals(isfinite(vals));

    if isempty(vals)
        ymin = -1;
        ymax = 1;
        yrng = 2;
        return;
    end

    ymin = min(vals);
    ymax = max(vals);
    yrng = ymax - ymin;

    if ~isfinite(yrng) || yrng <= 0
        base = max(max(abs(vals)), 1);
        ymin = min(vals) - 0.5 * base;
        ymax = max(vals) + 0.5 * base;
        yrng = ymax - ymin;
    end
end

% =========================================================================
% Draw one significance bracket
% =========================================================================
function dirsplit_add_sig_bracket(ax, x1, x2, y, txt, h)
    if nargin < 6 || ~isfinite(h) || h <= 0
        h = 0.05;
    end
    plot(ax, [x1 x1 x2 x2], [y y+h y+h y], 'k-', 'LineWidth', 1);
    text(ax, mean([x1 x2]), y + 1.15*h, txt, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontWeight', 'bold');
end

% =========================================================================
% Convert p-value to stars
% =========================================================================
function txt = dirsplit_p_to_stars(p)
    if ~isfinite(p)
        txt = 'n.s.';
    elseif p < 0.001
        txt = '***';
    elseif p < 0.01
        txt = '**';
    elseif p < 0.05
        txt = '*';
    else
        txt = 'n.s.';
    end
end

% =========================================================================
% Condition index helper
% Order:
%   1 grating-small-low
%   2 grating-small-high
%   3 grating-large-low
%   4 grating-large-high
%   5 plaid-small-low
%   6 plaid-small-high
%   7 plaid-large-low
%   8 plaid-large-high
% =========================================================================
function idx = dirsplit_cond_index(stimCode, sizeCode, contrastCode)
    idx = (stimCode - 1) * 4 + (sizeCode - 1) * 2 + contrastCode;
end

% =========================================================================
% SEM helper
% =========================================================================
function s = dirsplit_compute_sem(x)
    x = x(isfinite(x));
    if isempty(x)
        s = NaN;
    elseif numel(x) == 1
        s = 0;
    else
        s = std(x, 0) / sqrt(numel(x));
    end
end

% =========================================================================
% Safe p-value lookup
% =========================================================================
function p = dirsplit_pick_p(pVec, idx)
    if isempty(pVec) || numel(pVec) < idx
        p = NaN;
    else
        p = pVec(idx);
    end
end

% =========================================================================
% Create directory if needed
% =========================================================================
function dirsplit_ensure_dir(folderPath)
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end
end

% =========================================================================
% Sanitize file names
% =========================================================================
function out = dirsplit_sanitize_filename(in)
    out = regexprep(in, '[^a-zA-Z0-9_\-]', '_');
    out = regexprep(out, '_+', '_');
end

% =========================================================================
% Convert label to token for file names
% =========================================================================
function out = dirsplit_label_to_token(in)
    out = lower(in);
    out = strrep(out, ' ', '_');
    out = dirsplit_sanitize_filename(out);
end

% =========================================================================
% Format numeric value for titles and file names
% =========================================================================
function s = dirsplit_format_value(v)
    if ~isfinite(v)
        s = 'NaN';
        return;
    end

    if abs(v - round(v)) < 1e-10
        s = sprintf('%d', round(v));
    else
        s = sprintf('%.4g', v);
    end
end