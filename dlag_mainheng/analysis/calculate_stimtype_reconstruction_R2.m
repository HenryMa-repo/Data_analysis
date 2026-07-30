%% calculate_stimtype_reconstruction_R2.m
% Compute neuron-wise reconstruction R2 separately for grating and plaid.
%
% Run this script after data_reconstruction.m has added these fields to the
% corresponding bestmodel*.mat files:
%   seqEst.yRecon_use_all
%   seqEst.yRecon_use_across
%   seqEst.yRecon_use_within
%   seqEst.yRecon_use_feedforward
%   seqEst.yRecon_use_feedback
%
% Two model modes are supported:
%   1) data_condition = []
%        One all-condition model. Trials in seqEst are separated into
%        grating and plaid by:
%          seqEst.trialId
%          model_data_allruns{run}.condition_index_per_trial_full
%          model_data_allruns{run}.conditions_full(conditionId).stim_name
%
%   2) data_condition = 1:16, or another nonempty condition list
%        One condition-specific model per condition. Conditions are first
%        classified as grating or plaid from conditions_full, and all
%        selected conditions of the same stimulus type are pooled before
%        one R2 is computed.
%
% R2 definition:
%   R2 = 1 - RSS / TSS
%
% For each stimulus type, TSS is computed around the pooled per-neuron mean
% across all selected trials, time bins, and conditions belonging to that
% stimulus type. R2 is not averaged across conditions. Negative R2 values
% are retained.
%
% Unit IDs are read from the same model_data_allruns entry used to prepare
% the reconstruction data. For nan_trial_strategy == 6, data-content-
% specific unit-ID fields are required; otherwise standard probe unit-ID
% fields are required.
%
% Output variable:
%   stimtype_recon_R2
%
% Main output fields:
%   stimtype_recon_R2.grating.<reconstruction>.neuron_by_group
%   stimtype_recon_R2.plaid.<reconstruction>.neuron_by_group
%   stimtype_recon_R2.unit_ids_by_group
%   stimtype_recon_R2.yDims
%
% Output filename:
%   <data_content>_<model_mode>_stimtype_R2_all_across_within_ff_fb.mat
%
% Save location:
%   all_condition_model:
%       ./FA_Dlag_<data_content>/mat_results/runXXX/
%   condition_specific_models:
%       the folder containing this script

clc;
clear;

%% ------------------------------------------------------------------------
% User parameters
% -------------------------------------------------------------------------

data_content = 'demean_count_within_trial';
% Common options:
%   raw_count
%   raw_fr
%   z_within_trial
%   z_within_condition
%   z_across_conditions
%   demean_count_within_trial
%   demean_fr_within_trial
%   demean_pooledsd_within_condition

% [] uses one all-condition model.
% A nonempty vector uses condition-specific models and pools conditions of
% the same stimulus type.
data_condition = [];
% Example:
% data_condition = 1:16;

runIdx = 1;

% Metadata used to map trialId to stimulus type and to read group unit IDs.
dat_file = fullfile('.', 'model_data_allruns');
stim_tag = '_2[Gpl2_2c_2sz_400_2_200isi]';

% Group 1 corresponds to probe 0 and Group 2 corresponds to probe 1.
probe_ids = [0, 1];

%% ------------------------------------------------------------------------
% Reconstruction fields included in this analysis
% -------------------------------------------------------------------------

r2_specs = {
    'use_all',         'yRecon_use_all';
    'use_across',      'yRecon_use_across';
    'use_within',      'yRecon_use_within';
    'use_feedforward', 'yRecon_use_feedforward';
    'use_feedback',    'yRecon_use_feedback'
};

reconstruction_suffix = 'all_across_within_ff_fb';

%% ------------------------------------------------------------------------
% Main setup
% -------------------------------------------------------------------------

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

if isempty(data_condition)
    use_condition_mode = false;
    condition_list = [];
    model_mode = 'all_condition_model';
else
    use_condition_mode = true;
    condition_list = double(data_condition(:)');
    model_mode = 'condition_specific_models';

    if any(~isfinite(condition_list)) || ...
            any(condition_list < 1) || ...
            any(condition_list ~= round(condition_list))
        error('data_condition must contain positive integer condition IDs.');
    end

    if numel(unique(condition_list)) ~= numel(condition_list)
        error('data_condition contains duplicate condition IDs.');
    end
end

output_file_name = sprintf('%s_%s_stimtype_R2_%s.mat', ...
    data_content, model_mode, reconstruction_suffix);

fprintf('\n============================================================\n');
fprintf('Stimulus-type reconstruction R2\n');
fprintf('data_content : %s\n', data_content);
fprintf('model_mode   : %s\n', model_mode);
fprintf('runIdx       : %d\n', runIdx);
fprintf('stim_tag     : %s\n', stim_tag);
fprintf('dat_file     : %s\n', dat_file);
fprintf('============================================================\n');

%% ------------------------------------------------------------------------
% Load model metadata
% -------------------------------------------------------------------------

Sdata = loadMatFileFlexibleLocal(dat_file, 'model_data_allruns');
if ~isfield(Sdata, 'model_data_allruns')
    error('%s does not contain model_data_allruns.', dat_file);
end

model_data_allruns = Sdata.model_data_allruns;
metadata_run_idx = findRunByStimTagLocal(model_data_allruns, stim_tag);
metadata_run = model_data_allruns{metadata_run_idx};

requireFieldLocal(metadata_run, 'conditions_full', ...
    sprintf('model_data_allruns{%d}', metadata_run_idx));
requireFieldLocal(metadata_run, 'condition_index_per_trial_full', ...
    sprintf('model_data_allruns{%d}', metadata_run_idx));

conditions_full = metadata_run.conditions_full;
condition_index_per_trial_full = ...
    double(metadata_run.condition_index_per_trial_full(:));

if isempty(conditions_full)
    error('conditions_full is empty.');
end
if isempty(condition_index_per_trial_full)
    error('condition_index_per_trial_full is empty.');
end
if any(~isfinite(condition_index_per_trial_full)) || ...
        any(condition_index_per_trial_full < 1) || ...
        any(condition_index_per_trial_full ~= round(condition_index_per_trial_full))
    error('condition_index_per_trial_full contains invalid condition IDs.');
end
if any(condition_index_per_trial_full > numel(conditions_full))
    error('condition_index_per_trial_full references conditions outside conditions_full.');
end

validateAllConditionStimNamesLocal(conditions_full);

unit_ids_by_group = getModelUnitIdsByGroupLocal( ...
    metadata_run, data_content, probe_ids);

%% ------------------------------------------------------------------------
% Locate model files and build source list
% -------------------------------------------------------------------------

[source_specs, saveDir] = buildSourceSpecsLocal( ...
    data_content, condition_list, use_condition_mode, runIdx, ...
    conditions_full, scriptDir);

if ~exist(saveDir, 'dir')
    error('Output directory does not exist: %s', saveDir);
end

outFile = fullfile(saveDir, output_file_name);

fprintf('\nModel sources:\n');
for i = 1:numel(source_specs)
    if isempty(source_specs(i).condition_id)
        fprintf('  all-condition model: %s\n', source_specs(i).best_file);
    else
        fprintf('  condition %d (%s): %s\n', ...
            source_specs(i).condition_id, ...
            source_specs(i).stim_name, ...
            source_specs(i).best_file);
    end
end
fprintf('Output file: %s\n', outFile);

%% ------------------------------------------------------------------------
% Validate all best-model files and establish yDims
% -------------------------------------------------------------------------

yDims_ref = [];

for i = 1:numel(source_specs)
    S = load(source_specs(i).best_file, 'bestModel', 'res', 'seqEst');

    requiredVars = {'bestModel', 'res', 'seqEst'};
    for v = 1:numel(requiredVars)
        if ~isfield(S, requiredVars{v})
            error('File %s is missing variable %s.', ...
                source_specs(i).best_file, requiredVars{v});
        end
    end

    if ~isfield(S.res, 'estParams')
        error('File %s is missing res.estParams.', source_specs(i).best_file);
    end

    yDims = normalizeYDimsLocal(S.res.estParams);

    if isempty(yDims_ref)
        yDims_ref = yDims;
    elseif ~isequal(yDims_ref, yDims)
        error(['The model in %s has yDims %s, whereas the first model has ' ...
            'yDims %s. The neuron populations are not identical.'], ...
            source_specs(i).best_file, mat2str(yDims), mat2str(yDims_ref));
    end

    checkSeqEstFieldsLocal(S.seqEst, r2_specs, sum(yDims));

    if use_condition_mode
        verifyConditionSpecificTrialsLocal( ...
            S.seqEst, source_specs(i).condition_id, ...
            condition_index_per_trial_full, source_specs(i).best_file);
    else
        verifyAllConditionTrialsLocal( ...
            S.seqEst, condition_index_per_trial_full, ...
            conditions_full, source_specs(i).best_file);
    end
end

if numel(yDims_ref) ~= numel(probe_ids)
    error(['The reconstruction has %d groups, but probe_ids contains %d ' ...
        'entries.'], numel(yDims_ref), numel(probe_ids));
end

for g = 1:numel(yDims_ref)
    if numel(unit_ids_by_group{g}) ~= yDims_ref(g)
        error(['Group %d neuron count mismatch between model_data_allruns ' ...
            'and reconstruction: %d unit IDs versus yDims(%d) = %d.'], ...
            g, numel(unit_ids_by_group{g}), g, yDims_ref(g));
    end
end

%% ------------------------------------------------------------------------
% Compute grating and plaid R2 independently
% -------------------------------------------------------------------------

stim_types = {'grating', 'plaid'};

stimtype_recon_R2 = struct();
stimtype_recon_R2.data_content = data_content;
stimtype_recon_R2.model_mode = model_mode;
stimtype_recon_R2.data_condition = condition_list;
stimtype_recon_R2.runIdx = runIdx;
stimtype_recon_R2.stim_tag = stim_tag;
stimtype_recon_R2.metadata_run_idx = metadata_run_idx;
stimtype_recon_R2.dat_file = resolveMatFileLocal(dat_file);
stimtype_recon_R2.yDims = yDims_ref;
stimtype_recon_R2.probe_ids = probe_ids;
stimtype_recon_R2.unit_ids_by_group = unit_ids_by_group;
stimtype_recon_R2.reconstruction_specs = r2_specs;
stimtype_recon_R2.reconstruction_suffix = reconstruction_suffix;
stimtype_recon_R2.source_files = {source_specs.best_file};

for s = 1:numel(stim_types)
    stim_name = stim_types{s};

    fprintf('\n============================================================\n');
    fprintf('Computing pooled %s R2\n', stim_name);
    fprintf('============================================================\n');

    [stim_result, stim_info] = computeOneStimTypeR2Local( ...
        source_specs, stim_name, use_condition_mode, ...
        condition_index_per_trial_full, conditions_full, ...
        r2_specs, yDims_ref);

    stimtype_recon_R2.(stim_name) = stim_result;
    stimtype_recon_R2.(stim_name).info = stim_info;

    fprintf('%s: %d trials, %d concatenated time samples, conditions %s\n', ...
        stim_name, stim_info.n_trials, stim_info.n_time_samples, ...
        mat2str(stim_info.condition_ids));
end

%% ------------------------------------------------------------------------
% Save
% -------------------------------------------------------------------------

fprintf('\nSaving stimulus-type reconstruction R2:\n  %s\n', outFile);
save(outFile, 'stimtype_recon_R2', '-v7.3');
fprintf('Done.\n');

%% ========================================================================
% Local functions
% ========================================================================

function [source_specs, saveDir] = buildSourceSpecsLocal( ...
    data_content, condition_list, use_condition_mode, runIdx, ...
    conditions_full, scriptDir)

    source_specs = struct( ...
        'condition_id', {}, ...
        'stim_name', {}, ...
        'best_file', {});

    if use_condition_mode
        for i = 1:numel(condition_list)
            condID = condition_list(i);

            if condID > numel(conditions_full)
                error('Condition ID %d exceeds conditions_full length %d.', ...
                    condID, numel(conditions_full));
            end

            stim_name = getConditionStimNameLocal(conditions_full, condID);

            baseDir = ['./FA_Dlag_', data_content, ...
                '_condition', num2str(condID)];
            tempfname = sprintf('%s/mat_results/run%03d', baseDir, runIdx);
            bestFile = findOneFileLocal(tempfname, 'bestmodel*', true);

            source_specs(end+1).condition_id = condID; %#ok<AGROW>
            source_specs(end).stim_name = stim_name;
            source_specs(end).best_file = bestFile;
        end

        saveDir = scriptDir;
    else
        baseDir = ['./FA_Dlag_', data_content];
        tempfname = sprintf('%s/mat_results/run%03d', baseDir, runIdx);
        bestFile = findOneFileLocal(tempfname, 'bestmodel*', true);

        source_specs(1).condition_id = [];
        source_specs(1).stim_name = 'mixed';
        source_specs(1).best_file = bestFile;

        saveDir = tempfname;
    end
end

function [stim_result, stim_info] = computeOneStimTypeR2Local( ...
    source_specs, stim_name, use_condition_mode, ...
    condition_index_per_trial_full, conditions_full, ...
    r2_specs, yDims)

    yDim = sum(yDims);
    numGroups = numel(yDims);
    numSpecs = size(r2_specs, 1);

    sumY = cell(numSpecs, 1);
    countY = cell(numSpecs, 1);

    for r = 1:numSpecs
        sumY{r} = zeros(yDim, 1);
        countY{r} = zeros(yDim, 1);
    end

    nTrials = 0;
    nTimeSamples = 0;
    conditionIdsSeen = [];

    % First pass: pooled per-neuron means.
    for i = 1:numel(source_specs)
        S = load(source_specs(i).best_file, 'seqEst');

        [seqSel, condIdsSel] = selectSeqForStimLocal( ...
            S.seqEst, source_specs(i), stim_name, use_condition_mode, ...
            condition_index_per_trial_full, conditions_full);

        if isempty(seqSel)
            continue;
        end

        Ytrue = [seqSel.y];

        if size(Ytrue, 1) ~= yDim
            error('Ytrue in %s has %d rows; expected %d.', ...
                source_specs(i).best_file, size(Ytrue, 1), yDim);
        end

        nTrials = nTrials + numel(seqSel);
        nTimeSamples = nTimeSamples + size(Ytrue, 2);
        conditionIdsSeen = [conditionIdsSeen, condIdsSel(:)']; %#ok<AGROW>

        for r = 1:numSpecs
            fieldName = r2_specs{r, 2};
            Ypred = [seqSel.(fieldName)];

            if ~isequal(size(Ytrue), size(Ypred))
                error('%s size does not match y in %s.', ...
                    fieldName, source_specs(i).best_file);
            end

            valid = isfinite(Ytrue) & isfinite(Ypred);
            Ytmp = Ytrue;
            Ytmp(~valid) = 0;

            sumY{r} = sumY{r} + sum(Ytmp, 2);
            countY{r} = countY{r} + sum(valid, 2);
        end
    end

    if nTrials == 0
        error('No %s trials were found in the selected model source(s).', stim_name);
    end

    muY = cell(numSpecs, 1);
    for r = 1:numSpecs
        mu = sumY{r} ./ countY{r};
        mu(countY{r} == 0) = NaN;
        muY{r} = mu;
    end

    % Second pass: pooled RSS and TSS.
    RSS = cell(numSpecs, 1);
    TSS = cell(numSpecs, 1);

    for r = 1:numSpecs
        RSS{r} = zeros(yDim, 1);
        TSS{r} = zeros(yDim, 1);
    end

    for i = 1:numel(source_specs)
        S = load(source_specs(i).best_file, 'seqEst');

        seqSel = selectSeqForStimLocal( ...
            S.seqEst, source_specs(i), stim_name, use_condition_mode, ...
            condition_index_per_trial_full, conditions_full);

        if isempty(seqSel)
            continue;
        end

        Ytrue = [seqSel.y];

        for r = 1:numSpecs
            fieldName = r2_specs{r, 2};
            Ypred = [seqSel.(fieldName)];

            valid = isfinite(Ytrue) & isfinite(Ypred);

            D = Ytrue - repmat(muY{r}, 1, size(Ytrue, 2));
            E = Ytrue - Ypred;

            D(~valid) = 0;
            E(~valid) = 0;

            RSS{r} = RSS{r} + sum(E.^2, 2);
            TSS{r} = TSS{r} + sum(D.^2, 2);
        end
    end

    stim_result = struct();

    for r = 1:numSpecs
        r2Name = r2_specs{r, 1};
        rss = RSS{r};
        tss = TSS{r};

        neuronR2 = computeNeuronFromSumsLocal(rss, tss);

        stim_result.(r2Name).global_all = ...
            computeGlobalFromSumsLocal(rss, tss);
        stim_result.(r2Name).global_by_group = nan(1, numGroups);
        stim_result.(r2Name).neuron_by_group = cell(1, numGroups);
        stim_result.(r2Name).valid_sample_count_by_group = cell(1, numGroups);

        for g = 1:numGroups
            rows = getGroupRowsLocal(yDims, g);

            stim_result.(r2Name).global_by_group(g) = ...
                computeGlobalFromSumsLocal(rss(rows), tss(rows));
            stim_result.(r2Name).neuron_by_group{g} = neuronR2(rows);
            stim_result.(r2Name).valid_sample_count_by_group{g} = ...
                countY{r}(rows);
        end
    end

    stim_info = struct();
    stim_info.stim_name = stim_name;
    stim_info.condition_ids = unique(conditionIdsSeen, 'stable');
    stim_info.n_trials = nTrials;
    stim_info.n_time_samples = nTimeSamples;
end

function [seqSel, condIdsSel] = selectSeqForStimLocal( ...
    seqEst, source_spec, stim_name, use_condition_mode, ...
    condition_index_per_trial_full, conditions_full)

    if isempty(seqEst)
        seqSel = seqEst;
        condIdsSel = [];
        return;
    end

    if use_condition_mode
        if strcmp(source_spec.stim_name, stim_name)
            seqSel = seqEst;
            condIdsSel = repmat(source_spec.condition_id, 1, numel(seqSel));
        else
            seqSel = seqEst([]);
            condIdsSel = [];
        end
        return;
    end

    trialIds = getTrialIdsLocal(seqEst);
    condIds = condition_index_per_trial_full(trialIds);

    keep = false(size(condIds));
    for i = 1:numel(condIds)
        keep(i) = strcmp(getConditionStimNameLocal( ...
            conditions_full, condIds(i)), stim_name);
    end

    seqSel = seqEst(keep);
    condIdsSel = condIds(keep);
end

function verifyConditionSpecificTrialsLocal( ...
    seqEst, expectedCondition, condition_index_per_trial_full, sourceFile)

    trialIds = getTrialIdsLocal(seqEst);
    condIds = condition_index_per_trial_full(trialIds);

    if any(condIds ~= expectedCondition)
        bad = unique(condIds(condIds ~= expectedCondition));
        error(['Condition-specific model %s is expected to contain only ' ...
            'condition %d, but trialId metadata also maps to condition(s) %s.'], ...
            sourceFile, expectedCondition, mat2str(bad(:)'));
    end
end

function verifyAllConditionTrialsLocal( ...
    seqEst, condition_index_per_trial_full, conditions_full, sourceFile)

    trialIds = getTrialIdsLocal(seqEst);
    condIds = condition_index_per_trial_full(trialIds);

    for i = 1:numel(condIds)
        stim_name = getConditionStimNameLocal(conditions_full, condIds(i));
        if ~any(strcmp(stim_name, {'grating', 'plaid'}))
            error(['All-condition model %s contains trialId %d mapped to ' ...
                'unsupported stimulus type %s.'], ...
                sourceFile, trialIds(i), stim_name);
        end
    end
end

function trialIds = getTrialIdsLocal(seqEst)
    if isempty(seqEst)
        error('seqEst is empty.');
    end
    if ~isfield(seqEst, 'trialId')
        error(['seqEst is missing trialId. Stimulus type cannot be mapped ' ...
            'without trialId metadata.']);
    end

    trialIds = double(arrayfun(@(s) s.trialId, seqEst));
    trialIds = trialIds(:)';

    if any(~isfinite(trialIds)) || any(trialIds < 1) || ...
            any(trialIds ~= round(trialIds))
        error('seqEst.trialId contains invalid values.');
    end
end

function checkSeqEstFieldsLocal(seqEst, r2_specs, yDim)
    if isempty(seqEst)
        error('seqEst is empty.');
    end
    if ~isfield(seqEst, 'y')
        error('seqEst is missing field y.');
    end

    getTrialIdsLocal(seqEst);

    for n = 1:numel(seqEst)
        if isempty(seqEst(n).y)
            error('seqEst(%d).y is empty.', n);
        end
        if size(seqEst(n).y, 1) ~= yDim
            error('seqEst(%d).y has %d rows; expected %d.', ...
                n, size(seqEst(n).y, 1), yDim);
        end
    end

    for r = 1:size(r2_specs, 1)
        fieldName = r2_specs{r, 2};

        if ~isfield(seqEst, fieldName)
            error(['seqEst is missing reconstruction field %s. Run ' ...
                'data_reconstruction.m with the corresponding option first.'], ...
                fieldName);
        end

        for n = 1:numel(seqEst)
            if isempty(seqEst(n).(fieldName))
                error('seqEst(%d).%s is empty.', n, fieldName);
            end
            if ~isequal(size(seqEst(n).(fieldName)), size(seqEst(n).y))
                error('seqEst(%d).%s size does not match seqEst(%d).y.', ...
                    n, fieldName, n);
            end
        end
    end
end

function unit_ids_by_group = getModelUnitIdsByGroupLocal( ...
    metadata_run, data_content, probe_ids)

    requireFieldLocal(metadata_run, 'nan_trial_strategy', 'metadata run');
    nan_trial_strategy = double(metadata_run.nan_trial_strategy);

    unit_ids_by_group = cell(1, numel(probe_ids));

    for g = 1:numel(probe_ids)
        probe_id = probe_ids(g);

        if nan_trial_strategy == 6
            fieldName = sprintf('%s_probe%d_usedunit_ids', ...
                data_content, probe_id);
        else
            fieldName = sprintf('probe%d_usedunit_ids', probe_id);
        end

        if ~isfield(metadata_run, fieldName)
            error(['model_data_allruns entry is missing required unit-ID ' ...
                'field %s for data_content %s and nan_trial_strategy %g.'], ...
                fieldName, data_content, nan_trial_strategy);
        end

        ids = double(metadata_run.(fieldName));
        ids = ids(:);

        if isempty(ids)
            error('Unit-ID field %s is empty.', fieldName);
        end
        if any(~isfinite(ids))
            error('Unit-ID field %s contains nonfinite values.', fieldName);
        end
        if numel(unique(ids)) ~= numel(ids)
            error('Unit-ID field %s contains duplicate IDs.', fieldName);
        end

        unit_ids_by_group{g} = ids;
    end
end

function validateAllConditionStimNamesLocal(conditions_full)
    for condID = 1:numel(conditions_full)
        stim_name = getConditionStimNameLocal(conditions_full, condID);
        if ~any(strcmp(stim_name, {'grating', 'plaid'}))
            error('Unsupported stim_name in conditions_full(%d): %s', ...
                condID, stim_name);
        end
    end
end

function stim_name = getConditionStimNameLocal(conditions_full, condID)
    if condID < 1 || condID > numel(conditions_full)
        error('Condition ID %d is outside conditions_full.', condID);
    end

    if iscell(conditions_full)
        cond = conditions_full{condID};
    else
        cond = conditions_full(condID);
    end

    if ~isstruct(cond) || ~isfield(cond, 'stim_name')
        error('conditions_full condition %d is missing stim_name.', condID);
    end

    stim_name = lower(strtrim(char(cond.stim_name)));
end

function yDims = normalizeYDimsLocal(params)
    if ~isfield(params, 'yDims') || isempty(params.yDims)
        error('res.estParams is missing yDims.');
    end

    yDims = double(params.yDims(:)');

    if any(~isfinite(yDims)) || any(yDims <= 0) || ...
            any(yDims ~= round(yDims))
        error('res.estParams.yDims must contain positive integers.');
    end
end

function R2 = computeGlobalFromSumsLocal(rss, tss)
    rssTotal = sum(rss(:));
    tssTotal = sum(tss(:));

    if isfinite(rssTotal) && isfinite(tssTotal) && tssTotal > 0
        R2 = 1 - rssTotal / tssTotal;
    else
        R2 = NaN;
    end
end

function R2 = computeNeuronFromSumsLocal(rss, tss)
    R2 = nan(numel(rss), 1);
    valid = isfinite(rss) & isfinite(tss) & tss > 0;
    R2(valid) = 1 - rss(valid) ./ tss(valid);
end

function rows = getGroupRowsLocal(yDims, groupIdx)
    startIdx = sum(yDims(1:groupIdx-1)) + 1;
    rows = startIdx:(startIdx + yDims(groupIdx) - 1);
end

function run_idx = findRunByStimTagLocal(model_data_allruns, stim_tag)
    if ~iscell(model_data_allruns)
        error('model_data_allruns must be a cell array.');
    end

    hits = [];
    for i = 1:numel(model_data_allruns)
        entry = model_data_allruns{i};
        if isempty(entry) || ~isstruct(entry) || ~isfield(entry, 'stim_tag')
            continue;
        end
        if strcmp(entry.stim_tag, stim_tag)
            hits(end+1) = i; %#ok<AGROW>
        end
    end

    if isempty(hits)
        error('stim_tag not found in model_data_allruns: %s', stim_tag);
    end
    if numel(hits) > 1
        error('Duplicate stim_tag found in model_data_allruns: %s', stim_tag);
    end

    run_idx = hits;
end

function requireFieldLocal(S, fieldName, sourceName)
    if ~isfield(S, fieldName)
        error('%s is missing field %s.', sourceName, fieldName);
    end
end

function filePath = resolveMatFileLocal(fileBase)
    if exist(fileBase, 'file') == 2
        filePath = fileBase;
        return;
    end

    if ~endsWith(fileBase, '.mat')
        candidate = [fileBase, '.mat'];
        if exist(candidate, 'file') == 2
            filePath = candidate;
            return;
        end
    end

    error('MAT file not found: %s', fileBase);
end

function S = loadMatFileFlexibleLocal(fileBase, varargin)
    filePath = resolveMatFileLocal(fileBase);
    S = load(filePath, varargin{:});
end

function fname = findOneFileLocal(parentDir, pattern, mustExist)
    if ~exist(parentDir, 'dir')
        error('Directory not found: %s', parentDir);
    end

    files = dir(fullfile(parentDir, pattern));
    files = files(~[files.isdir]);

    if isempty(files)
        if mustExist
            error('No %s file found in %s.', pattern, parentDir);
        else
            fname = '';
            return;
        end
    end

    [~, idx] = sort([files.datenum], 'descend');
    files = files(idx);
    fname = fullfile(parentDir, files(1).name);
end
