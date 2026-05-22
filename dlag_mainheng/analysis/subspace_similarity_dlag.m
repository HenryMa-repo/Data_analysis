%% subspace_similarity_dlag.m
% Compute within-group DLAG subspace principal angles and subspace overlap.
%
% This version adds:
%   1) File names indicate latent selection mode:
%        subspace_similarity_all_latents.mat
%        subspace_similarity_DSL_<dsl_field>_filtered.mat
%   2) The saved SubspaceSim no longer contains SubspaceSim.meta,
%      use_dsl_filter, or dsl_field.
%   3) In condition mode, all requested conditions are summarized after the
%      loop. The summary mat keeps the same top-level structure as the
%      single-condition result:
%        SubspaceSim.classification
%        SubspaceSim.group
%      but scalar values are expanded across conditions.
%   4) Human-readable wide tables are printed and saved with two header rows.
%
% For every neural group, this computes:
%   1) across subspace vs within subspace
%   2) feedforward across subspace vs feedback across subspace
%
% Subspace overlap formula:
%   S(U,V) = 1 - ||(I - U*(U'*U)^(-1)*U')*V||_F / ||V||_F
% This is directional, so both S(A,B) and S(B,A) are computed.
%
% A captures B means B lies in A; equivalently, B projects little to A's
% null space.

clc;
clear;

% -------------------------------------------------------------------------
% User parameters
% -------------------------------------------------------------------------
data_content = 'z_across_conditions';
% options:
% raw_count, raw_fr, z_within_trial, z_within_condition,
% z_across_conditions, demean_count_within_trial, demean_fr_within_trial,
% demean_pooledsd_within_condition

data_condition = [1:16];
% [] for pooled all-condition mode, or e.g. 1:16 for condition mode

runIdx = 1;

% Whether to remove latents marked as DSL-remove by Latents_compare.m.
% If true, this script loads DSL_and_latent_category_stats.mat from each
% model folder and uses DSL.(dsl_field){groupIdx} as the keep mask.
use_dsl_filter = false;
dsl_field = 'logical';

% Stimulus metadata. Used only when data_condition is not empty, to make
% short condition labels such as G-S-L, G-L-H, P-S-L, etc.
dat_file = fullfile('.', 'model_data_allruns');
stim_tag = '_2[Gpl2_2c_2sz_400_2_200isi]';

% Save / print options
save_results = true;
print_single_model_verbose = false;
print_summary_tables = true;
save_summary_tables = true;

% Same outer location style as Latents_compare condition-mode outputs.
% Usually this script is run from the catfolder, so '.' is the catfolder.
summary_output_dir = '.';

% -------------------------------------------------------------------------
% Main loop setup
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

latentSelectionTag = makeLatentSelectionTagLocal(use_dsl_filter, dsl_field);
latentSelectionDisplay = makeLatentSelectionDisplayLocal(use_dsl_filter, dsl_field);

condition_full = [];
stim_abbrev = cell(1, numConditions);
conditionMap = [];

if use_condition_mode
    fprintf('Reading stimulus metadata from %s\n', dat_file);
    Sdata = load(dat_file, 'model_data_allruns');
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

AllSubspaceSim = cell(1, numConditions);
AllConditionInfo = struct([]);

% -------------------------------------------------------------------------
% Main loop: pooled all-condition mode or condition mode
% -------------------------------------------------------------------------
for cond_i = 1:numConditions
    if use_condition_mode
        this_condition = condition_list(cond_i);
        baseDir = ['./FA_Dlag_', data_content, '_condition', num2str(this_condition)];
        rowLabel = stim_abbrev{cond_i};
    else
        this_condition = [];
        baseDir = ['./FA_Dlag_', data_content];
        rowLabel = 'all';
    end

    tempfname = sprintf('%s/mat_results/run%03d', baseDir, runIdx);

    fprintf('\n============================================================\n');
    if isempty(this_condition)
        fprintf('DLAG subspace similarity: pooled all-condition mode\n');
    else
        fprintf('DLAG subspace similarity: condition %d (%s)\n', this_condition, rowLabel);
    end
    fprintf('Latent selection: %s\n', latentSelectionDisplay);
    fprintf('Reading from %s\n', tempfname);

    bestFile = findOneFileLocal(tempfname, 'bestmodel*', true);
    Sbest = load(bestFile);
    requiredBestVars = {'bestModel', 'res'};
    for v = 1:numel(requiredBestVars)
        if ~isfield(Sbest, requiredBestVars{v})
            error('File %s is missing variable %s.', bestFile, requiredBestVars{v});
        end
    end

    bestModel = Sbest.bestModel;
    res = Sbest.res;
    params = res.estParams;

    if isfield(Sbest, 'gp_params')
        gp_params = Sbest.gp_params;
    else
        gp_params = struct();
    end

    ambiguousIdxs = [];
    bootFile = findOneFileLocal(tempfname, 'bootstrapResults*', false);
    if ~isempty(bootFile)
        Sboot = load(bootFile);
        if isfield(Sboot, 'ambiguousIdxs')
            ambiguousIdxs = Sboot.ambiguousIdxs;
        end
    else
        warning('No bootstrapResults* file found in %s. Treating ambiguousIdxs as empty.', tempfname);
    end

    DSL = [];
    if use_dsl_filter
        dslFile = fullfile(tempfname, 'DSL_and_latent_category_stats.mat');
        if ~exist(dslFile, 'file')
            error(['use_dsl_filter=true, but %s was not found. ', ...
                   'Run Latents_compare.m first or set use_dsl_filter=false.'], dslFile);
        end
        Sdsl = load(dslFile, 'DSL');
        if ~isfield(Sdsl, 'DSL')
            error('%s does not contain DSL.', dslFile);
        end
        DSL = Sdsl.DSL;
    end

    opts = struct();
    opts.use_dsl_filter = use_dsl_filter;
    opts.dsl_field = dsl_field;

    SubspaceSim = computeDlagSubspaceSimilarityLocal(bestModel, params, gp_params, ambiguousIdxs, DSL, opts);

    if print_single_model_verbose
        printSubspaceSimilarityResultsLocal(SubspaceSim);
    end

    if save_results
        outFile = fullfile(tempfname, sprintf('subspace_similarity_%s.mat', latentSelectionTag));
        save(outFile, 'SubspaceSim');
        fprintf('Saved %s\n', outFile);
    end

    AllSubspaceSim{cond_i} = SubspaceSim;
    AllConditionInfo(cond_i).condition = this_condition;
    AllConditionInfo(cond_i).stim_abbrev = rowLabel;
    AllConditionInfo(cond_i).baseDir = baseDir;
    AllConditionInfo(cond_i).tempfname = tempfname;
end

% -------------------------------------------------------------------------
% Summary tables and condition-summary mat
% -------------------------------------------------------------------------
if use_condition_mode
    SubspaceSim = buildConditionSummarySubspaceSimLocal(AllSubspaceSim, condition_list, stim_abbrev);
    if save_results
        summaryMat = fullfile(summary_output_dir, sprintf( ...
            '%s_condition_mode_subspace_similarity_%s.mat', ...
            data_content, latentSelectionTag));
        save(summaryMat, 'SubspaceSim');
        fprintf('\nSaved condition-mode summary mat: %s\n', summaryMat);
    end

    if print_summary_tables || save_summary_tables
        for g = 1:numel(SubspaceSim.group)
            tableCell = buildConditionGroupTableCellLocal(SubspaceSim, g, stim_abbrev);
            titleStr = sprintf('%s condition mode subspace similarity - Group %d (%s)', ...
                data_content, g, latentSelectionDisplay);

            if print_summary_tables
                fprintf('\n============================================================\n');
                fprintf('%s\n', titleStr);
                printCellTableLocal(1, tableCell, [3 9 11]);
            end

            if save_summary_tables
                txtFile = fullfile(summary_output_dir, sprintf( ...
                    '%s_condition_mode_subspace_similarity_table_%s_group%d.txt', ...
                    data_content, latentSelectionTag, g));
                csvFile = fullfile(summary_output_dir, sprintf( ...
                    '%s_condition_mode_subspace_similarity_table_%s_group%d.csv', ...
                    data_content, latentSelectionTag, g));
                saveTableTextAndCsvLocal(tableCell, txtFile, csvFile, titleStr, [3 9 11]);
                fprintf('Saved summary table: %s\n', txtFile);
                fprintf('Saved summary table: %s\n', csvFile);
            end
        end
    end
else
    % Pooled all-condition mode: no extra summary mat is saved.
    % The single SubspaceSim mat above is the only mat output.
    if print_summary_tables || save_summary_tables
        SubspaceSim = AllSubspaceSim{1};
        tableCell = buildPooledGroupTableCellLocal(SubspaceSim);
        titleStr = sprintf('%s pooled all-condition mode subspace similarity (%s)', ...
            data_content, latentSelectionDisplay);

        if print_summary_tables
            fprintf('\n============================================================\n');
            fprintf('%s\n', titleStr);
            printCellTableLocal(1, tableCell, [3 9 11]);
        end

        if save_summary_tables && save_results
            pooledDir = AllConditionInfo(1).tempfname;
            txtFile = fullfile(pooledDir, sprintf('subspace_similarity_table_pooled_all_condition_mode_%s.txt', latentSelectionTag));
            csvFile = fullfile(pooledDir, sprintf('subspace_similarity_table_pooled_all_condition_mode_%s.csv', latentSelectionTag));
            saveTableTextAndCsvLocal(tableCell, txtFile, csvFile, titleStr, [3 9 11]);
        end
    end
end

%% ========================================================================
% Local functions
%% ========================================================================

function SubspaceSim = computeDlagSubspaceSimilarityLocal(bestModel, params, gp_params, ambiguousIdxs, DSL, opts)
    if ~isfield(params, 'C')
        error('params must contain loading matrix C.');
    end

    if isfield(bestModel, 'xDim_across')
        xDim_across = bestModel.xDim_across;
    elseif isfield(params, 'xDim_across')
        xDim_across = params.xDim_across;
    else
        error('Could not find xDim_across in bestModel or params.');
    end

    if isfield(bestModel, 'xDim_within')
        xDim_within = bestModel.xDim_within;
    elseif isfield(params, 'xDim_within')
        xDim_within = params.xDim_within;
    else
        error('Could not find xDim_within in bestModel or params.');
    end
    xDim_within = reshape(xDim_within, 1, []);

    if isfield(params, 'yDims')
        yDims = reshape(params.yDims, 1, []);
    else
        error('params must contain yDims.');
    end

    numGroups = numel(yDims);
    if numel(xDim_within) ~= numGroups
        error('Length of xDim_within (%d) must match number of groups (%d).', ...
            numel(xDim_within), numGroups);
    end

    localDims = xDim_across + xDim_within;

    if size(params.C, 2) ~= sum(localDims)
        error('params.C has %d columns, expected sum(xDim_across+xDim_within) = %d.', ...
            size(params.C, 2), sum(localDims));
    end
    if size(params.C, 1) ~= sum(yDims)
        error('params.C has %d rows, expected sum(yDims) = %d.', ...
            size(params.C, 1), sum(yDims));
    end

    latentClass = classifyDlagLatentsLocal(xDim_across, params, gp_params, ambiguousIdxs);

    SubspaceSim = struct();
    SubspaceSim.classification = latentClass;

    obsStart = cumsum([1, yDims(1:end-1)]);
    obsEnd = cumsum(yDims);
    latStart = cumsum([1, localDims(1:end-1)]);
    latEnd = cumsum(localDims);

    for g = 1:numGroups
        obsIdx = obsStart(g):obsEnd(g);
        latIdx = latStart(g):latEnd(g);
        Cg = params.C(obsIdx, latIdx);

        groupKeepMask = true(1, localDims(g));
        if opts.use_dsl_filter
            groupKeepMask = getDslKeepMaskLocal(DSL, opts.dsl_field, g, localDims(g));
        end

        acrossLocal = 1:xDim_across;
        withinLocal = (xDim_across+1):localDims(g);
        ffLocal = latentClass.feedforwardIdx;
        fbLocal = latentClass.feedbackIdx;

        acrossLocal = intersect(acrossLocal, find(groupKeepMask), 'stable');
        withinLocal = intersect(withinLocal, find(groupKeepMask), 'stable');
        ffLocal = intersect(ffLocal, find(groupKeepMask), 'stable');
        fbLocal = intersect(fbLocal, find(groupKeepMask), 'stable');

        SubspaceSim.group(g).name = sprintf('Group %d', g);

        pairSpecs = {
            'across_vs_within',        'Across',      'Within',   acrossLocal, withinLocal;
            'feedforward_vs_feedback', 'Feedforward', 'Feedback', ffLocal,     fbLocal
        };

        SubspaceSim.group(g).pairNames = pairSpecs(:, 1)';
        SubspaceSim.group(g).pair = cell(1, size(pairSpecs, 1));

        for p = 1:size(pairSpecs, 1)
            pairName = pairSpecs{p, 1};
            labelA = pairSpecs{p, 2};
            labelB = pairSpecs{p, 3};
            idxA = pairSpecs{p, 4};
            idxB = pairSpecs{p, 5};

            Araw = Cg(:, idxA);
            Braw = Cg(:, idxB);

            pairResult = compareTwoSubspacesLocal(Araw, Braw, labelA, labelB);
            pairResult.name = pairName;
            pairResult.labelA = labelA;
            pairResult.labelB = labelB;
            pairResult.idxA = idxA;
            pairResult.idxB = idxB;
            pairResult.rawDimA = size(Araw, 2);
            pairResult.rawDimB = size(Braw, 2);

            SubspaceSim.group(g).pair{p} = pairResult;
        end
    end
end

function pairResult = compareTwoSubspacesLocal(Araw, Braw, labelA, labelB)
    pairResult = struct();
    pairResult.status = 'ok';
    pairResult.warning = '';

    QA = dlagSvdBasisLocal(Araw);
    QB = dlagSvdBasisLocal(Braw);
    pairResult.basisA = QA;
    pairResult.basisB = QB;

    captureFieldAB = makeCaptureFieldNameLocal(labelA, labelB);
    captureFieldBA = makeCaptureFieldNameLocal(labelB, labelA);

    if isempty(QA) || isempty(QB) || size(QA, 2) == 0 || size(QB, 2) == 0
        pairResult.status = 'skipped_empty_subspace';
        pairResult.warning = sprintf('%s or %s has zero usable dimension.', labelA, labelB);
        pairResult.principal.cosine = [];
        pairResult.principal.angle_deg = [];
        pairResult.principal.first_angle_deg = NaN;
        pairResult.principal.last_angle_deg = NaN;
        pairResult.principal.median_angle_deg = NaN;
        pairResult.similarity.(captureFieldAB) = NaN;
        pairResult.similarity.(captureFieldBA) = NaN;
        pairResult.similarity.avg = NaN;
        return;
    end

    % Principal angles: singular values of QA' * QB are principal cosines.
    s = svd(QA' * QB, 'econ');
    s = min(max(s, 0), 1);
    thetaRad = acos(s);
    thetaDeg = thetaRad * 180 / pi;

    pairResult.principal.cosine = reshape(s, 1, []);
    pairResult.principal.angle_deg = reshape(thetaDeg, 1, []);
    pairResult.principal.first_angle_deg = thetaDeg(1);
    pairResult.principal.last_angle_deg = thetaDeg(end);
    pairResult.principal.median_angle_deg = median(thetaDeg);

    % Directional subspace overlap. Computed twice because the measure is
    % not symmetric.
    pairResult.similarity.(captureFieldAB) = directionalSubspaceOverlapLocal(QA, QB);
    pairResult.similarity.(captureFieldBA) = directionalSubspaceOverlapLocal(QB, QA);
    simVals = [pairResult.similarity.(captureFieldAB), pairResult.similarity.(captureFieldBA)];
    simVals = simVals(~isnan(simVals));
    if isempty(simVals)
        pairResult.similarity.avg = NaN;
    else
        pairResult.similarity.avg = mean(simVals);
    end
end

function fieldName = makeCaptureFieldNameLocal(labelA, labelB)
    % Make explicit overlap field names, e.g.:
    % across_captures_within, within_captures_across,
    % feedforward_captures_feedback, feedback_captures_feedforward.
    fieldName = sprintf('%s_captures_%s', normalizeCaptureLabelLocal(labelA), normalizeCaptureLabelLocal(labelB));
end

function label = normalizeCaptureLabelLocal(labelIn)
    label = lower(char(string(labelIn)));
    label = regexprep(label, '[^a-z0-9]+', '_');
    label = regexprep(label, '^_+|_+$', '');
    if isempty(label)
        label = 'subspace';
    end
end

function overlap = directionalSubspaceOverlapLocal(U, V)
    % S(U,V) = 1 - ||(I - P_U)V||_F / ||V||_F
    % P_U = U * inv(U' * U) * U'
    % U and V are already DLAG-style SVD bases, so U'*U should be close to I.

    if isempty(U) || isempty(V) || size(U, 2) == 0 || size(V, 2) == 0
        overlap = NaN;
        return;
    end

    if size(U, 1) ~= size(V, 1)
        error('U and V must live in the same observed neural space.');
    end

    denom = norm(V, 'fro');
    if denom <= eps
        overlap = NaN;
        return;
    end

    P_U = U * ((U' * U) \ U');
    residual = (eye(size(U, 1)) - P_U) * V;
    overlap = 1 - norm(residual, 'fro') / denom;

    if overlap < 0 && overlap > -1e-12
        overlap = 0;
    end
end

function Q = dlagSvdBasisLocal(L)
    % Build an orthonormal basis using the same SVD idea as DLAG's
    % orthogonalize.m and dominantModes_dlag.m.
    % This function intentionally does not call MATLAB orth().

    if isempty(L) || size(L, 2) == 0
        Q = zeros(size(L, 1), 0);
        return;
    end

    xDim = size(L, 2);

    if xDim == 1
        mag = sqrt(L' * L);
        if mag <= eps
            Q = zeros(size(L, 1), 0);
        else
            Q = L / mag;
        end
        return;
    end

    [UU, DD, ~] = svd(L, 'econ');
    s = diag(DD);
    if isempty(s) || max(s) <= eps
        Q = zeros(size(L, 1), 0);
        return;
    end

    r = min(xDim, size(UU, 2));
    Q = UU(:, 1:r);
end

function latentClass = classifyDlagLatentsLocal(xDim_across, params, gp_params, ambiguousIdxs)
    % Classify across latents into feedforward / feedback / ambiguous.
    % positive delay -> feedforward
    % negative delay -> feedback
    % zero / NaN / bootstrap-ambiguous -> ambiguous

    acrossDelay = [];
    if isstruct(gp_params) && isfield(gp_params, 'delays') && ~isempty(gp_params.delays)
        acrossDelay = reshape(gp_params.delays, 1, []);
    elseif isfield(params, 'DelayMatrix') && ~isempty(params.DelayMatrix)
        if size(params.DelayMatrix, 1) >= 2
            acrossDelay = params.DelayMatrix(2, :) - params.DelayMatrix(1, :);
        else
            acrossDelay = params.DelayMatrix(1, :);
        end
    end

    if isempty(acrossDelay)
        acrossDelay = nan(1, xDim_across);
    end

    if numel(acrossDelay) < xDim_across
        error('Delay vector has fewer entries than xDim_across.');
    end
    acrossDelay = acrossDelay(1:xDim_across);

    ambiguousIdxs = normalizeAmbiguousIdxsLocal(ambiguousIdxs, xDim_across);

    acrossIdx = 1:xDim_across;
    zeroOrNaNIdx = acrossIdx((acrossDelay == 0) | isnan(acrossDelay));
    ambiguousAll = unique([ambiguousIdxs, zeroOrNaNIdx]);

    ffIdx = find(acrossDelay > 0);
    fbIdx = find(acrossDelay < 0);
    ffIdx = setdiff(ffIdx, ambiguousAll, 'stable');
    fbIdx = setdiff(fbIdx, ambiguousAll, 'stable');

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

function ambiguousIdxs = normalizeAmbiguousIdxsLocal(ambiguousIdxs, xDim_across)
    if isempty(ambiguousIdxs)
        ambiguousIdxs = [];
        return;
    end

    if iscell(ambiguousIdxs)
        tmp = [];
        for ii = 1:numel(ambiguousIdxs)
            if ~isempty(ambiguousIdxs{ii})
                tmp = [tmp, reshape(ambiguousIdxs{ii}, 1, [])]; %#ok<AGROW>
            end
        end
        ambiguousIdxs = tmp;
    elseif islogical(ambiguousIdxs)
        ambiguousIdxs = find(ambiguousIdxs);
    else
        ambiguousIdxs = reshape(ambiguousIdxs, 1, []);
    end

    ambiguousIdxs = unique(ambiguousIdxs);
    ambiguousIdxs = ambiguousIdxs(ambiguousIdxs >= 1 & ambiguousIdxs <= xDim_across);
end

function keepMask = getDslKeepMaskLocal(DSL, fieldName, groupIdx, localDim)
    if isempty(DSL) || ~isstruct(DSL) || ~isfield(DSL, fieldName)
        error('DSL must contain field %s when use_dsl_filter=true.', fieldName);
    end

    fieldVal = DSL.(fieldName);
    if ~iscell(fieldVal) || numel(fieldVal) < groupIdx
        error('DSL.%s must be a cell array with one entry per group.', fieldName);
    end

    keepMask = reshape(fieldVal{groupIdx}, 1, []) ~= 0;
    if numel(keepMask) ~= localDim
        error('DSL.%s{%d} has length %d, expected %d.', ...
            fieldName, groupIdx, numel(keepMask), localDim);
    end
end

function Summary = buildConditionSummarySubspaceSimLocal(AllSubspaceSim, condition_list, stim_abbrev)
    nCond = numel(AllSubspaceSim);
    if nCond < 1
        error('AllSubspaceSim is empty.');
    end

    firstSim = AllSubspaceSim{1};
    numGroups = numel(firstSim.group);

    Summary = struct();
    Summary.classification = buildSummaryClassificationLocal(AllSubspaceSim, condition_list, stim_abbrev);

    for g = 1:numGroups
        Summary.group(g).name = sprintf('Group %d', g);
        Summary.group(g).pairNames = firstSim.group(g).pairNames;
        Summary.group(g).pair = cell(1, numel(firstSim.group(g).pair));

        for p = 1:numel(firstSim.group(g).pair)
            pr0 = firstSim.group(g).pair{p};
            prSum = struct();
            prSum.name = pr0.name;
            prSum.labelA = pr0.labelA;
            prSum.labelB = pr0.labelB;
            prSum.condition_id = condition_list;
            prSum.stim_abbrev = stim_abbrev;
            prSum.status = cell(1, nCond);
            prSum.warning = cell(1, nCond);
            prSum.idxA = cell(1, nCond);
            prSum.idxB = cell(1, nCond);
            prSum.rawDimA = nan(1, nCond);
            prSum.rawDimB = nan(1, nCond);
            prSum.basisA = cell(1, nCond);
            prSum.basisB = cell(1, nCond);
            prSum.principal.cosine = cell(1, nCond);
            prSum.principal.angle_deg = cell(1, nCond);
            prSum.principal.first_angle_deg = nan(1, nCond);
            prSum.principal.last_angle_deg = nan(1, nCond);
            prSum.principal.median_angle_deg = nan(1, nCond);
            simFields = fieldnames(pr0.similarity);
            for sf = 1:numel(simFields)
                prSum.similarity.(simFields{sf}) = nan(1, nCond);
            end

            for c = 1:nCond
                thisSim = AllSubspaceSim{c};
                if numel(thisSim.group) < g || numel(thisSim.group(g).pair) < p
                    error('Condition %d is missing group %d pair %d.', c, g, p);
                end
                pr = thisSim.group(g).pair{p};
                prSum.status{c} = pr.status;
                prSum.warning{c} = pr.warning;
                prSum.idxA{c} = pr.idxA;
                prSum.idxB{c} = pr.idxB;
                prSum.rawDimA(c) = pr.rawDimA;
                prSum.rawDimB(c) = pr.rawDimB;
                prSum.basisA{c} = pr.basisA;
                prSum.basisB{c} = pr.basisB;
                prSum.principal.cosine{c} = pr.principal.cosine;
                prSum.principal.angle_deg{c} = pr.principal.angle_deg;
                prSum.principal.first_angle_deg(c) = pr.principal.first_angle_deg;
                prSum.principal.last_angle_deg(c) = pr.principal.last_angle_deg;
                prSum.principal.median_angle_deg(c) = pr.principal.median_angle_deg;
                for sf = 1:numel(simFields)
                    simField = simFields{sf};
                    if isfield(pr.similarity, simField) && isnumeric(pr.similarity.(simField)) && isscalar(pr.similarity.(simField))
                        prSum.similarity.(simField)(c) = pr.similarity.(simField);
                    else
                        prSum.similarity.(simField)(c) = NaN;
                    end
                end
            end

            Summary.group(g).pair{p} = prSum;
        end
    end
end

function classification = buildSummaryClassificationLocal(AllSubspaceSim, condition_list, stim_abbrev)
    nCond = numel(AllSubspaceSim);
    class0 = AllSubspaceSim{1}.classification;

    classification = struct();
    if isfield(class0, 'categoryLabels')
        classification.categoryLabels = class0.categoryLabels;
    end
    classification.condition_id = condition_list;
    classification.stim_abbrev = stim_abbrev;

    fieldsToCollect = {'acrossDelay', 'acrossIdx', 'feedforwardIdx', 'feedbackIdx', 'ambiguousIdx'};
    for f = 1:numel(fieldsToCollect)
        fieldName = fieldsToCollect{f};
        classification.(fieldName) = cell(1, nCond);
        for c = 1:nCond
            if isfield(AllSubspaceSim{c}.classification, fieldName)
                classification.(fieldName){c} = AllSubspaceSim{c}.classification.(fieldName);
            else
                classification.(fieldName){c} = [];
            end
        end
    end
end

function tableCell = buildConditionGroupTableCellLocal(SubspaceSim, groupIdx, stim_abbrev)
    nRows = numel(stim_abbrev);
    tableCell = makeSummaryHeaderCellLocal('Stim');
    for r = 1:nRows
        rowCell = makeOneTableRowFromGroupLocal(SubspaceSim.group(groupIdx), r, stim_abbrev{r});
        tableCell(end+1, :) = rowCell; %#ok<AGROW>
    end
end

function tableCell = buildPooledGroupTableCellLocal(SubspaceSim)
    nGroups = numel(SubspaceSim.group);
    tableCell = makeSummaryHeaderCellLocal('Group');
    for g = 1:nGroups
        rowLabel = sprintf('Group %d', g);
        rowCell = makeOneTableRowFromGroupLocal(SubspaceSim.group(g), 1, rowLabel);
        tableCell(end+1, :) = rowCell; %#ok<AGROW>
    end
end

function tableCell = makeSummaryHeaderCellLocal(firstColName)
    header1 = {firstColName, 'Across dim', 'Within dim', ...
        'Across vs Within angle', '', '', ...
        'Across vs Within overlap', '', '', ...
        'FF dim', 'FB dim', ...
        'Feedforward vs Feedback angle', '', '', ...
        'Feedforward vs Feedback overlap', '', ''};

    header2 = {'', '', '', ...
        'First', 'Last', 'Median', ...
        'Across captures Within', 'Within captures Across', 'Avg', ...
        '', '', ...
        'First', 'Last', 'Median', ...
        'Feedforward captures Feedback', 'Feedback captures Feedforward', 'Avg'};

    tableCell = [header1; header2];
end

function rowCell = makeOneTableRowFromGroupLocal(groupStruct, idx, rowLabel)
    avw = getPairByNameLocal(groupStruct, 'across_vs_within');
    fvf = getPairByNameLocal(groupStruct, 'feedforward_vs_feedback');

    rowCell = {rowLabel, ...
        formatDimLocal(getVectorValueLocal(avw.rawDimA, idx)), ...
        formatDimLocal(getVectorValueLocal(avw.rawDimB, idx)), ...
        formatMetricLocal(getVectorValueLocal(avw.principal.first_angle_deg, idx)), ...
        formatMetricLocal(getVectorValueLocal(avw.principal.last_angle_deg, idx)), ...
        formatMetricLocal(getVectorValueLocal(avw.principal.median_angle_deg, idx)), ...
        formatMetricLocal(getCaptureValueLocal(avw, idx, 'Across', 'Within')), ...
        formatMetricLocal(getCaptureValueLocal(avw, idx, 'Within', 'Across')), ...
        formatMetricLocal(getVectorValueLocal(avw.similarity.avg, idx)), ...
        formatDimLocal(getVectorValueLocal(fvf.rawDimA, idx)), ...
        formatDimLocal(getVectorValueLocal(fvf.rawDimB, idx)), ...
        formatMetricLocal(getVectorValueLocal(fvf.principal.first_angle_deg, idx)), ...
        formatMetricLocal(getVectorValueLocal(fvf.principal.last_angle_deg, idx)), ...
        formatMetricLocal(getVectorValueLocal(fvf.principal.median_angle_deg, idx)), ...
        formatMetricLocal(getCaptureValueLocal(fvf, idx, 'Feedforward', 'Feedback')), ...
        formatMetricLocal(getCaptureValueLocal(fvf, idx, 'Feedback', 'Feedforward')), ...
        formatMetricLocal(getVectorValueLocal(fvf.similarity.avg, idx))};
end

function pr = getPairByNameLocal(groupStruct, pairName)
    pr = [];
    for p = 1:numel(groupStruct.pair)
        thisPr = groupStruct.pair{p};
        if isfield(thisPr, 'name') && strcmp(thisPr.name, pairName)
            pr = thisPr;
            return;
        end
    end
    error('Could not find pair %s in %s.', pairName, groupStruct.name);
end

function val = getCaptureValueLocal(pairStruct, idx, labelA, labelB)
    fieldName = makeCaptureFieldNameLocal(labelA, labelB);
    if isfield(pairStruct, 'similarity') && isfield(pairStruct.similarity, fieldName)
        val = getVectorValueLocal(pairStruct.similarity.(fieldName), idx);
    else
        val = NaN;
    end
end

function val = getVectorValueLocal(x, idx)
    if iscell(x)
        if idx <= numel(x)
            val = x{idx};
        else
            val = NaN;
        end
        return;
    end

    if isempty(x)
        val = NaN;
    elseif isscalar(x)
        val = x;
    elseif idx <= numel(x)
        val = x(idx);
    else
        val = NaN;
    end
end

function s = formatDimLocal(v)
    if isempty(v) || ~isnumeric(v) || ~isscalar(v) || ~isfinite(v)
        s = 'NaN';
    else
        s = sprintf('%d', round(v));
    end
end

function s = formatMetricLocal(v)
    if isempty(v) || ~isnumeric(v) || ~isscalar(v) || ~isfinite(v)
        s = 'NaN';
    else
        s = sprintf('%.4f', v);
    end
end

function printCellTableLocal(fid, tableCell, sepAfter)
    if nargin < 3
        sepAfter = [];
    end

    strCell = cell(size(tableCell));
    for r = 1:size(tableCell, 1)
        for c = 1:size(tableCell, 2)
            strCell{r, c} = char(string(tableCell{r, c}));
        end
    end

    widths = zeros(1, size(strCell, 2));
    for c = 1:size(strCell, 2)
        maxWidth = 1;
        for r = 1:size(strCell, 1)
            maxWidth = max(maxWidth, length(strCell{r, c}));
        end
        widths(c) = maxWidth;
    end

    totalWidth = sum(widths + 2) + 3 * numel(sepAfter);

    for r = 1:size(strCell, 1)
        for c = 1:size(strCell, 2)
            if c == 1
                fprintf(fid, '%-*s', widths(c), strCell{r, c});
            else
                fprintf(fid, '  %*s', widths(c), strCell{r, c});
            end
            if ismember(c, sepAfter)
                fprintf(fid, ' |');
            end
        end
        fprintf(fid, '\n');
        if r == 2
            fprintf(fid, '%s\n', repmat('-', 1, totalWidth));
        end
    end
end

function saveTableTextAndCsvLocal(tableCell, txtFile, csvFile, titleStr, sepAfter)
    fid = fopen(txtFile, 'w');
    if fid < 0
        error('Could not open %s for writing.', txtFile);
    end
    cleanupObj = onCleanup(@() fclose(fid));
    fprintf(fid, '%s\n', titleStr);
    printCellTableLocal(fid, tableCell, sepAfter);
    clear cleanupObj;

    writeCellCsvLocal(tableCell, csvFile);
end

function writeCellCsvLocal(C, filename)
    fid = fopen(filename, 'w');
    if fid < 0
        error('Could not open %s for writing.', filename);
    end
    cleanupObj = onCleanup(@() fclose(fid));

    for r = 1:size(C, 1)
        for c = 1:size(C, 2)
            txt = char(string(C{r, c}));
            txt = strrep(txt, '"', '""');
            if contains(txt, ',') || contains(txt, '"') || contains(txt, newline)
                txt = ['"', txt, '"']; %#ok<AGROW>
            end
            if c > 1
                fprintf(fid, ',');
            end
            fprintf(fid, '%s', txt);
        end
        fprintf(fid, '\n');
    end
    clear cleanupObj;
end

function printSubspaceSimilarityResultsLocal(SubspaceSim)
    for g = 1:numel(SubspaceSim.group)
        fprintf('\n%s\n', SubspaceSim.group(g).name);
        for p = 1:numel(SubspaceSim.group(g).pair)
            pr = SubspaceSim.group(g).pair{p};
            fprintf('  %s: %s vs %s\n', pr.name, pr.labelA, pr.labelB);
            fprintf('    dims: %s=%d, %s=%d\n', ...
                pr.labelA, pr.rawDimA, pr.labelB, pr.rawDimB);
            if strcmp(pr.status, 'ok')
                fprintf('    principal angle first / last / median deg = %.4f / %.4f / %.4f\n', ...
                    pr.principal.first_angle_deg, ...
                    pr.principal.last_angle_deg, ...
                    pr.principal.median_angle_deg);
                fprintf('    subspace overlap: %s captures %s = %.4f; %s captures %s = %.4f; avg = %.4f\n', ...
                    pr.labelA, pr.labelB, getCaptureValueLocal(pr, 1, pr.labelA, pr.labelB), ...
                    pr.labelB, pr.labelA, getCaptureValueLocal(pr, 1, pr.labelB, pr.labelA), ...
                    pr.similarity.avg);
            else
                fprintf('    skipped: %s\n', pr.warning);
            end
        end
    end
end

function fname = findOneFileLocal(parentDir, pattern, mustExist)
    if ~exist(parentDir, 'dir')
        error('Directory not found: %s', parentDir);
    end

    files = dir(fullfile(parentDir, pattern));
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

function tag = makeLatentSelectionTagLocal(use_dsl_filter, dsl_field)
    if use_dsl_filter
        fieldTag = char(string(dsl_field));
        fieldTag = regexprep(fieldTag, '[^A-Za-z0-9]+', '_');
        fieldTag = regexprep(fieldTag, '^_+|_+$', '');
        tag = sprintf('DSL_%s_filtered', fieldTag);
    else
        tag = 'all_latents';
    end
end

function displayName = makeLatentSelectionDisplayLocal(use_dsl_filter, dsl_field)
    if use_dsl_filter
        displayName = sprintf('DSL %s filtered', char(string(dsl_field)));
    else
        displayName = 'all latents';
    end
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
            error('Unsupported stim_name in condition_full(%d): %s', k, char(currStim));
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
            error('Stim %s does not have exactly 2 contrast levels.', char(stimLabels(s)));
        end
        contrastValuesByStim.(char(stimLabels(s))) = cvals;
    end

    dirVals = unique(effDirAll);
    dirVals = sort(dirVals(:)');
    if numel(dirVals) ~= 2
        error('Expected exactly 2 effective direction values in condition_full.');
    end

    stimDirLabels = {'stim_dir1', 'stim_dir2'};
    condLabels = { ...
        'grating-small-low', 'grating-small-high', ...
        'grating-large-low', 'grating-large-high', ...
        'plaid-small-low',   'plaid-small-high', ...
        'plaid-large-low',   'plaid-large-high'};
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
