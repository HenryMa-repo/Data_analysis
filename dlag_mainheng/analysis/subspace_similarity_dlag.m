%% subspace_similarity_dlag.m
% Compute within-group DLAG subspace principal angles and subspace overlap.
%
% Run this after plot_dlag_results.m has saved the bestmodel* file in
%   ./FA_Dlag_<data_content>/mat_results/run%03d
% or, for condition-specific models,
%   ./FA_Dlag_<data_content>_condition<condition>/mat_results/run%03d
%
% The program follows the same loading pattern as Latents_compare.m and uses
% DLAG-style SVD bases, following orthogonalize.m / dominantModes_dlag.m,
% instead of MATLAB orth().
%
% It computes, for every neural group:
%   1) across subspace vs within subspace
%   2) feedforward across subspace vs feedback across subspace
%
% For each pair, it stores:
%   - all principal cosines
%   - all principal angles in degrees
%   - median angle in degrees
%   - first and last principal angle in degrees
%   - directional subspace overlap A captures B and B captures A
%   - average of the two directional overlaps
%
% Subspace overlap formula:
%   S(U,V) = 1 - ||(I - U*(U'*U)^(-1)*U')*V||_F / ||V||_F
% This is directional, so the program computes both S(A,B) and S(B,A).
% A captures B means B lies in A; equivalently, B projects little to A's null space.

clc; clear;


% -------------------------------------------------------------------------
data_content = 'z_across_conditions';
% options:
% raw_count, raw_fr, z_within_trial, z_within_condition,
% z_across_conditions, demean_count_within_trial, demean_fr_within_trial,
% demean_pooledsd_within_condition

data_condition = [];     % [] for pooled all-condition model, or e.g. [1 2 3]
runIdx = 1;

% Whether to remove latents marked as DSL-remove by Latents_compare.m.
% Default is false because this program is intended to be runnable right
% after plot_dlag_results.m, before Latents_compare.m.
use_dsl_filter = false;
dsl_field = 'logical';    % if use_dsl_filter=true, use DSL.(dsl_field){groupIdx}

% Save / print options
save_results = true;
print_results = true;

% -------------------------------------------------------------------------
% Main loop: pooled mode or condition-specific mode
% -------------------------------------------------------------------------
if isempty(data_condition)
    condition_list = [];
    numConditions = 1;
else
    condition_list = data_condition(:)';
    numConditions = numel(condition_list);
end

for cond_i = 1:numConditions
    if isempty(condition_list)
        this_condition = [];
        baseDir = ['./FA_Dlag_', data_content];
    else
        this_condition = condition_list(cond_i);
        baseDir = ['./FA_Dlag_', data_content, '_condition', num2str(this_condition)];
    end

    tempfname = sprintf('%s/mat_results/run%03d', baseDir, runIdx);

    fprintf('\n============================================================\n');
    if isempty(this_condition)
        fprintf('DLAG subspace similarity: pooled model\n');
    else
        fprintf('DLAG subspace similarity: condition %d\n', this_condition);
    end
    fprintf('Reading from %s\n', tempfname);

    bestFile = findOneFile(tempfname, 'bestmodel*', true);
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
    bootFile = findOneFile(tempfname, 'bootstrapResults*', false);
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

    SubspaceSim = computeDlagSubspaceSimilarity(bestModel, params, gp_params, ambiguousIdxs, DSL, opts);

    if print_results
        printSubspaceSimilarityResults(SubspaceSim);
    end

    if save_results
        outFile = fullfile(tempfname, 'subspace_similarity_results.mat');
        save(outFile, 'SubspaceSim', 'use_dsl_filter', 'dsl_field');
        fprintf('Saved %s\n', outFile);
    end
end

%% ========================================================================
% Local functions
%% ========================================================================

function SubspaceSim = computeDlagSubspaceSimilarity(bestModel, params, gp_params, ambiguousIdxs, DSL, opts)

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
    SubspaceSim.meta.use_dsl_filter = opts.use_dsl_filter;
    SubspaceSim.meta.dsl_field = opts.dsl_field;
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

            % Use a cell array so skipped pairs do not need fake NaN fields.
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

    if isempty(QA) || isempty(QB) || size(QA, 2) == 0 || size(QB, 2) == 0
        pairResult.status = 'skipped_empty_subspace';
        pairResult.warning = sprintf('%s or %s has zero usable dimension.', labelA, labelB);
        return;
    end

    % Principal angles:
    % Singular values of QA' * QB are principal cosines.
    s = svd(QA' * QB, 'econ');

    % Numerical protection before acos. SVD values should be in [0, 1].
    s = min(max(s, 0), 1);

    thetaRad = acos(s);
    thetaDeg = thetaRad * 180 / pi;

    pairResult.principal.cosine = reshape(s, 1, []);
    pairResult.principal.angle_deg = reshape(thetaDeg, 1, []);
    pairResult.principal.first_angle_deg = thetaDeg(1);
    pairResult.principal.last_angle_deg = thetaDeg(end);
    pairResult.principal.median_angle_deg = median(thetaDeg);

    % Directional subspace overlap.
    % Computed twice because the measure is not symmetric.
    pairResult.similarity.A_captures_B = directionalSubspaceOverlapLocal(QA, QB);
    pairResult.similarity.B_captures_A = directionalSubspaceOverlapLocal(QB, QA);

    simVals = [pairResult.similarity.A_captures_B, pairResult.similarity.B_captures_A];
    simVals = simVals(~isnan(simVals));

    if isempty(simVals)
        pairResult.similarity.avg = NaN;
    else
        pairResult.similarity.avg = mean(simVals);
    end
end

function overlap = directionalSubspaceOverlapLocal(U, V)
    % S(U,V) = 1 - ||(I - P_U)V||_F / ||V||_F
    % P_U = U * inv(U' * U) * U'
    %
    % U and V are already DLAG-style SVD bases, so U'*U should be close to I.
    % The formula is kept in this explicit form to match the subspace
    % similarity definition.

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

    % Numerical protection against tiny negative values caused by floating-point error.
    if overlap < 0 && overlap > -1e-12
        overlap = 0;
    end
end

function Q = dlagSvdBasisLocal(L)
    % Build an orthonormal basis using the same SVD idea as DLAG's
    % orthogonalize.m and dominantModes_dlag.m.
    %
    % Important:
    %   - This function intentionally does not call MATLAB orth().
    %   - It does not estimate effective rank.
    %   - It takes the first xDim left singular vectors, where xDim=size(L,2).

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
    %
    % This follows the same idea as your Latents_compare logic:
    %   positive delay -> feedforward
    %   negative delay -> feedback
    %   zero / NaN / bootstrap-ambiguous -> ambiguous
    %
    % The exact biological interpretation of delay sign should match how
    % gp_params.delays was stored by plot_dlag_results.m.

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
    % Make ambiguousIdxs robust to empty, vector, logical mask, or cell input.

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

function printSubspaceSimilarityResults(SubspaceSim)

    for g = 1:numel(SubspaceSim.group)
        fprintf('\n%s\n', SubspaceSim.group(g).name);

        for p = 1:numel(SubspaceSim.group(g).pair)
            pr = SubspaceSim.group(g).pair{p};

            fprintf('  %s: %s vs %s\n', pr.name, pr.labelA, pr.labelB);
            fprintf('    dims: %s=%d, %s=%d\n', ...
                pr.labelA, pr.rawDimA, ...
                pr.labelB, pr.rawDimB);

            if strcmp(pr.status, 'ok')
                fprintf('    principal angle first / last / median deg = %.4f / %.4f / %.4f\n', ...
                    pr.principal.first_angle_deg, ...
                    pr.principal.last_angle_deg, ...
                    pr.principal.median_angle_deg);

                fprintf('    subspace overlap: %s captures %s = %.4f; %s captures %s = %.4f; avg = %.4f\n', ...
                    pr.labelA, pr.labelB, pr.similarity.A_captures_B, ...
                    pr.labelB, pr.labelA, pr.similarity.B_captures_A, ...
                    pr.similarity.avg);
            else
                fprintf('    skipped: %s\n', pr.warning);
            end
        end
    end
end

function fname = findOneFile(parentDir, pattern, mustExist)

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