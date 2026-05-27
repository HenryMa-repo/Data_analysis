%% data_reconstruction.m
% Post-hoc DLAG data reconstruction from already saved bestmodel*.mat files.
%
% Run this script after plot_dlag_results.m has saved the best-model mat file.
%
% This script:
%   1) Does not refit the model.
%   2) Does not rerun inference.
%   3) Uses the existing seqEst.xsm in the saved bestmodel*.mat file.
%   4) Adds reconstruction fields to seqEst.
%   5) Overwrites the original bestmodel*.mat with the updated seqEst.
%   6) Saves R2 results separately as reconstruction_R2.mat.
%
% Reconstruction fields added to each seqEst(n):
%   seqEst(n).yRecon_use_across
%   seqEst(n).yRecon_use_within
%   seqEst(n).yRecon_use_all
%   seqEst(n).yRecon_across_excl_within
%   seqEst(n).yRecon_within_excl_across
%
% R2 output:
%   reconstruction_R2.mat contains recon_R2, with fields:
%     recon_R2.use_across.global_all
%     recon_R2.use_across.global_by_group
%     recon_R2.use_across.neuron_by_group
%     recon_R2.use_within.global_all
%     recon_R2.use_within.global_by_group
%     recon_R2.use_within.neuron_by_group
%     recon_R2.use_all.global_all
%     recon_R2.use_all.global_by_group
%     recon_R2.use_all.neuron_by_group
%     recon_R2.across_excl_within.global_all
%     recon_R2.across_excl_within.global_by_group
%     recon_R2.across_excl_within.neuron_by_group
%     recon_R2.within_excl_across.global_all
%     recon_R2.within_excl_across.global_by_group
%     recon_R2.within_excl_across.neuron_by_group
%
% Notes:
%   - Reconstruction uses internally orthogonalized loading blocks.
%   - SVD uses svd(...,'econ').
%   - For nonzero blocks, dimensions are taken as xDim, not numerical rank.
%   - Completely zero loading blocks are treated as empty.
%   - Across/within 0-dimensional cases are allowed.
%   - R2 values are allowed to be negative and are not clipped.

clc;
clear;

%% ------------------------------------------------------------------------
% User parameters
% -------------------------------------------------------------------------

data_content = 'demean_count_within_trial';
% options usually include:
% raw_count, raw_fr, z_within_trial, z_within_condition,
% z_across_conditions, demean_count_within_trial, demean_fr_within_trial,
% demean_pooledsd_within_condition

data_condition = [1:16];
% [] for pooled all-condition mode, or e.g. 1:16 for condition mode.

runIdx = 1;

%% ------------------------------------------------------------------------
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

%% ------------------------------------------------------------------------
% Main loop
% -------------------------------------------------------------------------

for cond_i = 1:numConditions

    if use_condition_mode
        this_condition = condition_list(cond_i);
        baseDir = ['./FA_Dlag_', data_content, '_condition', num2str(this_condition)];
    else
        this_condition = [];
        baseDir = ['./FA_Dlag_', data_content];
    end

    tempfname = sprintf('%s/mat_results/run%03d', baseDir, runIdx);

    fprintf('\n============================================================\n');
    if isempty(this_condition)
        fprintf('DLAG data reconstruction: pooled all-condition mode\n');
    else
        fprintf('DLAG data reconstruction: condition %d\n', this_condition);
    end
    fprintf('Reading from %s\n', tempfname);

    bestFile = findOneFileLocal(tempfname, 'bestmodel*', true);
    fprintf('Loading best model: %s\n', bestFile);

    Sbest = load(bestFile);

    requiredBestVars = {'bestModel', 'res', 'seqEst'};
    for v = 1:numel(requiredBestVars)
        if ~isfield(Sbest, requiredBestVars{v})
            error('File %s is missing variable %s.', bestFile, requiredBestVars{v});
        end
    end

    bestModel = Sbest.bestModel;
    res = Sbest.res;
    seqEst = Sbest.seqEst;

    if ~isfield(res, 'estParams')
        error('File %s is missing res.estParams.', bestFile);
    end

    params = res.estParams;
    params = normalizeParamDimsLocal(params, bestModel);

    if ~isfield(seqEst, 'xsm')
        error('seqEst.xsm not found in %s. This script expects saved inference results.', bestFile);
    end

    if ~isfield(seqEst, 'y')
        error('seqEst.y not found in %s.', bestFile);
    end

    fprintf('Adding reconstruction fields to seqEst...\n');
    seqEst = addDlagReconstructionFieldsLocal(seqEst, params);

    fprintf('Computing reconstruction R2...\n');
    recon_R2 = computeReconstructionR2Local(seqEst, params.yDims);

    Sbest.seqEst = seqEst;

    fprintf('Overwriting best-model mat with augmented seqEst...\n');
    save(bestFile, '-struct', 'Sbest', '-v7.3');

    r2File = fullfile(fileparts(bestFile), 'reconstruction_R2.mat');
    fprintf('Saving R2 results: %s\n', r2File);
    save(r2File, 'recon_R2', '-v7.3');

    fprintf('Done.\n');
end

%% ========================================================================
% Local functions
% ========================================================================

function seqEst = addDlagReconstructionFieldsLocal(seqEst, params)
% Add the five requested reconstruction fields to each trial in seqEst.

params = normalizeParamDimsLocal(params, []);

yDims = params.yDims;
numGroups = numel(yDims);
yDim = sum(yDims);

xDim_across = params.xDim_across;
xDim_within = params.xDim_within;
localDims = xDim_across + xDim_within;
xDim_total = sum(localDims);

if size(params.C, 1) ~= yDim
    error('params.C has %d rows, expected sum(params.yDims) = %d.', ...
        size(params.C, 1), yDim);
end

if size(params.C, 2) ~= xDim_total
    error(['params.C has %d columns, expected sum(xDim_across + xDim_within) ', ...
        '= %d.'], size(params.C, 2), xDim_total);
end

if numel(params.d) ~= yDim
    error('params.d length %d does not match sum(params.yDims) = %d.', ...
        numel(params.d), yDim);
end

params.d = params.d(:);

blocks = precomputeReconstructionBlocksLocal(params);

for n = 1:numel(seqEst)

    if isempty(seqEst(n).xsm)
        error('seqEst(%d).xsm is empty.', n);
    end

    if isempty(seqEst(n).y)
        error('seqEst(%d).y is empty.', n);
    end

    if size(seqEst(n).xsm, 1) ~= xDim_total
        error('seqEst(%d).xsm has %d rows, expected %d.', ...
            n, size(seqEst(n).xsm, 1), xDim_total);
    end

    T = size(seqEst(n).xsm, 2);

    if size(seqEst(n).y, 1) ~= yDim
        error('seqEst(%d).y has %d rows, expected %d.', ...
            n, size(seqEst(n).y, 1), yDim);
    end

    if size(seqEst(n).y, 2) ~= T
        error('seqEst(%d).y and seqEst(%d).xsm have different time lengths.', n, n);
    end

    yBase = repmat(params.d, 1, T);

    yRecon_use_across = yBase;
    yRecon_use_within = yBase;
    yRecon_use_all = yBase;
    yRecon_across_excl_within = yBase;
    yRecon_within_excl_across = yBase;

    for groupIdx = 1:numGroups

        rows = blocks(groupIdx).obsIdx;
        d_g = repmat(params.d(rows), 1, T);

        X_across = seqEst(n).xsm(blocks(groupIdx).latIdx_across, :);
        X_within = seqEst(n).xsm(blocks(groupIdx).latIdx_within, :);
        X_all = seqEst(n).xsm(blocks(groupIdx).latIdx_all, :);

        Y_across = reconstructFromBlockLocal( ...
            blocks(groupIdx).Q_across, ...
            blocks(groupIdx).TT_across, ...
            X_across, ...
            numel(rows), T);

        Y_within = reconstructFromBlockLocal( ...
            blocks(groupIdx).Q_within, ...
            blocks(groupIdx).TT_within, ...
            X_within, ...
            numel(rows), T);

        Y_all = reconstructFromBlockLocal( ...
            blocks(groupIdx).Q_all, ...
            blocks(groupIdx).TT_all, ...
            X_all, ...
            numel(rows), T);

        Y_across_excl_within = projectOrthogonalComplementLocal( ...
            Y_across, blocks(groupIdx).Q_within);

        Y_within_excl_across = projectOrthogonalComplementLocal( ...
            Y_within, blocks(groupIdx).Q_across);

        yRecon_use_across(rows, :) = d_g + Y_across;
        yRecon_use_within(rows, :) = d_g + Y_within;
        yRecon_use_all(rows, :) = d_g + Y_all;
        yRecon_across_excl_within(rows, :) = d_g + Y_across_excl_within;
        yRecon_within_excl_across(rows, :) = d_g + Y_within_excl_across;
    end

    seqEst(n).yRecon_use_across = yRecon_use_across;
    seqEst(n).yRecon_use_within = yRecon_use_within;
    seqEst(n).yRecon_use_all = yRecon_use_all;
    seqEst(n).yRecon_across_excl_within = yRecon_across_excl_within;
    seqEst(n).yRecon_within_excl_across = yRecon_within_excl_across;
end
end

function blocks = precomputeReconstructionBlocksLocal(params)
% Precompute per-group loading bases and latent transforms.

yDims = params.yDims;
xDim_across = params.xDim_across;
xDim_within = params.xDim_within;

numGroups = numel(yDims);
localDims = xDim_across + xDim_within;

obsStart = cumsum([1, yDims(1:end-1)]);
obsEnd = cumsum(yDims);

latStart = cumsum([1, localDims(1:end-1)]);
latEnd = cumsum(localDims);

blocks = struct([]);

for groupIdx = 1:numGroups

    obsIdx = obsStart(groupIdx):obsEnd(groupIdx);
    latIdx = latStart(groupIdx):latEnd(groupIdx);

    localDim_g = localDims(groupIdx);

    acrossLocal = 1:xDim_across;
    withinLocal = (xDim_across + 1):localDim_g;

    latIdx_across = latIdx(acrossLocal);
    latIdx_within = latIdx(withinLocal);
    latIdx_all = latIdx;

    Cg = params.C(obsIdx, latIdx_all);

    C_across = Cg(:, acrossLocal);
    C_within = Cg(:, withinLocal);
    C_all = Cg;

    [Q_across, TT_across] = orthogonalizeLoadingBlockLocal(C_across);
    [Q_within, TT_within] = orthogonalizeLoadingBlockLocal(C_within);
    [Q_all, TT_all] = orthogonalizeLoadingBlockLocal(C_all);

    blocks(groupIdx).obsIdx = obsIdx;

    blocks(groupIdx).latIdx_across = latIdx_across;
    blocks(groupIdx).latIdx_within = latIdx_within;
    blocks(groupIdx).latIdx_all = latIdx_all;

    blocks(groupIdx).Q_across = Q_across;
    blocks(groupIdx).TT_across = TT_across;

    blocks(groupIdx).Q_within = Q_within;
    blocks(groupIdx).TT_within = TT_within;

    blocks(groupIdx).Q_all = Q_all;
    blocks(groupIdx).TT_all = TT_all;

    fprintf('  Group %d: xAcross=%d, xWithin=%d, basisAcross=%d, basisWithin=%d, basisAll=%d\n', ...
        groupIdx, xDim_across, xDim_within(groupIdx), ...
        size(Q_across, 2), size(Q_within, 2), size(Q_all, 2));
end
end

function Y = reconstructFromBlockLocal(Q, TT, X, yDim_group, T)
% Return Q * TT * X, with clean behavior for 0-dimensional blocks.

if isempty(Q) || isempty(TT) || isempty(X) || size(Q, 2) == 0
    Y = zeros(yDim_group, T);
    return;
end

Y = Q * (TT * X);
end

function Y_resid = projectOrthogonalComplementLocal(Y, Q_remove)
% Project Y into the orthogonal complement of span(Q_remove).
%
% Since Q_remove is an orthonormal SVD basis:
%   P = Q_remove * Q_remove'
%   Y_resid = (I - P) * Y
%           = Y - Q_remove * (Q_remove' * Y)
%
% This is the same null-space residual idea as the subspace-overlap code,
% but here the object being projected is reconstruction data rather than
% another loading subspace basis.

if isempty(Q_remove) || size(Q_remove, 2) == 0
    Y_resid = Y;
else
    Y_resid = Y - Q_remove * (Q_remove' * Y);
end
end

function [Q, TT] = orthogonalizeLoadingBlockLocal(L)
% Internal orthogonalization of one loading block.
%
% For L with size yDim x xDim, return Q and TT such that:
%   L * X == Q * (TT * X)
% up to numerical precision.
%
% For nonzero blocks, the basis dimension is xDim, matching the DLAG-style
% SVD basis convention. If the whole block is exactly/effectively zero, the
% contribution is treated as empty.
%
% svd(...,'econ') is used because these loading blocks are expected to be
% tall matrices.

[yDim, xDim] = size(L);

if xDim == 0
    Q = zeros(yDim, 0);
    TT = zeros(0, 0);
    return;
end

if xDim == 1
    mag = sqrt(L' * L);

    if mag <= eps
        Q = zeros(yDim, 0);
        TT = zeros(0, 1);
    else
        Q = L / mag;
        TT = mag;
    end

    return;
end

[U, S, V] = svd(L, 'econ');
s = diag(S);

if isempty(s) || max(s) <= eps
    Q = zeros(yDim, 0);
    TT = zeros(0, xDim);
    return;
end

if size(U, 2) < xDim || size(S, 1) < xDim || size(V, 2) < xDim
    error(['SVD returned fewer columns than xDim. This script expects tall ', ...
        'or at least not rank-shape-conflicting loading blocks. Got yDim=%d, xDim=%d.'], ...
        yDim, xDim);
end

Q = U(:, 1:xDim);
TT = S(1:xDim, 1:xDim) * V(:, 1:xDim)';
end

function recon_R2 = computeReconstructionR2Local(seqEst, yDims)
% Compute global, group-wise global, and neuron-wise R2 for each recon field.

yDims = reshape(yDims, 1, []);
numGroups = numel(yDims);

r2_specs = {
    'use_across',          'yRecon_use_across';
    'use_within',          'yRecon_use_within';
    'use_all',             'yRecon_use_all';
    'across_excl_within',  'yRecon_across_excl_within';
    'within_excl_across',  'yRecon_within_excl_across'};

Ytrue = [seqEst.y];

recon_R2 = struct();

for specIdx = 1:size(r2_specs, 1)

    r2Name = r2_specs{specIdx, 1};
    fieldName = r2_specs{specIdx, 2};

    if ~isfield(seqEst, fieldName)
        error('Field %s is missing from seqEst.', fieldName);
    end

    Ypred = [seqEst.(fieldName)];

    recon_R2.(r2Name).global_all = computeGlobalR2Local(Ytrue, Ypred);
    recon_R2.(r2Name).global_by_group = nan(1, numGroups);
    recon_R2.(r2Name).neuron_by_group = cell(1, numGroups);

    for groupIdx = 1:numGroups
        rows = getGroupRowsLocal(yDims, groupIdx);

        recon_R2.(r2Name).global_by_group(groupIdx) = ...
            computeGlobalR2Local(Ytrue(rows, :), Ypred(rows, :));

        recon_R2.(r2Name).neuron_by_group{groupIdx} = ...
            computeNeuronR2Local(Ytrue(rows, :), Ypred(rows, :));
    end
end
end

function R2 = computeGlobalR2Local(Ytrue, Ypred)
% Global R2 with per-neuron mean baseline.
% Negative values are allowed.

if ~isequal(size(Ytrue), size(Ypred))
    error('Ytrue and Ypred must have the same size.');
end

valid = isfinite(Ytrue) & isfinite(Ypred);
numValid = sum(valid, 2);

Ytmp = Ytrue;
Ytmp(~valid) = 0;

mu = sum(Ytmp, 2) ./ numValid;
mu(numValid == 0) = NaN;

D = Ytrue - repmat(mu, 1, size(Ytrue, 2));
E = Ytrue - Ypred;

D(~valid) = 0;
E(~valid) = 0;

RSS = sum(E(:).^2);
TSS = sum(D(:).^2);

if TSS > 0 && isfinite(TSS)
    R2 = 1 - RSS / TSS;
else
    R2 = NaN;
end
end

function R2 = computeNeuronR2Local(Ytrue, Ypred)
% Per-neuron R2 with each neuron centered by its own mean.
% Negative values are allowed.

if ~isequal(size(Ytrue), size(Ypred))
    error('Ytrue and Ypred must have the same size.');
end

numNeurons = size(Ytrue, 1);
R2 = nan(numNeurons, 1);

for i = 1:numNeurons

    valid = isfinite(Ytrue(i, :)) & isfinite(Ypred(i, :));

    if ~any(valid)
        continue;
    end

    yt = Ytrue(i, valid);
    yp = Ypred(i, valid);

    mu = mean(yt);

    RSS = sum((yt - yp).^2);
    TSS = sum((yt - mu).^2);

    if TSS > 0 && isfinite(TSS)
        R2(i) = 1 - RSS / TSS;
    end
end
end

function params = normalizeParamDimsLocal(params, bestModel)
% Normalize params.yDims, params.xDim_across, and params.xDim_within.
% If params lacks xDim fields, use bestModel.

if ~isfield(params, 'C')
    error('params must contain loading matrix C.');
end

if ~isfield(params, 'd')
    error('params must contain baseline vector d.');
end

if ~isfield(params, 'yDims')
    error('params must contain yDims.');
end

params.yDims = reshape(params.yDims, 1, []);
numGroups = numel(params.yDims);

if ~isfield(params, 'xDim_across') || isempty(params.xDim_across)
    if ~isempty(bestModel) && isfield(bestModel, 'xDim_across')
        params.xDim_across = bestModel.xDim_across;
    else
        params.xDim_across = 0;
    end
end

if isscalar(params.xDim_across)
    params.xDim_across = double(params.xDim_across);
else
    error('params.xDim_across must be scalar for standard DLAG.');
end

if ~isfield(params, 'xDim_within') || isempty(params.xDim_within)
    if ~isempty(bestModel) && isfield(bestModel, 'xDim_within')
        params.xDim_within = bestModel.xDim_within;
    else
        params.xDim_within = zeros(1, numGroups);
    end
end

params.xDim_within = reshape(params.xDim_within, 1, []);

if isscalar(params.xDim_within) && numGroups > 1
    params.xDim_within = repmat(params.xDim_within, 1, numGroups);
end

if numel(params.xDim_within) ~= numGroups
    error('xDim_within must have one entry per group.');
end

params.xDim_within = double(params.xDim_within);
end

function rows = getGroupRowsLocal(yDims, groupIdx)
% Row indices for one observed group.

startIdx = sum(yDims(1:groupIdx-1)) + 1;
rows = startIdx:(startIdx + yDims(groupIdx) - 1);
end

function fname = findOneFileLocal(parentDir, pattern, mustExist)
% Return the newest file matching pattern in parentDir.

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