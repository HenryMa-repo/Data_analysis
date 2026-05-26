%% pick_latents_by_svexp.m
% Pick DLAG latents by cumulative shared variance explained.
%
% For each neural group, latents are sorted by varexp.indiv{g} in descending
% order. The smallest number of latents whose cumulative shared variance
% reaches shared_varexp_threshold are labeled keep = 1. Remaining latents are
% labeled remove = 0.
%
% Two masks are saved:
%   SVExpFilter.rawlogical{g}: group-specific keep mask.
%       Across latent j may be kept in one group and removed in another.
%
%   SVExpFilter.logical{g}: across latents are intersected across groups.
%       Across latent j is kept in every group only if every group keeps it in
%       rawlogical. Within latents stay group-specific.
%
% This matches the subspace_similarity_dlag.m filter convention:
%   1 = keep, 0 = remove.

clc; clear;

% -------------------------------------------------------------------------
% User parameters
% -------------------------------------------------------------------------
data_content = 'raw_count';
% options:
% raw_count, raw_fr, z_within_trial, z_within_condition,
% z_across_conditions, demean_count_within_trial, demean_fr_within_trial,
% demean_pooledsd_within_condition

data_condition = [1:16];   % [] for pooled all-condition mode, or e.g. 1:16
runIdx = 1;

% Keep the smallest number of sorted latents that explain this fraction of
% group shared variance. 0.95 means keep latents up to 95% cumulative shared
% variance explained.
shared_varexp_threshold = 0.95;

save_results = true;


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

% -------------------------------------------------------------------------
% Main loop
% -------------------------------------------------------------------------
for cond_i = 1:numConditions
    if use_condition_mode
        this_condition = condition_list(cond_i);
        baseDir = ['./FA_Dlag_', data_content, '_condition', num2str(this_condition)];
        rowLabel = sprintf('condition%d', this_condition);
    else
        this_condition = [];
        baseDir = ['./FA_Dlag_', data_content];
        rowLabel = 'all';
    end

    tempfname = sprintf('%s/mat_results/run%03d', baseDir, runIdx);

    fprintf('\n============================================================\n');
    if isempty(this_condition)
        fprintf('Pick latents by shared variance explained: pooled all-condition mode\n');
    else
        fprintf('Pick latents by shared variance explained: condition %d\n', this_condition);
    end
    fprintf('Reading from %s\n', tempfname);

    bestFile = findOneFileLocal(tempfname, 'bestmodel*', true);
    Sbest = load(bestFile);

    requiredVars = {'bestModel', 'varexp'};
    for v = 1:numel(requiredVars)
        if ~isfield(Sbest, requiredVars{v})
            error('File %s is missing variable %s.', bestFile, requiredVars{v});
        end
    end

    bestModel = Sbest.bestModel;
    varexp = Sbest.varexp;

    if ~isfield(bestModel, 'xDim_across')
        error('bestModel is missing xDim_across.');
    end
    if ~isfield(bestModel, 'xDim_within')
        error('bestModel is missing xDim_within.');
    end

    xDim_across = bestModel.xDim_across;
    xDim_within = reshape(bestModel.xDim_within, 1, []);

    SVExpFilter = computeSVExpFilterLocal( ...
        varexp, xDim_across, xDim_within, shared_varexp_threshold);
    
    SVExpFilter.threshold=shared_varexp_threshold;


if save_results
    thresholdTag = thresholdTagLocal(shared_varexp_threshold);
    outFile = fullfile(tempfname, ['SVExpFilter_', thresholdTag, '.mat']);

    save(outFile, ...
        'SVExpFilter');

    fprintf('Saved %s\n', outFile);
end

end

%% ========================================================================
% Local functions
%% ========================================================================

function SVExpFilter = computeSVExpFilterLocal(varexp, xDim_across, xDim_within, threshold)

if ~isscalar(threshold) || ~isfinite(threshold) || threshold <= 0 || threshold > 1
    error('shared_varexp_threshold must be a finite scalar in (0, 1].');
end

if ~isstruct(varexp) || ~isfield(varexp, 'indiv')
    error('varexp must be a struct containing field indiv.');
end

if ~iscell(varexp.indiv)
    error('varexp.indiv must be a cell array with one entry per group.');
end

xDim_within = reshape(xDim_within, 1, []);
numGroups = numel(xDim_within);
localDims = xDim_across + xDim_within;

if numel(varexp.indiv) < numGroups
    error('varexp.indiv has %d groups, expected at least %d.', numel(varexp.indiv), numGroups);
end

SVExpFilter = struct();
SVExpFilter.indiv = cell(1, numGroups);
SVExpFilter.rawlogical = cell(1, numGroups);
SVExpFilter.logical = cell(1, numGroups);



for g = 1:numGroups
    v = reshape(varexp.indiv{g}, 1, []);
    if numel(v) ~= localDims(g)
        error('varexp.indiv{%d} has length %d, expected %d.', g, numel(v), localDims(g));
    end

    SVExpFilter.indiv{g} = v;

    vClean = v;
    vClean(~isfinite(vClean)) = 0;
    vClean(vClean < 0) = 0;

    [vSort, order] = sort(vClean, 'descend');
    totalVal = sum(vSort);

    keepMask = false(1, localDims(g));

    if totalVal > eps
        cumFrac = cumsum(vSort) ./ totalVal;
        M = find(cumFrac >= threshold, 1, 'first');
        if isempty(M)
            M = numel(vSort);
        end
        keepMask(order(1:M)) = true;

    end

    SVExpFilter.rawlogical{g} = double(keepMask);
    SVExpFilter.logical{g} = SVExpFilter.rawlogical{g};

end

% Across latents: intersect keep masks across neural groups.
% Within latents remain group-specific.
for j = 1:xDim_across
    keepAllGroups = true;
    for g = 1:numGroups
        keepAllGroups = keepAllGroups && (SVExpFilter.rawlogical{g}(j) ~= 0);
    end
    for g = 1:numGroups
        SVExpFilter.logical{g}(j) = double(keepAllGroups);
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



function tag = thresholdTagLocal(threshold)

% Convert numeric threshold to a filename-safe tag.
% Examples:
%   0.95  -> threshold0p95
%   0.975 -> threshold0p975
%   1     -> threshold1

tag = sprintf('threshold%.6g', threshold);
tag = strrep(tag, '.', 'p');
tag = strrep(tag, '-', 'm');

end
