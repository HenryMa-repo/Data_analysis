%% data_reconstruction.m
%
% Post-hoc DLAG data reconstruction from already saved bestmodel*.mat files.
%
% Run this script after plot_dlag_results.m has saved the best-model mat file.
%
% This script:
% 1) Does not refit the model.
% 2) Does not rerun inference.
% 3) Uses the existing seqEst.xsm in the saved bestmodel*.mat file.
% 4) Adds requested reconstruction fields to seqEst.
% 5) Can skip existing fields, so old reconstruction fields are not
%    overwritten unless overwrite_existing_recon_fields = true.
% 6) Overwrites the original bestmodel*.mat with the updated seqEst.
% 7) Saves R2 results separately as reconstruction_R2.mat.
%
% -------------------------------------------------------------------------
% Existing d / no-d reconstruction fields
% -------------------------------------------------------------------------
% seqEst(n).d
% seqEst(n).yRecon_use_across_no_d
% seqEst(n).yRecon_use_within_no_d
% seqEst(n).yRecon_use_all_no_d
%
% -------------------------------------------------------------------------
% Existing base reconstruction fields
% -------------------------------------------------------------------------
% seqEst(n).yRecon_use_across
% seqEst(n).yRecon_use_within
% seqEst(n).yRecon_use_all
% seqEst(n).yRecon_across_excl_within
% seqEst(n).yRecon_within_excl_across
%
% -------------------------------------------------------------------------
% Existing sampled-noise fields, if add_R_noise_reconstruction = true
% -------------------------------------------------------------------------
% seqEst(n).yRecon_use_across_with_R
% seqEst(n).yRecon_use_within_with_R
% seqEst(n).yRecon_use_all_with_R
% seqEst(n).yRecon_across_excl_within_with_R
% seqEst(n).yRecon_within_excl_across_with_R
%
% -------------------------------------------------------------------------
% Existing residual-preserving fields, if add_keep_resid_reconstruction = true
% -------------------------------------------------------------------------
% seqEst(n).yRecon_use_across_keep_resid
% seqEst(n).yRecon_use_within_keep_resid
% seqEst(n).yRecon_use_all_keep_resid
% seqEst(n).yRecon_across_excl_within_keep_resid
% seqEst(n).yRecon_within_excl_across_keep_resid
%
% -------------------------------------------------------------------------
% Existing directional reconstruction fields, if add_directional_reconstruction = true
% -------------------------------------------------------------------------
% seqEst(n).yRecon_use_feedback
% seqEst(n).yRecon_feedback_excl_within_ff_ambiguous
% seqEst(n).yRecon_feedback_excl_within
% seqEst(n).yRecon_feedback_excl_ff_ambiguous
%
% seqEst(n).yRecon_use_feedforward
% seqEst(n).yRecon_feedforward_excl_within_fb_ambiguous
% seqEst(n).yRecon_feedforward_excl_within
% seqEst(n).yRecon_feedforward_excl_fb_ambiguous
%
% -------------------------------------------------------------------------
% Directional + timescale reconstruction fields,
% if add_timescale_directional_reconstruction = true
% -------------------------------------------------------------------------
% Examples:
% seqEst(n).yRecon_use_ff_model_ts_v0_v120
% seqEst(n).yRecon_use_fb_model_ts_v0_v120
% seqEst(n).yRecon_use_ff_psd_ts_v0_v120
% seqEst(n).yRecon_use_fb_psd_ts_v0_v120
%
% These are use-style reconstructions:
%   d + selected FF or FB latent contribution
%
% Empty selections are saved as d-only numeric fields.
%
% The corresponding selected across-latent indices are saved in:
%   timescale_recon_info.(fieldName).selected_across_idx
%
% -------------------------------------------------------------------------
% Across / within + timescale reconstruction fields,
% if add_timescale_within_across_reconstruction = true
% -------------------------------------------------------------------------
% Examples:
% seqEst(n).yRecon_use_across_model_ts_v0_v120
% seqEst(n).yRecon_use_within_model_ts_v0_v120
% seqEst(n).yRecon_use_across_psd_ts_v0_v120
% seqEst(n).yRecon_use_within_psd_ts_v0_v120
%
% These are use-style reconstructions:
%   d + selected across or within latent contribution
%
% Empty selections are saved as d-only numeric fields.
%
% The corresponding selected indices are saved in:
%   timescale_recon_info.(fieldName).selected_across_idx
%   timescale_recon_info.(fieldName).selected_within_idx_by_group
%
% -------------------------------------------------------------------------
% R2 output
% -------------------------------------------------------------------------
% reconstruction_R2.mat contains recon_R2.
% For each supported reconstruction field that exists and is non-empty in
% every seqEst trial, recon_R2 has:
% .global_all
% .global_by_group
% .neuron_by_group

clc;
clear;

%% ------------------------------------------------------------------------
% User parameters
% -------------------------------------------------------------------------

data_content = 'raw_count';
% options usually include:
% raw_count, raw_fr, z_within_trial, z_within_condition,
% z_across_conditions, demean_count_within_trial, demean_fr_within_trial,
% demean_pooledsd_within_condition

data_condition = [];
% [] for pooled all-condition mode, or e.g. 1:16 for condition mode.

runIdx = 1;

%% ------------------------------------------------------------------------
% Reconstruction switches
% -------------------------------------------------------------------------

add_d_no_d_and_base_reconstruction = false;
add_R_noise_reconstruction = false;
add_keep_resid_reconstruction = false;

add_directional_reconstruction = false;

add_timescale_directional_reconstruction = false;
add_timescale_within_across_reconstruction = true;

% If false, existing seqEst fields will not be overwritten.
% Empty fields will still be filled.
overwrite_existing_recon_fields = false;

% Used only when add_R_noise_reconstruction = true.
use_fixed_noise_seed = false;
noise_seed = 1;

%% ------------------------------------------------------------------------
% Timescale reconstruction settings
% -------------------------------------------------------------------------

timescale_recon_sources = {'model-timescale'};
% options:
% {'model-timescale'}
% {'psd-timescale'}
% {'model-timescale', 'psd-timescale'}

timescale_ranges_ms = [
    0   120
    120 300
    300 Inf
];
% Left-closed and right-open:
% [0,120), [120,300), [300,Inf)

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

if add_R_noise_reconstruction
    if use_fixed_noise_seed
        rng(noise_seed, 'twister');
    else
        rng('shuffle');
    end
end

opts = struct();
opts.add_d_no_d_and_base_reconstruction = add_d_no_d_and_base_reconstruction;
opts.add_R_noise_reconstruction = add_R_noise_reconstruction;
opts.add_keep_resid_reconstruction = add_keep_resid_reconstruction;
opts.add_directional_reconstruction = add_directional_reconstruction;
opts.add_timescale_directional_reconstruction = add_timescale_directional_reconstruction;
opts.add_timescale_within_across_reconstruction = add_timescale_within_across_reconstruction;
opts.overwrite_existing_recon_fields = overwrite_existing_recon_fields;
opts.timescale_recon_specs = [];

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
        error(['seqEst.xsm not found in %s.\n', ...
               'This script expects saved inference results.'], bestFile);
    end

    if ~isfield(seqEst, 'y')
        error('seqEst.y not found in %s.', bestFile);
    end

    needDirectionalClass = ...
        opts.add_directional_reconstruction || ...
        opts.add_timescale_directional_reconstruction;

    needAnyTimescaleRecon = ...
        opts.add_timescale_directional_reconstruction || ...
        opts.add_timescale_within_across_reconstruction;

    needGpParams = ...
        needDirectionalClass || ...
        (needAnyTimescaleRecon && hasTimescaleSourceLocal(timescale_recon_sources, 'model-timescale'));

    latentClass = [];
    gp_params = [];

    if needGpParams
        gp_params = getGpParamsLocal(Sbest, bestFile);
    end

    if needDirectionalClass
        bootstrapFile = findOneFileLocal(tempfname, 'bootstrapResults*', true);
        fprintf('Loading bootstrap ambiguity file: %s\n', bootstrapFile);

        Sboot = load(bootstrapFile, 'ambiguousIdxs');

        if ~isfield(Sboot, 'ambiguousIdxs')
            error('%s is missing ambiguousIdxs.', bootstrapFile);
        end

        latentClass = classifyDlagLatentsLocal( ...
            params.xDim_across, gp_params, Sboot.ambiguousIdxs);

        printLatentClassificationLocal(latentClass);
    end

    if needAnyTimescaleRecon
        timescale_recon_specs = struct([]);

        if opts.add_timescale_directional_reconstruction
            timescale_directional_specs = buildTimescaleDirectionalSpecsLocal( ...
                params, ...
                latentClass, ...
                gp_params, ...
                tempfname, ...
                timescale_recon_sources, ...
                timescale_ranges_ms);

            timescale_recon_specs = appendStructArrayLocal( ...
                timescale_recon_specs, timescale_directional_specs);
        end

        if opts.add_timescale_within_across_reconstruction
            timescale_within_across_specs = buildTimescaleWithinAcrossSpecsLocal( ...
                params, ...
                gp_params, ...
                tempfname, ...
                timescale_recon_sources, ...
                timescale_ranges_ms);

            timescale_recon_specs = appendStructArrayLocal( ...
                timescale_recon_specs, timescale_within_across_specs);
        end

        opts.timescale_recon_specs = timescale_recon_specs;

        timescale_recon_info_new = buildTimescaleReconInfoLocal(timescale_recon_specs);

        if isfield(Sbest, 'timescale_recon_info') && isstruct(Sbest.timescale_recon_info)
            timescale_recon_info = mergeTimescaleReconInfoLocal( ...
                Sbest.timescale_recon_info, timescale_recon_info_new);
        else
            timescale_recon_info = timescale_recon_info_new;
        end

        fprintf('Timescale reconstruction fields to add:\n');

        for s = 1:numel(timescale_recon_specs)
            fprintf('  %s: %s\n', ...
                timescale_recon_specs(s).fieldName, ...
                formatTimescaleSpecSelectionLocal(timescale_recon_specs(s)));
        end
    else
        opts.timescale_recon_specs = [];

        if isfield(Sbest, 'timescale_recon_info') && isstruct(Sbest.timescale_recon_info)
            timescale_recon_info = Sbest.timescale_recon_info;
        else
            timescale_recon_info = struct();
        end
    end

    fprintf('Adding requested reconstruction fields to seqEst...\n');

    seqEst = addDlagReconstructionFieldsLocal( ...
        seqEst, ...
        params, ...
        latentClass, ...
        opts);

    fprintf('Computing reconstruction R2 from all existing supported reconstruction fields...\n');

    recon_R2 = computeReconstructionR2Local(seqEst, params.yDims);

    Sbest.seqEst = seqEst;

    if needAnyTimescaleRecon || ~isempty(fieldnames(timescale_recon_info))
        Sbest.timescale_recon_info = timescale_recon_info;
    end

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

function seqEst = addDlagReconstructionFieldsLocal(seqEst, params, latentClass, opts)

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

needStandardInMemory = ...
    opts.add_d_no_d_and_base_reconstruction || ...
    opts.add_R_noise_reconstruction || ...
    opts.add_keep_resid_reconstruction;

if opts.add_directional_reconstruction && isempty(latentClass)
    error('add_directional_reconstruction is true, but latentClass is empty.');
end

if opts.add_timescale_directional_reconstruction && isempty(latentClass)
    error('add_timescale_directional_reconstruction is true, but latentClass is empty.');
end

blocks = precomputeReconstructionBlocksLocal( ...
    params, latentClass, opts.add_directional_reconstruction);

needTimescaleRecon = ...
    opts.add_timescale_directional_reconstruction || ...
    opts.add_timescale_within_across_reconstruction;

if needTimescaleRecon
    timescaleBlocks = precomputeTimescaleBlocksLocal( ...
        params, opts.timescale_recon_specs);
else
    timescaleBlocks = [];
end

if opts.add_R_noise_reconstruction
    if ~isfield(params, 'R') || isempty(params.R)
        error('add_R_noise_reconstruction is true, but params.R is missing or empty.');
    end

    Rstd = buildDiagonalNoiseStdLocal(params.R, yDim);
else
    Rstd = [];
end

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

    if needStandardInMemory
        yRecon_d = yBase;

        yRecon_use_across_no_d = zeros(yDim, T);
        yRecon_use_within_no_d = zeros(yDim, T);
        yRecon_use_all_no_d = zeros(yDim, T);

        yRecon_use_across = yBase;
        yRecon_use_within = yBase;
        yRecon_use_all = yBase;

        yRecon_across_excl_within = yBase;
        yRecon_within_excl_across = yBase;
    end

    if opts.add_directional_reconstruction
        yRecon_use_feedback = yBase;
        yRecon_feedback_excl_within_ff_ambiguous = yBase;
        yRecon_feedback_excl_within = yBase;
        yRecon_feedback_excl_ff_ambiguous = yBase;

        yRecon_use_feedforward = yBase;
        yRecon_feedforward_excl_within_fb_ambiguous = yBase;
        yRecon_feedforward_excl_within = yBase;
        yRecon_feedforward_excl_fb_ambiguous = yBase;
    end

    if needTimescaleRecon
        nTsSpecs = numel(opts.timescale_recon_specs);
        yRecon_timescale = cell(1, nTsSpecs);

        for s = 1:nTsSpecs
            yRecon_timescale{s} = yBase;
        end
    end

    for groupIdx = 1:numGroups

        rows = blocks(groupIdx).obsIdx;
        d_g = repmat(params.d(rows), 1, T);

        if needStandardInMemory
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

            yRecon_use_across_no_d(rows, :) = Y_across;
            yRecon_use_within_no_d(rows, :) = Y_within;
            yRecon_use_all_no_d(rows, :) = Y_all;

            yRecon_use_across(rows, :) = d_g + Y_across;
            yRecon_use_within(rows, :) = d_g + Y_within;
            yRecon_use_all(rows, :) = d_g + Y_all;

            yRecon_across_excl_within(rows, :) = d_g + Y_across_excl_within;
            yRecon_within_excl_across(rows, :) = d_g + Y_within_excl_across;
        end

        if opts.add_directional_reconstruction
            X_feedback = seqEst(n).xsm(blocks(groupIdx).latIdx_feedback, :);
            X_feedforward = seqEst(n).xsm(blocks(groupIdx).latIdx_feedforward, :);

            Y_feedback = reconstructFromBlockLocal( ...
                blocks(groupIdx).Q_feedback, ...
                blocks(groupIdx).TT_feedback, ...
                X_feedback, ...
                numel(rows), T);

            Y_feedforward = reconstructFromBlockLocal( ...
                blocks(groupIdx).Q_feedforward, ...
                blocks(groupIdx).TT_feedforward, ...
                X_feedforward, ...
                numel(rows), T);

            Y_feedback_excl_within_ff_ambiguous = ...
                projectOrthogonalComplementLocal( ...
                Y_feedback, ...
                blocks(groupIdx).Q_remove_feedback_within_ff_ambiguous);

            Y_feedback_excl_within = ...
                projectOrthogonalComplementLocal( ...
                Y_feedback, ...
                blocks(groupIdx).Q_remove_feedback_within);

            Y_feedback_excl_ff_ambiguous = ...
                projectOrthogonalComplementLocal( ...
                Y_feedback, ...
                blocks(groupIdx).Q_remove_feedback_ff_ambiguous);

            Y_feedforward_excl_within_fb_ambiguous = ...
                projectOrthogonalComplementLocal( ...
                Y_feedforward, ...
                blocks(groupIdx).Q_remove_feedforward_within_fb_ambiguous);

            Y_feedforward_excl_within = ...
                projectOrthogonalComplementLocal( ...
                Y_feedforward, ...
                blocks(groupIdx).Q_remove_feedforward_within);

            Y_feedforward_excl_fb_ambiguous = ...
                projectOrthogonalComplementLocal( ...
                Y_feedforward, ...
                blocks(groupIdx).Q_remove_feedforward_fb_ambiguous);

            yRecon_use_feedback(rows, :) = d_g + Y_feedback;
            yRecon_feedback_excl_within_ff_ambiguous(rows, :) = ...
                d_g + Y_feedback_excl_within_ff_ambiguous;
            yRecon_feedback_excl_within(rows, :) = ...
                d_g + Y_feedback_excl_within;
            yRecon_feedback_excl_ff_ambiguous(rows, :) = ...
                d_g + Y_feedback_excl_ff_ambiguous;

            yRecon_use_feedforward(rows, :) = d_g + Y_feedforward;
            yRecon_feedforward_excl_within_fb_ambiguous(rows, :) = ...
                d_g + Y_feedforward_excl_within_fb_ambiguous;
            yRecon_feedforward_excl_within(rows, :) = ...
                d_g + Y_feedforward_excl_within;
            yRecon_feedforward_excl_fb_ambiguous(rows, :) = ...
                d_g + Y_feedforward_excl_fb_ambiguous;
        end

        if needTimescaleRecon
            for s = 1:numel(opts.timescale_recon_specs)
                X_selected = seqEst(n).xsm(timescaleBlocks(s, groupIdx).latIdx_selected, :);

                Y_selected = reconstructFromBlockLocal( ...
                    timescaleBlocks(s, groupIdx).Q_selected, ...
                    timescaleBlocks(s, groupIdx).TT_selected, ...
                    X_selected, ...
                    numel(rows), T);

                yRecon_timescale{s}(rows, :) = d_g + Y_selected;
            end
        end
    end

    if opts.add_d_no_d_and_base_reconstruction
        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'd', ...
            yRecon_d, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_use_across_no_d', ...
            yRecon_use_across_no_d, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_use_within_no_d', ...
            yRecon_use_within_no_d, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_use_all_no_d', ...
            yRecon_use_all_no_d, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_use_across', ...
            yRecon_use_across, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_use_within', ...
            yRecon_use_within, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_use_all', ...
            yRecon_use_all, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_across_excl_within', ...
            yRecon_across_excl_within, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_within_excl_across', ...
            yRecon_within_excl_across, opts.overwrite_existing_recon_fields);
    end

    if opts.add_R_noise_reconstruction
        noise_R = repmat(Rstd, 1, T) .* randn(yDim, T);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_use_across_with_R', ...
            yRecon_use_across + noise_R, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_use_within_with_R', ...
            yRecon_use_within + noise_R, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_use_all_with_R', ...
            yRecon_use_all + noise_R, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_across_excl_within_with_R', ...
            yRecon_across_excl_within + noise_R, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_within_excl_across_with_R', ...
            yRecon_within_excl_across + noise_R, opts.overwrite_existing_recon_fields);
    end

    if opts.add_keep_resid_reconstruction
        full_resid = seqEst(n).y - yRecon_use_all;

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_use_across_keep_resid', ...
            yRecon_use_across + full_resid, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_use_within_keep_resid', ...
            yRecon_use_within + full_resid, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_use_all_keep_resid', ...
            yRecon_use_all + full_resid, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_across_excl_within_keep_resid', ...
            yRecon_across_excl_within + full_resid, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_within_excl_across_keep_resid', ...
            yRecon_within_excl_across + full_resid, opts.overwrite_existing_recon_fields);
    end

    if opts.add_directional_reconstruction
        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_use_feedback', ...
            yRecon_use_feedback, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_feedback_excl_within_ff_ambiguous', ...
            yRecon_feedback_excl_within_ff_ambiguous, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_feedback_excl_within', ...
            yRecon_feedback_excl_within, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_feedback_excl_ff_ambiguous', ...
            yRecon_feedback_excl_ff_ambiguous, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_use_feedforward', ...
            yRecon_use_feedforward, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_feedforward_excl_within_fb_ambiguous', ...
            yRecon_feedforward_excl_within_fb_ambiguous, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_feedforward_excl_within', ...
            yRecon_feedforward_excl_within, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal(seqEst, n, 'yRecon_feedforward_excl_fb_ambiguous', ...
            yRecon_feedforward_excl_fb_ambiguous, opts.overwrite_existing_recon_fields);
    end

    if needTimescaleRecon
        for s = 1:numel(opts.timescale_recon_specs)
            fieldName = opts.timescale_recon_specs(s).fieldName;

            seqEst = maybeSetTrialFieldLocal(seqEst, n, fieldName, ...
                yRecon_timescale{s}, opts.overwrite_existing_recon_fields);
        end
    end
end
end

function seqEst = maybeSetTrialFieldLocal(seqEst, trialIdx, fieldName, value, overwriteExisting)

if overwriteExisting || ~isfield(seqEst, fieldName) || isempty(seqEst(trialIdx).(fieldName))
    seqEst(trialIdx).(fieldName) = value;
end
end

function blocks = precomputeReconstructionBlocksLocal(params, latentClass, addDirectional)

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

if addDirectional
    ffIdx = latentClass.feedforwardIdx(:)';
    fbIdx = latentClass.feedbackIdx(:)';
    ambiguousIdx = latentClass.ambiguousIdx(:)';

    if any(ffIdx < 1 | ffIdx > xDim_across)
        error('Feedforward indices are outside 1:xDim_across.');
    end

    if any(fbIdx < 1 | fbIdx > xDim_across)
        error('Feedback indices are outside 1:xDim_across.');
    end

    if any(ambiguousIdx < 1 | ambiguousIdx > xDim_across)
        error('Ambiguous indices are outside 1:xDim_across.');
    end
else
    ffIdx = [];
    fbIdx = [];
    ambiguousIdx = [];
end

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

    if addDirectional
        C_feedforward = Cg(:, ffIdx);
        C_feedback = Cg(:, fbIdx);
        C_ambiguous = Cg(:, ambiguousIdx);

        latIdx_feedforward = latIdx(ffIdx);
        latIdx_feedback = latIdx(fbIdx);
        latIdx_ambiguous = latIdx(ambiguousIdx);

        [Q_feedforward, TT_feedforward] = orthogonalizeLoadingBlockLocal(C_feedforward);
        [Q_feedback, TT_feedback] = orthogonalizeLoadingBlockLocal(C_feedback);

        C_remove_feedback_within_ff_ambiguous = [C_within, C_feedforward, C_ambiguous];
        C_remove_feedback_within = C_within;
        C_remove_feedback_ff_ambiguous = [C_feedforward, C_ambiguous];

        Q_remove_feedback_within_ff_ambiguous = ...
            orthonormalColumnSpaceLocal(C_remove_feedback_within_ff_ambiguous);
        Q_remove_feedback_within = ...
            orthonormalColumnSpaceLocal(C_remove_feedback_within);
        Q_remove_feedback_ff_ambiguous = ...
            orthonormalColumnSpaceLocal(C_remove_feedback_ff_ambiguous);

        C_remove_feedforward_within_fb_ambiguous = [C_within, C_feedback, C_ambiguous];
        C_remove_feedforward_within = C_within;
        C_remove_feedforward_fb_ambiguous = [C_feedback, C_ambiguous];

        Q_remove_feedforward_within_fb_ambiguous = ...
            orthonormalColumnSpaceLocal(C_remove_feedforward_within_fb_ambiguous);
        Q_remove_feedforward_within = ...
            orthonormalColumnSpaceLocal(C_remove_feedforward_within);
        Q_remove_feedforward_fb_ambiguous = ...
            orthonormalColumnSpaceLocal(C_remove_feedforward_fb_ambiguous);

        blocks(groupIdx).latIdx_feedforward = latIdx_feedforward;
        blocks(groupIdx).latIdx_feedback = latIdx_feedback;
        blocks(groupIdx).latIdx_ambiguous = latIdx_ambiguous;

        blocks(groupIdx).Q_feedforward = Q_feedforward;
        blocks(groupIdx).TT_feedforward = TT_feedforward;
        blocks(groupIdx).Q_feedback = Q_feedback;
        blocks(groupIdx).TT_feedback = TT_feedback;

        blocks(groupIdx).Q_remove_feedback_within_ff_ambiguous = ...
            Q_remove_feedback_within_ff_ambiguous;
        blocks(groupIdx).Q_remove_feedback_within = ...
            Q_remove_feedback_within;
        blocks(groupIdx).Q_remove_feedback_ff_ambiguous = ...
            Q_remove_feedback_ff_ambiguous;

        blocks(groupIdx).Q_remove_feedforward_within_fb_ambiguous = ...
            Q_remove_feedforward_within_fb_ambiguous;
        blocks(groupIdx).Q_remove_feedforward_within = ...
            Q_remove_feedforward_within;
        blocks(groupIdx).Q_remove_feedforward_fb_ambiguous = ...
            Q_remove_feedforward_fb_ambiguous;

        fprintf(['  Group %d: xAcross=%d, xWithin=%d, ', ...
            'FF=%d, FB=%d, Ambiguous=%d, ', ...
            'basisFF=%d, basisFB=%d, basisWithin=%d\n'], ...
            groupIdx, xDim_across, xDim_within(groupIdx), ...
            numel(ffIdx), numel(fbIdx), numel(ambiguousIdx), ...
            size(Q_feedforward, 2), size(Q_feedback, 2), size(Q_within, 2));
    else
        fprintf('  Group %d: xAcross=%d, xWithin=%d, basisAcross=%d, basisWithin=%d, basisAll=%d\n', ...
            groupIdx, xDim_across, xDim_within(groupIdx), ...
            size(Q_across, 2), size(Q_within, 2), size(Q_all, 2));
    end
end
end

function timescaleBlocks = precomputeTimescaleBlocksLocal(params, timescaleSpecs)

yDims = params.yDims;
xDim_across = params.xDim_across;
xDim_within = params.xDim_within;

numGroups = numel(yDims);
localDims = xDim_across + xDim_within;

obsStart = cumsum([1, yDims(1:end-1)]);
obsEnd = cumsum(yDims);

latStart = cumsum([1, localDims(1:end-1)]);
latEnd = cumsum(localDims);

nSpecs = numel(timescaleSpecs);

if nSpecs == 0
    timescaleBlocks = struct([]);
    return;
end

timescaleBlocks = repmat(struct( ...
    'obsIdx', [], ...
    'latIdx_selected', [], ...
    'Q_selected', [], ...
    'TT_selected', []), nSpecs, numGroups);

for s = 1:nSpecs
    kind = char(timescaleSpecs(s).recon_kind);

    for groupIdx = 1:numGroups
        obsIdx = obsStart(groupIdx):obsEnd(groupIdx);
        latIdx = latStart(groupIdx):latEnd(groupIdx);

        switch kind
            case {'ff', 'fb', 'across'}
                selectedAcrossIdx = timescaleSpecs(s).selected_across_idx(:)';

                if any(selectedAcrossIdx < 1 | selectedAcrossIdx > xDim_across)
                    error('Selected across indices are outside 1:xDim_across for field %s.', ...
                        timescaleSpecs(s).fieldName);
                end

                selectedLocalIdx = selectedAcrossIdx;

            case 'within'
                if ~isfield(timescaleSpecs(s), 'selected_within_idx_by_group') || ...
                        numel(timescaleSpecs(s).selected_within_idx_by_group) < groupIdx
                    error('Missing selected_within_idx_by_group for field %s.', ...
                        timescaleSpecs(s).fieldName);
                end

                selectedWithinIdx = timescaleSpecs(s).selected_within_idx_by_group{groupIdx};
                selectedWithinIdx = selectedWithinIdx(:)';

                if any(selectedWithinIdx < 1 | selectedWithinIdx > xDim_within(groupIdx))
                    error('Selected within indices are outside 1:xDim_within(groupIdx) for field %s, group %d.', ...
                        timescaleSpecs(s).fieldName, groupIdx);
                end

                selectedLocalIdx = xDim_across + selectedWithinIdx;

            otherwise
                error('Unknown timescale reconstruction kind: %s', kind);
        end

        latIdx_selected = latIdx(selectedLocalIdx);
        C_selected = params.C(obsIdx, latIdx_selected);

        [Q_selected, TT_selected] = orthogonalizeLoadingBlockLocal(C_selected);

        timescaleBlocks(s, groupIdx).obsIdx = obsIdx;
        timescaleBlocks(s, groupIdx).latIdx_selected = latIdx_selected;
        timescaleBlocks(s, groupIdx).Q_selected = Q_selected;
        timescaleBlocks(s, groupIdx).TT_selected = TT_selected;
    end
end
end

function specs = buildTimescaleDirectionalSpecsLocal( ...
    params, latentClass, gp_params, runDir, timescaleSources, timescaleRangesMs)

if isempty(latentClass)
    error('latentClass is required for timescale directional reconstruction.');
end

xDim_across = params.xDim_across;

if isempty(timescaleSources)
    specs = struct([]);
    return;
end

if ischar(timescaleSources) || isstring(timescaleSources)
    timescaleSources = {char(timescaleSources)};
end

if ~iscell(timescaleSources)
    error('timescale_recon_sources must be a character vector, string, or cell array.');
end

if isempty(timescaleRangesMs) || size(timescaleRangesMs, 2) ~= 2
    error('timescale_ranges_ms must be an N x 2 matrix.');
end

if any(~isfinite(timescaleRangesMs(:, 1)))
    error('Lower bounds in timescale_ranges_ms must be finite.');
end

if any(timescaleRangesMs(:, 2) <= timescaleRangesMs(:, 1))
    error('Each timescale range must satisfy upper > lower.');
end

specs = makeEmptyTimescaleSpecStructLocal();

for src_i = 1:numel(timescaleSources)
    sourceName = char(timescaleSources{src_i});
    sourceTag = makeTimescaleSourceTagLocal(sourceName);

    ts_across_ms = loadAcrossTimescaleLocal( ...
        sourceName, gp_params, runDir, xDim_across);

    for r = 1:size(timescaleRangesMs, 1)
        rangeMs = timescaleRangesMs(r, :);

        lower = rangeMs(1);
        upper = rangeMs(2);

        inRange = ts_across_ms >= lower & ts_across_ms < upper;
        inRange = reshape(inRange, 1, []);

        rangeTag = makeTimescaleRangeTagLocal(rangeMs);

        ffSelected = intersect(latentClass.feedforwardIdx(:)', find(inRange), 'stable');
        fbSelected = intersect(latentClass.feedbackIdx(:)', find(inRange), 'stable');

        fieldFF = sprintf('yRecon_use_ff_%s_ts_%s', sourceTag, rangeTag);
        fieldFB = sprintf('yRecon_use_fb_%s_ts_%s', sourceTag, rangeTag);

        specs(end+1) = makeOneTimescaleSpecLocal( ...
            fieldFF, 'ff', sourceName, sourceTag, rangeMs, rangeTag, ...
            ts_across_ms, ffSelected, {}, {}); %#ok<AGROW>

        specs(end+1) = makeOneTimescaleSpecLocal( ...
            fieldFB, 'fb', sourceName, sourceTag, rangeMs, rangeTag, ...
            ts_across_ms, fbSelected, {}, {}); %#ok<AGROW>
    end
end
end

function specs = buildTimescaleWithinAcrossSpecsLocal( ...
    params, gp_params, runDir, timescaleSources, timescaleRangesMs)

xDim_across = params.xDim_across;
xDim_within = params.xDim_within;
numGroups = numel(xDim_within);

if isempty(timescaleSources)
    specs = struct([]);
    return;
end

if ischar(timescaleSources) || isstring(timescaleSources)
    timescaleSources = {char(timescaleSources)};
end

if ~iscell(timescaleSources)
    error('timescale_recon_sources must be a character vector, string, or cell array.');
end

if isempty(timescaleRangesMs) || size(timescaleRangesMs, 2) ~= 2
    error('timescale_ranges_ms must be an N x 2 matrix.');
end

if any(~isfinite(timescaleRangesMs(:, 1)))
    error('Lower bounds in timescale_ranges_ms must be finite.');
end

if any(timescaleRangesMs(:, 2) <= timescaleRangesMs(:, 1))
    error('Each timescale range must satisfy upper > lower.');
end

specs = makeEmptyTimescaleSpecStructLocal();

for src_i = 1:numel(timescaleSources)
    sourceName = char(timescaleSources{src_i});
    sourceTag = makeTimescaleSourceTagLocal(sourceName);

    [ts_across_ms, ts_within_ms_by_group] = loadAcrossWithinTimescaleLocal( ...
        sourceName, gp_params, runDir, xDim_across, xDim_within);

    for r = 1:size(timescaleRangesMs, 1)
        rangeMs = timescaleRangesMs(r, :);

        lower = rangeMs(1);
        upper = rangeMs(2);

        rangeTag = makeTimescaleRangeTagLocal(rangeMs);

        acrossInRange = ts_across_ms >= lower & ts_across_ms < upper;
        acrossSelected = find(reshape(acrossInRange, 1, []));

        selectedWithinByGroup = cell(1, numGroups);

        for groupIdx = 1:numGroups
            ts_within_g = ts_within_ms_by_group{groupIdx};
            withinInRange = ts_within_g >= lower & ts_within_g < upper;
            selectedWithinByGroup{groupIdx} = find(reshape(withinInRange, 1, []));
        end

        fieldAcross = sprintf('yRecon_use_across_%s_ts_%s', sourceTag, rangeTag);
        fieldWithin = sprintf('yRecon_use_within_%s_ts_%s', sourceTag, rangeTag);

        specs(end+1) = makeOneTimescaleSpecLocal( ...
            fieldAcross, 'across', sourceName, sourceTag, rangeMs, rangeTag, ...
            ts_across_ms, acrossSelected, ts_within_ms_by_group, {}); %#ok<AGROW>

        specs(end+1) = makeOneTimescaleSpecLocal( ...
            fieldWithin, 'within', sourceName, sourceTag, rangeMs, rangeTag, ...
            ts_across_ms, [], ts_within_ms_by_group, selectedWithinByGroup); %#ok<AGROW>
    end
end
end

function specs = makeEmptyTimescaleSpecStructLocal()

specs = struct( ...
    'fieldName', {}, ...
    'r2Name', {}, ...
    'recon_kind', {}, ...
    'timescale_source', {}, ...
    'source_tag', {}, ...
    'range_ms', {}, ...
    'range_tag', {}, ...
    'timescale_across_ms', {}, ...
    'timescale_within_ms_by_group', {}, ...
    'selected_across_idx', {}, ...
    'selected_within_idx_by_group', {});
end

function spec = makeOneTimescaleSpecLocal( ...
    fieldName, reconKind, sourceName, sourceTag, rangeMs, rangeTag, ...
    tsAcrossMs, selectedAcrossIdx, tsWithinByGroup, selectedWithinByGroup)

if nargin < 9 || isempty(tsWithinByGroup)
    tsWithinByGroup = {};
end

if nargin < 10 || isempty(selectedWithinByGroup)
    selectedWithinByGroup = {};
end

spec = struct();
spec.fieldName = fieldName;
spec.r2Name = regexprep(fieldName, '^yRecon_', '');
spec.recon_kind = reconKind;
spec.timescale_source = sourceName;
spec.source_tag = sourceTag;
spec.range_ms = rangeMs;
spec.range_tag = rangeTag;
spec.timescale_across_ms = tsAcrossMs;
spec.timescale_within_ms_by_group = tsWithinByGroup;
spec.selected_across_idx = selectedAcrossIdx(:)';
spec.selected_within_idx_by_group = selectedWithinByGroup;
end

function sourceTag = makeTimescaleSourceTagLocal(sourceName)

switch char(sourceName)
    case 'model-timescale'
        sourceTag = 'model';

    case 'psd-timescale'
        sourceTag = 'psd';

    otherwise
        error('Unknown timescale source: %s', sourceName);
end
end

function ts_across_ms = loadAcrossTimescaleLocal(sourceName, gp_params, runDir, xDim_across)

switch char(sourceName)
    case 'model-timescale'
        if isempty(gp_params) || ~isfield(gp_params, 'tau_across')
            error('gp_params.tau_across is required for model-timescale reconstruction.');
        end

        ts_across_ms = reshape(double(gp_params.tau_across), 1, []);

    case 'psd-timescale'
        psdFile = fullfile(runDir, 'psd_timescale_stats.mat');

        if ~exist(psdFile, 'file')
            error(['psd_timescale_stats.mat not found:\n%s\n', ...
                   'Run compute_psd_timescale.m first, or remove psd-timescale from timescale_recon_sources.'], ...
                   psdFile);
        end

        Spsd = load(psdFile, 'PSD_timescale');

        if ~isfield(Spsd, 'PSD_timescale')
            error('PSD_timescale variable not found in:\n%s', psdFile);
        end

        PSD_timescale = Spsd.PSD_timescale;

        if isfield(PSD_timescale, 'across') && ...
                isfield(PSD_timescale.across, 'period_ms') && ...
                ~isempty(PSD_timescale.across.period_ms)

            ts_across_ms = reshape(double(PSD_timescale.across.period_ms), 1, []);

        elseif isfield(PSD_timescale, 'local') && ...
                numel(PSD_timescale.local) >= 1 && ...
                isfield(PSD_timescale.local(1), 'period_ms') && ...
                numel(PSD_timescale.local(1).period_ms) >= xDim_across

            ts_across_ms = reshape(double(PSD_timescale.local(1).period_ms(1:xDim_across)), 1, []);

        else
            error('Could not find PSD across period_ms in:\n%s', psdFile);
        end

    otherwise
        error('Unknown timescale source: %s', sourceName);
end

if numel(ts_across_ms) < xDim_across
    error('%s has fewer across timescales than xDim_across.', sourceName);
end

ts_across_ms = ts_across_ms(1:xDim_across);

if any(~isfinite(ts_across_ms))
    warning('%s contains non-finite across timescale values.', sourceName);
end
end

function [ts_across_ms, ts_within_ms_by_group] = loadAcrossWithinTimescaleLocal( ...
    sourceName, gp_params, runDir, xDim_across, xDim_within)

switch char(sourceName)
    case 'model-timescale'
        ts_across_ms = loadAcrossTimescaleLocal( ...
            sourceName, gp_params, runDir, xDim_across);

        ts_within_ms_by_group = loadModelWithinTimescaleLocal( ...
            gp_params, xDim_within);

    case 'psd-timescale'
        psdFile = fullfile(runDir, 'psd_timescale_stats.mat');

        if ~exist(psdFile, 'file')
            error(['psd_timescale_stats.mat not found:\n%s\n', ...
                   'Run compute_psd_timescale.m first, or remove psd-timescale from timescale_recon_sources.'], ...
                   psdFile);
        end

        Spsd = load(psdFile, 'PSD_timescale');

        if ~isfield(Spsd, 'PSD_timescale')
            error('PSD_timescale variable not found in:\n%s', psdFile);
        end

        PSD_timescale = Spsd.PSD_timescale;

        ts_across_ms = loadAcrossTimescaleLocal( ...
            sourceName, gp_params, runDir, xDim_across);

        ts_within_ms_by_group = loadPsdWithinTimescaleLocal( ...
            PSD_timescale, xDim_across, xDim_within, psdFile);

    otherwise
        error('Unknown timescale source: %s', sourceName);
end
end

function ts_within_ms_by_group = loadModelWithinTimescaleLocal(gp_params, xDim_within)

if isempty(gp_params) || ~isfield(gp_params, 'tau_within')
    error('gp_params.tau_within is required for model-timescale within reconstruction.');
end

tau_within = gp_params.tau_within;
xDim_within = reshape(xDim_within, 1, []);
numGroups = numel(xDim_within);
ts_within_ms_by_group = cell(1, numGroups);

if iscell(tau_within)
    if numel(tau_within) < numGroups
        error('gp_params.tau_within has fewer cells than number of groups.');
    end

    for groupIdx = 1:numGroups
        vals = reshape(double(tau_within{groupIdx}), 1, []);

        if numel(vals) < xDim_within(groupIdx)
            error('gp_params.tau_within{%d} has fewer values than xDim_within(%d).', ...
                groupIdx, groupIdx);
        end

        ts_within_ms_by_group{groupIdx} = vals(1:xDim_within(groupIdx));
    end

elseif isnumeric(tau_within)
    vals = double(tau_within);

    if isvector(vals)
        vals = reshape(vals, 1, []);

        if numel(vals) < sum(xDim_within)
            error('gp_params.tau_within vector has fewer values than sum(xDim_within).');
        end

        startIdx = 1;

        for groupIdx = 1:numGroups
            stopIdx = startIdx + xDim_within(groupIdx) - 1;
            ts_within_ms_by_group{groupIdx} = vals(startIdx:stopIdx);
            startIdx = stopIdx + 1;
        end

    elseif size(vals, 1) == numGroups
        for groupIdx = 1:numGroups
            rowVals = vals(groupIdx, :);

            if numel(rowVals) < xDim_within(groupIdx)
                error('gp_params.tau_within row %d has fewer values than xDim_within(%d).', ...
                    groupIdx, groupIdx);
            end

            ts_within_ms_by_group{groupIdx} = reshape(rowVals(1:xDim_within(groupIdx)), 1, []);
        end

    elseif size(vals, 2) == numGroups
        for groupIdx = 1:numGroups
            colVals = vals(:, groupIdx);

            if numel(colVals) < xDim_within(groupIdx)
                error('gp_params.tau_within column %d has fewer values than xDim_within(%d).', ...
                    groupIdx, groupIdx);
            end

            ts_within_ms_by_group{groupIdx} = reshape(colVals(1:xDim_within(groupIdx)), 1, []);
        end

    else
        error('Could not interpret numeric gp_params.tau_within.');
    end
else
    error('gp_params.tau_within must be a cell array or numeric array.');
end

for groupIdx = 1:numGroups
    if any(~isfinite(ts_within_ms_by_group{groupIdx}))
        warning('gp_params.tau_within for group %d contains non-finite values.', groupIdx);
    end
end
end

function ts_within_ms_by_group = loadPsdWithinTimescaleLocal( ...
    PSD_timescale, xDim_across, xDim_within, psdFile)

xDim_within = reshape(xDim_within, 1, []);
numGroups = numel(xDim_within);
ts_within_ms_by_group = cell(1, numGroups);

for groupIdx = 1:numGroups
    vals = [];

    if isfield(PSD_timescale, 'within') && ...
            numel(PSD_timescale.within) >= groupIdx && ...
            isfield(PSD_timescale.within(groupIdx), 'period_ms') && ...
            ~isempty(PSD_timescale.within(groupIdx).period_ms)

        vals = reshape(double(PSD_timescale.within(groupIdx).period_ms), 1, []);

    elseif isfield(PSD_timescale, 'local') && ...
            numel(PSD_timescale.local) >= groupIdx && ...
            isfield(PSD_timescale.local(groupIdx), 'period_ms') && ...
            numel(PSD_timescale.local(groupIdx).period_ms) >= xDim_across + xDim_within(groupIdx)

        localVals = reshape(double(PSD_timescale.local(groupIdx).period_ms), 1, []);
        vals = localVals((xDim_across + 1):(xDim_across + xDim_within(groupIdx)));
    end

    if isempty(vals)
        error('Could not find PSD within period_ms for group %d in:\n%s', ...
            groupIdx, psdFile);
    end

    if numel(vals) < xDim_within(groupIdx)
        error('PSD within period_ms for group %d has fewer values than xDim_within(%d).', ...
            groupIdx, groupIdx);
    end

    vals = vals(1:xDim_within(groupIdx));
    ts_within_ms_by_group{groupIdx} = vals;

    if any(~isfinite(vals))
        warning('PSD within period_ms for group %d contains non-finite values.', groupIdx);
    end
end
end

function rangeTag = makeTimescaleRangeTagLocal(rangeMs)

lowerTag = makeNumberTagLocal(rangeMs(1));
upperTag = makeNumberTagLocal(rangeMs(2));

rangeTag = [lowerTag, '_', upperTag];
end

function tag = makeNumberTagLocal(x)

if isinf(x)
    if x > 0
        tag = 'inf';
    else
        tag = 'neginf';
    end
    return;
end

tag = sprintf('%g', x);
tag = strrep(tag, '-', 'neg');
tag = strrep(tag, '.', 'p');

if isempty(tag)
    tag = 'x';
end

if ~isletter(tag(1))
    tag = ['v', tag];
end
end

function info = buildTimescaleReconInfoLocal(specs)

info = struct();

for s = 1:numel(specs)
    fieldName = specs(s).fieldName;

    if ~isempty(specs(s).selected_across_idx)
        info.(fieldName).selected_across_idx = specs(s).selected_across_idx;
    elseif any(strcmp(specs(s).recon_kind, {'ff', 'fb', 'across'}))
        info.(fieldName).selected_across_idx = [];
    end

    if ~isempty(specs(s).selected_within_idx_by_group)
        info.(fieldName).selected_within_idx_by_group = specs(s).selected_within_idx_by_group;
    end
end
end

function infoOut = mergeTimescaleReconInfoLocal(infoOld, infoNew)

infoOut = infoOld;

newFields = fieldnames(infoNew);

for i = 1:numel(newFields)
    f = newFields{i};
    infoOut.(f) = infoNew.(f);
end
end

function txt = formatTimescaleSpecSelectionLocal(spec)

kind = char(spec.recon_kind);

switch kind
    case {'ff', 'fb', 'across'}
        txt = sprintf('selected_across_idx = %s', ...
            mat2str(spec.selected_across_idx));

    case 'within'
        parts = cell(1, numel(spec.selected_within_idx_by_group));

        for groupIdx = 1:numel(spec.selected_within_idx_by_group)
            parts{groupIdx} = sprintf('G%d=%s', groupIdx, ...
                mat2str(spec.selected_within_idx_by_group{groupIdx}));
        end

        txt = ['selected_within_idx_by_group: ', strjoin(parts, ', ')];

    otherwise
        txt = 'unknown selection';
end
end

function Y = reconstructFromBlockLocal(Q, TT, X, yDim_group, T)

if isempty(Q) || isempty(TT) || isempty(X) || size(Q, 2) == 0
    Y = zeros(yDim_group, T);
    return;
end

Y = Q * (TT * X);
end

function Y_resid = projectOrthogonalComplementLocal(Y, Q_remove)

if isempty(Q_remove) || size(Q_remove, 2) == 0
    Y_resid = Y;
else
    Y_resid = Y - Q_remove * (Q_remove' * Y);
end
end

function [Q, TT] = orthogonalizeLoadingBlockLocal(L)

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
           'or at least not rank-shape-conflicting loading blocks. ', ...
           'Got yDim=%d, xDim=%d.'], yDim, xDim);
end

Q = U(:, 1:xDim);
TT = S(1:xDim, 1:xDim) * V(:, 1:xDim)';
end

function Q = orthonormalColumnSpaceLocal(L)

[yDim, nCols] = size(L);

if nCols == 0
    Q = zeros(yDim, 0);
    return;
end

if norm(L, 'fro') <= eps
    Q = zeros(yDim, 0);
    return;
end

[U, S, ~] = svd(L, 'econ');
s = diag(S);

if isempty(s) || max(s) <= eps
    Q = zeros(yDim, 0);
    return;
end

tol = max(size(L)) * eps(max(s));
r = sum(s > tol);

if r == 0
    Q = zeros(yDim, 0);
else
    Q = U(:, 1:r);
end
end

function Rstd = buildDiagonalNoiseStdLocal(R, yDim)

if ~isnumeric(R) || ~ismatrix(R)
    error('params.R must be a numeric matrix.');
end

if ~isequal(size(R), [yDim, yDim])
    error('params.R size is %d x %d, expected %d x %d.', ...
        size(R, 1), size(R, 2), yDim, yDim);
end

R = double(R);

if any(~isfinite(R(:)))
    error('params.R contains non-finite values.');
end

rvar = diag(R);
offdiag = R - diag(rvar);
tol = 1e-10 * max(1, max(abs(rvar)));

if max(abs(offdiag(:))) > tol
    error('params.R is expected to be diagonal, but off-diagonal entries are nonzero.');
end

if any(rvar < -tol)
    error('params.R has negative diagonal variance values.');
end

rvar(rvar < 0) = 0;
Rstd = sqrt(rvar);
Rstd = Rstd(:);
end

function recon_R2 = computeReconstructionR2Local(seqEst, yDims)

yDims = reshape(yDims, 1, []);
numGroups = numel(yDims);

r2_specs_all = getAllReconstructionR2SpecsLocal(seqEst);
r2_specs = selectExistingReconstructionSpecsLocal(seqEst, r2_specs_all);

if isempty(r2_specs)
    error('No usable reconstruction fields were found in seqEst.');
end

fprintf('  R2 will be computed for %d reconstruction fields:\n', size(r2_specs, 1));

for i = 1:size(r2_specs, 1)
    fprintf('    %s\n', r2_specs{i, 1});
end

Ytrue = [seqEst.y];

recon_R2 = struct();

for specIdx = 1:size(r2_specs, 1)
    r2Name = r2_specs{specIdx, 1};
    fieldName = r2_specs{specIdx, 2};

    Ypred = [seqEst.(fieldName)];

    if ~isequal(size(Ytrue), size(Ypred))
        error('Field %s has size mismatch with seqEst.y.', fieldName);
    end

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

function r2_specs = getAllReconstructionR2SpecsLocal(seqEst)

r2_specs = {
    'd_only', ...
    'd';

    'use_across_no_d', ...
    'yRecon_use_across_no_d';

    'use_within_no_d', ...
    'yRecon_use_within_no_d';

    'use_all_no_d', ...
    'yRecon_use_all_no_d';

    'use_across', ...
    'yRecon_use_across';

    'use_within', ...
    'yRecon_use_within';

    'use_all', ...
    'yRecon_use_all';

    'across_excl_within', ...
    'yRecon_across_excl_within';

    'within_excl_across', ...
    'yRecon_within_excl_across';

    'use_across_with_R', ...
    'yRecon_use_across_with_R';

    'use_within_with_R', ...
    'yRecon_use_within_with_R';

    'use_all_with_R', ...
    'yRecon_use_all_with_R';

    'across_excl_within_with_R', ...
    'yRecon_across_excl_within_with_R';

    'within_excl_across_with_R', ...
    'yRecon_within_excl_across_with_R';

    'use_across_keep_resid', ...
    'yRecon_use_across_keep_resid';

    'use_within_keep_resid', ...
    'yRecon_use_within_keep_resid';

    'use_all_keep_resid', ...
    'yRecon_use_all_keep_resid';

    'across_excl_within_keep_resid', ...
    'yRecon_across_excl_within_keep_resid';

    'within_excl_across_keep_resid', ...
    'yRecon_within_excl_across_keep_resid';

    'use_feedback', ...
    'yRecon_use_feedback';

    'feedback_excl_within_ff_ambiguous', ...
    'yRecon_feedback_excl_within_ff_ambiguous';

    'feedback_excl_within', ...
    'yRecon_feedback_excl_within';

    'feedback_excl_ff_ambiguous', ...
    'yRecon_feedback_excl_ff_ambiguous';

    'use_feedforward', ...
    'yRecon_use_feedforward';

    'feedforward_excl_within_fb_ambiguous', ...
    'yRecon_feedforward_excl_within_fb_ambiguous';

    'feedforward_excl_within', ...
    'yRecon_feedforward_excl_within';

    'feedforward_excl_fb_ambiguous', ...
    'yRecon_feedforward_excl_fb_ambiguous'
    };

if ~isempty(seqEst)
    allFields = fieldnames(seqEst);

    tsMask = startsWithLocal(allFields, 'yRecon_use_ff_') | ...
             startsWithLocal(allFields, 'yRecon_use_fb_') | ...
             startsWithLocal(allFields, 'yRecon_use_across_') | ...
             startsWithLocal(allFields, 'yRecon_use_within_');

    tsFields = allFields(tsMask);

    for i = 1:numel(tsFields)
        fieldName = tsFields{i};

        if contains(fieldName, '_ts_')
            r2Name = regexprep(fieldName, '^yRecon_', '');
            r2_specs(end+1, :) = {r2Name, fieldName}; %#ok<AGROW>
        end
    end
end

r2_specs = uniqueRowsStableLocal(r2_specs);
end

function mask = startsWithLocal(strs, prefix)

mask = false(size(strs));

for i = 1:numel(strs)
    s = strs{i};
    mask(i) = numel(s) >= numel(prefix) && strcmp(s(1:numel(prefix)), prefix);
end
end

function C = uniqueRowsStableLocal(C)

if isempty(C)
    return;
end

keep = true(size(C, 1), 1);

for i = 1:size(C, 1)
    for j = 1:(i-1)
        if strcmp(C{i, 1}, C{j, 1}) && strcmp(C{i, 2}, C{j, 2})
            keep(i) = false;
            break;
        end
    end
end

C = C(keep, :);
end

function r2_specs = selectExistingReconstructionSpecsLocal(seqEst, r2_specs_all)

keep = false(size(r2_specs_all, 1), 1);

for s = 1:size(r2_specs_all, 1)
    fieldName = r2_specs_all{s, 2};

    if ~isfield(seqEst, fieldName)
        continue;
    end

    allNonEmpty = true;

    for n = 1:numel(seqEst)
        if isempty(seqEst(n).(fieldName))
            allNonEmpty = false;
            break;
        end
    end

    if allNonEmpty
        keep(s) = true;
    else
        warning('Skipping R2 for %s because at least one trial is empty.', fieldName);
    end
end

r2_specs = r2_specs_all(keep, :);
end

function R2 = computeGlobalR2Local(Ytrue, Ypred)

if ~isequal(size(Ytrue), size(Ypred))
    error('Ytrue and Ypred must have the same size.');
end

valid = isfinite(Ytrue) & isfinite(Ypred);

numValid = sum(valid, 2);
Ytmp = Ytrue;
Ytmp(~valid) = 0;

mu = nan(size(Ytrue, 1), 1);
hasValid = numValid > 0;
mu(hasValid) = sum(Ytmp(hasValid, :), 2) ./ numValid(hasValid);

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

function latentClass = classifyDlagLatentsLocal(xDim_across, gp_params, ambiguousIdxs)

if ~isfield(gp_params, 'delays')
    error('gp_params must contain field delays.');
end

if ~isscalar(xDim_across) || xDim_across < 0 || mod(xDim_across, 1) ~= 0
    error('xDim_across must be a nonnegative integer scalar.');
end

acrossDelay = reshape(gp_params.delays, 1, []);

if numel(acrossDelay) < xDim_across
    error('gp_params.delays has fewer entries than xDim_across.');
end

acrossDelay = acrossDelay(1:xDim_across);

if isempty(ambiguousIdxs)
    ambiguousIdxs = [];
elseif islogical(ambiguousIdxs)
    ambiguousIdxs = find(ambiguousIdxs);
end

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

function printLatentClassificationLocal(latentClass)

fprintf('Directional across-latent classification:\n');
fprintf('  xDim_across: %d\n', numel(latentClass.acrossIdx));
fprintf('  Feedforward indices: %s\n', mat2str(latentClass.feedforwardIdx));
fprintf('  Feedback indices: %s\n', mat2str(latentClass.feedbackIdx));
fprintf('  Ambiguous indices: %s\n', mat2str(latentClass.ambiguousIdx));
fprintf('  Across delays: %s\n', mat2str(latentClass.acrossDelay));
end

function gp_params = getGpParamsLocal(Sbest, bestFile)

if isfield(Sbest, 'gp_params') && ~isempty(Sbest.gp_params)
    gp_params = Sbest.gp_params;
    return;
end

if isfield(Sbest, 'res') && isfield(Sbest.res, 'estParams') && ...
        isfield(Sbest.res.estParams, 'gp_params') && ...
        ~isempty(Sbest.res.estParams.gp_params)
    gp_params = Sbest.res.estParams.gp_params;
    return;
end

if isfield(Sbest, 'bestModel') && isfield(Sbest.bestModel, 'gp_params') && ...
        ~isempty(Sbest.bestModel.gp_params)
    gp_params = Sbest.bestModel.gp_params;
    return;
end

error('gp_params not found in %s.', bestFile);
end

function params = normalizeParamDimsLocal(params, bestModel)

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

startIdx = sum(yDims(1:groupIdx-1)) + 1;
rows = startIdx:(startIdx + yDims(groupIdx) - 1);
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

function tf = hasTimescaleSourceLocal(timescaleSources, targetSource)

if isempty(timescaleSources)
    tf = false;
    return;
end

if ischar(timescaleSources) || isstring(timescaleSources)
    timescaleSources = {char(timescaleSources)};
end

tf = false;

for i = 1:numel(timescaleSources)
    if strcmp(char(timescaleSources{i}), targetSource)
        tf = true;
        return;
    end
end
end

function out = appendStructArrayLocal(out, add)

if isempty(add)
    return;
end

if isempty(out)
    out = add;
else
    out = [out, add];
end
end