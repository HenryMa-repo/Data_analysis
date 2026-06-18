%% data_reconstruction.m
% Post-hoc DLAG data reconstruction from already saved bestmodel*.mat files.
%
% Run this script after plot_dlag_results.m has saved the best-model mat file.
%
% This script:
%   1) Does not refit the model.
%   2) Does not rerun inference.
%   3) Uses the existing seqEst.xsm in the saved bestmodel*.mat file.
%   4) Adds requested reconstruction fields to seqEst.
%   5) Can skip existing fields, so old reconstruction fields are not
%      overwritten unless overwrite_existing_recon_fields = true.
%   6) Overwrites the original bestmodel*.mat with the updated seqEst.
%   7) Saves R2 results separately as reconstruction_R2.mat.
%
% -------------------------------------------------------------------------
% Existing d / no-d reconstruction fields
% -------------------------------------------------------------------------
%   seqEst(n).d
%   seqEst(n).yRecon_use_across_no_d
%   seqEst(n).yRecon_use_within_no_d
%   seqEst(n).yRecon_use_all_no_d
%
% d:
%   Constant observation mean term, repeated across time:
%       repmat(params.d, 1, T)
%
% yRecon_use_*_no_d:
%   Latent-only reconstruction:
%       selected latent contribution only, without d
%
% -------------------------------------------------------------------------
% Existing base reconstruction fields
% -------------------------------------------------------------------------
%   seqEst(n).yRecon_use_across
%   seqEst(n).yRecon_use_within
%   seqEst(n).yRecon_use_all
%   seqEst(n).yRecon_across_excl_within
%   seqEst(n).yRecon_within_excl_across
%
% -------------------------------------------------------------------------
% Existing sampled-noise fields, if add_R_noise_reconstruction = true
% -------------------------------------------------------------------------
%   seqEst(n).yRecon_use_across_with_R
%   seqEst(n).yRecon_use_within_with_R
%   seqEst(n).yRecon_use_all_with_R
%   seqEst(n).yRecon_across_excl_within_with_R
%   seqEst(n).yRecon_within_excl_across_with_R
%
% -------------------------------------------------------------------------
% Existing residual-preserving fields, if add_keep_resid_reconstruction = true
% -------------------------------------------------------------------------
%   seqEst(n).yRecon_use_across_keep_resid
%   seqEst(n).yRecon_use_within_keep_resid
%   seqEst(n).yRecon_use_all_keep_resid
%   seqEst(n).yRecon_across_excl_within_keep_resid
%   seqEst(n).yRecon_within_excl_across_keep_resid
%
% -------------------------------------------------------------------------
% New directional reconstruction fields, if add_directional_reconstruction = true
% -------------------------------------------------------------------------
%   seqEst(n).yRecon_use_feedback
%   seqEst(n).yRecon_feedback_excl_within_ff_ambiguous
%   seqEst(n).yRecon_feedback_excl_within
%   seqEst(n).yRecon_feedback_excl_ff_ambiguous
%
%   seqEst(n).yRecon_use_feedforward
%   seqEst(n).yRecon_feedforward_excl_within_fb_ambiguous
%   seqEst(n).yRecon_feedforward_excl_within
%   seqEst(n).yRecon_feedforward_excl_fb_ambiguous
%
% Directional definitions:
%
%   use_feedback:
%       d + feedback signal
%
%   feedback_excl_within_ff_ambiguous:
%       d + feedback signal projected away from the combined remove subspace:
%       within + feedforward + ambiguous
%
%   feedback_excl_within:
%       d + feedback signal projected away from within only
%
%   feedback_excl_ff_ambiguous:
%       d + feedback signal projected away from the combined remove subspace:
%       feedforward + ambiguous
%
%   use_feedforward:
%       d + feedforward signal
%
%   feedforward_excl_within_fb_ambiguous:
%       d + feedforward signal projected away from the combined remove subspace:
%       within + feedback + ambiguous
%
%   feedforward_excl_within:
%       d + feedforward signal projected away from within only
%
%   feedforward_excl_fb_ambiguous:
%       d + feedforward signal projected away from the combined remove subspace:
%       feedback + ambiguous
%
% Important subspace rule:
%   For every exclusion case, the remove subspace is built by first
%   concatenating the relevant original loading blocks, then orthogonalizing
%   that combined loading block. Do not concatenate already-orthogonalized
%   Q matrices directly.
%
% Important directional note:
%   Directional fields are NOT expanded with _with_R or _keep_resid.
%
% R2 output:
%   reconstruction_R2.mat contains recon_R2.
%   For each supported reconstruction field that exists and is non-empty in
%   every seqEst trial, recon_R2 has:
%       .global_all
%       .global_by_group
%       .neuron_by_group
%
% Notes:
%   - Reconstruction uses internally orthogonalized loading blocks.
%   - SVD uses svd(...,'econ').
%   - For nonzero target reconstruction blocks, dimensions are taken as xDim,
%     not numerical rank, matching your previous reconstruction convention.
%   - For remove subspaces used only for projection, the actual column-space
%     rank is used. This is intentional because a remove subspace only needs
%     a projector basis.
%   - Completely zero loading blocks are treated as empty.
%   - Across/within/feedforward/feedback/ambiguous zero-dimensional cases
%     are allowed.
%   - R2 values are allowed to be negative and are not clipped.
%   - params.R is assumed to be a diagonal yDim x yDim matrix.

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

data_condition = [1:16];
% [] for pooled all-condition mode, or e.g. 1:16 for condition mode.

runIdx = 1;

%% ------------------------------------------------------------------------
% Reconstruction switches
% -------------------------------------------------------------------------
% Current recommended use:
%   If the bestmodel*.mat files already contain the previous d/no-d/base/
%   with_R/keep_resid fields, keep the defaults below. This will append only
%   the 8 new directional FF/FB fields and recompute reconstruction_R2.mat.
%
% For a fresh bestmodel without previous reconstruction fields, set:
%   add_d_no_d_and_base_reconstruction = true;
%   add_R_noise_reconstruction = true;       % optional
%   add_keep_resid_reconstruction = true;    % optional

add_d_no_d_and_base_reconstruction = false;
add_R_noise_reconstruction = false;
add_keep_resid_reconstruction = false;

add_directional_reconstruction = true;

% If false, existing seqEst fields will not be overwritten.
% Empty fields will still be filled.
overwrite_existing_recon_fields = false;

% Used only when add_R_noise_reconstruction = true.
% If use_fixed_noise_seed = false, MATLAB uses rng('shuffle'), so different
% script runs generate different random noise.
use_fixed_noise_seed = false;
noise_seed = 1;

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
opts.overwrite_existing_recon_fields = overwrite_existing_recon_fields;

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

    % ---------------------------------------------------------------------
    % Directional latent classification
    % ---------------------------------------------------------------------
    % Feedforward / feedback / ambiguous classification follows the same
    % logic as Latents_compare.m:
    %   delay > 0  -> feedforward
    %   delay < 0  -> feedback
    %   delay == 0, NaN delay, or bootstrap ambiguous -> ambiguous
    %
    % ambiguous latents are removed from feedforward and feedback.

    latentClass = [];

    if opts.add_directional_reconstruction
        gp_params = getGpParamsLocal(Sbest, bestFile);

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

    fprintf('Adding requested reconstruction fields to seqEst...\n');

    seqEst = addDlagReconstructionFieldsLocal( ...
        seqEst, ...
        params, ...
        latentClass, ...
        opts);

    fprintf('Computing reconstruction R2 from all existing supported reconstruction fields...\n');
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

function seqEst = addDlagReconstructionFieldsLocal(seqEst, params, latentClass, opts)
% Add requested reconstruction fields to each trial in seqEst.
%
% This function can append only new fields without overwriting old fields.
% This is useful when the bestmodel*.mat file already contains previous
% reconstruction fields and you only want to add new directional fields.

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

blocks = precomputeReconstructionBlocksLocal( ...
    params, latentClass, opts.add_directional_reconstruction);

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

    % ---------------------------------------------------------------------
    % Standard d/no-d/base reconstruction matrices.
    %
    % These may be needed even if not written, because _with_R and
    % _keep_resid depend on the standard base reconstructions.
    % ---------------------------------------------------------------------

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

    % ---------------------------------------------------------------------
    % Directional reconstruction matrices.
    % These are the 8 new fields. They are base/noiseless only.
    % ---------------------------------------------------------------------

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

    for groupIdx = 1:numGroups

        rows = blocks(groupIdx).obsIdx;
        d_g = repmat(params.d(rows), 1, T);

        % -----------------------------------------------------------------
        % Standard across / within / all reconstructions
        % -----------------------------------------------------------------

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

            % Latent-only, without d.
            yRecon_use_across_no_d(rows, :) = Y_across;
            yRecon_use_within_no_d(rows, :) = Y_within;
            yRecon_use_all_no_d(rows, :) = Y_all;

            % Original-style, with d.
            yRecon_use_across(rows, :) = d_g + Y_across;
            yRecon_use_within(rows, :) = d_g + Y_within;
            yRecon_use_all(rows, :) = d_g + Y_all;
            yRecon_across_excl_within(rows, :) = d_g + Y_across_excl_within;
            yRecon_within_excl_across(rows, :) = d_g + Y_within_excl_across;
        end

        % -----------------------------------------------------------------
        % Directional feedback / feedforward reconstructions
        % -----------------------------------------------------------------

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

            yRecon_use_feedback(rows, :) = ...
                d_g + Y_feedback;

            yRecon_feedback_excl_within_ff_ambiguous(rows, :) = ...
                d_g + Y_feedback_excl_within_ff_ambiguous;

            yRecon_feedback_excl_within(rows, :) = ...
                d_g + Y_feedback_excl_within;

            yRecon_feedback_excl_ff_ambiguous(rows, :) = ...
                d_g + Y_feedback_excl_ff_ambiguous;

            yRecon_use_feedforward(rows, :) = ...
                d_g + Y_feedforward;

            yRecon_feedforward_excl_within_fb_ambiguous(rows, :) = ...
                d_g + Y_feedforward_excl_within_fb_ambiguous;

            yRecon_feedforward_excl_within(rows, :) = ...
                d_g + Y_feedforward_excl_within;

            yRecon_feedforward_excl_fb_ambiguous(rows, :) = ...
                d_g + Y_feedforward_excl_fb_ambiguous;
        end
    end

    % ---------------------------------------------------------------------
    % Save d-only and base noiseless reconstructions, if requested.
    % ---------------------------------------------------------------------

    if opts.add_d_no_d_and_base_reconstruction
        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'd', ...
            yRecon_d, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_use_across_no_d', ...
            yRecon_use_across_no_d, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_use_within_no_d', ...
            yRecon_use_within_no_d, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_use_all_no_d', ...
            yRecon_use_all_no_d, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_use_across', ...
            yRecon_use_across, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_use_within', ...
            yRecon_use_within, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_use_all', ...
            yRecon_use_all, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_across_excl_within', ...
            yRecon_across_excl_within, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_within_excl_across', ...
            yRecon_within_excl_across, opts.overwrite_existing_recon_fields);
    end

    % ---------------------------------------------------------------------
    % Method A: sampled observation noise from params.R.
    %
    % The same sampled noise matrix is added to all five original-style
    % reconstructions within the same trial.
    %
    % Directional fields are intentionally NOT expanded with _with_R.
    % ---------------------------------------------------------------------

    if opts.add_R_noise_reconstruction
        noise_R = repmat(Rstd, 1, T) .* randn(yDim, T);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_use_across_with_R', ...
            yRecon_use_across + noise_R, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_use_within_with_R', ...
            yRecon_use_within + noise_R, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_use_all_with_R', ...
            yRecon_use_all + noise_R, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_across_excl_within_with_R', ...
            yRecon_across_excl_within + noise_R, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_within_excl_across_with_R', ...
            yRecon_within_excl_across + noise_R, opts.overwrite_existing_recon_fields);
    end

    % ---------------------------------------------------------------------
    % Method B: keep original full-model residual.
    %
    % residual = original data - full reconstruction.
    % Therefore:
    %   yRecon_use_all_keep_resid == seqEst(n).y
    %
    % Directional fields are intentionally NOT expanded with _keep_resid.
    % ---------------------------------------------------------------------

    if opts.add_keep_resid_reconstruction
        full_resid = seqEst(n).y - yRecon_use_all;

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_use_across_keep_resid', ...
            yRecon_use_across + full_resid, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_use_within_keep_resid', ...
            yRecon_use_within + full_resid, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_use_all_keep_resid', ...
            yRecon_use_all + full_resid, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_across_excl_within_keep_resid', ...
            yRecon_across_excl_within + full_resid, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_within_excl_across_keep_resid', ...
            yRecon_within_excl_across + full_resid, opts.overwrite_existing_recon_fields);
    end

    % ---------------------------------------------------------------------
    % New directional fields.
    % No _with_R and no _keep_resid are generated for these fields.
    % ---------------------------------------------------------------------

    if opts.add_directional_reconstruction

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_use_feedback', ...
            yRecon_use_feedback, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_feedback_excl_within_ff_ambiguous', ...
            yRecon_feedback_excl_within_ff_ambiguous, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_feedback_excl_within', ...
            yRecon_feedback_excl_within, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_feedback_excl_ff_ambiguous', ...
            yRecon_feedback_excl_ff_ambiguous, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_use_feedforward', ...
            yRecon_use_feedforward, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_feedforward_excl_within_fb_ambiguous', ...
            yRecon_feedforward_excl_within_fb_ambiguous, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_feedforward_excl_within', ...
            yRecon_feedforward_excl_within, opts.overwrite_existing_recon_fields);

        seqEst = maybeSetTrialFieldLocal( ...
            seqEst, n, 'yRecon_feedforward_excl_fb_ambiguous', ...
            yRecon_feedforward_excl_fb_ambiguous, opts.overwrite_existing_recon_fields);
    end
end
end

function seqEst = maybeSetTrialFieldLocal(seqEst, trialIdx, fieldName, value, overwriteExisting)
% Set seqEst(trialIdx).(fieldName) only when:
%   1) overwriteExisting is true, or
%   2) the field does not exist, or
%   3) the field exists but this trial's value is empty.
%
% This prevents repeated runs from overwriting previously generated fields.

if overwriteExisting || ~isfield(seqEst, fieldName) || isempty(seqEst(trialIdx).(fieldName))
    seqEst(trialIdx).(fieldName) = value;
end
end

function blocks = precomputeReconstructionBlocksLocal(params, latentClass, addDirectional)
% Precompute per-group loading bases and latent transforms.
%
% For standard reconstructions:
%   across, within, all are orthogonalized separately.
%
% For directional reconstructions:
%   feedback and feedforward target blocks are orthogonalized separately.
%
% For exclusion reconstructions:
%   each remove subspace is built by concatenating the original loading
%   blocks first, then orthogonalizing that combined block into a projector
%   basis.

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

        [Q_feedforward, TT_feedforward] = ...
            orthogonalizeLoadingBlockLocal(C_feedforward);

        [Q_feedback, TT_feedback] = ...
            orthogonalizeLoadingBlockLocal(C_feedback);

        % -------------------------------------------------------------
        % Feedback exclusion remove spaces.
        % For each case, concatenate original loading blocks first, then
        % orthogonalize the combined remove block into a projector basis.
        % -------------------------------------------------------------

        C_remove_feedback_within_ff_ambiguous = ...
            [C_within, C_feedforward, C_ambiguous];

        C_remove_feedback_within = ...
            C_within;

        C_remove_feedback_ff_ambiguous = ...
            [C_feedforward, C_ambiguous];

        Q_remove_feedback_within_ff_ambiguous = ...
            orthonormalColumnSpaceLocal(C_remove_feedback_within_ff_ambiguous);

        Q_remove_feedback_within = ...
            orthonormalColumnSpaceLocal(C_remove_feedback_within);

        Q_remove_feedback_ff_ambiguous = ...
            orthonormalColumnSpaceLocal(C_remove_feedback_ff_ambiguous);

        % -------------------------------------------------------------
        % Feedforward exclusion remove spaces.
        % -------------------------------------------------------------

        C_remove_feedforward_within_fb_ambiguous = ...
            [C_within, C_feedback, C_ambiguous];

        C_remove_feedforward_within = ...
            C_within;

        C_remove_feedforward_fb_ambiguous = ...
            [C_feedback, C_ambiguous];

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
% Since Q_remove is an orthonormal basis:
%   P = Q_remove * Q_remove'
%   Y_resid = (I - P) * Y
%           = Y - Q_remove * (Q_remove' * Y)

if isempty(Q_remove) || size(Q_remove, 2) == 0
    Y_resid = Y;
else
    Y_resid = Y - Q_remove * (Q_remove' * Y);
end
end

function [Q, TT] = orthogonalizeLoadingBlockLocal(L)
% Internal orthogonalization of one reconstruction loading block.
%
% For L with size yDim x xDim, return Q and TT such that:
%   L * X == Q * (TT * X)
% up to numerical precision.
%
% For nonzero reconstruction blocks, the basis dimension is xDim, matching
% your earlier reconstruction convention. If the whole block is exactly or
% effectively zero, the contribution is treated as empty.
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
        'or at least not rank-shape-conflicting loading blocks. ', ...
        'Got yDim=%d, xDim=%d.'], yDim, xDim);
end

Q = U(:, 1:xDim);
TT = S(1:xDim, 1:xDim) * V(:, 1:xDim)';
end

function Q = orthonormalColumnSpaceLocal(L)
% Return an orthonormal basis for the column span of L.
%
% This is used only for remove-subspace projection. Unlike reconstruction
% blocks, a remove subspace only needs the actual column span. Therefore
% rank-based basis selection is appropriate here.

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
% Extract observation-noise standard deviations from DLAG params.R.
%
% In standard DLAG, params.R is a yDim x yDim diagonal observation-noise
% covariance matrix.
%
% This function checks that format and returns:
%   Rstd = sqrt(diag(R))

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

% Clip tiny negative numerical values to zero.
rvar(rvar < 0) = 0;

Rstd = sqrt(rvar);
Rstd = Rstd(:);
end

function recon_R2 = computeReconstructionR2Local(seqEst, yDims)
% Compute global, group-wise global, and neuron-wise R2 for each supported
% reconstruction field currently present in seqEst.
%
% This function dynamically detects reconstruction fields. Therefore, if
% old fields already exist and this run only adds directional fields, R2 is
% recomputed for both old and new fields.

yDims = reshape(yDims, 1, []);
numGroups = numel(yDims);

r2_specs_all = getAllReconstructionR2SpecsLocal();
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

function r2_specs = getAllReconstructionR2SpecsLocal()
% All reconstruction fields currently supported by this script.

r2_specs = {
    % -------------------------------------------------------------
    % d / no-d diagnostic fields
    % -------------------------------------------------------------
    'd_only', ...
        'd';

    'use_across_no_d', ...
        'yRecon_use_across_no_d';

    'use_within_no_d', ...
        'yRecon_use_within_no_d';

    'use_all_no_d', ...
        'yRecon_use_all_no_d';

    % -------------------------------------------------------------
    % Original base fields
    % -------------------------------------------------------------
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

    % -------------------------------------------------------------
    % Sampled-noise fields
    % -------------------------------------------------------------
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

    % -------------------------------------------------------------
    % Residual-preserving fields
    % -------------------------------------------------------------
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

    % -------------------------------------------------------------
    % New feedback directional fields
    % -------------------------------------------------------------
    'use_feedback', ...
        'yRecon_use_feedback';

    'feedback_excl_within_ff_ambiguous', ...
        'yRecon_feedback_excl_within_ff_ambiguous';

    'feedback_excl_within', ...
        'yRecon_feedback_excl_within';

    'feedback_excl_ff_ambiguous', ...
        'yRecon_feedback_excl_ff_ambiguous';

    % -------------------------------------------------------------
    % New feedforward directional fields
    % -------------------------------------------------------------
    'use_feedforward', ...
        'yRecon_use_feedforward';

    'feedforward_excl_within_fb_ambiguous', ...
        'yRecon_feedforward_excl_within_fb_ambiguous';

    'feedforward_excl_within', ...
        'yRecon_feedforward_excl_within';

    'feedforward_excl_fb_ambiguous', ...
        'yRecon_feedforward_excl_fb_ambiguous'
};
end

function r2_specs = selectExistingReconstructionSpecsLocal(seqEst, r2_specs_all)
% Keep only reconstruction fields that exist and are non-empty for every
% trial. This allows the script to recompute R2 for whatever fields are
% currently present in seqEst.

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
% Global R2 with per-neuron mean baseline.
% Negative values are allowed.

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

function latentClass = classifyDlagLatentsLocal(xDim_across, gp_params, ambiguousIdxs)
% Classify across latents into feedforward, feedback, and ambiguous.
%
% Mirrors Latents_compare.m:
%   delay > 0           -> feedforward
%   delay < 0           -> feedback
%   bootstrap ambiguous -> ambiguous
%   delay == 0 or NaN   -> ambiguous
%
% Ambiguous latents are removed from feedforward and feedback.

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
% Print directional latent classification to command window.

fprintf('Directional across-latent classification:\n');
fprintf('  xDim_across: %d\n', numel(latentClass.acrossIdx));
fprintf('  Feedforward indices: %s\n', mat2str(latentClass.feedforwardIdx));
fprintf('  Feedback indices: %s\n', mat2str(latentClass.feedbackIdx));
fprintf('  Ambiguous indices: %s\n', mat2str(latentClass.ambiguousIdx));
fprintf('  Across delays: %s\n', mat2str(latentClass.acrossDelay));
end

function gp_params = getGpParamsLocal(Sbest, bestFile)
% Get gp_params from the loaded bestmodel file.
%
% The user's Latents_compare.m loads gp_params directly from bestmodel*.
% This helper allows a few fallback locations, but errors if none exist.

if isfield(Sbest, 'gp_params') && ~isempty(Sbest.gp_params)
    gp_params = Sbest.gp_params;
    return;
end

if isfield(Sbest, 'res') && isfield(Sbest.res, 'estParams') && ...
        isfield(Sbest.res.estParams, 'gp_params') && ~isempty(Sbest.res.estParams.gp_params)
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