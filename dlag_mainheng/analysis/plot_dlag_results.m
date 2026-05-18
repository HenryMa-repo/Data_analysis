% ==================
clc;clear
% Synthetic data generated from a DLAG model
% dat_file = 'I:\np_data\RafiL001p0120_g1\catgt_RafiL001p0120_g1/model_data_allruns';
dat_file = 'I:\np_data\RafiL001p0121_g1\catgt_RafiL001p0121_g1\model_data_allruns';
fprintf('Reading from %s \n',dat_file);
load(dat_file);
stim_tag = '_2[Gpl2_2c_2sz_400_2_200isi]';
data_content = 'raw_count';  
data_condtion=[];
% options:
% raw_count, raw_fr, z_within_trial, z_within_condition, 
% z_across_conditions, demean_count_within_trial, demean_fr_within_trial, demean_pooledsd_within_condition

% Extract all stim tags
all_run_tags = get_all_run_tags(model_data_allruns);


wanted_tag =stim_tag;


% Find the requested run by exact stim_tag match.
run_idx = find(strcmp(all_run_tags, wanted_tag));

if isempty(run_idx)
    error('Requested stim_tag not found %s', wanted_tag);
end

if numel(run_idx) > 1
    error('Duplicate stim_tag found %s', wanted_tag);
end


% -------------------------------------------------------------------------
% Get the selected run entry
% -------------------------------------------------------------------------
this_run = model_data_allruns{run_idx};

% Check whether requested data field exists
if ~isfield(this_run, data_content)
    error('Field %s not found in model_data_allruns{%d}.', data_content, run_idx);
end

% Bin width in ms
binWidth = this_run.bin_size * 1000;  % Sample period / spike count bin width, in units of time (e.g., ms)

% -------------------------------------------------------------------------
% Get yDims
% If nan_trial_strategy == 6, groupd is data-field-specific
% Otherwise use the shared groupd
% -------------------------------------------------------------------------
if isfield(this_run, 'nan_trial_strategy') && this_run.nan_trial_strategy == 6
    groupd_field = sprintf('%s_groupd', data_content);

    if ~isfield(this_run, groupd_field)
        error('For nan_trial_strategy == 6, field %s is missing.', groupd_field);
    end

    yDims = this_run.(groupd_field);  % Number of observed features (neurons) in each group (area)
else
    if ~isfield(this_run, 'groupd')
        error('Field groupd is missing in model_data_allruns{%d}.', run_idx);
    end

    yDims = this_run.groupd;  % Number of observed features (neurons) in each group (area)
end

xDims = {0:yDims(1)-1, 0:yDims(2)-1}; % Sweep over these dimensionalities

runIdx = 1;

if isempty(data_condtion)
    cond_list = [];  % no condition mode
else
    cond_list = data_condtion(:)'; % ensure row vector
end

% ================== LOOP START ==================
if isempty(cond_list)
    loop_range = 1;  % run once
else
    loop_range = cond_list;
end

for cond_i = 1:length(loop_range)

    if isempty(cond_list)
        baseDir = ['./FA_Dlag_',data_content];
        seqTrue = this_run.(data_content);
        cond_tag = '';
    else
        cond_val = loop_range(cond_i);
        baseDir = ['./FA_Dlag_',data_content,'_condition',num2str(cond_val)];
        seqTrue = this_run.([data_content,'_by_condition'])(cond_val).trials;
        cond_tag = ['_cond', num2str(cond_val)];
    end

rGroups = [1 2];          % For performance evaluation, we can regress group 2's activity with group 1

% If parallelize is true, all cross-validation folds and bootstrap samples
% will be analyzed in parallel using Matlab's parfor construct. 
% If you have access to multiple cores, this provides significant speedup.
parallelize = false;
numWorkers = 2;      % Adjust this to your computer's specs

tempfname = sprintf('%s/mat_results/run%03d', baseDir, runIdx);
%% Inspect full cross-validation results
[cvResults, bestModels] = getCrossValResults_fa(runIdx, 'baseDir', baseDir);

% Plot cross-validated performance vs estimated dimensionality
plotPerfvsDim_fa(cvResults, ...
                 'bestModels', bestModels);

savefig([tempfname,'/FA_cv_results.fig'])
exportgraphics(gcf, [tempfname,'/FA_cv_results.png'])


% Collect the optimal total dimensionality for each group.
numGroups = length(yDims);

%% here we are not use bestmodel, we use bestmodel with full data (cvResults.estParams), do LL' decomposition, to get lowest d can explain 0.95 shared var
%% code modified  decomposition and choose optd modified  CrossValFa and FactorAnalysisModelSelect(semedo paper)
csve=cell(length(yDims),1); %cumulative shared variance explained by the latent dimensions
 for groupIdx = 1:numGroups; 
% best factor number
if xDims{groupIdx}(bestModels(groupIdx)) == 0
    csve{groupIdx} = NaN;
else

    L = cvResults{groupIdx}(bestModels(groupIdx)).estParams.C;
    d = sort( eig( L*L' ), 'descend' );
    
    temp_cvLoss = (cumsum(d)/sum(d))';

    if xDims{groupIdx}(1) == 0
        temp_cvLoss = [0 temp_cvLoss];
        csve{groupIdx} = temp_cvLoss(xDims{groupIdx}+ 1);
    else
        csve{groupIdx} = temp_cvLoss(xDims{groupIdx});
    end
end

VAR_TRESHOLD = .95;  

if isnan(csve{groupIdx})
    xDim_opt_fa(groupIdx) = 0;
else
    xDim_opt_fa(groupIdx) = xDims{groupIdx}( find( csve{groupIdx} > VAR_TRESHOLD, 1 ) );   %d shared dimensionality 
end

 end

plotCSVEvsDim(csve, xDims, xDim_opt_fa);

savefig([tempfname,'/FA_cv_95dresults.fig'])
exportgraphics(gcf, [tempfname,'/FA_cv_95dresults.png'])


xDims_grid = construct_xDimsGrid(xDim_opt_fa);
%% Inspect cross-validation results
% Retrieve cross-validated results for all models in the results directory
[cvResults, ~] = getCrossValResults_dlag(runIdx, 'baseDir', baseDir);

% Plot a variety of performance metrics among the candidate models.
plotPerfvsDim_dlag(cvResults, 'xDims_grid', xDims_grid);

savefig([tempfname,'/Dlag_cv_results.fig'])
exportgraphics(gcf, [tempfname,'/Dlag_cv_results.png'])

% Select the model with the optimal number among candidates
bestModel = getNumAcrossDim_dlag(cvResults, xDims_grid);

%% =============================================================
% 3a) Evaluate how significantly each set of across-group delays
%     deviates from zero.
%  =============================================================

% Retrieve the best DLAG model
xDim_across = bestModel.xDim_across;
xDim_within = bestModel.xDim_within;
res = getModel_dlag(runIdx, xDim_across, xDim_within, 'baseDir', baseDir);

% Save all bootstrap results to a file
boot_fname = generate_inference_fname_dlag(runIdx, ...
                                           'bootstrapResults', ...
                                           xDim_across, ...
                                           xDim_within, ...
                                           'baseDir',baseDir);
numBootstrap = 100; % Number of bootstrap samples (the more the better)
delaySig = bootstrapDelaySignificance(seqTrue, ...
                                      res.estParams, ...
                                      numBootstrap, ...
                                      'parallelize', parallelize, ...
                                      'numWorkers', numWorkers);
% Label each delay as ambiguous (0) or unambiguous (1)
alpha = 0.05; % Significance level
ambiguousIdxs = find(delaySig >= alpha);
fprintf('Indexes of ambiguous delays: %s\n', num2str(ambiguousIdxs)); 


% Visualize non-zero and statistically ambiguous delays
gp_params=plotGPparams_dlag(res.estParams, binWidth, rGroups, ...
                  'plotAcross', true, ...
                  'plotWithin', false, ...
                  'units', 'ms', ...
                  'sig', delaySig, ...
                  'alpha', alpha);

if isempty(gp_params.DelayMatrix)
    gp_params.delays=[];
else
 gp_params.delays = gp_params.DelayMatrix(rGroups(2),:) ...
           - gp_params.DelayMatrix(rGroups(1),:);
end

savefig([tempfname,'/Delay_timescale_point_es_results.fig'])
exportgraphics(gcf, [tempfname,'/Delay_timescale_point_es_results.png'])

%% ================================================================
% 3b) Construct bootstrapped confidence intervals for latent delays
%     and timescales.
%  ================================================================

alpha = 0.05; % Construct (1-alpha) confidence intervals
bootParams = bootstrapGPparams(seqTrue, ...
                               res.estParams, ...
                               binWidth, ...
                               numBootstrap, ...
                               'alpha', alpha, ...
                               'parallelize', parallelize, ...
                               'numWorkers', numWorkers, ...
                               'segLength', Inf, ...
                               'tolLL', 1e-4, ...
                               'maxIters', 10);
save(boot_fname, 'delaySig','alpha',"ambiguousIdxs",'bootParams',"numBootstrap");
plotBootstrapGPparams_dlag(res.estParams, bootParams, binWidth, rGroups,...
                           'overlayParams', false);

savefig([tempfname,'/Delay_timescale_interval_es_results.fig'])
exportgraphics(gcf, [tempfname,'/Delay_timescale_interval_es_results.png'])

%%
% covType = 'rbf';          % Type of GP covariance kernel ('rbf' or 'sg')
% rGroups = [1 2];          % For performance evaluation, we can regress group 2's activity with group 1
% startTau = 2*binWidth;    % Initial timescale, in the same units of time as binWidth
% segLength = 25;           % Largest trial segment length, in no. of time points
% init_method = 'static';   % Initialize DLAG with fitted pCCA parameters
% learnDelays = true;       % Set to false if you want to fix delays at their initial value
% 
% 
% tolLL = 1e-8;             % Log-likelihood convergence tolerance
% freqLL = 1;               % Check for data log-likelihood convergence every freqLL EM iterations
% freqParam = 10;           % Store intermediate delay and timescale estimates every freqParam EM iterations
% minVarFrac = 0.001;       % Private noise variances will not be allowed to go below this value
% verbose = true;           % Toggle printed progress updates
% randomSeed = 0;           % Seed the random number generator, for reproducibility
% 
% % Plot training progress of various quantities. These plots can help with
% % troubleshooting, if necessary.
% plotFittingProgress(res, ...
%                     'freqLL', freqLL, ...
%                     'freqParam', freqParam, ...
%                     'units', 'ms');
% 


%% ========================================================
% Plot estimated latents
[seqEst, ~] = exactInferenceWithLL_dlag(seqTrue, res.estParams);
plotDimsVsTime_dlag(seqEst, 'xsm', res.estParams, res.binWidth, ...
                  'nPlotMax', 20, ...
                  'plotSingle', true, ...
                  'plotMean', false, ...
                  'units', []);

savefig([tempfname,'/Raw_latents_res.fig'])
exportgraphics(gcf, [tempfname,'/Raw_latents_res.png'])

% 1c) Visually scale latent trajectories by various metrics
% =========================================================

% Scale by variance explained
total = false; % true: denominator is total variance; else shared variance
[varexp, domexp] = computeVarExp_dlag(res.estParams, total);
[seqEst, sortParams] = scaleByVarExp(seqEst, ...
                                     res.estParams, ...
                                     varexp.indiv, ...
                                     'sortDims', true);
plotDimsVsTime_dlag(seqEst, 'xve', sortParams, res.binWidth, ...
                  'nPlotMax', 20, ...
                  'plotSingle', true, ...
                  'plotMean', false, ...
                  'units', []);

% Scale across-area latents by zero-delay correlation,  correlative mode
popcorr = computePopCorr_dlag(res.estParams);
[seqEst, sortParams] = scaleByCorr(seqEst, ...
                                   res.estParams, ...
                                   popcorr.lcorr, ...
                                   rGroups, ...
                                   'sortDims', true);
plotDimsVsTime_dlag(seqEst, 'xcorre', sortParams, res.binWidth, ...
                  'nPlotMax', 100, ...
                  'plotSingle', true, ...
                  'plotMean', false, ...
                  'units', []);
% correlative mode is not dominent mode.


% Scale across-area latents by zero-delay covariance
popcov = computePopCov_dlag(res.estParams);
[seqEst, sortParams] = scaleByCov(seqEst, ...
                                   res.estParams, ...
                                   popcov.lcov, ...
                                   rGroups, ...
                                   'sortDims', true);
plotDimsVsTime_dlag(seqEst, 'xcve', sortParams, res.binWidth, ...
                  'nPlotMax', 20, ...
                  'plotSingle', true, ...
                  'plotMean', false, ...
                  'units', []);



%% ================================================
% 1d) Project latents onto different types of modes
% =================================================

% Project latents onto dominant modes
seqEst = dominantProjection_dlag(seqEst, res.estParams, ...
                                 'includeAcross', true, ...
                                 'includeWithin', true);
plotDimsVsTime(seqEst, 'xdom', res.binWidth, ...
               'nPlotMax', 20, ...
               'nCol', xDim_across + max(xDim_within), ...
               'plotSingle', true, ...
               'plotMean', false, ...
               'units', 'ms');

savefig([tempfname,'/Projection_dom_modes_res.fig'])
exportgraphics(gcf, [tempfname,'/Projection_dom_modes_res.png'])

% Project latents onto zero-delay correlative modes
seqEst = correlativeProjection_dlag(seqEst, res.estParams,...
                                    'orth', true);
plotDimsVsTime(seqEst, 'xcorr', res.binWidth, ...
               'nPlotMax', 20, ...
               'nCol', xDim_across, ...
               'plotSingle', true, ...
               'plotMean', false, ...
               'units', 'ms'); 

savefig([tempfname,'/Projection_corr_modes_res.fig'])
exportgraphics(gcf, [tempfname,'/Projection_corr_modes_res.png'])

% Project latents onto zero-delay covariant modes
seqEst = covariantProjection_dlag(seqEst, res.estParams);
plotDimsVsTime(seqEst, 'xcov', res.binWidth, ...
               'nPlotMax', 20, ...
               'nCol', xDim_across, ...
               'plotSingle', true, ...
               'plotMean', false, ...
               'units', 'ms'); 

savefig([tempfname,'/Projection_cov_modes_res.fig'])
exportgraphics(gcf, [tempfname,'/Projection_cov_modes_res.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project latents onto zero-delay predictive modes
seqEst = predictiveProjection_dlag(seqEst,res.estParams, ...
                                   'orth', false, ...
                                   'groupIdxs', rGroups);
% Note that the order of rows corrsponds to the order of rGroups, where
% rGroups(1) is the source group, and rGroups(2) is the target group.
plotDimsVsTime(seqEst, 'xpred', res.binWidth, ...
               'nPlotMax', 20, ...
               'nCol', xDim_across, ...
               'plotSingle', true, ...
               'plotMean', false, ...
               'units', 'ms');

savefig([tempfname,'/Projection_pred_modes_res.fig'])
exportgraphics(gcf, [tempfname,'/Projection_pred_modes_res.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Related to dominant/covariant modes, inspect the spectra of DLAG loading 
% matrices. These give the amount of shared variance / cross-covariance 
% explained by each mode.
cutoffPC = 0.95;
d_shared = findSharedDimCutoff_dlag(res.estParams, cutoffPC, 'plotSpec', true)

savefig([tempfname,'/varexp_acrosscovexp_dom_covariant_mode.fig'])
exportgraphics(gcf, [tempfname,'/varexp_acrosscovexp_dom_covariant_mode.png'])
% Save best model details  to a file
bestm_fname = generate_inference_fname_dlag(runIdx, ...
                                           'bestmodel', ...
                                           xDim_across, ...
                                           xDim_within, ...
                                           'baseDir',baseDir);

save(bestm_fname,"bestModel","res","seqEst","varexp","domexp","popcorr","popcov","gp_params","cutoffPC","d_shared")
end % ================== LOOP END ==================


% % % % % % % % % % % % % % % % % % % % %% =============================================================================================
% % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % %% will be used in the future
% % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % Visualize the top three modes in 3D space
% % % % % % % % % % % % % % % % % % % % xspec = 'xdom';  % Set to any of the mode types above
% % % % % % % % % % % % % % % % % % % % plotTraj(seqEst, xspec, ...
% % % % % % % % % % % % % % % % % % % %          'dimsToPlot', 9:11, ...
% % % % % % % % % % % % % % % % % % % %          'nPlotMax', 20, ...
% % % % % % % % % % % % % % % % % % % %          'plotSingle', true, ...
% % % % % % % % % % % % % % % % % % % %          'plotMean', false);     
% % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % %% ==========================================
% % % % % % % % % % % % % % % % % % % % % 1e) Denoise observations using a DLAG model
% % % % % % % % % % % % % % % % % % % % % ===========================================
% % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % Denoise observations
% % % % % % % % % % % % % % % % % % % % [seqEst, ~, ~, ~, ~, ~] = denoise_dlag(seqEst, res.estParams);
% % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % Compare PSTHs of raw observations to PSTHs of denoised observations
% % % % % % % % % % % % % % % % % % % % psth_raw = get_psth(seqEst, 'spec', 'y');
% % % % % % % % % % % % % % % % % % % % spec = sprintf('yDenoisedOrth%02d', sum(xDim_across + xDim_within));
% % % % % % % % % % % % % % % % % % % % psth_denoised = get_psth(seqEst, 'spec', spec);
% % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % Raster plots
% % % % % % % % % % % % % % % % % % % % plotSeqRaster(psth_raw, res.binWidth, 'units', 'ms');
% % % % % % % % % % % % % % % % % % % % plotSeqRaster(psth_denoised, res.binWidth, 'units', 'ms');
% % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % Heat maps
% % % % % % % % % % % % % % % % % % % % figure;
% % % % % % % % % % % % % % % % % % % % hold on;
% % % % % % % % % % % % % % % % % % % % imagesc(flipud(psth_raw));
% % % % % % % % % % % % % % % % % % % % colormap('pink');
% % % % % % % % % % % % % % % % % % % % colorbar;
% % % % % % % % % % % % % % % % % % % % axis square;
% % % % % % % % % % % % % % % % % % % % xlabel('Time (ms)');
% % % % % % % % % % % % % % % % % % % % ylabel('Neurons');
% % % % % % % % % % % % % % % % % % % % title(sprintf('PSTHs, raw'));
% % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % figure;
% % % % % % % % % % % % % % % % % % % % hold on;
% % % % % % % % % % % % % % % % % % % % imagesc(flipud(psth_denoised));
% % % % % % % % % % % % % % % % % % % % colormap('pink');
% % % % % % % % % % % % % % % % % % % % colorbar;
% % % % % % % % % % % % % % % % % % % % axis square;
% % % % % % % % % % % % % % % % % % % % xlabel('Time (ms)');
% % % % % % % % % % % % % % % % % % % % ylabel('Neurons');
% % % % % % % % % % % % % % % % % % % % title(sprintf('PSTHs, denoised'));
% % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % 


















function all_tags = get_all_run_tags(model_data_allruns)
%% =========================================================================
% get_all_run_tags
%
% Purpose:
%   Extract all stim_tag values from model_data_allruns.
%
% Input:
%   unit_condition_metrics : cell array
%       One entry per run.
%
% Output:
%   all_tags : cell array of char
%       Run labels for all runs in this ksDir.
% =========================================================================

all_tags = cell(numel(model_data_allruns), 1);

for j = 1:numel(model_data_allruns)
    if ~isfield(model_data_allruns{j}, 'stim_tag')
        error('stim_tag missing in unit_condition_metrics{%d}.', j);
    end
    all_tags{j} = model_data_allruns{j}.stim_tag;
end

end