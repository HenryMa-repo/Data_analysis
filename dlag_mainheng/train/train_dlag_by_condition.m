%% =================
% 0a) Load demo data 
% ==================
clc;clear
% Synthetic data generated from a DLAG model
% dat_file = 'I:\np_data\RafiL001p0120_g1\catgt_RafiL001p0120_g1/model_data_allruns';
dat_file = 'I:\np_data\RafiL001p0120_g1\catgt_RafiL001p0120_g1\model_data_allruns';
fprintf('Reading from %s \n',dat_file);
load(dat_file);
stim_tag = '_2[Gpl2_2c_2sz_400_2_200isi]';
data_content = 'demean_count_within_trial';  
data_condtion=[2:16];
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

%% =========================
% 0b) Set up parallelization
% ==========================

% If parallelize is true, all cross-validation folds and bootstrap samples
% will be analyzed in parallel using Matlab's parfor construct. 
% If you have access to multiple cores, this provides significant speedup.
parallelize = true;
numWorkers = 16;      % Adjust this to your computer's specs  or try 32?

for i=1:numel(data_condtion)
% Load sequence data according to data_content
seqTrue = this_run.([data_content,'_by_condition'])(data_condtion(i)).trials;

%% =======================
% 1a) Fitting a DLAG model
% ========================

% Let's explicitly define all of the optional arguments, for 
% the sake of demonstration:

baseDir = ['./FA_Dlag_',data_content,'_condition',num2str(data_condtion(i))];        % Base directory where results will be saved
overwriteExisting = true; % Control whether existing results files are overwritten
saveData = false;         % Set to true to save train and test data (not recommended)
method = 'dlag';          % For now this is the only option, but that may change in the near future
% binWidth = 20;            % Sample period / spike count bin width, in units of time (e.g., ms)
% 
% 
% yDims = [10 10];          % Number of observed features (neurons) in each group (area)

covType = 'rbf';          % Type of GP covariance kernel ('rbf' or 'sg')
rGroups = [1 2];          % For performance evaluation, we can regress group 2's activity with group 1
startTau = 2*binWidth;    % Initial timescale, in the same units of time as binWidth
segLength = 20;           % Largest trial segment length, in no. of time points
init_method = 'static';   % Initialize DLAG with fitted pCCA parameters
learnDelays = true;       % Set to false if you want to fix delays at their initial value


tolLL = 1e-8;             % Log-likelihood convergence tolerance
freqLL = 1;               % Check for data log-likelihood convergence every freqLL EM iterations
freqParam = 10;           % Store intermediate delay and timescale estimates every freqParam EM iterations
minVarFrac = 0.001;       % Private noise variances will not be allowed to go below this value
verbose = true;           % Toggle printed progress updates
randomSeed = 0;           % Seed the random number generator, for reproducibility



%% ============================================================
% 2a) Cross-validate FA models to estimate total dimensionality 
%     (across+within) in each group.
%  ============================================================

% Change other input arguments as appropriate
      % Results will be saved in baseDir/mat_results/runXXX/,  
                          % where XXX is runIdx. Use a new runIdx for each dataset.
runIdx = 1;
numFolds = 4;
xDims = {0:yDims(1)-1, 0:yDims(2)-1}; % Sweep over these dimensionalities

fit_fa(runIdx, seqTrue, ...
       'baseDir', baseDir, ...
       'binWidth', binWidth, ...
       'numFolds', numFolds, ...
       'xDims', xDims, ...
       'yDims', yDims, ...
       'parallelize', parallelize, ...
       'randomSeed', randomSeed, ...
       'numWorkers', numWorkers, ...
       'overwriteExisting', overwriteExisting, ...
       'saveData', saveData);

%% Inspect full cross-validation results
[cvResults, bestModels] = getCrossValResults_fa(runIdx, 'baseDir', baseDir);

% Plot cross-validated performance vs estimated dimensionality
plotPerfvsDim_fa(cvResults, ...
                 'bestModels', bestModels);
               
% % % % % Collect the optimal total dimensionality for each group.
numGroups = length(yDims);
xDim_total_fa = nan(1,numGroups);
for groupIdx = 1:numGroups
    xDim_total_fa(groupIdx) = cvResults{groupIdx}(bestModels(groupIdx)).xDim;
end

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
%% ===============================================================
% 2b) Cross-validate DLAG models whose within- and across-group
%     dimensionalities are constrained to sum to the FA estimates.
%  ===============================================================

% Change other input arguments as appropriate
runIdx = 1;
numFolds = 4;
maxIters = 1000; % Limit EM iterations during cross-validation for speedup
fitAll = false; % Don't fit a model to all train data
% Determine DLAG models that satisfy the FA constraints
% xDims_grid = construct_xDimsGrid(xDim_total_fa);
%use opt d, 0.95% csve
xDims_grid = construct_xDimsGrid(xDim_opt_fa);

fit_dlag(runIdx, seqTrue, ...
         'baseDir', baseDir, ...
         'method', method, ...
         'binWidth', binWidth, ...
         'numFolds', numFolds, ...
         'fitAll', fitAll, ...
         'xDims_grid', xDims_grid, ...
         'yDims', yDims, ...
         'covType', covType, ...
         'rGroups', rGroups,...
         'startTau', startTau, ...
         'segLength', segLength, ...
         'init_method', init_method, ...
         'learnDelays', learnDelays, ...
         'maxIters', maxIters, ...
         'tolLL', tolLL, ...
         'freqLL', freqLL, ...
         'freqParam', freqParam, ...
         'minVarFrac', minVarFrac, ...
         'parallelize', parallelize, ...
         'randomSeed', randomSeed, ...
         'numWorkers', numWorkers, ...
         'overwriteExisting', overwriteExisting, ...
         'saveData', saveData);
     
%% Inspect cross-validation results
% Retrieve cross-validated results for all models in the results directory
[cvResults, ~] = getCrossValResults_dlag(runIdx, 'baseDir', baseDir);

% Plot a variety of performance metrics among the candidate models.
plotPerfvsDim_dlag(cvResults, 'xDims_grid', xDims_grid);

% Select the model with the optimal number among candidates
bestModel = getNumAcrossDim_dlag(cvResults, xDims_grid);

%% ====================================================================
% 2c) Fully train the optimal model selected in Section 2b, assuming EM 
%     iterations were limited during cross-validation.
% =====================================================================

% Change input arguments as appropriate
numFolds = 0;
xDims_across = bestModel.xDim_across;
xDims_within = num2cell(bestModel.xDim_within);
maxIters = 5e3;       % Set to even higher, if desired.

fit_dlag(runIdx, seqTrue, ...
         'baseDir', baseDir, ...
         'method', method, ...
         'binWidth', binWidth, ...
         'numFolds', numFolds, ...
         'xDims_across', xDims_across, ...
         'xDims_within', xDims_within, ...
         'yDims', yDims, ...
         'covType', covType, ...
         'rGroups', rGroups,...
         'startTau', startTau, ...
         'segLength', segLength, ...
         'init_method', init_method, ...
         'learnDelays', learnDelays, ...
         'maxIters', maxIters, ...
         'tolLL', tolLL, ...
         'freqLL', freqLL, ...
         'freqParam', freqParam, ...
         'minVarFrac', minVarFrac, ...
         'parallelize', false, ... % Only relevant for cross-validation
         'verbose', verbose, ...
         'randomSeed', randomSeed, ...
         'overwriteExisting', overwriteExisting, ...
         'saveData', saveData);                

end







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

