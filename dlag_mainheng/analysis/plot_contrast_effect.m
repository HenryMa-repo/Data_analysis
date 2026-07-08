%% =========================================================================
% plot_contrast_effect
%
% Purpose:
% For one selected stim_tag/run:
% 1. Load seqEst from DLAG bestmodel.
% 2. Use selected seqEst field:
%       response_data = seqEst.(analysis_field)
% 3. Pool trials by relative contrast level:
%       low contrast vs high contrast
%    Low/high are defined within each stim_name by sorting contrast values.
%    Then all low-contrast trials are pooled together, and all high-contrast
%    trials are pooled together.
%    Other condition features, such as size, direction, stim_name, etc.,
%    are not used for further subdivision.
% 4. For each trial, compute response within each analysis time window:
%       trial_response = sum(seqEst(t).(analysis_field)(:, bin_start:bin_end), 2)
% 5. Average trial_response across all low-contrast trials and all
%    high-contrast trials.
% 6. Compute contrast effect metrics for each cell:
%       delta_HL = H - L
%       CMI      = (H - L) ./ (abs(H) + abs(L))
% 7. Plot low vs high response separately for each group.
% 8. Save result mat and figures.
%
% Time-window mode:
% - binsize is read from this_run.bin_size, in seconds.
% - nbin is read from seqEst(1).T.
% - User parameters analysis_window and sliding_step are also in seconds.
% - If the window covers the full trial and only one window is generated,
%   filenames are exactly the old filenames, with no window tag.
% - If multiple windows are generated, filenames get:
%       _wn0_02_st0_02_w01
%       _wn0_02_st0_02_w02
%   No bin range is added to the filename.
%
% Model modes:
% - data_condition = []
%       Use pooled all-condition model:
%       ./FA_Dlag_<data_content>/mat_results/runXXX/bestmodel*
%       Save outputs into that model's runDir.
%
% - data_condition = 1:16, or any condition index list
%       Use condition-specific models:
%       ./FA_Dlag_<data_content>_conditionN/mat_results/runXXX/bestmodel*
%       Because contrast effect is across conditions, this mode first pools
%       responses from all requested condition-specific models by contrast,
%       then computes one contrast effect result.
%       Save outputs into this script's folder.
%
% Multiple analysis fields:
% - analysis_fields can contain one or multiple seqEst fields.
% - Each field is processed independently.
% - After each field is finished, close all figures.
%
% Loading strategy:
% - all_condition_model:
%       Load the pooled bestmodel only once, then read all analysis_fields
%       from the loaded seqEst.
%
% - condition_specific_models:
%       Load one condition model at a time to avoid keeping all condition
%       models in memory. For each loaded condition, read all analysis_fields
%       before clearing that condition model.
% =========================================================================

clc;
clear;

%% ----------------------- User parameters -----------------------

data_content = 'raw_count';
% options:
% raw_count, raw_fr, z_within_trial, z_within_condition,
% z_across_conditions, demean_count_within_trial,
% demean_fr_within_trial, demean_pooledsd_within_condition

% [] means pooled all-condition model.
% Non-empty means condition-specific models, e.g. 1:16.
data_condition = [];

runIdx = 1;

% Time-window parameters, in seconds.
% Example:
% binsize = 0.02, nbin = 20
% analysis_window = 0.2, sliding_step = 0.2
% gives windows bin 1:10 and bin 11:20.
%
% If analysis_window = nbin * binsize and only one window is generated,
% filenames are kept exactly the same as the old static version.
analysis_window = 0.2;
sliding_step    = 0.2;

% Fields to analyze from seqEst.
%
% Example:
% analysis_fields = {'y'};
% analysis_fields = {'y', 'yRecon_use_all', 'yRecon_use_across'};
%
% Original data:
%   y
%
% Base noiseless reconstruction fields after data_reconstruction.m:
%   d
%   yRecon_use_across
%   yRecon_use_within
%   yRecon_use_all
%   yRecon_use_across_no_d
%   yRecon_use_within_no_d
%   yRecon_use_all_no_d
%   yRecon_across_excl_within
%   yRecon_within_excl_across
%
% Feedback / feedforward fields, if generated:
%   yRecon_use_feedback
%   yRecon_feedback_excl_within_ff_ambiguous
%   yRecon_feedback_excl_within
%   yRecon_feedback_excl_ff_ambiguous
%
%   yRecon_use_feedforward
%   yRecon_feedforward_excl_within_fb_ambiguous
%   yRecon_feedforward_excl_within
%   yRecon_feedforward_excl_fb_ambiguous
%
% Reconstruction fields with R noise, if generated:
%   yRecon_use_across_with_R
%   yRecon_use_within_with_R
%   yRecon_use_all_with_R
%   yRecon_across_excl_within_with_R
%   yRecon_within_excl_across_with_R
%
% Reconstruction fields keeping residual, if generated:
%   yRecon_use_across_keep_resid
%   yRecon_use_within_keep_resid
%   yRecon_use_all_keep_resid
%   yRecon_across_excl_within_keep_resid
%   yRecon_within_excl_across_keep_resid

analysis_fields = { ...
    'y', ...
    'yRecon_use_across', ...
    'yRecon_use_within', ...
    'yRecon_use_all', ...
    'yRecon_across_excl_within', ...
    'yRecon_within_excl_across', ...
    'yRecon_use_feedback', ...
    'yRecon_feedback_excl_within_ff_ambiguous', ...
    'yRecon_feedback_excl_within', ...
    'yRecon_feedback_excl_ff_ambiguous', ...
    'yRecon_use_feedforward', ...
    'yRecon_feedforward_excl_within_fb_ambiguous', ...
    'yRecon_feedforward_excl_within', ...
    'yRecon_feedforward_excl_fb_ambiguous'};

dat_file = './model_data_allruns.mat';

stim_tag = '_2[Gpl2_2c_2sz_400_2_200isi]';

group_names = {};

% Save options.
save_fig = true;
save_png = true;
save_matlab_fig = true;
save_result_mat = true;

% Plot options.
plot_fullrange = false;
plot_brokenaxis = true;
plot_metric_hist = true;

% Response and metric cleanup.
% Finite values with abs(value) < tolerance are forced to exactly 0.
% response_zero_tolerance is applied before metric calculation.
% metric_zero_tolerance is applied after metric calculation.
response_zero_tolerance = 1e-10;
metric_zero_tolerance = 1e-10;

% Figure options.
fig_position = [100 100 1800 700];
marker_size = 40;
marker_face_alpha = 0.45;
marker_edge_alpha = 1;

% Broken-axis options, copied in spirit from plot_size_effect.m.
break_start_prctile = 98.0;
broken_axis_trigger_ratio = 1;
tail_display_frac = 0.08;
break_gap_frac = 0.03;

%% ----------------------- Resolve script folder and inputs -----------------------

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

analysis_fields = normalize_analysis_fields(analysis_fields);
nAnalysisFields = numel(analysis_fields);

if isempty(data_condition)
    model_mode = 'all_condition_model';
    use_condition_specific_models = false;
    condition_list = [];
else
    model_mode = 'condition_specific_models';
    use_condition_specific_models = true;
    condition_list = validate_condition_list(data_condition);
end

safe_data_content = sanitize_filename(data_content);

%% ----------------------- Load model data -----------------------

dat_file = ensure_mat_file(dat_file);

if ~isfile(dat_file)
    error('dat_file does not exist: %s', dat_file);
end

fprintf('Reading model data from:\n  %s\n', dat_file);
S_model = load(dat_file);

if ~isfield(S_model, 'model_data_allruns')
    error('model_data_allruns not found in %s.', dat_file);
end

model_data_allruns = S_model.model_data_allruns;

all_run_tags = get_all_run_tags_local(model_data_allruns);
run_data_idx = find(strcmp(all_run_tags, stim_tag));

if isempty(run_data_idx)
    error('Requested stim_tag not found: %s', stim_tag);
end

if numel(run_data_idx) > 1
    error('Duplicate stim_tag found: %s', stim_tag);
end

this_run = get_run_from_model_data(model_data_allruns, run_data_idx);

fprintf('Selected model_data_allruns index: %d\n', run_data_idx);
fprintf('Selected stim_tag: %s\n', this_run.stim_tag);

if isfield(this_run, 'conditions_full')
    conditions_full = this_run.conditions_full;
elseif isfield(this_run, 'condition_full')
    conditions_full = this_run.condition_full;
else
    error('Neither conditions_full nor condition_full was found in selected run.');
end

if isempty(conditions_full)
    error('conditions_full is empty.');
end

if ~isfield(conditions_full, 'contrast')
    error('Condition field "contrast" not found in conditions_full.');
end

if ~isfield(conditions_full, 'stim_name')
    error(['Condition field "stim_name" not found in conditions_full. ', ...
        'This script defines low/high contrast within each stim_name.']);
end

binsize = get_binsize_from_run_local(this_run);

fprintf('Using binsize from this_run.bin_size: %.12g s\n', binsize);

%% ----------------------- Determine low and high contrast -----------------------

contrast_info = build_contrast_pair_info_local(conditions_full);

fprintf('\nContrast grouping:\n');
fprintf('Low and high contrast are defined within each stim_name, then pooled.\n');

fprintf('\nContrast values by stim_name:\n');
disp(contrast_info.contrastValuesByStim);

fprintf('Low-contrast condition indices:\n');
disp(contrast_info.low_condition_indices);

fprintf('High-contrast condition indices:\n');
disp(contrast_info.high_condition_indices);

if isempty(contrast_info.low_condition_indices)
    error('No low-contrast conditions found.');
end

if isempty(contrast_info.high_condition_indices)
    error('No high-contrast conditions found.');
end

requested_condition_list = condition_list;

if use_condition_specific_models
    if any(condition_list < 1) || any(condition_list > numel(conditions_full))
        error('data_condition contains indices outside 1:%d.', numel(conditions_full));
    end

    fprintf('\nRequested condition-specific models:\n');
    disp(requested_condition_list);
else
    fprintf('\nUsing pooled all-condition model.\n');
end

fprintf('\nAnalysis fields:\n');
disp(analysis_fields(:));

%% ----------------------- Initialize per-field accumulation cache -----------------------

field_cache = initialize_field_cache(analysis_fields);

groupd = [];
nGroups = [];
group_names_this = {};
group_row_ranges = {};
output_dir = '';

time_window_info = [];

%% ----------------------- Load model(s) once per needed model and accumulate all fields -----------------------

if ~use_condition_specific_models

    %% ----------------------- All-condition pooled model: load once -----------------------

    this_condition = [];
    [baseDir, runDir, bestmodel_file] = resolve_bestmodel_from_training_settings( ...
        data_content, this_condition, runIdx);

    output_dir = runDir;

    fprintf('\nUsing baseDir:\n  %s\n', baseDir);
    fprintf('Using runDir:\n  %s\n', runDir);
    fprintf('Loading bestmodel once for all analysis fields:\n  %s\n', bestmodel_file);

    [S_best, seqEst] = load_bestmodel_once(bestmodel_file);

    nbin = get_nbin_from_seqEst_local(seqEst);
    time_window_info = build_time_window_info_local( ...
        binsize, nbin, analysis_window, sliding_step);

    print_time_window_info_local(time_window_info);

    trial_ids = extract_trial_ids(seqEst);
    condition_index_seq = map_seq_trials_to_conditions( ...
        this_run, conditions_full, trial_ids);

    if any(condition_index_seq < 1) || any(condition_index_seq > numel(conditions_full))
        error('condition_index_seq contains invalid condition indices.');
    end

    trial_contrast_code = contrast_info.condition_contrast_code(condition_index_seq);

    low_trial_mask = trial_contrast_code == 1;
    high_trial_mask = trial_contrast_code == 2;

    low_trial_indices_in_seqEst = find(low_trial_mask);
    high_trial_indices_in_seqEst = find(high_trial_mask);

    for f = 1:nAnalysisFields
        analysis_field = analysis_fields{f};

        fprintf('\nReading seqEst.%s from loaded all-condition model\n', analysis_field);

        [nUnits_this, nTimeBins_this] = validate_seq_field_shape(seqEst, analysis_field);

        if f == 1
            groupd = get_groupd(S_best, this_run, data_content, nUnits_this);
            groupd = groupd(:)';

            if sum(groupd) ~= nUnits_this
                error('sum(groupd) = %d, but seqEst.%s has %d units.', ...
                    sum(groupd), analysis_field, nUnits_this);
            end

            nGroups = numel(groupd);
            group_names_this = normalize_group_names(group_names, nGroups, this_run);
            [~, ~, group_row_ranges] = build_group_index(groupd);
        else
            if nUnits_this ~= sum(groupd)
                error('seqEst.%s has %d units, expected %d from groupd.', ...
                    analysis_field, nUnits_this, sum(groupd));
            end
        end

        trial_response_by_window = compute_trial_response_from_seqEst( ...
            seqEst, analysis_field, nUnits_this, time_window_info);

        field_cache(f).nUnits = nUnits_this;
        field_cache(f).nTimeBins = nTimeBins_this;
        field_cache(f).low_trial_response_by_window = cell(1, time_window_info.nWindows);
        field_cache(f).high_trial_response_by_window = cell(1, time_window_info.nWindows);

        for w = 1:time_window_info.nWindows
            field_cache(f).low_trial_response_by_window{w} = ...
                trial_response_by_window{w}(:, low_trial_mask);
            field_cache(f).high_trial_response_by_window{w} = ...
                trial_response_by_window{w}(:, high_trial_mask);
        end

        field_cache(f).low_trial_condition_index = condition_index_seq(low_trial_mask)';
        field_cache(f).high_trial_condition_index = condition_index_seq(high_trial_mask)';
        field_cache(f).low_trial_indices_in_seqEst = low_trial_indices_in_seqEst;
        field_cache(f).high_trial_indices_in_seqEst = high_trial_indices_in_seqEst;

        field_cache(f).model_source(1).model_mode = model_mode;
        field_cache(f).model_source(1).condition = [];
        field_cache(f).model_source(1).baseDir = baseDir;
        field_cache(f).model_source(1).runDir = runDir;
        field_cache(f).model_source(1).bestmodel_file = bestmodel_file;
        field_cache(f).model_source(1).nTrials = numel(seqEst);
    end

    clear S_best seqEst;

else

    %% ----------------------- Condition-specific models: load one condition at a time -----------------------

    output_dir = scriptDir;
    first_condition_loaded = false;

    for cc = 1:numel(condition_list)
        this_condition = condition_list(cc);

        this_condition_contrast_code = contrast_info.condition_contrast_code(this_condition);
        this_condition_contrast_label = contrast_info.condition_contrast_label{this_condition};
        this_condition_contrast_value = contrast_info.condition_contrast_value(this_condition);
        this_condition_stim_name = contrast_info.condition_stim_name{this_condition};

        [baseDir, runDir, bestmodel_file] = resolve_bestmodel_from_training_settings( ...
            data_content, this_condition, runIdx);

        fprintf('\nCondition %d/%d: condition %d, contrast = %g, label = %s, stim_name = %s\n', ...
            cc, numel(condition_list), this_condition, ...
            this_condition_contrast_value, this_condition_contrast_label, ...
            this_condition_stim_name);

        fprintf('Using baseDir:\n  %s\n', baseDir);
        fprintf('Using runDir:\n  %s\n', runDir);
        fprintf('Loading bestmodel once for all analysis fields:\n  %s\n', bestmodel_file);

        [S_best, seqEst] = load_bestmodel_once(bestmodel_file);

        if ~first_condition_loaded
            nbin = get_nbin_from_seqEst_local(seqEst);
            time_window_info = build_time_window_info_local( ...
                binsize, nbin, analysis_window, sliding_step);

            print_time_window_info_local(time_window_info);
        end

        for f = 1:nAnalysisFields
            analysis_field = analysis_fields{f};

            fprintf('  Reading seqEst.%s\n', analysis_field);

            [nUnits_this, nTimeBins_this] = validate_seq_field_shape(seqEst, analysis_field);

            if ~first_condition_loaded && f == 1
                groupd = get_groupd(S_best, this_run, data_content, nUnits_this);
                groupd = groupd(:)';

                if sum(groupd) ~= nUnits_this
                    error('sum(groupd) = %d, but seqEst.%s has %d units.', ...
                        sum(groupd), analysis_field, nUnits_this);
                end

                nGroups = numel(groupd);
                group_names_this = normalize_group_names(group_names, nGroups, this_run);
                [~, ~, group_row_ranges] = build_group_index(groupd);
            else
                if nUnits_this ~= sum(groupd)
                    error('Condition %d seqEst.%s has %d units, expected %d from groupd.', ...
                        this_condition, analysis_field, nUnits_this, sum(groupd));
                end

                groupd_this = get_groupd(S_best, this_run, data_content, nUnits_this);
                groupd_this = groupd_this(:)';

                if ~isequal(groupd_this, groupd)
                    error('Condition %d has groupd %s, expected %s.', ...
                        this_condition, mat2str(groupd_this), mat2str(groupd));
                end
            end

            if isempty(field_cache(f).nUnits)
                field_cache(f).nUnits = nUnits_this;
                field_cache(f).nTimeBins = nTimeBins_this;
                field_cache(f).low_trial_response_by_window = cell(1, time_window_info.nWindows);
                field_cache(f).high_trial_response_by_window = cell(1, time_window_info.nWindows);

                for w = 1:time_window_info.nWindows
                    field_cache(f).low_trial_response_by_window{w} = [];
                    field_cache(f).high_trial_response_by_window{w} = [];
                end
            else
                if field_cache(f).nUnits ~= nUnits_this || ...
                        field_cache(f).nTimeBins ~= nTimeBins_this
                    error(['Condition %d has seqEst.%s size %d x %d, ', ...
                        'but previous conditions used %d x %d.'], ...
                        this_condition, analysis_field, nUnits_this, nTimeBins_this, ...
                        field_cache(f).nUnits, field_cache(f).nTimeBins);
                end
            end

            trial_response_by_window = compute_trial_response_from_seqEst( ...
                seqEst, analysis_field, nUnits_this, time_window_info);

            n_trials_this = size(trial_response_by_window{1}, 2);

            if this_condition_contrast_code == 1
                for w = 1:time_window_info.nWindows
                    field_cache(f).low_trial_response_by_window{w} = [ ...
                        field_cache(f).low_trial_response_by_window{w}, ...
                        trial_response_by_window{w}]; %#ok<AGROW>
                end

                field_cache(f).low_trial_condition_index = [ ...
                    field_cache(f).low_trial_condition_index, ...
                    repmat(this_condition, 1, n_trials_this)]; %#ok<AGROW>

            elseif this_condition_contrast_code == 2
                for w = 1:time_window_info.nWindows
                    field_cache(f).high_trial_response_by_window{w} = [ ...
                        field_cache(f).high_trial_response_by_window{w}, ...
                        trial_response_by_window{w}]; %#ok<AGROW>
                end

                field_cache(f).high_trial_condition_index = [ ...
                    field_cache(f).high_trial_condition_index, ...
                    repmat(this_condition, 1, n_trials_this)]; %#ok<AGROW>

            else
                error('Condition %d has unexpected contrast code %g.', ...
                    this_condition, this_condition_contrast_code);
            end

            field_cache(f).model_source(cc).model_mode = model_mode; %#ok<SAGROW>
            field_cache(f).model_source(cc).condition = this_condition;
            field_cache(f).model_source(cc).condition_contrast_code = this_condition_contrast_code;
            field_cache(f).model_source(cc).condition_contrast_label = this_condition_contrast_label;
            field_cache(f).model_source(cc).condition_contrast_value = this_condition_contrast_value;
            field_cache(f).model_source(cc).condition_stim_name = this_condition_stim_name;
            field_cache(f).model_source(cc).baseDir = baseDir;
            field_cache(f).model_source(cc).runDir = runDir;
            field_cache(f).model_source(cc).bestmodel_file = bestmodel_file;
            field_cache(f).model_source(cc).nTrials = n_trials_this;
        end

        first_condition_loaded = true;

        clear S_best seqEst;
    end
end

if isempty(groupd)
    error('groupd was not initialized. No model data were processed.');
end

%% ----------------------- Process each analysis field and each time window -----------------------

for analysisFieldIdx = 1:nAnalysisFields
    analysis_field = field_cache(analysisFieldIdx).analysis_field;
    safe_field = sanitize_filename(analysis_field);

    fprintf('\n============================================================\n');
    fprintf('Processing accumulated analysis field %d/%d: seqEst.%s\n', ...
        analysisFieldIdx, nAnalysisFields, analysis_field);
    fprintf('============================================================\n');

    for windowIdx = 1:time_window_info.nWindows

        window_suffix = time_window_info.file_suffix{windowIdx};
        window_label = time_window_info.label{windowIdx};

        contrast_effect_base_name = sprintf('%s_%s_contrast_effect_%s%s', ...
            safe_data_content, model_mode, safe_field, window_suffix);

        lowvshigh_fullrange_base_name = sprintf('%s_%s_lowvshigh_%s_fullrange%s', ...
            safe_data_content, model_mode, safe_field, window_suffix);

        lowvshigh_brokenaxis_base_name = sprintf('%s_%s_lowvshigh_%s_brokenaxis%s', ...
            safe_data_content, model_mode, safe_field, window_suffix);

        fprintf('\n------------------------------------------------------------\n');
        fprintf('Time window %d/%d: %s\n', ...
            windowIdx, time_window_info.nWindows, window_label);
        fprintf('------------------------------------------------------------\n');

        low_trial_response = ...
            field_cache(analysisFieldIdx).low_trial_response_by_window{windowIdx};
        high_trial_response = ...
            field_cache(analysisFieldIdx).high_trial_response_by_window{windowIdx};

        low_trial_condition_index = field_cache(analysisFieldIdx).low_trial_condition_index;
        high_trial_condition_index = field_cache(analysisFieldIdx).high_trial_condition_index;
        low_trial_indices_in_seqEst = field_cache(analysisFieldIdx).low_trial_indices_in_seqEst;
        high_trial_indices_in_seqEst = field_cache(analysisFieldIdx).high_trial_indices_in_seqEst;
        model_source = field_cache(analysisFieldIdx).model_source;

        nUnits = field_cache(analysisFieldIdx).nUnits;
        nTimeBins = field_cache(analysisFieldIdx).nTimeBins;

        n_low_trials = size(low_trial_response, 2);
        n_high_trials = size(high_trial_response, 2);

        if n_low_trials == 0
            error('No trials found for low contrast.');
        end

        if n_high_trials == 0
            error('No trials found for high contrast.');
        end

        fprintf('\nNumber of pooled trials:\n');
        fprintf('  low contrast: %d trials\n', n_low_trials);
        fprintf('  high contrast: %d trials\n', n_high_trials);

        fprintf('\nOutput folder:\n  %s\n', output_dir);

        %% ----------------------- Compute pooled responses by contrast -----------------------

        low_response = mean(low_trial_response, 2, 'omitnan');
        high_response = mean(high_trial_response, 2, 'omitnan');

        % Remove tiny numerical residuals in pooled responses before metric calculation.
        low_response = force_small_metric_values_to_zero(low_response, response_zero_tolerance);
        high_response = force_small_metric_values_to_zero(high_response, response_zero_tolerance);

        low_response_std = std(low_trial_response, 0, 2, 'omitnan');
        high_response_std = std(high_trial_response, 0, 2, 'omitnan');

        low_response_sem = low_response_std ./ sqrt(sum(isfinite(low_trial_response), 2));
        high_response_sem = high_response_std ./ sqrt(sum(isfinite(high_trial_response), 2));

        if any(isnan(low_response)) || any(isnan(high_response))
            warning('NaN found in computed low_response or high_response. Please inspect seqEst.%s.', ...
                analysis_field);
        end

        %% ----------------------- Compute metrics -----------------------

        valid_delta_mask = isfinite(low_response) & isfinite(high_response);

        CMI_denominator = abs(high_response) + abs(low_response);
        valid_CMI_mask = valid_delta_mask & CMI_denominator > response_zero_tolerance;

        delta_HL = nan(nUnits, 1);
        CMI = nan(nUnits, 1);

        delta_HL(valid_delta_mask) = ...
            high_response(valid_delta_mask) - low_response(valid_delta_mask);

        CMI(valid_CMI_mask) = ...
            delta_HL(valid_CMI_mask) ./ CMI_denominator(valid_CMI_mask);

        % Remove tiny numerical residuals in metrics.
        % NaN and Inf are preserved.
        delta_HL = force_small_metric_values_to_zero(delta_HL, metric_zero_tolerance);
        CMI = force_small_metric_values_to_zero(CMI, metric_zero_tolerance);

        effect = struct();
        effect.delta_HL = delta_HL;
        effect.CMI = CMI;

        valid_effect_mask = struct();
        valid_effect_mask.delta_HL = valid_delta_mask;
        valid_effect_mask.CMI = valid_CMI_mask;

        %% ----------------------- Save mat result -----------------------

        contrast_effect_result = struct();

        contrast_effect_result.data_content = data_content;
        contrast_effect_result.model_mode = model_mode;
        contrast_effect_result.data_condition = data_condition;
        contrast_effect_result.requested_condition_list = requested_condition_list;
        contrast_effect_result.condition_list = condition_list;

        contrast_effect_result.runIdx = runIdx;
        contrast_effect_result.stim_tag = stim_tag;
        contrast_effect_result.dat_file = dat_file;
        contrast_effect_result.output_dir = output_dir;

        contrast_effect_result.analysis_field = analysis_field;
        contrast_effect_result.analysis_fields = analysis_fields;
        contrast_effect_result.response_per_trial = ...
            'sum across selected time-window bins';

        contrast_effect_result.binsize = binsize;
        contrast_effect_result.nbin = time_window_info.nbin;
        contrast_effect_result.analysis_window = analysis_window;
        contrast_effect_result.sliding_step = sliding_step;
        contrast_effect_result.window_nbin = time_window_info.window_nbin;
        contrast_effect_result.step_nbin = time_window_info.step_nbin;
        contrast_effect_result.nWindows = time_window_info.nWindows;
        contrast_effect_result.time_window_index = windowIdx;
        contrast_effect_result.bin_start = time_window_info.bin_start(windowIdx);
        contrast_effect_result.bin_end = time_window_info.bin_end(windowIdx);
        contrast_effect_result.time_start = time_window_info.time_start(windowIdx);
        contrast_effect_result.time_end = time_window_info.time_end(windowIdx);
        contrast_effect_result.time_window_label = window_label;
        contrast_effect_result.time_file_suffix = window_suffix;
        contrast_effect_result.time_parameter_file_tag = time_window_info.parameter_file_tag;
        contrast_effect_result.is_full_trial_single_window = ...
            time_window_info.is_full_trial_single_window;

        contrast_effect_result.contrast_grouping_rule = ...
            ['Low/high contrast are defined within each stim_name by sorting ', ...
            'that stim_name''s contrast values. Then all low conditions and all ', ...
            'high conditions are pooled. Other condition features are not used ', ...
            'for subdivision.'];

        contrast_effect_result.contrastValuesByStim = contrast_info.contrastValuesByStim;
        contrast_effect_result.stimLabelsForContrast = contrast_info.stimLabels;

        contrast_effect_result.condition_contrast_code = contrast_info.condition_contrast_code;
        contrast_effect_result.condition_contrast_label = contrast_info.condition_contrast_label;
        contrast_effect_result.condition_contrast_value = contrast_info.condition_contrast_value;
        contrast_effect_result.condition_stim_name = contrast_info.condition_stim_name;

        contrast_effect_result.low_condition_indices = contrast_info.low_condition_indices;
        contrast_effect_result.high_condition_indices = contrast_info.high_condition_indices;

        contrast_effect_result.low_trial_indices_in_seqEst = low_trial_indices_in_seqEst;
        contrast_effect_result.high_trial_indices_in_seqEst = high_trial_indices_in_seqEst;
        contrast_effect_result.low_trial_condition_index = low_trial_condition_index;
        contrast_effect_result.high_trial_condition_index = high_trial_condition_index;

        contrast_effect_result.n_low_trials = n_low_trials;
        contrast_effect_result.n_high_trials = n_high_trials;

        contrast_effect_result.nUnits = nUnits;
        contrast_effect_result.nTimeBins = nTimeBins;
        contrast_effect_result.groupd = groupd;
        contrast_effect_result.group_names = group_names_this;
        contrast_effect_result.model_source = model_source;

        contrast_effect_result.metric_formulas.delta_HL = 'H - L';
        contrast_effect_result.metric_formulas.CMI = '(H - L) ./ (abs(H) + abs(L))';

        contrast_effect_result.response_zero_tolerance = response_zero_tolerance;
        contrast_effect_result.metric_zero_tolerance = metric_zero_tolerance;
        contrast_effect_result.zero_tolerance_rule = ...
            ['Finite low_response and high_response values with abs(value) < ', ...
            'response_zero_tolerance are forced to exactly 0 before metric ', ...
            'calculation. Finite metric values with abs(value) < ', ...
            'metric_zero_tolerance are forced to exactly 0 after metric calculation. ', ...
            'No epsilon is added to CMI. CMI is computed only when ', ...
            'abs(high_response) + abs(low_response) > response_zero_tolerance.'];

        contrast_effect_result.low_response = low_response;
        contrast_effect_result.high_response = high_response;
        contrast_effect_result.low_response_std = low_response_std;
        contrast_effect_result.high_response_std = high_response_std;
        contrast_effect_result.low_response_sem = low_response_sem;
        contrast_effect_result.high_response_sem = high_response_sem;

        contrast_effect_result.delta_HL = delta_HL;
        contrast_effect_result.CMI = CMI;
        contrast_effect_result.valid_effect_mask = valid_effect_mask;

        contrast_effect_result.low_trial_response = low_trial_response;
        contrast_effect_result.high_trial_response = high_trial_response;

        for g = 1:nGroups
            rows = group_row_ranges{g};

            contrast_effect_result.group(g).group_name = group_names_this{g};
            contrast_effect_result.group(g).group_index = g;
            contrast_effect_result.group(g).nUnits = numel(rows);

            contrast_effect_result.group(g).low_response = low_response(rows);
            contrast_effect_result.group(g).high_response = high_response(rows);
            contrast_effect_result.group(g).low_response_std = low_response_std(rows);
            contrast_effect_result.group(g).high_response_std = high_response_std(rows);
            contrast_effect_result.group(g).low_response_sem = low_response_sem(rows);
            contrast_effect_result.group(g).high_response_sem = high_response_sem(rows);

            contrast_effect_result.group(g).delta_HL = delta_HL(rows);
            contrast_effect_result.group(g).CMI = CMI(rows);
        end

        output_mat = fullfile(output_dir, sprintf('%s.mat', contrast_effect_base_name));

        if save_result_mat
            save(output_mat, 'contrast_effect_result', '-v7.3');
            fprintf('\nSaved result mat:\n  %s\n', output_mat);
        end

        %% ----------------------- Plot 1: full-range linear axis -----------------------

        if plot_fullrange
            hfig_full = figure( ...
                'Color', 'w', ...
                'Name', 'Low vs high contrast response by group, full range', ...
                'Position', fig_position);

            plot_contrast_effect_groups_fullrange( ...
                low_response, high_response, groupd, ...
                stim_tag, analysis_field, model_mode, window_label, ...
                marker_size, marker_face_alpha, marker_edge_alpha);

            if save_fig
                save_current_figure_local(hfig_full, output_dir, ...
                    lowvshigh_fullrange_base_name, ...
                    save_png, save_matlab_fig);
            end
        end

        %% ----------------------- Plot 2: clean broken-axis display -----------------------

        if plot_brokenaxis
            hfig_broken = figure( ...
                'Color', 'w', ...
                'Name', 'Low vs high contrast response by group, clean broken axis', ...
                'Position', fig_position);

            plot_contrast_effect_groups_clean_brokenaxis( ...
                low_response, high_response, groupd, ...
                stim_tag, analysis_field, model_mode, window_label, ...
                marker_size, marker_face_alpha, marker_edge_alpha, ...
                break_start_prctile, broken_axis_trigger_ratio, ...
                tail_display_frac, break_gap_frac);

            if save_fig
                save_current_figure_local(hfig_broken, output_dir, ...
                    lowvshigh_brokenaxis_base_name, ...
                    save_png, save_matlab_fig);
            end
        end

        %% ----------------------- Plot 3: metric histograms -----------------------

        if plot_metric_hist
            hfig_hist = figure( ...
                'Color', 'w', ...
                'Name', 'Contrast effect metric histograms by group', ...
                'Position', fig_position);

            plot_metric_histograms_by_group( ...
                effect, valid_effect_mask, groupd, ...
                stim_tag, analysis_field, model_mode, window_label);

            if save_fig
                save_current_figure_local(hfig_hist, output_dir, ...
                    contrast_effect_base_name, ...
                    save_png, save_matlab_fig);
            end
        end

        fprintf('\nFinished seqEst.%s, %s\n', analysis_field, window_label);

        close all;
    end
end

fprintf('\nDone.\n');

%% =========================================================================
% Local functions
% =========================================================================

function fname = ensure_mat_file(fname)
if isstring(fname)
    fname = char(fname);
end

if isfile(fname)
    return;
end

[p, n, e] = fileparts(fname);

if isempty(e)
    fname = fullfile(p, [n, '.mat']);
end
end

function analysis_fields = normalize_analysis_fields(analysis_fields)
if ischar(analysis_fields)
    analysis_fields = {analysis_fields};
elseif isstring(analysis_fields)
    analysis_fields = cellstr(analysis_fields(:))';
elseif iscell(analysis_fields)
    analysis_fields = reshape(analysis_fields, 1, []);

    for i = 1:numel(analysis_fields)
        if isstring(analysis_fields{i})
            analysis_fields{i} = char(analysis_fields{i});
        end

        if ~ischar(analysis_fields{i})
            error('Each entry in analysis_fields must be a char or string.');
        end
    end
else
    error('analysis_fields must be a char, string, string array, or cell array of chars.');
end

empty_mask = cellfun(@isempty, analysis_fields);

if any(empty_mask)
    error('analysis_fields contains empty entries.');
end

if numel(unique(analysis_fields, 'stable')) ~= numel(analysis_fields)
    error('analysis_fields contains duplicate fields.');
end
end

function contrast_info = build_contrast_pair_info_local(conditions_full)
% Build per-condition relative contrast labels.
%
% Important:
% - Low/high contrast is defined within each stim_name.
% - This avoids using the global min/max contrast, which can be wrong
%   when different stim_name values have different numeric contrast levels.
% - After labeling, this script pools all low conditions and all high
%   conditions across stim_name.

nCond = numel(conditions_full);

required_fields = {'stim_name', 'contrast'};

for f = 1:numel(required_fields)
    fn = required_fields{f};

    if ~isfield(conditions_full, fn)
        error('conditions_full.%s is required for contrast effect analysis.', fn);
    end
end

stimNameAll = strings(nCond, 1);
contrastAll = nan(nCond, 1);

for k = 1:nCond
    stimNameAll(k) = lower(string(conditions_full(k).stim_name));
    contrastAll(k) = conditions_full(k).contrast;

    if ~isfinite(contrastAll(k))
        error('conditions_full(%d).contrast is not finite.', k);
    end
end

allStim = unique(stimNameAll, 'stable');
allStim = lower(allStim);

% Keep the grating/plaid order when both are present, following the
% condition-summary convention used in the size-effect script.
if all(ismember(["grating", "plaid"], allStim))
    stimLabels = ["grating", "plaid"];
else
    stimLabels = allStim(:)';
end

contrastValuesByStim = struct();
condition_contrast_code = nan(1, nCond);
condition_contrast_label = cell(1, nCond);
condition_contrast_value = nan(1, nCond);
condition_stim_name = cell(1, nCond);

for s = 1:numel(stimLabels)
    stim = stimLabels(s);
    idx = (stimNameAll == stim);

    if ~any(idx)
        continue;
    end

    cvals = unique(contrastAll(idx));
    cvals = sort(cvals(:)');

    if numel(cvals) ~= 2
        error(['Stim %s does not have exactly 2 contrast levels. ', ...
            'Found %d levels: %s.'], ...
            char(stim), numel(cvals), mat2str(cvals));
    end

    contrastValuesByStim.(char(stim)) = cvals;

    condIdx = find(idx);

    for jj = 1:numel(condIdx)
        condID = condIdx(jj);
        currContrast = contrastAll(condID);

        contrastCode = find(abs(cvals - currContrast) < 1e-10, 1);

        if isempty(contrastCode)
            error('Could not map condition %d contrast value %g within stim %s.', ...
                condID, currContrast, char(stim));
        end

        condition_contrast_code(condID) = contrastCode;
        condition_contrast_label{condID} = ternary_label_local(contrastCode, 'low', 'high');
        condition_contrast_value(condID) = currContrast;
        condition_stim_name{condID} = char(stim);
    end
end

if any(isnan(condition_contrast_code))
    missingID = find(isnan(condition_contrast_code), 1);
    error('Could not assign contrast code for condition %d.', missingID);
end

contrast_info = struct();
contrast_info.condition_contrast_code = condition_contrast_code;
contrast_info.condition_contrast_label = condition_contrast_label;
contrast_info.condition_contrast_value = condition_contrast_value;
contrast_info.condition_stim_name = condition_stim_name;
contrast_info.low_condition_indices = find(condition_contrast_code == 1);
contrast_info.high_condition_indices = find(condition_contrast_code == 2);
contrast_info.contrastValuesByStim = contrastValuesByStim;
contrast_info.stimLabels = cellstr(stimLabels);
end

function label = ternary_label_local(code, label1, label2)
if code == 1
    label = label1;
elseif code == 2
    label = label2;
else
    error('ternary_label_local expects code 1 or 2.');
end
end

function cache = initialize_field_cache(analysis_fields)
nFields = numel(analysis_fields);

cache = repmat(struct( ...
    'analysis_field', '', ...
    'low_trial_response_by_window', {{}}, ...
    'high_trial_response_by_window', {{}}, ...
    'low_trial_condition_index', [], ...
    'high_trial_condition_index', [], ...
    'low_trial_indices_in_seqEst', [], ...
    'high_trial_indices_in_seqEst', [], ...
    'model_source', struct([]), ...
    'nUnits', [], ...
    'nTimeBins', []), 1, nFields);

for f = 1:nFields
    cache(f).analysis_field = analysis_fields{f};
end
end

function condition_list = validate_condition_list(data_condition)
condition_list = data_condition;

if isstring(condition_list)
    condition_list = str2double(condition_list);
end

if ~isnumeric(condition_list)
    error('data_condition must be [] or a numeric vector of condition indices.');
end

condition_list = reshape(condition_list, 1, []);

if isempty(condition_list)
    return;
end

if any(~isfinite(condition_list)) || any(mod(condition_list, 1) ~= 0)
    error('data_condition must contain finite integer condition indices.');
end

if any(condition_list < 1)
    error('data_condition must contain positive condition indices.');
end

if numel(unique(condition_list, 'stable')) ~= numel(condition_list)
    error('data_condition contains duplicate condition indices.');
end
end

function binsize = get_binsize_from_run_local(this_run)
if isfield(this_run, 'bin_size')
    binsize = this_run.bin_size;
else
    error('this_run.bin_size was not found.');
end

if ~isnumeric(binsize) || ~isscalar(binsize) || ~isfinite(binsize) || binsize <= 0
    error('this_run.bin_size must be a finite positive scalar in seconds.');
end
end

function nbin = get_nbin_from_seqEst_local(seqEst)
if isempty(seqEst)
    error('seqEst is empty.');
end

if ~isfield(seqEst, 'T')
    error('seqEst(1).T was not found.');
end

nbin = seqEst(1).T;

if ~isnumeric(nbin) || ~isscalar(nbin) || ~isfinite(nbin) || nbin <= 0 || mod(nbin, 1) ~= 0
    error('seqEst(1).T must be a finite positive integer.');
end

nbin = double(nbin);
end

function time_window_info = build_time_window_info_local( ...
    binsize, nbin, analysis_window, sliding_step)

if ~isnumeric(analysis_window) || ~isscalar(analysis_window) || ...
        ~isfinite(analysis_window) || analysis_window <= 0
    error('analysis_window must be a finite positive scalar in seconds.');
end

if ~isnumeric(sliding_step) || ~isscalar(sliding_step) || ...
        ~isfinite(sliding_step) || sliding_step <= 0
    error('sliding_step must be a finite positive scalar in seconds.');
end

window_nbin_raw = analysis_window / binsize;
step_nbin_raw = sliding_step / binsize;

window_nbin = round(window_nbin_raw);
step_nbin = round(step_nbin_raw);

if abs(window_nbin_raw - window_nbin) > 1e-9 || ...
        abs(step_nbin_raw - step_nbin) > 1e-9
    error(['analysis_window and sliding_step must be integer multiples ', ...
        'of this_run.bin_size. binsize = %.12g, analysis_window = %.12g, sliding_step = %.12g.'], ...
        binsize, analysis_window, sliding_step);
end

if window_nbin > nbin
    error('analysis_window is longer than the trial. window_nbin = %d, nbin = %d.', ...
        window_nbin, nbin);
end

bin_start = 1:step_nbin:(nbin - window_nbin + 1);
bin_end = bin_start + window_nbin - 1;

if isempty(bin_start)
    error('No valid time windows were generated.');
end

nWindows = numel(bin_start);

time_start = (bin_start - 1) * binsize;
time_end = bin_end * binsize;

parameter_file_tag = sprintf('_wn%s_st%s', ...
    format_seconds_for_filename_local(analysis_window), ...
    format_seconds_for_filename_local(sliding_step));

is_full_trial_single_window = ...
    nWindows == 1 && bin_start(1) == 1 && bin_end(1) == nbin;

file_suffix = cell(1, nWindows);
label = cell(1, nWindows);

for w = 1:nWindows
    if is_full_trial_single_window
        file_suffix{w} = '';
    else
        file_suffix{w} = sprintf('%s_w%02d', parameter_file_tag, w);
    end

    label{w} = sprintf('window %d/%d, %.6g-%.6g s, bins %d-%d', ...
        w, nWindows, time_start(w), time_end(w), bin_start(w), bin_end(w));
end

time_window_info = struct();
time_window_info.binsize = binsize;
time_window_info.nbin = nbin;
time_window_info.analysis_window = analysis_window;
time_window_info.sliding_step = sliding_step;
time_window_info.window_nbin = window_nbin;
time_window_info.step_nbin = step_nbin;
time_window_info.nWindows = nWindows;
time_window_info.bin_start = bin_start;
time_window_info.bin_end = bin_end;
time_window_info.time_start = time_start;
time_window_info.time_end = time_end;
time_window_info.parameter_file_tag = parameter_file_tag;
time_window_info.file_suffix = file_suffix;
time_window_info.label = label;
time_window_info.is_full_trial_single_window = is_full_trial_single_window;
end

function print_time_window_info_local(time_window_info)
fprintf('\nTime-window analysis settings:\n');
fprintf('  binsize         = %.12g s\n', time_window_info.binsize);
fprintf('  nbin            = %d\n', time_window_info.nbin);
fprintf('  analysis_window = %.12g s = %d bins\n', ...
    time_window_info.analysis_window, time_window_info.window_nbin);
fprintf('  sliding_step    = %.12g s = %d bins\n', ...
    time_window_info.sliding_step, time_window_info.step_nbin);
fprintf('  nWindows        = %d\n', time_window_info.nWindows);

if time_window_info.is_full_trial_single_window
    fprintf('  Full-trial single window detected: old filenames will be used.\n');
else
    fprintf('  Filename parameter tag: %s\n', time_window_info.parameter_file_tag);
end

for w = 1:time_window_info.nWindows
    fprintf('  %s, file_suffix = %s\n', ...
        time_window_info.label{w}, time_window_info.file_suffix{w});
end
end

function s = format_seconds_for_filename_local(x)
s = sprintf('%.12g', x);
s = strrep(s, '.', '_');
s = strrep(s, '-', 'm');
s = regexprep(s, '[^a-zA-Z0-9_]', '_');
s = regexprep(s, '_+', '_');
end

function run_data = get_run_from_model_data(model_data_allruns, idx)
if iscell(model_data_allruns)
    run_data = model_data_allruns{idx};
else
    run_data = model_data_allruns(idx);
end
end

function all_tags = get_all_run_tags_local(model_data_allruns)
all_tags = cell(numel(model_data_allruns), 1);

for j = 1:numel(model_data_allruns)
    run_data = get_run_from_model_data(model_data_allruns, j);

    if ~isfield(run_data, 'stim_tag')
        error('stim_tag missing in model_data_allruns entry %d.', j);
    end

    all_tags{j} = run_data.stim_tag;
end
end

function [baseDir, runDir, bestmodel_file] = resolve_bestmodel_from_training_settings( ...
    data_content, this_condition, runIdx)

if isempty(this_condition)
    baseDir = ['./FA_Dlag_', data_content];
else
    baseDir = ['./FA_Dlag_', data_content, '_condition', num2str(this_condition)];
end

runDir = sprintf('%s/mat_results/run%03d', baseDir, runIdx);

files = dir(fullfile(runDir, 'bestmodel*'));

if isempty(files)
    error('No bestmodel* file found in %s', runDir);
end

[~, newestIdx] = max([files.datenum]);
bestmodel_file = fullfile(runDir, files(newestIdx).name);
end

function [S_best, seqEst] = load_bestmodel_once(bestmodel_file)
S_best = load(bestmodel_file);

if ~isfield(S_best, 'seqEst')
    error('seqEst not found in bestmodel file: %s', bestmodel_file);
end

seqEst = S_best.seqEst;

if isempty(seqEst)
    error('seqEst is empty in bestmodel file: %s', bestmodel_file);
end
end

function [nUnits, nTimeBins] = validate_seq_field_shape(seqEst, analysis_field)
if ~isfield(seqEst, analysis_field)
    error('Field seqEst.%s not found. Choose another analysis_field.', ...
        analysis_field);
end

y0 = seqEst(1).(analysis_field);

if ~isnumeric(y0) || ndims(y0) ~= 2
    error('seqEst(1).%s must be a numeric nUnit x nTimeBin matrix.', ...
        analysis_field);
end

[nUnits, nTimeBins] = size(y0);

for t = 1:numel(seqEst)
    if ~isfield(seqEst(t), analysis_field)
        error('seqEst(%d).%s is missing.', t, analysis_field);
    end

    y = seqEst(t).(analysis_field);

    if ~isnumeric(y) || ndims(y) ~= 2
        error('seqEst(%d).%s must be a numeric nUnit x nTimeBin matrix.', ...
            t, analysis_field);
    end

    if size(y, 1) ~= nUnits || size(y, 2) ~= nTimeBins
        error('seqEst(%d).%s size mismatch. Expected %d x %d, got %d x %d.', ...
            t, analysis_field, nUnits, nTimeBins, size(y, 1), size(y, 2));
    end
end
end

function trial_response_by_window = compute_trial_response_from_seqEst( ...
    seqEst, analysis_field, nUnits, time_window_info)

nTrials = numel(seqEst);
nWindows = time_window_info.nWindows;

trial_response_by_window = cell(1, nWindows);

for w = 1:nWindows
    trial_response_by_window{w} = nan(nUnits, nTrials);
end

for t = 1:nTrials
    y = seqEst(t).(analysis_field);

    if size(y, 1) ~= nUnits
        error('Unit number mismatch in seqEst trial %d.', t);
    end

    for w = 1:nWindows
        bin_start = time_window_info.bin_start(w);
        bin_end = time_window_info.bin_end(w);

        trial_response_by_window{w}(:, t) = ...
            sum(y(:, bin_start:bin_end), 2, 'omitnan');
    end
end
end

function groupd = get_groupd(S_best, this_run, data_content, nUnits)
if isfield(S_best, 'res') && ...
        isfield(S_best.res, 'estParams') && ...
        isfield(S_best.res.estParams, 'yDims')
    groupd = S_best.res.estParams.yDims;

elseif isfield(S_best, 'bestModel') && ...
        isfield(S_best.bestModel, 'estParams') && ...
        isfield(S_best.bestModel.estParams, 'yDims')
    groupd = S_best.bestModel.estParams.yDims;

elseif isfield(S_best, 'bestModel') && ...
        isfield(S_best.bestModel, 'yDims')
    groupd = S_best.bestModel.yDims;

else
    groupd_field = sprintf('%s_groupd', data_content);

    if isfield(this_run, 'nan_trial_strategy') && ...
            this_run.nan_trial_strategy == 6 && ...
            isfield(this_run, groupd_field)
        groupd = this_run.(groupd_field);
    elseif isfield(this_run, 'groupd')
        groupd = this_run.groupd;
    else
        warning('Could not find groupd/yDims. Treating all cells as one group.');
        groupd = nUnits;
    end
end

groupd = groupd(:)';
end

function group_names = normalize_group_names(group_names, nGroups, this_run)
if isempty(group_names)
    candidate_fields = {'group_names', 'area_names', 'group_name', 'areas'};

    for f = 1:numel(candidate_fields)
        fn = candidate_fields{f};

        if isfield(this_run, fn)
            candidate = this_run.(fn);

            if isstring(candidate)
                candidate = cellstr(candidate);
            end

            if iscell(candidate) && numel(candidate) == nGroups
                group_names = candidate;
                break;
            end
        end
    end
end

if isempty(group_names)
    group_names = cell(1, nGroups);

    for g = 1:nGroups
        group_names{g} = sprintf('Group %d', g);
    end
else
    if isstring(group_names)
        group_names = cellstr(group_names);
    end

    if ~iscell(group_names) || numel(group_names) ~= nGroups
        error('group_names must be empty or a cell array with nGroups entries.');
    end

    for g = 1:nGroups
        group_names{g} = char(group_names{g});
    end
end
end

function [group_index_all, unit_index_within_group_all, group_row_ranges] = build_group_index(groupd)
nUnits = sum(groupd);

group_index_all = nan(nUnits, 1);
unit_index_within_group_all = nan(nUnits, 1);
group_row_ranges = cell(1, numel(groupd));

group_start = 1;

for g = 1:numel(groupd)
    group_end = group_start + groupd(g) - 1;
    rows = group_start:group_end;

    group_row_ranges{g} = rows;
    group_index_all(rows) = g;
    unit_index_within_group_all(rows) = (1:groupd(g))';

    group_start = group_end + 1;
end
end

function trial_ids = extract_trial_ids(seqEst)
if isfield(seqEst, 'trialId')
    trial_ids = nan(numel(seqEst), 1);

    for t = 1:numel(seqEst)
        trial_ids(t) = seqEst(t).trialId;
    end
else
    trial_ids = (1:numel(seqEst))';
end

trial_ids = trial_ids(:);
end

function condition_index_seq = map_seq_trials_to_conditions(this_run, conditions_full, trial_ids)
trial_ids = trial_ids(:);

if isfield(this_run, 'condition_index_per_trial_full')
    condition_index_per_trial_full = this_run.condition_index_per_trial_full(:);

    if any(trial_ids < 1) || any(trial_ids > numel(condition_index_per_trial_full))
        error('Some seqEst trialId values are outside condition_index_per_trial_full range.');
    end

    condition_index_seq = condition_index_per_trial_full(trial_ids);
    return;
end

if ~isfield(conditions_full, 'trial_indices')
    error(['Cannot map seqEst trials to conditions. Need either ', ...
        'this_run.condition_index_per_trial_full or conditions_full(k).trial_indices.']);
end

maxTrialIndex = 0;

for c = 1:numel(conditions_full)
    idx = conditions_full(c).trial_indices(:);

    if ~isempty(idx)
        maxTrialIndex = max(maxTrialIndex, max(idx));
    end
end

if maxTrialIndex < 1
    error('No valid trial indices found in conditions_full.trial_indices.');
end

trial_to_condition = nan(maxTrialIndex, 1);

for c = 1:numel(conditions_full)
    idx = conditions_full(c).trial_indices(:);

    if isempty(idx)
        continue;
    end

    if any(idx < 1) || any(mod(idx, 1) ~= 0)
        error('conditions_full(%d).trial_indices contains invalid values.', c);
    end

    already_assigned = ~isnan(trial_to_condition(idx));

    if any(already_assigned)
        dupIdx = idx(find(already_assigned, 1));
        error('Trial index %d appears in multiple conditions.', dupIdx);
    end

    trial_to_condition(idx) = c;
end

if any(trial_ids < 1) || any(trial_ids > numel(trial_to_condition))
    error('Some seqEst trialId values are outside conditions_full.trial_indices range.');
end

condition_index_seq = trial_to_condition(trial_ids);

if any(isnan(condition_index_seq))
    missingTrial = trial_ids(find(isnan(condition_index_seq), 1));
    error('Could not map seqEst trialId %d to any condition.', missingTrial);
end
end

function plot_contrast_effect_groups_fullrange( ...
    low_response, high_response, groupd, ...
    stim_tag, analysis_field, model_mode, window_label, ...
    marker_size, marker_face_alpha, marker_edge_alpha)

nGroups = numel(groupd);
group_start = 1;

for g = 1:nGroups
    group_end = group_start + groupd(g) - 1;
    group_idx = group_start:group_end;

    x = low_response(group_idx);
    y = high_response(group_idx);

    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    subplot(1, nGroups, g);

    scatter(x, y, marker_size, 'filled', ...
        'MarkerFaceAlpha', marker_face_alpha, ...
        'MarkerEdgeAlpha', marker_edge_alpha);

    hold on;

    all_val = [x(:); y(:)];
    all_val = all_val(isfinite(all_val));

    if isempty(all_val)
        min_val = -1;
        max_val = 1;
    else
        min_val = min(all_val);
        max_val = max(all_val);

        if min_val == max_val
            min_val = min_val - 1;
            max_val = max_val + 1;
        end
    end

    pad = 0.05 * max(eps, max_val - min_val);
    axis_min = min_val - pad;
    axis_max = max_val + pad;

    plot([axis_min axis_max], [axis_min axis_max], 'k--', 'LineWidth', 1);

    axis square;
    xlim([axis_min axis_max]);
    ylim([axis_min axis_max]);

    xlabel('Low contrast: response in window');
    ylabel('High contrast: response in window');
    title(sprintf('Group %d, n = %d cells', g, groupd(g)), ...
        'Interpreter', 'none');

    grid off;
    box off;

    group_start = group_end + 1;
end

sgtitle(sprintf('%s, %s, seqEst.%s, %s', ...
    stim_tag, model_mode, analysis_field, window_label), ...
    'Interpreter', 'none');
end

function plot_contrast_effect_groups_clean_brokenaxis( ...
    low_response, high_response, groupd, ...
    stim_tag, analysis_field, model_mode, window_label, ...
    marker_size, marker_face_alpha, marker_edge_alpha, ...
    break_start_prctile, broken_axis_trigger_ratio, ...
    tail_display_frac, break_gap_frac)

nGroups = numel(groupd);
group_start = 1;

for g = 1:nGroups
    group_end = group_start + groupd(g) - 1;
    group_idx = group_start:group_end;

    x = low_response(group_idx);
    y = high_response(group_idx);

    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    subplot(1, nGroups, g);

    all_val = [x(:); y(:)];
    all_val = all_val(isfinite(all_val));

    if isempty(all_val)
        all_val = 1;
    end

    min_val = min(all_val);
    max_val = max(all_val);

    if min_val < 0
        scatter(x, y, marker_size, 'filled', ...
            'MarkerFaceAlpha', marker_face_alpha, ...
            'MarkerEdgeAlpha', marker_edge_alpha);

        hold on;

        if min_val == max_val
            min_val = min_val - 1;
            max_val = max_val + 1;
        end

        pad = 0.05 * max(eps, max_val - min_val);
        axis_min = min_val - pad;
        axis_max = max_val + pad;

        plot([axis_min axis_max], [axis_min axis_max], 'k--', 'LineWidth', 1);

        axis square;
        xlim([axis_min axis_max]);
        ylim([axis_min axis_max]);

        xlabel('Low contrast: response in window');
        ylabel('High contrast: response in window');
        title(sprintf('Group %d, n = %d cells, linear axis', g, groupd(g)), ...
            'Interpreter', 'none');

        grid off;
        box off;
    else
        if isempty(max_val) || ~isfinite(max_val) || max_val <= 0
            max_val = 1;
        end

        break_start = prctile(all_val, break_start_prctile);

        if isempty(break_start) || ~isfinite(break_start) || break_start <= 0
            break_start = max_val;
        end

        use_broken_axis = max_val > break_start * broken_axis_trigger_ratio;

        if ~use_broken_axis
            scatter(x, y, marker_size, 'filled', ...
                'MarkerFaceAlpha', marker_face_alpha, ...
                'MarkerEdgeAlpha', marker_edge_alpha);

            hold on;

            plot([0 max_val], [0 max_val], 'k--', 'LineWidth', 1);

            axis square;
            xlim([0 max_val * 1.05]);
            ylim([0 max_val * 1.05]);

            raw_ticks = choose_prebreak_integer_ticks_local(max_val);

            xticks(raw_ticks);
            yticks(raw_ticks);
            xticklabels(compose_integer_tick_labels_local(raw_ticks));
            yticklabels(compose_integer_tick_labels_local(raw_ticks));

            xlabel('Low contrast: response in window');
            ylabel('High contrast: response in window');
            title(sprintf('Group %d, n = %d cells, linear axis', g, groupd(g)), ...
                'Interpreter', 'none');

            grid off;
            box off;
        else
            high_values = all_val(all_val > break_start);
            break_end = min(high_values);

            if isempty(break_end) || ~isfinite(break_end) || break_end <= break_start
                break_end = break_start;
            end

            tail_display_len = max(eps, break_start * tail_display_frac);
            break_gap = max(eps, break_start * break_gap_frac);

            x_plot = broken_axis_transform_local( ...
                x, break_start, break_end, max_val, tail_display_len, break_gap);
            y_plot = broken_axis_transform_local( ...
                y, break_start, break_end, max_val, tail_display_len, break_gap);

            scatter(x_plot, y_plot, marker_size, 'filled', ...
                'MarkerFaceAlpha', marker_face_alpha, ...
                'MarkerEdgeAlpha', marker_edge_alpha);

            hold on;

            low_ref_raw = linspace(0, break_start, 100);
            low_ref_plot = broken_axis_transform_local( ...
                low_ref_raw, break_start, break_end, max_val, tail_display_len, break_gap);
            plot(low_ref_plot, low_ref_plot, 'k--', 'LineWidth', 1);

            high_ref_raw = linspace(break_end, max_val, 100);
            high_ref_plot = broken_axis_transform_local( ...
                high_ref_raw, break_start, break_end, max_val, tail_display_len, break_gap);
            plot(high_ref_plot, high_ref_plot, 'k--', 'LineWidth', 1);

            display_max = broken_axis_transform_local( ...
                max_val, break_start, break_end, max_val, tail_display_len, break_gap);

            axis square;
            xlim([0 display_max * 1.05]);
            ylim([0 display_max * 1.05]);

            raw_ticks = choose_prebreak_integer_ticks_local(break_start);
            plot_ticks = broken_axis_transform_local( ...
                raw_ticks, break_start, break_end, max_val, tail_display_len, break_gap);

            xticks(plot_ticks);
            yticks(plot_ticks);
            xticklabels(compose_integer_tick_labels_local(raw_ticks));
            yticklabels(compose_integer_tick_labels_local(raw_ticks));

            xlabel('Low contrast: response in window');
            ylabel('High contrast: response in window');
            title(sprintf('Group %d, n = %d cells', g, groupd(g)), ...
                'Interpreter', 'none');

            grid off;
            box off;

            draw_axis_break_marks_local(gca, break_start, break_gap, display_max);
        end
    end

    group_start = group_end + 1;
end

sgtitle(sprintf('%s, %s, seqEst.%s, %s', ...
    stim_tag, model_mode, analysis_field, window_label), ...
    'Interpreter', 'none');
end

function plot_metric_histograms_by_group(effect, valid_effect_mask, groupd, stim_tag, analysis_field, model_mode, window_label)
metric_names = {'delta_HL', 'CMI'};
nGroups = numel(groupd);

for m = 1:numel(metric_names)
    metric_name = metric_names{m};
    metric_values = effect.(metric_name);
    metric_valid = valid_effect_mask.(metric_name);

    for g = 1:nGroups
        plot_idx = (m - 1) * nGroups + g;
        subplot(numel(metric_names), nGroups, plot_idx);

        group_start = sum(groupd(1:g-1)) + 1;
        group_end = sum(groupd(1:g));
        group_idx = group_start:group_end;

        vals = metric_values(group_idx);
        valid = metric_valid(group_idx) & isfinite(vals);
        vals = vals(valid);

        histogram(vals, 30);

        grid off;
        box off;

        xlabel(metric_label(metric_name), 'Interpreter', 'none');
        ylabel('Cell count');

        if isempty(vals)
            title(sprintf('Group %d, n = 0', g), 'Interpreter', 'none');
        else
            title(sprintf('Group %d, n = %d, median = %.3f', ...
                g, numel(vals), median(vals, 'omitnan')), ...
                'Interpreter', 'none');
        end
    end
end

sgtitle(sprintf('%s, %s, seqEst.%s, contrast effect metrics, %s', ...
    stim_tag, model_mode, analysis_field, window_label), ...
    'Interpreter', 'none');
end

function v_plot = broken_axis_transform_local( ...
    v, break_start, break_end, max_val, tail_display_len, break_gap)

v_plot = v;

if max_val <= break_end
    return;
end

high_mask = v >= break_end;
middle_mask = v > break_start & v < break_end;

v_plot(middle_mask) = break_start;
v_plot(high_mask) = break_start + break_gap + ...
    (v(high_mask) - break_end) ./ max(eps, max_val - break_end) .* tail_display_len;
end

function raw_ticks = choose_prebreak_integer_ticks_local(prebreak_max)
if isempty(prebreak_max) || ~isfinite(prebreak_max) || prebreak_max <= 0
    raw_ticks = 0;
    return;
end

max_tick = floor(prebreak_max);

if max_tick < 1
    raw_ticks = 0;
    return;
end

target_intervals = 5;
rough_step = max_tick / target_intervals;
step = choose_nice_integer_step_local(rough_step);

last_tick = floor(max_tick / step) * step;
raw_ticks = 0:step:last_tick;

if numel(raw_ticks) < 3 && max_tick >= 2
    step = 1;
    raw_ticks = 0:step:max_tick;
end
end

function step = choose_nice_integer_step_local(rough_step)
if rough_step <= 1
    step = 1;
    return;
end

exponent = floor(log10(rough_step));
base = 10 ^ exponent;
candidates = [1 2 5 10] * base;
candidates = candidates(candidates >= rough_step);

if isempty(candidates)
    step = 10 * base;
else
    step = candidates(1);
end

step = max(1, round(step));
end

function labels = compose_integer_tick_labels_local(raw_ticks)
labels = cell(size(raw_ticks));

for i = 1:numel(raw_ticks)
    labels{i} = sprintf('%d', round(raw_ticks(i)));
end
end

function draw_axis_break_marks_local(ax, break_start, break_gap, display_max)
axes(ax); %#ok<LAXES>
hold(ax, 'on');

x_break_center = break_start + break_gap / 2;
y_break_center = break_start + break_gap / 2;

slash_dx = display_max * 0.012;
slash_dy = display_max * 0.025;
offset = display_max * 0.018;

x0 = ax.XLim(1);
y0 = ax.YLim(1);

plot(ax, ...
    [x_break_center - slash_dx, x_break_center + slash_dx], ...
    [y0 + slash_dy, y0 - slash_dy], ...
    'k-', 'LineWidth', 1.2, 'Clipping', 'off');

plot(ax, ...
    [x_break_center - slash_dx + offset, x_break_center + slash_dx + offset], ...
    [y0 + slash_dy, y0 - slash_dy], ...
    'k-', 'LineWidth', 1.2, 'Clipping', 'off');

plot(ax, ...
    [x0 + slash_dx, x0 - slash_dx], ...
    [y_break_center - slash_dy, y_break_center + slash_dy], ...
    'k-', 'LineWidth', 1.2, 'Clipping', 'off');

plot(ax, ...
    [x0 + slash_dx, x0 - slash_dx], ...
    [y_break_center - slash_dy + offset, y_break_center + slash_dy + offset], ...
    'k-', 'LineWidth', 1.2, 'Clipping', 'off');
end

function save_current_figure_local(hfig, fig_out_dir, base_name, save_png, save_matlab_fig)
if ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

if save_png
    png_file = fullfile(fig_out_dir, sprintf('%s.png', base_name));
    exportgraphics(hfig, png_file, 'Resolution', 300);
    fprintf('\nSaved PNG figure:\n  %s\n', png_file);
end

if save_matlab_fig
    fig_file = fullfile(fig_out_dir, sprintf('%s.fig', base_name));
    savefig(hfig, fig_file);
    fprintf('Saved MATLAB figure:\n  %s\n', fig_file);
end
end

function label = metric_label(metric_name)
switch metric_name
    case 'delta_HL'
        label = 'delta H-L = H - L';
    case 'CMI'
        label = 'CMI = (H - L) / (abs(H) + abs(L))';
    otherwise
        label = metric_name;
end
end

function x = force_small_metric_values_to_zero(x, tol)
if isempty(tol)
    return;
end

if ~isscalar(tol) || ~isnumeric(tol) || ~isfinite(tol) || tol < 0
    error('zero tolerance must be a finite nonnegative scalar.');
end

finite_small_mask = isfinite(x) & abs(x) < tol;
x(finite_small_mask) = 0;
end

function name = sanitize_filename(name)
if isstring(name)
    name = char(name);
end

name = regexprep(name, '[^a-zA-Z0-9_\-]', '_');
name = regexprep(name, '_+', '_');
name = strtrim(name);

if isempty(name)
    name = 'unnamed';
end
end