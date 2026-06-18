%% =========================================================================
% plot_size_effect
%
% Purpose:
% For one selected stim_tag/run:
% 1. Load seqEst from DLAG bestmodel.
% 2. Use selected seqEst field:
%       response_data = seqEst.(analysis_field)
% 3. Optionally filter conditions/trials by relative contrast level, then
%    pool trials by condition size:
%       small size vs large size
% 4. For each trial, compute total response across all time bins:
%       trial_response = sum(seqEst(t).(analysis_field), 2)
% 5. Average trial_response across all small-size trials and large-size trials.
% 6. Compute size effect metrics for each cell:
%       classic_SI  = (S - L) ./ S
%       delta_SL    = S - L
%       S_norm_diff = (S - L) ./ abs(S)
% 7. Plot small vs large response separately for each group.
% 8. Save result mat and figures.
%
% Contrast filter:
% - pick_contrast = 0
%       Use all contrasts, exactly matching the original behavior.
% - pick_contrast = 1
%       Use only low-contrast conditions. Low is defined within each
%       stim_name by sorting that stim_name's contrast values.
% - pick_contrast = 2
%       Use only high-contrast conditions. High is defined within each
%       stim_name by sorting that stim_name's contrast values.
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
%       Because size effect is across conditions, this mode first pools
%       responses from all requested condition-specific models by size,
%       then computes one size effect result.
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

% Contrast filter for size-effect pooling.
%   0 = use all contrasts, same as old behavior.
%   1 = use low contrast only.
%   2 = use high contrast only.
% Low/high are defined within each stim_name, not by global min/max contrast.
pick_contrast = 1;

runIdx = 1;

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
%   seqEst(n).d
%   seqEst(n).yRecon_use_across
%   seqEst(n).yRecon_use_within
%   seqEst(n).yRecon_use_all
%   seqEst(n).yRecon_use_across_no_d
%   seqEst(n).yRecon_use_within_no_d
%   seqEst(n).yRecon_use_all_no_d
%   seqEst(n).yRecon_across_excl_within
%   seqEst(n).yRecon_within_excl_across
%
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
    'd'...
     'y', ...
    'yRecon_use_across', ...
    'yRecon_use_within', ...
    'yRecon_use_all', ...
    'yRecon_across_excl_within', ...
    'yRecon_within_excl_across'};

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

% Broken-axis options, copied in spirit from plot_size_rawcount_scatter.m.
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
pick_contrast = validate_pick_contrast_local(pick_contrast);

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

if ~isfield(conditions_full, 'size')
    error('Condition field "size" not found in conditions_full.');
end

%% ----------------------- Determine contrast filter -----------------------

contrast_filter = build_contrast_filter_info_local(conditions_full, pick_contrast);

fprintf('\nContrast filter mode: %s\n', contrast_filter.label);
if pick_contrast ~= 0
    fprintf('Selected contrast condition indices:\n');
    disp(contrast_filter.selected_condition_indices);
    fprintf('Contrast values by stim_name:\n');
    disp(contrast_filter.contrastValuesByStim);
end

%% ----------------------- Determine small and large size -----------------------

size_values = [conditions_full.size];
unique_sizes = unique(size_values);

fprintf('\nUnique size values found in this run:\n');
disp(unique_sizes);

if numel(unique_sizes) ~= 2
    error(['Expected exactly 2 unique size values for small vs large, ', ...
        'but found %d. Please inspect conditions_full.size.'], numel(unique_sizes));
end

small_size = min(unique_sizes);
large_size = max(unique_sizes);

small_condition_indices = find(size_values == small_size & contrast_filter.condition_mask);
large_condition_indices = find(size_values == large_size & contrast_filter.condition_mask);

fprintf('Using small_size = %g\n', small_size);
fprintf('Using large_size = %g\n', large_size);

fprintf('Small-size condition indices:\n');
disp(small_condition_indices);

fprintf('Large-size condition indices:\n');
disp(large_condition_indices);

if isempty(small_condition_indices)
    error('No small-size conditions remain after applying pick_contrast = %d.', pick_contrast);
end

if isempty(large_condition_indices)
    error('No large-size conditions remain after applying pick_contrast = %d.', pick_contrast);
end

requested_condition_list = condition_list;

if use_condition_specific_models
    if any(condition_list < 1) || any(condition_list > numel(conditions_full))
        error('data_condition contains indices outside 1:%d.', numel(conditions_full));
    end

    if pick_contrast ~= 0
        condition_list = condition_list(contrast_filter.condition_mask(condition_list));

        if isempty(condition_list)
            error(['No requested condition-specific models remain after applying ', ...
                'pick_contrast = %d. Requested conditions were: %s.'], ...
                pick_contrast, mat2str(requested_condition_list));
        end
    end

    fprintf('\nRequested condition-specific models before contrast filtering:\n');
    disp(requested_condition_list);
    fprintf('Using condition-specific models after contrast filtering:\n');
    disp(condition_list);
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

    trial_ids = extract_trial_ids(seqEst);
    condition_index_seq = map_seq_trials_to_conditions( ...
        this_run, conditions_full, trial_ids);

    if any(condition_index_seq < 1) || any(condition_index_seq > numel(conditions_full))
        error('condition_index_seq contains invalid condition indices.');
    end

    trial_size = size_values(condition_index_seq);
    trial_contrast_mask = contrast_filter.condition_mask(condition_index_seq);

    small_trial_mask = trial_size == small_size & trial_contrast_mask;
    large_trial_mask = trial_size == large_size & trial_contrast_mask;

    small_trial_indices_in_seqEst = find(small_trial_mask);
    large_trial_indices_in_seqEst = find(large_trial_mask);

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

        trial_response = compute_trial_response_from_seqEst(seqEst, analysis_field, nUnits_this);

        field_cache(f).nUnits = nUnits_this;
        field_cache(f).nTimeBins = nTimeBins_this;
        field_cache(f).small_trial_response = trial_response(:, small_trial_mask);
        field_cache(f).large_trial_response = trial_response(:, large_trial_mask);
        field_cache(f).small_trial_condition_index = condition_index_seq(small_trial_mask)';
        field_cache(f).large_trial_condition_index = condition_index_seq(large_trial_mask)';
        field_cache(f).small_trial_indices_in_seqEst = small_trial_indices_in_seqEst;
        field_cache(f).large_trial_indices_in_seqEst = large_trial_indices_in_seqEst;

        field_cache(f).model_source(1).model_mode = model_mode;
        field_cache(f).model_source(1).condition = [];
        field_cache(f).model_source(1).baseDir = baseDir;
        field_cache(f).model_source(1).runDir = runDir;
        field_cache(f).model_source(1).bestmodel_file = bestmodel_file;
        field_cache(f).model_source(1).nTrials = numel(seqEst);
        field_cache(f).model_source(1).pick_contrast = pick_contrast;
        field_cache(f).model_source(1).contrast_label = contrast_filter.label;
    end

    clear S_best seqEst;

else

    %% ----------------------- Condition-specific models: load one condition at a time -----------------------

    output_dir = scriptDir;
    first_condition_loaded = false;

    for cc = 1:numel(condition_list)

        this_condition = condition_list(cc);
        this_condition_size = size_values(this_condition);

        [baseDir, runDir, bestmodel_file] = resolve_bestmodel_from_training_settings( ...
            data_content, this_condition, runIdx);

        fprintf('\nCondition %d/%d: condition %d, size = %g\n', ...
            cc, numel(condition_list), this_condition, this_condition_size);
        fprintf('Using baseDir:\n  %s\n', baseDir);
        fprintf('Using runDir:\n  %s\n', runDir);
        fprintf('Loading bestmodel once for all analysis fields:\n  %s\n', bestmodel_file);

        [S_best, seqEst] = load_bestmodel_once(bestmodel_file);

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
            else
                if field_cache(f).nUnits ~= nUnits_this || field_cache(f).nTimeBins ~= nTimeBins_this
                    error(['Condition %d has seqEst.%s size %d x %d, ', ...
                        'but previous conditions used %d x %d.'], ...
                        this_condition, analysis_field, nUnits_this, nTimeBins_this, ...
                        field_cache(f).nUnits, field_cache(f).nTimeBins);
                end
            end

            trial_response = compute_trial_response_from_seqEst(seqEst, analysis_field, nUnits_this);
            n_trials_this = size(trial_response, 2);

            if this_condition_size == small_size
                field_cache(f).small_trial_response = [ ...
                    field_cache(f).small_trial_response, trial_response]; %#ok<AGROW>
                field_cache(f).small_trial_condition_index = [ ...
                    field_cache(f).small_trial_condition_index, ...
                    repmat(this_condition, 1, n_trials_this)]; %#ok<AGROW>

            elseif this_condition_size == large_size
                field_cache(f).large_trial_response = [ ...
                    field_cache(f).large_trial_response, trial_response]; %#ok<AGROW>
                field_cache(f).large_trial_condition_index = [ ...
                    field_cache(f).large_trial_condition_index, ...
                    repmat(this_condition, 1, n_trials_this)]; %#ok<AGROW>

            else
                error('Condition %d has unexpected size value %g.', ...
                    this_condition, this_condition_size);
            end

            field_cache(f).model_source(cc).model_mode = model_mode; %#ok<SAGROW>
            field_cache(f).model_source(cc).condition = this_condition;
            field_cache(f).model_source(cc).condition_size = this_condition_size;
            field_cache(f).model_source(cc).condition_contrast_code = contrast_filter.condition_contrast_code(this_condition);
            field_cache(f).model_source(cc).condition_contrast_label = contrast_filter.condition_contrast_label{this_condition};
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

%% ----------------------- Process each analysis field from accumulated data -----------------------

for analysisFieldIdx = 1:nAnalysisFields

    analysis_field = field_cache(analysisFieldIdx).analysis_field;
    safe_field = sanitize_filename(analysis_field);

    size_effect_base_name = sprintf('%s_%s_size_effect_%s%s', ...
        safe_data_content, model_mode, safe_field, contrast_filter.file_suffix);

    svsl_fullrange_base_name = sprintf('%s_%s_svsl_%s_fullrange%s', ...
        safe_data_content, model_mode, safe_field, contrast_filter.file_suffix);

    svsl_brokenaxis_base_name = sprintf('%s_%s_svsl_%s_brokenaxis%s', ...
        safe_data_content, model_mode, safe_field, contrast_filter.file_suffix);

    fprintf('\n============================================================\n');
    fprintf('Processing accumulated analysis field %d/%d: seqEst.%s\n', ...
        analysisFieldIdx, nAnalysisFields, analysis_field);
    fprintf('============================================================\n');

    small_trial_response = field_cache(analysisFieldIdx).small_trial_response;
    large_trial_response = field_cache(analysisFieldIdx).large_trial_response;

    small_trial_condition_index = field_cache(analysisFieldIdx).small_trial_condition_index;
    large_trial_condition_index = field_cache(analysisFieldIdx).large_trial_condition_index;

    small_trial_indices_in_seqEst = field_cache(analysisFieldIdx).small_trial_indices_in_seqEst;
    large_trial_indices_in_seqEst = field_cache(analysisFieldIdx).large_trial_indices_in_seqEst;

    model_source = field_cache(analysisFieldIdx).model_source;

    nUnits = field_cache(analysisFieldIdx).nUnits;
    nTimeBins = field_cache(analysisFieldIdx).nTimeBins;

    n_small_trials = size(small_trial_response, 2);
    n_large_trials = size(large_trial_response, 2);

    if n_small_trials == 0
        error('No trials found for small_size = %g.', small_size);
    end

    if n_large_trials == 0
        error('No trials found for large_size = %g.', large_size);
    end

    fprintf('\nNumber of pooled trials:\n');
    fprintf('  small size %g: %d trials\n', small_size, n_small_trials);
    fprintf('  large size %g: %d trials\n', large_size, n_large_trials);

    fprintf('\nOutput folder:\n  %s\n', output_dir);

    %% ----------------------- Compute pooled responses by size -----------------------

    S_response = mean(small_trial_response, 2, 'omitnan');
    L_response = mean(large_trial_response, 2, 'omitnan');

    % Remove tiny numerical residuals in pooled responses before metric calculation.
    % This is important because S_response is used as denominator.
    S_response = force_small_metric_values_to_zero(S_response, response_zero_tolerance);
    L_response = force_small_metric_values_to_zero(L_response, response_zero_tolerance);

    S_response_std = std(small_trial_response, 0, 2, 'omitnan');
    L_response_std = std(large_trial_response, 0, 2, 'omitnan');

    S_response_sem = S_response_std ./ sqrt(sum(isfinite(small_trial_response), 2));
    L_response_sem = L_response_std ./ sqrt(sum(isfinite(large_trial_response), 2));

    if any(isnan(S_response)) || any(isnan(L_response))
        warning('NaN found in computed S_response or L_response. Please inspect seqEst.%s.', analysis_field);
    end

    %% ----------------------- Compute metrics -----------------------

    valid_delta_mask = isfinite(S_response) & isfinite(L_response);
    valid_denominator_mask = valid_delta_mask & abs(S_response) > response_zero_tolerance;

    classic_SI = nan(nUnits, 1);
    delta_SL = nan(nUnits, 1);
    S_norm_diff = nan(nUnits, 1);

    delta_SL(valid_delta_mask) = S_response(valid_delta_mask) - L_response(valid_delta_mask);
    classic_SI(valid_denominator_mask) = delta_SL(valid_denominator_mask) ./ S_response(valid_denominator_mask);
    S_norm_diff(valid_denominator_mask) = delta_SL(valid_denominator_mask) ./ abs(S_response(valid_denominator_mask));

    % Remove tiny numerical residuals in metrics.
    % NaN and Inf are preserved.
    delta_SL = force_small_metric_values_to_zero(delta_SL, metric_zero_tolerance);
    classic_SI = force_small_metric_values_to_zero(classic_SI, metric_zero_tolerance);
    S_norm_diff = force_small_metric_values_to_zero(S_norm_diff, metric_zero_tolerance);

    effect = struct();
    effect.classic_SI = classic_SI;
    effect.delta_SL = delta_SL;
    effect.S_norm_diff = S_norm_diff;

    valid_effect_mask = struct();
    valid_effect_mask.classic_SI = valid_denominator_mask;
    valid_effect_mask.delta_SL = valid_delta_mask;
    valid_effect_mask.S_norm_diff = valid_denominator_mask;

    valid_metric_names = fieldnames(effect); %#ok<NASGU>

    %% ----------------------- Save mat result -----------------------

    size_effect_result = struct();

    size_effect_result.data_content = data_content;
    size_effect_result.model_mode = model_mode;
    size_effect_result.data_condition = data_condition;
    size_effect_result.requested_condition_list = requested_condition_list;
    size_effect_result.condition_list = condition_list;

    size_effect_result.pick_contrast = pick_contrast;
    size_effect_result.contrast_label = contrast_filter.label;
    size_effect_result.contrast_file_suffix = contrast_filter.file_suffix;
    size_effect_result.contrast_condition_mask = contrast_filter.condition_mask;
    size_effect_result.selected_contrast_condition_indices = contrast_filter.selected_condition_indices;
    size_effect_result.condition_contrast_code = contrast_filter.condition_contrast_code;
    size_effect_result.condition_contrast_label = contrast_filter.condition_contrast_label;
    size_effect_result.condition_contrast_value = contrast_filter.condition_contrast_value;
    size_effect_result.condition_stim_name = contrast_filter.condition_stim_name;
    size_effect_result.contrastValuesByStim = contrast_filter.contrastValuesByStim;
    size_effect_result.stimLabelsForContrast = contrast_filter.stimLabels;
    size_effect_result.runIdx = runIdx;
    size_effect_result.stim_tag = stim_tag;
    size_effect_result.dat_file = dat_file;
    size_effect_result.output_dir = output_dir;
    size_effect_result.analysis_field = analysis_field;
    size_effect_result.analysis_fields = analysis_fields;

    size_effect_result.response_per_trial = 'sum across all time bins';

    size_effect_result.small_size = small_size;
    size_effect_result.large_size = large_size;
    size_effect_result.small_condition_indices = small_condition_indices;
    size_effect_result.large_condition_indices = large_condition_indices;

    size_effect_result.small_trial_indices_in_seqEst = small_trial_indices_in_seqEst;
    size_effect_result.large_trial_indices_in_seqEst = large_trial_indices_in_seqEst;
    size_effect_result.small_trial_condition_index = small_trial_condition_index;
    size_effect_result.large_trial_condition_index = large_trial_condition_index;

    size_effect_result.n_small_trials = n_small_trials;
    size_effect_result.n_large_trials = n_large_trials;
    size_effect_result.nUnits = nUnits;
    size_effect_result.nTimeBins = nTimeBins;
    size_effect_result.groupd = groupd;
    size_effect_result.group_names = group_names_this;
    size_effect_result.model_source = model_source;

    size_effect_result.metric_formulas.classic_SI = '(S - L) ./ S';
    size_effect_result.metric_formulas.delta_SL = 'S - L';
    size_effect_result.metric_formulas.S_norm_diff = '(S - L) ./ abs(S)';

    size_effect_result.response_zero_tolerance = response_zero_tolerance;
    size_effect_result.metric_zero_tolerance = metric_zero_tolerance;
    size_effect_result.zero_tolerance_rule = ...
        ['Finite S_response and L_response values with abs(value) < response_zero_tolerance ', ...
         'are forced to exactly 0 before metric calculation. Finite metric values with ', ...
         'abs(value) < metric_zero_tolerance are forced to exactly 0 after metric calculation.'];

    size_effect_result.S_response = S_response;
    size_effect_result.L_response = L_response;
    size_effect_result.S_response_std = S_response_std;
    size_effect_result.L_response_std = L_response_std;
    size_effect_result.S_response_sem = S_response_sem;
    size_effect_result.L_response_sem = L_response_sem;

    size_effect_result.classic_SI = classic_SI;
    size_effect_result.delta_SL = delta_SL;
    size_effect_result.S_norm_diff = S_norm_diff;
    size_effect_result.valid_effect_mask = valid_effect_mask;

    size_effect_result.small_trial_response = small_trial_response;
    size_effect_result.large_trial_response = large_trial_response;

    for g = 1:nGroups
        rows = group_row_ranges{g};

        size_effect_result.group(g).group_name = group_names_this{g};
        size_effect_result.group(g).group_index = g;
        size_effect_result.group(g).nUnits = numel(rows);
        size_effect_result.group(g).S_response = S_response(rows);
        size_effect_result.group(g).L_response = L_response(rows);
        size_effect_result.group(g).S_response_std = S_response_std(rows);
        size_effect_result.group(g).L_response_std = L_response_std(rows);
        size_effect_result.group(g).S_response_sem = S_response_sem(rows);
        size_effect_result.group(g).L_response_sem = L_response_sem(rows);
        size_effect_result.group(g).classic_SI = classic_SI(rows);
        size_effect_result.group(g).delta_SL = delta_SL(rows);
        size_effect_result.group(g).S_norm_diff = S_norm_diff(rows);
    end

    output_mat = fullfile(output_dir, sprintf('%s.mat', size_effect_base_name));

    if save_result_mat
        save(output_mat, 'size_effect_result', '-v7.3');
        fprintf('\nSaved result mat:\n  %s\n', output_mat);
    end

    %% ----------------------- Plot 1: full-range linear axis -----------------------

    if plot_fullrange
        hfig_full = figure('Color', 'w', ...
            'Name', 'Small vs large response by group, full range', ...
            'Position', fig_position);

        plot_size_effect_groups_fullrange( ...
            S_response, L_response, groupd, ...
            small_size, large_size, stim_tag, analysis_field, model_mode, ...
            marker_size, marker_face_alpha, marker_edge_alpha);

        if save_fig
            save_current_figure_local(hfig_full, output_dir, ...
                svsl_fullrange_base_name, ...
                save_png, save_matlab_fig);
        end
    end

    %% ----------------------- Plot 2: clean broken-axis display -----------------------

    if plot_brokenaxis
        hfig_broken = figure('Color', 'w', ...
            'Name', 'Small vs large response by group, clean broken axis', ...
            'Position', fig_position);

        plot_size_effect_groups_clean_brokenaxis( ...
            S_response, L_response, groupd, ...
            small_size, large_size, stim_tag, analysis_field, model_mode, ...
            marker_size, marker_face_alpha, marker_edge_alpha, ...
            break_start_prctile, broken_axis_trigger_ratio, ...
            tail_display_frac, break_gap_frac);

        if save_fig
            save_current_figure_local(hfig_broken, output_dir, ...
                svsl_brokenaxis_base_name, ...
                save_png, save_matlab_fig);
        end
    end

    %% ----------------------- Plot 3: metric histograms -----------------------

    if plot_metric_hist
        hfig_hist = figure('Color', 'w', ...
            'Name', 'Size effect metric histograms by group', ...
            'Position', fig_position);

        plot_metric_histograms_by_group( ...
            effect, valid_effect_mask, groupd, ...
            stim_tag, analysis_field, model_mode);

        if save_fig
            save_current_figure_local(hfig_hist, output_dir, ...
                size_effect_base_name, ...
                save_png, save_matlab_fig);
        end
    end

    fprintf('\nFinished seqEst.%s\n', analysis_field);

    close all;

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

function pick_contrast = validate_pick_contrast_local(pick_contrast)
    if isstring(pick_contrast)
        pick_contrast = str2double(pick_contrast);
    end

    if ~isnumeric(pick_contrast) || ~isscalar(pick_contrast) || ~isfinite(pick_contrast)
        error('pick_contrast must be 0, 1, or 2.');
    end

    if mod(pick_contrast, 1) ~= 0 || ~ismember(pick_contrast, [0 1 2])
        error('pick_contrast must be 0, 1, or 2.');
    end

    pick_contrast = double(pick_contrast);
end

function contrast_filter = build_contrast_filter_info_local(conditions_full, pick_contrast)
% Build per-condition relative contrast labels.
%
% Important:
% - Low/high contrast is defined within each stim_name.
% - This avoids using the global min/max contrast, which is wrong when
%   grating and plaid have different numeric contrast levels.
% - pick_contrast = 0 keeps all conditions and does not require contrast
%   fields, preserving the old behavior.

    nCond = numel(conditions_full);

    contrast_filter = struct();
    contrast_filter.pick_contrast = pick_contrast;
    contrast_filter.condition_mask = true(1, nCond);
    contrast_filter.selected_condition_indices = 1:nCond;
    contrast_filter.condition_contrast_code = nan(1, nCond);
    contrast_filter.condition_contrast_label = repmat({'all'}, 1, nCond);
    contrast_filter.condition_contrast_value = nan(1, nCond);
    contrast_filter.condition_stim_name = repmat({''}, 1, nCond);
    contrast_filter.contrastValuesByStim = struct();
    contrast_filter.stimLabels = {};

    if pick_contrast == 0
        contrast_filter.label = 'all_contrasts';
        contrast_filter.file_suffix = '';
        return;
    elseif pick_contrast == 1
        contrast_filter.label = 'low_contrast';
        contrast_filter.file_suffix = '_low_contrast';
    elseif pick_contrast == 2
        contrast_filter.label = 'high_contrast';
        contrast_filter.file_suffix = '_high_contrast';
    else
        error('pick_contrast must be 0, 1, or 2.');
    end

    required_fields = {'stim_name', 'contrast'};
    for f = 1:numel(required_fields)
        fn = required_fields{f};
        if ~isfield(conditions_full, fn)
            error(['pick_contrast = %d requires conditions_full.%s, ', ...
                'but that field is missing.'], pick_contrast, fn);
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
    % condition-summary convention used in Latents_compare.m.
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
                'Found %d levels: %s.'], char(stim), numel(cvals), mat2str(cvals));
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

    condition_mask = condition_contrast_code == pick_contrast;

    contrast_filter.condition_mask = condition_mask;
    contrast_filter.selected_condition_indices = find(condition_mask);
    contrast_filter.condition_contrast_code = condition_contrast_code;
    contrast_filter.condition_contrast_label = condition_contrast_label;
    contrast_filter.condition_contrast_value = condition_contrast_value;
    contrast_filter.condition_stim_name = condition_stim_name;
    contrast_filter.contrastValuesByStim = contrastValuesByStim;
    contrast_filter.stimLabels = cellstr(stimLabels);
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
        'small_trial_response', [], ...
        'large_trial_response', [], ...
        'small_trial_condition_index', [], ...
        'large_trial_condition_index', [], ...
        'small_trial_indices_in_seqEst', [], ...
        'large_trial_indices_in_seqEst', [], ...
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
        error('seqEst(1).%s must be a numeric nUnit x nTimeBin matrix.', analysis_field);
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

function trial_response = compute_trial_response_from_seqEst(seqEst, analysis_field, nUnits)
    nTrials = numel(seqEst);
    trial_response = nan(nUnits, nTrials);

    for t = 1:nTrials
        y = seqEst(t).(analysis_field);

        if size(y, 1) ~= nUnits
            error('Unit number mismatch in seqEst trial %d.', t);
        end

        trial_response(:, t) = sum(y, 2, 'omitnan');
    end
end

function groupd = get_groupd(S_best, this_run, data_content, nUnits)
    if isfield(S_best, 'res') && isfield(S_best.res, 'estParams') && isfield(S_best.res.estParams, 'yDims')
        groupd = S_best.res.estParams.yDims;

    elseif isfield(S_best, 'bestModel') && isfield(S_best.bestModel, 'estParams') && isfield(S_best.bestModel.estParams, 'yDims')
        groupd = S_best.bestModel.estParams.yDims;

    elseif isfield(S_best, 'bestModel') && isfield(S_best.bestModel, 'yDims')
        groupd = S_best.bestModel.yDims;

    else
        groupd_field = sprintf('%s_groupd', data_content);

        if isfield(this_run, 'nan_trial_strategy') && this_run.nan_trial_strategy == 6 ...
                && isfield(this_run, groupd_field)
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

function plot_size_effect_groups_fullrange( ...
    S_response, L_response, groupd, ...
    small_size, large_size, stim_tag, analysis_field, model_mode, ...
    marker_size, marker_face_alpha, marker_edge_alpha)

    nGroups = numel(groupd);
    group_start = 1;

    for g = 1:nGroups
        group_end = group_start + groupd(g) - 1;
        group_idx = group_start:group_end;

        x = S_response(group_idx);
        y = L_response(group_idx);

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

        xlabel(sprintf('Small size %g: total response per trial', small_size));
        ylabel(sprintf('Large size %g: total response per trial', large_size));
        title(sprintf('Group %d, n = %d cells', g, groupd(g)), ...
            'Interpreter', 'none');

        grid off;
        box off;

        group_start = group_end + 1;
    end

    sgtitle(sprintf('%s, %s, seqEst.%s, sum over all time bins', ...
        stim_tag, model_mode, analysis_field), ...
        'Interpreter', 'none');
end

function plot_size_effect_groups_clean_brokenaxis( ...
    S_response, L_response, groupd, ...
    small_size, large_size, stim_tag, analysis_field, model_mode, ...
    marker_size, marker_face_alpha, marker_edge_alpha, ...
    break_start_prctile, broken_axis_trigger_ratio, ...
    tail_display_frac, break_gap_frac)

    nGroups = numel(groupd);
    group_start = 1;

    for g = 1:nGroups
        group_end = group_start + groupd(g) - 1;
        group_idx = group_start:group_end;

        x = S_response(group_idx);
        y = L_response(group_idx);

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

            xlabel(sprintf('Small size %g: total response per trial', small_size));
            ylabel(sprintf('Large size %g: total response per trial', large_size));
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

                xlabel(sprintf('Small size %g: total response per trial', small_size));
                ylabel(sprintf('Large size %g: total response per trial', large_size));
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

                xlabel(sprintf('Small size %g: total response per trial', small_size));
                ylabel(sprintf('Large size %g: total response per trial', large_size));
                title(sprintf('Group %d, n = %d cells', g, groupd(g)), ...
                    'Interpreter', 'none');

                grid off;
                box off;

                draw_axis_break_marks_local(gca, break_start, break_gap, display_max);
            end
        end

        group_start = group_end + 1;
    end

    sgtitle(sprintf('%s, %s, seqEst.%s, sum over all time bins', ...
        stim_tag, model_mode, analysis_field), ...
        'Interpreter', 'none');
end

function plot_metric_histograms_by_group(effect, valid_effect_mask, groupd, stim_tag, analysis_field, model_mode)
    metric_names = {'classic_SI', 'delta_SL', 'S_norm_diff'};
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

    sgtitle(sprintf('%s, %s, seqEst.%s, size effect metrics', ...
        stim_tag, model_mode, analysis_field), ...
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
        case 'classic_SI'
            label = 'classic SI = (S - L) / S';

        case 'delta_SL'
            label = 'delta S-L = S - L';

        case 'S_norm_diff'
            label = 'S norm diff = (S - L) / abs(S)';

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