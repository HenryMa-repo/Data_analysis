%% =========================================================================
% GPL_analysis
%
% Purpose:
% Analyze a grating/plaid direction-tuning run after first restricting the
% analysis to neurons selected in a separate target run.
%
% Main steps for each probe/kilosort folder:
% 1. Read target-selected unit IDs from model_data_allruns.mat.
% 2. Find gpl_stim_tag in bined_data_allruns.mat.
% 3. Keep only target-selected units that also exist in the GPL run.
% 4. Average raw_fr across bins to obtain one mean firing rate per trial.
% 5. Compute grating and plaid direction tuning curves (mean, SD, SEM).
% 6. Optionally load a matching spontaneous/baseline data file.
% 7. Compute grating DI and vector-based OSI, with optional baseline subtraction.
% 8. Build pattern and component predictions for the plaid tuning curve.
% 9. Compute ordinary correlations, partial correlations, Fisher Z values,
%    pattern index PI = Zp - Zc, and pattern/component classification.
% 10. Save numerical results and all figures as both PNG and FIG files.
%
% Stimulus-design handling:
% - The number and spacing of grating/plaid directions are read from data.
% - The plaid component separation is detected automatically from the
%   stored plaid dir1 and dir2 values. Swapping dir1 and dir2 is allowed.
% - For detected separation alpha, the component predictor is
%       C(theta) = G(theta - alpha/2) + G(theta + alpha/2)
% - Pattern/component predictions and DI require exact measured grating
%   directions. This script never interpolates missing directions; it
%   raises an error instead.
%
% Required input in each kilosort folder:
% - bined_data_allruns.mat
%
% Required input in the shared CatGT folder:
% - model_data_allruns.mat
%
% Main output in each kilosort folder:
% - unit_gpl_results.mat
% - Group1/Group2_Grating_tuning_SEM-or-SD.png and .fig
% - Group1/Group2_Plaid_tuning_SEM-or-SD.png and .fig
% - Group1/Group2_DI_OSI_distribution.png and .fig
% - Group1/Group2_Pattern_component_analysis.png and .fig
%
% Notes:
% 1. Spontaneous subtraction is optional. The baseline file is selected by
%    stim_tag, analysis_window, and bin_size through baseline_data_index.mat.
% 2. raw_fr is unit x trial x bin. The trial response is the mean raw_fr
%    across bins, which is the average firing rate across the full window.
% 3. SD and SEM are always both saved. error_method only controls the
%    shaded error band in the tuning figures and the tuning filenames.
% 4. In spontaneous-subtracted mode, DI values above 1.6 are grouped into
%    one overflow bar labeled >1.6 in the distribution figure.
% 5. DI and OSI distribution panels show both mean and median.
% 6. Probe/group mapping follows the order below:
%       probe 0 -> Group1
%       probe 1 -> Group2
% =========================================================================

clc;
clear;

%% ----------------------- User parameters -----------------------
root_folder = 'I:\np_data';
runName = 'RafiL001p0120';
runind = 1;                         % run index after -g

probes = [0, 1];                    % probe indices after -prb
probe_group_names = {'Group1', 'Group2'};  % same order as probes

% GPL run used for grating/plaid tuning and pattern/component analysis.
gpl_stim_tag = '[dir12_gpl_2_200isi_fixedphase]';

% Whether DI, OSI, and the component predictor use spontaneous/baseline
% firing rates. Raw grating/plaid tuning figures remain on the raw FR scale.
subtract_spontaneous = true;

% Requested baseline version. The program requires one unique index entry
% matching gpl_stim_tag, this window, and this bin size.
baseline_analysis_window = [-0.1 0];
baseline_bin_size = 0.1;
baseline_index_filename = 'baseline_data_index.mat';

% Separate run used only to define the neuron whitelist.
target_stim_tag = '_2[Gpl2_2c_2sz_400_2_200isi]';

% Leave empty to use model_data_allruns.mat in the shared CatGT folder.
dat_file = '';

% Error band shown in grating/plaid tuning figures.
% Options: 'SEM' or 'SD'. Both are saved in unit_gpl_results.mat.
error_method = 'SD';

% Pattern/component classification threshold.
pattern_z_threshold = 1.28;

% Distribution binning and bar appearance.
selectivity_bin_edges = 0:0.1:1;
distribution_bar_width = 0.82;      % < 1 leaves a gap between bars

% Output MAT filename in each kilosort folder.
gpl_mat_name = 'unit_gpl_results.mat';

% Figure visibility while the script runs: 'on' or 'off'.
figure_visibility = 'on';

%% ----------------------- Validate user parameters -----------------------
if numel(probes) ~= numel(probe_group_names)
    error('probe_group_names must have one entry for each value in probes.');
end

error_method = upper(char(error_method));
if ~any(strcmp(error_method, {'SEM', 'SD'}))
    error('error_method must be ''SEM'' or ''SD''.');
end

if ~isscalar(pattern_z_threshold) || ...
        ~isnumeric(pattern_z_threshold) || ...
        ~isfinite(pattern_z_threshold) || ...
        pattern_z_threshold <= 0
    error('pattern_z_threshold must be a positive finite scalar.');
end

if ~isscalar(distribution_bar_width) || ...
        ~isnumeric(distribution_bar_width) || ...
        ~isfinite(distribution_bar_width) || ...
        distribution_bar_width <= 0 || distribution_bar_width >= 1
    error('distribution_bar_width must be a finite scalar in (0, 1).');
end

if ~islogical(subtract_spontaneous) || ~isscalar(subtract_spontaneous)
    error('subtract_spontaneous must be one logical scalar.');
end

if subtract_spontaneous
    baseline_analysis_window = double(baseline_analysis_window(:)');
    baseline_bin_size = double(baseline_bin_size);

    if numel(baseline_analysis_window) ~= 2 || ...
            any(~isfinite(baseline_analysis_window)) || ...
            baseline_analysis_window(2) <= baseline_analysis_window(1)
        error(['baseline_analysis_window must be [start end] with two ' ...
            'finite values and end > start.']);
    end

    if baseline_analysis_window(2) > 0
        error('baseline_analysis_window must end at or before time zero.');
    end

    if ~isscalar(baseline_bin_size) || ...
            ~isfinite(baseline_bin_size) || baseline_bin_size <= 0
        error('baseline_bin_size must be a positive finite scalar.');
    end

    if isempty(baseline_index_filename)
        error('baseline_index_filename cannot be empty.');
    end
end

if ~any(strcmp(lower(char(figure_visibility)), {'on', 'off'}))
    error('figure_visibility must be ''on'' or ''off''.');
end

%% ----------------------- Build shared session paths -----------------------
run_g = sprintf('%s_g%d', runName, runind);
destDir = fullfile(root_folder, run_g);
cat_folder = fullfile(destDir, ['catgt_' run_g]);

if isempty(dat_file)
    dat_file = fullfile(cat_folder, 'model_data_allruns.mat');
end

fprintf('destDir        : %s\n', destDir);
fprintf('cat_folder     : %s\n', cat_folder);
fprintf('GPL stim tag   : %s\n', gpl_stim_tag);
fprintf('Target stim tag: %s\n', target_stim_tag);
fprintf('Model data file: %s\n', dat_file);
fprintf('Error method   : %s\n', error_method);
fprintf('Spon subtract  : %d\n', subtract_spontaneous);

if subtract_spontaneous
    fprintf('Baseline window: [%g %g] s\n', ...
        baseline_analysis_window(1), baseline_analysis_window(2));
    fprintf('Baseline bin   : %g s\n', baseline_bin_size);
end

%% ----------------------- Load target-selected unit IDs -----------------------
if ~isfile(dat_file)
    error('Missing model data file: %s', dat_file);
end

M = load(dat_file, 'model_data_allruns');
if ~isfield(M, 'model_data_allruns')
    error('model_data_allruns not found in %s', dat_file);
end

model_data_allruns = M.model_data_allruns;
target_run_idx = find_cell_run_by_stim_tag( ...
    model_data_allruns, target_stim_tag, 'model_data_allruns');
target_run = model_data_allruns{target_run_idx};

%% ----------------------- Process each probe folder -----------------------
for ip = 1:numel(probes)

    thisProbe = probes(ip);
    group_name = char(probe_group_names{ip});
    imecStr = sprintf('imec%d', thisProbe);
    probe_folder = fullfile(cat_folder, [run_g '_' imecStr]);

    fprintf('\n============================================================\n');
    fprintf('Processing probe %d -> %s\n', thisProbe, group_name);
    fprintf('probe_folder: %s\n', probe_folder);
    fprintf('============================================================\n');

    if ~isfolder(probe_folder)
        warning('probe_folder does not exist, skipping probe %d: %s', ...
            thisProbe, probe_folder);
        continue;
    end

    [target_unit_ids, target_unit_field] = ...
        get_selected_unit_ids_for_probe(target_run, thisProbe);

    fprintf('Target unit field: %s\n', target_unit_field);
    fprintf('Target-selected units: %d\n', numel(target_unit_ids));

    %% ----------------------- Find kilosort folders -----------------------
    d = dir(fullfile(probe_folder, 'kilosort*'));
    d = d([d.isdir]);

    if isempty(d)
        warning('No kilosort* folders found under probe %d: %s', ...
            thisProbe, probe_folder);
        continue;
    end

    [~, idx] = sort(lower({d.name}));
    d = d(idx);

    fprintf('Found %d kilosort folder(s) under probe %d.\n', ...
        numel(d), thisProbe);

    %% ----------------------- Process each kilosort folder -----------------------
    for i = 1:numel(d)

        ksDir = fullfile(d(i).folder, d(i).name);
        fprintf('\nProcessing probe %d, ksDir: %s\n', thisProbe, ksDir);

        try
            bined_file = fullfile(ksDir, 'bined_data_allruns.mat');
            if ~isfile(bined_file)
                error('Missing file: %s', bined_file);
            end

            S = load(bined_file, 'bined_data_allruns');
            if ~isfield(S, 'bined_data_allruns')
                error('bined_data_allruns not found in %s', bined_file);
            end

            bined_data_allruns = S.bined_data_allruns;
            gpl_run_idx = find_cell_run_by_stim_tag( ...
                bined_data_allruns, gpl_stim_tag, 'bined_data_allruns');
            gpl_data = bined_data_allruns{gpl_run_idx};

            %% ----------------------- Load requested baseline version -----------------------
            if subtract_spontaneous
                [baseline_data, baseline_source] = ...
                    load_matching_baseline_data( ...
                    ksDir, ...
                    baseline_index_filename, ...
                    gpl_stim_tag, ...
                    baseline_analysis_window, ...
                    baseline_bin_size);
            else
                baseline_data = [];
                baseline_source = struct();
            end

            %% ----------------------- Compute all GPL results -----------------------
            gpl_results = compute_gpl_results( ...
                gpl_data, ...
                target_unit_ids, ...
                target_unit_field, ...
                thisProbe, ...
                group_name, ...
                ksDir, ...
                dat_file, ...
                target_run_idx, ...
                gpl_run_idx, ...
                gpl_stim_tag, ...
                target_stim_tag, ...
                error_method, ...
                pattern_z_threshold, ...
                distribution_bar_width, ...
                subtract_spontaneous, ...
                baseline_data, ...
                baseline_source);

            %% ----------------------- Build output filenames -----------------------
            if subtract_spontaneous
                output_suffix = '_sponsub';
            else
                output_suffix = '';
            end

            output_mat_name = append_suffix_before_extension( ...
                gpl_mat_name, output_suffix);

            % Tuning figures always show raw firing-rate tuning, so their
            % filenames are shared by raw and spontaneous-subtracted modes.
            base_grating = sprintf('%s_Grating_tuning_%s', ...
                group_name, error_method);
            base_plaid = sprintf('%s_Plaid_tuning_%s', ...
                group_name, error_method);

            % Selectivity and pattern/component results change with baseline
            % subtraction, so these names receive _sponsub when enabled.
            base_selectivity = sprintf('%s_DI_OSI_distribution%s', ...
                group_name, output_suffix);
            base_pattern = sprintf('%s_Pattern_component_analysis%s', ...
                group_name, output_suffix);

            gpl_results.output_suffix = output_suffix;
            gpl_results.files = struct();
            gpl_results.files.mat = fullfile(ksDir, output_mat_name);

            gpl_results.files.grating_png = fullfile(ksDir, ...
                [base_grating '.png']);
            gpl_results.files.grating_fig = fullfile(ksDir, ...
                [base_grating '.fig']);

            gpl_results.files.plaid_png = fullfile(ksDir, ...
                [base_plaid '.png']);
            gpl_results.files.plaid_fig = fullfile(ksDir, ...
                [base_plaid '.fig']);

            gpl_results.files.selectivity_png = fullfile(ksDir, ...
                [base_selectivity '.png']);
            gpl_results.files.selectivity_fig = fullfile(ksDir, ...
                [base_selectivity '.fig']);

            gpl_results.files.pattern_component_png = fullfile(ksDir, ...
                [base_pattern '.png']);
            gpl_results.files.pattern_component_fig = fullfile(ksDir, ...
                [base_pattern '.fig']);

            % Save once before plotting so numerical results remain available
            % even if a later figure export fails.
            save(gpl_results.files.mat, 'gpl_results');

            %% ----------------------- Plot grating tuning -----------------------
            fig = plot_direction_tuning( ...
                gpl_results.grating.directions, ...
                gpl_results.grating.mean_fr, ...
                select_error_matrix(gpl_results.grating, error_method), ...
                gpl_results.plot.depth_desc_order, ...
                sprintf('%s grating tuning (%s)', group_name, error_method), ...
                figure_visibility);

            save_figure_pair(fig, ...
                gpl_results.files.grating_png, ...
                gpl_results.files.grating_fig);
            close(fig);

            %% ----------------------- Plot plaid tuning -----------------------
            fig = plot_direction_tuning( ...
                gpl_results.plaid.directions, ...
                gpl_results.plaid.mean_fr, ...
                select_error_matrix(gpl_results.plaid, error_method), ...
                gpl_results.plot.depth_desc_order, ...
                sprintf('%s plaid tuning (%s)', group_name, error_method), ...
                figure_visibility);

            save_figure_pair(fig, ...
                gpl_results.files.plaid_png, ...
                gpl_results.files.plaid_fig);
            close(fig);

            %% ----------------------- Plot DI and OSI distributions -----------------------
            fig = plot_di_osi_distributions( ...
                gpl_results.grating.DI, ...
                gpl_results.grating.OSI, ...
                selectivity_bin_edges, ...
                distribution_bar_width, ...
                group_name, ...
                figure_visibility, ...
                subtract_spontaneous);

            save_figure_pair(fig, ...
                gpl_results.files.selectivity_png, ...
                gpl_results.files.selectivity_fig);
            close(fig);

            %% ----------------------- Plot pattern/component analysis -----------------------
            fig = plot_pattern_component_analysis( ...
                gpl_results.plaid.Zp, ...
                gpl_results.plaid.Zc, ...
                gpl_results.plaid.PI, ...
                gpl_results.plaid.is_pattern, ...
                gpl_results.plaid.is_component, ...
                gpl_results.plaid.is_unclassified, ...
                pattern_z_threshold, ...
                distribution_bar_width, ...
                group_name, ...
                figure_visibility);

            save_figure_pair(fig, ...
                gpl_results.files.pattern_component_png, ...
                gpl_results.files.pattern_component_fig);
            close(fig);

            % Save again after all expected output paths have been recorded.
            save(gpl_results.files.mat, 'gpl_results');

            fprintf('Saved GPL analysis outputs:\n');
            fprintf('  %s\n', gpl_results.files.mat);
            fprintf('  %s\n', gpl_results.files.grating_png);
            fprintf('  %s\n', gpl_results.files.grating_fig);
            fprintf('  %s\n', gpl_results.files.plaid_png);
            fprintf('  %s\n', gpl_results.files.plaid_fig);
            fprintf('  %s\n', gpl_results.files.selectivity_png);
            fprintf('  %s\n', gpl_results.files.selectivity_fig);
            fprintf('  %s\n', gpl_results.files.pattern_component_png);
            fprintf('  %s\n', gpl_results.files.pattern_component_fig);

        catch ME
            fprintf(2, 'Error in probe %d, ksDir %s\n', thisProbe, ksDir);
            fprintf(2, '%s\n', getReport(ME, 'extended', ...
                'hyperlinks', 'off'));
        end
    end
end

fprintf('\nDone.\n');

%% ======================= Local functions =======================

function run_idx = find_cell_run_by_stim_tag(cell_of_structs, stim_tag, source_name)
% Find exactly one cell entry whose .stim_tag matches stim_tag.

    if ~iscell(cell_of_structs)
        error('%s must be a cell array.', source_name);
    end

    run_idx = [];
    for r = 1:numel(cell_of_structs)
        entry = cell_of_structs{r};
        if isempty(entry) || ~isstruct(entry) || ~isfield(entry, 'stim_tag')
            continue;
        end
        if strcmp(entry.stim_tag, stim_tag)
            run_idx(end+1) = r; %#ok<AGROW>
        end
    end

    if isempty(run_idx)
        error('stim_tag not found in %s: %s', source_name, stim_tag);
    end
    if numel(run_idx) > 1
        error('Duplicate stim_tag found in %s: %s', source_name, stim_tag);
    end
end

function [selected_unit_ids, field_name] = ...
        get_selected_unit_ids_for_probe(target_run, probe_id)
% Read the neuron whitelist for one probe from the target model run.
%
% Normal model_data_prepar outputs use:
%   probe0_usedunit_ids, probe1_usedunit_ids, ...
%
% For nan_trial_strategy == 6, unit IDs are data-field specific. Since this
% analysis uses raw firing rate, the fallback field is:
%   raw_fr_probe0_usedunit_ids, raw_fr_probe1_usedunit_ids, ...

    standard_field = sprintf('probe%d_usedunit_ids', probe_id);
    rawfr_field = sprintf('raw_fr_probe%d_usedunit_ids', probe_id);

    if isfield(target_run, standard_field)
        field_name = standard_field;
    elseif isfield(target_run, rawfr_field)
        field_name = rawfr_field;
    else
        error(['Target run contains neither %s nor %s. ' ...
            'Cannot determine selected units for probe %d.'], ...
            standard_field, rawfr_field, probe_id);
    end

    selected_unit_ids = double(target_run.(field_name));
    selected_unit_ids = selected_unit_ids(:);

    if isempty(selected_unit_ids)
        error('Target-selected unit list is empty in field %s.', field_name);
    end
    if any(~isfinite(selected_unit_ids))
        error('Target-selected unit IDs in %s contain nonfinite values.', ...
            field_name);
    end
    if numel(unique(selected_unit_ids)) ~= numel(selected_unit_ids)
        error('Target-selected unit IDs in %s contain duplicates.', field_name);
    end
end

function gpl_results = compute_gpl_results( ...
        gpl_data, target_unit_ids, target_unit_field, ...
        probe_id, group_name, ksDir, dat_file, ...
        target_run_idx, gpl_run_idx, ...
        gpl_stim_tag, target_stim_tag, error_method, ...
        pattern_z_threshold, distribution_bar_width, ...
        subtract_spontaneous, baseline_data, baseline_source)
% Compute tuning, selectivity, and pattern/component metrics for one probe.

    validate_gpl_bined_data(gpl_data);

    gpl_all_unit_ids = double(gpl_data.unit_ids(:));
    if numel(unique(gpl_all_unit_ids)) ~= numel(gpl_all_unit_ids)
        error('gpl_data.unit_ids contains duplicates.');
    end

    % Keep target units in target-list order, and map them to GPL data rows.
    [in_gpl, gpl_unit_index_all] = ismember(target_unit_ids, gpl_all_unit_ids);
    target_unit_index = find(in_gpl);
    gpl_unit_index = gpl_unit_index_all(in_gpl);
    used_unit_ids = target_unit_ids(in_gpl);
    missing_target_unit_ids = target_unit_ids(~in_gpl);

    if isempty(used_unit_ids)
        error(['No target-selected units for probe %d were found in the ' ...
            'GPL run of %s.'], probe_id, ksDir);
    end

    if ~isempty(missing_target_unit_ids)
        warning(['%d/%d target-selected units were not found in the GPL run ' ...
            'for probe %d. Example missing unit ID: %g'], ...
            numel(missing_target_unit_ids), numel(target_unit_ids), ...
            probe_id, missing_target_unit_ids(1));
    end

    raw_fr_all = double(gpl_data.raw_fr);
    nTrial = size(raw_fr_all, 2);

    [selected_raw_fr_bin, response_fr_trial] = ...
        select_and_average_raw_fr( ...
        raw_fr_all, gpl_unit_index, nTrial, 'GPL stimulus raw_fr');

    % Trial-wise stimulus identity reconstructed from conditions.
    [stim_name_trial, grating_dir_trial, plaid_dir_trial, ...
        dir1_trial, dir2_trial] = get_gpl_trial_information( ...
        gpl_data.conditions, ...
        gpl_data.condition_index_per_trial, ...
        nTrial);

    is_grating = strcmpi(stim_name_trial, 'grating');
    is_plaid = strcmpi(stim_name_trial, 'plaid');

    if ~any(is_grating)
        error('No grating trials found in GPL conditions.');
    end
    if ~any(is_plaid)
        error('No plaid trials found in GPL conditions.');
    end

    grating_directions = unique_angles_deg( ...
        grating_dir_trial(is_grating & isfinite(grating_dir_trial)));
    plaid_directions = unique_angles_deg( ...
        plaid_dir_trial(is_plaid & isfinite(plaid_dir_trial)));

    if isempty(grating_directions)
        error('No finite grating directions were found.');
    end
    if isempty(plaid_directions)
        error('No finite plaid directions were found.');
    end

    % Detect one common component separation from plaid dir1/dir2.
    % The calculation is invariant to swapping dir1 and dir2.
    [plaid_component_separation_deg, ...
        plaid_component_separation_per_trial] = ...
        detect_plaid_component_separation( ...
        dir1_trial, dir2_trial, is_plaid);

    % Raw stimulus-period tuning curves.
    [g_mean, g_sd, g_sem, g_n_trials, g_n_valid_trials] = ...
        compute_direction_tuning( ...
        response_fr_trial, grating_dir_trial, is_grating, ...
        grating_directions);

    [p_mean, p_sd, p_sem, p_n_trials, p_n_valid_trials] = ...
        compute_direction_tuning( ...
        response_fr_trial, plaid_dir_trial, is_plaid, ...
        plaid_directions);

    % Optional baseline alignment and baseline-subtracted trial responses.
    if subtract_spontaneous
        baseline = prepare_baseline_for_gpl( ...
            baseline_data, baseline_source, gpl_data, used_unit_ids, ...
            stim_name_trial, grating_dir_trial, plaid_dir_trial, ...
            dir1_trial, dir2_trial, nTrial);

        baseline_fr_trial = baseline.response_fr_trial;
        response_fr_trial_sponsub = response_fr_trial - baseline_fr_trial;

        [g_baseline_mean, g_baseline_sd, g_baseline_sem, ...
            g_baseline_n_trials, g_baseline_n_valid_trials] = ...
            compute_direction_tuning( ...
            baseline_fr_trial, grating_dir_trial, is_grating, ...
            grating_directions);

        [p_baseline_mean, p_baseline_sd, p_baseline_sem, ...
            p_baseline_n_trials, p_baseline_n_valid_trials] = ...
            compute_direction_tuning( ...
            baseline_fr_trial, plaid_dir_trial, is_plaid, ...
            plaid_directions);

        [g_sponsub_mean, g_sponsub_sd, g_sponsub_sem, ...
            g_sponsub_n_trials, g_sponsub_n_valid_trials] = ...
            compute_direction_tuning( ...
            response_fr_trial_sponsub, ...
            grating_dir_trial, is_grating, grating_directions);

        [p_sponsub_mean, p_sponsub_sd, p_sponsub_sem, ...
            p_sponsub_n_trials, p_sponsub_n_valid_trials] = ...
            compute_direction_tuning( ...
            response_fr_trial_sponsub, ...
            plaid_dir_trial, is_plaid, plaid_directions);

        selectivity_mean_fr = g_sponsub_mean;
    else
        baseline = struct();
        baseline_fr_trial = [];
        response_fr_trial_sponsub = [];

        g_baseline_mean = [];
        g_baseline_sd = [];
        g_baseline_sem = [];
        g_baseline_n_trials = [];
        g_baseline_n_valid_trials = [];

        p_baseline_mean = [];
        p_baseline_sd = [];
        p_baseline_sem = [];
        p_baseline_n_trials = [];
        p_baseline_n_valid_trials = [];

        g_sponsub_mean = [];
        g_sponsub_sd = [];
        g_sponsub_sem = [];
        g_sponsub_n_trials = [];
        g_sponsub_n_valid_trials = [];

        p_sponsub_mean = [];
        p_sponsub_sd = [];
        p_sponsub_sem = [];
        p_sponsub_n_trials = [];
        p_sponsub_n_valid_trials = [];

        selectivity_mean_fr = g_mean;
    end

    % DI preferred direction and OSI are based on selectivity_mean_fr.
    % In spontaneous-subtracted mode this is G(theta)-B(theta).
    grating_selectivity = compute_grating_selectivity( ...
        selectivity_mean_fr, grating_directions, subtract_spontaneous);

    % D and P remain raw. Only C uses baseline when requested.
    plaid_analysis = compute_pattern_component_analysis( ...
        g_mean, grating_directions, ...
        p_mean, plaid_directions, ...
        plaid_component_separation_deg, ...
        pattern_z_threshold, ...
        subtract_spontaneous, ...
        g_baseline_mean);

    nUsed = numel(used_unit_ids);
    unit_depth_um_all = get_optional_unit_vector( ...
        gpl_data, 'unit_depth_um', numel(gpl_all_unit_ids));
    unit_channel_all = get_optional_unit_vector( ...
        gpl_data, 'unit_channel', numel(gpl_all_unit_ids));

    unit_depth_um = unit_depth_um_all(gpl_unit_index);
    unit_channel = unit_channel_all(gpl_unit_index);
    depth_desc_order = get_depth_desc_order(unit_depth_um, nUsed);

    % ----------------------- Pack output structure -----------------------
    gpl_results = struct();

    gpl_results.probe_id = probe_id;
    gpl_results.group_name = group_name;
    gpl_results.ksDir = ksDir;

    gpl_results.gpl_stim_tag = gpl_stim_tag;
    gpl_results.target_stim_tag = target_stim_tag;
    gpl_results.model_data_file = dat_file;
    gpl_results.target_run_index = target_run_idx;
    gpl_results.gpl_run_index = gpl_run_idx;
    gpl_results.target_unit_id_field = target_unit_field;

    gpl_results.analysis_window = double(gpl_data.analysis_window);
    gpl_results.window_duration_s = ...
        double(gpl_data.analysis_window(2) - gpl_data.analysis_window(1));
    gpl_results.bin_size = double(gpl_data.bin_size);
    gpl_results.bin_edges = double(gpl_data.bin_edges);
    gpl_results.bin_centers = double(gpl_data.bin_centers);

    gpl_results.error_method_for_plot = error_method;
    gpl_results.plaid_component_separation_deg = ...
        plaid_component_separation_deg;
    gpl_results.pattern_z_threshold = pattern_z_threshold;
    gpl_results.baseline_subtracted = logical(subtract_spontaneous);

    % Unit identity and cross-run alignment.
    gpl_results.target_unit_ids = target_unit_ids(:);
    gpl_results.gpl_all_unit_ids = gpl_all_unit_ids(:);
    gpl_results.used_unit_ids = used_unit_ids(:);
    gpl_results.missing_target_unit_ids = missing_target_unit_ids(:);
    gpl_results.used_neuron_index = (1:nUsed)';
    gpl_results.target_unit_index = target_unit_index(:);
    gpl_results.gpl_unit_index = gpl_unit_index(:);
    gpl_results.unit_depth_um = unit_depth_um(:);
    gpl_results.unit_channel = unit_channel(:);

    % Trial-level stimulus-period data after target-unit filtering.
    gpl_results.raw_fr_bin = selected_raw_fr_bin;
    gpl_results.response_fr_trial = response_fr_trial;
    gpl_results.response_fr_trial_sponsub = response_fr_trial_sponsub;
    gpl_results.condition_index_per_trial = ...
        double(gpl_data.condition_index_per_trial(:));
    gpl_results.conditions = gpl_data.conditions;
    gpl_results.stim_name_trial = stim_name_trial;
    gpl_results.grating_dir_trial = grating_dir_trial;
    gpl_results.plaid_dir_trial = plaid_dir_trial;
    gpl_results.dir1_trial = dir1_trial;
    gpl_results.dir2_trial = dir2_trial;
    gpl_results.plaid_component_separation_per_trial = ...
        plaid_component_separation_per_trial;
    gpl_results.is_grating_trial = is_grating;
    gpl_results.is_plaid_trial = is_plaid;

    % Baseline source, alignment, and trial-level data.
    gpl_results.baseline = baseline;
    if subtract_spontaneous
        gpl_results.baseline.response_fr_trial = baseline_fr_trial;
    end

    % Grating results.
    gpl_results.grating = struct();
    gpl_results.grating.directions = grating_directions;
    gpl_results.grating.mean_fr = g_mean;
    gpl_results.grating.sd_fr = g_sd;
    gpl_results.grating.sem_fr = g_sem;
    gpl_results.grating.n_trials = g_n_trials;
    gpl_results.grating.n_valid_trials = g_n_valid_trials;

    gpl_results.grating.baseline_mean_fr = g_baseline_mean;
    gpl_results.grating.baseline_sd_fr = g_baseline_sd;
    gpl_results.grating.baseline_sem_fr = g_baseline_sem;
    gpl_results.grating.baseline_n_trials = g_baseline_n_trials;
    gpl_results.grating.baseline_n_valid_trials = ...
        g_baseline_n_valid_trials;

    gpl_results.grating.sponsub_mean_fr = g_sponsub_mean;
    gpl_results.grating.sponsub_sd_fr = g_sponsub_sd;
    gpl_results.grating.sponsub_sem_fr = g_sponsub_sem;
    gpl_results.grating.sponsub_n_trials = g_sponsub_n_trials;
    gpl_results.grating.sponsub_n_valid_trials = ...
        g_sponsub_n_valid_trials;
    gpl_results.grating.selectivity_input_mean_fr = selectivity_mean_fr;
    gpl_results.grating.selectivity_input_is_sponsub = ...
        logical(subtract_spontaneous);

    selectivity_fields = fieldnames(grating_selectivity);
    for k = 1:numel(selectivity_fields)
        f = selectivity_fields{k};
        gpl_results.grating.(f) = grating_selectivity.(f);
    end

    % Plaid results.
    gpl_results.plaid = struct();
    gpl_results.plaid.directions = plaid_directions;
    gpl_results.plaid.mean_fr = p_mean;
    gpl_results.plaid.sd_fr = p_sd;
    gpl_results.plaid.sem_fr = p_sem;
    gpl_results.plaid.n_trials = p_n_trials;
    gpl_results.plaid.n_valid_trials = p_n_valid_trials;

    gpl_results.plaid.baseline_mean_fr = p_baseline_mean;
    gpl_results.plaid.baseline_sd_fr = p_baseline_sd;
    gpl_results.plaid.baseline_sem_fr = p_baseline_sem;
    gpl_results.plaid.baseline_n_trials = p_baseline_n_trials;
    gpl_results.plaid.baseline_n_valid_trials = ...
        p_baseline_n_valid_trials;

    gpl_results.plaid.sponsub_mean_fr = p_sponsub_mean;
    gpl_results.plaid.sponsub_sd_fr = p_sponsub_sd;
    gpl_results.plaid.sponsub_sem_fr = p_sponsub_sem;
    gpl_results.plaid.sponsub_n_trials = p_sponsub_n_trials;
    gpl_results.plaid.sponsub_n_valid_trials = ...
        p_sponsub_n_valid_trials;

    plaid_fields = fieldnames(plaid_analysis);
    for k = 1:numel(plaid_fields)
        f = plaid_fields{k};
        gpl_results.plaid.(f) = plaid_analysis.(f);
    end

    gpl_results.plot = struct();
    gpl_results.plot.depth_desc_order = depth_desc_order(:);
    gpl_results.plot.distribution_bar_width = distribution_bar_width;
    gpl_results.plot.sponsub_DI_regular_max = 1.6;
end

function [baseline_data, source] = load_matching_baseline_data( ...
        ksDir, index_filename, stim_tag, requested_window, requested_bin_size)
% Select one baseline file by stim_tag, analysis_window, and bin_size.

    index_file = fullfile(ksDir, index_filename);
    if ~isfile(index_file)
        error('Missing baseline index file: %s', index_file);
    end

    S = load(index_file, 'baseline_data_index');
    if ~isfield(S, 'baseline_data_index') || ...
            ~isstruct(S.baseline_data_index)
        error(['baseline_data_index is missing or invalid in baseline ' ...
            'index file: %s'], index_file);
    end

    baseline_data_index = S.baseline_data_index;
    required_fields = {'stim_tag', 'analysis_window', 'bin_size', 'filename'};

    for f = 1:numel(required_fields)
        if ~isfield(baseline_data_index, required_fields{f})
            error('baseline_data_index is missing field %s in %s.', ...
                required_fields{f}, index_file);
        end
    end

    match = false(1, numel(baseline_data_index));

    for k = 1:numel(baseline_data_index)
        same_tag = strcmp( ...
            char(string(baseline_data_index(k).stim_tag)), stim_tag);

        same_window = same_numeric_vector_local( ...
            baseline_data_index(k).analysis_window, requested_window);

        same_bin = numeric_scalar_equal_local( ...
            baseline_data_index(k).bin_size, requested_bin_size);

        match(k) = same_tag && same_window && same_bin;
    end

    if ~any(match)
        error(['No baseline index entry matches the requested GPL baseline.\n' ...
            'stim_tag: %s\nanalysis_window: [%g %g]\nbin_size: %g\n' ...
            'Index file: %s'], ...
            stim_tag, requested_window(1), requested_window(2), ...
            requested_bin_size, index_file);
    end

    if sum(match) > 1
        error(['Multiple baseline index entries match stim_tag %s, window ' ...
            '[%g %g], and bin_size %g in %s.'], ...
            stim_tag, requested_window(1), requested_window(2), ...
            requested_bin_size, index_file);
    end

    entry_index = find(match, 1);
    entry = baseline_data_index(entry_index);
    stored_filename = char(string(entry.filename));

    if isfile(stored_filename)
        baseline_file = stored_filename;
    else
        baseline_file = fullfile(ksDir, stored_filename);
    end

    if ~isfile(baseline_file)
        error('Baseline data file referenced by the index is missing: %s', ...
            baseline_file);
    end

    B = load(baseline_file, 'baseline_data');
    if ~isfield(B, 'baseline_data') || ~isstruct(B.baseline_data)
        error('baseline_data is missing or invalid in %s.', baseline_file);
    end

    baseline_data = B.baseline_data;

    validate_loaded_baseline_identity( ...
        baseline_data, stim_tag, requested_window, requested_bin_size, ...
        baseline_file);

    source = struct();
    source.index_file = index_file;
    source.data_file = baseline_file;
    source.index_entry_number = entry_index;
    source.index_entry = entry;
    source.requested_stim_tag = stim_tag;
    source.requested_analysis_window = double(requested_window(:)');
    source.requested_bin_size = double(requested_bin_size);
end

function validate_loaded_baseline_identity( ...
        baseline_data, stim_tag, requested_window, requested_bin_size, ...
        baseline_file)
% Recheck the baseline file itself rather than trusting only the index.

    required_fields = {'stim_tag', 'analysis_window', 'bin_size'};

    for k = 1:numel(required_fields)
        if ~isfield(baseline_data, required_fields{k})
            error('Baseline file %s is missing field %s.', ...
                baseline_file, required_fields{k});
        end
    end

    if ~strcmp(char(string(baseline_data.stim_tag)), stim_tag)
        error(['Baseline stim_tag in file does not match the request.\n' ...
            'File: %s\nExpected: %s\nFound: %s'], ...
            baseline_file, stim_tag, char(string(baseline_data.stim_tag)));
    end

    if ~same_numeric_vector_local( ...
            baseline_data.analysis_window, requested_window)
        error(['Baseline analysis_window in file does not match the request. ' ...
            'File: %s'], baseline_file);
    end

    if ~numeric_scalar_equal_local( ...
            baseline_data.bin_size, requested_bin_size)
        error(['Baseline bin_size in file does not match the request. ' ...
            'File: %s'], baseline_file);
    end
end

function baseline = prepare_baseline_for_gpl( ...
        baseline_data, baseline_source, gpl_data, used_unit_ids, ...
        stim_name_trial, grating_dir_trial, plaid_dir_trial, ...
        dir1_trial, dir2_trial, nTrial)
% Validate and align one baseline file to the selected GPL units and trials.

    required_fields = { ...
        'stim_tag', ...
        'unit_ids', ...
        'analysis_window', ...
        'bin_size', ...
        'condition_index_per_trial', ...
        'conditions', ...
        'raw_fr'};

    for k = 1:numel(required_fields)
        if ~isfield(baseline_data, required_fields{k})
            error('Loaded baseline_data is missing field %s.', ...
                required_fields{k});
        end
    end

    if ~strcmp(char(string(baseline_data.stim_tag)), ...
            char(string(gpl_data.stim_tag)))
        error('Baseline and GPL stimulus tags do not match.');
    end

    baseline_condition_index = ...
        double(baseline_data.condition_index_per_trial(:));

    gpl_condition_index = ...
        double(gpl_data.condition_index_per_trial(:));

    if numel(baseline_condition_index) ~= nTrial
        error(['Baseline trial count from condition_index_per_trial (%d) ' ...
            'does not match GPL trial count (%d).'], ...
            numel(baseline_condition_index), nTrial);
    end

    same_condition_index = same_value_vector_local( ...
        baseline_condition_index, gpl_condition_index);

    if ~all(same_condition_index)
        mismatch = find(~same_condition_index, 1);
        error(['Baseline and GPL condition_index_per_trial differ. ' ...
            'First mismatch is trial %d.'], mismatch);
    end

    [base_stim_name, base_grating_dir, base_plaid_dir, ...
        base_dir1, base_dir2] = ...
        get_gpl_trial_information( ...
        baseline_data.conditions, ...
        baseline_condition_index, ...
        nTrial);

    assert_same_trial_text( ...
        base_stim_name, stim_name_trial, 'stim_name');
    assert_same_trial_angles( ...
        base_grating_dir, grating_dir_trial, 'grating_dir');
    assert_same_trial_angles( ...
        base_plaid_dir, plaid_dir_trial, 'plaid_dir');
    assert_same_trial_angles( ...
        base_dir1, dir1_trial, 'dir1');
    assert_same_trial_angles( ...
        base_dir2, dir2_trial, 'dir2');

    baseline_all_unit_ids = double(baseline_data.unit_ids(:));

    if numel(unique(baseline_all_unit_ids)) ~= ...
            numel(baseline_all_unit_ids)
        error('baseline_data.unit_ids contains duplicates.');
    end

    [found, baseline_unit_index] = ...
        ismember(used_unit_ids, baseline_all_unit_ids);

    if any(~found)
        missing = used_unit_ids(~found);
        error(['%d selected GPL units are missing from the baseline file. ' ...
            'Example missing unit ID: %g'], ...
            numel(missing), missing(1));
    end

    baseline_raw_fr_all = double(baseline_data.raw_fr);

    if size(baseline_raw_fr_all, 1) ~= numel(baseline_all_unit_ids)
        error(['First dimension of baseline raw_fr does not match the ' ...
            'number of baseline unit IDs.']);
    end

    if size(baseline_raw_fr_all, 2) ~= nTrial
        error(['Second dimension of baseline raw_fr (%d) does not match ' ...
            'GPL trial count (%d).'], ...
            size(baseline_raw_fr_all, 2), nTrial);
    end

    [selected_raw_fr_bin, response_fr_trial] = ...
        select_and_average_raw_fr( ...
        baseline_raw_fr_all, ...
        baseline_unit_index, ...
        nTrial, ...
        'baseline raw_fr');

    baseline = struct();
    baseline.stim_tag = char(string(baseline_data.stim_tag));
    baseline.analysis_window = ...
        double(baseline_data.analysis_window(:)');
    baseline.bin_size = double(baseline_data.bin_size);

    baseline.index_file = baseline_source.index_file;
    baseline.data_file = baseline_source.data_file;
    baseline.index_entry_number = ...
        baseline_source.index_entry_number;
    baseline.index_entry = baseline_source.index_entry;

    baseline.all_unit_ids = baseline_all_unit_ids;
    baseline.unit_index = baseline_unit_index(:);
    baseline.used_unit_ids = used_unit_ids(:);

    baseline.raw_fr_bin = selected_raw_fr_bin;
    baseline.response_fr_trial = response_fr_trial;

    baseline.condition_index_per_trial = baseline_condition_index;
    baseline.conditions = baseline_data.conditions;
end

function [selected_raw_fr_bin, response_fr_trial] = ...
        select_and_average_raw_fr( ...
        raw_fr_all, unit_index, nTrial, label_text)
% Select unit rows and average equal-duration firing-rate bins per trial.

    if ~(ismatrix(raw_fr_all) || ndims(raw_fr_all) == 3)
        error('%s must be unit x trial or unit x trial x bin.', label_text);
    end

    if size(raw_fr_all, 2) ~= nTrial
        error('%s trial dimension does not match nTrial.', label_text);
    end

    if ismatrix(raw_fr_all)
        selected_raw_fr_bin = raw_fr_all(unit_index, :);
        response_fr_trial = selected_raw_fr_bin;
    else
        selected_raw_fr_bin = raw_fr_all(unit_index, :, :);
        response_fr_trial = mean(selected_raw_fr_bin, 3, 'omitnan');
        response_fr_trial = reshape( ...
            response_fr_trial, ...
            [numel(unit_index), nTrial]);
    end
end

function assert_same_trial_text(a, b, field_name)
% Require identical trial-wise text labels, ignoring letter case.

    a = string(a(:));
    b = string(b(:));

    if numel(a) ~= numel(b)
        error('Trial-wise %s vectors have different lengths.', field_name);
    end

    same = strcmpi(strtrim(a), strtrim(b));

    if any(~same)
        idx = find(~same, 1);
        error(['Baseline and GPL trial-wise %s differ at trial %d. ' ...
            'Baseline: %s; GPL: %s'], ...
            field_name, idx, char(a(idx)), char(b(idx)));
    end
end

function assert_same_trial_angles(a, b, field_name)
% Require identical trial-wise angles, treating paired NaNs as equal.

    a = double(a(:));
    b = double(b(:));

    if numel(a) ~= numel(b)
        error('Trial-wise %s vectors have different lengths.', field_name);
    end

    both_nan = isnan(a) & isnan(b);
    both_finite_equal = angle_equal(a, b);
    same = both_nan | both_finite_equal;

    if any(~same)
        idx = find(~same, 1);
        error(['Baseline and GPL trial-wise %s differ at trial %d. ' ...
            'Baseline: %.12g; GPL: %.12g'], ...
            field_name, idx, a(idx), b(idx));
    end
end

function tf = same_value_vector_local(a, b)
% Element-wise equality helper used for mismatch localization.

    a = double(a(:));
    b = double(b(:));

    if numel(a) ~= numel(b)
        tf = false(max(numel(a), numel(b)), 1);
        return;
    end

    tf = (a == b) | (isnan(a) & isnan(b));
end

function tf = same_numeric_vector_local(a, b)
% Compare two finite numeric vectors with a small tolerance.

    a = double(a(:));
    b = double(b(:));

    if numel(a) ~= numel(b)
        tf = false;
        return;
    end

    if isempty(a)
        tf = true;
        return;
    end

    if any(~isfinite(a)) || any(~isfinite(b))
        tf = isequaln(a, b);
        return;
    end

    tf = all(abs(a - b) <= 1e-12);
end

function tf = numeric_scalar_equal_local(a, b)
% Compare two finite numeric scalars with a small tolerance.

    tf = isnumeric(a) && isnumeric(b) && ...
        isscalar(a) && isscalar(b) && ...
        isfinite(a) && isfinite(b) && ...
        abs(double(a) - double(b)) <= 1e-12;
end

function filename_out = ...
        append_suffix_before_extension(filename_in, suffix)
% Add a suffix immediately before a filename extension.

    [folder_part, base_name, extension] = ...
        fileparts(char(filename_in));

    filename_only = [base_name char(suffix) extension];

    if isempty(folder_part)
        filename_out = filename_only;
    else
        filename_out = fullfile(folder_part, filename_only);
    end
end

function validate_gpl_bined_data(gpl_data)
% Validate the fields required from one bined_data_allruns GPL entry.

    required_fields = { ...
        'stim_tag', ...
        'unit_ids', ...
        'unit_depth_um', ...
        'unit_channel', ...
        'analysis_window', ...
        'bin_size', ...
        'bin_edges', ...
        'bin_centers', ...
        'condition_index_per_trial', ...
        'conditions', ...
        'raw_fr'};

    for k = 1:numel(required_fields)
        f = required_fields{k};
        if ~isfield(gpl_data, f)
            error('Required field missing from GPL bined data: %s', f);
        end
    end

    if ~(ismatrix(gpl_data.raw_fr) || ndims(gpl_data.raw_fr) == 3)
        error('gpl_data.raw_fr must be unit x trial or unit x trial x bin.');
    end

    if ~isstruct(gpl_data.conditions)
        error('gpl_data.conditions must be a struct array.');
    end

    if ~isnumeric(gpl_data.condition_index_per_trial)
        error('gpl_data.condition_index_per_trial must be numeric.');
    end

    if numel(gpl_data.analysis_window) ~= 2 || ...
            gpl_data.analysis_window(2) <= gpl_data.analysis_window(1)
        error('gpl_data.analysis_window must be [start end] with end > start.');
    end
end

function v = get_optional_unit_vector(S, field_name, nUnit)
% Read an optional unit-level vector; return NaN if invalid or absent.

    if isfield(S, field_name)
        v = double(S.(field_name));
        v = v(:);
        if numel(v) ~= nUnit
            warning('%s length does not match nUnit. Filling with NaN.', ...
                field_name);
            v = nan(nUnit, 1);
        end
    else
        warning('%s is absent. Filling with NaN.', field_name);
        v = nan(nUnit, 1);
    end
end

function [stim_name_trial, grating_dir_trial, plaid_dir_trial, ...
        dir1_trial, dir2_trial] = get_gpl_trial_information( ...
        conditions, condition_index_per_trial, nTrial)
% Reconstruct trial-wise GPL labels from the condition structure.
%
% Required condition fields:
%   stim_name, grating_dir, plaid_dir, dir1, dir2
%
% dir1 and dir2 are used only to detect the common plaid component
% separation. Their order is not important.

    condition_index_per_trial = double(condition_index_per_trial(:));
    if numel(condition_index_per_trial) ~= nTrial
        error(['condition_index_per_trial length (%d) does not match ' ...
            'the GPL trial count (%d).'], ...
            numel(condition_index_per_trial), nTrial);
    end

    stim_field = find_field_case_insensitive(conditions, {'stim_name'});
    gdir_field = find_field_case_insensitive(conditions, {'grating_dir'});
    pdir_field = find_field_case_insensitive(conditions, {'plaid_dir'});
    dir1_field = find_field_case_insensitive(conditions, {'dir1'});
    dir2_field = find_field_case_insensitive(conditions, {'dir2'});

    if isempty(stim_field)
        error('No stim_name field found in GPL conditions.');
    end
    if isempty(gdir_field)
        error('No grating_dir field found in GPL conditions.');
    end
    if isempty(pdir_field)
        error('No plaid_dir field found in GPL conditions.');
    end
    if isempty(dir1_field)
        error('No dir1 field found in GPL conditions.');
    end
    if isempty(dir2_field)
        error('No dir2 field found in GPL conditions.');
    end

    stim_name_trial = strings(nTrial, 1);
    grating_dir_trial = nan(nTrial, 1);
    plaid_dir_trial = nan(nTrial, 1);
    dir1_trial = nan(nTrial, 1);
    dir2_trial = nan(nTrial, 1);

    nCond = numel(conditions);
    for t = 1:nTrial
        c = condition_index_per_trial(t);
        if ~isfinite(c) || c ~= round(c) || c < 1 || c > nCond
            error('Invalid condition index at trial %d: %g', t, c);
        end

        stim_name_trial(t) = get_text_scalar_from_struct( ...
            conditions(c), stim_field);
        grating_dir_trial(t) = get_numeric_scalar_from_struct( ...
            conditions(c), gdir_field);
        plaid_dir_trial(t) = get_numeric_scalar_from_struct( ...
            conditions(c), pdir_field);
        dir1_trial(t) = get_numeric_scalar_from_struct( ...
            conditions(c), dir1_field);
        dir2_trial(t) = get_numeric_scalar_from_struct( ...
            conditions(c), dir2_field);
    end

    % Preserve non-integer direction designs. Only canonicalize to [0, 360).
    grating_dir_trial = canonical_angle_deg(grating_dir_trial);
    plaid_dir_trial = canonical_angle_deg(plaid_dir_trial);
    dir1_trial = canonical_angle_deg(dir1_trial);
    dir2_trial = canonical_angle_deg(dir2_trial);
end

function field_name = find_field_case_insensitive(S, candidates)
% Case-insensitive field-name lookup.

    field_name = '';
    if isempty(S) || ~isstruct(S)
        return;
    end

    fn = fieldnames(S);
    fn_lower = lower(fn);

    for k = 1:numel(candidates)
        idx = find(strcmp(fn_lower, lower(candidates{k})), 1);
        if ~isempty(idx)
            field_name = fn{idx};
            return;
        end
    end
end

function x = get_numeric_scalar_from_struct(S, field_name)
% Extract one numeric scalar, allowing NaN values.

    val = S.(field_name);

    if isnumeric(val) || islogical(val)
        if isempty(val) || numel(val) ~= 1
            error('Field %s is not a scalar.', field_name);
        end
        x = double(val);
        return;
    end

    if isstring(val)
        if isempty(val) || numel(val) ~= 1
            error('Field %s is not a scalar string.', field_name);
        end
        x = str2double(val);
        return;
    end

    if ischar(val)
        x = str2double(val);
        return;
    end

    error('Field %s has unsupported type.', field_name);
end

function txt = get_text_scalar_from_struct(S, field_name)
% Extract one text scalar and return it as a string scalar.

    val = S.(field_name);

    if isstring(val)
        if numel(val) ~= 1
            error('Field %s is not a scalar string.', field_name);
        end
        txt = val;
    elseif ischar(val)
        txt = string(val);
    elseif iscell(val)
        if numel(val) ~= 1
            error('Field %s is not a scalar cell.', field_name);
        end
        txt = string(val{1});
    elseif iscategorical(val)
        if numel(val) ~= 1
            error('Field %s is not a scalar categorical value.', field_name);
        end
        txt = string(val);
    else
        error('Field %s has unsupported text type.', field_name);
    end

    txt = strtrim(txt);
end

function a = canonical_angle_deg(a, doRound)
% Convert finite angles in degrees to canonical [0, 360).

    if nargin < 2
        doRound = false;
    end

    if doRound
        finite_mask = isfinite(a);
        a(finite_mask) = round(a(finite_mask));
    end

    a = mod(a, 360);
    tol = 1e-10;
    a(abs(a) < tol) = 0;
    a(abs(a - 360) < tol) = 0;
end

function directions = unique_angles_deg(a)
% Return sorted unique canonical directions without assuming fixed spacing.

    a = canonical_angle_deg(double(a(:)));
    a = sort(a(isfinite(a)));

    if isempty(a)
        directions = zeros(1, 0);
        return;
    end

    directions = a(1);
    for k = 2:numel(a)
        if ~angle_equal(a(k), directions(end))
            directions(end+1) = a(k); %#ok<AGROW>
        end
    end

    directions = directions(:)';
end

function [separation_deg, separation_per_trial] = ...
        detect_plaid_component_separation(dir1_trial, dir2_trial, is_plaid)
% Detect one common minimum circular separation between plaid components.
%
% This is invariant to whether the two component directions are stored as
% (dir1, dir2) or (dir2, dir1).

    dir1_trial = double(dir1_trial(:));
    dir2_trial = double(dir2_trial(:));
    is_plaid = logical(is_plaid(:));

    if numel(dir1_trial) ~= numel(dir2_trial) || ...
            numel(dir1_trial) ~= numel(is_plaid)
        error('dir1, dir2, and is_plaid must have the same length.');
    end

    if any(is_plaid & (~isfinite(dir1_trial) | ~isfinite(dir2_trial)))
        bad_trial = find(is_plaid & ...
            (~isfinite(dir1_trial) | ~isfinite(dir2_trial)), 1);
        error(['Cannot detect plaid component separation because plaid ' ...
            'trial %d has nonfinite dir1 or dir2.'], bad_trial);
    end

    separation_per_trial = nan(numel(is_plaid), 1);
    plaid_delta = mod(dir1_trial(is_plaid) - dir2_trial(is_plaid) + 180, ...
        360) - 180;
    plaid_separation = abs(plaid_delta);
    separation_per_trial(is_plaid) = plaid_separation;

    if isempty(plaid_separation)
        error('No plaid trials are available for component-separation detection.');
    end

    tol = 1e-8;
    if any(plaid_separation <= tol)
        bad_trial_rel = find(plaid_separation <= tol, 1);
        plaid_trials = find(is_plaid);
        error(['Invalid plaid component separation at trial %d: dir1 and ' ...
            'dir2 are identical.'], plaid_trials(bad_trial_rel));
    end

    if any(plaid_separation >= 180 - tol)
        bad_trial_rel = find(plaid_separation >= 180 - tol, 1);
        plaid_trials = find(is_plaid);
        error(['Invalid or ambiguous plaid component separation at trial %d: ' ...
            'the minimum circular separation is %.12g degrees.'], ...
            plaid_trials(bad_trial_rel), plaid_separation(bad_trial_rel));
    end

    separation_deg = median(plaid_separation);
    inconsistent = abs(plaid_separation - separation_deg) > tol;
    if any(inconsistent)
        unique_sep = unique(round(plaid_separation * 1e8) / 1e8);
        error(['Plaid component separation is not constant across trials. ' ...
            'Detected separations (deg): %s'], mat2str(unique_sep(:)'));
    end

    % Use the common measured value rather than a user-entered parameter.
    % The median is robust to tiny floating-point residue.
end

function [mean_fr, sd_fr, sem_fr, n_trials, n_valid_trials] = ...
        compute_direction_tuning( ...
        response_fr_trial, direction_trial, trial_type_mask, directions)
% Compute mean, SD, and SEM for each direction.
%
% response_fr_trial: unit x trial
% direction_trial   : trial x 1
% trial_type_mask   : trial x 1
% directions        : 1 x nDirection

    [nUnit, nTrial] = size(response_fr_trial);
    direction_trial = direction_trial(:);
    trial_type_mask = logical(trial_type_mask(:));

    if numel(direction_trial) ~= nTrial || numel(trial_type_mask) ~= nTrial
        error('Trial metadata length does not match response_fr_trial.');
    end

    nDir = numel(directions);
    mean_fr = nan(nUnit, nDir);
    sd_fr = nan(nUnit, nDir);
    sem_fr = nan(nUnit, nDir);
    n_trials = zeros(1, nDir);
    n_valid_trials = zeros(nUnit, nDir);

    for d = 1:nDir
        take = trial_type_mask & angle_equal(direction_trial, directions(d));
        n_trials(d) = sum(take);

        if n_trials(d) == 0
            continue;
        end

        X = response_fr_trial(:, take);
        mean_fr(:, d) = mean(X, 2, 'omitnan');
        sd_fr(:, d) = std(X, 0, 2, 'omitnan');
        n_valid_trials(:, d) = sum(isfinite(X), 2);

        denom = sqrt(n_valid_trials(:, d));
        sem_fr(:, d) = sd_fr(:, d) ./ denom;
        sem_fr(n_valid_trials(:, d) == 0, d) = NaN;
    end
end

function tf = angle_equal(a, b)
% Circular angle equality with a small numerical tolerance.

    tol = 1e-8;
    delta = mod(a - b + 180, 360) - 180;
    tf = isfinite(a) & isfinite(b) & abs(delta) < tol;
end

function out = compute_grating_selectivity( ...
        mean_fr, directions, use_absolute_osi_denominator)
% Compute DI and vector-based OSI from one grating mean tuning curve.
%
% In raw mode:
%   DI  = 1 - Rnull/Ropt
%   OSI = |sum R(theta) exp(i*2*theta)| / sum R(theta)
%
% In spontaneous-subtracted mode, mean_fr is G(theta)-B(theta):
%   preferred direction is selected from G-B;
%   DI is valid only when Ropt > 0;
%   OSI denominator is sum(abs(G-B)).
%
% OSI is calculated directly from all measured directions. Opposite
% directions are not explicitly averaged or merged before the vector sum.
%
% DI never uses interpolation. Every measured grating direction must have
% an exactly measured direction 180 degrees away; otherwise this function
% raises an error before processing neurons.

    [nUnit, nDir] = size(mean_fr);
    directions = directions(:)';

    if numel(directions) ~= nDir
        error('Direction count does not match grating tuning columns.');
    end

    use_absolute_osi_denominator = ...
        logical(use_absolute_osi_denominator);

    % Precompute and strictly validate the opposite-direction mapping.
    opposite_index = nan(1, nDir);

    for d = 1:nDir
        required_null_dir = mod(directions(d) + 180, 360);
        idx = find(angle_equal(directions, required_null_dir));

        if isempty(idx)
            error(['Cannot calculate DI because grating direction %.12g deg ' ...
                'does not have a measured opposite direction %.12g deg. ' ...
                'Interpolation is disabled.'], ...
                directions(d), required_null_dir);
        end

        if numel(idx) > 1
            error(['Cannot calculate DI because opposite direction %.12g deg ' ...
                'matches multiple grating-direction columns.'], ...
                required_null_dir);
        end

        opposite_index(d) = idx;
    end

    preferred_direction = nan(nUnit, 1);
    null_direction = nan(nUnit, 1);

    preferred_direction_index = nan(nUnit, 1);
    null_direction_index = nan(nUnit, 1);

    Ropt_fr = nan(nUnit, 1);
    Rnull_fr = nan(nUnit, 1);

    DI = nan(nUnit, 1);
    DI_valid = false(nUnit, 1);
    DI_excluded_nonpositive_opt = false(nUnit, 1);

    OSI = nan(nUnit, 1);
    preferred_orientation = nan(nUnit, 1);

    orientation_vector = ...
        complex(nan(nUnit, 1), nan(nUnit, 1));

    OSI_denominator = nan(nUnit, 1);
    OSI_valid = false(nUnit, 1);

    for u = 1:nUnit
        R = mean_fr(u, :);
        finite_mask = isfinite(R) & isfinite(directions);

        %% Direction index
        if any(finite_mask)
            finite_idx = find(finite_mask);

            [rmax, local_idx] = max(R(finite_mask));
            pref_idx = finite_idx(local_idx);  % first maximum if tied
            null_idx = opposite_index(pref_idx);

            preferred_direction(u) = directions(pref_idx);
            null_direction(u) = directions(null_idx);

            preferred_direction_index(u) = pref_idx;
            null_direction_index(u) = null_idx;

            Ropt_fr(u) = rmax;

            if isfinite(R(null_idx))
                Rnull_fr(u) = R(null_idx);
            end

            if isfinite(rmax) && ...
                    rmax > 0 && ...
                    isfinite(R(null_idx))

                DI(u) = 1 - R(null_idx) / rmax;
                DI_valid(u) = true;

            elseif isfinite(rmax) && rmax <= 0

                DI_excluded_nonpositive_opt(u) = true;
            end
        end

        %% Orientation selectivity index
        if any(finite_mask)
            Ruse = R(finite_mask);
            theta = directions(finite_mask);

            if use_absolute_osi_denominator
                denominator = sum(abs(Ruse));
            else
                denominator = sum(Ruse);
            end

            OSI_denominator(u) = denominator;

            if isfinite(denominator) && denominator > 0
                vec = sum( ...
                    Ruse .* exp(1i * 2 * deg2rad(theta)));

                OSI(u) = abs(vec) / denominator;
                orientation_vector(u) = vec;
                OSI_valid(u) = true;

                if abs(vec) > 0
                    preferred_orientation(u) = ...
                        mod(rad2deg(angle(vec)) / 2, 180);
                end
            end
        end
    end

    DI_finite = DI(isfinite(DI));
    OSI_finite = OSI(isfinite(OSI));

    out = struct();

    out.preferred_direction = preferred_direction;
    out.null_direction = null_direction;
    out.preferred_direction_index = preferred_direction_index;
    out.null_direction_index = null_direction_index;
    out.opposite_direction_index = opposite_index;

    out.Ropt_fr = Ropt_fr;
    out.Rnull_fr = Rnull_fr;

    out.DI = DI;
    out.DI_valid = DI_valid;
    out.DI_excluded_nonpositive_opt = ...
        DI_excluded_nonpositive_opt;
    out.DI_mean = mean(DI_finite, 'omitnan');
    out.DI_median = median(DI_finite, 'omitnan');
    out.DI_n_valid = numel(DI_finite);
    out.DI_n_above_1p6 = sum(DI_finite > 1.6);

    out.preferred_orientation = preferred_orientation;
    out.orientation_vector = orientation_vector;
    out.OSI_denominator = OSI_denominator;

    out.OSI_denominator_uses_absolute_response = ...
        use_absolute_osi_denominator;

    out.OSI = OSI;
    out.OSI_valid = OSI_valid;
    out.OSI_mean = mean(OSI_finite, 'omitnan');
    out.OSI_median = median(OSI_finite, 'omitnan');
    out.OSI_n_valid = numel(OSI_finite);
end

function out = compute_pattern_component_analysis( ...
        grating_mean_fr, grating_directions, ...
        plaid_mean_fr, plaid_directions, ...
        component_separation_deg, z_threshold, ...
        subtract_spontaneous, grating_baseline_mean_fr)
% Build pattern/component predictors and compute correlation-based metrics.
%
% The actual plaid response and pattern prediction remain raw:
%   D(theta) = measured raw plaid firing rate
%   P(theta) = raw G(theta)
%
% Without baseline subtraction:
%   C(theta) = G1 + G2
%
% With baseline subtraction:
%   C(theta) = (G1-B1) + (G2-B2) + (B1+B2)/2
%            = G1 + G2 - (B1+B2)/2
%
% No interpolation is used.

    [nUnit, nGdir] = size(grating_mean_fr);
    [nUnitP, nPdir] = size(plaid_mean_fr);

    if nUnitP ~= nUnit
        error('Grating and plaid tuning matrices have different unit counts.');
    end

    grating_directions = grating_directions(:)';
    plaid_directions = plaid_directions(:)';

    if numel(grating_directions) ~= nGdir
        error('grating_directions does not match grating tuning columns.');
    end

    if numel(plaid_directions) ~= nPdir
        error('plaid_directions does not match plaid tuning columns.');
    end

    subtract_spontaneous = logical(subtract_spontaneous);

    if subtract_spontaneous
        if ~isequal(size(grating_baseline_mean_fr), size(grating_mean_fr))
            error(['grating_baseline_mean_fr must match grating_mean_fr ' ...
                'when spontaneous subtraction is enabled.']);
        end
    end

    if ~isscalar(component_separation_deg) || ...
            ~isfinite(component_separation_deg) || ...
            component_separation_deg <= 0 || ...
            component_separation_deg >= 180

        error(['component_separation_deg must be a finite scalar strictly ' ...
            'between 0 and 180 degrees.']);
    end

    half_sep = component_separation_deg / 2;

    component_direction_1 = ...
        mod(plaid_directions - half_sep, 360);

    component_direction_2 = ...
        mod(plaid_directions + half_sep, 360);

    pattern_prediction = nan(nUnit, nPdir);
    component_prediction = nan(nUnit, nPdir);

    component_raw_response_1 = nan(nUnit, nPdir);
    component_raw_response_2 = nan(nUnit, nPdir);

    component_baseline_1 = nan(nUnit, nPdir);
    component_baseline_2 = nan(nUnit, nPdir);

    component_evoked_response_1 = nan(nUnit, nPdir);
    component_evoked_response_2 = nan(nUnit, nPdir);

    component_baseline_addback = nan(nUnit, nPdir);

    pattern_grating_index = nan(1, nPdir);
    component1_grating_index = nan(1, nPdir);
    component2_grating_index = nan(1, nPdir);

    % Strictly require every direction used by P(theta) and C(theta) to
    % have been measured in the grating run. No interpolation is allowed.
    for d = 1:nPdir
        pattern_idx = ...
            require_unique_direction_index( ...
            grating_directions, ...
            plaid_directions(d), ...
            sprintf( ...
            'pattern prediction for plaid direction %.12g deg', ...
            plaid_directions(d)));

        comp1_idx = ...
            require_unique_direction_index( ...
            grating_directions, ...
            component_direction_1(d), ...
            sprintf([ ...
            'component prediction for plaid direction %.12g deg ' ...
            '(theta - alpha/2)'], ...
            plaid_directions(d)));

        comp2_idx = ...
            require_unique_direction_index( ...
            grating_directions, ...
            component_direction_2(d), ...
            sprintf([ ...
            'component prediction for plaid direction %.12g deg ' ...
            '(theta + alpha/2)'], ...
            plaid_directions(d)));

        pattern_grating_index(d) = pattern_idx;
        component1_grating_index(d) = comp1_idx;
        component2_grating_index(d) = comp2_idx;

        pattern_prediction(:, d) = ...
            grating_mean_fr(:, pattern_idx);

        G1 = grating_mean_fr(:, comp1_idx);
        G2 = grating_mean_fr(:, comp2_idx);

        component_raw_response_1(:, d) = G1;
        component_raw_response_2(:, d) = G2;

        if subtract_spontaneous
            B1 = grating_baseline_mean_fr(:, comp1_idx);
            B2 = grating_baseline_mean_fr(:, comp2_idx);

            addback = (B1 + B2) / 2;

            component_baseline_1(:, d) = B1;
            component_baseline_2(:, d) = B2;

            component_evoked_response_1(:, d) = G1 - B1;
            component_evoked_response_2(:, d) = G2 - B2;

            component_baseline_addback(:, d) = addback;

            component_prediction(:, d) = ...
                (G1 - B1) + ...
                (G2 - B2) + ...
                addback;
        else
            component_prediction(:, d) = G1 + G2;
        end
    end

    rp = nan(nUnit, 1);
    rc = nan(nUnit, 1);
    rpc = nan(nUnit, 1);

    Rp = nan(nUnit, 1);
    Rc = nan(nUnit, 1);

    Zp = nan(nUnit, 1);
    Zc = nan(nUnit, 1);
    PI = nan(nUnit, 1);

    n_valid_directions = zeros(nUnit, 1);
    partial_corr_valid = false(nUnit, 1);
    fisher_z_clipped = false(nUnit, 1);

    for u = 1:nUnit
        D = plaid_mean_fr(u, :);
        P = pattern_prediction(u, :);
        C = component_prediction(u, :);

        valid = isfinite(D) & isfinite(P) & isfinite(C);

        nValid = sum(valid);
        n_valid_directions(u) = nValid;

        if nValid < 4
            continue;
        end

        Dv = D(valid);
        Pv = P(valid);
        Cv = C(valid);

        rp(u) = pearson_corr_safe(Dv, Pv);
        rc(u) = pearson_corr_safe(Dv, Cv);
        rpc(u) = pearson_corr_safe(Pv, Cv);

        if ~all(isfinite([rp(u), rc(u), rpc(u)]))
            continue;
        end

        denom_p = ...
            (1 - rc(u)^2) * ...
            (1 - rpc(u)^2);

        denom_c = ...
            (1 - rp(u)^2) * ...
            (1 - rpc(u)^2);

        if denom_p <= 0 || denom_c <= 0
            continue;
        end

        Rp(u) = ...
            (rp(u) - rc(u) * rpc(u)) / ...
            sqrt(denom_p);

        Rc(u) = ...
            (rc(u) - rp(u) * rpc(u)) / ...
            sqrt(denom_c);

        % Remove tiny numerical excursions outside the correlation range.
        Rp(u) = min(max(Rp(u), -1), 1);
        Rc(u) = min(max(Rc(u), -1), 1);

        partial_corr_valid(u) = true;

        % atanh is infinite at exactly +/-1. Clip only for the Fisher-Z
        % calculation so a numerically perfect fit remains plottable.
        clip_eps = 1e-12;

        Rp_for_z = ...
            min(max(Rp(u), -1 + clip_eps), 1 - clip_eps);

        Rc_for_z = ...
            min(max(Rc(u), -1 + clip_eps), 1 - clip_eps);

        fisher_z_clipped(u) = ...
            (Rp_for_z ~= Rp(u)) || ...
            (Rc_for_z ~= Rc(u));

        z_scale = sqrt(nValid - 3);

        Zp(u) = atanh(Rp_for_z) * z_scale;
        Zc(u) = atanh(Rc_for_z) * z_scale;
        PI(u) = Zp(u) - Zc(u);
    end

    valid_classification = ...
        isfinite(Zp) & ...
        isfinite(Zc) & ...
        isfinite(PI);

    is_pattern = ...
        valid_classification & ...
        (Zp > z_threshold) & ...
        ((Zp - Zc) > z_threshold);

    is_component = ...
        valid_classification & ...
        (Zc > z_threshold) & ...
        ((Zc - Zp) > z_threshold);

    is_unclassified = ...
        valid_classification & ...
        ~is_pattern & ...
        ~is_component;

    is_invalid = ~valid_classification;

    class_label = strings(nUnit, 1);
    class_label(is_pattern) = 'pattern';
    class_label(is_component) = 'component';
    class_label(is_unclassified) = 'unclassified';
    class_label(is_invalid) = 'invalid';

    out = struct();

    out.actual_response = plaid_mean_fr;
    out.actual_response_is_raw = true;

    out.pattern_prediction = pattern_prediction;
    out.pattern_prediction_is_raw = true;

    out.component_prediction = component_prediction;
    out.component_prediction_uses_baseline = ...
        subtract_spontaneous;

    out.component_separation_deg = ...
        component_separation_deg;

    out.half_component_separation_deg = half_sep;
    out.interpolation_used = false;
    out.n_total_directions = nPdir;

    out.component_direction_1 = component_direction_1;
    out.component_direction_2 = component_direction_2;

    out.pattern_grating_index = pattern_grating_index;
    out.component1_grating_index = component1_grating_index;
    out.component2_grating_index = component2_grating_index;

    out.component_raw_response_1 = component_raw_response_1;
    out.component_raw_response_2 = component_raw_response_2;

    out.component_baseline_1 = component_baseline_1;
    out.component_baseline_2 = component_baseline_2;

    out.component_evoked_response_1 = ...
        component_evoked_response_1;

    out.component_evoked_response_2 = ...
        component_evoked_response_2;

    out.component_baseline_addback = ...
        component_baseline_addback;

    out.rp = rp;
    out.rc = rc;
    out.rpc = rpc;

    out.Rp = Rp;
    out.Rc = Rc;

    out.Zp = Zp;
    out.Zc = Zc;
    out.PI = PI;

    out.n_valid_directions = n_valid_directions;
    out.partial_corr_valid = partial_corr_valid;
    out.fisher_z_clipped = fisher_z_clipped;

    out.pattern_z_threshold = z_threshold;

    out.is_pattern = is_pattern;
    out.is_component = is_component;
    out.is_unclassified = is_unclassified;
    out.is_invalid = is_invalid;
    out.class_label = class_label;
end

function idx = require_unique_direction_index( ...
        measured_directions, required_direction, context_text)
% Return the unique measured-direction index required by one calculation.
% Missing directions are errors because interpolation is intentionally off.

    measured_directions = measured_directions(:)';
    required_direction = mod(required_direction, 360);

    match = find(angle_equal( ...
        measured_directions, required_direction));

    if isempty(match)
        error(['Cannot construct %s because required grating direction ' ...
            '%.12g deg was not measured. Interpolation is disabled.'], ...
            context_text, required_direction);
    end

    if numel(match) > 1
        error(['Cannot construct %s because required grating direction ' ...
            '%.12g deg matches multiple measured direction columns.'], ...
            context_text, required_direction);
    end

    idx = match;
end

function r = pearson_corr_safe(x, y)
% Pearson correlation with explicit checks for constant vectors.

    x = double(x(:));
    y = double(y(:));

    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    if numel(x) < 2
        r = NaN;
        return;
    end

    x = x - mean(x);
    y = y - mean(y);

    denom = sqrt(sum(x.^2) * sum(y.^2));

    if denom <= 0 || ~isfinite(denom)
        r = NaN;
        return;
    end

    r = sum(x .* y) / denom;
    r = min(max(r, -1), 1);
end

function order = get_depth_desc_order(unit_depth_um, nUnit)
% Sort by descending depth; ties retain original unit order; NaNs last.

    unit_depth_um = unit_depth_um(:);

    if numel(unit_depth_um) ~= nUnit
        unit_depth_um = nan(nUnit, 1);
    end

    unit_index = (1:nUnit)';
    finite_depth = isfinite(unit_depth_um);

    finite_idx = unit_index(finite_depth);
    finite_depth_val = unit_depth_um(finite_depth);

    if isempty(finite_idx)
        order = unit_index;
        return;
    end

    sort_table = [-finite_depth_val, finite_idx];
    [~, sort_pos] = sortrows(sort_table, [1 2]);

    finite_order = finite_idx(sort_pos);
    nan_order = unit_index(~finite_depth);

    order = [finite_order; nan_order];
end

function err_matrix = select_error_matrix(tuning_struct, error_method)
% Return either SEM or SD for plotting.

    switch upper(error_method)
        case 'SEM'
            err_matrix = tuning_struct.sem_fr;

        case 'SD'
            err_matrix = tuning_struct.sd_fr;

        otherwise
            error('Unknown error_method: %s', error_method);
    end
end

function fig = plot_direction_tuning( ...
        directions, mean_fr, err_fr, plot_order, ...
        figure_title, visibility)
% Plot one compact direction-tuning subplot per neuron.
%
% Each neuron uses its own y-axis range, determined from that neuron's
% mean tuning curve and the displayed mean +/- SD or SEM error band.
%
% The lower limit is not forced to zero. A small padding is added above
% and below the displayed data so that the curve and error band are not
% clipped by the subplot boundaries.

    [nUnit, nDir] = size(mean_fr);
    directions = directions(:)';

    if numel(directions) ~= nDir
        error('Direction count does not match tuning matrix columns.');
    end

    if ~isequal(size(err_fr), size(mean_fr))
        error('Error matrix size does not match mean tuning matrix size.');
    end

    if numel(plot_order) ~= nUnit
        error('plot_order length does not match neuron count.');
    end

    nCols = ceil(sqrt(nUnit));
    nRows = ceil(nUnit / nCols);

    fig_width = max(900, 180 * nCols);
    fig_height = max(700, 145 * nRows);

    fig = figure( ...
        'Color', 'w', ...
        'Visible', visibility, ...
        'Position', [80 80 fig_width fig_height]);

    sgtitle(figure_title, 'Interpreter', 'none');

    % Embedded vertical padding for every individual neuron.
    ylim_padding_fraction = 0.10;

    for ii = 1:nUnit

        u = plot_order(ii);

        ax = subplot(nRows, nCols, ii, 'Parent', fig);
        hold(ax, 'on');

        r_mean = mean_fr(u, :);
        r_err = err_fr(u, :);

        % Do not truncate the lower error-band boundary at zero.
        r_upper = r_mean + r_err;
        r_lower = r_mean - r_err;

        %% Shaded SD or SEM region
        valid_band = ...
            isfinite(directions) & ...
            isfinite(r_upper) & ...
            isfinite(r_lower);

        if sum(valid_band) >= 2

            x_band = directions(valid_band);
            upper_use = r_upper(valid_band);
            lower_use = r_lower(valid_band);

            hBand = fill( ...
                ax, ...
                [x_band fliplr(x_band)], ...
                [upper_use fliplr(lower_use)], ...
                [0.75 0.75 0.75], ...
                'EdgeColor', 'none');

            uistack(hBand, 'bottom');
        end

        %% Mean tuning curve
        plot( ...
            ax, ...
            directions, ...
            r_mean, ...
            'Color', [0 0 0], ...
            'LineWidth', 1);

        xlim(ax, [0 360]);

        %% Set an independent y-axis range for this neuron
        local_values = [ ...
            r_mean(:); ...
            r_lower(:); ...
            r_upper(:)];

        local_values = local_values(isfinite(local_values));

        if isempty(local_values)

            local_ylim = [-0.5 0.5];

        else

            y_min = min(local_values);
            y_max = max(local_values);
            y_range = y_max - y_min;

            % Handle completely flat tuning curves.
            if y_range <= eps(max(abs([y_min, y_max, 1])))

                reference_scale = max(abs([y_min, y_max]));

                if reference_scale == 0
                    extra = 0.5;
                else
                    extra = 0.10 * reference_scale;
                end

            else

                extra = ylim_padding_fraction * y_range;

            end

            local_ylim = [ ...
                y_min - extra, ...
                y_max + extra];

            % Final protection against an invalid zero-width range.
            if local_ylim(2) <= local_ylim(1)
                local_ylim = [y_min - 0.5, y_max + 0.5];
            end
        end

        ylim(ax, local_ylim);

        box(ax, 'off');

        set( ...
            ax, ...
            'XTick', [], ...
            'YTick', []);
    end
end

function fig = plot_di_osi_distributions( ...
        DI, OSI, bin_edges, bar_width, ...
        group_name, visibility, baseline_subtracted)
% Plot DI and OSI probability distributions in a two-row figure.
% Bars are gray, separated, and outlined in black.
%
% In spontaneous-subtracted mode, regular DI bins cover 0 to 1.6.
% Every finite DI value greater than 1.6 is combined into one overflow bar
% labeled >1.6. Mean and median are calculated from all finite values,
% including the values represented by the overflow bar.
%
% Mean is marked by a solid vertical line. Median is marked by a dashed
% vertical line. Their numerical values are printed in each panel.

    DI = double(DI(:));
    OSI = double(OSI(:));

    DI = DI(isfinite(DI));
    OSI = OSI(isfinite(OSI));

    baseline_subtracted = logical(baseline_subtracted);
    bin_edges = double(bin_edges(:)');

    if numel(bin_edges) < 2 || ...
            any(~isfinite(bin_edges)) || ...
            any(diff(bin_edges) <= 0)
        error('bin_edges must be a strictly increasing finite vector.');
    end

    fig = figure( ...
        'Color', 'w', ...
        'Visible', visibility, ...
        'Position', [120 80 760 850]);

    tl = tiledlayout( ...
        fig, ...
        2, 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    %% DI distribution
    ax1 = nexttile(tl, 1);

    if baseline_subtracted
        di_regular_max = 1.6;

        [di_overflow_center, di_display_right] = ...
            plot_probability_bars_with_upper_overflow( ...
            ax1, ...
            DI, ...
            bin_edges, ...
            bar_width, ...
            di_regular_max);

        add_mean_median_markers( ...
            ax1, ...
            DI, ...
            0, ...
            di_display_right, ...
            di_regular_max, ...
            di_overflow_center);
    else
        plot_probability_bars( ...
            ax1, ...
            DI, ...
            bin_edges, ...
            bar_width);

        xlim(ax1, [bin_edges(1), bin_edges(end)]);

        add_mean_median_markers( ...
            ax1, ...
            DI, ...
            bin_edges(1), ...
            bin_edges(end), ...
            [], ...
            []);
    end

    xlabel(ax1, 'Direction index (DI)');
    ylabel(ax1, 'Proportion of cells');

    title( ...
        ax1, ...
        sprintf('DI distribution, n = %d', numel(DI)));

    box(ax1, 'off');

    %% OSI distribution
    ax2 = nexttile(tl, 2);

    plot_probability_bars( ...
        ax2, ...
        OSI, ...
        bin_edges, ...
        bar_width);

    xlim(ax2, [bin_edges(1), bin_edges(end)]);

    add_mean_median_markers( ...
        ax2, ...
        OSI, ...
        bin_edges(1), ...
        bin_edges(end), ...
        [], ...
        []);

    xlabel(ax2, 'Orientation selectivity index (OSI)');
    ylabel(ax2, 'Proportion of cells');

    title( ...
        ax2, ...
        sprintf('OSI distribution, n = %d', numel(OSI)));

    box(ax2, 'off');

    title( ...
        tl, ...
        sprintf('%s grating selectivity', group_name), ...
        'Interpreter', 'none');
end

function [overflow_center, display_right] = ...
        plot_probability_bars_with_upper_overflow( ...
        ax, values, base_bin_edges, bar_width, upper_limit)
% Plot regular probability bars through upper_limit and one overflow bar.
%
% The overflow bar is separated from the regular bins and represents all
% finite values strictly greater than upper_limit.

    values = double(values(:));
    values = values(isfinite(values));

    base_bin_edges = double(base_bin_edges(:)');
    bin_widths = diff(base_bin_edges);
    bin_width = median(bin_widths);

    if ~isfinite(bin_width) || bin_width <= 0 || ...
            any(abs(bin_widths - bin_width) > 1e-10)
        error(['base_bin_edges must have a constant positive spacing for ' ...
            'the DI overflow histogram.']);
    end

    if abs(base_bin_edges(1)) > 1e-10
        error('DI overflow histogram requires base_bin_edges to begin at 0.');
    end

    n_regular_bins_float = upper_limit / bin_width;
    n_regular_bins = round(n_regular_bins_float);

    if abs(n_regular_bins_float - n_regular_bins) > 1e-10
        error(['DI upper limit %.12g is not an integer multiple of the ' ...
            'histogram bin width %.12g.'], upper_limit, bin_width);
    end

    regular_edges = (0:n_regular_bins) * bin_width;
    regular_edges(end) = upper_limit;
    regular_centers = ...
        (regular_edges(1:end-1) + regular_edges(2:end)) / 2;

    regular_counts = histcounts(values, regular_edges);

    % DI should not be negative because Ropt is the maximal response, but
    % include any unexpected negative finite values in the first bar so the
    % displayed probabilities still sum to one.
    regular_counts(1) = regular_counts(1) + sum(values < regular_edges(1));

    overflow_count = sum(values > upper_limit);

    if isempty(values)
        regular_probability = zeros(size(regular_counts));
        overflow_probability = 0;
    else
        regular_probability = regular_counts / numel(values);
        overflow_probability = overflow_count / numel(values);
    end

    hold(ax, 'on');

    % Leave one empty-bin-width gap before the overflow category.
    overflow_center = upper_limit + 1.5 * bin_width;

    all_centers = [regular_centers, overflow_center];
    all_probability = [regular_probability, overflow_probability];

    bar( ...
        ax, ...
        all_centers, ...
        all_probability, ...
        bar_width, ...
        'FaceColor', [0.55 0.55 0.55], ...
        'EdgeColor', [0 0 0], ...
        'LineWidth', 0.8);

    display_right = overflow_center + 0.6 * bin_width;
    xlim(ax, [0, display_right]);

    numeric_ticks = [0, 0.4, 0.8, 1.2, upper_limit];
    numeric_ticks = numeric_ticks( ...
        numeric_ticks >= 0 & numeric_ticks <= upper_limit);

    tick_values = [numeric_ticks, overflow_center];
    tick_labels = arrayfun( ...
        @(x) sprintf('%g', x), ...
        numeric_ticks, ...
        'UniformOutput', false);

    tick_labels{end + 1} = sprintf('>%g', upper_limit);

    set( ...
        ax, ...
        'XTick', tick_values, ...
        'XTickLabel', tick_labels);
end

function add_mean_median_markers( ...
        ax, values, display_left, display_right, ...
        overflow_limit, overflow_center)
% Mark and print the mean and median of all finite values.
%
% For an overflow histogram, a center value above overflow_limit is drawn
% at overflow_center, while the text reports its true untruncated value.

    values = double(values(:));
    values = values(isfinite(values));

    if isempty(values)
        return;
    end

    mean_value = mean(values);
    median_value = median(values);

    mean_plot_x = center_value_to_display_position( ...
        mean_value, display_left, display_right, ...
        overflow_limit, overflow_center);

    median_plot_x = center_value_to_display_position( ...
        median_value, display_left, display_right, ...
        overflow_limit, overflow_center);

    xline( ...
        ax, ...
        mean_plot_x, ...
        '-', ...
        'Color', [0 0 0], ...
        'LineWidth', 1.3, ...
        'HandleVisibility', 'off');

    xline( ...
        ax, ...
        median_plot_x, ...
        '--', ...
        'Color', [0 0 0], ...
        'LineWidth', 1.3, ...
        'HandleVisibility', 'off');

    text( ...
        ax, ...
        0.98, ...
        0.96, ...
        sprintf(['Mean = %.3f (solid)\n' ...
                 'Median = %.3f (dashed)'], ...
                 mean_value, median_value), ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 10, ...
        'Interpreter', 'none');
end

function x_display = center_value_to_display_position( ...
        x_value, display_left, display_right, ...
        overflow_limit, overflow_center)
% Map a mean/median value to its visible x position.

    if ~isempty(overflow_limit) && ...
            ~isempty(overflow_center) && ...
            x_value > overflow_limit

        x_display = overflow_center;
        return;
    end

    x_display = min(max(x_value, display_left), display_right);
end

function fig = plot_pattern_component_analysis( ...
        Zp, Zc, PI, ...
        is_pattern, is_component, is_unclassified, ...
        z_threshold, bar_width, group_name, visibility)
% Plot Zp versus Zc and the PI probability distribution.
%
% Classification boundaries are shown in the scatter but are deliberately
% excluded from the legend. The PI panel contains no reference lines.

    valid_z = isfinite(Zp) & isfinite(Zc);
    valid_pi = isfinite(PI);

    all_z = [Zp(valid_z); Zc(valid_z)];

    if isempty(all_z)
        z_limits = [-1, 4];
    else
        z_lo = min(-1, floor(min(all_z) - 0.5));
        z_hi = max(4, ceil(max(all_z) + 0.5));

        if z_hi <= z_lo
            z_hi = z_lo + 1;
        end

        z_limits = [z_lo, z_hi];
    end

    fig = figure( ...
        'Color', 'w', ...
        'Visible', visibility, ...
        'Position', [100 180 1150 500]);

    tl = tiledlayout( ...
        fig, ...
        1, 2, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    %% Zp-Zc scatter
    ax1 = nexttile(tl, 1);
    hold(ax1, 'on');

    scatter( ...
        ax1, ...
        Zc(is_unclassified), ...
        Zp(is_unclassified), ...
        42, ...
        'o', ...
        'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', [0.25 0.25 0.25], ...
        'LineWidth', 1, ...
        'DisplayName', 'Unclassified');

    scatter( ...
        ax1, ...
        Zc(is_component), ...
        Zp(is_component), ...
        42, ...
        'o', ...
        'MarkerFaceColor', [0.55 0.55 0.55], ...
        'MarkerEdgeColor', [0 0 0], ...
        'LineWidth', 0.8, ...
        'DisplayName', 'Component');

    scatter( ...
        ax1, ...
        Zc(is_pattern), ...
        Zp(is_pattern), ...
        42, ...
        'o', ...
        'MarkerFaceColor', [0 0 0], ...
        'MarkerEdgeColor', [0 0 0], ...
        'DisplayName', 'Pattern');

    plot_z_classification_boundaries( ...
        ax1, ...
        z_limits, ...
        z_threshold);

    xlim(ax1, z_limits);
    ylim(ax1, z_limits);
    axis(ax1, 'square');

    xlabel(ax1, 'Component correlation (Z_c)');
    ylabel(ax1, 'Pattern correlation (Z_p)');

    title( ...
        ax1, ...
        sprintf( ...
        'Pattern/component classification, n = %d', ...
        sum(valid_z)));

    legend(ax1, 'Location', 'best');
    box(ax1, 'off');

    %% PI distribution
    ax2 = nexttile(tl, 2);

    PI_valid = PI(valid_pi);

    if isempty(PI_valid)
        pi_limit = 2;
    else
        pi_limit = max( ...
            2, ...
            ceil(max(abs(PI_valid)) + 0.5));
    end

    pi_edges = linspace( ...
        -pi_limit, ...
        pi_limit, ...
        13);

    plot_probability_bars( ...
        ax2, ...
        PI_valid, ...
        pi_edges, ...
        bar_width);

    xlim(ax2, [-pi_limit, pi_limit]);

    xlabel(ax2, 'Pattern index (Z_p - Z_c)');
    ylabel(ax2, 'Proportion of cells');

    title( ...
        ax2, ...
        sprintf( ...
        'PI distribution, n = %d', ...
        numel(PI_valid)));

    box(ax2, 'off');

    title( ...
        tl, ...
        sprintf( ...
        '%s plaid pattern/component analysis', ...
        group_name), ...
        'Interpreter', 'none');
end

function plot_probability_bars( ...
        ax, values, bin_edges, bar_width)
% Draw a probability histogram using bar so adjacent bins have visible gaps.

    values = double(values(:));
    values = values(isfinite(values));

    bin_edges = double(bin_edges(:)');

    if numel(bin_edges) < 2 || ...
            any(~isfinite(bin_edges)) || ...
            any(diff(bin_edges) <= 0)
        error('bin_edges must be a strictly increasing finite vector.');
    end

    counts = histcounts(values, bin_edges);

    if isempty(values)
        probability = zeros(size(counts));
    else
        probability = counts / numel(values);
    end

    centers = ...
        (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

    bar( ...
        ax, ...
        centers, ...
        probability, ...
        bar_width, ...
        'FaceColor', [0.55 0.55 0.55], ...
        'EdgeColor', [0 0 0], ...
        'LineWidth', 0.8);
end

function plot_z_classification_boundaries( ...
        ax, z_limits, threshold)
% Draw the conventional piecewise Z-threshold classification boundaries.
%
% Pattern region requires:
%   Zp > threshold and Zp - Zc > threshold
%
% Component region requires:
%   Zc > threshold and Zc - Zp > threshold
%
% Boundary lines are not included in the legend.

    lo = z_limits(1);
    hi = z_limits(2);

    line_args = { ...
        'LineStyle', '--', ...
        'Color', [0 0 0], ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off'};

    % Pattern boundary:
    % horizontal for Zc <= 0, diagonal for Zc >= 0.
    if lo < 0 && ...
            threshold >= lo && ...
            threshold <= hi

        plot( ...
            ax, ...
            [lo, min(0, hi)], ...
            [threshold, threshold], ...
            line_args{:});
    end

    x_start = max(0, lo);
    x_end = min(hi, hi - threshold);

    if x_end > x_start

        x = linspace( ...
            x_start, ...
            x_end, ...
            200);

        plot( ...
            ax, ...
            x, ...
            x + threshold, ...
            line_args{:});
    end

    % Component boundary:
    % vertical for Zp <= 0, diagonal for Zp >= 0.
    if lo < 0 && ...
            threshold >= lo && ...
            threshold <= hi

        plot( ...
            ax, ...
            [threshold, threshold], ...
            [lo, min(0, hi)], ...
            line_args{:});
    end

    x_start = max(threshold, lo);
    x_end = hi;

    if x_end > x_start

        x = linspace( ...
            x_start, ...
            x_end, ...
            200);

        y = x - threshold;

        keep = ...
            y >= lo & ...
            y <= hi;

        plot( ...
            ax, ...
            x(keep), ...
            y(keep), ...
            line_args{:});
    end
end

function save_figure_pair(fig, png_file, fig_file)
% Save one MATLAB figure as both a high-resolution PNG and a FIG file.

    savefig(fig, fig_file);
    print(fig, png_file, '-dpng', '-r300');
end