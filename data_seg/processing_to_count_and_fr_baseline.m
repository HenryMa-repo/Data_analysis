%% =========================================================================
% processing_to_count_and_fr_baseline
%
% Purpose:
% For selected run(s) identified by stim_tag, generate baseline spike-count
% and firing-rate data using a pre-stimulus analysis window.
%
% This script is based on processing_to_count_and_fr.m, but is specialized
% for baseline data:
%
% 1. Each requested run may use its own negative-time baseline window and
%    its own bin size.
%
% 2. Only the following response arrays are calculated and saved:
%       raw_count
%       raw_fr
%
%    No z-scored or demeaned response arrays are calculated.
%
% 3. Each run is saved in a separate MAT file. The filename contains:
%       stim_tag
%       analysis_window
%       bin_size
%
%    Therefore, baseline files generated with different windows or bin
%    sizes do not overwrite one another.
%
% 4. Each kilosort folder also contains:
%       baseline_data_index.mat
%
%    baseline_data_index maps stim_tag, analysis_window, and bin_size to
%    the corresponding baseline filename. Downstream programs should use
%    this index rather than parsing filenames.
%
% Required input in each kilosort folder:
% - spike_unit_time_trial.mat
% - unit_condition_metrics.mat
% - cluster_info.tsv
%
% Output for each requested run:
% - baseline_data_<stim-tag>_win_<start>_to_<end>_bin_<size>.mat
%
% Output index:
% - baseline_data_index.mat
%
% Baseline-file variable:
% - baseline_data
%
% Index-file variable:
% - baseline_data_index
%
% baseline_data fields:
% .data_type
% .stim_tag
% .source_run_index
% .unit_ids
% .unit_depth_um
% .unit_channel
% .analysis_window
% .bin_size
% .bin_edges
% .bin_centers
% .condition_fields
% .condition_index_per_trial
% .conditions
% .raw_count
% .raw_fr
%
% Notes:
% 1. raw_count and raw_fr are stored as:
%       unit x trial x bin
%
% 2. Time zero is stimulus onset. Baseline analysis windows must end at or
%    before zero.
%
% 3. The baseline window must be contained within the prestimulus coverage
%    stored in spike_unit_time_trial.mat.
%
% 4. The analysis-window duration must be an integer multiple of bin_size.
%
% 5. If the same stim_tag/window/bin-size combination is rerun, that exact
%    baseline version is updated. Other baseline versions are preserved.
% =========================================================================

clc;
clear;

addpath(genpath(fullfile('.', 'expo_tools')));
addpath(genpath(fullfile('.', 'utils')));

%% ----------------------- User parameters -----------------------
root_folder = 'I:\np_data';
runName = 'RafiL001p0120';
runind = 1;          % run index after -g
probes = [0,1];      % probe indices after -prb

% -------------------------------------------------------------------------
% One entry per requested run.
%
% Time zero is stimulus onset. For baseline data, analysis_window should
% lie before stimulus onset and therefore normally ends at zero.
%
% Each requested run may use a different analysis_window and bin_size.
% -------------------------------------------------------------------------
run_specs = struct([]);

run_specs(1).stim_tag = '[RFG_coarse2dg_99_4_150isi]';
run_specs(1).analysis_window = [-0.1 0];
run_specs(1).bin_size = 0.1;

run_specs(2).stim_tag = '[dir12_gpl_2_200isi_fixedphase]';
run_specs(2).analysis_window = [-0.1 0];
run_specs(2).bin_size = 0.1;

run_specs(3).stim_tag = '_2[Gpl2_2c_2sz_400_2_200isi]';
run_specs(3).analysis_window = [-0.1 0];
run_specs(3).bin_size = 0.1;

% Name of the index file saved in each kilosort folder.
baseline_index_filename = 'baseline_data_index.mat';

%% ----------------------- Validate run specifications -----------------------
if isempty(run_specs)
    error('run_specs is empty.');
end

for r = 1:numel(run_specs)

    if ~isfield(run_specs(r), 'stim_tag') || ...
            isempty(run_specs(r).stim_tag)
        error('run_specs(%d).stim_tag is missing or empty.', r);
    end

    if ~isfield(run_specs(r), 'analysis_window')
        error('run_specs(%d).analysis_window is missing.', r);
    end

    if ~isfield(run_specs(r), 'bin_size')
        error('run_specs(%d).bin_size is missing.', r);
    end

    analysis_window = double(run_specs(r).analysis_window);
    bin_size = double(run_specs(r).bin_size);

    if numel(analysis_window) ~= 2 || ...
            any(~isfinite(analysis_window))
        error(['run_specs(%d).analysis_window must contain two finite ' ...
            'numeric values.'], r);
    end

    analysis_window = analysis_window(:)';

    if analysis_window(2) <= analysis_window(1)
        error(['run_specs(%d).analysis_window must satisfy ' ...
            'end > start.'], r);
    end

    if analysis_window(1) >= 0
        error(['run_specs(%d).analysis_window does not begin before ' ...
            'stimulus onset.'], r);
    end

    if analysis_window(2) > 0
        error(['run_specs(%d).analysis_window ends after stimulus onset. ' ...
            'Baseline windows must end at or before zero.'], r);
    end

    if ~isscalar(bin_size) || ...
            ~isfinite(bin_size) || ...
            bin_size <= 0
        error('run_specs(%d).bin_size must be a positive finite scalar.', r);
    end

    % Also verifies that the window duration is an integer multiple of the
    % requested bin size.
    make_bin_definition(analysis_window, bin_size);

    run_specs(r).analysis_window = analysis_window;
    run_specs(r).bin_size = bin_size;
end

% Prevent duplicate specifications within the current run_specs array.
for r1 = 1:numel(run_specs)
    for r2 = (r1 + 1):numel(run_specs)

        same_tag = strcmp( ...
            run_specs(r1).stim_tag, ...
            run_specs(r2).stim_tag);

        same_window = same_numeric_vector_local( ...
            run_specs(r1).analysis_window, ...
            run_specs(r2).analysis_window);

        same_bin = abs( ...
            run_specs(r1).bin_size - ...
            run_specs(r2).bin_size) <= 1e-12;

        if same_tag && same_window && same_bin
            error(['run_specs(%d) and run_specs(%d) specify the same ' ...
                'stim_tag, analysis_window, and bin_size.'], r1, r2);
        end
    end
end

%% ----------------------- Build shared session paths -----------------------
run_g = sprintf('%s_g%d', runName, runind);
destDir = fullfile(root_folder, run_g);

fprintf('destDir: %s\n', destDir);

%% ----------------------- Process each probe folder -----------------------
for ip = 1:numel(probes)

    thisProbe = probes(ip);
    imecStr = sprintf('imec%d', thisProbe);

    probe_folder = fullfile( ...
        destDir, ...
        ['catgt_' run_g], ...
        [run_g '_' imecStr]);

    fprintf('\n============================================================\n');
    fprintf('Processing probe %d\n', thisProbe);
    fprintf('probe_folder: %s\n', probe_folder);
    fprintf('============================================================\n');

    if ~isfolder(probe_folder)
        warning( ...
            'probe_folder does not exist, skipping probe %d: %s', ...
            thisProbe, probe_folder);
        continue;
    end

    %% ----------------------- Find kilosort folders -----------------------
    d = dir(fullfile(probe_folder, 'kilosort*'));
    d = d([d.isdir]);

    if isempty(d)
        warning( ...
            'No kilosort* folders found under probe %d: %s', ...
            thisProbe, probe_folder);
        continue;
    end

    [~, idx] = sort(lower({d.name}));
    d = d(idx);

    fprintf( ...
        'Found %d kilosort folder(s) under probe %d.\n', ...
        numel(d), thisProbe);

    %% ----------------------- Process each kilosort folder -----------------------
    for i = 1:numel(d)

        ksDir = fullfile(d(i).folder, d(i).name);

        fprintf('\nProcessing probe %d, ksDir: %s\n', ...
            thisProbe, ksDir);

        try
            %% ----------------------- Load required files -----------------------
            trial_file = fullfile( ...
                ksDir, ...
                'spike_unit_time_trial.mat');

            if ~isfile(trial_file)
                error('Missing file: %s', trial_file);
            end

            S = load( ...
                trial_file, ...
                'spike_unit_time_trial', ...
                'prestim_t', ...
                'poststim_t');

            if ~isfield(S, 'spike_unit_time_trial')
                error( ...
                    'spike_unit_time_trial not found in %s', ...
                    trial_file);
            end

            if ~isfield(S, 'prestim_t')
                error('prestim_t not found in %s', trial_file);
            end

            if ~isfield(S, 'poststim_t')
                error('poststim_t not found in %s', trial_file);
            end

            spike_unit_time_trial = S.spike_unit_time_trial;
            prestim_t = double(S.prestim_t);

            % poststim_t is loaded and checked to preserve the input
            % requirements of processing_to_count_and_fr.m. It is not
            % needed for a strictly prestimulus baseline window.
            poststim_t = double(S.poststim_t); %#ok<NASGU>

            if ~isscalar(prestim_t) || ...
                    ~isfinite(prestim_t) || ...
                    prestim_t < 0
                error('prestim_t must be a finite nonnegative scalar.');
            end

            cond_file = fullfile( ...
                ksDir, ...
                'unit_condition_metrics.mat');

            if ~isfile(cond_file)
                error('Missing file: %s', cond_file);
            end

            Sc = load( ...
                cond_file, ...
                'unit_condition_metrics');

            if ~isfield(Sc, 'unit_condition_metrics')
                error( ...
                    'unit_condition_metrics not found in %s', ...
                    cond_file);
            end

            unit_condition_metrics = Sc.unit_condition_metrics;

            if numel(spike_unit_time_trial) ~= ...
                    numel(unit_condition_metrics)

                error([ ...
                    'Number of runs mismatch between ' ...
                    'spike_unit_time_trial (%d) and ' ...
                    'unit_condition_metrics (%d) in ksDir %s'], ...
                    numel(spike_unit_time_trial), ...
                    numel(unit_condition_metrics), ...
                    ksDir);
            end

            %% ----------------------- Load Phy unit metadata -----------------------
            % Read curated unit-level metadata from cluster_info.tsv.
            % The metadata are later aligned separately to the unit_ids
            % stored for each run.
            unit_meta = load_unit_meta_from_phy(ksDir);

            %% ----------------------- Load or initialize baseline index -----------------------
            baseline_index_file = fullfile( ...
                ksDir, ...
                baseline_index_filename);

            baseline_data_index = load_existing_baseline_index( ...
                baseline_index_file);

            %% ----------------------- Find all available run tags -----------------------
            all_run_tags = get_all_run_tags( ...
                unit_condition_metrics);

            %% ----------------------- Process requested runs -----------------------
            for r = 1:numel(run_specs)

                wanted_tag = run_specs(r).stim_tag;
                analysis_window = run_specs(r).analysis_window;
                bin_size = run_specs(r).bin_size;

                run_idx = find(strcmp( ...
                    all_run_tags, ...
                    wanted_tag));

                if isempty(run_idx)
                    error( ...
                        'Requested stim_tag not found in ksDir %s: %s', ...
                        ksDir, wanted_tag);
                end

                if numel(run_idx) > 1
                    error([ ...
                        'Duplicate stim_tag found in ' ...
                        'unit_condition_metrics for ksDir %s: %s'], ...
                        ksDir, wanted_tag);
                end

                fprintf('\n  Baseline run %d/%d\n', ...
                    r, numel(run_specs));
                fprintf('    stim_tag        : %s\n', wanted_tag);
                fprintf('    analysis_window : [%g %g] s\n', ...
                    analysis_window(1), analysis_window(2));
                fprintf('    bin_size        : %g s\n', bin_size);

                this_run_trials = ...
                    spike_unit_time_trial{run_idx};

                this_run_metrics = ...
                    unit_condition_metrics{run_idx};

                %% ----------------------- Calculate baseline data -----------------------
                baseline_data = process_one_run_baseline( ...
                    this_run_trials, ...
                    this_run_metrics, ...
                    analysis_window, ...
                    bin_size, ...
                    prestim_t, ...
                    run_idx, ...
                    wanted_tag, ...
                    unit_meta);

                %% ----------------------- Create version-specific filename -----------------------
                baseline_filename = make_baseline_filename( ...
                    wanted_tag, ...
                    analysis_window, ...
                    bin_size);

                baseline_file = fullfile( ...
                    ksDir, ...
                    baseline_filename);

                %% ----------------------- Create/update index entry -----------------------
                new_index_entry = make_baseline_index_entry( ...
                    baseline_data, ...
                    baseline_filename);

                % Update the in-memory index before saving the data file.
                % This also checks whether the generated filename collides
                % with a different baseline specification.
                updated_baseline_data_index = ...
                    upsert_baseline_index( ...
                    baseline_data_index, ...
                    new_index_entry);

                %% ----------------------- Save this run's baseline data -----------------------
                save( ...
                    baseline_file, ...
                    'baseline_data');

                %% ----------------------- Save updated index -----------------------
                baseline_data_index = ...
                    updated_baseline_data_index;

                save( ...
                    baseline_index_file, ...
                    'baseline_data_index');

                fprintf('    Saved baseline file:\n');
                fprintf('      %s\n', baseline_file);
                fprintf('    Updated baseline index:\n');
                fprintf('      %s\n', baseline_index_file);
            end

        catch ME
            fprintf(2, ...
                'Error in probe %d, ksDir %s\n', ...
                thisProbe, ksDir);

            fprintf(2, '%s\n', ...
                getReport(ME, 'extended', ...
                'hyperlinks', 'off'));
        end
    end
end

fprintf('\nDone.\n');

%% ======================= Local functions =======================

function unit_meta = load_unit_meta_from_phy(ksDir)
%% =========================================================================
% load_unit_meta_from_phy
%
% Purpose:
% Read Phy cluster_info.tsv from one kilosort folder and extract unit-level
% metadata needed by this script.
%
% Expected useful columns:
% - cluster_id
% - depth
% - ch
%
% Output:
% unit_meta.unit_ids
% unit_meta.depth_um
% unit_meta.channel
%
% If cluster_info.tsv or an optional metadata column is unavailable, the
% corresponding output is empty or NaN and a warning is printed.
% =========================================================================

    unit_meta = struct();
    unit_meta.unit_ids = [];
    unit_meta.depth_um = [];
    unit_meta.channel = [];

    cluster_info_file = fullfile( ...
        ksDir, ...
        'cluster_info.tsv');

    if ~isfile(cluster_info_file)
        warning([ ...
            'cluster_info.tsv not found in ksDir. ' ...
            'Unit depth/channel will be NaN: %s'], ...
            ksDir);
        return;
    end

    try
        opts = detectImportOptions( ...
            cluster_info_file, ...
            'FileType', 'text', ...
            'Delimiter', '\t');

        T = readtable( ...
            cluster_info_file, ...
            opts);

    catch
        T = readtable( ...
            cluster_info_file, ...
            'FileType', 'text', ...
            'Delimiter', '\t');
    end

    varnames = T.Properties.VariableNames;

    id_col = find_table_column( ...
        varnames, ...
        {'cluster_id', 'clusterid', 'id'});

    depth_col = find_table_column( ...
        varnames, ...
        {'depth', 'depth_um', 'y', 'ypos', 'y_pos'});

    ch_col = find_table_column( ...
        varnames, ...
        {'ch', 'channel', 'best_channel', 'peak_channel'});

    if isempty(id_col)
        warning([ ...
            'No cluster_id column found in cluster_info.tsv. ' ...
            'Unit depth/channel will be NaN: %s'], ...
            cluster_info_file);
        return;
    end

    unit_ids = table_column_to_numeric( ...
        T, id_col);

    if isempty(depth_col)
        warning([ ...
            'No depth column found in cluster_info.tsv. ' ...
            'unit_depth_um will be NaN: %s'], ...
            cluster_info_file);

        depth_um = nan(size(unit_ids));
    else
        depth_um = table_column_to_numeric( ...
            T, depth_col);
    end

    if isempty(ch_col)
        warning([ ...
            'No ch/channel column found in cluster_info.tsv. ' ...
            'unit_channel will be NaN: %s'], ...
            cluster_info_file);

        channel = nan(size(unit_ids));
    else
        channel = table_column_to_numeric( ...
            T, ch_col);
    end

    unit_meta.unit_ids = unit_ids(:);
    unit_meta.depth_um = depth_um(:);
    unit_meta.channel = channel(:);

    if numel(unit_meta.depth_um) ~= ...
            numel(unit_meta.unit_ids)

        warning([ ...
            'Depth column length mismatch in cluster_info.tsv. ' ...
            'unit_depth_um will be NaN: %s'], ...
            cluster_info_file);

        unit_meta.depth_um = ...
            nan(size(unit_meta.unit_ids));
    end

    if numel(unit_meta.channel) ~= ...
            numel(unit_meta.unit_ids)

        warning([ ...
            'Channel column length mismatch in cluster_info.tsv. ' ...
            'unit_channel will be NaN: %s'], ...
            cluster_info_file);

        unit_meta.channel = ...
            nan(size(unit_meta.unit_ids));
    end
end

function colname = find_table_column(varnames, candidates)
% Case-insensitive search for one table-column name.

    colname = '';

    if isempty(varnames)
        return;
    end

    var_lower = lower(varnames);

    for k = 1:numel(candidates)

        cand = lower(candidates{k});
        idx = find(strcmp(var_lower, cand), 1);

        if ~isempty(idx)
            colname = varnames{idx};
            return;
        end
    end
end

function x = table_column_to_numeric(T, colname)
% Convert one table column to a numeric column vector.

    raw = T.(colname);

    if isnumeric(raw)
        x = double(raw);

    elseif iscell(raw)
        x = str2double(string(raw));

    elseif isstring(raw)
        x = str2double(raw);

    elseif iscategorical(raw)
        x = str2double(string(raw));

    elseif ischar(raw)
        x = str2double(cellstr(raw));

    else
        try
            x = double(raw);
        catch
            x = str2double(string(raw));
        end
    end

    x = x(:);
end

function out = process_one_run_baseline( ...
        this_run_trials, ...
        this_run_metrics, ...
        analysis_window, ...
        bin_size, ...
        prestim_t, ...
        run_idx, ...
        stim_tag, ...
        unit_meta)
%% =========================================================================
% process_one_run_baseline
%
% Purpose:
% For one selected run, convert prestimulus spike times into baseline
% raw-count and raw-firing-rate matrices.
%
% No normalized or demeaned arrays are calculated.
% =========================================================================

    if ~isfield(this_run_metrics, 'unit_ids')
        error( ...
            'unit_ids missing in unit_condition_metrics{%d}.', ...
            run_idx);
    end

    if ~isfield( ...
            this_run_metrics, ...
            'condition_index_per_trial')

        error([ ...
            'condition_index_per_trial missing in ' ...
            'unit_condition_metrics{%d}.'], ...
            run_idx);
    end

    if ~isfield(this_run_metrics, 'conditions')
        error( ...
            'conditions missing in unit_condition_metrics{%d}.', ...
            run_idx);
    end

    if ~isfield(this_run_metrics, 'condition_fields')
        error( ...
            'condition_fields missing in unit_condition_metrics{%d}.', ...
            run_idx);
    end

    unit_ids = ...
        this_run_metrics.unit_ids(:);

    condition_index_per_trial = ...
        this_run_metrics.condition_index_per_trial(:);

    condition_fields = ...
        this_run_metrics.condition_fields(:);

    conditions_identity = extract_condition_identity( ...
        this_run_metrics.conditions, ...
        condition_fields);

    nUnit = numel(unit_ids);
    nTrial = numel(this_run_trials);

    if numel(condition_index_per_trial) ~= nTrial
        error([ ...
            'condition_index_per_trial length (%d) does not match ' ...
            'number of trials (%d) for run %s'], ...
            numel(condition_index_per_trial), ...
            nTrial, ...
            stim_tag);
    end

    [nBin, bin_edges, bin_centers] = ...
        make_bin_definition( ...
        analysis_window, ...
        bin_size);

    validate_baseline_window( ...
        analysis_window, ...
        prestim_t, ...
        stim_tag);

    %% ----------------------- Align Phy unit metadata -----------------------
    [unit_depth_um, unit_channel] = ...
        align_unit_meta_to_unit_ids( ...
        unit_ids, ...
        unit_meta, ...
        stim_tag);

    %% ----------------------- Build raw count matrix -----------------------
    % Dimensions:
    %   unit x trial x bin
    raw_count = zeros( ...
        nUnit, ...
        nTrial, ...
        nBin);

    for t = 1:nTrial

        tr = this_run_trials{t};

        if isempty(tr)
            error( ...
                'Encountered an empty trial matrix at trial %d for run %s.', ...
                t, stim_tag);
        end

        if ~isnumeric(tr) || size(tr, 2) < 2
            error([ ...
                'Trial %d for run %s must be a numeric matrix with ' ...
                'at least two columns: unit ID and spike time.'], ...
                t, stim_tag);
        end

        for u = 1:nUnit

            uid = unit_ids(u);

            spk_t = tr( ...
                tr(:,1) == uid, ...
                2);

            if isempty(spk_t)

                c = zeros(1, nBin);

            else

                spk_keep = spk_t( ...
                    spk_t >= analysis_window(1) & ...
                    spk_t <= analysis_window(2));

                c = histcounts( ...
                    spk_keep, ...
                    bin_edges);
            end

            raw_count(u, t, :) = ...
                reshape(c, [1 1 nBin]);
        end
    end

    raw_count = double(raw_count);

    % Firing rate in Hz for each individual bin.
    raw_fr = raw_count ./ bin_size;

    %% ----------------------- Pack output -----------------------
    out = struct();

    out.data_type = 'baseline';
    out.stim_tag = stim_tag;
    out.source_run_index = run_idx;

    out.unit_ids = unit_ids;

    % Same order as out.unit_ids.
    out.unit_depth_um = unit_depth_um;
    out.unit_channel = unit_channel;

    out.analysis_window = analysis_window;
    out.bin_size = bin_size;
    out.bin_edges = bin_edges;
    out.bin_centers = bin_centers;

    out.condition_fields = condition_fields;
    out.condition_index_per_trial = ...
        condition_index_per_trial;

    out.conditions = conditions_identity;

    out.raw_count = raw_count;
    out.raw_fr = raw_fr;
end

function [unit_depth_um, unit_channel] = ...
        align_unit_meta_to_unit_ids( ...
        unit_ids, unit_meta, stim_tag)
% Align metadata from cluster_info.tsv to the unit IDs used in one run.

    unit_ids = unit_ids(:);

    unit_depth_um = nan(size(unit_ids));
    unit_channel = nan(size(unit_ids));

    if isempty(unit_meta) || ...
            ~isfield(unit_meta, 'unit_ids') || ...
            isempty(unit_meta.unit_ids)

        warning([ ...
            'No unit metadata available for run %s. ' ...
            'unit_depth_um/unit_channel will be NaN.'], ...
            stim_tag);
        return;
    end

    meta_ids = unit_meta.unit_ids(:);

    [tf, loc] = ismember( ...
        unit_ids, ...
        meta_ids);

    if any(tf)
        unit_depth_um(tf) = ...
            unit_meta.depth_um(loc(tf));

        unit_channel(tf) = ...
            unit_meta.channel(loc(tf));
    end

    if any(~tf)

        missing_ids = unit_ids(~tf);

        warning([ ...
            '%d/%d unit_ids in run %s were not found in ' ...
            'cluster_info.tsv. Example missing unit_id: %g'], ...
            numel(missing_ids), ...
            numel(unit_ids), ...
            stim_tag, ...
            missing_ids(1));
    end
end

function [nBin, bin_edges, bin_centers] = ...
        make_bin_definition( ...
        analysis_window, bin_size)
% Construct bin edges and centers for one analysis window.

    analysis_window = double(analysis_window(:)');

    if numel(analysis_window) ~= 2 || ...
            any(~isfinite(analysis_window))
        error( ...
            'analysis_window must contain two finite values.');
    end

    if ~isscalar(bin_size) || ...
            ~isfinite(bin_size) || ...
            bin_size <= 0
        error( ...
            'bin_size must be a positive finite scalar.');
    end

    t0 = analysis_window(1);
    t1 = analysis_window(2);

    span = t1 - t0;

    if span <= 0
        error( ...
            'analysis_window must satisfy end > start.');
    end

    nBin_float = span / bin_size;
    nBin = round(nBin_float);

    tol = 1e-10;

    if abs(nBin_float - nBin) > tol
        error([ ...
            'analysis_window span (%.12g s) is not an integer ' ...
            'multiple of bin_size (%.12g s).'], ...
            span, bin_size);
    end

    if nBin < 1
        error( ...
            'The requested analysis window produces fewer than one bin.');
    end

    bin_edges = ...
        t0 + (0:nBin) * bin_size;

    % Force the final edge to equal the requested endpoint exactly.
    bin_edges(end) = t1;

    bin_centers = ...
        bin_edges(1:end-1) + bin_size / 2;
end

function validate_baseline_window( ...
        analysis_window, prestim_t, stim_tag)
% Validate that the requested window is completely prestimulus and lies
% inside the available prestimulus coverage.

    analysis_window = double(analysis_window(:)');

    if numel(analysis_window) ~= 2 || ...
            any(~isfinite(analysis_window))
        error( ...
            'Invalid baseline analysis_window for run %s.', ...
            stim_tag);
    end

    if analysis_window(2) <= analysis_window(1)
        error([ ...
            'Baseline analysis_window must satisfy end > start ' ...
            'for run %s.'], ...
            stim_tag);
    end

    if analysis_window(1) >= 0
        error([ ...
            'Baseline analysis_window start %.6g is not before ' ...
            'stimulus onset for run %s.'], ...
            analysis_window(1), stim_tag);
    end

    if analysis_window(2) > 0
        error([ ...
            'Baseline analysis_window end %.6g is after stimulus ' ...
            'onset for run %s.'], ...
            analysis_window(2), stim_tag);
    end

    if analysis_window(1) < -prestim_t
        error([ ...
            'Baseline analysis_window start %.6g is earlier than ' ...
            'available prestimulus coverage ' ...
            '(-prestim_t = %.6g) for run %s.'], ...
            analysis_window(1), ...
            -prestim_t, ...
            stim_tag);
    end
end

function all_tags = get_all_run_tags( ...
        unit_condition_metrics)
% Extract stim_tag from every entry in unit_condition_metrics.

    all_tags = cell( ...
        numel(unit_condition_metrics), ...
        1);

    for j = 1:numel(unit_condition_metrics)

        if ~isfield( ...
                unit_condition_metrics{j}, ...
                'stim_tag')

            error( ...
                'stim_tag missing in unit_condition_metrics{%d}.', ...
                j);
        end

        all_tags{j} = ...
            unit_condition_metrics{j}.stim_tag;
    end
end

function conditions_identity = ...
        extract_condition_identity( ...
        conditions_in, condition_fields)
% Preserve condition identity and trial membership without retaining
% condition-derived response metrics.

    nCond = numel(conditions_in);

    conditions_identity = ...
        repmat(struct(), nCond, 1);

    for c = 1:nCond

        if ~isfield(conditions_in(c), 'trial_indices')
            error( ...
                'trial_indices is missing in conditions(%d).', ...
                c);
        end

        conditions_identity(c).trial_indices = ...
            conditions_in(c).trial_indices(:);

        for k = 1:numel(condition_fields)

            f = condition_fields{k};

            if isfield(conditions_in(c), f)

                conditions_identity(c).(f) = ...
                    conditions_in(c).(f);

            else

                error( ...
                    'Condition field %s is missing in conditions(%d).', ...
                    f, c);
            end
        end
    end
end

function filename = make_baseline_filename( ...
        stim_tag, analysis_window, bin_size)
% Create a readable, filesystem-safe filename containing the run identity,
% baseline window, and bin size.

    stim_file_tag = ...
        sanitize_stim_tag_for_filename(stim_tag);

    start_tag = ...
        make_number_tag_local(analysis_window(1));

    end_tag = ...
        make_number_tag_local(analysis_window(2));

    bin_tag = ...
        make_number_tag_local(bin_size);

    filename = sprintf( ...
        'baseline_data_%s_win_%s_to_%s_bin_%s.mat', ...
        stim_file_tag, ...
        start_tag, ...
        end_tag, ...
        bin_tag);
end

function tag = sanitize_stim_tag_for_filename(stim_tag)
% Convert a stim_tag into a readable filename component.
%
% The exact original stim_tag remains stored inside baseline_data and
% baseline_data_index. This sanitized form is only used in the filename.

    tag = char(string(stim_tag));
    tag = strtrim(tag);

    % Replace every group of non-ASCII alphanumeric characters with "_".
    tag = regexprep( ...
        tag, ...
        '[^A-Za-z0-9]+', ...
        '_');

    % Collapse repeated underscores.
    tag = regexprep( ...
        tag, ...
        '_+', ...
        '_');

    % Remove leading and trailing underscores.
    tag = regexprep( ...
        tag, ...
        '^_+|_+$', ...
        '');

    if isempty(tag)
        tag = 'stim';
    end
end

function tag = make_number_tag_local(x)
% Convert a numeric scalar into a compact filename-safe text tag.
%
% Examples:
%   -0.2 -> neg0p2
%    0   -> 0
%    0.05 -> 0p05

    if ~isscalar(x) || ~isnumeric(x) || ~isfinite(x)
        error( ...
            'Filename numeric tag requires one finite numeric scalar.');
    end

    tag = sprintf('%.12g', double(x));

    tag = strrep(tag, '-', 'neg');
    tag = strrep(tag, '+', 'pos');
    tag = strrep(tag, '.', 'p');

    if isempty(tag)
        tag = '0';
    end
end

function baseline_data_index = ...
        load_existing_baseline_index(index_file)
% Load an existing index so older baseline versions are preserved.
% Return an empty canonical index when the file does not yet exist.

    baseline_data_index = empty_baseline_index();

    if ~isfile(index_file)
        return;
    end

    S = load( ...
        index_file, ...
        'baseline_data_index');

    if ~isfield(S, 'baseline_data_index')
        error([ ...
            'The existing index file does not contain ' ...
            'baseline_data_index: %s'], ...
            index_file);
    end

    loaded_index = S.baseline_data_index;

    if ~isstruct(loaded_index)
        error([ ...
            'baseline_data_index in %s is not a struct array.'], ...
            index_file);
    end

    required_fields = fieldnames( ...
        empty_baseline_index());

    for k = 1:numel(required_fields)

        if ~isfield(loaded_index, required_fields{k})
            error([ ...
                'Existing baseline_data_index in %s is missing ' ...
                'required field: %s'], ...
                index_file, ...
                required_fields{k});
        end
    end

    baseline_data_index = loaded_index(:)';
end

function index = empty_baseline_index()
% Return an empty baseline index with a fixed field definition.

    index = struct( ...
        'stim_tag', {}, ...
        'analysis_window', {}, ...
        'bin_size', {}, ...
        'filename', {}, ...
        'source_run_index', {}, ...
        'n_units', {}, ...
        'n_trials', {}, ...
        'n_bins', {});
end

function entry = make_baseline_index_entry( ...
        baseline_data, baseline_filename)
% Create one index entry from a baseline_data structure.

    entry = struct();

    entry.stim_tag = ...
        char(string(baseline_data.stim_tag));

    entry.analysis_window = ...
        double(baseline_data.analysis_window(:)');

    entry.bin_size = ...
        double(baseline_data.bin_size);

    entry.filename = ...
        char(string(baseline_filename));

    entry.source_run_index = ...
        double(baseline_data.source_run_index);

    entry.n_units = ...
        size(baseline_data.raw_count, 1);

    entry.n_trials = ...
        size(baseline_data.raw_count, 2);

    entry.n_bins = ...
        size(baseline_data.raw_count, 3);
end

function baseline_data_index = ...
        upsert_baseline_index( ...
        baseline_data_index, new_entry)
% Add a new baseline version to the index, or replace the entry describing
% the exact same stim_tag, analysis_window, and bin_size.
%
% Different windows or bin sizes are preserved as separate entries.
% A filename collision between different specifications raises an error.

    if isempty(baseline_data_index)

        baseline_data_index = new_entry;
        return;
    end

    baseline_data_index = ...
        baseline_data_index(:)';

    nEntry = numel(baseline_data_index);
    exact_match = false(1, nEntry);
    filename_match = false(1, nEntry);

    for k = 1:nEntry

        same_tag = strcmp( ...
            baseline_data_index(k).stim_tag, ...
            new_entry.stim_tag);

        same_window = same_numeric_vector_local( ...
            baseline_data_index(k).analysis_window, ...
            new_entry.analysis_window);

        same_bin = abs( ...
            double(baseline_data_index(k).bin_size) - ...
            double(new_entry.bin_size)) <= 1e-12;

        exact_match(k) = ...
            same_tag && same_window && same_bin;

        filename_match(k) = strcmp( ...
            baseline_data_index(k).filename, ...
            new_entry.filename);
    end

    if sum(exact_match) > 1
        error([ ...
            'baseline_data_index contains duplicate entries for ' ...
            'stim_tag %s, window [%g %g], bin_size %g.'], ...
            new_entry.stim_tag, ...
            new_entry.analysis_window(1), ...
            new_entry.analysis_window(2), ...
            new_entry.bin_size);
    end

    % A filename may be reused only when it represents the exact same
    % baseline specification.
    conflicting_filename = ...
        filename_match & ~exact_match;

    if any(conflicting_filename)

        k = find(conflicting_filename, 1);

        error([ ...
            'Generated baseline filename collides with a different ' ...
            'index entry.\nFilename: %s\nExisting stim_tag: %s\n' ...
            'New stim_tag: %s'], ...
            new_entry.filename, ...
            baseline_data_index(k).stim_tag, ...
            new_entry.stim_tag);
    end

    if any(exact_match)

        % Rerunning the exact same version updates its metadata and file.
        baseline_data_index(exact_match) = ...
            new_entry;

    else

        % A new window/bin-size version is appended and older versions stay.
        baseline_data_index(end + 1) = ...
            new_entry;
    end
end

function tf = same_numeric_vector_local(a, b)
% Compare two numeric vectors using a small floating-point tolerance.

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