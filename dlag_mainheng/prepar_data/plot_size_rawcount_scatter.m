%% =========================================================================
% plot_size_rawcount_scatter
%
% Purpose:
% For one selected stim_tag/run in model_data_allruns.mat:
% 1. Use raw_count_by_condition.
% 2. Pool trials by condition size: small vs large.
% 3. For each trial, compute total raw spike count across time bins:
%       total_count = sum(trial.y, 2)
%    where trial.y is nUnit x nTimeBin.
% 4. Average total_count across trials within each size.
% 5. Plot large vs small response separately for each group/probe.
%
% Each point = one cell.
% x-axis = small size mean total spike count per trial.
% y-axis = large size mean total spike count per trial.
%
% This version saves/displays two figures:
%   1) Full-range linear axis.
%   2) Clean broken-axis display.
% =========================================================================

clc; clear;

%% ----------------------- User parameters -----------------------
dat_file = 'I:\np_data\RafiL001p0121_g1\catgt_RafiL001p0121_g1\model_data_allruns.mat';
stim_tag = '_2[Gpl2_2c_2sz_400_2_200isi]';

data_content = 'raw_count';
by_condition_field = sprintf('%s_by_condition', data_content);

% Save options.
% save_fig is the master switch for saving files.
save_fig = true;
fig_out_dir = pwd;
save_png = true;
save_matlab_fig = false;

% Figure appearance.
fig_position = [100 100 1800 700];
marker_size = 40;
marker_face_alpha = 0.45;
marker_edge_alpha = 1;

% Broken-axis options.
% break_start is estimated from this percentile of values in each group.
break_start_prctile = 98.0;

% Broken axis is used only when max_val > break_start * this ratio.
broken_axis_trigger_ratio = 1;

% Display length assigned to the high-value tail, relative to break_start.
tail_display_frac = 0.08;

% Display gap around the broken part, relative to break_start.
break_gap_frac = 0.03;

%% ----------------------- Load data -----------------------
fprintf('Reading from %s\n', dat_file);
S = load(dat_file);

if ~isfield(S, 'model_data_allruns')
    error('model_data_allruns not found in %s.', dat_file);
end

model_data_allruns = S.model_data_allruns;

%% ----------------------- Select run by stim_tag -----------------------
all_run_tags = get_all_run_tags_local(model_data_allruns);
run_idx = find(strcmp(all_run_tags, stim_tag));

if isempty(run_idx)
    error('Requested stim_tag not found: %s', stim_tag);
end

if numel(run_idx) > 1
    error('Duplicate stim_tag found: %s', stim_tag);
end

this_run = model_data_allruns{run_idx};

fprintf('Selected run index: %d\n', run_idx);
fprintf('Selected stim_tag: %s\n', this_run.stim_tag);

%% ----------------------- Check fields -----------------------
if ~isfield(this_run, by_condition_field)
    error('Field %s not found in selected run.', by_condition_field);
end

if ~isfield(this_run, 'groupd')
    error('Field groupd not found in selected run.');
end

cond_data = this_run.(by_condition_field);
groupd = this_run.groupd(:)';

if isempty(cond_data)
    error('%s is empty.', by_condition_field);
end

if ~isfield(cond_data, 'size')
    error('Condition field "size" not found in %s.', by_condition_field);
end

if ~isfield(cond_data, 'trials')
    error('Field "trials" not found in %s condition structs.', by_condition_field);
end

%% ----------------------- Determine small and large size -----------------------
size_values = [cond_data.size];
unique_sizes = unique(size_values);

fprintf('\nUnique size values found in this run:\n');
disp(unique_sizes);

if numel(unique_sizes) ~= 2
    error(['Expected exactly 2 unique size values for small vs large, ', ...
           'but found %d. Please inspect cond_data.size.'], numel(unique_sizes));
end

small_size = min(unique_sizes);
large_size = max(unique_sizes);

fprintf('Using small_size = %g\n', small_size);
fprintf('Using large_size = %g\n', large_size);

%% ----------------------- Compute response per cell -----------------------
nUnits = get_n_units_from_condition_trials(cond_data);

if sum(groupd) ~= nUnits
    error('sum(groupd) = %d, but raw_count has %d units.', sum(groupd), nUnits);
end

[small_mean_count, n_small_trials] = compute_mean_total_count_for_size( ...
    cond_data, small_size, nUnits);

[large_mean_count, n_large_trials] = compute_mean_total_count_for_size( ...
    cond_data, large_size, nUnits);

fprintf('\nNumber of pooled trials:\n');
fprintf(' small size %g: %d trials\n', small_size, n_small_trials);
fprintf(' large size %g: %d trials\n', large_size, n_large_trials);

%% ----------------------- Optional sanity check -----------------------
if any(isnan(small_mean_count)) || any(isnan(large_mean_count))
    warning('NaN found in computed mean responses. Please inspect raw_count data.');
end

%% ----------------------- Prepare output -----------------------
if save_fig && ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

safe_tag = regexprep(stim_tag, '[^\w]', '_');

%% ----------------------- Plot 1: full-range linear axis -----------------------
hfig_full = figure('Color', 'w', ...
                   'Name', 'Small vs large raw count by group, full range', ...
                   'Position', fig_position);

plot_size_rawcount_groups_fullrange( ...
    small_mean_count, large_mean_count, groupd, ...
    small_size, large_size, stim_tag, ...
    marker_size, marker_face_alpha, marker_edge_alpha);

if save_fig
    save_current_figure_local(hfig_full, fig_out_dir, ...
        sprintf('small_vs_large_rawcount_%s_fullrange', safe_tag), ...
        save_png, save_matlab_fig);
end

%% ----------------------- Plot 2: clean broken-axis display -----------------------
hfig_broken = figure('Color', 'w', ...
                     'Name', 'Small vs large raw count by group, clean broken axis', ...
                     'Position', fig_position);

plot_size_rawcount_groups_clean_brokenaxis( ...
    small_mean_count, large_mean_count, groupd, ...
    small_size, large_size, stim_tag, ...
    marker_size, marker_face_alpha, marker_edge_alpha, ...
    break_start_prctile, broken_axis_trigger_ratio, ...
    tail_display_frac, break_gap_frac);

if save_fig
    save_current_figure_local(hfig_broken, fig_out_dir, ...
        sprintf('small_vs_large_rawcount_%s_brokenaxis', safe_tag), ...
        save_png, save_matlab_fig);
end

fprintf('\nDone.\n');

%% =========================================================================
% Local functions
% =========================================================================

function plot_size_rawcount_groups_fullrange( ...
    small_mean_count, large_mean_count, groupd, ...
    small_size, large_size, stim_tag, ...
    marker_size, marker_face_alpha, marker_edge_alpha)

    nGroups = numel(groupd);
    group_start = 1;

    for g = 1:nGroups
        group_end = group_start + groupd(g) - 1;
        group_idx = group_start:group_end;

        x = small_mean_count(group_idx);
        y = large_mean_count(group_idx);

        valid = isfinite(x) & isfinite(y);
        x = x(valid);
        y = y(valid);

        subplot(1, nGroups, g);

        scatter(x, y, marker_size, 'filled', ...
            'MarkerFaceAlpha', marker_face_alpha, ...
            'MarkerEdgeAlpha', marker_edge_alpha);
        hold on;

        max_val = max([x(:); y(:)]);
        if isempty(max_val) || ~isfinite(max_val) || max_val <= 0
            max_val = 1;
        end

        plot([0 max_val], [0 max_val], 'k--', 'LineWidth', 1);

        axis square;
        xlim([0 max_val * 1.05]);
        ylim([0 max_val * 1.05]);

        % Full-range figure also uses clean integer ticks.
        raw_ticks = choose_prebreak_integer_ticks_local(max_val);
        xticks(raw_ticks);
        yticks(raw_ticks);
        xticklabels(compose_integer_tick_labels_local(raw_ticks));
        yticklabels(compose_integer_tick_labels_local(raw_ticks));

        xlabel(sprintf('Small size %g: mean total spike count per trial', small_size));
        ylabel(sprintf('Large size %g: mean total spike count per trial', large_size));
        title(sprintf('Group %d, n = %d cells', g, groupd(g)), ...
            'Interpreter', 'none');

        grid off;
        box off;

        group_start = group_end + 1;
    end

    sgtitle(sprintf('%s raw_count total over time bins', stim_tag), ...
        'Interpreter', 'none');
end

function plot_size_rawcount_groups_clean_brokenaxis( ...
    small_mean_count, large_mean_count, groupd, ...
    small_size, large_size, stim_tag, ...
    marker_size, marker_face_alpha, marker_edge_alpha, ...
    break_start_prctile, broken_axis_trigger_ratio, ...
    tail_display_frac, break_gap_frac)

    nGroups = numel(groupd);
    group_start = 1;

    for g = 1:nGroups
        group_end = group_start + groupd(g) - 1;
        group_idx = group_start:group_end;

        x = small_mean_count(group_idx);
        y = large_mean_count(group_idx);

        valid = isfinite(x) & isfinite(y);
        x = x(valid);
        y = y(valid);

        subplot(1, nGroups, g);

        all_val = [x(:); y(:)];
        all_val = all_val(isfinite(all_val));

        if isempty(all_val)
            all_val = 1;
        end

        max_val = max(all_val);
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

            xlabel(sprintf('Small size %g: mean total spike count per trial', small_size));
            ylabel(sprintf('Large size %g: mean total spike count per trial', large_size));
            title(sprintf('Group %d, n = %d cells, linear axis', g, groupd(g)), ...
                'Interpreter', 'none');

            grid off;
            box off;

        else
            % break_end is the first actual data value above break_start.
            % Therefore, the skipped interval usually contains no real data point.
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

            % y = x reference line. It is split into two parts so that no line
            % is drawn across the skipped interval.
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

            % Only show clean integer ticks before the break.
            % No ticks or tick labels are shown after the break.
            raw_ticks = choose_prebreak_integer_ticks_local(break_start);
            plot_ticks = broken_axis_transform_local( ...
                raw_ticks, break_start, break_end, max_val, tail_display_len, break_gap);

            xticks(plot_ticks);
            yticks(plot_ticks);
            xticklabels(compose_integer_tick_labels_local(raw_ticks));
            yticklabels(compose_integer_tick_labels_local(raw_ticks));

            xlabel(sprintf('Small size %g: mean total spike count per trial', small_size));
            ylabel(sprintf('Large size %g: mean total spike count per trial', large_size));
            title(sprintf('Group %d, n = %d cells', g, groupd(g)), ...
                'Interpreter', 'none');

            grid off;
            box off;

            draw_axis_break_marks_local(gca, break_start, break_gap, display_max);
        end

        group_start = group_end + 1;
    end

    sgtitle(sprintf('%s raw_count total over time bins', stim_tag), ...
        'Interpreter', 'none');
end

function v_plot = broken_axis_transform_local( ...
    v, break_start, break_end, max_val, tail_display_len, break_gap)
% Map raw values to display coordinates for a broken-axis plot.
%
% Raw coordinate:
%   [0, break_start]        : displayed normally.
%   (break_start, break_end): skipped visually.
%   [break_end, max_val]    : compressed into a short tail segment.
%
% Display coordinate:
%   [0, break_start] is unchanged.
%   tail starts after a small visual gap.

    v_plot = v;

    if max_val <= break_end
        return;
    end

    high_mask = v >= break_end;
    middle_mask = v > break_start & v < break_end;

    % Values in the skipped region should be rare because break_end is chosen
    % as the first actual value above break_start. If they exist, place them
    % at the break boundary.
    v_plot(middle_mask) = break_start;

    v_plot(high_mask) = break_start + break_gap + ...
        (v(high_mask) - break_end) ./ max(eps, max_val - break_end) .* tail_display_len;
end

function raw_ticks = choose_prebreak_integer_ticks_local(prebreak_max)
% Choose clean integer ticks before the break.
%
% Example outputs:
%   0 5 10 15
%   0 10 20 30
%   0 20 40 60
%
% No ticks are generated after the break.

    if isempty(prebreak_max) || ~isfinite(prebreak_max) || prebreak_max <= 0
        raw_ticks = 0;
        return;
    end

    max_tick = floor(prebreak_max);

    if max_tick < 1
        raw_ticks = 0;
        return;
    end

    % Aim for about 4 to 6 intervals before the break.
    target_intervals = 5;
    rough_step = max_tick / target_intervals;
    step = choose_nice_integer_step_local(rough_step);

    last_tick = floor(max_tick / step) * step;
    raw_ticks = 0:step:last_tick;

    % If the range is small, use every integer.
    if numel(raw_ticks) < 3 && max_tick >= 2
        step = 1;
        raw_ticks = 0:step:max_tick;
    end
end

function step = choose_nice_integer_step_local(rough_step)
% Choose a nice integer step such as 1, 2, 5, 10, 20, 50, 100.

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
% Format tick labels as integers.

    labels = cell(size(raw_ticks));

    for i = 1:numel(raw_ticks)
        labels{i} = sprintf('%d', round(raw_ticks(i)));
    end
end

function draw_axis_break_marks_local(ax, break_start, break_gap, display_max)
% Draw simple double slash marks on x and y axes to indicate the break.

    axes(ax); %#ok<LAXES>
    hold(ax, 'on');

    x_break_center = break_start + break_gap / 2;
    y_break_center = break_start + break_gap / 2;

    slash_dx = display_max * 0.012;
    slash_dy = display_max * 0.025;
    offset = display_max * 0.018;

    x0 = ax.XLim(1);
    y0 = ax.YLim(1);

    % x-axis break mark: two slashes near bottom axis.
    plot(ax, ...
        [x_break_center - slash_dx, x_break_center + slash_dx], ...
        [y0 + slash_dy, y0 - slash_dy], ...
        'k-', 'LineWidth', 1.2, 'Clipping', 'off');

    plot(ax, ...
        [x_break_center - slash_dx + offset, x_break_center + slash_dx + offset], ...
        [y0 + slash_dy, y0 - slash_dy], ...
        'k-', 'LineWidth', 1.2, 'Clipping', 'off');

    % y-axis break mark: two slashes near left axis.
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
% Save current figure as PNG and/or MATLAB FIG.

    if save_png
        png_file = fullfile(fig_out_dir, sprintf('%s.png', base_name));
        exportgraphics(hfig, png_file, 'Resolution', 300);
        fprintf('\nSaved PNG figure:\n %s\n', png_file);
    end

    if save_matlab_fig
        fig_file = fullfile(fig_out_dir, sprintf('%s.fig', base_name));
        savefig(hfig, fig_file);
        fprintf('Saved MATLAB figure:\n %s\n', fig_file);
    end
end

function all_tags = get_all_run_tags_local(model_data_allruns)
    all_tags = cell(numel(model_data_allruns), 1);

    for j = 1:numel(model_data_allruns)
        if ~isfield(model_data_allruns{j}, 'stim_tag')
            error('stim_tag missing in model_data_allruns{%d}.', j);
        end

        all_tags{j} = model_data_allruns{j}.stim_tag;
    end
end

function nUnits = get_n_units_from_condition_trials(cond_data)
    nUnits = [];

    for c = 1:numel(cond_data)
        trials = cond_data(c).trials;

        if isempty(trials)
            continue;
        end

        if ~isfield(trials, 'y')
            error('Field y missing in cond_data(%d).trials.', c);
        end

        nUnits = size(trials(1).y, 1);
        return;
    end

    error('No non-empty trials found in cond_data.');
end

function [mean_count, nTrials] = compute_mean_total_count_for_size(cond_data, target_size, nUnits)
    sum_count = zeros(nUnits, 1);
    nTrials = 0;

    for c = 1:numel(cond_data)
        if cond_data(c).size ~= target_size
            continue;
        end

        trials = cond_data(c).trials;

        for t = 1:numel(trials)
            y = trials(t).y; % nUnit x nTimeBin

            if size(y, 1) ~= nUnits
                error('Unit number mismatch in condition %d trial %d.', c, t);
            end

            % Total raw spike count across time bins for each cell.
            trial_total_count = sum(y, 2);

            sum_count = sum_count + trial_total_count;
            nTrials = nTrials + 1;
        end
    end

    if nTrials == 0
        error('No trials found for size = %g.', target_size);
    end

    mean_count = sum_count ./ nTrials;
end
