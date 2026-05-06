%% =========================================================================
% plot_size_rawcount_scatter
%
% Purpose:
%   For one selected stim_tag/run in model_data_allruns.mat:
%   1. Use raw_count_by_condition.
%   2. Pool trials by condition size: small vs large.
%   3. For each trial, compute total raw spike count across time bins:
%          total_count = sum(trial.y, 2)
%      where trial.y is nUnit x nTimeBin.
%   4. Average total_count across trials within each size.
%   5. Plot large vs small response separately for each group/probe.
%
% Each point = one cell.
% x-axis = small size mean total spike count per trial.
% y-axis = large size mean total spike count per trial.
% =========================================================================

clc; clear;

%% ----------------------- User parameters -----------------------

dat_file = 'I:\np_data\RafiL001p0120_g1\catgt_RafiL001p0120_g1\model_data_allruns.mat';

stim_tag = '_2[Gpl2_2c_2sz_400_2_200isi]';

data_content = 'raw_count';
by_condition_field = sprintf('%s_by_condition', data_content);

save_fig = false;
fig_out_dir = pwd;

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
fprintf('  small size %g: %d trials\n', small_size, n_small_trials);
fprintf('  large size %g: %d trials\n', large_size, n_large_trials);

%% ----------------------- Optional sanity check for raw_count NaN -----------------------

if any(isnan(small_mean_count)) || any(isnan(large_mean_count))
    warning('NaN found in computed mean responses. Please inspect raw_count data.');
end

%% ----------------------- Plot by group -----------------------

if save_fig && ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

nGroups = numel(groupd);
group_start = 1;

figure('Color', 'w', 'Name', 'Small vs large raw count by group');

for g = 1:nGroups
    group_end = group_start + groupd(g) - 1;
    group_idx = group_start:group_end;

    x = small_mean_count(group_idx);
    y = large_mean_count(group_idx);

    subplot(1, nGroups, g);
    scatter(x, y, 40, 'filled');
    hold on;

    max_val = max([x(:); y(:)]);
    if isempty(max_val) || ~isfinite(max_val)
        max_val = 1;
    end

    plot([0 max_val], [0 max_val], 'k--', 'LineWidth', 1);

    axis square;
    xlim([0 max_val * 1.05]);
    ylim([0 max_val * 1.05]);

    xlabel(sprintf('Small size %g: mean total spike count per trial', small_size));
    ylabel(sprintf('Large size %g: mean total spike count per trial', large_size));

    title(sprintf('Group %d, n = %d cells', g, groupd(g)), 'Interpreter', 'none');

    grid on;
    box off;

    group_start = group_end + 1;
end

sgtitle(sprintf('%s raw_count total over time bins', stim_tag), ...
    'Interpreter', 'none');

%% ----------------------- Save figure if requested -----------------------

if save_fig
    safe_tag = regexprep(stim_tag, '[^\w]', '_');
    fig_file = fullfile(fig_out_dir, sprintf('small_vs_large_rawcount_%s.png', safe_tag));

    exportgraphics(gcf, fig_file, 'Resolution', 300);
    fprintf('\nSaved figure:\n  %s\n', fig_file);
end

fprintf('\nDone.\n');

%% =========================================================================
% Local functions
% =========================================================================

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

        y = trials(t).y;   % nUnit x nTimeBin

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