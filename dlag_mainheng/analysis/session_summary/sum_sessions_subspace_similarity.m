%% sum_sessions_subspace_similarity.m
% Pool the directional subspace-similarity results saved by
% subspace_similarity_dlag.m.
% Revision: 2026-08-07a
%
% One figure is made with this x-axis order:
%   Group 1: Across captures Within, Within captures Across,
%            Feedforward captures Feedback, Feedback captures Feedforward
%   Group 2: the same four directional comparisons
%
% all_condition_model:
%   Each session contributes one observation.
%
% condition_specific_models:
%   Each session-condition contributes one observation. Conditions are
%   pooled directly; no within-session averaging is performed first.

clc;
clear;

%% ========================== USER SETTINGS ==============================

root_dir = 'I:\np_data';

% -------------------------------------------------------------------------
% Data/model settings
% -------------------------------------------------------------------------

data_content = 'raw_count';
% Options usually include:
% raw_count, raw_fr, z_within_trial, z_within_condition,
% z_across_conditions, demean_count_within_trial,
% demean_fr_within_trial, demean_pooledsd_within_condition

model_mode = 'all_condition_model';
% Options:
% 'all_condition_model'
% 'condition_specific_models'

runIdx = 1;

% These options must match the options used in subspace_similarity_dlag.m.
use_dsl_filter = false;
dsl_field = 'logical';             % 'rawlogical' or 'logical'

use_svexp_filter = false;
svexp_field = 'logical';           % 'rawlogical' or 'logical'
shared_varexp_threshold = 0.95;

% -------------------------------------------------------------------------
% Plot settings
% -------------------------------------------------------------------------

plot_style = 'bar_points';
% Options:
% 'bar_points'    : mean bar + raw points + black median square
% 'violin'        : violin + raw points + median line
% 'points_median' : raw points + median line, without violin outline

% Whether to draw the mean for each plot style.
% Defaults preserve the behavior of the previous version:
% bar_points draws a mean bar; violin and points_median do not draw mean.
show_mean = struct( ...
    'bar_points', true, ...
    'violin', false, ...
    'points_median', false);

% Whether to draw the median for each plot style.
% In bar_points, median is a black square.
% In violin and points_median, median is a horizontal line.
show_median = struct( ...
    'bar_points', false, ...
    'violin', false, ...
    'points_median', false);

% Saving switches. save_figure controls both FIG and vector SVG output.
save_figure = true;
save_mat = true;

figure_visible = 'on';
close_after_save = false;

% Group 1 / Group 2 colors used in the existing analysis programs.
group_colors = [
    0.0000, 0.4470, 0.7410;
    0.8500, 0.3250, 0.0980
];

% Within-group categories are close; the two groups are farther apart.
within_group_spacing = 1.0;
between_group_spacing = 2.0;

y_limits = [0, 1];
figure_size_inches = [8.8, 4.2];

point_size = 24;
point_alpha = 0.72;

bar_width = 0.62;
bar_alpha = 0.88;

violin_width = 0.34;
violin_alpha = 0.22;
violin_point_jitter_width = 0.12;
median_line_half_width = 0.18;
mean_marker_size = 7;

% Legend locations are set separately because MATLAB's automatic 'best'
% location can place the legend on top of the bars in bar_points mode.
bar_points_legend_location = 'northeastoutside';
other_plot_legend_location = 'best';

%% ======================== VALIDATE SETTINGS ============================

model_mode = normalizeModelModeLocal(model_mode);
plot_style = lower(strtrim(char(plot_style)));

if ~isfolder(root_dir)
    error('root_dir does not exist: %s', root_dir);
end

if ~ismember(plot_style, {'bar_points', 'violin', 'points_median'})
    error(['plot_style must be ''bar_points'', ''violin'', or ' ...
           '''points_median''.']);
end

validatePlotStatisticOptionsLocal(show_mean, 'show_mean');
validatePlotStatisticOptionsLocal(show_median, 'show_median');

show_mean_current = show_mean.(plot_style);
show_median_current = show_median.(plot_style);

validateattributes(save_figure, {'logical', 'numeric'}, ...
    {'scalar'}, mfilename, 'save_figure');
validateattributes(save_mat, {'logical', 'numeric'}, ...
    {'scalar'}, mfilename, 'save_mat');

save_figure = logical(save_figure);
save_mat = logical(save_mat);

validateattributes(runIdx, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, mfilename, 'runIdx');

validateattributes(shared_varexp_threshold, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnegative'}, ...
    mfilename, 'shared_varexp_threshold');

validateattributes(group_colors, {'numeric'}, ...
    {'size', [2, 3], '>=', 0, '<=', 1}, mfilename, 'group_colors');

if ~(isnumeric(y_limits) && numel(y_limits) == 2 && ...
        all(isfinite(y_limits)) && y_limits(2) > y_limits(1))
    error('y_limits must be a finite increasing two-element vector.');
end
y_limits = reshape(y_limits, 1, 2);

latentSelectionTag = makeLatentSelectionTagLocal( ...
    use_dsl_filter, dsl_field, ...
    use_svexp_filter, svexp_field, shared_varexp_threshold);

model_tag = modelModeToFileTagLocal(model_mode);

if show_mean_current
    mean_tag = 'mean';
else
    mean_tag = 'nomean';
end

if show_median_current
    median_tag = 'median';
else
    median_tag = 'nomedian';
end

out_base = sprintf( ...
    '%s_%s_run%03d_subspace_similarity_%s_%s_%s_%s', ...
    data_content, model_tag, runIdx, latentSelectionTag, ...
    plot_style, mean_tag, median_tag);

fprintf('Root dir              : %s\n', root_dir);
fprintf('Model mode            : %s\n', model_mode);
fprintf('Latent selection tag  : %s\n', latentSelectionTag);
fprintf('Plot style             : %s\n', plot_style);
fprintf('Show mean              : %d\n', show_mean_current);
fprintf('Show median            : %d\n', show_median_current);
fprintf('Save figure            : %d\n', save_figure);
fprintf('Save mat               : %d\n', save_mat);
fprintf('Output base            : %s\n', out_base);

%% ======================= FIND SESSION FOLDERS ==========================

session_dirs = findCatgtSessionDirsLocal(root_dir);

if isempty(session_dirs)
    error('No catgt_* folders found one level below root_dir: %s', root_dir);
end

fprintf('Found %d candidate catgt_* session folders.\n', numel(session_dirs));

%% ======================== READ ALL SESSIONS ============================

category_keys = {
    'group1_across_captures_within', ...
    'group1_within_captures_across', ...
    'group1_feedforward_captures_feedback', ...
    'group1_feedback_captures_feedforward', ...
    'group2_across_captures_within', ...
    'group2_within_captures_across', ...
    'group2_feedforward_captures_feedback', ...
    'group2_feedback_captures_feedforward'
};

category_labels = {
    'Across captures Within', ...
    'Within captures Across', ...
    'FF captures FB', ...
    'FB captures FF', ...
    'Across captures Within', ...
    'Within captures Across', ...
    'FF captures FB', ...
    'FB captures FF'
};

category_group = [1, 1, 1, 1, 2, 2, 2, 2];

x_positions = [
    1, ...
    1 + within_group_spacing, ...
    1 + 2 * within_group_spacing, ...
    1 + 3 * within_group_spacing, ...
    1 + 3 * within_group_spacing + between_group_spacing, ...
    1 + 4 * within_group_spacing + between_group_spacing, ...
    1 + 5 * within_group_spacing + between_group_spacing, ...
    1 + 6 * within_group_spacing + between_group_spacing
];

num_categories = numel(category_keys);
values = nan(0, num_categories);

records = struct( ...
    'session_name', {}, ...
    'session_dir', {}, ...
    'input_file', {}, ...
    'condition_id', {}, ...
    'stim_abbrev', {}, ...
    'group1_across_captures_within', {}, ...
    'group1_within_captures_across', {}, ...
    'group1_feedforward_captures_feedback', {}, ...
    'group1_feedback_captures_feedforward', {}, ...
    'group2_across_captures_within', {}, ...
    'group2_within_captures_across', {}, ...
    'group2_feedforward_captures_feedback', {}, ...
    'group2_feedback_captures_feedforward', {});

skipped = struct('session_dir', {}, 'reason', {});

for si = 1:numel(session_dirs)

    session_dir = session_dirs{si};
    [~, session_name] = fileparts(session_dir);

    input_file = makeInputFileLocal( ...
        session_dir, data_content, model_mode, runIdx, latentSelectionTag);

    fprintf('[%d/%d] Reading %s\n', si, numel(session_dirs), session_dir);

    if ~isfile(input_file)
        reason = sprintf('Input file not found: %s', input_file);
        warning('%s', reason);
        skipped(end + 1) = struct( ...
            'session_dir', session_dir, 'reason', reason); %#ok<SAGROW>
        continue;
    end

    try
        S = load(input_file, 'SubspaceSim');

        if ~isfield(S, 'SubspaceSim')
            error('%s does not contain SubspaceSim.', input_file);
        end

        SubspaceSim = S.SubspaceSim;
        validateSubspaceSimLocal(SubspaceSim);

        if strcmp(model_mode, 'all_condition_model')

            row_values = readOnePooledModelLocal(SubspaceSim);

            rec = makeRecordLocal( ...
                session_name, session_dir, input_file, ...
                NaN, 'all', row_values);

            values(end + 1, :) = row_values; %#ok<SAGROW>
            records(end + 1) = rec; %#ok<SAGROW>

        else

            [session_values, condition_ids, stim_abbrev] = ...
                readOneConditionSummaryLocal(SubspaceSim);

            for ci = 1:size(session_values, 1)

                rec = makeRecordLocal( ...
                    session_name, session_dir, input_file, ...
                    condition_ids(ci), stim_abbrev{ci}, ...
                    session_values(ci, :));

                values(end + 1, :) = session_values(ci, :); %#ok<SAGROW>
                records(end + 1) = rec; %#ok<SAGROW>
            end
        end

    catch ME
        reason = ME.message;
        warning('Skipping %s:\n%s', session_dir, reason);
        skipped(end + 1) = struct( ...
            'session_dir', session_dir, 'reason', reason); %#ok<SAGROW>
    end
end

if isempty(records)
    error('No valid subspace-similarity records were collected.');
end

if ~any(isfinite(values(:)))
    error('All collected subspace-similarity values are NaN.');
end

%% ======================== BUILD SUMMARY ================================

values_by_category = cell(1, num_categories);
category_mean = nan(1, num_categories);
category_median = nan(1, num_categories);
category_n = zeros(1, num_categories);

for k = 1:num_categories
    v = values(:, k);
    v = v(isfinite(v));

    values_by_category{k} = v;
    category_n(k) = numel(v);

    if ~isempty(v)
        category_mean(k) = mean(v);
        category_median(k) = median(v);
    end
end

SubspaceSimilaritySummary = struct();

SubspaceSimilaritySummary.meta = struct();
SubspaceSimilaritySummary.meta.program = mfilename;
SubspaceSimilaritySummary.meta.root_dir = root_dir;
SubspaceSimilaritySummary.meta.data_content = data_content;
SubspaceSimilaritySummary.meta.model_mode = model_mode;
SubspaceSimilaritySummary.meta.runIdx = runIdx;
SubspaceSimilaritySummary.meta.use_dsl_filter = use_dsl_filter;
SubspaceSimilaritySummary.meta.dsl_field = dsl_field;
SubspaceSimilaritySummary.meta.use_svexp_filter = use_svexp_filter;
SubspaceSimilaritySummary.meta.svexp_field = svexp_field;
SubspaceSimilaritySummary.meta.shared_varexp_threshold = ...
    shared_varexp_threshold;
SubspaceSimilaritySummary.meta.latentSelectionTag = ...
    latentSelectionTag;
SubspaceSimilaritySummary.meta.plot_style = plot_style;
SubspaceSimilaritySummary.meta.show_mean = show_mean;
SubspaceSimilaritySummary.meta.show_median = show_median;
SubspaceSimilaritySummary.meta.show_mean_current = show_mean_current;
SubspaceSimilaritySummary.meta.show_median_current = show_median_current;
SubspaceSimilaritySummary.meta.save_figure = save_figure;
SubspaceSimilaritySummary.meta.save_mat = save_mat;

if strcmp(model_mode, 'condition_specific_models')
    SubspaceSimilaritySummary.meta.observation_unit = ...
        'one session-condition';
else
    SubspaceSimilaritySummary.meta.observation_unit = 'one session';
end

SubspaceSimilaritySummary.meta.num_session_folders_found = ...
    numel(session_dirs);
SubspaceSimilaritySummary.meta.num_session_folders_skipped = ...
    numel(skipped);
SubspaceSimilaritySummary.meta.num_observations = size(values, 1);

SubspaceSimilaritySummary.category_keys = category_keys;
SubspaceSimilaritySummary.category_labels = category_labels;
SubspaceSimilaritySummary.category_group = category_group;
SubspaceSimilaritySummary.x_positions = x_positions;

SubspaceSimilaritySummary.values = values;
SubspaceSimilaritySummary.values_by_category = values_by_category;
SubspaceSimilaritySummary.mean = category_mean;
SubspaceSimilaritySummary.median = category_median;
SubspaceSimilaritySummary.n = category_n;

SubspaceSimilaritySummary.records = records;
SubspaceSimilaritySummary.skipped = skipped;

mat_file = fullfile(root_dir, [out_base, '.mat']);
fig_file = fullfile(root_dir, [out_base, '.fig']);
svg_file = fullfile(root_dir, [out_base, '.svg']);

SubspaceSimilaritySummary.output_files = struct( ...
    'mat', mat_file, ...
    'fig', fig_file, ...
    'svg', svg_file);

if save_mat
    save(mat_file, 'SubspaceSimilaritySummary', '-v7.3');
    fprintf('Saved mat: %s\n', mat_file);
else
    fprintf('MAT saving is disabled.\n');
end

%% ============================== PLOT ===================================

fig = figure( ...
    'Visible', figure_visible, ...
    'Color', 'w', ...
    'Units', 'inches', ...
    'Position', [1, 1, figure_size_inches], ...
    'Renderer', 'painters');

ax = axes(fig);
hold(ax, 'on');

for k = 1:num_categories

    x = x_positions(k);
    g = category_group(k);
    this_color = group_colors(g, :);
    v = values_by_category{k};

    if isempty(v)
        continue;
    end

    switch plot_style

        case 'bar_points'

            if show_mean_current
                bar(ax, x, category_mean(k), bar_width, ...
                    'FaceColor', this_color, ...
                    'FaceAlpha', bar_alpha, ...
                    'EdgeColor', 'none', ...
                    'HandleVisibility', 'off');
            end

            drawPointsLocal( ...
                ax, repmat(x, size(v)), v, point_size, ...
                this_color, 'k', point_alpha);

            if show_median_current
                plot(ax, x, category_median(k), 'ks', ...
                    'MarkerSize', 7, ...
                    'MarkerFaceColor', 'k', ...
                    'HandleVisibility', 'off');
            end

        case {'violin', 'points_median'}

            if strcmp(plot_style, 'violin')
                [density, y_grid] = estimateDensityLocal(v, y_limits);

                if ~isempty(density)
                    density = density ./ max(density) .* violin_width;

                    patch(ax, ...
                        [x - density, fliplr(x + density)], ...
                        [y_grid, fliplr(y_grid)], ...
                        this_color, ...
                        'FaceAlpha', violin_alpha, ...
                        'EdgeColor', this_color, ...
                        'LineWidth', 1, ...
                        'HandleVisibility', 'off');
                end
            end

            x_jitter = deterministicJitterLocal( ...
                numel(v), violin_point_jitter_width, k);

            drawPointsLocal( ...
                ax, x + x_jitter, v, point_size, ...
                this_color, 'none', point_alpha);

            if show_mean_current
                plot(ax, x, category_mean(k), 'kd', ...
                    'MarkerSize', mean_marker_size, ...
                    'MarkerFaceColor', 'k', ...
                    'HandleVisibility', 'off');
            end

            if show_median_current
                plot(ax, ...
                    [x - median_line_half_width, ...
                     x + median_line_half_width], ...
                    [category_median(k), category_median(k)], ...
                    '-', ...
                    'Color', this_color, ...
                    'LineWidth', 2, ...
                    'HandleVisibility', 'off');
            end
    end
end

% Legend shows Group 1 and Group 2 colors only.
legend_handles = gobjects(2, 1);

for g = 1:2
    legend_handles(g) = plot(ax, nan, nan, 'o', ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', group_colors(g, :), ...
        'MarkerEdgeColor', group_colors(g, :), ...
        'LineStyle', 'none');
end

if strcmp(plot_style, 'bar_points')
    legend_location = bar_points_legend_location;
else
    legend_location = other_plot_legend_location;
end

legend(ax, legend_handles, {'Group 1', 'Group 2'}, ...
    'Location', legend_location, ...
    'Box', 'off', ...
    'FontName', 'Arial', ...
    'FontSize', 9);

x_margin = 0.55;
xlim(ax, [x_positions(1) - x_margin, x_positions(end) + x_margin]);
ylim(ax, y_limits);

set(ax, ...
    'XTick', x_positions, ...
    'XTickLabel', category_labels, ...
    'TickLabelInterpreter', 'none', ...
    'TickDir', 'out', ...
    'TickLength', [0.018, 0.018], ...
    'Box', 'off', ...
    'FontName', 'Arial', ...
    'FontSize', 9, ...
    'LineWidth', 0.9, ...
    'XColor', 'k', ...
    'YColor', 'k', ...
    'Layer', 'top');

xtickangle(ax, 25);
ylabel(ax, 'Subspace similarity', 'FontName', 'Arial');
grid(ax, 'off');

% No title: intended for direct poster/paper assembly.
set(findall(fig, '-property', 'FontName'), 'FontName', 'Arial');

if save_figure
    savefig(fig, fig_file);

    try
        exportgraphics(fig, svg_file, ...
            'ContentType', 'vector', ...
            'BackgroundColor', 'none');
    catch
        print(fig, svg_file, '-dsvg', '-painters');
    end

    fprintf('Saved fig: %s\n', fig_file);
    fprintf('Saved svg: %s\n', svg_file);
else
    fprintf('Figure saving is disabled.\n');
end
fprintf('N by category: [%s]\n', num2str(category_n));

if ~isempty(skipped)
    if save_mat
        fprintf(['Skipped %d session folder(s). Reasons are stored in ' ...
                 'the mat file.\n'], numel(skipped));
    else
        fprintf(['Skipped %d session folder(s). MAT saving is disabled, ' ...
                 'so inspect the warnings above for the reasons.\n'], ...
                numel(skipped));
    end
end

if close_after_save
    close(fig);
end

fprintf('Done.\n');

%% ======================================================================
% Local functions
% =======================================================================

function validatePlotStatisticOptionsLocal(options, variable_name)

    required_fields = {'bar_points', 'violin', 'points_median'};

    if ~(isstruct(options) && isscalar(options))
        error('%s must be a scalar struct.', variable_name);
    end

    for k = 1:numel(required_fields)
        field_name = required_fields{k};

        if ~isfield(options, field_name)
            error('%s is missing field %s.', variable_name, field_name);
        end

        value = options.(field_name);

        if ~((islogical(value) || isnumeric(value)) && isscalar(value) && ...
                isfinite(double(value)) && ismember(double(value), [0, 1]))
            error('%s.%s must be true or false.', ...
                variable_name, field_name);
        end

        options.(field_name) = logical(value); %#ok<NASGU>
    end
end

function model_mode = normalizeModelModeLocal(model_mode)

    model_mode = lower(strtrim(char(model_mode)));
    compact = regexprep(model_mode, '[^a-z0-9]', '');

    switch compact
        case {'allconditionmodel', 'allconditionsmodel', 'allcondition'}
            model_mode = 'all_condition_model';

        case {'conditionspecificmodels', 'conditionspecificmodel', ...
              'conditionmodels', 'conditionmode'}
            model_mode = 'condition_specific_models';

        otherwise
            error(['Unknown model_mode: %s. Use ''all_condition_model'' ' ...
                   'or ''condition_specific_models''.'], model_mode);
    end
end

function tag = modelModeToFileTagLocal(model_mode)

    switch normalizeModelModeLocal(model_mode)
        case 'all_condition_model'
            tag = 'all_condition_M';

        case 'condition_specific_models'
            tag = 'condition_specific_M';
    end
end

function tag = makeLatentSelectionTagLocal( ...
        use_dsl_filter, dsl_field, ...
        use_svexp_filter, svexp_field, shared_varexp_threshold)

    parts = {};

    if use_dsl_filter
        parts{end + 1} = sprintf( ...
            'DSL_%s_filtered', sanitizeTagLocal(dsl_field)); %#ok<AGROW>
    end

    if use_svexp_filter
        parts{end + 1} = sprintf( ...
            'SVExp_%s_%s_filtered', ...
            thresholdTagLocal(shared_varexp_threshold), ...
            sanitizeTagLocal(svexp_field)); %#ok<AGROW>
    end

    if isempty(parts)
        tag = 'all_latents';
    else
        tag = strjoin(parts, '_and_');
    end
end

function tag = sanitizeTagLocal(x)

    tag = char(string(x));
    tag = regexprep(tag, '[^A-Za-z0-9]+', '_');
    tag = regexprep(tag, '^_+|_+$', '');

    if isempty(tag)
        tag = 'field';
    end
end

function tag = thresholdTagLocal(threshold)

    tag = sprintf('threshold%.6g', threshold);
    tag = strrep(tag, '.', 'p');
    tag = strrep(tag, '-', 'm');
end

function session_dirs = findCatgtSessionDirsLocal(root_dir)

    listing = dir(fullfile(root_dir, '*', 'catgt_*'));
    session_dirs = {};

    for k = 1:numel(listing)
        if ~listing(k).isdir
            continue;
        end

        session_dirs{end + 1, 1} = fullfile( ...
            listing(k).folder, listing(k).name); %#ok<AGROW>
    end

    if isempty(session_dirs)
        return;
    end

    session_dirs = unique(session_dirs, 'stable');
    lower_dirs = cellfun(@lower, session_dirs, 'UniformOutput', false);
    [~, order] = sort(lower_dirs);
    session_dirs = session_dirs(order);
end

function input_file = makeInputFileLocal( ...
        session_dir, data_content, model_mode, runIdx, latentSelectionTag)

    if strcmp(model_mode, 'all_condition_model')
        input_file = fullfile( ...
            session_dir, ...
            ['FA_Dlag_', data_content], ...
            'mat_results', ...
            sprintf('run%03d', runIdx), ...
            sprintf('subspace_similarity_%s.mat', latentSelectionTag));
    else
        input_file = fullfile( ...
            session_dir, ...
            sprintf('%s_condition_mode_subspace_similarity_%s.mat', ...
            data_content, latentSelectionTag));
    end
end

function validateSubspaceSimLocal(SubspaceSim)

    if ~isstruct(SubspaceSim) || ~isfield(SubspaceSim, 'group')
        error('SubspaceSim must be a struct containing group.');
    end

    if numel(SubspaceSim.group) < 2
        error('SubspaceSim.group must contain at least two groups.');
    end

    required_pairs = {'across_vs_within', 'feedforward_vs_feedback'};

    for g = 1:2
        for p = 1:numel(required_pairs)
            getPairByNameLocal(SubspaceSim, g, required_pairs{p});
        end
    end
end

function row_values = readOnePooledModelLocal(SubspaceSim)

    specs = similaritySpecsLocal();
    row_values = nan(1, size(specs, 1));

    for k = 1:size(specs, 1)
        row_values(k) = getScalarSimilarityLocal( ...
            SubspaceSim, specs{k, 1}, specs{k, 2}, specs{k, 3});
    end
end

function [values, condition_ids, stim_abbrev] = ...
        readOneConditionSummaryLocal(SubspaceSim)

    ref_pair = getPairByNameLocal( ...
        SubspaceSim, 1, 'across_vs_within');

    if ~isfield(ref_pair, 'condition_id')
        error('Condition-summary pair is missing condition_id.');
    end

    if ~isfield(ref_pair, 'stim_abbrev')
        error('Condition-summary pair is missing stim_abbrev.');
    end

    condition_ids = reshape(ref_pair.condition_id, [], 1);
    stim_abbrev = normalizeLabelCellLocal(ref_pair.stim_abbrev);

    nCond = numel(condition_ids);

    if numel(stim_abbrev) ~= nCond
        error('condition_id and stim_abbrev have different lengths.');
    end

    specs = similaritySpecsLocal();
    values = nan(nCond, size(specs, 1));

    for k = 1:size(specs, 1)
        g = specs{k, 1};
        pair_name = specs{k, 2};
        similarity_field = specs{k, 3};
        pair_result = getPairByNameLocal(SubspaceSim, g, pair_name);

        if ~isfield(pair_result, 'similarity') || ...
                ~isfield(pair_result.similarity, similarity_field)
            error('Group %d pair %s is missing similarity.%s.', ...
                g, pair_name, similarity_field);
        end

        v = pair_result.similarity.(similarity_field);

        if ~(isnumeric(v) && isreal(v))
            error('Group %d pair %s similarity.%s must be numeric.', ...
                g, pair_name, similarity_field);
        end

        v = reshape(double(v), [], 1);

        if numel(v) ~= nCond
            error(['Group %d pair %s similarity.%s has %d values; ' ...
                   'expected %d conditions.'], ...
                  g, pair_name, similarity_field, numel(v), nCond);
        end

        values(:, k) = v;
    end
end

function value = getScalarSimilarityLocal( ...
        SubspaceSim, group_idx, pair_name, similarity_field)

    pair_result = getPairByNameLocal( ...
        SubspaceSim, group_idx, pair_name);

    if ~isfield(pair_result, 'similarity') || ...
            ~isfield(pair_result.similarity, similarity_field)
        error('Group %d pair %s is missing similarity.%s.', ...
            group_idx, pair_name, similarity_field);
    end

    value = pair_result.similarity.(similarity_field);

    if ~(isnumeric(value) && isscalar(value) && isreal(value))
        error('Group %d pair %s similarity.%s must be a numeric scalar.', ...
            group_idx, pair_name, similarity_field);
    end

    value = double(value);
end

function specs = similaritySpecsLocal()

    specs = {
        1, 'across_vs_within', ...
           'across_captures_within';
        1, 'across_vs_within', ...
           'within_captures_across';
        1, 'feedforward_vs_feedback', ...
           'feedforward_captures_feedback';
        1, 'feedforward_vs_feedback', ...
           'feedback_captures_feedforward';
        2, 'across_vs_within', ...
           'across_captures_within';
        2, 'across_vs_within', ...
           'within_captures_across';
        2, 'feedforward_vs_feedback', ...
           'feedforward_captures_feedback';
        2, 'feedforward_vs_feedback', ...
           'feedback_captures_feedforward'
    };
end

function pair_result = getPairByNameLocal( ...
        SubspaceSim, group_idx, requested_pair_name)

    group_result = SubspaceSim.group(group_idx);

    if ~isfield(group_result, 'pair') || ~iscell(group_result.pair)
        error('SubspaceSim.group(%d).pair must be a cell array.', group_idx);
    end

    pair_idx = [];

    if isfield(group_result, 'pairNames')
        pair_names = cellstr(string(group_result.pairNames));
        pair_idx = find(strcmp(pair_names, requested_pair_name), 1);
    end

    if isempty(pair_idx)
        for p = 1:numel(group_result.pair)
            this_pair = group_result.pair{p};

            if isstruct(this_pair) && isfield(this_pair, 'name') && ...
                    strcmp(char(string(this_pair.name)), requested_pair_name)
                pair_idx = p;
                break;
            end
        end
    end

    if isempty(pair_idx)
        error('Group %d is missing pair %s.', ...
            group_idx, requested_pair_name);
    end

    if pair_idx > numel(group_result.pair)
        error('Group %d pair index %d exceeds the pair cell array.', ...
            group_idx, pair_idx);
    end

    pair_result = group_result.pair{pair_idx};

    if ~isstruct(pair_result)
        error('SubspaceSim.group(%d).pair{%d} must be a struct.', ...
            group_idx, pair_idx);
    end
end

function labels = normalizeLabelCellLocal(labels_in)

    if isstring(labels_in)
        labels = cellstr(labels_in(:));
    elseif ischar(labels_in)
        labels = cellstr(labels_in);
    elseif iscell(labels_in)
        labels = cell(numel(labels_in), 1);

        for k = 1:numel(labels_in)
            labels{k} = char(string(labels_in{k}));
        end
    else
        error('stim_abbrev must be a char array, string array, or cell array.');
    end

    labels = reshape(labels, [], 1);
end

function rec = makeRecordLocal( ...
        session_name, session_dir, input_file, ...
        condition_id, stim_abbrev, row_values)

    rec = struct();
    rec.session_name = session_name;
    rec.session_dir = session_dir;
    rec.input_file = input_file;
    rec.condition_id = condition_id;
    rec.stim_abbrev = stim_abbrev;
    rec.group1_across_captures_within = row_values(1);
    rec.group1_within_captures_across = row_values(2);
    rec.group1_feedforward_captures_feedback = row_values(3);
    rec.group1_feedback_captures_feedforward = row_values(4);
    rec.group2_across_captures_within = row_values(5);
    rec.group2_within_captures_across = row_values(6);
    rec.group2_feedforward_captures_feedback = row_values(7);
    rec.group2_feedback_captures_feedforward = row_values(8);
end

function drawPointsLocal( ...
        ax, x, y, point_size, face_color, edge_color, point_alpha)

    try
        scatter(ax, x, y, point_size, ...
            'o', ...
            'MarkerFaceColor', face_color, ...
            'MarkerEdgeColor', edge_color, ...
            'MarkerFaceAlpha', point_alpha, ...
            'MarkerEdgeAlpha', point_alpha, ...
            'LineWidth', 0.45, ...
            'HandleVisibility', 'off');
    catch
        scatter(ax, x, y, point_size, ...
            'o', ...
            'MarkerFaceColor', face_color, ...
            'MarkerEdgeColor', edge_color, ...
            'LineWidth', 0.45, ...
            'HandleVisibility', 'off');
    end
end

function jitter = deterministicJitterLocal(n, half_width, seed_offset)

    if n <= 1 || half_width <= 0
        jitter = zeros(n, 1);
        return;
    end

    stream = RandStream('mt19937ar', 'Seed', 1000 + seed_offset);
    jitter = (2 .* rand(stream, n, 1) - 1) .* half_width;
end

function [density, y_grid] = estimateDensityLocal(values, y_limits)

    values = values(:);
    values = values(isfinite(values));

    density = [];
    y_grid = [];

    if numel(values) < 2 || max(values) == min(values)
        return;
    end

    if exist('ksdensity', 'file') == 2
        try
            [density, y_grid] = ksdensity(values, 'NumPoints', 100);
        catch
            density = [];
            y_grid = [];
        end
    end

    if isempty(density)
        numBins = min(20, max(5, round(sqrt(numel(values)))));
        [counts, edges] = histcounts(values, numBins, ...
            'Normalization', 'pdf');
        y_grid = edges(1:end-1) + diff(edges) ./ 2;
        density = counts;
    end

    keep = isfinite(density) & isfinite(y_grid) & ...
        y_grid >= y_limits(1) & y_grid <= y_limits(2);

    density = density(keep);
    y_grid = y_grid(keep);

    if isempty(density) || max(density) <= 0
        density = [];
        y_grid = [];
        return;
    end

    density = reshape(density, 1, []);
    y_grid = reshape(y_grid, 1, []);
end
