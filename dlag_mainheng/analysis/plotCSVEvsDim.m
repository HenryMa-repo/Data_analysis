
function plotCSVEvsDim(csve, xDims, optimalD)
%
% plotCSVEvsDim(csve, xDims, optimalD)
%
% Description:
%     Plot cumulative % shared variance explained versus dimensionality.
%     The plotting style follows plotPerfvsDim_fa.
%
% Arguments:
%
%     Required:
%
%     csve     -- cell array; each cell contains a vector of
%                 cumulative % shared variance explained values
%
%     xDims    -- cell array; same size as csve, each cell contains
%                 the corresponding dimensionality values
%
%     optimalD -- vector; length must equal numel(csve), indicating
%                 the shared dimensionality to highlight in each plot
%
% Outputs:
%     None. (But creates figures)
%

numGroups = numel(csve);

% basic checks
if ~iscell(csve) || ~iscell(xDims)
    error('csve and xDims must both be cell arrays.');
end

if numel(csve) ~= numel(xDims)
    error('csve and xDims must have the same number of cells.');
end

if numel(optimalD) ~= numGroups
    error('Length of optimalD must equal the number of cells in csve/xDims.');
end

figure;

for groupIdx = 1:numGroups
    y = csve{groupIdx};
    x = xDims{groupIdx};

    if numel(x) ~= numel(y)
        error('xDims{%d} and csve{%d} must have the same length.', groupIdx, groupIdx);
    end

    colors = generateColors(); % same style as original code

    subplot(1, numGroups, groupIdx);
    hold on;

    xlabel('dimensionality');
    ylabel('cumulative % shared variance explained for matrix LL^T');

    % main line
    plot(x, y, 'o-', ...
        'color', colors.grays{1}, ...
        'MarkerFaceColor', colors.grays{1}, ...
        'linewidth', 1.5);

    % find and mark optimalD
    idx = find(x == optimalD(groupIdx), 1);

    if ~isempty(idx)
        legendEntries = plot(x(idx), y(idx), 'p', ...
            'color', colors.reds{4}, ...
            'markerfacecolor', colors.reds{4}, ...
            'markersize', 10);

        legendLabels = 'shared dimensionality';
        legend(legendEntries, legendLabels, 'Location', 'southeast');
    else
        warning('optimalD(%d) = %g is not found in xDims{%d}.', ...
            groupIdx, optimalD(groupIdx), groupIdx);
    end

    hold off;
end