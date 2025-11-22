%Axon analysis and quantification code associated with Gutman-Wei et al.,
%2025 - function for analysis of horizontal segment and layer-specific data
%Function for computing descriptive data on a subset of a cell tracing, 
%to be applied to vertically-segmented portions of an axon tracing
%Calculates total length, branch points, end points contained within the
%data passed in. Calculates 2-d convex hull in the x-z plane for points
%passed in. Plots overlaid tracings of all cells within group for each
%group, viewed as a projection onto x-z plane. Makes plots and statistical
%comparisons of all metrics. Operates on a cell array of tracings, multiple
%tracings can be processed at once.

%Inputs: 
% allCells: cell array of data from each cell, including points, lengths, 
% point type, etc. Each element is a structure corresponding to one
% reconstruction
% subset: cell array of logical vectors, with each element being a cell
% (indexed same as allCells). Each element of logical vector denotes
% whether the corresponding point is to be included in the subset or not
% group: defines the group that each cell belongs to (usually age or Kir2.1
% status). 
% ageGroups: names of the groups, number of elements matches number of
% unique values in group
% figureTitle: string describing what subset is being plotted, will be
% inserted in all figure titles

%Outputs:
% horzData: structure with fields that contain numeric arrays summarizing 
% data from all cells. Totals of length, branch points, end points, convex
% hull are reported, as aare the groups associated with each cell.
% plots: plots of total length, branch points, end points, and horizontal
% (x-z) convex hulls of the subset are plotted across groups with
% statistical comparisons. In addition, two views of the subsets en face
% are plotted. One shows all axon points of overlaid cells, the other shows
% only branch points and end points of overlaid cells. Number of cells
% overlaid is determined by the group with the smallest number of cells
% (N). The N cells closest to the mean subset length in the other groups
% are plotted along with all cellss from the smallest groups.


function [horzData] = HorizontalSegmentAnalysis(allCells, subset, group, ageGroups, figureTitle)
%% Loop through and find all metrics
subsetData = {nan(size(allCells)), nan(size(allCells)),...
    nan(size(allCells)), nan(size(allCells))};
for i = 1:length(allCells)
    subsetData{1}(i) = sum(allCells{i}.lengths(subset{i}));
    subsetData{2}(i) = sum(allCells{i}.pointType(subset{i}, 1));
    subsetData{3}(i) = sum(allCells{i}.pointType(subset{i}, 2));
    if nnz(allCells{i}.points(subset{i}, 1)) > 2
        [~, subsetData{4}(i)] = convhull(allCells{i}.points(subset{i}, 1), allCells{i}.points(subset{i},3)); 
    else
        subsetData{4}(i) = 0;
    end
end

%% Make structure
horzData = struct('length', subsetData{1}, 'branchpoints',...
    subsetData{2}, 'endpoints', subsetData{3}, 'convexHull', subsetData{4},...
    'group', group);

%% Plot summary data by group
figure('Name', figureTitle)
ylabels = {'Length (um)', 'Number of branch points', ...
    'Number of endpoints', 'Convex hull area (um^2)'};
for i = 1:length(subsetData)
    subplot(1,4,i)
    plotBetweenGroups(subsetData{i}, group, ageGroups)
    ylabel(ylabels{i})
    figQuality(gcf, gca, [12 4])
end

%% Even out numbers for subset data if groups have uneven number
%Comment out if plotting projections of all cells. Takes the group with the
%lowest number of cells (N), selects the N cells closest to the mean for
%display from the rest of the groups.
%Find group counts
groupInfo = [unique(group), zeros(size(unique(group)))]; % group index, count
for i = 1:size(groupInfo, 1)
    groupInfo(i, 2) = nnz(group == groupInfo(i,1));
end
groupN = min(groupInfo(:,2)); %number of cells to display for each group

%Alternatively, can manually override and specify 10:
%groupN = 10;
%disp('10 cells used for en face views - manually specified');

cellsToPlot = false(size(allCells));

%Filter by lengths closest to the group mean
for i = 1:length(unique(group))
    tempGroupIndex = find(group == i);
    tempDev = subsetData{1}(tempGroupIndex) - mean(subsetData{1}(tempGroupIndex)); %deviation from mean of length only!
    [~, keepIndex] = mink(abs(tempDev), groupN); %Find the minimum N cells to keep, closest to mean
    cellsToPlot(tempGroupIndex(keepIndex)) = true; %Make logical indexing array of cells to discard
end

%% Make x-z projection plots 
%Plot figure of axon overlays separated by group
figure('Name', ['Axon overlays ' figureTitle])
%count points per group
for i = 1:length(allCells)
    if cellsToPlot(i)
        subplot(1, numel(unique(group)), group(i));
        hold on
        scatter(allCells{i}.points(subset{i}, 1), allCells{i}.points(subset{i}, 3), 7, [0.5 0.5 0.5], 'filled');
    end
end

%Find axes limits
xlimits = [0 0]; %Initialize x-limits for plots
zlimits = [0 0]; %Initialize z-limits for plots
for i = 1:numel(ageGroups)
    subplot(1, numel(ageGroups), i)
    tempXLim = xlim;
    tempZLim = ylim;
    
    if tempXLim(1) < xlimits(1)
        xlimits(1) = tempXLim(1);
    end
    if tempXLim(2) > xlimits(2)
        xlimits(2) = tempXLim(2);
    end
    
    if tempZLim(1) < zlimits(1)
        zlimits(1) = tempZLim(1);
    end
    if tempZLim(2) > zlimits(2)
        zlimits(2) = tempZLim(2);
    end
end

%unify x and z limits to make it square
if zlimits(2) > xlimits(2)
    xlimits(2) = zlimits(2);
end
if zlimits(1) < xlimits(1)
    xlimits(1) = zlimits(1);
end

%Apply limits and formatting to axon plot
for i = 1:numel(ageGroups)
    subplot(1, numel(ageGroups), i)
    xlim(xlimits)
    ylim(xlimits)
    xlabel('Tangential axis 1 (um)')
    ylabel('Tangential axis 2 (um)')
    title([ageGroups{i}, ' (' num2str(nnz(group(cellsToPlot) == i)), ' cells)']);
    figQuality(gcf, gca, [numel(ageGroups) * 6 5])
end

%Plot branch and end points only
figure('Name', ['Branch and end points ' figureTitle])
for i = 1:length(allCells)
    if cellsToPlot(i)
        subplot(1, numel(unique(group)), group(i));
        hold on
        tempBranch = subset{i} & allCells{i}.pointType(:, 1); 
        scatter(allCells{i}.points(tempBranch, 1), allCells{i}.points(tempBranch, 3), 20, 'g', 'filled');
        tempEnd = subset{i} & allCells{i}.pointType(:, 2);  
        scatter(allCells{i}.points(tempEnd, 1), allCells{i}.points(tempEnd, 3), 20, 'm', 'filled');
    end
end

%Apply limits and formatting to branch/end point plot
for i = 1:numel(ageGroups)
    subplot(1, numel(ageGroups), i)
    xlim(xlimits)
    ylim(xlimits)
    xlabel('Tangential axis 1 (um)')
    ylabel('Tangential axis 2 (um)')
    title([ageGroups{i}, ' (' num2str(nnz(group(cellsToPlot) == i)), ' cells)']);
    figQuality(gcf, gca, [numel(ageGroups) * 6 5])
end

end

%% Function for plotting
function plotBetweenGroups(data, group, labels)
%makes plot but need to separately specify the subplot used, run
%fiqQuality afterwards, and add titles/axis legends as desired

if numel(unique(group)) > 2 %K-W test followed by pairwise comparisons of consecutive groups
    %nonparametric 1-factor test across ages
    [pAll,~, ~] = kruskalwallis(data, group, 'off');
    pPairwise = nan(length(unique(group)) - 1, 1);
    %pairwise rank sum with adjacent ages
    for i = 1:length(pPairwise)
        pPairwise(i) = ranksum(data(group == i), data(group == i+1));
    end
elseif numel(unique(group)) == 2 %Pairwise M-W test only
    pPairwise = ranksum(data(group == 1), data(group == 2));
elseif numel(unique(group)) < 2
    disp('not enough groups')
    return
end

hold on
for j = 1:length(unique(group))
    tempData = data(group == j);
    scatter(group(group == j), tempData, 20, [0.5 0.5 0.5], 'filled', 'jitter', 'on')
    errorbar(j, mean(tempData), std(tempData)/sqrt(numel(tempData)), ...
        'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'CapSize', 0,...
        'MarkerSize', 8); 
end
xlim([0,(length(labels) + 1)]);
xticks(1:length(labels));
xticklabels(labels);
ylimits = ylim;
ylimits = [ylimits(1), ylimits(2)*1.25];
ylim(ylimits)
for i = 1:length(pPairwise)
    if rem(i, 2) == 0
        yLoc = ylimits(2) * 0.8;
    else
        yLoc = ylimits(2) * 0.85;
    end
    plot((i:i+1), yLoc * [1 1], '-k', 'LineWidth', 1)
    text(i + 0.1, yLoc * 1.03, ['p = ' num2str(round(pPairwise(i), 3))]);
end
    
if numel(unique(group)) > 2 %put overall K-W p-value as x-axis label if more than 2 groups
    xlabel(['KW p: ' num2str(round(pAll, 4))])
    %disp(pAll)
end
%disp(pPairwise)
    
end
