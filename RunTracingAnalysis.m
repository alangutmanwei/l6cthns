%Axon analysis and quantification code associated with Gutman-Wei et al.,
%2025 - main script
%Originally written in MATLAB 2019a

%Updated version of tracing analysis, which can be used for all tracing
%comparisons of axons and dendrites. This main script reads tracing files
%formatted in SWC file format and saved as txt files with an additional
%column at rightmost that specifies the layer location of each point. 
% Requires a metadata spreadsheet that contains the following:
% 1: File names of each tracing ("filename_txt")
% 2: Experimental group of each tracing ("Condition")
% 3: x, y, z coordinates, in the same frame of reference of the swc tracing
% file, of the white matter and pia intersections of a radially oriented
% vector that passes through the cell soma. This vector will be used to
% rotate tracings such that the y-axes of tracings are aligned with the
% wm-pia axis. The locations of these are defined manually on images.


%Processes axon tracing by removing short branches, removing the white
%matter axon, and rotating the reconstruction to be radially oriented.
%Calculates point locations relative to a scaled grid originating at the
%cell soma, and also calculates scaled heights of points within the 
%reconstruction. Compares length, branch points, end points,
%and volume of cortex filled between groups. Allows for segmenting
%horizontal portions of the cell (by depth or by layer) and comparing the
%subsets between groups using the same metrics. Includes functions at the
%end of the script to compare data between groups, to construct average
%heatmaps for each group/condition, and to plot depth profiles of axon
%length and convex hull area

%User input determines if script performs statistical comparisons between
%groups, or if it displays and stores visualizations of the cells and ends
%before performing statistics (user will be prompted for directory)

%Parameters to set:
%gridProp: size of scaled grid to use for binning axons to quantify axon
%distribution across cortex. As a proportion of cortical depth
%minLength: minimum length, in microns of branches to be considered
%gridType: string input, 'wm' if grid starts at the wm and goes to pia, or
%'cell' if the grid starts at 0 at the soma and goes as far as cell extends

%Calls functions stored in separate files:
%figQuality.m - spruces up figure to preferred style. Only works for 2D
%ProcessReconstruction.m - trims short branches below specified threshold
%(10um), trims off white matter processes, rotates the tracing to
%be aligned to a wm-pia vector specified by metadata file
%GenerateGridData.m - Bins axon length, branch points, and end points into
%grids based on user-specified grid size. Returns data structure that can
%be used to make scaled depth-distribution plots and scaled heatmaps
%HorizontalSegmentAnalysis.m - Takes a horizontal slice of the tracing,
%specified using a logical index array of included vs. excluded points of 
%the tracing, and computes characteristics of that slice. E.g. axon in top
%or bottom halves of cortex, or in layer X (if file has layer annotations)

%Inputs, set differently for each set of tracings:
% Note: relative paths are used for 'metadata' and 'path'
% metadata file: spreadsheet with file names and pia and wm locations, and
% the age/condition (e.g. 'P4', or 'Kir'). should include columns called
% piaX, piaY, piaZ, wmX, wmY, wmZ, and condition
% path: location of tracing swc files
% ageGroups: Name of each group or condition, to use for labeling plots
% compLabel: x-axis title, describing the groups being compared
% nameVar: variable name in the metadata table to use for the filename
% extension: whatever needs to be appended to the end of the filename.
% usually '.txt', '.swc', or '_Annotated.txt' for annotated tracings
% Tracings: located in path specified, must be formatted in swc pattern (7
% columns), with first column as index, last column as parent index, and
% columns 3-5 as XYZ values for points. Layer-annotated files are
% similarly formatted, but with an eighth column denoting layers as such:
% 2-> L1; 3-> L2/3; 4-> L4; 5-> L5; 6->L6 and 7-> WM. SWC headers are not
% necessarily, and headers with discontinuous lines tend to generate errors
% - remove header if in doubt.

%Outputs: 
% If user specifies to generate figures and not perform statistics, script
% will plot the tracings of axons/dendrites in 3D, and generate PNG files
% of 2D projections in XY and ZY planes (where Y is the WM-Pia axis).
% Figure characteristics are optimized to be generated when connected to a
% large desktop monitor, and are sized for printing
% Otherwise
% First plot is of total axon length and total branch number. An axon
% branch is defined as any segment of axon starting with soma or
% branchpoint, ending with a branchpoint or endpoint, with no branchpoints
% or endpoints in between.
% XY, ZY, and XZ projected heatmaps of neuronal process distribution across
% cortex are plotted, scaled to cortical size changes between conditions or
% across development. Uses integrated function display2DHeatmaps
% Curves of axon length and horizontal convex hull area are plotted across
% depth bins, scaled to cortical thickness. Uses integrated function 
% plotTracingsByDepth. Includes both average and individual curve versions
% Plots maximum height in cortex across conditions, at 100th, 90th, and
% 80th percentile of tracing.

% Queries user to stop code if desired before plotting subset of tracing
% Subsets are: upper half, lower half, lower 30%, layer 6, layer 5, and
% layer 4. Plots summary data of length, branch points,
% end points, and convex hull area contained in subset, then plots
% horizontal en face views (XZ plane) of the subset of axons. 2 versions of
% en face views shown, 1 with only axon length and the other with only axon
% branch points and end points. Same number of cells is plotted for each
% group for en face views, with the cells closest to the mean length
% plotted for each group if numbers are uneven.

%All summary data plots between conditions use the integrated function 
%plotBetweenGroups

clear
close all

%% Define files
metadata = readtable('Location/MetadataFileName.xlsx');
path = 'Location/';
disp('Groups 1 and 2 Comparison:')
ageGroups = {'Group1' 'Group2'}; 
compLabel = 'VariableCompared'; 
nameVar = 'filename_txt'; 
extension = '.txt';

%% Set parameters and find groups
gridProp = 0.05; %size of grid, defined as proportion of pia-WM distance
minLength = 10; %minimum length, in microns of branches to be considered
gridType = 'wm'; %WM-pia ('wm') or cell-centered grid ('cell')

group = zeros(height(metadata), 1); %group (age/condition) the cell belongs to
nCells = zeros(size(ageGroups)); %number of cells total
%count number of cells per age group
for i = 1:length(ageGroups)
   for j = 1:height(metadata)
      if contains(metadata.condition{j}, ageGroups{i})
          nCells(i) = nCells(i) + 1;
          group(j) = i;
      end
   end
end

%% Get cell data
allCells = cell(sum(nCells), 1); %data structures holding processed reconstruction
totalStats = nan(sum(nCells), 2); %length, number of branches

%Ask user if only doing example tracing figures and not analysis
outputFigures = logical(input('Only save figures and skip analysis?'));

%If so, set the file saving directory
if outputFigures
    savePath = uigetdir; %Initialize in current working directory
end

%Loop through all cells, process reconstructions and get summary data, save
%tracing to directory if commanded by user
for i = 1:length(allCells)
    tempName = [path metadata{i, nameVar}{1}, extension]; %add suffixes to file name
    disp(metadata{i, nameVar}{1}) %Print cell being processed
    [allCells{i}] = ProcessReconstruction(tempName, ...
        [metadata.piaX(i), metadata.piaY(i), metadata.piaZ(i)], ...
        [metadata.wmX(i), metadata.wmY(i), metadata.wmZ(i)], minLength, ...
        gridProp, outputFigures); 
    if outputFigures %Below for formatting figures to save
        title([metadata.condition{i}, ': ', metadata.filename_txt{i}, ' (',... 
            num2str(round(metadata.depth(i), 1)), '%)'], 'Interpreter', 'None');
        set(gcf, 'Position', [100 100 850 1100])
        xlimits = xlim;
        ylimits = ylim;
        totalInfo = text(xlimits(1)+20, ylimits(1)+20, ['length(um): ',...
            num2str(allCells{i}.totLength), ', branches: ',...
            num2str(allCells{i}.totBranch + allCells{i}.totEnd)]);
        axis equal
        %Optional comment zone to generate figures without saving: 
        %
        saveas(gcf, [savePath, '/', ageGroups{group(i)}, '_',...
            metadata.filename_txt{i}, '_XY']);
        print([savePath, '/', ageGroups{group(i)}, '_',...
            metadata.filename_txt{i}, '_XY'], ...
            '-dpng')
        view([1,0,0])
        delete(totalInfo)
        axis equal
        set(gcf, 'Position', [100 100 1100 850])
        set(gcf, 'PaperOrientation', 'landscape')
        print([savePath, '/', ageGroups{group(i)}, '_',...
            metadata.filename_txt{i}, '_ZY'], ...
            '-dpng')
        saveas(gcf, [savePath, '/', ageGroups{group(i)}, '_',...
            metadata.filename_txt{i}, '_ZY']);
        close all
        %}
    end
    totalStats(i,1) = allCells{i}.totLength; %transfer total length to array
    %Calculate branch number as branch points + end points
    totalStats(i,2) = allCells{i}.totBranch + allCells{i}.totEnd; 
end

%% Break if only saving figures
if outputFigures, return; end

%% Define WM-Pia reference grid and calculate gridData values
%Generate a grid starting at the white matter boundary, spaced using
%gridProp intervals 
wmGrid = 0:gridProp:1; %Cuts out points above and below the wm-pia line specified!
wmGridIndex = cell(size(allCells));
% clasify each point according to position relative to white matter. Points
% below WM remain NaN
for i = 1:length(wmGridIndex)
    wmGridIndex{i} = nan(size(allCells{i}.pointHeight));
    for j = 1:length(wmGrid)-1
        tempIndex = (allCells{i}.pointHeight >= wmGrid(j)) &...
            (allCells{i}.pointHeight < wmGrid(j+1));
        wmGridIndex{i}(tempIndex) = j - 0.5;
    end
end

%Bin tracings by grid using external function. Note that the wmGrid is
%passed in but NOT necessarily used, depending on gridType (only if 'wm'). 
gridData = cell(size(allCells));
for i = 1:length(wmGridIndex)
    [gridData{i}] = GenerateGridData(allCells{i}.pointGridIndex,...
        allCells{i}.lengths, allCells{i}.pointType,...
        allCells{i}.points, wmGridIndex{i}, gridType);
end

%% Plot totals
%First main plot - compare total statistics across conditions. 
figure %2 subplots: length, branches
ylabels = {'Length (um)', 'Number of branches'};
for i = 1:size(totalStats, 2)
    subplot(1,4,i)
    plotBetweenGroups(totalStats(:,i), group, ageGroups)
    ylimTemp = ylim;
    if i == 3
        ylim([1 ylimTemp(2)]); %Modify tortuosity lower bound to be 1
    end
    ylabel(ylabels{i})
    figQuality(gcf, gca, [6 4])
end

%% Normalize grid data outputs for plotting in heatmaps
% normalize all distributions by total length, branchpoints, endpoints
gridLengthNorm = cell(size(allCells));
gridBranchNorm = cell(size(allCells));
gridEndNorm = cell(size(allCells));
%normalize all arrays by the corresponding total value for the cell
for i = 1:length(gridLengthNorm) %Assume all measurement cell arrays are same dimensions
    gridLengthNorm{i} = gridData{i}.gridLength/allCells{i}.totLength;
    gridBranchNorm{i} = gridData{i}.gridBranch/allCells{i}.totBranch;
    gridEndNorm{i} = gridData{i}.gridEnd/allCells{i}.totEnd;
end

%% Plot heatmaps 
%Plot average heatmaps
ranges = display2DHeatmaps(gridData, gridLengthNorm, gridBranchNorm, gridEndNorm, group, ageGroups, 0.01);

%% Plot average curves of all cells across scaled depth
if contains(gridType, 'wm')
    plotTracingsByDepth(allCells, gridData, wmGridIndex,...
        [ranges{2}(1), round(1/gridProp - 0.5, 1)], group, ageGroups)
    for i = 1:2
        subplot(1,2,i)
        ylabel('Normalized distance from WM');
    end
elseif contains(gridType, 'cell')
    cellGridIndex = cell(size(allCells));
    for i = 1:length(allCells)
        cellGridIndex{i} = allCells{i}.pointGridIndex(:,2);
    end
    plotTracingsByDepth(allCells, gridData, cellGridIndex,...
        ranges{2}, group, ageGroups)

    for i = 1:2
        subplot(1,2,i)
        ylabel('Number of grids from cell body location');
    end
end

%% Ask user whether to continue to plotting subsets
stopEarly = logical(input('Do you want to stop early?'));
if stopEarly, return; end
% Following sections compare horizontal segments of tracing, including all
% data points within a certain vertical percentage of the cortex, or all
% data points within a certain layer (defined by layer annotations in the
% original tracing files - see Python code)

%% Horizontally subset based on cortical thickness proportions:
% Define upper half subset
upperHalfSubset = cell(size(allCells));
for i = 1:length(allCells)
    upperHalfSubset{i} = allCells{i}.pointHeight >= 0.5; %inclusive
end
HorizontalSegmentAnalysis(allCells, upperHalfSubset, group, ageGroups, '(upper half)');
% Define lower half subset
upperHalfSubset = cell(size(allCells));
for i = 1:length(allCells)
    upperHalfSubset{i} = allCells{i}.pointHeight < 0.5; %exclusive
end
HorizontalSegmentAnalysis(allCells, upperHalfSubset, group, ageGroups, '(lower half)');

%% Horizontally segment based on layer annotations:
% Layer 6
l6Subset = cell(size(allCells));
for i = 1:length(allCells)
    l6Subset{i} = allCells{i}.layerID >= 6; 
end
l6SegmentStats = HorizontalSegmentAnalysis(allCells, l6Subset, group, ageGroups, '(L6)');
% Layer 5
l5Subset = cell(size(allCells));
for i = 1:length(allCells)
    l5Subset{i} = allCells{i}.layerID == 5;
end
l5SegmentStats = HorizontalSegmentAnalysis(allCells, l5Subset, group, ageGroups, '(L5)');
% Layer 4
l4Subset = cell(size(allCells));
for i = 1:length(allCells)
    l4Subset{i} = allCells{i}.layerID == 4;
end
l4SegmentStats = HorizontalSegmentAnalysis(allCells, l4Subset, group, ageGroups, '(L4)');


%% Function to plot all summary data between groups
% Function compares summary data between groups. 
%makes plot but need to separately specify the subplot used, run
%fiqQuality afterwards, and add titles/axis legends as desired

%inputs:
%data: summary data to be plotted as a single vector
%group: integers assigning data to different groups
%labels: labels for each group, numbered in order of the group index 

%outputs (no variables):
%plot of the data across groups, in ascending order of group index, with
%individual points, means with SEM, and non-parametric statistical tests
%between groups. 

function plotBetweenGroups(data, group, labels)

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
ylimits = [0, ylimits(2)*1.25];
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
end
    
%Whether to display p-values
%
if exist('pAll', 'var'), disp(pAll); end
disp(pPairwise)
%}
end

%% Function to generate average heatmaps
%Plots average heatmaps of 2-D projections of cells, in X-Y, Y-Z, and X-Z
%planes where Y is the wm-pia axis. Uses black-to-white heatmaps, and
%boosts any non-zero value by a factor set by the input. Heatmap colobar
%is accurate even with the color boost

%Inputs
% gridData: cell array of axon length, branchpoint, and endpoint data 
% binned by grid for all cells, including grid parameters. Output from
% GenerateGridData.m
% gridLengthNorm: cell array of normalized lengths in each grid for each
% cell, where each element is a 3D numeric array. total 1 for each cell
% gridBranchNorm: cell array of normalized branch points in each grid for 
% each cell, where each element is a 3D numeric array. total 1 for each
% gridEndNorm: cell array of normalized end points in each grid for each
% cell, where each element is a 3D numeric array. total 1 for each
% group: array of integers specifying which group/condition each cell is
% ageGroups: labels for the groups, in order of ascending integer value
% colorBoost: number added to all heatmap values that are nonzero, to
% increase contrast. is a unitless proportion because all heatmaps are
% normalized before averaging. 

function [ranges] = display2DHeatmaps(gridData, gridLengthNorm, gridBranchNorm, gridEndNorm, group, ageGroups, colorBoost)
%% Unify single cells into average heatmaps for each age
%Note: very important that grids with no contained axon measurements are 0,
%not NaN
%find parameters for generalized cell-centered heatmaps
ranges = {[0 0], [1 0], [0 0]}; %x,y,z order
for j = 1:length(gridLengthNorm)
    if min(gridData{j}.gridX) < ranges{1}(1)
        ranges{1}(1) = min(gridData{j}.gridX);
    end
    if max(gridData{j}.gridX) > ranges{1}(2)
        ranges{1}(2) = max(gridData{j}.gridX);
    end
    if min(gridData{j}.gridY) < ranges{2}(1)
        ranges{2}(1) = min(gridData{j}.gridY);
    end
    if max(gridData{j}.gridY) > ranges{2}(2)
        ranges{2}(2) = max(gridData{j}.gridY);
    end
    if min(gridData{j}.gridZ) < ranges{3}(1)
        ranges{3}(1) = min(gridData{j}.gridZ);
    end
    if max(gridData{j}.gridZ) > ranges{3}(2)
        ranges{3}(2) = max(gridData{j}.gridZ);
    end
end

%Make index arrays for the three spatial dimensions
alignIndex = cell(size(ranges)); 
xzLim = max(abs([ranges{1} ranges{3}]));
alignIndex{1} = (-xzLim):1:xzLim; %Scale X and Z axis to be square/symmetric
alignIndex{2} = ranges{2}(1):1:ranges{2}(2); %Define y-axis using ranges
alignIndex{3} = (-xzLim):1:xzLim; %Scale X and Z axis to be square/symmetric

%Below: alternative option to not make X and Z square
%{
for i = 1:length(ranges)
    alignIndex{i} = ranges{i}(1):1:ranges{i}(2);
end
%}

%% convert grid data input to a nx3 cell array containing the three point types
grid3 = [gridLengthNorm, gridBranchNorm, gridEndNorm]; %Store 3 variables together

%initialize maps of the right dimensions: note that y-values are dim.1 
baseMap = zeros(length(alignIndex{1}), length(alignIndex{2}),...
    length(alignIndex{3}));
%data structures for storing age-specific group maps of each variable
alignedHeatmaps = cell(length(ageGroups), size(grid3, 2)); 

indivMaps = cell(size(grid3,1), 1); %individual aligned heatmaps for length only

%% populate wm-aligned heatmaps
for i = 1:size(grid3, 2) %loop through measurements
    for j = 1:length(ageGroups)
        alignedHeatmaps{j,i} = baseMap;
        for k = 1:size(grid3, 1) % loop through cells
            if group(k) == j
                tempMap = baseMap;
                tempRange = {[min(gridData{k}.gridX), max(gridData{k}.gridX)],...
                    [min(gridData{k}.gridY), max(gridData{k}.gridY)],...
                    [min(gridData{k}.gridZ), max(gridData{k}.gridZ)]};
                for l = 1:length(tempRange) %match range to elements in 
                    for m = 1:length(tempRange{l})
                        tempRange{l}(m) = find(tempRange{l}(m) == alignIndex{l});
                    end
                end
                %Populate temporary map with grid data from individual cell
                tempMap(tempRange{1}(1):tempRange{1}(2),...
                    tempRange{2}(1):tempRange{2}(2),...
                    tempRange{3}(1):tempRange{3}(2)) = grid3{k, i}; 
                if nnz(isnan(tempMap)) > 0 %check for NaN values
                    disp(k)
                end
                %add cell to group average
                alignedHeatmaps{j,i} = alignedHeatmaps{j,i} + tempMap;
                if i == 1 %Store aligned individual maps for length for plotting projections w/ errors later
                    indivMaps{k} = tempMap; %normalized within cell
                end
            end
        end
        %Average within group by dividing by number of cells
        alignedHeatmaps{j,i} = alignedHeatmaps{j,i}/nnz(group == j);
    end
end

%% plot xy, zy, and xz-projected heatmaps for each age group 
measureNames = {' Length (um)', ' Branchpoints', ' Endpoints'}; %for titles

%wm-aligned: nearly exactly same procedure as cell-centered
xy = cell(size(alignedHeatmaps)); %Currently only plotting XY
yz = cell(size(alignedHeatmaps));
xz = cell(size(alignedHeatmaps));
for i = 1:size(alignedHeatmaps,1)
    for j = 1:size(alignedHeatmaps, 2)
        xy{i,j} = sum(alignedHeatmaps{i,j},3);
        for k = 1:numel(xy{i,j})
            if xy{i,j}(k) > 0 %boost value if nonzero
                xy{i,j}(k) = xy{i,j}(k) + colorBoost;
            end
        end
        yz{i,j} = sum(alignedHeatmaps{i,j},1);
        yz{i,j} = permute(yz{i,j}, [3 2 1]); %Leaves y as second dimension
        for k = 1:numel(yz{i,j})
            if yz{i,j}(k) > 0 %boost value if nonzero
                yz{i,j}(k) = yz{i,j}(k) + colorBoost;
            end
        end
        
        xz{i,j} = sum(alignedHeatmaps{i,j},2);
        xz{i,j} = permute(xz{i,j}, [1 3 2]);
        for k = 1:numel(xz{i,j})
            if xz{i,j}(k) > 0 %boost value if nonzero
                xz{i,j}(k) = xz{i,j}(k) + colorBoost;
            end
        end
    end
end

figure('Name', 'xy projections')
%Find x and y labels
hmXlabels = sign(alignIndex{1}) + fix(alignIndex{1});
hmXlabels(~any(abs(hmXlabels') == [1 5 10 15 20], 2)) = NaN;
hmYlabels = string(alignIndex{2});
hmYlabels(1) = 'WM';
hmYlabels(length(hmYlabels)) = 'Pia';
for i = 2:(length(hmYlabels)-1)
    if rem(alignIndex{2}(i),5) ~= 0
        hmYlabels(i) = '';
    end
end
hmYlabels = flip(hmYlabels);
%Make heatmaps
for i = 1:length(ageGroups)
    for j = 1:3 %Loop through variables
        subplot(length(ageGroups), 3, (i-1)*3 + j)
        %Flips y values so positive is upwards
        h = heatmap(alignIndex{1}, flip(alignIndex{2}), ...
            flip(xy{i,j}', 1), 'FontSize', 6);     
        h.GridVisible = 'off';
        h.XDisplayLabels = hmXlabels;
        h.YDisplayLabels = hmYlabels;
        title([ageGroups{i} measureNames{j}]);
        set(gcf, 'Position', [50 50 900 250*numel(ageGroups)])
        xlabel('tangential axis 1')
        ylabel('pia-wm axis')
        %correct tick labels for boost value
        tempAxis = struct(gca);
        tempCB = tempAxis.Colorbar;
        tempTicks = tempCB.Ticks - colorBoost;
        tempTicks(tempTicks < 0) = 0;
        tempCB.TickLabels = tempTicks;
    end
end
colormap(gray);

figure('Name', 'zy projections')
%Find z labels
hmZlabels = sign(alignIndex{3}) + fix(alignIndex{3});
hmZlabels(~any(abs(hmZlabels') == [1 5 10 15 20], 2)) = NaN;

%Make heatmaps
for i = 1:length(ageGroups)
    for j = 1:3 %Loop through variables
        subplot(length(ageGroups), 3, (i-1)*3 + j)
        %Flips y values so positive is upwards
        h = heatmap(alignIndex{3}, flip(alignIndex{2}), ...
            flip(yz{i,j}', 1), 'FontSize', 6);     
        h.GridVisible = 'off';
        h.XDisplayLabels = hmZlabels;
        h.YDisplayLabels = hmYlabels;
        title([ageGroups{i} measureNames{j}]);
        set(gcf, 'Position', [50 50 900 250*numel(ageGroups)])
        xlabel('tangential axis 2')
        ylabel('pia-wm axis')
        %correct tick labels for boost value
        tempAxis = struct(gca);
        tempCB = tempAxis.Colorbar;
        tempTicks = tempCB.Ticks - colorBoost;
        tempTicks(tempTicks < 0) = 0;
        tempCB.TickLabels = tempTicks;
    end
end
colormap(gray);
%}

figure('Name', 'xz projections')

%Make heatmaps
for i = 1:length(ageGroups)
    for j = 1:3 %Loop through variables
        subplot(length(ageGroups), 3, (i-1)*3 + j)
        %Flips y values so positive is upwards
        h = heatmap(alignIndex{1}, alignIndex{3}, ...
            xz{i,j}', 'FontSize', 6);     
        h.GridVisible = 'off';
        h.XDisplayLabels = hmXlabels;
        h.YDisplayLabels = hmZlabels;
        title([ageGroups{i} measureNames{j}]);
        set(gcf, 'Position', [50 50 900 250*numel(ageGroups)])
        xlabel('tangential axis 1')
        ylabel('tangential axis 2')
        %correct tick labels for boost value
        tempAxis = struct(gca);
        tempCB = tempAxis.Colorbar;
        tempTicks = tempCB.Ticks - colorBoost;
        tempTicks(tempTicks < 0) = 0;
        tempCB.TickLabels = tempTicks;
    end
end
colormap(gray);
%}

end

%% Function to plot length and convex hull by depth
% Plots axon length and horizontal convex hull area as a function of
% cortical depth. 

%Inputs
% allCells: cell array of all the cells in the dataset, output from
% ProcessReconstruction.m.
% gridData: cell array of data binned by grid for each cell, output from
% GenerateGridData.m
% yGridIndex: cell array of 1 vertical vector per cell, specifying the y
% grid index that the corresponding point in the tracing belongs to. Can be
% zeroed on wm boundary or centered on the cell. 
% yRange: range of y-values in the yGridIndex
% group: array of integers specifying which group/condition each cell is
% ageGroups: labels for the groups, in order of ascending integer value

function plotTracingsByDepth(allCells, gridData, yGridIndex, yRange, group, ageGroups)
colors = {'c' 'm' 'g' 'r' 'b'}; %excess to use for up to 5  groups

%% Find values for all cells
yValues = yRange(1):1:yRange(2);
binnedLength = zeros(length(yValues), length(allCells)); %zero if no data
binnedCVH = zeros(length(yValues), length(allCells));%zero if no data

for j = 1:length(allCells)
    for i = 1:length(yValues)
        %Add relevant lengths
        binnedLength(i,j) = sum(allCells{j}.lengths(yGridIndex{j} == yValues(i)));
        %transfer correct convex hull area
        if nnz(gridData{j}.gridY == yValues(i)) > 0
            binnedCVH(i,j) = gridData{j}.cvhArea(gridData{j}.gridY == yValues(i));
        end
    end
end

%% Plot individual cells, colored by group
binnedData = {binnedLength, binnedCVH};
xlabels = {'Axon length (um)', 'Convex hull area (um^2)'};
figure
for i = 1:length(binnedData)
    subplot(1,2,i)
    hold on
    for j = 1:size(binnedData{i}, 2)
        plot(binnedData{i}(:,j), yValues, '-', 'Color', colors{group(j)})
    end
    xlabel(xlabels{i})
    figQuality(gcf, gca, [6 6])
end

%% Plot averages, colored by group
figure
for i = 1:length(binnedData)
    subplot(1,2,i)
    hold on
    for j = 1:length(ageGroups)
        tempData = binnedData{i}(:, group == j);
        errorbar(mean(tempData, 2), yValues,...
            std(tempData, [], 2)/sqrt(size(tempData, 2)), 'horizontal',...
            '-o', 'Color', colors{j}, 'LineWidth', 1.5, 'CapSize', 0, ...
            'MarkerFaceColor', colors{j});
    end
    legend(ageGroups)
    xlabel(xlabels{i})
    figQuality(gcf, gca, [6 6])
end

end