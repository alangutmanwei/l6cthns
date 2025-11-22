%Axon analysis and quantification code associated with Gutman-Wei et al.,
%2025 - function for initial processing of reconstruction
% Returns reconstruction with WM axon and small branches removed, and a 
% series of indices matching each xyz value of reconstruction with point 
% type, height, layer identity, and other information about the tracing.

% Requires the SWC file for the cell to be accessible from
%the local directory (relative path used). User toggles graphics
%windows depending on the processing speed desired.

%Single SWC file is processed at one time. File is read as a matrix, then
%branch points and end points are found by querying swc file for number of
%offspring points. >1 offspring = branch point, 0 offspring = endpoint.
%Lengths are calculated as euclidean distance between consecutive points in
%the tracing, and indexed at the second point in each such distance vector.
%Total length for each branch is calculated, and branches with lengths <
%user threshold (usually 10um) are discarded. The points associated with
%these branches are deleted from the tracing matrix, and then branchpoints
%and endpoint are re-identified after this culling. Cell is then rotated
%such that the vector from WM to Pia, specified by user input, is parallel
%to the y-axis. This is accomplished via 2 sets of rotations about the
%z-axis, then the x-axis, and leaving the cell at an arbitrary rotation
%about the y-axis. After rotation, any processes below the white matter are
%removed. A scaled grid is created according to the user's
%specified grid size, and all points in the tracing are assigned to a grid
%element. The grid is centered on the cell soma overall. Finally, the cell
%is plotted as a figure with all axon points shown, and branch points end
%points, cell soma, primary axon, wm location, and pia location plotted in
%different colors. 

%Inputs: 
%filename: relative path to SWC file (can have .swc or  .txt suffix)
%piaXYZ: 3-element vector of [x,y,z] of pia intersection point with axis of
%cell (obtained manually and pulled from metadata spreadsheet)
%wmXYZ: same as piaXYZ for white matter intersection point of cell axis
%minlength: minimum length of branches before they are removed from
%analysis, typically 10um
%gridProp: Size of grid elements, as a fraction of the pia-WM distance,
%typically 0.05 (5% of cortical thickness)
%figuresOn: true = plot cell after rotation and removal of small branches
%and WM axon, with branchpoints, endpoints, soma, and primary axon colored.
%false to report data much more quickly

%Outputs: single structure tracingStats containing:
% totLength: total length of tracing after removing WM axon, and removing
% short branches
% totBranch: total branch point number of tracing after removing WM axon 
% and removing short branches
% totEnd: total end point number of tracing after removing WM axon 
% and removing short branches, typically totBranch + 1
% cortThickness: thickness of cortex in this tracing, based on wm-pia
% vector specified
% points: 3-column matrix with x,y,z values of the full reconstruction after
% all processing
% lengths: array with the length value assigned to each xyz point in the
% reconstruction. Length is zero for primary initiation point and value for
% each voxel-length segment is assigned to the endpoint of the segment
% pointType: 2-column matrix of logicals, first column true/false for
% branchpoint, second for endpoint
% pointHeight: height of each point as a proportion of cortical thickness,
% with zero being white matter boundary
% pointGridIndex:index value of the grid that each point is in, given as
% the index value of the center of the grid in 3 dimensions. this grid is
% automatically centered such that the soma/axon initiation site is 0
% layerID: column vector indicating the layer each point is assigned to
% based on layer annotations included in the modified SWC file. If SWC file
% does not contain an 8th column of layer annotation, returns NaN
% trimmedSWC: all data in SWC format after points in WM or in branches
% <10um are trimmed
% wmPia: position of the wm and pia after rotation, relative to cell soma
% as origin (top row wm, bottom row pia)
% soma: original soma location before rotation 

%% Inputs and setup
function [tracingStats] = ProcessReconstruction(fileName, piaXYZ, wmXYZ, minLength, gridProp, figuresOn)
%Read file in
data = readmatrix(fileName, 'FileType', 'text');
%Find x, y, and z: columns assumed by typical output from SNT
x = data(:,3);
y = data(:,4);
z = data(:,5);


%% Reject if more than 1 primary path
if nnz(data(:,7) == -1) > 1
    disp('more than 1 primary path')
    tracingStats = NaN;
    return
end

%% Zero on soma
somaPosition = [x(1), y(1), z(1)];
x = x - somaPosition(1);
y = y - somaPosition(2);
z = z - somaPosition(3);

%% Find segment lengths
%all distances between points are scaled in microns, but points are not
%necessarily all the same distance apart - calculates distances in xyz
lengths = zeros(size(data, 1), 1);
for i = 2:size(data, 1) %Skip entry 1 which has gap size 0
    parentLocation = find(data(:,1) == data(i,7));
    lengths(i) = pdist2(data(i,3:5), data(parentLocation,3:5));
end

%% Find initial set of branchpoints and endpoints, pre minimum length filter
branchPoints = false(size(lengths)); %True if entry is a branchpoint
endPoints = false(size(lengths)); %True if entry is an endpoint
for i = 1:size(data,1)
    pointIndex = data(i,1); %temporary
    numOffspring = nnz(data(:,7) == pointIndex); %Number of offspring points from this point
    branchPoints(i) = numOffspring > 1; %True for all points with more than one offpsring
    endPoints(i) = numOffspring == 0; %True for all points that have no offspring
end

%% Segment out branches, pre-minimum length filter
branches = nan(nnz(branchPoints) + nnz(endPoints), 2); %start and end indices(columns) of each branch(rows)
branchPointIndices = find(branchPoints); %indices of all branch points
%First branch is primary axon to first branchpoint
branches(1,:) = [1 branchPointIndices(1)]; 
branchIndex = 2; %Initialize: start by filling in second branch, 1st is primary
%Rest of branches
for i = 1:length(branchPointIndices) %Loop through all points where branches start
    startPoint = branchPointIndices(i);
    offspring = find(data(:,7) == startPoint); %Find offspring (usually 2)
    branchEnds = nan(size(offspring)); %initialize storage for ends of branches
    for j = 1:length(offspring) %loop through all branches (usually 2)
        k = offspring(j); %counter that moves along branch 
        while ~branchPoints(k) && ~endPoints(k) %stops at next branch or end point
            k = k + 1;
        end
        branchEnds(j) = k; %Stores branch end points
        branches(branchIndex + (j-1), :) = [offspring(j) branchEnds(j)]; %Save branch start and end points
    end
    branchIndex = branchIndex + nnz(offspring); %increment the branch index
end

%% Remove branches < threshold
%Set path length interval for calculating tortuosity
branchLengths = nan(size(branches, 1), 1); %total lengths of each branch
for i = 1:size(branches, 1) %Loop through all branches
    branchLengths(i) = sum(lengths(branches(i,1):branches(i,2)));
end

branchesToRemove = branchLengths < minLength;
%Make sure only short branches that are endings are removed
for i = 1:length(branchesToRemove)
   if branchesToRemove(i)
       if ~endPoints(branches(i,2))
           branchesToRemove(i) = false;
       end
   end
end

%% Redefine x, y, and z values after removing branches 
pointsToRemove = false(size(x));
for i = 1:length(branchesToRemove)
   if branchesToRemove(i)
      pointsToRemove(branches(i,1):branches(i,2)) = true; 
   end
end

%Remove points associated with too-short branches
x(pointsToRemove) = [];
y(pointsToRemove) = [];
z(pointsToRemove) = [];
lengths(pointsToRemove) = []; %lengths are not redefined - old array removed
data(pointsToRemove, :) = [];

disp(['number of short branches removed: ' num2str(nnz(branchesToRemove))]);
%Note: does not cull branches or remake branch cell array
%% Find branchpoints and endpoints after filtering small branches
branchPoints = false(size(x)); %True if entry is a branchpoint
endPoints = false(size(x)); %True if entry is a an endpoint
for i = 1:size(data,1)
    pointIndex = data(i,1); %temporary
    numOffspring = nnz(data(:,7) == pointIndex); %Number of offspring points from this point
    branchPoints(i) = numOffspring > 1; %True for all points with more than one offpsring
    endPoints(i) = numOffspring == 0; %True for all points that have no offspring
end

%% Rotate cell based on soma-pia vector
%Find rotation vector
%piaVector = piaXYZ - somaPosition; %rotation vector from soma to pia
piaVector = piaXYZ - wmXYZ; %rotation vector from wm to pia
piaXYZ = piaXYZ - somaPosition;
wmXYZ = wmXYZ - somaPosition;

%Find angle for rotation about x (relative to y axis)
xRot = -atan(piaVector(3)/piaVector(2));
%Make rotation matrix
xRM = [1 0 0; 0 cos(xRot) -sin(xRot); 0 sin(xRot) cos(xRot)];

%Rotate all points 
[x, y, z, piaRotate, wmRotate] = rotatePoints(x, y, z, xRM, piaXYZ, wmXYZ);

%Repeat for z axis rotation
%Find angle for rotation about z, but after x rotation is complete
piaVector = piaRotate - wmRotate; %Redefine z rotation needed AFTER x rotation
zRot = atan(piaVector(1)/piaVector(2));
zRM = [cos(zRot) -sin(zRot) 0; sin(zRot) cos(zRot) 0; 0 0 1];
[x, y, z, piaRotate, wmRotate] = rotatePoints(x, y, z, zRM, piaRotate, wmRotate);

%% Find white-matter region of axon and mark for removal
%Finds the point where the primary axon enters white matter, and marks all
%offspring as wm axon. Used for removing from pointsContained array
%elements later in function
wmAxon = false(size(y));
primaryAxon = false(size(y));
primaryAxon(1:find(endPoints, 1)) = true; %Assumes primary axon is first PATH traced
wmAxon(primaryAxon & (y < wmRotate(2))) = true; %find points of primary axon below wm
%loop through points and extend WM status to all downstream nodes - assumes
%nodes are all in appropriate order
for i = 2:length(wmAxon) %Assumes element 1, soma, is not in WM
    if ~wmAxon(i)
        if wmAxon(data(:,1) == data(i,7))
            wmAxon(i) = true;
        end
    end
end

%print length trimmed
disp(['wm axon length trimmed: ' num2str(sum(lengths(wmAxon)))]);
disp(['wm branch points trimmed: ' num2str(sum(branchPoints(wmAxon)))]);


%% Remove WM axon
%Remove points associated with too-short branches
x(wmAxon) = [];
y(wmAxon) = [];
z(wmAxon) = [];
lengths(wmAxon) = []; %lengths are not redefined - old array removed
branchPoints(wmAxon) = [];
endPoints(wmAxon) = [];
data(wmAxon, :) = [];

%Set the truncated primary axon point to be endpoint
primaryAxon(wmAxon) = [];
endPoints(find(primaryAxon, 1, 'last')) = true;
%Note: does not cull branches or remake branch cell array 

%% Set grid size
cortexThickness = pdist2(piaRotate, wmRotate);
gridSize = cortexThickness * gridProp;

%% Generate grid and shape objects - CELL-CENTERED
% MATLAB alphashapes are shape objects that can be compared to points in 3D
% space to determine which points are contained within. Alphashapes are
% created for a 3D grid encompassing the cell and all annotated processes
% in file
maxima = [max(x), max(y), max(z)];
minima = [min(x), min(y), min(z)];
maxima = ceil(maxima/gridSize) * gridSize; %Force endpoints to be multiples of grid size
minima = floor(minima/gridSize) * gridSize; %assumes soma is not at extreme of cell
[gridX, gridY, gridZ] = meshgrid(minima(1):gridSize:maxima(1), ...
    minima(2):gridSize:maxima(2), minima(3):gridSize:maxima(3)); %Grid vertices

gridShapes = cell(size(gridX) - 1); %subtract one for fencepost, this is the number of grids
gridIndex = {nan(size(gridX) - 1), nan(size(gridX) - 1), nan(size(gridX) - 1)}; %index of grid distance from origin in each dimension, in units of grids, x,y,z
vertices = cell(size(gridX) - 1); %store cell array of shape vertices for plotting later

% Create alphashapes for all of the points
for k = 1:size(gridShapes,3)
    for j = 1:size(gridShapes, 2)
       for i = 1:size(gridShapes, 1)
           vertices{i,j,k} = [gridX(i,j,k), gridY(i,j,k), gridZ(i,j,k); ...
               gridX(i + 1, j, k), gridY(i + 1, j, k), gridZ(i + 1, j, k);...
               gridX(i, j + 1, k), gridY(i, j + 1, k), gridZ(i, j + 1, k);...
               gridX(i, j, k + 1), gridY(i, j, k + 1), gridZ(i, j, k + 1);...
               gridX(i + 1, j, k + 1), gridY(i + 1, j, k + 1), gridZ(i + 1, j, k + 1);...
               gridX(i + 1, j + 1, k), gridY(i + 1, j + 1, k), gridZ(i + 1, j + 1, k);...
               gridX(i, j + 1, k + 1), gridY(i, j + 1, k + 1), gridZ(i, j + 1, k + 1);...
               gridX(i + 1, j + 1, k + 1), gridY(i + 1, j + 1, k + 1), gridZ(i + 1, j + 1, k + 1)];
           tempCenter = round(mean(vertices{i,j,k})./gridSize, 1); 
           for l = 1:3
               gridIndex{l}(i,j,k) = tempCenter(l); 
           end
           gridShapes{i,j,k} = alphaShape(vertices{i,j,k}(:,1), ...
               vertices{i,j,k}(:,2), vertices{i,j,k}(:,3), Inf); %Inf generates convex hull
       end
    end
end

%% Count number of points and length within each shape
pointsContained = cell(size(gridShapes)); %Indices of points within the shape
pointIndex = nan(length(x), 3);
for k = 1:size(gridShapes, 3)
   for j = 1:size(gridShapes, 2)
      for i = 1:size(gridShapes, 1)
          pointsContained{i,j,k} = inShape(gridShapes{i,j,k}, x, y, z);
          for l = 1:3
            pointIndex(pointsContained{i,j,k}, l) = gridIndex{l}(i,j,k);
          end
      end
   end
end


%% Layer info
if size(data, 2) == 8
    layerID = data(:,8);
else
    layerID = NaN;
end

%% Plot cell
wmPia = [wmRotate; piaRotate];
if figuresOn
    figure
    hold on
    plotX = x(~(branchPoints));
    plotY = y(~(branchPoints));
    plotZ = z(~(branchPoints));
        %Plot cell with color coded points
    scatter3(0, 0, 0, 50, 'r', 'filled')
    scatter3(x(branchPoints), y(branchPoints), ...
        z(branchPoints), 30, 'k', 'filled')
    %Plot wm/pia boundaries
    scatter3(wmPia(:,1), wmPia(:,2), wmPia(:,3), 50, 'y', 'filled');
    if ~isnan(layerID)
        plotLayer = layerID(~branchPoints);
        restLayers = (plotLayer < 4 | plotLayer > 6);
        scatter3(plotX(plotLayer == 6), plotY(plotLayer == 6),...
            plotZ(plotLayer == 6), 10, 'g', 'filled')
        scatter3(plotX(plotLayer == 5), plotY(plotLayer == 5),...
            plotZ(plotLayer == 5), 10, 'm', 'filled')
        scatter3(plotX(plotLayer == 4), plotY(plotLayer == 4),...
            plotZ(plotLayer == 4), 10, 'c', 'filled')
        scatter3(plotX(restLayers), plotY(restLayers),...
            plotZ(restLayers), 10, [0.5 0.5 0.5], 'filled')
        legend({'soma', 'branch points', 'wm/pia', 'L6', 'L5', 'L4', 'other axon'})
    else
        scatter3(plotX, plotY, plotZ, 10, [0.5 0.5 0.5], 'filled')
        scatter3(x(primaryAxon), y(primaryAxon), z(primaryAxon), 10, 'b', 'filled')
        legend({'soma', 'branch points', 'wm/pia', 'axon', 'primary axon'})
    end
    xlabel('tangential axis 1')
    ylabel('wm-pia axis')
    zlabel('tangential axis 2')
    xticks(unique(gridX))
    yticks(unique(gridY))
    zticks(unique(gridZ))
    ylimits = ylim;
    ylim([ylimits(1), piaRotate(2)]);
    grid on
    set(gcf, 'Position', [100 100 600 600])
end

%% Make structure
tracingStats = struct('totLength', sum(lengths), 'totBranch', ...
    sum(branchPoints), 'totEnd', sum(endPoints), 'cortThickness',...
    cortexThickness, 'points', [x y z], 'lengths', lengths,...
    'pointType',[branchPoints endPoints], 'pointHeight', (y - wmRotate(2))/cortexThickness,...
    'pointGridIndex', pointIndex, 'layerID', layerID, 'trimmedSWC', data, 'wmPia', wmPia,...
    'soma', somaPosition);
end

%% Function for rotating points 
function [x, y, z, piaRotate, wmRotate] = rotatePoints(x, y, z, rotMatrix, piaXYZ, wmXYZ)
tempPoints = nan(length(x), 3);
for i = 1:length(x)
    tempPoints(i,:) = (rotMatrix * [x(i); y(i); z(i)])'; 
end
x = tempPoints(:,1);
y = tempPoints(:,2);
z = tempPoints(:,3);
%Rotate pia and wm points
piaRotate = piaXYZ;
piaRotate = (rotMatrix * piaRotate')';
wmRotate = wmXYZ;
wmRotate = (rotMatrix * wmRotate')';
end