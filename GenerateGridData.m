%Axon analysis and quantification code associated with Gutman-Wei et al.,
%2025 - function to calculate tracing statistics with a grid overlaid
%Generates quantification of length, branchpoints, and endpoints in each
%grid for cell input, as well as convex hulls at each level specified. The
%location of grid vertices and the assignment of each point within the
%tracing is completed after running ProcessReconstruction. This function
%reads the grid indices of the points for a tracing, and uses the length,
%point ID (branch/end), and grid index relative to white matter-pia axis to
%quantify the amount of axon length and number of axon branch points and
%end points contained within each grid unit. User can specify whether data
%from cell-centered or wm-aligned grids are reported. Processes one cell at
%a time. 

%Inputs:
% pointGridIndex: 3-column array, with each row corresponding to a point in
% the tracing. Values indicate the x(1), y(2), and z(3) grids the point
% belongs to as the x-, y-, or z-value of the center of that grid. y is
% cell-centered by default, but can be replaced with wm grid values in the
% function if specified by user
% lengths: distance between individual points of the tracing for all points
% pointType: 2-column logical indicating whether a point is a branch point
% (column 1) or an end point (column 2)
% points: 3-column array with x, y, and z values of all points
% wmGridIndex: grid location in the y dimension of all points in the
% tracing, with a grid indexed beginning from 0.5 (close to wm). Should be
% the same height as all other tracing data (1 row for each point), and
% calculated by user from pointHeights (see RunTracingAnalysis.m)
% type: string input of 'wm' or 'cell'. Determines whether the y-grid data
% is calculated from a grid centered on the cell, or from a grid scaled
% from wm to pia. 

%Outputs (no plots):
% gridData: structure with fields that contain 3-D arrays that quantify the
% total length, branch points, and end points in each grid. Also contains
% the 2D convex hull in the X-Z dimension of points in each horizontal
% slice (i.e. all points with the same y-grid). Reports the grid X, Y, and
% Z values contained within the tracing, and reports points excluded from
% the grid data (only those not within the wm-pia grid if using 'wm')

function [gridData] = GenerateGridData(pointGridIndex, lengths, pointType,...
    points, wmGridIndex, type)
%% Set y-grid according to type
if contains(type, 'wm')
    pointGridIndex(:,2) = wmGridIndex;
    pointsExcluded = isnan(wmGridIndex);
    if nnz(pointsExcluded) == 0
        pointsExcluded = 0;
    end
elseif contains(type, 'cell')
    pointsExcluded = 0;
else
    disp('incorrect selection of grid type')
    return
end

%% find grids
gridX = unique(pointGridIndex(:,1));
gridY = unique(pointGridIndex(:,2));
gridZ = unique(pointGridIndex(:,3));
%Remove nans, which are points outside of 0-1 range in normalized height
%(outside of the y-range spanned by the rotation vector points)
gridY(isnan(gridY)) = [];

gridLength = nan(length(gridX), length(gridY), length(gridZ));
gridBranch = nan(length(gridX), length(gridY), length(gridZ));
gridEnd = nan(length(gridX), length(gridY), length(gridZ));
for i = 1:length(gridX)
    for j = 1:length(gridY)
        for k = 1:length(gridZ)
            tempIndex = (pointGridIndex(:,1) == gridX(i)) &...
                (pointGridIndex(:,2) == gridY(j)) & (pointGridIndex(:,3) == gridZ(k));
            gridLength(i,j,k) = sum(lengths(tempIndex));
            gridBranch(i,j,k) = sum(pointType(tempIndex, 1));
            gridEnd(i,j,k) = sum(pointType(tempIndex, 2));
        end
    end
end

%% find hulls
%Following 2 sections create 2d convex hulls of all points in the cell at
%each level of y used in the grid (convex hulls are in the x-z plane), and
%plots them in 3d. Uses unscaled grid with um grid size used throughout
%the script, and shows horizontal extend of the cell binned by cortical
%depth. Also creates 2d angle plots of axon length extending in each
%direction of X-Z plane at each y level
yHulls = cell(size(gridY)); %cell array containing vertices of 2-D convex hulls at each level of y
gridCVH = zeros(size(gridY)); %convex hull area for each level of y
for i = 1:length(yHulls) %Looping through y values
    tempX = points(pointGridIndex(:,2) == gridY(i), 1);
    tempZ = points(pointGridIndex(:,2) == gridY(i), 2);
    if nnz(tempX) > 2
        [hullPoints, gridCVH(i)] = convhull(tempX, tempZ); %2d convex hull using only x and z values of points contained
        yHulls{i} = [tempX(hullPoints), tempZ(hullPoints)]; %stores x and z values as 2-column matrix
    end
end

%% Make structure
gridData = struct('gridLength', gridLength, 'gridBranch', gridBranch,...
    'gridEnd', gridEnd, 'cvhArea', gridCVH, 'gridX', gridX, 'gridY',...
    gridY, 'gridZ', gridZ, 'pointsExcluded', pointsExcluded); 

end



