% classify_edge_points.m
addpath('C:\Users\aashr\ISSM-Windows-MATLAB\bin');
load('data_pinns_Amery.mat');

% Toggle plotting on/off
doPlot = false;  % set to false to skip plotting

% Inputs
thickness_nc = 'ice_thickness.nc';
T = readtable('DomainOutlinePoints.csv', 'Delimiter', ' ', 'ReadVariableNames', true);

% Grid and mask
tx = double(ncread(thickness_nc, 'x'));
ty = double(ncread(thickness_nc, 'y'));
mesh = double(ncread(thickness_nc, 'mask'));

% Points (transform from table)
xgl = -T.y + 2.15e6;
ygl = -T.x + 2.55e6;

% Ensure xgrid, ygrid are column vectors
xgrid = tx(:);
ygrid = ty(:);

% Subregion limits (same as before)
xmin = -1146390;
xmax = -175452;
ymin = -2279238;
ymax = -1663947;
xlo = min(xmin,xmax); xhi = max(xmin,xmax);
ylo = min(ymin,ymax); yhi = max(ymin,ymax);

ix = find(xgrid >= xlo & xgrid <= xhi);
iy = find(ygrid >= ylo & ygrid <= yhi);

Z = mesh(iy, ix);
x_sub = xgrid(ix);
y_sub = ygrid(iy);

% Get subregion mask value at nearest node for each point
nPoints = numel(xgl);
maskVals_sub = NaN(nPoints,1);
ix_sub_all = ix;
iy_sub_all = iy;
mapX = NaN(numel(xgrid),1); mapX(ix_sub_all) = 1:numel(ix_sub_all);
mapY = NaN(numel(ygrid),1); mapY(iy_sub_all) = 1:numel(iy_sub_all);

% Precompute nearest grid indices for each point (global)
ix_nearest = zeros(nPoints,1);
iy_nearest = zeros(nPoints,1);
for k = 1:nPoints
    [~, ix_k] = min(abs(xgrid - xgl(k)));
    [~, iy_k] = min(abs(ygrid - ygl(k)));
    ix_nearest(k) = ix_k;
    iy_nearest(k) = iy_k;
    lx = mapX(ix_k); ly = mapY(iy_k);
    if ~isnan(lx) && ~isnan(ly)
        maskVals_sub(k) = Z(ly, lx);
    end
end

% Check neighborhood neighborhood (using grid indices) for any ocean cells (mask == 0)
windowSize = 12;
halfEven = windowSize/2;                % 3
offsets = -(halfEven-1) : halfEven;     % [-2 -1 0 1 2 3]

isOceanNearbyneighborhood = false(nPoints,1);
numOceanCells = zeros(nPoints,1);

nx = numel(xgrid);
ny = numel(ygrid);

for k = 1:nPoints
    ix_k = ix_nearest(k);
    iy_k = iy_nearest(k);

    ix_range = ix_k + offsets;
    iy_range = iy_k + offsets;
    ix_range = ix_range(ix_range >= 1 & ix_range <= nx);
    iy_range = iy_range(iy_range >= 1 & iy_range <= ny);

    if isempty(ix_range) || isempty(iy_range)
        continue
    end

    neighborhood = mesh(iy_range, ix_range);   % rows->y, cols->x
    nOcean = sum(neighborhood(:) == 0);
    numOceanCells(k) = nOcean;
    isOceanNearbyneighborhood(k) = (nOcean > 0);
end

% Plot subregion and mark points with ocean nearby as crosses
if doPlot
    figure
    imagesc(x_sub, y_sub, Z);
    axis equal
    set(gca, 'YDir', 'normal')
    hold on

    % points with ocean in neighborhood -> black 'x'
    idxOcean = find(isOceanNearbyneighborhood);
    if ~isempty(idxOcean)
        plot(xgl(idxOcean), ygl(idxOcean), 'kx', 'MarkerSize', 8, 'LineWidth', 1.5);
    end

    % other valid subregion points -> red filled circles
    idxValid = find(~isnan(maskVals_sub) & ~isOceanNearbyneighborhood);
    if ~isempty(idxValid)
        plot(xgl(idxValid), ygl(idxValid), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    end

    xlabel('X (m)')
    ylabel('Y (m)')
    title('Mask Subregion — points with ocean in neighborhood neighborhood marked with X')
    legend({'Ocean within neighborhood','Other points'}, 'Location', 'bestoutside')
    hold off
end

% assume xgl, ygl, and isOceanNearbyneighborhood are in workspace
inRange = logical(isOceanNearbyneighborhood(:));   % true = within range (ocean nearby)
n = numel(xgl);

% Prepare columns
transformed_x = xgl(:);
transformed_y = ygl(:);

% grounding_line = 0 when inRange, 1 otherwise
grounding_line = double(~inRange);   % ~inRange -> 1 means not in range => grounding_line=1
grounding_line(inRange) = 0;         % ensure explicitness

% ice_edge = 1 when inRange, 0 otherwise
ice_edge = double(inRange);
ice_edge(~inRange) = 0;

% Build table and write to CSV
Tout = table(transformed_x, transformed_y, grounding_line, ice_edge, ...
    'VariableNames', {'transformed_x','transformed_y','grounding_line','ice_edge'});

writetable(Tout, 'points_export.csv');
% disp(Tout)
