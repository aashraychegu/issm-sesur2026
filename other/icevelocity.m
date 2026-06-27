clear;
clc;
addpath('C:\Users\aashr\ISSM-Windows-MATLAB\bin');
load('data/data_pinns_Amery.mat')
% Create mesh / triangulation if needed (this creates a new md.mesh)
md = bamg(model, 'domain', 'data/DomainOutline.exp', 'hmin', 2500);
md.mesh;

ncfile = 'antarctica_ice_velocity_1995-2001_450m_v01.1.nc';

% --- Inspect NetCDF contents and try to pick coordinate/data vars ---
info = ncinfo(ncfile);
varnames = {info.Variables.Name};

crs = ncread(ncfile, 'crs');
x = ncread(ncfile, 'x');
y = ncread(ncfile, 'y');
ud = ncread(ncfile, 'VX');
vd = ncread(ncfile, 'VY');

[X, Y] = meshgrid(x, y);

% Flatten arrays
xv = X(:);
yv = Y(:);
udv = ud(:);
vdv = vd(:);

% Filter out points where ud or vd are NaN (keep only valid velocity points)
validMask = ~isnan(udv) & ~isnan(vdv);

% Optionally, if you meant to remove points that ARE NaN, use the same mask.
% If you wanted the opposite, invert mask: validMask = isnan(udv) | isnan(vdv);

xv = xv(validMask);
yv = yv(validMask);
udv = udv(validMask);
vdv = vdv(validMask);

% Now lengths must match
n = numel(xv);
ptsTable = table(xv, yv, udv, vdv, 'VariableNames', {'x', 'y', 'ud', 'vd'});

% --- Build triangulation from ISSM mesh ---
% md.mesh.elements is expected as Nx3 (node indices) and md.mesh.x, y as column vectors.
% Ensure column orientation:
nodes = [-1 * md.mesh.y(:) + 2.15e6,-1 * md.mesh.x(:)+2.55e6];
elems = md.mesh.elements; % typically 3 x NT or NT x 3 - ensure it's NT x 3
if size(elems,1) == 3 && size(elems,2) ~= 3
    elems = elems'; % make NT x 3
end
% Create triangulation object
TR = triangulation(elems, nodes);

% Use pointLocation to get triangle index for each point (NaN if outside)
triIdx = TR.pointLocation([ptsTable.x, ptsTable.y]);

% Logical mask of points inside mesh
insideMask = ~isnan(triIdx);

% Subset table to points inside mesh
pointsInsideTable = ptsTable(insideMask, :);

% Optional: show counts
fprintf('Total points: %d\n', height(ptsTable));
fprintf('Points inside mesh: %d\n', height(pointsInsideTable));

% Plot mesh and every 100th sampled grid point
figure('Name','Mesh and Sampled Velocity Points','Color',[1 1 1]);
hold on;
axis equal;
box on;

% Plot mesh edges (use triplot on node coordinates)

% Choose which points to sample. Use pointsInsideTable (valid points inside mesh).
pts = ptsTable;


% Define clipping bounding box
xlim_clip = [-971171, -398629];
ylim_clip = [-2331535, -1667448];

% Logical mask for points inside the clipping box
regionMask = pts.x >= xlim_clip(1) & pts.x <= xlim_clip(2) & ...
             pts.y >= ylim_clip(1) & pts.y <= ylim_clip(2);

% Keep only points in region
pts = pts(regionMask, :);

% Sample every 100th point from those in-region points
step = 1;
if isempty(pts)
    warning('No points found inside the specified region.');
    idx = [];
else
    idxAll = 1:height(pts);
    idx = idxAll(1:step:end);
end

% compute speed for color (use pts after filtering)
speed = hypot(pts.ud, pts.vd);

% Scatter sampled points colored by speed (no marker border)
scatter(pts.x(idx), pts.y(idx), 2, speed(idx), 'filled', 'MarkerEdgeColor', 'none');

% Set plot limits to exactly the requested bbox
xlim(xlim_clip);
ylim(ylim_clip);
ylabel('Y');
title('Mesh and Every 100th Velocity Point');


triplot(TR, 'Color', [1 0 0], 'LineWidth', 1); % light gray mesh
colormap(parula);
colorbar()

% Compute bounding box of mesh nodes
xmin = min(nodes(:,1));
xmax = max(nodes(:,1));
ymin = min(nodes(:,2));
ymax = max(nodes(:,2));

% Add a small padding (5% of range)
px = 0.05 * (xmax - xmin);
py = 0.05 * (ymax - ymin);

xlim([xmin-px, xmax+px]);
ylim([ymin-py, ymax+py]);

axis equal tight;
drawnow;

xlabel('X');
title('Mesh and Every 100th Velocity Point');
hold off;
