clear;
clc;
addpath('C:\Users\aashr\ISSM-Windows-MATLAB\bin');
load('data/data_pinns_Amery.mat')

% Create mesh / triangulation if needed (this creates a new md.mesh)
md = bamg(model, 'domain', 'data/DomainOutline.exp', 'hmin', 2500);
md.mesh;

% NetCDF file with thickness
ncfile = 'ice_thickness.nc';

% Read coordinates and thickness
x = ncread(ncfile, 'x');
y = ncread(ncfile, 'y');
h = ncread(ncfile, 'thickness');

[X, Y] = meshgrid(x, y);

% Define the clipping bounding box (icevelocity window)
xlim_clip = [-971171, -398629];
ylim_clip = [-2331535, -1667448];

% Logical mask for points inside the clipping box (flattened)
xv = X(:);
yv = Y(:);
boxMask = xv >= xlim_clip(1) & xv <= xlim_clip(2) & ...
          yv >= ylim_clip(1) & yv <= ylim_clip(2);

% Reshape mask to grid size and trim thickness to box
boxMaskGrid = reshape(boxMask, size(X));
h_clipped = nan(size(h));          % initialize with NaNs outside box
h_clipped(boxMaskGrid) = h(boxMaskGrid);

% --- Build triangulation from ISSM mesh (keep transform as before) ---
nodes = [-1 * md.mesh.y(:) + 2.15e6, -1 * md.mesh.x(:) + 2.55e6];
elems = md.mesh.elements;
if size(elems,1) == 3 && size(elems,2) ~= 3
    elems = elems'; % ensure NT x 3
end
TR = triangulation(elems, nodes);

% Plot clipped thickness and overlay mesh
figure('Name','Ice Thickness (Clipped) with Mesh','Color',[1 1 1]);
p = pcolor(X, Y, h_clipped);
shading interp;
colormap(parula);
colorbar;
hold on;

% Plot mesh edges (triangulation) - keep full mesh overlay (axis limits will clip display)
triplot(TR, 'Color', [0.8 0.2 0.2], 'LineWidth', 0.8);

% Set plot limits to the clipping bounding box
xlim(xlim_clip);
ylim(ylim_clip);

axis equal tight;
xlabel('X');
ylabel('Y');
title('Ice Thickness (m) Clipped To Box with Mesh Overlay');

hold off;
drawnow;
