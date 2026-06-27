clc; clear;
addpath('C:\Users\aashr\ISSM-Windows-MATLAB\bin');

doPlot = false;
load('data/data_pinns_Amery.mat');
md = model;

c1 = 2.15e6;
c2 = 2.55e6;

thickness_nc = 'data/ice_thickness.nc';
velocity_nc  = 'data/antarctica_ice_velocity_1995-2001_450m_v01.1.nc';
md = bamg(md, 'domain', 'data/DomainOutline.exp', 'hmin', 2500);

tx = double(ncread(thickness_nc, 'x'));
ty = double(ncread(thickness_nc, 'y'));
hd = double(ncread(thickness_nc, 'thickness'));
mask = double(ncread(thickness_nc, 'mask'));

vx = double(ncread(velocity_nc, 'x'));
vy = double(ncread(velocity_nc, 'y'));
ud = double(ncread(velocity_nc, 'VX'));
vd = double(ncread(velocity_nc, 'VY'));

%% Transform mesh to the right locations.
md = transformMesh(md,c1,c2);

%% Remesh near the edges
ice_shelf_idx = computeIceShelf(md, tx, ty, mask, 6);

nv = md.mesh.numberofvertices; % Number of vertices

% Initialize ice_levelset
md.mask.ice_levelset = -1 * ones(nv, 1);

% Set ice_levelset for ice shelf points
md.mask.ice_levelset(ice_shelf_idx) = 1;

% Reinitialize level set to create a signed distance field
md.mask.ice_levelset = reinitializelevelset(md, md.mask.ice_levelset);

% Set default resolution
hVertices = 2500 * ones(md.mesh.numberofvertices, 1);

% Refine elements near the ice front
hVertices(abs(md.mask.ice_levelset) < 5e3) = 500;

% Remesh with BAMG
md = bamg(md, 'hVertices', hVertices, 'hmin', 400, 'hmax', 2500);

ice_shelf_idx = computeIceShelf(md, tx, ty, mask, 6);

%% set ice shelf edge conditions
% Initialize ice_levelset
md.mask.ice_levelset = -1 * ones(nv, 1);
% Set ice_levelset for ice shelf points
md.mask.ice_levelset(ice_shelf_idx) = 1;

%% interpolate
md = transformMesh(md,c1,c2);
[hd, ud, vd, original_nans] = interpolateToNodes(md, tx, ty, hd, vx, vy, ud,vd);
md = setupModel(md, hd, ud,vd,ice_shelf_idx,3);

md.inversion.iscontrol = 0;
md = solve(md, 'Stressbalance');

plotmodel(md, 'axis#all', 'tight', ...
    'data', md.results.StressbalanceSolution.MaterialsRheologyBbar, ...
    'caxis', [1.3 1.9] * 1e8, 'title', 'Inferred B', ...
    'data', md.results.StressbalanceSolution.Vel, 'title', 'Modeled velocities');