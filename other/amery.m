clc;clear;
classify_edge_points;
addpath('C:\Users\aashr\ISSM-Windows-MATLAB\bin');

% Toggle plotting on/off
doPlot = false;

% Load ISSM model
load('data_pinns_Amery.mat');
md = model;

% -------------------------
% Read NetCDF files (only necessary data)
% -------------------------
thickness_nc = 'ice_thickness.nc';
velocity_nc  = 'antarctica_ice_velocity_1995-2001_450m_v01.1.nc';
md = bamg(md, 'domain', 'DomainOutline.exp', 'hmin', 2500);

tx = double(ncread(thickness_nc, 'x')); % Explicit double
ty = double(ncread(thickness_nc, 'y')); % Explicit double
hd = double(ncread(thickness_nc, 'thickness')); % Explicit double

vx = double(ncread(velocity_nc, 'x')); % Explicit double
vy = double(ncread(velocity_nc, 'y')); % Explicit double
ud = double(ncread(velocity_nc, 'VX')); % Explicit double
vd = double(ncread(velocity_nc, 'VY')); % Explicit double

% -------------------------
% Interpolate to mesh nodes
% -------------------------
[nodes_trans, hd_on_nodes, ud_on_nodes, vd_on_nodes, original_nan_mask] = ...
    interpolateToNodes(tx, ty, hd, vx, vy, ud, vd, md, 2.15e6, 2.55e6);

% -------------------------
% Read edge points table
% -------------------------
edge_types = readtable('points_export.csv');

% -------------------------
% Match Points Exactly To Mesh Nodes And Save Indices
% -------------------------
% Extract point coordinates (explicit double)

% -------------------------
% Match Points To Mesh Nodes Using pdist2 And Save Indices
% -------------------------
% Extract point coordinates (explicit double)
px = double(edge_types.transformed_x(:));
py = double(edge_types.transformed_y(:));

% Mesh node coordinates (explicit double)
mx = double(nodes_trans(:,1));
my = double(nodes_trans(:,2));

nPts = numel(px);
nNodes = numel(mx);

% Build arrays
ptsXY   = [px, py];
nodesXY = [mx, my];

% Match points to mesh nodes using pdist2
[dist, idx] = pdist2(nodesXY, ptsXY, 'fasteuclidean', 'Smallest', 1);

% idx is K-by-nPts for 'Smallest',K. With K==1 it is 1-by-nPts -> make it a column vector
idx = idx(:);

% Ensure table columns are numeric/logical and same length as points
isGroundingPoint = logical(edge_types.grounding_line(:));   % 1 = grounding
isIceEdgePoint   = logical(edge_types.ice_edge(:));         % 1 = ice edge

% Map point-level flags to node indices via idx
grounding_idx = unique(idx(isGroundingPoint & isfinite(idx)));
iceedge_idx   = unique(idx(isIceEdgePoint   & isfinite(idx)));

% Make them column vectors of doubles (safe for ismember etc.)
grounding_idx = double(grounding_idx(:));
iceedge_idx   = double(iceedge_idx(:));


% Create mesh
nv = md.mesh.numberofvertices;
ne = md.mesh.numberofelements;
md.mask.ice_levelset = -1 * ones(nv,1);
size(iceedge_idx)
md.mask.ice_levelset(iceedge_idx) = 1;
md.mask.ice_levelset = reinitializelevelset(md, md.mask.ice_levelset); % turns ice_levelset into a signed distance
hVertices = 2500 * ones(md.mesh.numberofvertices, 1); % default 2.5km resolution for each node 
hVertices(abs(md.mask.ice_levelset)<5e3) = 500; % refine elements near the ice front to 500m resolution
md = bamg(md, 'hVertices', hVertices, 'hmin', 400, 'hmax', 2500);

[nodes_trans, hd_on_nodes, ud_on_nodes, vd_on_nodes, original_nan_mask] = ...
    interpolateToNodes(tx, ty, hd, vx, vy, ud, vd, md, 2.15e6, 2.55e6);
nv = md.mesh.numberofvertices;
ne = md.mesh.numberofelements;

% Plotting (optional)
if doPlot
    elems = double(md.mesh.elements); % Explicit double
    if size(elems, 1) == 3 && size(elems, 2) ~= 3
        elems = elems';
    end
    raw_x = double(edge_types.transformed_x(:));
    raw_y = double(edge_types.transformed_y(:));
    plotInterpolatedFields(elems, nodes_trans, hd_on_nodes, sqrt(ud_on_nodes.^2 + vd_on_nodes.^2), ...
        original_nan_mask, grounding_idx, iceedge_idx);
end



% -------------------------
% Model setup
% -------------------------
m = md.materials;
md.geometry.thickness = hd_on_nodes;
md.geometry.surface  = hd_on_nodes * m.rho_ice / m.rho_water;
md.geometry.base     = hd_on_nodes * (m.rho_ice - m.rho_water) / m.rho_water;
md = setmask(md, '', '');
md = setflowequation(md, 'SSA', 'all');
md.cluster = generic('np', 2);
md.materials.rheology_B = 1.8e8 * ones(md.mesh.numberofelements, 1);
md.friction.coefficient = zeros(md.mesh.numberofvertices, 1);
md.miscellaneous.name = 'waow';

md.initialization.vx = ud_on_nodes;
md.initialization.vy = vd_on_nodes;
md.initialization.vel = sqrt(ud_on_nodes.^2 + vd_on_nodes.^2);

% Inversion setup
md.inversion = m1qn3inversion();
md.inversion.vx_obs = ud_on_nodes;
md.inversion.vy_obs = vd_on_nodes;
md.inversion.vel_obs = sqrt(ud_on_nodes.^2 + vd_on_nodes.^2);

% Inversion parameters
md.inversion.iscontrol = 0;
md.inversion.control_parameters = {'MaterialsRheologyBbar'};
md.inversion.maxsteps = 20;
md.inversion.cost_functions = [101, 502];
md.inversion.cost_functions_coefficients = ones(md.mesh.numberofvertices, 2);
md.inversion.min_parameters = cuffey(273) * ones(md.mesh.numberofelements, 1);
md.inversion.max_parameters = cuffey(200) * ones(md.mesh.numberofelements, 1);

% Initialize arrays (pre-allocate)
md.stressbalance.spcvx = ud_on_nodes;
md.stressbalance.spcvx(iceedge_idx) = NaN;
md.stressbalance.spcvy = vd_on_nodes;
md.stressbalance.spcvy(iceedge_idx) = NaN;
md.stressbalance.spcvz = nan(nv, 1);
md.stressbalance.loadingforce = zeros(nv, 3);
md.stressbalance.referential = nan(nv, 6);
md.friction.p = zeros(ne, 1);
md.friction.q = zeros(ne, 1);
md.mask.ocean_levelset = -1 * ones(nv, 1);


% Solve
md.verbose = verbose(0);
md = solve(md, 'Stressbalance');
plotmodel(md, 'axis#all', 'tight', ...
    'data', md.results.StressbalanceSolution.MaterialsRheologyBbar, 'caxis', [1.3 1.9]*10^8, 'title', 'inferred B', ...
    'data', md.results.StressbalanceSolution.Vel, 'title', 'modeled velocities');

% -------------------------
% Interpolation function
% -------------------------
function [nodes_trans, hd_nodes, ud_nodes_filled, vd_nodes_filled, original_nan_mask] = ...
        interpolateToNodes(tx, ty, hd, vx, vy, ud, vd, md, c1, c2)
    % Explicit double for all inputs
    tx = double(tx);
    ty = double(ty);
    hd = double(hd);
    vx = double(vx);
    vy = double(vy);
    ud = double(ud);
    vd = double(vd);

    % Sort and reorder data
    [ty_s, iy] = sort(ty, 'ascend');
    [tx_s, ix] = sort(tx, 'ascend');
    hd_s = hd(iy, ix);

    [vy_s, iyv] = sort(vy, 'ascend');
    [vx_s, ixv] = sort(vx, 'ascend');
    ud_s = ud(iyv, ixv);
    vd_s = vd(iyv, ixv);

    % Build gridded interpolants
    F_hd = griddedInterpolant({ty_s, tx_s}, hd_s, 'linear', 'nearest');
    F_ud = griddedInterpolant({vy_s, vx_s}, ud_s, 'linear', 'nearest');
    F_vd = griddedInterpolant({vy_s, vx_s}, vd_s, 'linear', 'nearest');

    % Transform mesh nodes (explicit double)
    nodes_trans = [-double(md.mesh.y(:)) + c1, -double(md.mesh.x(:)) + c2];

    % Query interpolants
    yq = nodes_trans(:, 2);
    xq = nodes_trans(:, 1);
    hd_nodes = F_hd(yq, xq);
    ud_nodes = F_ud(yq, xq);
    vd_nodes = F_vd(yq, xq);

    % Record original NaN mask
    original_nan_mask = ~(isfinite(ud_nodes) & isfinite(vd_nodes));

    % Fill NaNs
    valid_any = isfinite(ud_nodes) | isfinite(vd_nodes);
    ud_nodes_filled = ud_nodes;
    vd_nodes_filled = vd_nodes;

    if any(~valid_any)
        F_ud_nodes = scatteredInterpolant(xq(valid_any), yq(valid_any), double(ud_nodes(valid_any)), 'nearest', 'nearest');
        F_vd_nodes = scatteredInterpolant(xq(valid_any), yq(valid_any), double(vd_nodes(valid_any)), 'nearest', 'nearest');

        nan_ud = ~isfinite(ud_nodes_filled);
        if any(nan_ud)
            ud_nodes_filled(nan_ud) = F_ud_nodes(xq(nan_ud), yq(nan_ud));
        end

        nan_vd = ~isfinite(vd_nodes_filled);
        if any(nan_vd)
            vd_nodes_filled(nan_vd) = F_vd_nodes(xq(nan_vd), yq(nan_vd));
        end
    end
end

% -------------------------
% Plotting function (updated)
% -------------------------
function plotInterpolatedFields(elems, nodes_trans, hd_nodes, speed_nodes, original_nan_mask, grounding_idx, iceedge_idx, raw_x, raw_y)
    TR = triangulation(double(elems), double(nodes_trans));
    figure('Name', 'Interpolated Fields on Transformed Mesh', 'Color', [1 1 1], 'Position', [100 100 1200 450]);
    
    % Subplot 1: thickness (1/3 width)
    subplot(1, 3, 1);
    hold on; axis equal; box on;
    triplot(TR, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6);
    scatter(double(nodes_trans(:,1)), double(nodes_trans(:,2)), 14, double(hd_nodes), 'filled', 'MarkerEdgeColor', 'none');
    colormap(gca, parula);
    colorbar;
    title('Thickness (hd) interpolated on nodes');
    xlabel('X (transformed)'); ylabel('Y (transformed)');
    hold off;
    
    % Subplot 2: speed with region coloring (2/3 width)
    subplot(1, 3, [2 3]);
    hold on; axis equal; box on;
    triplot(TR, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6);
    
    % Create masks for different node types
    valid_mask = isfinite(double(speed_nodes));
    all_nodes = (1:size(nodes_trans, 1))';
    
    % Normal nodes: valid, not grounding, not ice edge
    normal_mask = valid_mask & ~ismember(all_nodes, double(grounding_idx)) & ~ismember(all_nodes, double(iceedge_idx));
    
    % Grounding line nodes
    grounding_mask = ismember(all_nodes, double(grounding_idx)) & valid_mask;
    
    % Ice edge nodes
    iceedge_mask = ismember(all_nodes, double(iceedge_idx)) & valid_mask;
    
    % Plot normal nodes with speed color
    scatter(double(nodes_trans(normal_mask,1)), double(nodes_trans(normal_mask,2)), 16, double(speed_nodes(normal_mask)), 'filled', 'MarkerEdgeColor', 'none');
    
    % Plot grounding line nodes as red X
    scatter(double(nodes_trans(grounding_mask,1)), double(nodes_trans(grounding_mask,2)), 50, 'r', 'x', 'LineWidth', 1.2);
    
    % Plot ice edge nodes as orange X
    scatter(double(nodes_trans(iceedge_mask,1)), double(nodes_trans(iceedge_mask,2)), 50, [1 0.65 0], 'x', 'LineWidth', 1.2);
    
    if nargin >= 9 && ~isempty(raw_x) && ~isempty(raw_y)
        % Plot raw points as small circles with a cyan fill and black edge
        scatter(double(raw_x), double(raw_y), 36, ...
            'MarkerEdgeColor', [1 1 0], 'MarkerFaceColor', [0 1 1], ...
            'LineWidth', 0.6, 'DisplayName', 'Raw CSV points');
    end

    % Add legend
    legend('Mesh', 'Normal Nodes (speed)', 'Grounding Line (red X)', 'Ice Edge (orange X)', 'Raw CSV points', 'Location', 'bestoutside');
    
    colormap(gca, parula);
    hcb = colorbar;
    hcb.Label.String = 'Speed (units)';
    title('Velocity Speed (grounding=red X, ice edge=orange X)');
    xlabel('X (transformed)'); ylabel('Y (transformed)');
    
    xmin = min(double(nodes_trans(:,1))); xmax = max(double(nodes_trans(:,1)));
    ymin = min(double(nodes_trans(:,2))); ymax = max(double(nodes_trans(:,2)));
    padx = 0.02*(xmax-xmin); pady = 0.02*(ymax-ymin);
    xlim([xmin-padx, xmax+padx]);
    ylim([ymin-pady, ymax+pady]);
    
    hold off;
end