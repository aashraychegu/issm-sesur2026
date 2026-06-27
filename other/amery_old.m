% amery.m
addpath('C:\Users\aashr\ISSM-Windows-MATLAB\bin');

% Toggle plotting on/off
doPlot = false;

% Load ISSM model and any workspace variables
load('data_pinns_Amery.mat');
md = model;

% Create mesh / triangulation if needed
md = bamg(model, 'domain', 'DomainOutline.exp', 'hmin', 2500);
md.mesh; % ensure mesh is created/filled

% -------------------------
% Read NetCDF files (thickness and velocity)
% -------------------------
thickness_nc = 'ice_thickness.nc';
velocity_nc  = 'antarctica_ice_velocity_1995-2001_450m_v01.1.nc';

tx = double(ncread(thickness_nc, 'x'));         % tx: x coords for thickness file
ty = double(ncread(thickness_nc, 'y'));         % ty: y coords for thickness file
hd = double(ncread(thickness_nc, 'thickness')); % hd: [ny, nx] (rows=ty, cols=tx)
mesh = double(ncread(thickness_nc, 'mask'));

vx = double(ncread(velocity_nc, 'x'));          % vx: x coords for velocity file
vy = double(ncread(velocity_nc, 'y'));          % vy: y coords for velocity file
ud = double(ncread(velocity_nc, 'VX'));         % ud: [ny_v, nx_v]
vd = double(ncread(velocity_nc, 'VY'));         % vd: [ny_v, nx_v]

% Transformation constants (used elsewhere)
c1 = 2.15e6;
c2 = 2.55e6;

% -------------------------
% Interpolate to mesh nodes (returns transformed nodes and filled fields)
% -------------------------
[nodes_trans, hd_on_nodes, ud_on_nodes, vd_on_nodes, original_nan_mask] = ...
    interpolateToNodes(tx, ty, hd, vx, vy, ud, vd, md, c1, c2);

% compute speed for later use
speed_nodes = sqrt(ud_on_nodes.^2 + vd_on_nodes.^2);

% -------------------------
% Plotting (separate function)
% -------------------------
if doPlot
    elems = md.mesh.elements;
    if size(elems,1) == 3 && size(elems,2) ~= 3
        elems = elems';
    end
    plotInterpolatedFields(elems, nodes_trans, hd_on_nodes, speed_nodes, original_nan_mask);
end



% -------------------------
% Continue with original workflow (use interpolated values where appropriate)
% -------------------------
md.geometry.thickness = hd_on_nodes;
md.geometry.surface = hd_on_nodes * 917 / 1028;
md.geometry.base = hd_on_nodes * -111 / 1028;
md = setmask(md,'','');
md = setflowequation(md,'SSA','all');
md.cluster = generic('np',2);
md.materials.rheology_B = 1.8e8 * ones(md.mesh.numberofelements, 1);
md.friction.coefficient = zeros(md.mesh.numberofvertices,1);
md.miscellaneous.name = 'wow';

% results of previous run are taken as observations
md.inversion = m1qn3inversion();
md.inversion.vx_obs = ud_on_nodes;
md.inversion.vy_obs = vd_on_nodes;
md.inversion.vel_obs = sqrt(ud_on_nodes.^2 + vd_on_nodes.^2);

edge_types = readtable('points_export.csv');

maxsteps = 20;
md.inversion.iscontrol = 1;
md.inversion.control_parameters = {'MaterialsRheologyBbar'};
md.inversion.maxsteps = maxsteps;
md.inversion.cost_functions = [101,502];
md.inversion.cost_functions_coefficients = ones(md.mesh.numberofvertices,2);
md.inversion.min_parameters = cuffey(273)*ones(md.mesh.numberofelements,1);
md.inversion.max_parameters = cuffey(200)*ones(md.mesh.numberofelements,1);
md.stressbalance.spcvx = nan*ones(md.mesh.numberofvertices, 1);
md.stressbalance.spcvy = nan*ones(md.mesh.numberofvertices, 1);
md.stressbalance.spcvz = nan*ones(md.mesh.numberofvertices, 1);
md.stressbalance.loadingforce = zeros(md.mesh.numberofvertices, 3);
md.stressbalance.referential = nan*ones(md.mesh.numberofvertices, 6);
md.friction.p = zeros(md.mesh.numberofelements,1);
md.friction.q = zeros(md.mesh.numberofelements,1);
md.mask.ocean_levelset = zeros(md.mesh.numberofvertices,1);


% Go solve!
md.verbose = verbose(0);
md = solve(md,'Stressbalance');
plotmodel(md,'axis#all','tight','data',md.results.StressbalanceSolution.MaterialsRheologyBbar,'caxis',[1.3 1.9]*10^8,'title','inferred B',...
    'data',md.results.StressbalanceSolution.Vel,'title','modeled velocities');

function [nodes_trans, hd_nodes_filled, ud_nodes_filled, vd_nodes_filled, original_nan_mask] = ...
        interpolateToNodes(tx, ty, hd, vx, vy, ud, vd, md, c1, c2)

    % Convert to double (defensive)
    tx = double(tx); ty = double(ty); hd = double(hd);
    vx = double(vx); vy = double(vy); ud = double(ud); vd = double(vd);

    % Sort grid vectors ascending and reorder data to match {y,x} ordering
    % hd: rows correspond to ty, columns to tx
    [ty_s, iy] = sort(ty, 'ascend');
    [tx_s, ix] = sort(tx, 'ascend');
    hd_s = hd(iy, ix);

    % ud/vd: rows correspond to vy, columns to vx
    [vy_s, iyv] = sort(vy, 'ascend');
    [vx_s, ixv] = sort(vx, 'ascend');
    ud_s = ud(iyv, ixv);
    vd_s = vd(iyv, ixv);

    % Build gridded interpolants using {y, x} ordering. Use 'none' extrap so we can detect missing values.
    F_hd = griddedInterpolant({ty_s, tx_s}, hd_s, 'linear', 'nearest');
    F_ud = griddedInterpolant({vy_s, vx_s}, ud_s, 'linear', 'nearest');
    F_vd = griddedInterpolant({vy_s, vx_s}, vd_s, 'linear', 'nearest');

    % Transform mesh nodes (only nodes)
    mesh_x = md.mesh.x(:);
    mesh_y = md.mesh.y(:);
    nodes_trans = [-mesh_y + c1, -mesh_x + c2]; % columns [x_trans, y_trans]

    % Query interpolants (note griddedInterpolant signature is F(y,x))
    yq = nodes_trans(:,2);
    xq = nodes_trans(:,1);

    hd_nodes = F_hd(yq, xq);
    ud_nodes = F_ud(yq, xq);
    vd_nodes = F_vd(yq, xq);

    % Record which nodes were NaN originally (before any fill)
    original_nan_mask = ~(isfinite(ud_nodes) & isfinite(vd_nodes));

    % Fill NaNs in ud/vd using nearest-neighbor from any node with a valid component
    xn = xq; yn = yq;
    valid_u = isfinite(ud_nodes);
    valid_v = isfinite(vd_nodes);
    valid_any = valid_u | valid_v;

    ud_nodes_filled = ud_nodes;
    vd_nodes_filled = vd_nodes;

    if any(~valid_any)
        % Use nearest neighbor scattered interpolant built from nodes that have any valid component.
        F_ud_nodes = scatteredInterpolant(xn(valid_any), yn(valid_any), ud_nodes(valid_any), 'nearest', 'nearest');
        F_vd_nodes = scatteredInterpolant(xn(valid_any), yn(valid_any), vd_nodes(valid_any), 'nearest', 'nearest');

        nan_ud = ~isfinite(ud_nodes_filled);
        if any(nan_ud)
            ud_nodes_filled(nan_ud) = F_ud_nodes(xn(nan_ud), yn(nan_ud));
        end

        nan_vd = ~isfinite(vd_nodes_filled);
        if any(nan_vd)
            vd_nodes_filled(nan_vd) = F_vd_nodes(xn(nan_vd), yn(nan_vd));
        end
    end

    hd_nodes_filled = hd_nodes;
end

function plotInterpolatedFields(elems, nodes_trans, hd_nodes, speed_nodes, original_nan_mask)
    % elems: connectivity NTx3
    % nodes_trans: Nx2 [x_trans, y_trans]
    % hd_nodes, speed_nodes: N x 1
    % original_nan_mask: N x 1 logical marking nodes originally NaN before fill

    TR = triangulation(elems, nodes_trans);

    figure('Name','Interpolated Fields on Transformed Mesh','Color',[1 1 1], 'Position',[100 100 1000 450]);

    % Subplot 1: thickness
    subplot(1,2,1);
    hold on; axis equal; box on;
    triplot(TR, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6);
    scatter(nodes_trans(:,1), nodes_trans(:,2), 14, hd_nodes, 'filled', 'MarkerEdgeColor','none');
    colormap(gca,parula);
    colorbar;
    title('Thickness (hd) interpolated on nodes');
    xlabel('X (transformed)'); ylabel('Y (transformed)');
    hold off;

    % Subplot 2: speed with original NaNs highlighted
    subplot(1,2,2);
    hold on; axis equal; box on;
    triplot(TR, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6);

    valid_mask = isfinite(speed_nodes);
    scatter(nodes_trans(valid_mask,1), nodes_trans(valid_mask,2), 16, speed_nodes(valid_mask), 'filled', 'MarkerEdgeColor','none');

    if any(original_nan_mask)
        scatter(nodes_trans(original_nan_mask,1), nodes_trans(original_nan_mask,2), 32, 'r', 'x', 'LineWidth', 1.2);
    end

    colormap(gca,parula);
    hcb = colorbar;
    hcb.Label.String = 'Speed (units)';
    title('Velocity Speed (original NaNs in red)');
    xlabel('X (transformed)'); ylabel('Y (transformed)');

    xmin = min(nodes_trans(:,1)); xmax = max(nodes_trans(:,1));
    ymin = min(nodes_trans(:,2)); ymax = max(nodes_trans(:,2));
    padx = 0.02*(xmax-xmin); pady = 0.02*(ymax-ymin);
    xlim([xmin-padx, xmax+padx]);
    ylim([ymin-pady, ymax+pady]);

    hold off;
end
