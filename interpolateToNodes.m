function [hd_nodes, ud_nodes_filled, vd_nodes_filled, original_nan_mask] = ...
    interpolateToNodes(md, tx, ty, hd, vx, vy, ud, vd)
    % INTERPOLATETONODES - Interpolate thickness, velocity, and mask data to mesh nodes.
    %
    % Inputs:
    %   md:        ISSM model structure containing mesh information.
    %   tx, ty:    1D arrays of x and y coordinates for the input grids.
    %   hd:        2D array of ice thickness.
    %   vx, vy:    1D arrays of x and y coordinates for velocity grids.
    %   ud, vd:    2D arrays of x and y velocity components.
    %
    % Outputs:
    %   nodes:             Coordinates of mesh nodes (Nx2 array).
    %   hd_nodes:          Interpolated ice thickness at mesh nodes.
    %   ud_nodes_filled:   Interpolated and filled x-velocity at mesh nodes.
    %   vd_nodes_filled:   Interpolated and filled y-velocity at mesh nodes.
    %   original_nan_mask: Logical array indicating original NaN locations.

    % Explicit double for all inputs
    

    % Extract mesh node coordinates
    nodes = [md.mesh.x(:), md.mesh.y(:)];

    % Sort and reorder thickness data
    [ty_s, iy] = sort(ty, 'ascend');
    [tx_s, ix] = sort(tx, 'ascend');
    hd_s = hd(iy, ix);

    % Sort and reorder velocity data
    [vy_s, iyv] = sort(vy, 'ascend');
    [vx_s, ixv] = sort(vx, 'ascend');
    ud_s = ud(iyv, ixv);
    vd_s = vd(iyv, ixv);

    % Build gridded interpolants
    F_hd = griddedInterpolant({ty_s, tx_s}, hd_s, 'linear', 'nearest');
    F_ud = griddedInterpolant({vy_s, vx_s}, ud_s, 'linear', 'nearest');
    F_vd = griddedInterpolant({vy_s, vx_s}, vd_s, 'linear', 'nearest');
    
    % Query interpolants at mesh nodes
    yq = nodes(:, 2);
    xq = nodes(:, 1);
    hd_nodes = F_hd(yq, xq);
    ud_nodes = F_ud(yq, xq);
    vd_nodes = F_vd(yq, xq);
    % size(hd_nodes), size(ud_nodes), size(vd_nodes)
    % Record original NaN mask
    original_nan_mask = ~(isfinite(ud_nodes) & isfinite(vd_nodes));

    % Fill NaNs using nearest-neighbor interpolation
    valid_any = isfinite(ud_nodes) | isfinite(vd_nodes);
    ud_nodes_filled = ud_nodes;
    vd_nodes_filled = vd_nodes;

    if any(~valid_any)
        % Use valid points for scattered interpolation
        valid_xq = xq(valid_any);
        valid_yq = yq(valid_any);
        valid_ud = ud_nodes(valid_any);
        valid_vd = vd_nodes(valid_any);

        % Create scattered interpolants for filling NaNs
        F_ud_nodes = scatteredInterpolant(valid_xq, valid_yq, valid_ud, 'nearest', 'nearest');
        F_vd_nodes = scatteredInterpolant(valid_xq, valid_yq, valid_vd, 'nearest', 'nearest');

        % Fill NaNs in ud_nodes_filled
        nan_ud = ~isfinite(ud_nodes_filled);
        if any(nan_ud)
            ud_nodes_filled(nan_ud) = F_ud_nodes(xq(nan_ud), yq(nan_ud));
        end

        % Fill NaNs in vd_nodes_filled
        nan_vd = ~isfinite(vd_nodes_filled);
        if any(nan_vd)
            vd_nodes_filled(nan_vd) = F_vd_nodes(xq(nan_vd), yq(nan_vd));
        end
    end
end