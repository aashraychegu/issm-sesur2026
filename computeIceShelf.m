function ice_shelf_idx = computeIceShelf(md, tx, ty, mask, radius)
    % COMPUTEICESHELF - Identify ice shelf points in a mesh based on proximity to ocean (mask == 0).
    %
    % Inputs:
    %   md:        ISSM model structure containing mesh information.
    %   tx:        1D array of x-coordinates for the mask grid.
    %   ty:        1D array of y-coordinates for the mask grid.
    %   mask:      2D array of mask values (0 = ocean, 1 = ice).
    %   radius:    Radius of the square grid to check for ocean proximity (default: 6).
    %
    % Output:
    %   ice_shelf_idx: indices of points on the ice shelf

    % Default radius
    if nargin < 5
        radius = 6;
    end

    % Extract mesh coordinates
    x_edge = md.mesh.x;
    y_edge = md.mesh.y;

    % Reshape tx and ty as column vectors
    tx_col = tx(:);  % Nx1
    ty_col = ty(:);  % Nx1

    % Reshape x_edge and y_edge as row vectors
    x_edge_row = x_edge.';
    y_edge_row = y_edge.';

    % Compute absolute differences (NxM)
    diff_x = abs(tx_col - x_edge_row);
    diff_y = abs(ty_col - y_edge_row);

    % Find indices of minimum differences (1xM)
    [~, idx_x] = min(diff_x, [], 1);
    [~, idx_y] = min(diff_y, [], 1);

    % Transpose to match edge point order (Mx1)
    idx_x = idx_x.';
    idx_y = idx_y.';

    % Initialize ice_shelf
    ice_shelf = zeros(size(x_edge));

    % Loop over each mesh point
    for i = 1:length(x_edge)
        % Get the closest (x, y) indices in the mask grid
        cx = idx_x(i);
        cy = idx_y(i);

        % Define the square grid boundaries (clamped to mask dimensions)
        x_min = max(1, cx - radius);
        x_max = min(size(mask, 2), cx + radius);
        y_min = max(1, cy - radius);
        y_max = min(size(mask, 1), cy + radius);

        % Extract the submatrix
        submask = mask(y_min:y_max, x_min:x_max);

        % Check if any value in the submatrix is 0 (ocean)
        if any(submask(:) == 0)
            ice_shelf(i) = 1;
        end
    end
    ice_shelf_idx = find(ice_shelf == 1);
end