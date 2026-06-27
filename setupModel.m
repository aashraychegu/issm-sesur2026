function md = setupModel(md, hd_on_nodes, ud_on_nodes, vd_on_nodes, iceedge_idx, verbose_level)
    % SETUPMODEL - Set up and solve an ISSM model with inversion for rheology parameter B.
    %
    % Inputs:
    %   md:             ISSM model structure.
    %   hd_on_nodes:    Ice thickness interpolated to mesh nodes.
    %   ud_on_nodes:    X-velocity interpolated to mesh nodes.
    %   vd_on_nodes:    Y-velocity interpolated to mesh nodes.
    %   iceedge_idx:    Indices of ice shelf edge nodes.
    %   verbose_level:  Verbosity level for ISSM (default: 0).
    %
    % Output:
    %   md:             Updated ISSM model structure after solving.

    % Set default verbosity level
    if nargin < 6
        verbose_level = 0;
    end

    % Extract materials and mesh properties
    m = md.materials;
    nv = md.mesh.numberofvertices;
    ne = md.mesh.numberofelements;

    % -------------------------
    % Model setup
    % -------------------------
    md.geometry.thickness = hd_on_nodes;
    md.geometry.surface = hd_on_nodes * m.rho_ice / m.rho_water;
    md.geometry.base = hd_on_nodes * (m.rho_ice - m.rho_water) / m.rho_water;

    % Reset mask
    md = setmask(md, '', '');

    % Set flow equation to SSA for all elements
    md = setflowequation(md, 'SSA', 'all');

    % Cluster configuration
    md.cluster = generic('np', 2);

    % Set rheology parameter B
    md.materials.rheology_B = 1.8e8 * ones(ne, 1);

    % Set friction coefficient to zero
    md.friction.coefficient = zeros(nv, 1);

    % Set model name
    md.miscellaneous.name = 'waow';

    % Initialize velocities
    md.initialization.vx = ud_on_nodes;
    md.initialization.vy = vd_on_nodes;
    md.initialization.vel = sqrt(ud_on_nodes.^2 + vd_on_nodes.^2);

    md.inversion = m1qn3inversion();
    md.inversion.vx_obs = ud_on_nodes;
    md.inversion.vy_obs = vd_on_nodes;
    md.inversion.vel_obs = sqrt(ud_on_nodes.^2 + vd_on_nodes.^2);

    % Inversion parameters
    md.inversion.control_parameters = {'MaterialsRheologyBbar'};
    md.inversion.maxsteps = 20;
    md.inversion.cost_functions = [101, 502];
    md.inversion.cost_functions_coefficients = ones(nv, 2);
    md.inversion.min_parameters = cuffey(273) * ones(ne, 1);
    md.inversion.max_parameters = cuffey(200) * ones(ne, 1);

    % -------------------------
    % Initialize arrays
    % -------------------------
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

    md.verbose = verbose(verbose_level);
end