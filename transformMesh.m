function md = transformMesh(md, c1, c2)
% TRANSFORMMESH  Apply a coordinate transform to an ISSM mesh
%
%   md = TRANSFORMMESH(md, c1, c2) replaces the mesh node coordinates in
%   the ISSM model structure `md` by performing the transform
%       x_new = -y_old + c1
%       y_new = -x_old + c2
%   Inputs:
%     md  - ISSM model struct with fields md.mesh.x and md.mesh.y
%     c1  - scalar offset applied to transformed x coordinates
%     c2  - scalar offset applied to transformed y coordinates
%
%   Output:
%     md  - model struct with updated md.mesh.x and md.mesh.y
%
%   Example:
%     md = transformMesh(md, 2.15e6, 2.55e6);

% Read original node coordinates and apply transform
o_mx = double(md.mesh.x(:));
o_my = double(md.mesh.y(:));

md.mesh.x = -o_my + double(c1);
md.mesh.y = -o_mx + double(c2);
end