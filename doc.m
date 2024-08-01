%% torus_homeo_sfc_isotropic_splg
%
% Function to isotropically sample a given torus-homeomorphic surface.
%
% Author : nicolas.douillet9 (at) gmail.com, 2017-2024.
%
%% Syntax
%
% torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z);
%
% torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z, range_u, range_v);
%
% torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z, range_u, range_v, option_random_sampling);
%
%% Description
%
% torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z)
% generates a tricolon vector / [60*120 3] matrix of X, Y, and Z coordinates
% of points sampling the -torus-homeomorphic- surface defined by
% function handles fctn_x, fctn_y, and fctn_z.
%
% torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z, range_u, range_v)
% generates range_u(1,3)*range_v(1,3) samples located in the area
% [min(u), max(u)] x [min(v) max(v)] = [range_u(1,1), range_u(1,2)] x [range_v(1,1) range_v(1,2)]
%
% torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z, range_u, range_v, option_random_sampling)
% randoms the sampling if option_random_sampling = true/1,
% else -option_random_sampling = false/0- sampling is uniform.
%
%% See also
%
% <https://fr.mathworks.com/matlabcentral/fileexchange/64284-sphere-homeomorphic-surface-quasi-isotropic-sampling sphere_homeo_sfc_isotropic_splg> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/85173-mesh-generation-toolbox?s_tid=srchtitle mesh_generation_toolbox> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/69212-geoid geoid> |
% <https://fr.mathworks.com/help/matlab/ref/rand.html rand> | 
% <https://fr.mathworks.com/help/matlab/ref/mesh.html mesh> | 
% <https://fr.mathworks.com/help/matlab/ref/trimesh.html trimesh>
%
%% Input arguments
%
% - fctn_x : function handle in x direction, in spherical coordinates, assumed overloaded for vectors and matrices
%
% - fctn_y : function handle in y direction, in spherical coordinates, assumed overloaded for vectors and matrices
%
% - fctn_z : function handle in z direction, in spherical coordinates, assumed overloaded for vectors and matrices
%
% - range_u : real row vector double, u parameter vector of type : [min(u), max(u), nb_samples_u].
%
% - range_v : real row vector double, v parameter vector of type : [min(v), max(v), nb_samples_v].
%
% - option_random_sampling : logical *true (1) / false (0).
%
%% Output arguments
%
%        [| | |]
% - M = [X Y Z], real matrix double of data. Size(M) = [nb_samples_u*nb_samples_v 3].
%        [| | |]
%
% - u : real matrix double, the sampling matrix / grid in u direction. Size(u) = [nb_samples_u,nb_samples_v].
%
% - v : real matrix double, the sampling matrix / grid in v direction. Size(v) = [nb_samples_u,nb_samples_v].
%
%        [|  |  | ]
% - T = [i1 i2 i3], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3],
%        [|  |  | ]
%
%                   with nb_triangles = nb_samples_u*(nb_samples_v-2).
%                   T is relevant only in the case option_random_sampling = false/0.
%
%% Example #1
% Isotropic random sampling
a = 1;
b = 1;
Rho = 8;
r = 2;

fctn_x = @(u,v)a*(Rho*cos(v)+r*sin(u).*cos(v));
fctn_y = @(u,v)b*(Rho*sin(v)+r*sin(u).*sin(v));
fctn_z = @(u,v)r*cos(u);

range_u = [0 2*pi 40]; % latitude angle for a torus
range_v = [0 2*pi 40]; % longitude angle for a torus

[M1, u1, v1] = torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z, range_u, range_v);
figure;
plot3(M1(:,1), M1(:,2), M1(:,3), 'b.'), hold on;
axis equal, axis tight;
colormap([0 0 1]);

%% Example #2
% Isotropic regular sampling
[M2, u2, v2, T] = torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z, range_u, range_v, 0);
TRI = triangulation(T, M2(:,1), M2(:,2), M2(:,3));

figure;
trimesh(TRI), hold on;
axis equal, axis tight;
colormap([0 0 1]);