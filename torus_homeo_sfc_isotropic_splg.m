function [M, u, v, T] = torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z, range_u, range_v, option_random_sampling);
%% torus_homeo_sfc_isotropic_splg : function to isotropically sample a given torus-homeomorphic surface.
%
% Author & support : nicolas.douillet (at) free.fr, 2017-2020.
%
%
% Syntax
%
% torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z);
% torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z, range_u, range_v);
% torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z, range_u, range_v,
% option_random_sampling);
%
%
% Description
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
% randoms the sampling if option_random_sampling = 1/true,
% else -option_random_sampling = false/0- sampling is uniform.
%
%
% See also : RAND, MESH, TRIMESH
%
%
% Input arguments
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
%
% Output arguments
%
%       [| | |]
% - M = [X Y Z], real matrix double of data. Size(M) = [nb_samples_u*nb_samples_v 3].
%       [| | |]
%
% - u : real matrix double, the sampling matrix / grid in u direction. Size(u) = [nb_samples_u,nb_samples_v].
%
% - v : real matrix double, the sampling matrix / grid in v direction. Size(v) = [nb_samples_u,nb_samples_v].
%
%       [|  |  | ]
% - T = [i1 i2 i3], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3],
%       [|  |  | ]
%
%                   with nb_triangles = nb_samples_u*(nb_samples_v-2).
%                   T is relevant only in the case option_random_sampling = false/0.
%
%
% Example
%
% Sampling a torus(Rho,r) centered on the origin
%
% a = 1;
% b = 1;
% Rho = 8;
% r = 2;
% 
% fctn_x = @(u,v)a*(Rho*cos(v)+r*sin(u).*cos(v));
% fctn_y = @(u,v)b*(Rho*sin(v)+r*sin(u).*sin(v));
% fctn_z = @(u,v)r*cos(u);
% 
% range_u = [0 2*pi 40]; % latitude angle for a torus
% range_v = [0 2*pi 40]; % longitude angle for a torus
%
% [M1, u1, v1] = torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z, range_u, range_v);
% [M2, u2, v2, T] = torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z, range_u, range_v, 0);
% TRI = triangulation(T, M2(:,1), M2(:,2), M2(:,3));
% figure;
% subplot(1,2,1);
% plot3(M1(:,1), M1(:,2), M1(:,3), 'b.'), hold on;
% axis equal, axis square, axis tight;
% colormap([0 0 1]);
% subplot(1,2,2);
% trimesh(TRI), hold on;
% axis equal, axis tight;
% colormap([0 0 1]);


%% Input parsing
assert(nargin >= 3, 'Not enough input arguments.');
assert(nargin < 7, 'Too many input arguments.');

if nargin < 6
    
    option_random_sampling = 1;
    disp('Random sampling option selected by default.');
    
    if nargin < 5
        
        range_v = [0 2*pi 120];
        disp('range_v missing ; default vector is [0 pi 120]');
        
        if nargin < 4
            
            range_u = [0 pi 60];
            disp('range_u missing ; default vector is [0 pi 60]');
            
        end
        
    end
    
end

% Check function handles
assert(isa(fctn_x,'function_handle'), 'fctn_x is not a function handle');
assert(isa(fctn_y,'function_handle'), 'fctn_y is not a function handle');
assert(isa(fctn_z,'function_handle'), 'fctn_z is not a function handle');

% Check (u, v) parameter ranges
assert(numel(range_u) >= 2, 'No range_u max value ; range_u must contain at least two elements.');

if size(range_u,2) < 2
    range_u = range_u';
end

if (size(range_u,2) < 3 || range_u(1,3) < 1)
    range_u = cat(2, range_u, 120);
    warning('range_u : missing or incorrect value for the number of samples ; set to 120 by default');
else
    assert(isreal(range_u(1,3)) && range_u(1,3) >= 1 && rem(range_u(1,3),1) == 0, 'range_u number of elements must be an integer.');        
end

u_min = range_u(1,1);
u_max = range_u(1,2);
nb_samples_u = range_u(1,3);
assert(u_max >= u_min, 'Error range_u : must be u_min <= u_max.');


assert(numel(range_v) >= 2, 'No range_v max value ; range_v must contain at least two elements.');

if size(range_v,2) < 2
    range_v = range_v';
end

if(size(range_v,2) < 3 || range_v(1,3) < 1)
    range_v = cat(2, range_v, 120);
    disp('range_v number of samples missing or incorrect value ; set to 120 by default');
else
    assert(isreal(range_v(1,3)) && range_v(1,3) >= 1 && rem(range_v(1,3),1) == 0, 'range_v number of elements must be an integer.');        
end


%% Body
v_min = range_v(1,1);
v_max = range_v(1,2);
nb_samples_v = range_v(1,3);
assert(v_max >= v_min, 'Error range_v : must be v_min <= v_max.');

% Surface sampling
[~, u_step] = meshgrid(linspace(0,1,nb_samples_v),linspace(0,1,nb_samples_u));
u_rand = rand(nb_samples_u, nb_samples_v); 

u = option_random_sampling*u_rand + (1-option_random_sampling)*u_step;
u = u_min + (u_max-u_min)*u;

v_step = meshgrid(linspace(0,1,nb_samples_v),linspace(0,1,nb_samples_u));
v_rand = rand(nb_samples_u, nb_samples_v);

v = option_random_sampling*v_rand + (1-option_random_sampling)*v_step;
v = v_min + (v_max - v_min)*v;

X = fctn_x(u,v);
Y = fctn_y(u,v);
Z = fctn_z(u,v);

% Triplet indices for mesh facets
T = zeros(nb_samples_u*nb_samples_v, 3);
row_idx = 1;
i = 1;

while(i < nb_samples_v)
    
    j = 1;
    
    while(j < nb_samples_u)       
        
        T(row_idx, :)   = [(i-1)*nb_samples_u+j (i-1)*nb_samples_u+j+1 i*nb_samples_u+j];
        row_idx = row_idx + 1;
        T(row_idx, :)   = [i*nb_samples_u+j i*nb_samples_u+j+1 (i-1)*nb_samples_u+j+1];
        row_idx = row_idx + 1;
        
        j = j + 1;
        
    end
    
    % begin-end loop junction
    T(row_idx, :)   = [(i-1)*nb_samples_u+j (i-1)*nb_samples_u+1 i*nb_samples_u+j];
    row_idx = row_idx + 1;
    T(row_idx, :)   = [i*nb_samples_u+j i*nb_samples_u (i-1)*nb_samples_u+1];
    row_idx = row_idx + 1;
    
    i = i + 1;
end

% Format X, Y, Z data into M matrix
assert(numel(X) == numel(Y) && numel(X) == numel(Z) && numel(Y) == numel(Z), 'X, Y, and Z numbers of elements mismatch.');
assert(size(X,1) == size(Y,1) && size(X,1) == size(Z,1) && size(Y,1) == size(Z,1), 'X, Y, and Z dimensions mismatch.');

if size(X,2) > 1
    X = X(:);
    Y = Y(:);
    Z = Z(:);
end

M = [X Y Z];


end % torus_homeo_sfc_isotropic_splg