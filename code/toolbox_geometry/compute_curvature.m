function [Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(V,F,options)

% compute_curvature - compute principal curvature directions and values
%
%   [Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(V,F,options);
%
%   Umin is the direction of minimum curvature
%   Umax is the direction of maximum curvature
%   Cmin is the minimum curvature
%   Cmax is the maximum curvature
%   Cmean=(Cmin+Cmax)/2
%   Cgauss=Cmin*Cmax
%   Normal is the normal to the surface
%
%   options.curvature_smoothing controls the size of the ring used for
%       averaging the curvature tensor.
%
%   The algorithm is detailed in
%       David Cohen-Steiner and Jean-Marie Morvan.
%       Restricted Delaunay triangulations and normal cycle.
%       In Proc. 19th Annual ACM Symposium on Computational Geometry,
%       pages 237-246, 2003.
%   and also in
%       Pierre Alliez, David Cohen-Steiner, Olivier Devillers, Bruno Le?vy, and Mathieu Desbrun.
%       Anisotropic Polygonal Remeshing.
%       ACM Transactions on Graphics, 2003.
%       Note: SIGGRAPH '2003 Conference Proceedings
%
%   Copyright (c) 2007 Gabriel Peyre

orient = 1;

options.null = 0;
naver = getoptions(options, 'curvature_smoothing', 3);
verb = getoptions(options, 'verb', 1);

[V,F] = check_face_vertex(V,F);

n = size(V,2);
m = size(F,2);

% associate each edge to a pair of faces
Af = -triangulation2adjacency(F);

i = [F(1,:) F(2,:) F(3,:)];
j = [F(2,:) F(3,:) F(1,:)];
s = [1:m 1:m 1:m];
Af = sparse(i,j,s,n,n);

%% PATCH %%%
Af(Af>m) = 0;
%% END PATCH %%%

[i,j,s1] = find(Af);     % direct link
[i,j,s2] = find(Af');    % reverse link

I = find( (s1>0) + (s2>0) == 2 );

% links edge->faces
E = [s1(I) s2(I)];
i = i(I); j = j(I);
% only directed edges
I = find(i<j);
E = E(I,:);
i = i(I); j = j(I);
ne = length(i); % number of directed edges

% normalized edge
e = V(:,j) - V(:,i);
d = sqrt(sum(e.^2,1));
e = e ./ repmat(d,3,1);
% avoid too large numerics
% d = d./mean(d);

% normals to faces
[tmp,normal] = compute_normal(V,F);

% inner product of normals
dp = sum( normal(:,E(:,1)) .* normal(:,E(:,2)), 1 );
% angle un-signed
beta = acos(clamp(dp,-1,1));
% sign
cp = crossp( normal(:,E(:,1))', normal(:,E(:,2))' )';
si = orient * sign( sum( cp.*e,1 ) );
% angle signed
beta = beta .* si;
% tensors
T = zeros(3,3,ne);
for x=1:3
    for y=1:x
        T(x,y,:) = reshape( e(x,:).*e(y,:), 1,1,ne );
        T(y,x,:) = T(x,y,:);
    end
end

T = T.*repmat( reshape(d.*beta,1,1,ne), [3,3,1] );
% T = T.*repmat( reshape(beta,1,1,ne), [3,3,1] );

% do pooling on vertices
Tv = zeros(3,3,n);
w = zeros(1,1,n);
for k=1:ne
%    progressbar(k,ne);
    Tv(:,:,i(k)) = Tv(:,:,i(k)) + T(:,:,k);
    Tv(:,:,j(k)) = Tv(:,:,j(k)) + T(:,:,k);
    w(:,:,i(k)) = w(:,:,i(k)) + 1;
    w(:,:,j(k)) = w(:,:,j(k)) + 1;
end
w(w<eps) = 1;
Tv = Tv; % ./repmat(w,[3,3,1]);

%%
% Compute area around each V

% area of each face
a = V(:,F(3,:)) - V(:,F(1,:));
b = V(:,F(2,:)) - V(:,F(1,:));
ab = crossp(a',b')';
Af = sqrt(sum(ab.^2))/2;
% area of each vertex
m = size(F,2);
U = sparse( [1:m, 1:m, 1:m], [F(1,:) F(2,:) F(3,:)], [Af,Af,Af] );
Av = full(sum(U,1));
% normalize
Tv = Tv ./ repmat( reshape(Av,[1 1 n]), [3 3 1] );

% do averaging to smooth the field
options.niter_averaging = naver;
for x=1:3
    for y=1:3
        a = Tv(x,y,:);
        a = perform_mesh_smoothing(F,V,a(:),options);
        Tv(x,y,:) = reshape( a, 1,1,n );
    end
end

% extract eigenvectors and eigenvalues
U = zeros(3,3,n);
D = zeros(3,n);
for k=1:n
    if verb==1
        progressbar(k,n);
    end
    [u,d] = eig(Tv(:,:,k));
    d = real(diag(d));
    % sort acording to [norma,min curv, max curv]
    [tmp,I] = sort(abs(d));
    D(:,k) = d(I);
    U(:,:,k) = real(u(:,I));
end

Umin = squeeze(U(:,3,:));
Umax = squeeze(U(:,2,:));
Cmin = D(2,:)';
Cmax = D(3,:)';
Normal = squeeze(U(:,1,:));
Cmean = (Cmin+Cmax)/2;
Cgauss = Cmin.*Cmax;

% enforce than min<max
I = find(Cmin>Cmax);
Cmin1 = Cmin; Umin1 = Umin;
Cmin(I) = Cmax(I); Cmax(I) = Cmin1(I);
Umin(:,I) = Umax(:,I); Umax(:,I) = Umin1(:,I);

% try to re-orient the normals
normal = compute_normal(V,F);
s = sign( sum(Normal.*normal,1) );
Normal = Normal .* repmat(s, 3,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = crossp(x,y)
% x and y are (m,3) dimensional
z = x;
z(:,1) = x(:,2).*y(:,3) - x(:,3).*y(:,2);
z(:,2) = x(:,3).*y(:,1) - x(:,1).*y(:,3);
z(:,3) = x(:,1).*y(:,2) - x(:,2).*y(:,1);

end

function A = triangulation2adjacency(face,vertex)

% triangulation2adjacency - compute the adjacency matrix
%   of a given triangulation.
%
%   A = triangulation2adjacency(face);
% or for getting a weighted graph
%   A = triangulation2adjacency(face,vertex);
%
%   Copyright (c) 2005 Gabriel Peyre


[tmp,face] = check_face_vertex([],face);
f = double(face)';

A = sparse([f(:,1); f(:,1); f(:,2); f(:,2); f(:,3); f(:,3)], ...
           [f(:,2); f(:,3); f(:,1); f(:,3); f(:,1); f(:,2)], ...
           1.0);
% avoid double links
A = double(A>0);

end

function y = clamp(x,a,b)

% clamp - clamp a value
%
%   y = clamp(x,a,b);
%
% Default is [a,b]=[0,1].
%
%   Copyright (c) 2004 Gabriel Peyre

if nargin<2
    a = 0;
end
if nargin<3
    b = 1;
end

if iscell(x)
    for i=1:length(x)
        y{i} = clamp(x{i},a,b);
    end
    return;
end

y = max(x,a);
y = min(y,b);

end

function f = perform_mesh_smoothing(face,vertex,f,options)

% perform_mesh_smoothing - smooth a function defined on a mesh by averaging
%
%   f = perform_mesh_smoothing(face,vertex,f,options);
%
%   Smooth a function f on a width of options.niter_averaging vertices.
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
naver = getoptions(options, 'niter_averaging', 1);
type = getoptions(options, 'averaging_type', 'combinatorial');

if nargin<3
    f = [];
end
if isempty(f)
    f = vertex;
end
if size(f,1)<size(f,2)
    f = f';
end
[vertex,face] = check_face_vertex(vertex,face);

if size(f,2)>1
    for i=1:size(f,2)
        f(:,i) = perform_mesh_smoothing(face,vertex,f(:,i),options);
    end
    return;
end

n = max(face(:));

% compute normalized averaging matrix
if strcmp(type, 'combinatorial')
    %add diagonal
    W = triangulation2adjacency(face) + speye(n);
    D = spdiags(full(sum(W,2).^(-1)),0,n,n);
    W = D*W;
else
    options.normalize=1;
    W = compute_mesh_weight(vertex,face,type,options);
end

% do averaging to smooth the field
for k=1:naver
    f = W*f;
end

end
