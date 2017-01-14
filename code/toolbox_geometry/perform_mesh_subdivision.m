function [f1,face1] = perform_mesh_subdivision(f, face, nsub, options)

% perform_mesh_subdivision - perfrom a mesh sub-division
%
%   [f1,face1] = perform_mesh_subdivision(f, face, nsub, options);
%
%   face is a (3,nface) matrix of original face adjacency
%   face1 is the new matrix after subdivision
%   f is a (d,nvert) matrix containing the value f(:,i) of a function
%       at vertex i on the original mesh. One should have
%           nvert=max(face(:))
%       (can be multi dimensional like point position in R^3, d=3)
%   f1 is the value of the function on the subdivided mesh.
%
%   options.sub_type is the kind of subvision applied:
%       'linear4': 1:4 tolopoligical subivision with linear interpolation
%       'linear3': 1:3 tolopoligical subivision with linear interpolation
%       'loop': 1:4 tolopoligical subivision with loop interpolation
%       'butterfly': 1:4 tolopoligical subivision with linear interpolation
%       'sqrt3': 1:3 topological subdivision with sqrt(3) interpolation
%          (dual scheme).
%       'spherical4': 1:4 tolopoligical subivision with linear
%           interpolation and projection of f on the sphere
%       'spherical3': 1:3 tolopoligical subivision with linear
%           interpolation and projection of f on the sphere
%
%   An excellent reference for mesh subdivision is
%       Subdivision for Modeling and Animation,
%       SIGGRAPH 2000 Course notes.
%       http://mrl.nyu.edu/publications/subdiv-course2000/
%
%   The sqrt(3) subdivision is explained in
%       \sqrt{3}-subdivision, Leif Kobbelt
%       Proc. of SIGGRAPH 2000
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
if nargin<2
    error('Not enough arguments');
end
if nargin==2
    nsub=1;
end

sub_type = getoptions(options, 'sub_type', '1:4');
spherical = getoptions(options, 'spherical', 0);
sanity_check = getoptions(options, 'sanity_check', 1);

switch lower(sub_type)
    case 'linear3'
        interpolation = 'linear';
        topology = 3;
    case 'linear4'
        interpolation = 'linear';
        topology = 4;
    case 'loop'
        interpolation = 'loop';
        topology = 4;
    case 'butterfly'
        interpolation = 'butterfly';
        topology = 4;
    case 'sqrt3';
        interpolation = 'sqrt3';
        topology = 3;
    case 'spherical3'
        interpolation = 'linear';
        topology = 3;
        spherical = 1;
    case 'spherical4'
        interpolation = 'linear';
        topology = 4;
        spherical = 1;
    case '1:3'
        interpolation = 'linear';
        topology = 3;
    case '1:4'
        interpolation = 'linear';
        topology = 4;
end

if nsub==0
    f1 = f;
    face1 = face;
    return;
end

if nsub>1
    % special case for multi-subdivision
    f1 = f;
    face1 = face;
    for i = 1:nsub
         [f1,face1] = perform_mesh_subdivision(f1,face1,1, options);
    end
    return;
end


if size(f,1)>size(f,2) && sanity_check
    f=f';
end
if size(face,1)>size(face,2) && sanity_check
    face=face';
end

m = size(face,2);
n = size(f,2);

verb = getoptions(options, 'verb', n>500);
loop_weigths = getoptions(options, 'loop_weigths', 1);

if topology==3
    f1 = ( f(:,face(1,:)) + f(:,face(2,:)) + f(:,face(3,:)))/3;
    f1 = cat(2, f, f1 );
    %%%%%% 1:3 subdivision %%%%%
    switch interpolation
        case 'linear'
            face1 = cat(2, ...
                [face(1,:); face(2,:); n+(1:m)], ...
                [face(2,:); face(3,:); n+(1:m)], ...
                [face(3,:); face(1,:); n+(1:m)] );
        case 'sqrt3'
            face1 = [];
            edge = compute_edges(face);
            ne = size(edge,2);
            e2f = compute_edge_face_ring(face);
            face1 = [];
            % create faces
            for i=1:ne
                if verb
                    progressbar(i,n+ne);
                end
                v1 = edge(1,i); v2 = edge(2,i);
                F1 = e2f(v1,v2); F2 = e2f(v2,v1);
                if min(F1,F2)<0
                    % special case
                    face1(:,end+1) = [v1 v2 n+max(F1,F2)];
                else
                    face1(:,end+1) = [v1 n+F1 n+F2];
                    face1(:,end+1) = [v2 n+F2 n+F1];
                end
            end
            % move old vertices
            vring0 = compute_vertex_ring(face);
            for k=1:n
                if verb
                    progressbar(k+ne,n+ne);
                end
                m = length(vring0{k});
               	beta = (4-2*cos(2*pi/m))/(9*m);         % warren weights
                f1(:,k) = f(:,k)*(1-m*beta) + beta*sum(f(:,vring0{k}),2);
            end

        otherwise
            error('Unknown scheme for 1:3 subdivision');
    end
else
    %%%%%% 1:4 subdivision %%%%%
    i = [face(1,:) face(2,:) face(3,:) face(2,:) face(3,:) face(1,:)];
    j = [face(2,:) face(3,:) face(1,:) face(1,:) face(2,:) face(3,:)];
    I = find(i<j);
    i = i(I); j = j(I);
    [tmp,I] = unique(i + 1234567*j);
    i = i(I); j = j(I);
    ne = length(i); % number of edges
    s = n+(1:ne);

    A = sparse([i;j],[j;i],[s;s],n,n);

    % first face
    v12 = full( A( face(1,:) + (face(2,:)-1)*n ) );
    v23 = full( A( face(2,:) + (face(3,:)-1)*n ) );
    v31 = full( A( face(3,:) + (face(1,:)-1)*n ) );

    face1 = [   cat(1,face(1,:),v12,v31),...
        cat(1,face(2,:),v23,v12),...
        cat(1,face(3,:),v31,v23),...
        cat(1,v12,v23,v31)   ];


    switch interpolation
        case 'linear'
            % add new vertices at the edges center
            f1 = [f, (f(:,i)+f(:,j))/2 ];

        case 'butterfly'

            global vring e2f fring facej;
            vring = compute_vertex_ring(face1);
            e2f = compute_edge_face_ring(face);
            fring = compute_face_ring(face);
            facej = face;
            f1 = zeros(size(f,1),n+ne);
            f1(:,1:n) = f;
            for k=n+1:n+ne
                if verb
                    progressbar(k-n,ne);
                end
                [e,v,g] = compute_butterfly_neighbors(k, n);
                f1(:,k) = 1/2*sum(f(:,e),2) + 1/8*sum(f(:,v),2) - 1/16*sum(f(:,g),2);
            end

        case 'loop'

            global vring e2f fring facej;
            vring = compute_vertex_ring(face1);
            vring0 = compute_vertex_ring(face);
            e2f = compute_edge_face_ring(face);
            fring = compute_face_ring(face);
            facej = face;
            f1 = zeros(size(f,1),n+ne);
            f1(:,1:n) = f;
            % move old vertices
            for k=1:n
                if verb
                    progressbar(k,n+ne);
                end
                m = length(vring0{k});
                if loop_weigths==1
                    beta = 1/m*( 5/8 - (3/8+1/4*cos(2*pi/m))^2 );   % loop original construction
                else
                    beta = 3/(8*m);         % warren weights
                end
                f1(:,k) = f(:,k)*(1-m*beta) + beta*sum(f(:,vring0{k}),2);
            end
            % move new vertices
            for k=n+1:n+ne
                if verb
                    progressbar(k,n+ne);
                end
                [e,v] = compute_butterfly_neighbors(k, n);
                f1(:,k) = 3/8*sum(f(:,e),2) + 1/8*sum(f(:,v),2);
            end

        otherwise
            error('Unknown scheme for 1:3 subdivision');
    end
end

if spherical
    % project on the sphere
    d = sqrt( sum(f1.^2,1) );
    d(d<eps)=1;
    f1 = f1 ./ repmat( d, [size(f,1) 1]);
end

end


function vring = compute_vertex_ring(face)

% compute_vertex_ring - compute the 1 ring of each vertex in a triangulation.
%
%   vring = compute_vertex_ring(face);
%
%   vring{i} is the set of vertices that are adjacent
%   to vertex i.
%
%   Copyright (c) 2004 Gabriel Peyre

[tmp,face] = check_face_vertex([],face);

nverts = max(max(face));

A = triangulation2adjacency(face);
[i,j,s] = find(sparse(A));

% create empty cell array
vring{nverts} = [];

for m = 1:length(i)
    vring{i(m)}(end+1) = j(m);
end

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

return;


nvert = max(max(face));
nface = size(face,1);
A = spalloc(nvert,nvert,3*nface);

for i=1:nface
    for k=1:3
        kk = mod(k,3)+1;
        if nargin<2
            A(face(i,k),face(i,kk)) = 1;
        else
            v = vertex(:,face(i,k))-vertex(:,face(i,kk));
            A(face(i,k),face(i,kk)) = sqrt( sum(v.^2) );    % euclidean distance
        end
    end
end
% make sure that all edges are symmetric
A = max(A,A');

end

function A = compute_edge_face_ring(face)

% compute_edge_face_ring - compute faces adjacent to each edge
%
%   e2f = compute_edge_face_ring(face);
%
%   e2f(i,j) and e2f(j,i) are the number of the two faces adjacent to
%   edge (i,j).
%
%   Copyright (c) 2007 Gabriel Peyre


[tmp,face] = check_face_vertex([],face);

n = max(face(:));
m = size(face,2);
i = [face(1,:) face(2,:) face(3,:)];
j = [face(2,:) face(3,:) face(1,:)];
s = [1:m 1:m 1:m];

% first without duplicate
[tmp,I] = unique( i+(max(i)+1)*j );
% remaining items
J = setdiff(1:length(s), I);

% flip the duplicates
i1 = [i(I) j(J)];
j1 = [j(I) i(J)];
s = [s(I) s(J)];

% remove doublons
[tmp,I] = unique( i1+(max(i1)+1)*j1 );
i1 = i1(I); j1 = j1(I); s = s(I);

A = sparse(i1,j1,s,n,n);


% add missing points
I = find( A'~=0 );
I = I( A(I)==0 );
A( I ) = -1;

end


function fring = compute_face_ring(face)

% compute_face_ring - compute the 1 ring of each face in a triangulation.
%
%   fring = compute_face_ring(face);
%
%   fring{i} is the set of faces that are adjacent
%   to face i.
%
%   Copyright (c) 2004 Gabriel Peyre

% the code assumes that faces is of size (3,nface)
[tmp,face] = check_face_vertex([],face);

nface = size(face,2);
nvert = max(max(face));

A = compute_edge_face_ring(face);
[i,j,s1] = find(A);     % direct link
[i,j,s2] = find(A');    % reverse link

I = find(i<j);
s1 = s1(I); s2 = s2(I);

fring{nface} = [];
for k=1:length(s1)
    if s1(k)>0 && s2(k)>0
        fring{s1(k)}(end+1) = s2(k);
        fring{s2(k)}(end+1) = s1(k);
    end
end

return;

for i=1:nface
    face(i,:) = sort(face(i,:));
end

edges = zeros(3*nface,2);

edges( 1:nface, : ) = face(:,1:2);
edges( (nface+1):2*nface, : ) = face(:,2:3);
edges( (2*nface+1):3*nface, : ) = face(:,[1,3]);

for i=1:nface
    fring{i} = [];
end

for i=1:3*nface
    e = edges(i,:);
    I = find(edges(:,1)==e(1) & edges(:,2)==e(2));
    if length(I)==2
        f1 = mod( I(1)-1, nface)+1;
        f2 = mod( I(2)-1, nface)+1;
        fring{f1} = [fring{f1}, f2];
        fring{f2} = [fring{f2}, f1];
        edges(I,:) = rand(length(I),2);
    end
end

return;

return;
h = waitbar(0,'Computing Face 1-ring');
for i=1:nface
    waitbar(i/nface);
    fring{i} = [];
    for j=1:3
        j1 = j;
        j2 = mod(j,3)+1;
        v1= face(i,j1);
        v2= face(i,j2);
        % test if another face share the same vertices
        f = [];
        for i1=1:nface
            for a=1:3
                if face(i1,a)==v1 && i1~=i
                    for b=1:3
                        if face(i1,b)==v2
                            % add to the ring
                            fring{i} = [fring{i}, i1];
                        end
                    end
                end
            end
        end

    end
end
close(h);

end

function [e,v,g] = compute_butterfly_neighbors(k, nj)

%   compute_butterfly_neighbors - compute local neighbors of a vertex
%
%   [e,v,g] = compute_butterfly_neighbors(k, nj);
%
%   This is for internal use.
%
% e are the 2 direct edge neighbors
% v are the 2 indirect neighbors
% g are the fare neighbors
%
%   You need to provide:
%       for e: vring, e2f
%       for v: fring
%       for g: facej
%
%   Copyright (c) 2007 Gabriel Peyre

global vring e2f fring facej;

% find the 2 edges in the fine subdivition
vr = vring{k};
I = find(vr<=nj);
e = vr(I);
% find the coarse faces associated to the edge e
f = [e2f(e(1),e(2)) e2f(e(2),e(1))];
% symmetrize for boundary faces ...
f(f==-1) = f(3 - find(f==-1));

if nargout>1
    F1 = mysetdiff(fring{f(1)}, f(2));
    F2 = mysetdiff(fring{f(2)}, f(1));
    % symmetrize for boundary faces
    F1 = [F1 repmat(f(1), 1, 3-length(F1))];
    F2 = [F2 repmat(f(2), 1, 3-length(F2))];
    v = [ mysetdiff( facej(:,f(1)), e )   mysetdiff( facej(:,f(2)), e ) ];
    if nargout>2
        d = [v, e];
        g = [setdiff( facej(:,F1(1)),d ), ...
            mysetdiff( facej(:,F1(2)),d ), ...
            mysetdiff( facej(:,F2(1)),d ), ...
            mysetdiff( facej(:,F2(2)),d ) ];
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = mysetdiff(a,b)
% removed in a entries equal to entries in b
for s=b
    a(a==s) = [];
end
end
