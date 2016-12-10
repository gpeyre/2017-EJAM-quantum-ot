function X = extract_direction_field( V, F, d0, d1, x )
% function X = extract_direction_field( V, F, d0, d1, x )
%
% INPUT:
%    V  - 3x|V| matrix of vertex positions
%    F  - |F|x3 matrix of faces as 1-based indices into vertex list
%    d0 - sparse |E|x|V| matrix representing discrete exterior derivative on 0-simplices
%    d1 - sparse |F|x|E| matrix representing discrete exterior derivative on 1-simplices
%
% OUTPUT:
%    X - 3x|F| matrix specifying a unit vector in each face
%

% allocate space for the vector field (one 3x1 vector per face)
X = zeros( 3, size(F,1) );

% initialize a record of faces that have been visited (1=visited; 0=not visited)
visited = zeros( size(F,1), 1 );

% choose an arbitrary initial face and unit vector
root = 1;
w = V(:,F(root,1)) - V(:,F(root,2));
%% w = randn(3,1);  % GP
w = w / norm(w);

% enqueue the root face
X(:,root) = w;
Q = [ root ];
visited( root ) = 1;

% transport the initial vector to all other triangles
while( length( Q ) > 0 )
    % dequeue face i
    i = Q( end );
    Q = Q( 1:end-1 );
    
    % visit the neighbors of face i
    for j = triangle_neighbors( i, d1 )
        
        % ignore neighbors we've already visited
        if( ~visited( j ))
            
            % transport the tangent vector at face i to face j
            k = sharedEdge( i, j, d1 );
            X(:,j) = transport( X(:,i), i, j, k, V, F, x );
            
            % enqueue face j
            Q = [ j, Q ];
            visited( j ) = 1;
        end
    end
end
end

function k = sharedEdge( i, j, d1 )
% function k = sharedEdge( i, j, d1 )
%
% Returns the index of the edge shared by faces i and j (or the empty
% set if i and j do not share a face).
%
% INPUT:
%    i - index of first face
%    j - index of second face
%    d1 - sparse |F|x|E| matrix representing discrete exterior derivative on 1-simplices

k = intersect( find(d1(i,:)), find(d1(j,:)) );

end

function I = triangle_neighbors( i, d1 )
% function I = triangle_neighbors( i, d1 )
%
% Finds the indices of all triangles that share a single edge with triangle i.
%
% INPUT:
%    i  - triangle index
%    d1 - sparse |F|x|E| matrix representing discrete exterior derivative on 1-simplices
%
% OUTPUT:
%    N - 1xn vector containing indices of n neighboring triangles
%

E = find( d1( i, : ));
[R,I] = find( d1( :, E )');
I = unique( I );
I = setdiff( I, [i] )';

end

function w = transport( w0, i, j, k, V, F, x )
% function w = transport( w0, i, j, k, V, F, x )
%
% Parallel transport the unit tangent vector w0 from face i to
% face j using the connection specified by x.
%
% INPUT:
%    w0 - 3x1 vector of coordinates for tangent vector
%    i  - index of source triangle
%    j  - index of destination triangle
%    k  - index of shared edge
%    V  - 3x|V| matrix of vertex positions
%    F  - |F|x3 matrix of faces as 1-based indices into vertex list
%    x  - |E|x1 vector giving deviation from Levi-Civita on each dual edge
%
% OUTPUT:
%    w - 3x1 vector of coordinates for transported vector
%

% Grab vertices a, b, c, and d of the two adjacent
% triangles i and j according to the following labels:
%
%                         b
%                        /|\
%                       / | \
%                      /  |  \
%                     /   |   \
%                    c  i | j  d
%                     \   |   /
%                      \  |  /
%                       \ | /
%                        \|/
%                         a

I = sort( intersect( F(i,:), F(j,:) ));
a = V( :, I(1) );
b = V( :, I(2) );
c = V( :, setdiff( F(i,:), I ) );
d = V( :, setdiff( F(j,:), I ) );

% compute the change of basis between triangles i and j
Ei = orthogonalize( b-a, c-a );
Ej = orthogonalize( b-a, b-d );

% also compute the rotation in the plane by the connection angle x(k)
R = [ cos(x(k)), -sin(x(k)), 0; ...
    sin(x(k)),  cos(x(k)), 0; ...
    0,             0, 1 ];

% compose these maps to compute the parallel transport of the
% vector w0 from triangle i to triangle j via the trivial connection
w = Ej*R*inv(Ei)*w0;

end

function E = orthogonalize( u, v )
% returns a 3x3 orthonormal matrix whose first column is parallel
% to u and whose final column is orthogonal to u and v
e1 = u;
e1 = e1/norm(e1);

e2 = v;
e2 = e2-(e2'*e1)*e1;
e2 = e2/norm(e2);

e3 = cross( e1, e2 );
e3 = e3/norm(e3);

E = [ e1, e2, e3 ];
end

