  function plot_direction_field( V, F, X, S )
% function plot_direction_field( V, F, X, S )
%
% Plots a tangent vector field on a mesh.  Also plots singular vertices.
%
% INPUT:
%    V - 3x|V| matrix of vertex positions
%    F - |F|x3 matrix of faces as 1-based indices into edge list
%    X - 3x|F| matrix specifying a unit vector in each face
%    S  - Nx2 vector specifying N singularities as (vertex ID, index) pairs
%

   figure(1);
   hold on;

   % plot mesh
   colormap gray;
   trimesh( F, V(1,:), V(2,:), V(3,:), zeros(length(V),1) );

   % plot vector field
   eps = 1e-3;
   for i = 1:size(F,1)
      c = barycenter( i, F, V );
      r = inradius( i, F, V );
      N = normal( i, F, V );
      a = c + r*X(:,i) + eps*N;
      b = c - r*X(:,i) + eps*N;
      plot3( [a(1),b(1)], ...
             [a(2),b(2)], ...
             [a(3),b(3)], ...
             'Color', [0 0 1], ...
             'LineWidth', 2 );
   end

   % plot singularities
   for i = 1:size(S,1)
      j = S(i,1);
      H = plot3( V(1,j), V(2,j), V(3,j), 'r.' );
      set( H, 'MarkerSize', 40 );
   end

   axis equal;
   axis off;
   view([-180,0]);
   zoom(1);
   rotate3d;
   hold off;

end

function c = barycenter( i, F, V )
% computes the barycenter of face i
   c = (sum(V(:,F(i,:))')/3)';
end

function r = inradius( i, F, V )
% computes the inradius of face i
   a = V( :, F(i,1) );
   b = V( :, F(i,2) );
   c = V( :, F(i,3) );

   u = norm( a-b );
   v = norm( b-c );
   w = norm( c-a );

   r = (1/2)*sqrt( ((u+v-w)*(w+u-v)*(v+w-u)) / (u+v+w) );
end

function N = normal( i, F, V )
% computes the normal of face i
   a = V( :, F(i,1) );
   b = V( :, F(i,2) );
   c = V( :, F(i,3) );

   N = cross( b-a, c-a );
   N = N / norm( N );
end

