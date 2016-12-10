  function K = get_gaussian_curvature( V, F )
% function K = get_gaussian_curvature( V, F )
%
% Computes the standard discretization of Gaussian curvature at
% mesh vertices, i.e., 2pi minus the sum of the angles around each
% vertex.
%
% INPUT:
%    V - 3x|V| matrix of vertex positions
%    F - |F|x3 matrix of faces as 1-based indices into vertex list
%
% OUTPUT:
%    K - |V|x1 vector of per-vertex Gaussian curvatures
%

   % grab the number of vertices and faces
   nV = size( V, 2 );
   nF = size( F, 1 );

   % initialize the vector of curvatures (initially 2pi at every vertex)
   K = zeros( nV, 1 ) + 2*pi;

   % iterate over faces, subtracting the contribution of each angle
   % from the Gaussian curvature of the corresponding vertex
   for f = F'
      for i = 1:3
         % get unit vectors u,v along two of the edges
         u = V(:,f(2))-V(:,f(1)); u = u/norm(u);
         v = V(:,f(3))-V(:,f(1)); v = v/norm(v);

         % subtract the angle between u and v from the total
         theta = acos(u'*v);
         K(f(1)) = K(f(1)) - theta;

         % continue to the next pair of edges
         f = circshift(f,1);
      end
   end
end

