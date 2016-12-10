  function x = compute_trivial_connection( d0, d1, K, S )
% function x = compute_trivial_connection( d0, d1, K, S )
%
% Computes a set of connection angles x associated with each dual edge of a mesh that
% describe a trivial connection with the specified singularities.  These angles are
% expressed relative to the discrete Levi-Civita connection, i.e., ``translation''
% across the shared edge between two faces.
%
% INPUT:
%    d0 - sparse |E|x|V| matrix representing discrete exterior derivative on 0-simplices
%    d1 - sparse |F|x|E| matrix representing discrete exterior derivative on 1-simplices
%    K  - |V|x1 vector of per-vertex Gaussian curvatures
%    S  - Nx2 vector specifying N singularities as (vertex ID, index) pairs -- indices must sum to 2
%
% OUTPUT:
%    x - |E|x1 vector of angles encoding deviations from the Levi-Civita connection
%

   % for surfaces of genus 0 (i.e., topological spheres), our constraint matrix
   % is simply the transpose of d0, the discrete exterior derivative on 0-forms
   A = d0';

   % in the genus-0 case, the right hand side of our system is the Riemannian
   % holonomy around vertices (i.e., Gaussian curvature) minus the target
   % curvature at singular vertices
   b = K;
   for s = S'
      b(s(1)) = K(s(1)) - 2*pi*s(2);
   end
   
   % first, compute *any* solution x to the underconstrained system Ax=-b
   x = A \ -b;

   % next, remove the null-space component of x to get the solution
   % to x with smallest l2 norm -- this is done by applying the
   % projection x - d1'*inv(d1*d1')*d1*x, where d1 is the discrete
   % exterior derivative on 1-forms
   L = d1*d1';
   y = L \ (d1*x);
   x = x - d1'*y;

end

