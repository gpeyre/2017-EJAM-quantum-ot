  function [d0,d1] = build_exterior_derivatives( V, F )
% function [d0,d1] = build_exterior_derivatives( V, F )
%
% Builds the discrete exterior derivative operators on 0- and 1-forms.
%
% INPUT:
%    V - 3x|V| matrix of vertex positions
%    F - |F|x3 matrix of faces as 1-based indices into edge list
%    e2v - |e2v|x2 matrix of edges as 1-based indices into vertex list
%    f2e - |f2e|x3 matrix of faces as 1-based indices into edge list
%
% OUTPUT:
%    d0 - sparse |e2v|x|V| matrix representing discrete exterior derivative on 0-simplices
%    d1 - sparse |f2e|x|e2v| matrix representing discrete exterior derivative on 1-simplices
%

   % construct a temporary map from edge indices to vertex indices (e2v),
   % and from face indices to edge indices (f2e).  The latter map gives
   % a positive or negative sign to edge indices depending on their
   % orientation relative to the orientation of the corresponding
   % face -- positive for consistent orientation, negative otherwise
   e2v = zeros( 0, 2 );
   f2e = zeros( size(F,1), 3 );
   edgeIndex = sparse( size(V,2), size(V,2) );
   nEdges = 0;

   for i = 1:size( F, 1 )
      for j = 0:2
         v1 = F( i, 1+j );
         v2 = F( i, 1+mod(j+1,3) );
         if( v2 < v1 )
            tmp = v1;
            v1 = v2;
            v2 = tmp;
            sgn = -1;
         else
            sgn = 1;
         end

         if( edgeIndex( v1, v2 ) == 0 )
            nEdges = nEdges+1;
            edgeIndex( v1, v2 ) = nEdges;
            e2v = [ e2v; [ v1, v2 ]];
         end

         f2e( i, 1+j ) = sgn*edgeIndex( v1, v2 );
      end
   end

   % now that we have vertex-edge and edge-face adjacency
   % information, build the exterior derivative matrices
   d0 = sparse( size(e2v,1), size(V,2) );
   d1 = sparse( size(f2e,1), size(e2v,1) );
   
   % d0
   for i = 1:size(e2v,1)
      d0( i, e2v( i, 1 )) = -1;
      d0( i, e2v( i, 2 )) =  1;
   end
   
   % d1
   for i = 1:size(f2e,1)
      for j = 1:3
         d1( i, abs( f2e( i, j ))) = sign( f2e( i, j ));
      end
   end
end

