  function [V,F] = read_mesh( filename )
% function [V,F] = read_mesh( filename )
%
% Loads a mesh in Wavefront OBJ format.
%
% INPUT:
%    filename - path to mesh file
%
% OUTPUT:
%    V - 3x|V| matrix of vertex positions
%    F - |F|x3 matrix of faces as 1-based indices into edge list
%

   fp = fopen( filename, 'r' );
   if( fp == -1 )
      disp( sprintf( 'Error: could not read mesh file "%s"\n', filename ));
      return;
   end

   V = zeros( 3, 0 );
   F = zeros( 0, 3 );
   while( ~feof( fp ))
       line = fscanf( fp, '%c %f %f %f' );
       label = line( 1 );
       values = line( 2:4 );

       if( label == 'v' ) % vertex
          V = [ V, values(1:3) ];
       end

       if( label == 'f' ) % face
          F = [ F; values' ];
       end
   end
   fclose( fp );

end

