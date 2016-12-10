%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
% run.m
% Keenan Crane
% June 15, 2010
%
% This program computes smooth direction fields on meshes using the method
% described in Crane et al, "Trivial Connections on Discrete Surfaces."
% The primary purpose of this code is to illustrate how the algorithm works;
% it is likely not suitable for large meshes.  Furthermore, for simplicity,
% this version of the code works only on genus-0 surfaces without boundary,
% i.e., topological spheres.
%
% This script is the main routine, and should be called from the MATLAB
% command line by simply typing 'run'.  A mesh file can be specified by
% editing the lines under INPUT below.  Meshes should be in Wavefront OBJ
% format with ONLY lines of the form
%
%    v x y z
%
% or
%
%    f i1 i2 i3
%
% i.e., vertex coordinates or triangles (no comments, texture coordinates,
% etc.).  An example file can be found in 'test.obj'.
%
% Singularities are specified via (vertex ID, index) pairs (also under INPUT)
% where the vertex IDs should correspond to 1-based indices used in the
% corresponding OBJ file.  Indices are used to specify the behavior of
% singularities at these vertices, and must sum to 2, which is the Euler
% characteristic of a topological sphere.
%
% ----------
%
% This code was developed for an introductory undergraduate course on discrete
% differential geometry, and follows three prior assignments:
%
%   Topological Invariants of Discrete Surfaces
%   http://www.cs.caltech.edu/~keenan/ddg_exercises/ddg_hw1.pdf
%
%   Mesh Smoothing
%   http://www.cs.caltech.edu/~keenan/ddg_exercises/ddg_hw2.pdf
%
%   Vector Field Decomposition
%   http://www.cs.caltech.edu/~keenan/ddg_exercises/ddg_hw3.pdf
%
% The material above builds all of the necessary concepts (both code and math)
% from the ground up, and should provide a concrete understanding of most of
% the ideas used in this code -- the remaining ideas are discussed in the paper
%
%    Crane, Desbrun, Schröder, "Trivial Connections on Discrete Surfaces"
%    SGP 2010 / Computer Graphics Forum
%    http://www.cs.caltech.edu/~keenan/pdf/connections.pdf
%
% ----------
%
% Finally, please do not hesitate to send bugs or
% questions to keenan@cs.caltech.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&

%==== INPUT ====================================================================

input_filename = 'test.obj';

% list of singularities as (vertex ID, index) pairs
% NOTE: indices *must* sum to 2!
% -------------------------------------------------
S = [ 1, 2 ];                    % one singularity of index +2
% S = [ 1, 1; 300, 1 ];            % two singularities of index +1
% S = [ 1, 1; 300, -1; 500, 2 ];   % three singularities: +1, -1, +2
% S = [ 1, .5; 300, -.5; 500, 2 ]; % three singularities: +1/2, -1/2, +2


%==== MAIN PROGRAM =============================================================

stdout = 1;

fprintf( stdout, 'Reading mesh...\n' );
[V,F] = read_mesh( input_filename );

fprintf( stdout, 'Building exterior derivatives...\n' );
[d0,d1] = build_exterior_derivatives( V, F );

fprintf( stdout, 'Getting Riemannian holonomy...\n' );
K = get_gaussian_curvature( V, F );

fprintf( stdout, 'Computing trivial connection...\n' );
x = compute_trivial_connection( d0, d1, K, S );

fprintf( stdout, 'Extracting parallel direction field...\n' );
X = extract_direction_field( V, F, d0, d1, x );

fprintf( stdout, 'Displaying result...\n' );
plot_direction_field( V, F, X, S );

fprintf( stdout, 'Done.\n' );

