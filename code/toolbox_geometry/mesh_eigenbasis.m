function U = mesh_eigenbasis(V,F, options)

% mesh_eigenbasis - compute eigenbasis for tangent/normal space
%
%   U = mesh_eigenbasis(V,F, basis_mode);
%
%   U(:,1,i) and U(:,2,i) are the tangents
%   U(:,3,i) is the normal
%
%   Use the code of Keenan Crane.
%
%   S = options.singularity: list of singularities as (vertex ID, index) pairs
%       NOTE: indices *must* sum to 2!
%   S = [ 1, 2 ];                    % one singularity of index +2
%   S = [ 1, 1; 300, 1 ];            % two singularities of index +1
%   S = [ 1, 1; 300, -1; 500, 2 ];   % three singularities: +1, -1, +2
%   S = [ 1, .5; 300, -.5; 500, 2 ]; % three singularities: +1/2, -1/2, +2

options.null = 0;
basis_mode = getoptions(options, 'basis_mode', 'connection');
S = getoptions(options, 'singularity', [ 1, 2 ]); 

N = size(V,2);
crossp = @(a,b)[a(:,2).*b(:,3)-a(:,3).*b(:,2), -a(:,1).*b(:,3)+a(:,3).*b(:,1), a(:,1).*b(:,2)-a(:,2).*b(:,1)];
normalize3 = @(x)x ./ repmat(sqrt(sum(x.^2,2)), [1 3]);
resh = @(x)reshape(x, [3 1 N]);
U0 = compute_normal(V,F); % normal

switch basis_mode
    case 'random'
        U1 = normalize3( crossp( U0', randn(N,3) ) )';
    case 'connection'
        %%% using Keenan Crane's method and code %%%              
        % Building exterior derivatives...
        [d0,d1] = build_exterior_derivatives( V, F' );
        % Getting Riemannian holonomy...
        K = get_gaussian_curvature( V, F' );
        % Computing trivial connection...
        warning off;
        x = compute_trivial_connection( d0, d1, K, S );
        warning on;
        % Extracting parallel direction field...
        Z = extract_direction_field( V, F', d0, d1, x );
        % display
        % plot_direction_field( V, F', U1, S );
        % interpolate on vertices
        U1 = zeros(3,N);
        for i=1:size(F,2)
            f = F(:,i);
            for j=1:3
                U1(:,f(j)) = U1(:,f(j)) + Z(:,i);
            end
        end
        % orthogonalize with respect to normal
        U1 = normalize3( crossp( U1', U0' ) )';
    otherwise
        error('Unknown');
end

U2 = normalize3( crossp( U1', U0' ) )';
U = cat(2, resh(U1), resh(U2), resh(U0));


end