% Builds the coordinate system and a simple near-field point light 
% condition for twoshot_render()

% S size of image
% P position of point light (e.g. [0.0 0.0 1.0] for above "origin")
% fov field of view, 1 is fine
function viewing = twoshot_pointlight(S, P, fov)
    a = S(1)/S(2);  % aspect
    
    % Eye vector field
    E = cat(3, fov*repmat(-linspace(-1, 1, S(2)), [S(1) 1]), ...
               fov*repmat(-linspace(a, -a, S(1))', [1 S(2)]), ...
               ones(S));
    E = E ./ repmat(sqrt(sum(E.^2,3)), [1 1 3]);    % normalize
    
    % Light direction vector field
    L = cat(3, P(1)+repmat(-linspace(-1, 1, S(2)), [S(1) 1]), ...
               P(2)+repmat(-linspace(a, -a, S(1))', [1 S(2)]), ...
               P(3)*ones(S));
    
    I = 1./sum(L.^2, 3);    % inverse square distance
    L = L .* repmat(sqrt(I), [1 1 3]);  % normalize

    viewing = struct;
    viewing.eye = E;
    viewing.light = L;
    viewing.invsq = I;
    viewing.illum = ones([S 3]);
end

