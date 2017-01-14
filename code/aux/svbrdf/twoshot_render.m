% This is ported back from the C++ renderer (the one used by the optimizer
% assumes eye position and light position coincide, so it's slightly
% simpler)

function img = twoshot_render(mat, viewing)
    S = size(mat.alb_diff);
    S = S(1:2);

    % Flatten the images into horizontal lists of pixels, for convenient
    % bulk computations
    N = reshape(permute(mat.normal, [3 1 2]), 3, []);
    L = reshape(permute(viewing.light, [3 1 2]), 3, []);
    E = reshape(permute(viewing.eye, [3 1 2]), 3, []);
    
    % Halfway
    H = L+E;
    H = H ./ repmat(sqrt(sum(H.^2,1)), [3 1]);
    H = reshape(permute(H, [3 1 2]), 3, []);

    % Rotation to local BRDF coordinates
    Nmat = rot_from(N);
    
    Hn = mat_mul_bulk3(Nmat, H);
    
    Hnn = Hn ./ repmat(Hn(3,:), [3 1]);
    
    M = reshape(permute(mat.gloss, [3 4 1 2]), 2, 2, []);
    
    HnnW = mat_mul_bulk2(M, Hnn(1:2,:));
    HnnWD = sqrt(sum(HnnW.*Hnn(1:2,:)));

    Fs = exp(-HnnWD.^mat.alpha);

    F0 = 0.04;
	fres = F0 + (1-F0)*(1.0-max(0,(sum(H.*L,1)))).^5;
    Fs = (Fs .* fres) / F0; % Fresnel
					
	Fs = Fs ./ sum(H.*L,1); % from Brady et al.
    
    Scol = reshape(permute(mat.alb_spec, [3 1 2]), 3, []);
    Dcol = reshape(permute(mat.alb_diff, [3 1 2]), 3, []);
    
    dshade = max(0,sum(N .* L, 1));
    
    img = repmat(dshade .* viewing.invsq(:)', [3 1]) .* (Scol .* repmat(Fs, [3 1]) + Dcol);
    img = viewing.illum .* permute(reshape(img, [3 S]), [2 3 1]);
end

% helper
function R = rot_from(a)
    n = size(a,2);
    R = zeros(3,3,n);
    b = a ./ repmat(sqrt(a(3,:)+1), [3 1]);
    
    R(1,1,:) = -b(1,:).^2 + 1;
    R(1,2,:) = -b(1,:).*b(2,:);
    R(2,1,:) = R(1,2,:);
    R(2,2,:) = -b(2,:).^2 + 1;
    R(3,3,:) = -b(1,:).^2 - b(2,:).^2 + 1;
    R(1,3,:) = -a(1,:);
    R(3,1,:) = a(1,:);
    R(2,3,:) = -a(2,:);
    R(3,2,:) = a(2,:);
end

