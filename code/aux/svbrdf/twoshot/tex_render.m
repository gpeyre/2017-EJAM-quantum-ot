% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function img = tex_render(cX, cY, X, g)
    g(4) = exp(g(4)) + 0.5;
    
    % always: horizontal is points, vertical is variables    
    np = size(cX,2);
    nv = size(X,1);
    
    if size(X,2) == 1
        X = repmat(X, [1 np]);
    end
    
    D2 = cY.^2 + cX.^2;
    inv_sq = 1 ./ (1 + D2);
    inv_sq = inv_sq .* exp(-D2/(exp(g(3))+0.3)); % bake the flash emission here for convenience
   
    % Eye/light/half-vectors
    E = [-cX; -cY; ones(1,np)];
    
    E = E ./ repmat(sqrt(sum(E.^2,1)), [3 1]);
    
    % Normal vectors
    N = [X(1:2,:); ones(1,np)];
    N = N ./ repmat(sqrt(sum(N.^2,1)), [3 1]);
    Nmat = rot_from(N);
    
    En = mat_mul_bulk3(Nmat,E);
    
    % En is our eye vector in coordinates where the normal is [0 0 1]'.
    
    diffuse = repmat(max(0,En(3,:)), [3 1]) .* exp(X(3:5,:));

    
    % Enn(1:2) are the coordinates on the unit distance plane perpendicular
    % to the normal
    Enn = En ./ repmat(En(3,:), [3 1]);
    
    M = expm_bulk(X(7:9,:));
    M(1,1,:) = M(1,1,:) + 0.001;
    M(2,2,:) = M(2,2,:) + 0.001;
    M = mat_inv_bulk2(M);

    EnnW = mat_mul_bulk2(M, Enn(1:2,:));
    EnnWD = sqrt(sum(EnnW.*Enn(1:2,:)));
%    detM = vec(M(1,1,:).*M(2,2,:) - M(1,2,:).^2)';

    Fs = exp(-EnnWD.^g(4));

    YI = [1 0 1.13983; 1 -0.39465 -0.5806; 1 2.03211 0];
    Scol =  max(0, YI* ([1;g(1);g(2)]) *exp(X(6,:)));
    
    img = repmat(max(0,En(3,:)), [3 1]) .* Scol .* repmat(Fs, [3 1]) + diffuse;

    img = img .* repmat(inv_sq .* (En(3,:)>0), [3 1]);
end



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

