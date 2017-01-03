function Q = tensor_inv(P)

% tensor_inv - explicit inverse of a tensor field
%
%   Q = tensor_inv(P)
%
%   Q(:,:,i,j) = inv( P(:,:,i,j) )
%   
%   Copyright (c) 2016 Gabriel Peyre

d = size(P,1); 
switch d
    case 1
        Q = 1./P;
    case 2 
        % det
        m = P(1,1,:,:).*P(2,2,:,:) - P(2,1,:,:).*P(1,2,:,:);
        Q = P; 
        Q(1,1,:,:) = P(2,2,:,:)./m;
        Q(2,2,:,:) = P(1,1,:,:)./m;
        Q(1,2,:,:) = -P(1,2,:,:)./m;
        Q(2,1,:,:) = -P(2,1,:,:)./m;
    case 3
        error('Need to be implemented using tensor_eigenrecomp');
end

end
