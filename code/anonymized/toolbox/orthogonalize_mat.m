function X = orthogonalize_mat(X)

% orthogonalize_mat - spectral orthogonalization 
%
%   X = orthogonalize_mat(X);
%
%   Copyright (c) 2017 [anonymized]

[U,S,V] = svd(X);
X = U*V';

end 