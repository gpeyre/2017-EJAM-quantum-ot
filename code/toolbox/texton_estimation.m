function [Ts,f0,T,S,U] = texton_estimation(f, n, options)

% texton_estimation - estimate covariance and texton
%
%   [Ts,f0,T,S,U]  = texton_estimation(f, n, options); 
%
%   Ts(:,:,i,j) is the sqtm of T(:,:,i,j), i.e. the texton.
%
%   calling synth_func(0) generate each time a new random texture.
%   calling synth_func(k) with k>0 generate a fixed new texture. 
%
%   n is the size of the synthesized texture. 
%
%   Copyright (c) 2016 Gabriel Peyre

% number of extracted window
estimator_type = getoptions(options, 'estimator_type', 'periodic');
q = getoptions(options, 'samples', 200);

d = size(f,3);
n0 = size(f);

% window for estimation
w = sin( pi*(0:n-1)'/n ).^2; w = w*w';

T = zeros(d,d,n,n);
for it=1:q
    progressbar(it,q);
    a = 1+ floor( [rand*(n0(1)-n) rand*(n0(1)-n)] );
    F = f(a(1):a(1)+n-1, a(2):a(2)+n-1,:); 
    switch estimator_type
        case 'window'
            h = repmat(w, [1 1 d]) .* F;
        case 'periodic'
            h = periodic_comp(F);
    end
    h = h - repmat( mean(mean(F,1),2), [n n 1] );
    H = fft2(h);
    for i=1:d
        for j=1:d
            T(i,j,:,:) = T(i,j,:,:) + reshape(H(:,:,i) .* conj(H(:,:,j)), [1 1 n n]);
        end
    end
end
T = T/q;

% extract eigenvalues
S = zeros(n,n,d); U = [];
for i=1:n
    for j=1:n
        [U(:,:,i,j),s] = eig(T(:,:,i,j)); s = diag(s);
        S(i,j,:) = reshape( s, [1 1 d] ); 
    end
end

% compute sqrtm of the covariance
warning off;
for i=1:n
    for j=1:n
        Ts(:,:,i,j) = sqrtm(T(:,:,i,j));
    end
end
warning on;

f0 = mean(mean(f,1),2);


end
