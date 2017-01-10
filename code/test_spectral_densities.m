%%
% Test for matrix-valued densities.
% Need to implement a windowing estimator. 

addpath('toolbox/');
addpath('toolbox_quantum/');

d = 3; 
N = 128; 

X = randn(N,3);

% outer product
Xf = fft(X);
% Xf = Xf(1,:);
Cf = tensor_mult( reshape(permute(Xf, [2 1]), [3 1 N]), reshape(Xf', [1 3 N]) );
% Cf = abs( fft(C, [], 3) ).^2;

U = [];S = [];
for i=1:N
    [U(:,:,i), s] = eig(Cf(:,:,i)); 
    s = abs( real( diag(s) ) );
    [s,I] = sort(s); 
    S(:,i) = s; 
    U(:,:,i) = U(:,I,i);
end

