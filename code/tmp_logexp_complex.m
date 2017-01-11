%%
% test for expM/logM of complex matrices

global logexp_fast_mode;
logexp_fast_mode = 2;

d = 2;
A = randn(d)+1i*randn(d);
X = A*A';

logm(X)
logM(X)

expm(logm(X))
expM(logm(X))