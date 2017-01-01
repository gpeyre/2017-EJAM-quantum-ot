function [P,err] = quantum_interp(Q,R, t, options)

% quantum_interp - interpolation (without transport) between two tensors
%
%   [P,err] = quantum_interp(Q,R, t, options);
%
% ***** For options.barymode=='hard' *****
%
%   P = expM( (1-t)*logM(Q) + t*logM(R) )
%
% ***** For options.barymode=='soft' *****
%
% P solves:
%   min_P (1-t)*D(P,Q)^2 + t*D(P,R)^2
% where
%   D(P,Q)^2 = tr( P+Q - 2*m(P,Q) )
% where the log-geometric-mean is
%   m(P,Q) = expM( logM(P)/2+logM(Q)/2 )
%
%   Copyright (c) Gabriel Peyre 2016

options.null = 0;
niter = getoptions(options, 'niter', 100);
barymode = getoptions(options, 'barymode', 'hard');
d = size(Q,1);

if length(t)>1
    P = [];
    for k=1:length(t)
        [P(:,:,k),err] = quantum_interp(Q,R, t(k), options);
    end
    return;
end

err = [];
lQ = logM(Q); lR = logM(R);

if strcmp(barymode, 'hard')
    P =  expM( (1-t)*lQ + t*lR );
    return;
end

A = zeros(d);
for it=1:niter
    A1 = 2*logM( (1-t)*expM(A/2+lQ/2) + t*expM(A/2+lR/2) ) - A;
    err(it) = norm(A-A1, 1);
    A = A1;
end
P = expM(A);

end