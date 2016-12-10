function [a,e] = eigen_remaper(u,v,direc)

% eigen_remaper - compute (ansiotropy,energy) from eigenvalues.
%
%   [a,e] = eigen_remaper(u,v,+1)
%   [u,v] = eigen_remaper(u,v,-1)
%
% If direc==1, compute (ansiotropy,energy) from a pair u>v as
%    a = (u-v)./(u+v);
%    e = sqrt(u+v);
% If direct==-1, inverse the computation.
%
%   a is clost to 1 for high anisotropy, and close to 0 for low anisotropy.
%
%   Copyright (c) Gabriel Peyre

if direc==+1
    a = (u-v)./(u+v);
    e = sqrt(u+v);
else
    a = u; e = v;
    u = (1+a).*e.^2/2;
    v = (1-a).*e.^2/2;
    a = u; e = v; 
end


end