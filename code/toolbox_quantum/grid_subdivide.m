function [i,j,c] = grid_subdivide(i,j,n,d, nsub)

% grid_subdivide - subdivide a sparse grid
%
%   [i,j,c] = grid_subdivide(i,j,n, d);
%
%   d is the dimension
%   (i(k),j(k))_k are set of coupling location in dimension d, 
%       the output is a set of 2^(2*d) coupling points.
%
%   Copyright (c) 2016 Gabriel Peyre


if d~=2
    error('Only implemented in 2D');
end

if nsub>1
    for k=1:nsub
        [i,j,c] = grid_subdivide(i,j,n,d, 1);
    end
    return;
end


% 2D indexing
[i1,i2] = ind2sub([n n], i);
[j1,j2] = ind2sub([n n], j);
I1 = []; I2 = []; J1 = []; J2 = [];
for a1=0:1
    for a2=0:1
        for b1=0:1
            for b2=0:1
                I1 = [I1;2*i1-1+a1];
                I2 = [I2;2*i2-1+a2];
                J1 = [J1;2*j1-1+b1];
                J2 = [J2;2*j2-1+b2];                
            end
        end
    end
end
n1 = 2*n;
i = sub2ind([n1 n1], I1, I2);
j = sub2ind([n1 n1], J1, J2);
cIJ = (I1-J1).^2/n1^2 + (I2-J2).^2/n1^2;
c = sparse(i,j,cIJ,n1*n1,n1*n1);

end