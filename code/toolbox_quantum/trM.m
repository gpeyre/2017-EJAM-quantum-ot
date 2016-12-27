function y = trM(x, do_squeeze)

% trM - trace of a tensor
%
%   y = trM(x, do_squeeze=0);
%
%   y(1,1,:,:) = sum_i x(i,i,:,:)
%
%   squeeze() is applied to y if do_squeeze=1.
%
%   Copyright (c) Gabriel Peyre

if nargin<2
    do_squeeze = 0;
end

y = zeros(1,1,size(x,3),size(x,4));
for i=1:size(x,1)
    y = y+x(i,i,:,:);
end

if do_squeeze==1
    y = squeeze(y);
end

end