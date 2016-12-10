function y = trM(x)

% trM - trace of a tensor
%
%   y = trM(x);
%
%   y(1,1,:,:) = sum_i x(i,i,:,:)
%
%   Copyright (c) Gabriel Peyre

y = zeros(1,1,size(x,3),size(x,4));
for i=1:size(x,1)
    y = y+x(i,i,:,:);
end

end