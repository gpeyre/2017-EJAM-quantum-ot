function D = tensor_id(a, d)

% tensor_id - diagonal tensor creation
%
%   D = tensor_id(a, d)
%   D(:,:,i) = a(i)*eye(d,d)
%
%   Copyright (c) 2016 Gabriel Peyre


D = zeros(d,d,size(a,1),size(a,2));
a = reshape(a, [1 1 size(a)]);
for k=1:d   
    D(k,k,:,:) = a;
end

end