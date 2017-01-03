function Q = tensor_transp(T)

% tensor_transp - transpose a tensor field
%
%   Q = tensor_transp(T)
%
%   Q(:,:,i,j) = Q(:,:,i,j)'
%
%   Copyright (c) 2016 Gabriel Peyre

Q = zeros(size(T,2), size(T,1), size(T,3), size(T,4));

for i=1:size(T,1)
    for j=1:size(T,2)
        Q(j,i,:,:) = T(i,j,:,:);
    end
end

end
