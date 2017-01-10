% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function B = mat_inv_bulk2(A)
    %{
    tic
    y = zeros(size(x));
    for i = 1:size(x,2)
        y(:,i) = A(:,:,i)*x(:,i);
    end
    toc
    %}

    %tic
     
    detA = (A(1,1,:).*A(2,2,:) - A(1,2,:).*A(2,1,:));
    B = zeros(size(A));
    B(1,1,:) = A(2,2,:);
    B(2,1,:) = -A(2,1,:);
    B(1,2,:) = -A(1,2,:);
    B(2,2,:) = A(1,1,:);
    B = B ./ repmat(detA, [2 2 1]);
    %toc
end