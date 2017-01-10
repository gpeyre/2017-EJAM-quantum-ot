% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function y = mat_mul_bulk3(A,x)
    %{
    % debug verification:
    tic
    y = zeros(size(x));
    for i = 1:size(x,2)
        y(:,i) = A(:,:,i)*x(:,i);
    end
    toc
    %}

    %tic
    if size(x,2) == 1
        y = A*x;
        return
    end
    A = permute(A, [3 1 2]);
    x = x';
    y = [vec(A(:,1,1)).*x(:,1) + vec(A(:,1,2)).*x(:,2) + vec(A(:,1,3)).*x(:,3) ...
         vec(A(:,2,1)).*x(:,1) + vec(A(:,2,2)).*x(:,2) + vec(A(:,2,3)).*x(:,3) ...
         vec(A(:,3,1)).*x(:,1) + vec(A(:,3,2)).*x(:,2) + vec(A(:,3,3)).*x(:,3)]';
    %toc
end