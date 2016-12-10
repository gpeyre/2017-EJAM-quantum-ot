function A = tensor_operation(A, op, diag_op)

if nargin<3
    diag_op = 0;
end

% use faster expression
slow_mode = 0;

if diag_op==0
    for i=1:size(A,3)
        for j=1:size(A,4)
            A(:,:,i,j) = op(A(:,:,i,j));
        end
    end
else
    switch size(A,1)
        case 2
            [e1,e2,l1,l2] = tensor_eigendecomp(A);
            A = tensor_eigenrecomp(e1,e2,op(l1),op(l2));
        case 3
            n1 = size(A,3); n2 = size(A,4); 
            A = reshape(A, [3 3 n1*n2]);
            %
            [D,U] = Eigens3x3(A);
            for k=1:3
                D(k,1,:) = op(D(k,1,:));
            end
            % U*diag(D)*U'
            tprod = @(a,b)[...
                a(1,1,:).*b(1,1,:), a(1,1,:).*b(2,1,:), a(1,1,:).*b(3,1,:); ...
                a(2,1,:).*b(1,1,:), a(2,1,:).*b(2,1,:), a(2,1,:).*b(3,1,:); ...
                a(3,1,:).*b(1,1,:), a(3,1,:).*b(2,1,:), a(3,1,:).*b(3,1,:)];
            re = @(d)repmat(d, [3 1 1]);
            A = tprod(re(D(1,1,:)).*U(:,1,:), U(:,1,:)) + ...
                tprod(re(D(2,1,:)).*U(:,2,:), U(:,2,:)) + ...
                tprod(re(D(3,1,:)).*U(:,3,:), U(:,3,:));
            %
            A = reshape(A, [3 3 n1 n2]);
            
        otherwise
            error('Only 2D and 3D array supported');            
    end
end
    
end
