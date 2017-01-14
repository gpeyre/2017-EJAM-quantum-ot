% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function R = tex_reverse_transport(fdata, transport, solution)
    STEP = 2;

    % Get the exemplar variables
    l = solution.step{STEP}.vars.l;
	V = size(l,2);      % number of variables
    M = sqrt(size(l,1));    % exemplar pixel size

    % Size of the final image
    S = M*fdata.par.B;

    % Derivative matrices
    alpha = 0.75;
    Dx1 = conv2d_mat([-1 1], [M M], 'replicate'); 
    Dy1 = conv2d_mat([1;-1], [M M], 'replicate');
    
    % Pixel values and the derivatives, concatenated (this is what will be
    % transported)
    l1 = cat(3, l, Dx1*l, Dy1*l);
    
    B = fdata.par.B;        % number of blocks
    R = cell(fdata.par.B);
    
    % Transport each tile in turn, collect them to a big cell array
    for i = 1:B(1)
        i
        for j = 1:B(2)
            r = zeros(M^2,V);
            % Get the reverse transport mapping
            P = double(transport.P_m2b(:,:,i,j));
            % Clean it up a bit
            P = tex_filter_transport(P,4);
            % Then just do the permutation
            r = l1(P, :,:);
            r = reshape(r, [M M 3*V]);
            R{i,j} = r;
            
        end
    end
    
    % Fuse the cell array to a big image
    R = cell2mat(R);
    I2 = zeros([prod(S) V]);
    
    % For simultaneously enforcing the integrability of the gradient
    % reconstructed normals (doesn't matter that much really)
	Curl = [conv2d_mat([-1;1], S, 'replicate') conv2d_mat([-1 1], S, 'replicate')];
    
    b = cat(1, ...
        alpha*reshape(R(:,:,1:V), [], V), ...
        reshape(R(:,:,V+1:2*V), [], V), ...
        reshape(R(:,:,2*V+1:3*V), [], V));

    Dx2 = conv2d_mat([-1 1], S, 'replicate'); 
    Dy2 = conv2d_mat([1;-1], S, 'replicate');
    
    % The gradient operator for the full-size image (to be inverted)
    A = [alpha*speye(prod(S)); Dx2; Dy2];
    
    Z = sparse(size(A,1), size(A,2));

    % Augment with integrability prior
    AC = [A Z;
          Z A;
          5*Curl];
	bc = [b(:,1); b(:,2); zeros(prod(S),1)];

    ACTAC = AC'*AC;
    ACTbc = AC'*bc;

    % Reconstruct the normals from their transported gradients
    %MIC = ichol(ACTAC,struct('michol','on'));
    temp = pcg(ACTAC, ACTbc, 1e-5,200);%,MIC,MIC');
    
    I2(:,1:2) = reshape(temp, [], 2);
    
    % Same for the rest of variables (no curl constraint here)
    ATA = A'*A;
    ATb = A'*b;
    
    MIC = ichol(ATA,struct('michol','on'));
    
    for v = 3:V
        I2(:,v) = pcg(ATA, ATb(:,v), 1e-5,200,MIC,MIC');
    end
    
%    R0 = R;
    R = reshape(I2, [S V]);
    
%    imagesc(R0(:,:,1)-R(:,:,1))    % If we want to compare before/after
    
    
end