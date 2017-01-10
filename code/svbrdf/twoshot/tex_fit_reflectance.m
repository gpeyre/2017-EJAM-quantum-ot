% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function [vars,R] = tex_fit_reflectance(recons, fdata, coords, initial, maxiters)

    % global constant stuff for the optimizer (e.g. input data, useful
    % matrices, settings, etc)
    glob = struct;
    glob.S0 = size(recons);
        
    nth = 1;
    glob.recons = recons(:,:,:,1:nth:end,1:nth:end);

    glob.S0 = size(glob.recons);
    
    % Collapse the 2D image and 2D block dims to linear lists, because
    % every residual will rather operate on such form. So now the indices
    % are [pixels channels blocks], e.g. [64*64 3 16*16] = [4096 3 256]
    glob.S = [prod(glob.S0(1:2)) glob.S0(3) prod(glob.S0(4:5))];
    glob.recons = reshape(glob.recons, glob.S);

    glob.initial = glob.recons;
    
    glob.DjImg = sparse_clip_boundary(glob.S0(1:2), [0 1]) * ...
        conv2d_mat([-1 1], glob.S0(1:2), 'replicate');
    glob.DiImg = sparse_clip_boundary(glob.S0(1:2), [1 0]) * ...
        conv2d_mat([1;-1], glob.S0(1:2), 'replicate');
 
    glob.n_var = 9; % number of variables per pixel
    
    Z = sparse(glob.S(1), glob.S(1));
    I = speye(glob.S(1), glob.S(1));
    
    % Rotate a vector field by 90 degrees (for curl)
    glob.R = [...
        Z I;
        -I Z;
        ];
    
    glob.Grad = [conv2d_mat([-1 1], glob.S0(1:2), 'replicate'); conv2d_mat([1;-1], glob.S0(1:2), 'replicate')];
    glob.Curl = [conv2d_mat([-1;1], glob.S0(1:2), 'replicate') conv2d_mat([-1 1], glob.S0(1:2), 'replicate')];
    
    
    % Optimization variables
    vars = struct;

    % Initial guess is the current reconstruction
    if nargin < 4 || numel(initial) == 0
        % if there's no initial guess supplied, make one
        vars.g = [0 0 -1 0.4];

        nn = glob.S(1);

        RGB_guess = log(mean(reshape(fdata.img_double.flash,[],3),1));
        vars.l = repmat([0 0  RGB_guess -1 -3 -3 0.01], [nn 1]);    % XXX -4 for -2 in spec size
        vars.l = vars.l + rand(size(vars.l))*0.01;
    else
        % ... otherwise use the supplied one
        vars = initial;
    end
    
    % Rendering coordinate system
    glob.cX = coords(:,:,1);
    glob.cY = coords(:,:,2);
    
    glob.vars_l_size = size(vars.l);
    glob.vars_g_size = size(vars.g);
    
    rng('default');
    rng(123425);

    % For historical reasons we call the full tree rebuild for every
    % iteration, it's not a big deal but a bit wasteful...
    
    % Form the residual function for L-M
    residual = @(Y) pass_nargout(@(X) eval_tree(build_tree_full(glob), glob, X), @flatten_tree, unpack_vars(glob,Y));
    outfun = @(x,o,s) render_show(glob, x,o); % function that shows the current result between iterations
    options = optimoptions('lsqnonlin','Jacobian','on','Algorithm','levenberg-marquardt','OutputFcn', outfun, 'MaxIter', maxiters, 'InitDamping', 1);%,'ScaleProblem','Jacobian');
    % Optimize
    x = lsqnonlin(residual, pack_vars(vars), [],[],options);

    % Done, x now has the solution. Re-render the lumitexel for use in H-B.
    vars = unpack_vars(glob, x);
    tree_render = build_tree_render(glob, false);
    R = eval_tree(tree_render, glob, vars, []);
    R = reshape(R{1}, glob.S0);
        
    
end

function tree = build_tree_full(glob)

    tree = struct;
    tree.F = @tree_nop;

    smooth_w = 4*[0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]';

    % Data fidelity residual
    tree.child{1} = build_tree_render(glob, true);
    % Normal map curl residual
    tree.child{2} = build_tree_normals(glob);
    % Smoothness residual
    tree.child{3} = build_tree_smooth_spatial(glob, smooth_w*0.05);

    % Pointwise prior expectations
    Q = struct;
    Q.g = [];
    Q.l = repmat(-[0 0 0 0 0 -6 -6 -6 0], [glob.S(1) 1]);

    % Prior weights
    w = [1 1 2 2 2 4 0.5 0.5 0.5]';
    W = spdiags(w, 0, numel(w), numel(w));

    % Pointwise prior residual
    tree.child{4}.F =  @(a,b,c) tree_chain(a,b,c, ...
                { ...
                @tree_strip_globals, ...
                @(x,y,z) tree_add(x,y,z, Q), ...
                @(x,y,z) tree_permute(x,y,z, [2 1]), ...
                @(x,y,z) tree_matrix_prod(x,y,z, W), ...
                @(x,y,z) tree_mul_scalar(x,y,z, 0.01), ...
                });
end


function tree = build_tree_render(glob, resid)
    tree = struct;
    Z = struct;
    Z.g = [];
    Z.l = -glob.initial;
	
    
    if resid == true
        % If we asked for a comparison against input data...
    	tree.F = ...
                @(a,b,c) tree_chain(a,b,c, ...
                    { ...
                    @tree_render, ...
                    @tree_strip_globals, ...
                    @(x,y,z) tree_clamp(x,y,z, [-100 1]), ...
                    @(x,y,z) tree_add(x,y,z,Z), ...
                    @(x,y,z) tree_huber(x,y,z, 0.1), ...
                    @(x,y,z) tree_mul_scalar(x,y,z, sqrt(150/glob.S(3))), ...
                    });
    else
        % ... but if we just wanted the picture:
    	tree.F = ...
                @(a,b,c) tree_chain(a,b,c, ...
                    { ...
                    @tree_render, ...
                    @tree_strip_globals, ...
                    @(x,y,z) tree_clamp(x,y,z, [0 1]), ...
                    });
    end
end


function tree = build_tree_normals(glob)
    tree = struct;
	
	tree.F = ...
        @(a,b,c) tree_chain(a,b,c, ...
        { ...
            ... % Get the first two components, flattened
            @(x,y,z) tree_extract_head(x,y,z, glob.S(1)*2), ...
            ... % and compute the curl.
            @(x,y,z) tree_matrix_prod(x,y,z, glob.Curl), ...
            @(x,y,z) tree_mul_scalar(x,y,z, 5), ...
        });
end


function tree = build_tree_smooth_spatial(glob, w)
    tree = struct;
	
    W = spdiags(w, 0, numel(w), numel(w));
    
	tree.F = ...
        @(a,b,c) tree_chain(a,b,c, ...
        { ...
            @(x,y,z) tree_matrix_prod(x,y,z, [glob.DiImg; glob.DjImg]), ...
            @(x,y,z) tree_permute(x,y,z, [2 1]), ...
            @(x,y,z) tree_matrix_prod(x,y,z, W), ...
            @(x,y,z) tree_huber(x,y,z, 0.1), ...
            @(x,y,z) tree_mul_scalar(x,y,z, sqrt(96/glob.S0(1))), ... 
        });
end


% Display some intermediate data during optimization (not a very good
% selection of data at that, but whatever!)
function stop = render_show(glob,x,opv)
    stop = false;
    
    disp('iter done.');
    
    if exist('_break', 'file')
        stop = true;
    end
    
    
    clf;
  
    vars = unpack_vars(glob, x);

    vars.g

    
    tree = build_tree_render(glob, false);
    R = eval_tree(tree, glob, vars, []);

    R = reshape(R{1}, glob.S0);

    % Show a few different lighting directions...
    n = 1;
    is = [1 round(size(R,4)/2) min(size(R,4),size(R,5))];
    for i = 1:3
        for j = 1:3
            subplot(6,3,n);
            
            imagec(R(:,:,:,is(i),is(j)));
            
            n = n+1;
        end
    end
    
    % Show the rendered lumitexels at a few different locations...
    for i = 1:3
        for j = 1:3
            subplot(6,3,n);
            
            imagec(permute(squeeze(R(3*i+32,3*j+32,:,:,:)), [2 3 1]));
            
            
            n = n+1;
        end
    end
    
    % Show some kind of an encoding of the current variables
    L = reshape(vars.l, [glob.S0(1:2) size(vars.l, 2)]);
    colormap gray;
    subplot(6,3,10);
    imagec(exp(L(:,:,3:5)));
    subplot(6,3,11);
    imagesc(exp(L(:,:,6))*64); % 9
    subplot(6,3,12);
    imagesc(L(:,:,7)); % 11
    subplot(6,3,13);
    imagesc(L(:,:,8));
    subplot(6,3,14);
    imagesc(L(:,:,9));
    
    subplot(6,3,18);
    size(vars.l)
    V = reshape(vars.l, [glob.S0(1:2) size(vars.l,2)]);
    imagesc(max(0,min(1,cat(3, 0.5+V(:,:,1)-mean(vec(V(:,:,1))), 0.5+V(:,:,2)-mean(vec(V(:,:,2))), 1+0*V(:,:,1)))));
    
    drawnow;

end


%% Hierarchical evalo-differentiator

% This is the system that drives the residual and derivative evaluation.
% eval_tree traverses the hierarchical residual tree and returns a
% hierarchical list of residuals, and Jacobians if such an output parameter
% is requested (if not, derivative computations are automatically skipped).

% See the setup functions above to see how the hierarchical residuals are
% constructed. This is basically similar to forward mode automatic 
% differentation in massive bulk matrix form, allowing for the kind of 
% structured sparsity we need. The system is hard-coded with the assumption
% that there is a multidimensional array of bulk variables (x.l), typically
% per-pixel, and a vector of global variables (x.g). This results in an
% arrow-head structure for J'J. It is assumed that no global variable ever
% depends on local variables, as this wasn't required in our algorithm (and
% also has the potential to generate _very_ much fill in J'J). The
% corresponding dependencies are encoded in J.ll, J.lg and J.gg, which are
% essentially blocks of J with lower-left missing.

% The user must supply the chain rule Jacobian
% multiplication routines by hand, and there is obviously no verification
% that the supplied derivatives match the functions themselves; see
% example implemtations of multiple basic operations below.

% Note that the Jacobians are accumulated in transpose matrices, as
% horizontal concatenation of sparse matrices is very much significantly 
% faster than vertical.

function [x_,J_] = eval_tree(tree, glob, x, varargin)
    if (nargout == 1)
        %% Jacobian not needed.
        
        % Evaluate the branch function at current position
        disp(sprintf('Tree eval: \t%s', func2str(tree.F)));
        xx = tree.F(glob, x, []);
        
        if isfield(tree,'child') == 0 || numel(tree.child) == 0
            % If this was a leaf, just return the value, flattened
            x_ = {xx.l; xx.g};
        else
            % Otherwise pass the value to all children and concatenate the
            % return values
            x_ = cell(numel(tree.child),1);
            for b = 1:numel(tree.child)
                x_{b} = eval_tree(tree.child{b}, glob, xx);
            end
        end
    else
        %% Jacobian needed.
        J = [];
        
        % Use the J if one passed; otherwise make a seed identity
        if numel(varargin) > 0
            J = varargin{1};
        else
            J = make_J(x);
        end
            
        % Evaluate and get the Jacobian too.
        disp(sprintf('Tree eval J: \t%s', func2str(tree.F)));

        [xx, JJ] = tree.F(glob, x, J);
        
        % NOTE! F must return J_F * J, not just J_F! This is because
        % sometimes J_F's action is easier and  more efficient to express 
        % in non-matrix form (e.g. permutation or row choice).
        
        if isfield(tree,'child') == 0 || numel(tree.child) == 0
            % If it's a leaf, return the result
            
            % Concatenate the locals and globals
            x_ = {xx.l; xx.g};
            % Open up the block-structured J, so it will be in suitable
            % form for vertical concatenation; no function will receive
            % this J as a parameter anymore.
            % NOTE: transposed for catenation efficiency
            J_ = {[JJ.ll JJ.lg]';
                [sparse(size(JJ.gg,1),size(JJ.ll,2)) JJ.gg]'}';
        else
            % Otherwise, concatenate results of all children
            
            J_ = cell(1,numel(tree.child)); % (transposed!)
            x_ = cell(numel(tree.child),1); 
            for b = 1:numel(tree.child)
                [xx_, JJ_] = eval_tree(tree.child{b}, glob, xx, JJ);
                x_{b} = xx_;
                J_{1,b} = JJ_; % (transposed!)
            end
        end
    end
end

% The hierarchical list of residuals and Jacobians is collapsed to a 
% matrix and a vector with this function.
function [x_,J_] = flatten_tree(x, varargin)
    %disp('flattening tree.');
    if nargout == 1
        %% Without Jacobians
        % Special case: empty leaf
        if numel(x) == 0
        %    x_ = [];
        %    return;
            warning('empty leaf');  % not sure if it matters
        end

        if ~iscell(x{1})
            % If this is a leaf, combine all into one.
            x_ = x;
            for i = 1:numel(x_)
                x_{i} = x_{i}(:);   % flatten the matrix
            end
            x_ = cell2mat(x_);
            return;
        else
            % If it's not a leaf, combine results from children
            x_ = cell(numel(x),1);
            for i = 1:numel(x)
                x_{i} = flatten_tree(x{i});
            end
            x_ = cell2mat(x_);
        end
        
    else
        %% With Jacobians
        % Special case: empty leaf
        if numel(x) == 0
            warning('empty leaf');  % not sure if it matters
        end

        J = varargin{1};
        
        if ~iscell(x{1})
            % If this is a leaf, combine all into one.
            x_ = x;
            for i = 1:numel(x_)
                x_{i} = x_{i}(:);   % flatten the matrix
            end
            x_ = cell2mat(x_);
            %J
            J_ = cell2mat(J);
            return;
        else
            % If it's not a leaf, combine results from children
            x_ = cell(numel(x),1);  % (transposed!)
            J_ = cell(1,numel(x));
            for i = 1:numel(x)
                [x_{i}, J_{i}] = flatten_tree(x{i},J{i});
            end
            x_ = cell2mat(x_);
            J_ = cell2mat(J_);
        end
        
    end
end

% A helper that gets around some Matlab input/output inference problem
% we had...
function [x_,J_] = pass_nargout(f1, f2, x)
    if nargout == 1
        x_ = f2(f1(x));
    else
        [xx,JJ] = f1(x);
        [x_,J_] = f2(xx,JJ);
        J_ = J_';
    end 
end


% A convenience function for constructing linear chains of
% compositions.
function [x_, J_] = tree_chain(glob, x, J, funs)
    % Recursively evaluate all the function handles in funs, starting from
    % the first (at the deepest recursion level) and feeding its output 
    % onward.
    if nargout == 1
        % Plain evaluation
        xx = x;
        if numel(funs) > 1
            xx = tree_chain(glob, x, J, funs(1:end-1));
        end
        %disp(sprintf('\tchain eval: \t%s', func2str(funs{end})));
        x_ = funs{end}(glob, xx, J);
    else
        % Evaluation with Jacobians
        xx = x;
        JJ = J;
        if numel(funs) > 1
            [xx, JJ] = tree_chain(glob, x, J, funs(1:end-1));
        end
        %disp(sprintf('\tchain eval J: \t%s', func2str(funs{end})));
        [x_, J_] = funs{end}(glob, xx, JJ);
    end
end


% Helper for multiplication of two block matrices of form
% [ ll  lg
%    0  gg ]
% (where lg means "derivatives of local vars wrt. global vars, etc.; 
% note that we don't allow globals to depend on locals)
function J = Jmul(K, L)
    J.ll = K.ll * L.ll;
    J.lg = K.ll * L.lg + K.lg * L.gg;
    J.gg = K.gg * L.gg;
end

% Initial seed (identity matrix)
function J = make_J(vars)
    J = struct;
    Nl = numel(vars.l);
    Ng = numel(vars.g);
    J.ll = spdiags(ones(Nl,1), 0, Nl, Nl);
    J.lg = zeros(Nl,Ng);
    J.gg = eye(Ng);
end


%% General tree utility component functions (not all are used currently)

function [x_, J_] = tree_nop(glob, x, J)
    x_ = x;
    if nargout > 1
        J_ = J;
    end
end

function [x_, J_] = tree_add(glob, x, J, y)
    x_ = x;
    x_.l = x_.l + y.l;
    x_.g = x_.g + y.g;
    
    if nargout > 1
        J_ = J;
    end
end

function [x_, J_] = tree_add_locals(glob, x, J, y)
    x_ = x;
    x_.l = x_.l + y;
    
    if nargout > 1
        J_ = J;
    end
end

function [x_, J_] = tree_mul_scalar(glob, x, J, s)
    x_ = x;
    x_.l = s * x_.l;
    x_.g = s * x_.g;
    
    if nargout > 1
        J_.ll = s * J.ll;
        J_.lg = s * J.lg;
        J_.gg = s * J.gg;
    end
end


function [x_, J_] = tree_permute(glob, x, J, p)
    s = size(x.l);
    
    x_ = x;
    
    x_.l = permute(x.l, p);
    
    if nargout > 1
        J_ = J;
        o = reshape(1:prod(s), s);
        o = vec(permute(o, p));
        
        J_.ll = J.ll(o,:);
        J_.lg = J.lg(o,:);
    end
end

function [x_, J_] = tree_flatten_left_dims(glob, x, J)
    s = [size(x.l) 1];  % Add a dummy-1 so we can also handle n*2 matrices
    
    x_ = x;
    x_.l = reshape(x.l, [s(1)*s(2) s(3:end)]);
    
    if nargout > 1
        J_ = J;
    end
end


function [x_, J_] = tree_flatten_all(glob, x, J)
    x_ = x;
    x_.l = x.l(:);
    
    if nargout > 1
        J_ = J;
    end
end

function [x_, J_] = tree_extract_slice_rightdim(glob, x, J, pos)
    s = size(x.l);
    
    x_ = x;

    % Number of elements in the all-but-right dimensions
    len = prod(s(1:end-1));
    
    x_.l = x.l(len*(pos-1)+1:len*pos);
    x_.l = reshape(x_.l, s(1:end-1));
    
    if nargout > 1
        J_ = J;
        
        J_.ll = J.ll(len*(pos-1)+1:len*pos,:);
        J_.lg = J.lg(len*(pos-1)+1:len*pos,:);
    end
end

function [x_, J_] = tree_sort(glob, x, J)
    x_ = x;

    [x_.l, ord] = sort(x.l(:));
        
    if nargout > 1
        J_ = J;
        
        J_.ll = J.ll(ord,:);
        J_.lg = J.lg(ord,:);
    end
end



function [x_,J_] = tree_matrix_prod(glob, x, J, A)
    
    x_ = x;
    x_.l = A * x.l;

    if nargout > 1
        J_ = J;
        
        DA = kron(speye(size(x.l,2)), A);
        J_.ll = DA * J.ll;
        J_.lg = DA * J.lg;

    end
end

function [x_, J_] = tree_reshape(glob, x, J, r)
    x_ = x;
    x_.l = reshape(x.l, r);
    
    if nargout > 1
        % The Jacobian cannot maintain multiple dimensions (due to 
        % sparsity), so this is a nop.
        J_ = J;
    end
end


function [x_, J_] = tree_extract_head(glob, x, J, n)
    x_ = x;
    x_.l = x.l((1:n)');
        
    if nargout > 1
        J_ = J;
        J_.ll = J.ll((1:n)',:);
        J_.lg = J.lg((1:n)',:);
    end
end

% The purpose of this often-called node is that we generally want to only
% use the globals as "parameters" for the render function; if they are not
% deleted, their values will be included as residuals, resulting in an 
% accidental zero-mean prior.
function [x_,J_] = tree_strip_globals(glob, x, J)
    x_.l = x.l;     % Pass the locals...
    x_.g = [];      % but delete the locals
    
    if nargout > 1
        J_.ll = J.ll;
        % Note that the dependencies of the locals on original global
        % remain.
        J_.lg = J.lg;
        J_.gg = zeros(0,size(x.g,1)); % Make a zero matrix of proper width
    end    
end


function [x_, J_] = tree_clamp(glob, x, J, range)
    x_ = x;

    to_min = x.l < range(1);
    to_max = x.l > range(2);
    
    x_.l(to_min) = range(1);
    x_.l(to_max) = range(2);
    
    if nargout > 1
        J_ = J;
        J_.ll(to_min(:),:) = 0;
        J_.ll(to_max(:),:) = 0;
        J_.lg(to_min(:),:) = 0;
        J_.lg(to_max(:),:) = 0;
    end
end


function [x_, J_] = tree_mul_pointwise(glob, x, J, ml, mg)
    x_ = x;
    x_.l = ml .* x_.l;
    x_.g = mg .* x_.g;
    
    if nargout > 1
        ML = spdiags(ml(:), 0, numel(ml), numel(ml));
        J_.ll = ML * J.ll;
        J_.lg = ML * J.lg;
        J_.gg = diag(mg) * J.gg;
    end
end


function [x_,J_] = tree_huber(glob, x, J, a)
    x_.l = hub(x.l, a);
    x_.g = hub(x.g, a);
    
    if nargout > 1
        nl = numel(x.l); %size(x.l,1);
        hub_Jl = spdiags(hub_J(x.l(:),a),0,nl,nl);
        J_.ll = hub_Jl * J.ll;
        J_.lg = hub_Jl * J.lg;
        
        J_.gg = diag(hub_J(x.g(:),a)) * J.gg;
    end
end


function x_ = hub(x,a)
    x_ = sqrt(2) * a * sign(x) .* sqrt(sqrt(1+x.^2/a^2)-1);
end

function J_ = hub_J(x,a)
    denom = a * sign(x) .* sqrt(2*x.^2/a^2+2).*sqrt(sqrt(x.^2/a^2+1)-1);
    J_ = x ./ denom;
    J_(denom == 0) = 1;
end


function [x_,J_] = tree_exp(glob, x, J, idx)
    x_ = x;
    J_ = J;

    bidx = zeros(size(x_.l));
    bidx(:,idx) = 1;
    idxes = bidx(:) ~= 0;
    
    x_.l(idxes) = exp(x_.l(idxes));
    
    if nargout > 1
        je = ones(numel(bidx),1);
        je(idxes) = x_.l(idxes);
        JE = spdiags(je, 0, numel(je), numel(je));
        J_.ll = JE * J.ll;
        J_.lg = JE * J.lg;        
    end
end


%%

% Tree node that takes in the per-pixel parameters and gives out the
% corresponding lumitexel arrays. If Jacobian is requested, it is computed
% with finite differences in bulk, and composed into the suitable 
% block-matrix form (as there are no inter-pixel dependencies in 
% rendering).
function [x_, J_] = tree_render(glob, x, J)
    % Render at current position
    x_ = struct;
    
    x_.l = tree_render_all(glob, x);
    
    x_.g = x.g; % Pass the globals through
  
    % Sparse diag helper
    spdiag = @(x) spdiags(x, 0, numel(x), numel(x));

    % Need to differentiate?
    if nargout > 1
       e = 0.000001;  % FD epsilon       
       J_ = struct;
       
       % Loop through the pointwise local variables and perturb each in
       % turn. 
       J_.ll = sparse(numel(x_.l), 0);%sparse(0,0);
       ll_cell = cell(1, size(x.l,2));
       for v = 1:size(x.l,2)
           disp(sprintf('Evaluating finite differences of local variable %i', v));
           x2 = x;
           x2.l(:,v) = x2.l(:,v) + e;

           % Evaluate finite difference...
           d = (tree_render_all(glob, x2) - x_.l) / e;
           
           % Shuffle the derivatives into their proper order...
           D = spconvert([(1:numel(d))', mod((0:numel(d)-1)', size(d,1))+1, d(:)]);
                      
           % Store the newest block in the Jacobian
           ll_cell{v} = D;
           
       end
       J_.ll = cell2mat(ll_cell);
    
       % Evaluate finite differences for each global variable, and cat them
       % to the end of J.
       J_.gg = eye(size(x.g,1));
       J_.lg = sparse(numel(x_.l), 0);
       lg_cell = cell(1, numel(x.g));
       for v = 1:numel(x.g)
           
           disp(sprintf('Evaluating finite differences of global variable %i', v));
           x2 = x;
           x2.g(v) = x2.g(v) + e;
           d = (tree_render_all(glob, x2) - x_.l) / e;
           %J_.lg = [J_.lg d(:)];
           lg_cell{v} = d(:);
       end
       J_.lg = cell2mat(lg_cell);

        J_ = Jmul(J_,J);
    end
   
end

% Call the actual rendering model
function R = tree_render_impl(glob, X, g)
    R = tex_render(vec(glob.cX)', vec(glob.cY)', X', g)';
end

function l = tree_render_all(glob, x)
    l = zeros(glob.S);
    for pid = 1:glob.S(1)
        % Extract the parameters for this particular pixel and render its
        % microfacet distribution
        l(pid,:,:) = tree_render_impl(glob, x.l(pid,:), x.g)';        
    end
end


%% 

% Functions for going back and forth with the "human-readable" .l and .g
% arrays, and the flattened list that L-M likes to work with...
function vars = unpack_vars(glob, x)
    vars.l = reshape(x(1:prod(glob.vars_l_size)), glob.vars_l_size);
    vars.g = [];
    if numel(x) > prod(glob.vars_l_size)
        vars.g = x(prod(glob.vars_l_size)+1:end);
    end
end

function x = pack_vars(vars)
    x = [vars.l(:);vars.g(:)];
end

