function op = load_derivative_2d()

% load_derivative_2d - load 2D image derivatives
%
%   op = load_derivative_2d();
%
%   op.Hess implements hessian.
%
%   Copyright (c) 2016 Gabriel Peyre


% hesian derivatives
op.dXf = @(f)f([2:end end],:) - f;
op.dXb = @(f)f - f([1 1:end-1],:);
op.dXc = @(f)( f([2:end end],:) - f([1 1:end-1],:) )/2;
op.dXX = @(f)op.dXf(op.dXb(f));
op.dYc = @(f)op.dXc(f')';
op.dYY = @(f)op.dXX(f')';
op.dXY = @(f)op.dXc(op.dYc(f));
op.Hess = @(f)permute( reshape(cat(3, op.dXX(f), op.dXY(f), op.dXY(f), op.dYY(f)), [size(f,1) size(f,1) 2 2]), [3 4 1 2]);

end