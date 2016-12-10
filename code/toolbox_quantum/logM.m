function y = logM(x)

% logM - matrix logarithm
%
%   x = logM(x);
%
%   You can set the global variable logexp_fast_mode:
%       - logexp_fast_mode=0 : matlab expm (slow)
%       - logexp_fast_mode=1 : explicit diagonalization (faster, only for 2x2 matrices)
%
%   Copyright (c) 2016 Gabriel Peyre

global logexp_fast_mode;
if isempty(logexp_fast_mode)
    logexp_fast_mode = 1;
end

if size(x,1)==1
    % 1D case
    y = log(x);
    return;
end

switch logexp_fast_mode
    case 0
        y = tensor_operation(x, @logm, 0);
    case {1,2}
        y = tensor_operation(x, @log, 1);
    case 3        % pbm?
        % explicit 2x2
        a = x(1,1,:,:); b = x(2,2,:,:); c = x(1,2,:,:);
        % eigenvalue
        D = sqrt( (a-b).^2 + 4*c.^2 ) / 2;
        T = (a+b) / 2;
        S = sign(a-b);
        % apply op to eigenvalues
        u = log(T+S.*D); v = log(T-S.*D);
        % angle of eigenvectors
        % m = 2*c ./ (a-b);
        % phi = atan(m); % 2*theta
        % phi(isnan(theta)) = 0; % update this
        % Cphi = 1 ./ sqrt(1+m.^2);
        % Sphi = m ./ sqrt(1+m.^2);
        Cphi = abs(a-b)./(2*D);
        Sphi = 2*c.*S./(2*D);
        C2 = (Cphi+1)/2; CS = Sphi/2;
        % remap
        uv = u-v;
        y(1,1,:,:) = uv.*C2+v;
        y(2,2,:,:) = -uv.*C2+u;
        y(1,2,:,:) = -uv.*CS;
        y(2,1,:,:) = y(1,2,:,:); % JUSTIN CHANGED RHS FROM x TO y
    case 4
        y = tensorLog2x2(x);
    otherwise
        error('Unavailable');
end

% yy = tensorLog2x2(x);
% err = norm(y(:)-yy(:),'fro');
% fprintf('C++ error for log:  %g\n',err);
% if err > 1e-9
%     save x.mat
%     error('log failed!');
% end

end