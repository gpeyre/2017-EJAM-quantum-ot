function x = expM(x)

% expM - matrix exponential
%
%   x = expM(x, use_fast);
%
%   You can set the global variable logexp_fast_mode:
%       - logexp_fast_mode=0 : matlab expm (slow)
%       - logexp_fast_mode=1 : explicit diagonalization (faster, only for 2x2 matrices)
%       - logexp_fast_mode=2 : explicit formula (even faster, only for 2x2 matrices)
%
%   Copyright (c) 2016 Gabriel Peyre

global logexp_fast_mode;
if isempty(logexp_fast_mode)
    logexp_fast_mode = 1;
end

if size(x,1)==1
    % 1D case
    x = exp(x);
    return;
end

switch logexp_fast_mode
    case 0
        % slow
        x = tensor_operation(x, @expm, 0);
    case 1
        % diagonalization
        x = tensor_operation(x, @exp, 1);
    case 2
        if size(x,1)~=2
            error('Works only for 2x2 matrices');
        end
        % explicit 2x2
        a = x(1,1,:,:); b = x(1,2,:,:); d = x(2,2,:,:);
        S = exp((a+d)/2);
        D = 1/2 * sqrt( (a-d).^2 + 4*abs(b).^2 );
        CD = cosh(D); SD = sinh(D)./D;
        L = (a-d) .* SD / 2;
        x(1,1,:,:) = S .* ( CD + L );
        x(2,2,:,:) = S .* ( CD - L );
        x(1,2,:,:) = S .* b .* SD;
        x(2,1,:,:) = conj(x(1,2,:,:));
    case 4
        x = tensorExp2x2(x);
    otherwise
        error('Unavailable');
end

% fprintf('C++ error for exp:  %g\n',norm(x(:)-xx(:),'fro'));
% save x.mat

end
