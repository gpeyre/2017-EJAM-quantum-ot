function [op,func,mat] = load_graph_operators(F, t, epsilon, options)

% load_graph_operators - load various graph-based operators
%
%   [op,func,mat] = load_graph_operators(F, t, epsilon);
%
%   F should be of size (3,P) where P is #faces.
%
%   t is the step-size to apply Vardhan formula.
%   epsilon is a safety parameter when inverting distances.
%
%   Copyright (c) 2015 Gabriel Peyre

options.null = 0;

n = max(F(:));

% stores parameters
mat.t = t;
mat.epsilon = epsilon;
mat.faces = F;

% compute edges
E = [F([1 2],:) F([2 3],:) F([3 1],:)];
E = unique([E E(2:-1:1,:)]', 'rows')';
mat.E = E;
p = size(E,2);

% gradient operator
mat.nabla = sparse( [1:p 1:p], [E(1,:) E(2,:)], [ones(1,p) -ones(1,p)] );

%% 
% Load scalar functions.

func.w  = getoptions(options, 'func_w', []);
func.wD = getoptions(options, 'func_wD', []);
func.h  = getoptions(options, 'func_h', []);
func.hD = getoptions(options, 'func_hD', []);

if isempty(func.w) || isempty(func.wD)
    func.w = @(s)1./(s+epsilon^2);
    func.wD = @(s)-1./(s+epsilon^2).^2; % derivative
end
if isempty(func.h) || isempty(func.hD)
    func.h = @(u)-1/t * log( u );
    func.hD = @(u)-1./(t*u); % derivative
end

% load operators functions
op.S = @(x)sum( (mat.nabla*x).^2, 2 );
op.W = @(S)full( sparse( E(1,:),E(2,:), func.w(S), n,n ) );
op.L = @(W)diag(sum(W,2))-W;
op.D = @(D)func.h(D);
op.H = @(L)inv(eye(n)+t*L);

% load matrices functions
mat.L = @(x)op.L(op.W(op.S(x)));
mat.D = @(x)op.D(op.H(mat.L(x)));

end