function M = perform_blurring(M, sigma, options)

% perform_blurring - gaussian blurs an image
%
%   M = perform_blurring(M, sigma, options);
%
%   M is the original data
%   sigma is the width of blurs (in pixels)
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;

if iscell(M)
    for i=1:length(M)
        M{i} = perform_blurring(M{i}, sigma, options);
    end
    return;
end

if sigma==0
    return;
end

if size(M,3)>1
    for i=1:size(M,3)
        for j=1:size(M,4)
            M(:,:,i,j) = perform_blurring(M(:,:,i,j), sigma, options);
        end
    end
    return;
end

n = max(size(M));

eta = 4;
p = round((sigma*eta)/2)*2+1;
p = min(p,round(n/2)*2-1);

A = [1 1];
if size(M,1)==1 || size(M,2)==1
    A = 1; % 1D.
end
h = compute_gaussian_filter(p*A,sigma/(4*n),n*A);
M = perform_convolution(M, h, options);

end


function f = compute_gaussian_filter(n,s,N);

% compute_gaussian_filter - compute a 1D or 2D Gaussian filter.
%
%   f = compute_gaussian_filter(n,s,N);
%
%   'n' is the size of the filter, odd for no phase in the filter.
%       (if too small it will alterate the filter).
%       use n=[n1,n2] for a 2D filter
%   's' is the standard deviation of the filter.
%   'N' is the size of the big signal/image (supposed to lie in [0,1] or [0,1]x[0,1]).
%       use N=[N1,N2] for a 2D filter.
%
%   The equation (in 1D) is
%       f[k] = exp( -(x(k)^2/(2*s^2)) );
%   where x span [-1/2,1/2].
%
%   The filter is normalised so that it sums to 1.
%
%   Copyright (c) 2004 Gabriel Peyre

nd = 1;
if length(n)>1 & n(2)>1
    nd = 2;
end

if nd==2 && length(s)==1
    s = [s s];
end
if nd==2 && length(N)==1
    N = [N N];
end

if nd==1
    f = build_gaussian_filter_1d(n,s,N);
else
    f = build_gaussian_filter_2d(n,s,N);
end

end

% build_gaussian_filter_2d - compute a 2D Gaussian filter.
%
%   f = build_gaussian_filter_2d(n,s,N);
%
%   'n' is the size of the filter, odd for no phase in the filter.
%       (if too small it will alterate the filter).
%   's' is the standard deviation of the filter.
%   'N' is the size of the big image (supposed to lie in [0,1]x[0,1]).
%
%   The filter is normalised so that it sums to 1.
%
%   Copyright (c) 2004 Gabriel Peyre

function f = build_gaussian_filter_2d(n,s,N)

if nargin<2
    error('Not enough arguments.');
end
if nargin<3
    N = n;
end

if length(N)==1 || N(1)==1
    N = N(:); N = [N N];
end
if length(s)==1 || s(1)==1
    s = s(:); s = [s s];
end

if s<=0
    f = zeros(n);
    f(round((n-1)/2),round((n-1)/2)) = 1;
    return;
end

x = ( (0:n(1)-1)-(n(1)-1)/2 )/(N(1)-1);
y = ( (0:n(2)-1)-(n(2)-1)/2 )/(N(2)-1);
[Y,X] = meshgrid(y,x);
f = exp( -(X.^2/ (2*s(1)^2)) - (Y.^2/ (2*s(2)^2)) );
f = f / sum(f(:));

end

% build_gaussian_filter_1d - compute a Gaussian filter.
%
%   f = build_gaussian_filter_1d(n,s,N);
%
%   Copyright (c) 2004 Gabriel Peyre

function f = build_gaussian_filter_1d(n,s,N)

if nargin<2
    error('Not enough arguments.');
end
if nargin<3
    N = n;
end

if s<=0
    f = zeros(n,1);
    f(round((n-1)/2)) = 1;
    return;
end

x = ( (0:n-1)-(n-1)/2 )/(N-1);
f = exp( -x.^2/(2*s^2) );
f = f / sum(f(:));

end


function y = perform_convolution(x,h, bound)

% perform_convolution - compute convolution with centered filter.
%
%   y = perform_convolution(x,h,bound);
%
%   The filter 'h' is centred at 0 for odd
%   length of the filter, and at 1/2 otherwise.
%
%   This works either for 1D or 2D convolution.
%   For 2D the matrix have to be square.
%
%   'bound' is either 'per' (periodic extension)
%   or 'sym' (symmetric extension).
%
%   Copyright (c) 2004 Gabriel Peyre

if nargin<3
    bound = 'sym';
end

if isstruct(bound)
    bound = getoptions(bound, 'bound', 'sym');
end

if not(strcmp(bound, 'sym')) && not(strcmp(bound, 'per'))
    error('bound should be sym or per');
end

if iscell(x)
    for i=1:length(x)
        y{i} = perform_convolution(x{i},h, bound);
    end
    return;
end

if size(x,3)>1 && size(x,3)<4
    % for color images
    y = x;
    for i=1:size(x,3)
        y(:,:,i) = perform_convolution(x(:,:,i),h, bound);
    end
    return;
end

if size(x,3)>1
    error('Not yet implemented for 3D array, use smooth3 instead.');
end

n = size(x);
p = size(h);

bound = lower(bound);

nd = ndims(x);
if size(x,1)==1 || size(x,2)==1
    nd = 1;
end
if nd==1
    n = length(x);
    p = length(h);
end

if strcmp(bound, 'sym')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % symmetric boundary conditions

    d1 = floor( p/2 );  % padding before
    d2 = p-d1-1;            % padding after

    if nd==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% 1D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x = x(:); h = h(:);
        xx = [ x(d1:-1:1); x; x(end:-1:end-d2+1) ];
        y = conv(xx,h);
        y = y(p:end-p+1);
    elseif nd==2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % double symmetry
        xx = x;
        xx = [ xx(d1(1):-1:1,:); xx; xx(end:-1:end-d2(1)+1,:) ];
        xx = [ xx(:,d1(2):-1:1), xx, xx(:,end:-1:end-d2(2)+1) ];

        y = conv2(xx,h);
        y = y( (2*d1(1)+1):(2*d1(1)+n(1)), (2*d1(2)+1):(2*d1(2)+n(2)) );
    end

else

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % periodic boundary conditions

    if p>n
        error('h filter should be shorter than x.');
    end
    d = floor((p-1)/2);
    if nd==1
        x = x(:); h = h(:);
        h = [ h(d+1:end); zeros(n-p,1); h(1:d) ];
        y = real( ifft( fft(x).*fft(h) ) );
    else
        h = [ h(d(1)+1:end,:); zeros( n(1)-p(1),p(2) ); h( 1:d(1),: ) ];
        h = [ h(:,d(2)+1:end), zeros(n(1),n(2)-p(2)), h(:,1:d(2)) ];
        y = real( ifft2( fft2(x).*fft2(h) ) );
    end

end

end
