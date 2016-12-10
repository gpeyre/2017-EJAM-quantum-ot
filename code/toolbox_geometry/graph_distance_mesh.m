function D = graph_distance_mesh(V,F)

% graph_distance_mesh - geodesic distance on the graph of edges
%
%   D = graph_distance_mesh(V,F);
%
%   Copyright (c) 2016 Gabriel Peyre

N = size(V,2);

I = []; J = []; d = [];
for i=1:3
    j = mod(i,3)+1;
    I = [I F(i,:)]; 
    J = [J F(j,:)];
    X = V(:,F(i,:)); Y = V(:,F(j,:));
    d = [d sqrt(sum((X-Y).^2))];    
end

A = sparse(I,J,d,N,N);
A = full(max(A,A'));
A = A + speye(N,N)*1e-20;
A(A==0) = Inf;
D = FastFloyd(A);

end