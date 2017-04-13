d = 3;
S = @(x)x*x';
D = @(P,Q)trace( P+Q-2*expm(logm(P)/2+logm(Q)/2) )^(1/2);
k = 0;
p={}; q={}; r={};v=[];
while true
    k = k+1;
    P = S(randn(d));
    Q = S(randn(d));
    R = S(randn(d));
    v(end+1) = ( D(P,Q)+D(Q,R) )/D(P,R) - 1;
    if v(end)<0
        warning(['Problem:' num2str(k)]);
        p{end+1}=P; q{end+1}=Q; r{end+1}=R;
    end
    
end

return;

A = @(P,Q)sqrtm( sqrtm(P)*Q*sqrtm(P) );