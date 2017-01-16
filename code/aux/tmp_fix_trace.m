
%%
% Fix distribution of traces to some input --- NOT USED.

if 0
    % reference histogram for eigenvalues
    sigma = .1;
    H = linspace(vmin,1,N).^2;
    H = vmin + exp(-linspace(-1,1,N).^2/(2*sigma^2));
    for k=1:m
        [e1,e2,l1,l2] = tensor_eigendecomp(nu{k});
        l1 = max(l1,vmin); l2 = max(l2,vmin);
        %
        Tr = l1+l2;
        Tr1 = hist_eq(Tr, H);
        l1 = l1 .* Tr1 ./ Tr; l2 = l2 .* Tr1 ./ Tr;
        %
        nu{k} = tensor_eigenrecomp(e1,e2,l1,l2);
    end
    
    for s=1:2
        [~,~,L1{s},L2{s}] = tensor_eigendecomp(mu{s});
        L1{s} = sort(L1{s});
        L2{s} = sort(L2{s});
    end
    for k=1:m
        t = (k-1)/(m-1);
        [e1,e2,l1,l2] = tensor_eigendecomp(nu{k});
        l1 = hist_eq(l1, (1-t)*L1{1} + t*L1{2});
        l2 = hist_eq(l2, (1-t)*L2{1} + t*L2{2});
        nu{k} = tensor_eigenrecomp(e1,e2,l1,l2);
    end
end