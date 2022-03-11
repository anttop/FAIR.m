function L=operator_norm(K,Kadj,n)
%computes the spectral norm of an operator via a potent method
    if length(n)==1
        u = randn(n,1); u = u/norm(u);
    else
        u = n;
        u = u/norm(u);
    end
    e = [];
    for i=1:30
        v = K(u);
        e(end+1) = sum(v(:).*v(:));
        u = v/norm(v(:));
        v= Kadj(u);
        e(end+1) = sum(v(:).*v(:));
        u = v/norm(v(:));
    end
    L = sqrt(e(end));