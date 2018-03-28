function [ X,T ] = DC(D,rho,T0)

[U,S,V] = svd(D,'econ');

for t = 1:1
    
    [ X,T1 ] = DCInner(S,rho,T0,U,V);
    err = norm(T1-T0,'fro')/norm(T0,'fro');
    if err < 1e-3
        break
    end
    T0 = T1;
end
T = T1;
end


function [ X,t ] = DCInner(S,rho,J,U,V)
lambda=1/2/rho;
S0 = diag(S);
grad=1./(1+J.^2);
t=max(S0-lambda*grad,0);
X=U*diag(t)*V';
end
