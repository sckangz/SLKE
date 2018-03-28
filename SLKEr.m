function [result]= SLKEr( K,D,y,gamma,mu)
%initialization
n=length(unique(y));
maxIter=500;
nn=length(y);
Y1=zeros(nn);
Y2=zeros(nn);
Z=eye(nn);
%Z=rand(n);
W=Z;

%main function
for i=1:maxIter
    Zold=Z;
J=(D*W*W'*D'+mu*eye(nn))\(mu*Z+Y1+D*W*K');
% J=J-diag(diag(J));
 J(find(J<0))=0;
W=(D'*J*J'*D+mu*eye(nn))\(mu*Z+Y2+D'*J*K);
% W=W-diag(diag(W));
 W(find(W<0))=0;  
 
H=(W-Y1/mu+J-Y2/mu)/2;
%Z=max(abs(H)-gamma/(mu*2),0).*sign(H);
[U S V] = svd(H, 'econ');
       
        diagS = diag(S);
        svp = length(find(diagS > gamma/(2*mu)));
        diagS = max(0,diagS - gamma/(2*mu));
        
        if svp < 0.5 %svp = 0
            svp = 1;
        end
        Z = U(:,1:svp)*diag(diagS(1:svp))*V(:,1:svp)';

Z(find(Z<0))=0;  

Y1=Y1+mu*(Z-J);
Y2=Y2+mu*(Z-W);

 mu=mu*1.25;

if((i>5)&(norm(Z-Zold,'fro') < norm(Zold,'fro') * 1e-7))  
        break
    end
end

L=(Z+Z')/2;
actual_ids = spectral_clustering(L, n);

result=ClusteringMeasure(actual_ids ,y);