function [error]= errormin(Y1,X,W,lambda,mu,type);

switch type
    case 1
%         g=Y1/mu+(X-X*W);
%         epsilon=lambda/mu;
%         [m,n]=size(g);
%         for i=1:m
%             for j=1:n
%                 if g(i,j)>epsilon
%                     E(i,j)=g(i,j)-episilon;
%                 elseif g(i,j)<-epsilon
%                     E(i,j)=g(i,j)+episilon;
%                 else
%                     E(i,j)=0;
%                 end
%             end
%         end
D = Y1/mu+(X-X*W);
E = zeros(size(D));
epsilon = lambda/mu;
DD = abs(D)-epsilon;
DD2 = DD.*sign(D);
ID = abs(D)>epsilon;
E(ID) = DD2(ID);
    case 2
        E=(Y1+mu*(X-X*W))/(2*lambda+mu); 
    case 21
        alpha=lambda/mu;
        G=X-X*W+Y1/mu;
        G1 = sqrt(sum(G.^2,1));
G1(G1==0) = alpha;
G2 = (G1-alpha)./G1;
E = G*diag((G1>alpha).*G2);

%         alpha=mu/lambda;
%         Q=X-X*W+Y1/mu;
%         [n,m]=size(Q);
%         for i=1:m
%             t=sum(Q(:,i).*Q(:,i));
%             if t>alpha
%             E(:,i)=(t-alpha)/t*Q(:,i);
%             else
%                 E(:,i)=zeros(n,1);
%             end
%         end
%             case 0
%         alpha=mu/lambda;
%         Q=X-X*W+Y1/mu;
%         [n,m]=size(Q);
%         for i=1:m
%             t=sum(Q(:,i).*Q(:,i));
%             if t>alpha
%             E(:,i)=(t-alpha)/t*Q(:,i);
%             else
%                 E(:,i)=zeros(n,1);
%             end
%         end
end
error=E;