alpha=[ .0001  10 ];
beta=[.001  10];
% gamma=[1e-4  1e-3  .01];
%m=length(lambda)*length(mu);
load COIL20_1440n_1024d_20c;
load('COIL20_1440n_1024d_20c_kernel/COIL20_1440n_1024d_20c_kernel_gaussian_50_post_Sample-Scale.mat');
warning off
addpath('/root/permanent/qpc.txt')
best=zeros(1,3);
for i=1:length(alpha)
    for j=1:length(beta)
%       for ij=1:length(gamma);
fprintf('params %12.6f%12.6f\n',alpha(i),beta(j))
 % [res] =localmulticluster(X,K,y,alpha(i),beta(j),gamma(ij))
  [result]=unifiedcluster(K,y,alpha(i),beta(j))
% if (res(1)>best(1))
%     best=res;
%     alp=alpha(i);bta=beta(j);gam=gamma(ij);
%      save('C:\Users\cs2143\Desktop\Zhao\kernelclustering\unifiedres\coil.mat','res','alp','bta','gam')
% end
dlmwrite('/root/permanent/result/COIL5/gaussian50.txt',[alpha(i),beta(j),result],'-append','delimiter','\t','newline','pc');
%     end
    end
end
