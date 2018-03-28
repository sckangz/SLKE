load YALE_165n_1024d_15c_uni;
load('data/YALE_165n_1024d_15c_zscore_uni_kernel/YALE_165n_1024d_15c_zscore_uni_kernel_gaussian_0.01_post_Sample-Scale.mat');
D=K;
%load('data/YALE_165n_1024d_15c_zscore_uni_kernel/YALE_165n_1024d_15c_zscore_uni_kernel_polyplus_2_post_Sample-Scale.mat')
para1=[1e-5 1e-3 .10 10 100 200 1000];
para2=[.001 .1  10 50 100 200 300 400];
for ij=1:length(para1)
gamma=para1(ij);
for iji=1:length(para2)
mu=para2(iji);
s=y;
[result]= SLKEs( K,D,y,gamma,mu)
dlmwrite('/Users/apple/Desktop/Yale.txt',[lambda,mu,result],'-append','delimiter','\t','newline','pc');
    
end
end
