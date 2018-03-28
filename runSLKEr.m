load tr11_414n_6429d_9c_tfidf_uni;
load('data/tr11_414n_6429d_9c_tfidf_uni_kernel/tr11_414n_6429d_9c_tfidf_uni_kernel_polynomial_2_post_Sample-Scale.mat');
D=K;
para1=[1e-5 1e-3 .10 10 100 200 1000];
para2=[.001 .1  10 50 100 200 300 400];

for ij=1:length(para1)
gamma=para1(ij);
for iji=1:length(para2)
mu=para2(iji);
s=y;

fprintf('params%12.6f%12.6f\n',gamma,mu)
[result]= SLKEs( K,D,y,gamma,mu)



end
end
