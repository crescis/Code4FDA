function[RMSP,RMSC,Y_full,Y_s1,real_func,alpha_s1] = s1(N,I,M_pc,I0)
N = 10;I=20; I0 = 20; M_pc = 5;rng(1022);
[A_origin,A,Y,Y_full,real_func,Y_pca] = datagenerate(N,I,M_pc,I0);
%parameter setting
n = N; p = I; K = M_pc;
alpha1 = 0.1; alpha2 = 1;
Z = [];
for l = 1:n
    Z = blkdiag(Z,A_origin(l,:));
end
Q_origin = buildQ(A_origin);
Q_origin = sparse(Q_origin);
[m,~] = size(Q_origin); %number of edge

%parameter setting
n = N; p = I; K = M_pc;
%algorithm parameter setting
sppara.tau = 1.618;
sppara.sigma = 0.3;
sppara.tol = 10^-3;
sppara.kmax = 5000;
sppara.k = 0;
sppara.alpha1 = alpha1;
sppara.alpha2 = alpha2;
%spADMM
tic
[B,iter] = spADMM(Q_origin,Z,Y,sppara,n,p,m,K);
toc
Y_pre = Z*B*Y_pca.vectors(1:50,1:M_pc)';
subplot(1,2,1); plot(Y_full'); subplot(1,2,2); plot(Y_pre')
plot(Y_full')
axis([0,50,-12,10]);
plot(Y_pre')
axis([0,50,-12,10]);
norm(Y_full-Y_pre,'fro')/(size(Y_full,1)*size(Y_full,2))

