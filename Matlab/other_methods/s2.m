function[RMSP,RMSC,Y_full,Y_s2,real_func,alpha_s2] = s2(N,I,M_pc,I0)
%[A_origin,A,b,Y_full,real_func,Y_pca] = datagenerate(N,I,M_pc,I0);
[A_origin,A,b,Y_full,real_func,Y_pca] = datagenerate_s2(N,I,M_pc,I0);
[A_final,Q,Q_new,d,algpara,T_inv,x0,y0,z0] = datagenerate_2(A_origin,M_pc,I,N);
A_final = blockdiag(A_origin,M_pc);
A_final = sparse(A_final);
Q = buildQ(A_origin);
Q = sparse(Q);
[m,~] = size(Q); %number of edge
%[Q_new] = buildQ_new(Q,m,I*M_pc);
Q_new = kron(Q,speye(M_pc*I));
%parameter setting
d = N*M_pc*I; %dimension of x
eigval = max(eig(Q'*Q));
%algorithm parameter setting
algpara.x_tol = 10^(-3);
algpara.z_tol = 10^(-3);
algpara.rho = 0.5;
algpara.alpha = algpara.rho*eigval;
algpara.beta = algpara.rho;
algpara.iter = 10000;
algpara.lambda = 10;
algpara.lambda1 = 0; %control lasso term
algpara.lambda2 = 10;                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0; %control sum-of-norms regularizer term
[algpara.N_alg,algpara.d_alg] = size(Q_new);
%initial value
T = (1/algpara.alpha)*(A_final'*A_final) + speye(d); %sparse form of matrix helps a lot
T_inv = inv(T);
x0 = rand(algpara.d_alg,1);
z0 = rand(algpara.N_alg,1);
y0 = rand(algpara.N_alg,1);
%ADLPMM
[x,y,z,record] = ADLPMM_2(x0,y0,z0,Q_new,A_final,b,T_inv,algpara);
rng(427)
[A_origin,A,b,Y_full,real_func,Y_pca] = datagenerate_s2(N,I,M_pc,I0);
[A_final,Q,Q_new,d,algpara,T_inv,x0,y0,z0] = datagenerate_2(A_origin,M_pc,I,N);
Y_s2 = reshape(A_final*x,M_pc,N)'*Y_pca.vectors(1:50,1:M_pc)'; %every column of pca.vectors represent one principle function's value on time points
result_s1 = reshape(x,M_pc*I,N)';
temp = [];
alpha_s2 = [];
    for j = (1:(N/3))
        pc_s1 = reshape(result_s1(j,:),M_pc,I)';
        alpha_s1 = pc_s1*Y_pca.vectors(1:50,1:M_pc)';
        alpha_s2 = [alpha_s2,alpha_s1];
        temp(j) = (norm(alpha_s1-real_func(:,1:50),'fro')/(size(real_func,1)*size(real_func,2)))^0.5;
    end
    for j = ((N/3+1):(N/3)*2)
        pc_s1 = reshape(result_s1(j,:),M_pc,I)';
        alpha_s1 = pc_s1*Y_pca.vectors(1:50,1:M_pc)';
        alpha_s2 = [alpha_s2,alpha_s1];
        temp(j) = (norm(alpha_s1-real_func(:,51:100),'fro')/(size(real_func,1)*size(real_func,2)))^0.5;
    end
    for j = (((N/3)*2+1):N)
        pc_s1 = reshape(result_s1(j,:),M_pc,I)';
        alpha_s1 = pc_s1*Y_pca.vectors(1:50,1:M_pc)';
        alpha_s2 = [alpha_s2,alpha_s1];
        temp(j) = (norm(alpha_s1-real_func(:,101:150),'fro')/(size(real_func,1)*size(real_func,2)))^0.5;
    end
real_func = reshape(real_func',50,3*I)';
RMSP = (norm(Y_s2-Y_full,'fro')/(size(Y_full,1)*size(Y_full,2)))^0.5 %root-mean-squared prediction error
RMSC = mean(temp)