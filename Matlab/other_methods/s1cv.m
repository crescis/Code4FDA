function[RMSP,RMSC,Y_full,Y_s1,real_func,alpha_s1,RMSP_series,RMSC_series] = s1cv(N,I,M_pc,I0,lambda_series,lambda1)
[A_origin,A,b,Y_full,real_func,Y_pca] = datagenerate(N,I,M_pc,I0);
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
algpara.lambda1 = lambda1; %control lasso term
RMSP_series = [];
RMSC_series = [];
k = numel(lambda_series);
alpha_s = [];
for i = (1:k)
    algpara.lambda2 = lambda_series(i); %control sum-of-norms regularizer term
    %ADLPMM
    [x,y,z,record] = ADLPMM_2(x0,y0,z0,Q_new,A_final,b,T_inv,algpara);
    Y_s1 = reshape(A_final*x,M_pc,N)'*Y_pca.vectors(1:50,1:M_pc)'; %every column of pca.vectors represent one principle function's value on time points
    result_s1 = reshape(x,M_pc*I,N)';
    temp = [];
    for j = (1:N)
        pc_s1 = reshape(result_s1(j,:),M_pc,I)';
        alpha_s1 = pc_s1*Y_pca.vectors(1:50,1:M_pc)';
        temp(j) = (norm(alpha_s1-real_func,'fro')/(size(real_func,1)*size(real_func,2)))^0.5;
        alpha_s = [alpha_s,alpha_s1];
    end
    RMSP = (norm(Y_s1-Y_full,'fro')/(size(Y_full,1)*size(Y_full,2)))^0.5; %root-mean-squared prediction error
    RMSC = mean(temp);
    RMSP_series(i) = RMSP;
    RMSC_series(i) = RMSC;
    alpha_s1 = alpha_s;
end

[RMSP,b] = min(RMSP_series);
RMSP
RMSC = RMSC_series(b)

