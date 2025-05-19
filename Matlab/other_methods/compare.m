%model 0: main model-calculate the value on pc_alpha with fpca and with Q
%via ADLPMM
base0_result = reshape(x,M_pc*I,N)';
%after run the main 'test.m' file
base0_pc = reshape(base0_result(1,:),M_pc,I)';
base0_alpha = base0_pc*Y_pca.vectors(1:50,1:M_pc)';
subplot(1,2,1),plot(base0_alpha')
subplot(1,2,2),plot(real_func')
norm((base0_alpha-real_func),'fro')
base0_y_pc = reshape(A_final*x,M_pc,N)';
base0_y = base0_y_pc*Y_pca.vectors(1:50,1:M_pc)'; %every column of pca.vectors represent one principle function's value on time points
subplot(1,2,1),plot(base0_y')
subplot(1,2,2),plot(Y_full')
norm(base0_y-Y_full,'fro')/(size(Y_full,1)*size(Y_full,2))




%model 1: calculate the value on time points without fpca and without Q
%using lsqr
base1_X = [];
X = A_origin;
M = 50;
for i = 1:size(X,1)
    temp1 = [];
    for j = 1:M
        temp1 = blkdiag(temp1,X(i,:));
    end
    base1_X((50*(i-1)+1):(50*i),:) = temp1;
end
base1_b = reshape(Y_full',N*M,1);
tic
base1_alpha = lsqr(base1_X,base1_b);
toc
base1_alpha_f = reshape(base1_alpha,M,I)';
subplot(1,2,1),plot(real_func')
subplot(1,2,2),plot(base1_alpha_f')
base1_y = X*base1_alpha_f;
subplot(1,2,1),plot(base1_y')
subplot(1,2,2),plot(Y_full')
base1_eps = base1_y - Y_full;
norm(base1_eps,'fro')/(size(Y_full,1)*size(Y_full,2))


%model 2: calculate the value on time points without fpca and without Q
%using lasso via ADLPMM
base2_I = eye(size(base1_X,2));
[algpara.N_alg,algpara.d_alg] = size(base2_I);
%initial value
T = (1/algpara.alpha)*(base1_X'*base1_X) + speye(size(base1_X,2)); %sparse form of matrix helps a lot
T_inv = inv(T);
x0 = rand(algpara.d_alg,1);
z0 = rand(algpara.N_alg,1);
y0 = rand(algpara.N_alg,1);
algpara.lambda = 10;
[x,y,z,record] = ADLPMM(x0,y0,z0,base2_I,base1_X,base1_b,T_inv,algpara);
base3_alpha_pc = x;
base4_alpha_f = reshape(base3_alpha_pc,50,20)';
subplot(1,2,1),plot(real_func')
subplot(1,2,2),plot(base4_alpha_f')
base4_y = X*base4_alpha_f;
subplot(1,2,1),plot(base4_y')
subplot(1,2,2),plot(Y_full')
base2_eps = base4_y - Y_full;
norm(base2_eps,'fro')


%model3: calculate the value on pc_alpha with fpca and without Q
%using lsqr
base3_alpha_pc_origin = lsqr(A,b);
base3_alpha_pc = reshape(base3_alpha_pc_origin,M_pc,I)';
base3_alpha = base3_alpha_pc*Y_pca.vectors(1:50,1:M_pc)';
subplot(1,2,1),plot(base3_alpha_pc')
subplot(1,2,2),plot(alpha_1')
norm(base3_alpha-alpha_1,'fro')
base3_y = A_origin*base3_alpha_pc*Y_pca.vectors(1:50,1:M_pc)';
subplot(1,2,1),plot(base3_y')
subplot(1,2,2),plot(Y_full')
norm(base3_y-Y_full,'fro')

%model4: calculate the value on pc_alpha with fpca and without Q
%using lasso via ADLPMM
base4_I = eye(size(A,2));
[algpara.N_alg,algpara.d_alg] = size(base4_I);
%initial value
T = (1/algpara.alpha)*(A'*A) + speye(size(A,2)); %sparse form of matrix helps a lot
T_inv = inv(T);
x0 = rand(algpara.d_alg,1);
z0 = rand(algpara.N_alg,1);
y0 = rand(algpara.N_alg,1);
algpara.lambda = 3;
[x,y,z,record] = ADLPMM(x0,y0,z0,base4_I,A,b,T_inv,algpara);
base3_alpha_pc = reshape(x,M_pc,I)';
base4_alpha = base3_alpha_pc*Y_pca.vectors(1:50,1:M_pc)';
subplot(1,2,1),plot(real_func')
subplot(1,2,2),plot(base4_alpha')
norm(base4_alpha-alpha_1,'fro')
base4_y_pc = reshape(A*x,M_pc,N)';
base4_y = base4_y_pc*Y_pca.vectors(1:50,1:M_pc)';
subplot(1,2,1),plot(base4_y')
subplot(1,2,2),plot(Y_full')
norm(base4_y - Y_full,'fro')






%baseline 2: comparision with FSL
% set parameters for FSL
nlam_base=100;
BIC_para = 1;
lamratio=0.001;
% Call on FSL and AFSL
[history]=AFSL(Y_f, X,M, M_pc, nlam_base, BIC_para, lamratio)
subplot(1,2,1),plot(real_func')
subplot(1,2,2),plot(history.alpha_hat_origin_FSL')
base1_eps = X*history.alpha_hat_origin_FSL - Y_full;
norm(base1_eps,'fro')