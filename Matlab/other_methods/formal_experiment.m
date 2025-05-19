%--------------------------------------------------------------------------
%  Simulation Driver for the Manuscript Experiments
%
%  The concrete implementations of all algorithms referenced in the paper live
%  in the `./algorithm/` directory.
%--------------------------------------------------------------------------

clear all
N = 10;
I = 5;
I0 = I; 
M_pc = 5;
rng(428)
[A_origin,A,b,Y_full,real_func,Y_pca] = datagenerate(N,I,M_pc,I0);
%GPR toolbox is needed for data generation, and the plot function should be
%adjusted to the builtin function.
[A_final,Q,Q_new,d,algpara,T_inv,x0,y0,z0] = datagenerate_2(A_origin,M_pc,I,N);
algpara.lambda = 3;
tic 
[x,y,z,record] = ADLPMM(x0,y0,z0,Q_new,A_final,b,T_inv,algpara);
toc
base0_y = reshape(A_final*x,M_pc,N)'*Y_pca.vectors(1:50,1:M_pc)'; %every column of pca.vectors represent one principle function's value on time points
norm(base0_y-Y_full,'fro')/(size(Y_full,1)*size(Y_full,2))
algpara.lambda1 = 2;%l1 term
algpara.lambda2 = 1;%sum of norms term
tic
[x,y,z,record] = ADLPMM_2(x0,y0,z0,Q_new,A_final,b,T_inv,algpara);
toc
base0_y = reshape(A_final*x,M_pc,N)'*Y_pca.vectors(1:50,1:M_pc)'; %every column of pca.vectors represent one principle function's value on time points
norm(base0_y-Y_full,'fro')/(size(Y_full,1)*size(Y_full,2))



AFSL_pre = A_origin*history.alpha_hat_origin_FSL-Y_full;
norm(AFSL_pre(:,5:50),'fro')/(size(Y_full,1)*size(Y_full,2))