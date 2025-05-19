%%%Setting 1: low dimension
%1. algpara.lambda2 = 10; i.e. for every obdervation, coefficient is the
%same, use ADMMLP_2
%n = 15,50,100,200,500 p=n/5
N = 15; I = 20;I0 = I;M_pc = 5;rng(428)
tic
[RMSP,RMSC,Y_full,Y_s1,real_func,alpha_s1] = s1(N,I,M_pc,I0);
toc
% subplot(1,2,1); plot(Y_full'); subplot(1,2,2); plot(Y_s1')
% subplot(1,2,1); plot(real_func'); subplot(1,2,2); plot(alpha_s1')
% plot(Y_full(1:4,:)')
%2. use cv to choose the best value of lambda2 (RMSP as the criteria)
N = 500; I = 20;I0 = I; M_pc = 5;rng(428)
lambda_series = linspace(0,0.1,11)
tic
[RMSP,RMSC,Y_full,Y_s1,real_func,alpha_s1,RMSP_series,RMSC_series] = s1cv(N,I,M_pc,I0,lambda_series);
toc
RMSP_series
RMSC_series
subplot(1,2,1); plot(Y_full'); subplot(1,2,2); plot(Y_s1')
subplot(1,2,1); plot(real_func'); subplot(1,2,2); plot(alpha_s1')

%%345 is common among all settings, just remember to change the
%%datagenerate function and RMSC
%3. use ls to calculate the value on time points without fpca and without Q
N = 500; I = 20;I0 = I; M_pc = 5;rng(428)
tic
[RMSP,RMSC,Y_full,Y_s1lr,real_func,alpha_s1lr] = s1lr(N,I,M_pc,I0);
toc

%4. fsl use ADMMLP_2
N = 500; I = 20;I0 = I; M_pc = 5;rng(428)
tic
[RMSP,RMSC,Y_full,Y_s1lr,real_func,alpha_s1lr] = s1fsl(N,I,M_pc,I0);
toc

%5. use lasso to calculate the value on time points without fpca and without Q
N = 500; I = 20;I0 = I; M_pc = 5;rng(428)
tic
[RMSP,RMSC,Y_full,Y_s1lasso,real_func,alpha_s1lasso] = s1lasso(N,I,M_pc,I0);
toc
%%%Setting 2: cluster data - 3 clusters
%1. algpara.lambda2 = 10; i.e. for every obdervation, coefficient is the
%same, use ADMMLP_2
%n = 15,50,100,200,500 p=n/5
N = 375; I = 15;I0 = I; M_pc = 5;
rng(428)
tic
[RMSP,RMSC,Y_full,Y_s2,real_func,alpha_s2] = s2(N,I,M_pc,I0);
toc
subplot(1,2,1); plot(Y_full'); subplot(1,2,2); plot(Y_s2')
subplot(1,2,1),plot(alpha_s2');subplot(1,2,2),plot(real_func')
a = reshape(alpha_s2',50,450)';
%2.  use cv to choose the best value of lambda2 (RMSP as the criteria)
N = 375; I = 15;I0 = I;M_pc = 5;rng(428);
lambda_series = linspace(0,0.1,11);
tic
[RMSP,RMSC,Y_full,Y_s2,real_func,alpha_s2,RMSP_series,RMSC_series] = s2cv(N,I,M_pc,I0,0,0);   
toc
RMSP_series
RMSC_series

a = reshape(alpha_s2',50,1800)';
plot(a')
%%3,4,5 don't consider group structure, just the same as setting 1,change
%%the data generation function and RMSC
%3. use ls to calculate the value on time points without fpca and without Q
N = 375; I = 15;I0 = I; M_pc = 5;rng(428)
tic
[RMSP,RMSC,Y_full,Y_s1lr,real_func,alpha_s1lr] = s1lr(N,I,M_pc,I0);
toc

%4. fsl use ADMMLP_2
N = 375; I = 15;I0 = I; M_pc = 5;rng(428)
tic
[RMSP,RMSC,Y_full,Y_s1lr,real_func,alpha_s1lr] = s1fsl(N,I,M_pc,I0);
toc
%5. use lasso to calculate the value on time points without fpca and without Q
N = 30; I = 15;I0 = I; M_pc = 5;rng(428)
tic
[RMSP,RMSC,Y_full,Y_s1lasso,real_func,alpha_s1lasso] = s1lasso(N,I,M_pc,I0);
toc

%%%Setting 3: high dimensional
%1. algpara.lambda2 = 10; i.e. for every observation, coefficient is the
%same, use ADMMLP_2
%n = 15,30,60,90,120 p=2*n
N = 15; I = 120;I0 = 15;M_pc = 5;rng(428)
tic
[RMSP,RMSC,Y_full,Y_s1,real_func,alpha_s1] = s1(N,I,M_pc,I0);
toc
real_func_S3s1 = [real_func_S3s1,real_func]; 
alpha_S3s1 = [alpha_S3s1,alpha_s1];
subplot(1,2,1); plot(Y_full'); subplot(1,2,2); plot(Y_s1')
subplot(1,2,1); plot(real_func'); subplot(1,2,2); plot(alpha_s1')

%2. use cv to choose the best value of lambda2 (RMSP as the criteria)
N = 15; I = 120;I0 = 15; M_pc = 5;rng(428)
lambda_series = 0.01;
tic
[RMSP,RMSC,Y_full,Y_s1,real_func,alpha_s1,RMSP_series,RMSC_series] = s2cv(N,I,M_pc,I0,lambda_series,0.01);
toc
real_func_S3s2 = [real_func_S3s2,real_func]; 
alpha_S3s2 = [alpha_S3s2,alpha_s1];
RMSP_series
RMSC_series
RMSP_matrix = [RMSP_matrix;RMSP_series];
RMSC_matrix = [RMSC_matrix;RMSC_series];
subplot(1,2,1); plot(Y_full'); subplot(1,2,2); plot(Y_s1')
subplot(1,2,1); plot(real_func'); subplot(1,2,2); plot(alpha_s1')
a = reshape(alpha_s1',50,1800)';
plot(a')
for i = 1:size(a,1)
    temp = norm(a(i,:));
    if temp < 0.01
        a(i,:) = 0;
    end 
end

sum(sum(alpha_S3s2(1:15,:)<0.00000000001))
15*15750 - 131243
sum(sum(alpha_S3s2(16:120,:)<0.00000000001))
105*15750-848983

num = [0,15,30,60,90,120];
for i = 1:5
    TP2(i) = (15*50*num(i+1)-sum(sum(alpha_S3s2(1:15,(50*num(i)+1):50*num(i+1))<0.1)))/(15*50*num(i+1));
    FP2(i) = (105*50*num(i+1)-sum(sum(alpha_S3s2(16:120,(50*num(i)+1):50*num(i+1))<0.01)))/(105*50*num(i+1));
end
TP2
FP2


%3. use ls to calculate the value on time points without fpca and without Q
N = 120; I = 120;I0 = 15; M_pc = 5;rng(428)
tic
[RMSP,RMSC,Y_full,Y_s1lr,real_func,alpha_s1lr] = s1lr(N,I,M_pc,I0);
toc
real_func_S3s4 = [real_func_S3s4,real_func]; 
alpha_S3s4 = [alpha_S3s4,alpha_s1lr];



%4
N = 120; I = 120;I0 = 15; M_pc = 5;rng(428)
tic
[RMSP,RMSC,Y_full,Y_s1lr,real_func,alpha_s1lr] = s1fsl(N,I,M_pc,I0);
toc
real_func_S3s3 = [real_func_S3s3,real_func]; 
alpha_S3s3 = [alpha_S3s3,alpha_s1lr];
TP4 = [];FP4 = [];
for i = 1:5
    TP4(i) = (15*50-sum(sum(alpha_S3s4(1:15,(50*(i-1)+1):50*i)<0.1)))/(15*50);
    FP4(i) = (105*50-sum(sum(alpha_S3s4(16:120,(50*(i-1)+1):50*i)<0.1)))/(105*50);
end
TP4
FP4
%5. use lasso to calculate the value on time points without fpca and without Q
N = 120; I = 120;I0 = 15; M_pc = 5;rng(428)
tic
[RMSP,RMSC,Y_full,Y_s1lasso,real_func,alpha_s1lasso] = s1lasso(N,I,M_pc,I0);
toc
real_func_S3s5 = [real_func_S3s5,real_func]; 
alpha_S3s5 = [alpha_S3s5,alpha_s1lasso];
TP5 = [];FP5 = [];
for i = 1:5
    TP5(i) = (15*50-sum(sum(alpha_S3s5(1:15,(50*(i-1)+1):50*i)<0.1)))/(15*50);
    FP5(i) = (105*50-sum(sum(alpha_S3s5(16:120,(50*(i-1)+1):50*i)<0.1)))/(105*50);
end
TP5
FP5

nnz(alpha_S3s5(1:15,1:250))
15*250-nnz(alpha_S3s5(1:15,1:250))
nnz(alpha_S3s5(16:120,1:250))
105*250-nnz(alpha_S3s5(16:120,1:250))
%%%Setting 4: high dimensional & network
%1. algpara.lambda2 = 10; i.e. for every observation, coefficient is the
%same, use ADMMLP_2
%n = 15,30,60,90,120 p=2*n
N = 120; I = 120;I0 = 15;M_pc = 5;rng(428)
tic
[RMSP,RMSC,Y_full,Y_s1,real_func,alpha_s1] = s2(N,I,M_pc,I0);
toc
real_func_S4s1 = [real_func_S4s1,real_func]; 
alpha_S4s1 = [alpha_S4s1,alpha_s1];
subplot(1,2,1); plot(Y_full'); subplot(1,2,2); plot(Y_s1')
subplot(1,2,1); plot(real_func'); subplot(1,2,2); plot(alpha_s1')

%2. use cv to choose the best value of lambda2 (RMSP as the criteria)
N = 120; I = 120; I0 = 15; M_pc = 5;rng(428)
lambda_series = 0.1;
tic
[RMSP,RMSC,Y_full,Y_s1,real_func,alpha_s1,RMSP_series,RMSC_series] = s2cv(N,I,M_pc,I0,lambda_series,1);
toc
real_func_S4s2 = [real_func_S4s2,real_func]; 
alpha_S4s2 = [alpha_S4s2,alpha_s1];
RMSP_series
RMSC_series
RMSP_matrix = [RMSP_matrix;RMSP_series];
RMSC_matrix = [RMSC_matrix;RMSC_series];
subplot(1,2,1); plot(Y_full'); subplot(1,2,2); plot(Y_s1')
subplot(1,2,1); plot(real_func'); subplot(1,2,2); plot(alpha_s1')
a = reshape(alpha_s1',50,1800)';
plot(a')
for i = 1:size(a,1)
    temp = norm(a(i,:));
    if temp < 0.01
        a(i,:) = 0;
    end 
end

num = [0,15,30,60,90,120];
for i = 1:5
    TP2(i) = (15*50*num(i+1)-sum(sum(alpha_S4s2(1:15,(50*num(i)+1):50*num(i+1))<0.001)))/(15*50*num(i+1));
    FP2(i) = (105*50*num(i+1)-sum(sum(alpha_S4s2(16:120,(50*num(i)+1):50*num(i+1))<0.01)))/(105*50*num(i+1));
end
TP2
FP2

%3. use ls to calculate the value on time points without fpca and without Q
N = 120; I = 120;I0 = 15; M_pc = 5;rng(428)
tic
[RMSP,RMSC,Y_full,Y_s1lr,real_func,alpha_s1lr] = s1lr(N,I,M_pc,I0);
toc
real_func_S4s3 = [real_func_S4s3,real_func]; 
alpha_S4s3 = [alpha_S4s3,alpha_s1lr];

%4
N = 30; I = 120;I0 = 15; M_pc = 5;rng(428)
tic
[RMSP,RMSC,Y_full,Y_s1lr,real_func,alpha_s1lr] = s1fsl(N,I,M_pc,I0);
toc
real_func_S4s4 = [real_func_S4s4,real_func]; 
alpha_S4s4 = [alpha_S4s4,alpha_s1lr];

%5. use lasso to calculate the value on time points without fpca and without Q
% remember to change the datageneration_s2 twice
N = 120; I = 120;I0 = 15; M_pc = 5;rng(428)
tic
[RMSP,RMSC,Y_full,Y_s1lasso,real_func,alpha_s1lasso] = s1lasso(N,I,M_pc,I0);
toc
real_func_S4s5 = [real_func_S4s5,real_func]; 
alpha_S4s5 = [alpha_S4s5,alpha_s1lasso];