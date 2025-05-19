function[A_origin,A,Y,Y_full,real_func,Y_pca] = datagenerate(N,I,M_pc,I0)


%generalization of functional data
%rng(1010)

%%%%%%%%%%%%%%%%%
M=50;
rho = 0.75;
nbasis = min([50 M]);
T_domain = (0:(M-1))/(M-1); 

% Next we set up the covariates parameters for Matern process%

mu_X = zeros(I,1); 
Sig_X = ones(I,I); 
for i = 1:(I-1)
    for j = (i+1):I
        Sig_X(i,j) = rho^(j-i);
        Sig_X(j,i) = Sig_X(i,j);
    end
end 
% Here we set  signals parameters for our model %
nu_alpha = 2.5; 
range = 1/4;
variance = 1;
hyp = [log(range),log(variance)/2]; 
% we set errors parameters for our model%
nu_eps = 1.5;
mu_eps = zeros(M,1);
range1 = 0.25;  % the larger range1 is, Y_full and Y_full_noeps will be more similar
variance = 1;
hyp1 = [log(range1),log(variance)/2]; 
Sig_eps=covMaterniso(2*nu_eps,hyp1,T_domain');
Sig_eps2 = (Sig_eps + Sig_eps.')./2;
mu_alpha = zeros(M,1);
Sig_alpha = covMaterniso(2*nu_alpha,hyp,T_domain');

X = mvnrnd(mu_X,Sig_X,N);
X= zscore(X);
%rng(428)
alpha_1 = mvnrnd(mu_alpha,Sig_alpha,I); %true coefficient function in matrix form
alpha_1 = zscore(alpha_1);
eps = mvnrnd(mu_eps,Sig_eps2,N);
Y_full = X(:,1:I0)*alpha_1(1:I0,:)+eps; %X(:,I_X)*alpha_1+eps
%Y_full = X(:,:)*alpha_1+eps;
%Y_full = zscore(Y_full);
%Y_full_noeps = X(:,I_X)*alpha_1;
%subplot(1,2,1),plot(Y_full');
%subplot(1,2,2),plot(Y_full_noeps');

%try to add a trend to alpha
%alpha_trend=linspace(0,10,50);
%alpha_trend = repmat(alpha_trend,20,1)+alpha_1;
%alpha_1 = zscore(alpha_trend')';
%plot(alpha_1')

% If Y(t) is not a functional object, then the following steps are to do
% the FD conversion
Y_obs = zeros(N,M); 
T_obs = zeros(N,M); 
T_pos = zeros(N,M); 
Y_obs_cell = cell(1,N);
T_obs_cell = cell(1,N);
    for i = 1:N
        T_pos(i,:) = 1:M;
        T_obs(i,:) = T_domain(T_pos(i,:));
        Y_obs(i,:) = Y_full(i,T_pos(i,:)); 
        Y_obs_cell{i} = Y_obs(i,:);
        T_obs_cell{i} = T_obs(i,:);
    end
    Y_obs_vec = reshape(Y_obs',N*M,1);

bspline_basis = create_bspline_basis([0 1], nbasis, 4);
Y_f = data2fd(T_domain',Y_obs', bspline_basis,2,0.00001); %the fourth parameter 
% of data2fd is 0.00001 in the origin code, we keep it, when lambda equals
% to 0, the plot of functional version Y_f have a large deviation in edge
eval_pts = (0:(M-1))/(M-1);
Y_f_eval = eval_fd(eval_pts,Y_f)';
Y_f_eps = Y_full - Y_f_eval; %compare the difference between origin discrete data 
% Y_full and the functional version converted by data2fd
norm(Y_f_eps,'fro') %the error basically equals to 0 when lambda is set to 0
% subplot(1,2,1),plot(Y_f_eval')
% subplot(1,2,2),plot(Y_full')

%use fpca in ocean the larger M_pc is, the better fpca is, maybe we need to
%find a balance between the number of pc and computation power
Y_pca = fpca(Y_f);
Y_pc = proj(Y_f,Y_pca);
Y_pc = Y_pc(:,1:M_pc);
Y_f_reco = reco(Y_pca,Y_pc,M_pc); %functional object recovered by fpca
%in the reco function, fd is basically the same as data2fd
Y_pca_eval = eval_fd(eval_pts,Y_f_reco)'; %matrix object evaulated by the  
% functional object recovered by fpca
feps = Y_full-Y_pca_eval; %pca result vs origin data
eeps = Y_f_eval-Y_pca_eval; %pca result vs functional evaulated data
norm(feps,'fro');
norm(eeps,'fro');
% subplot(1,2,1),plot(Y_pca_eval')
% subplot(1,2,2),plot(Y_full')





%convert by FPCA, the error caused by fpca in this method is too much huge,
%the reason is not clear yet
% [N,I]=size(X); Y_pca = pca_fd(Y_f,M_pc); Y_scores = Y_pca.harmscr;
% FPCA_func_tmp = Y_pca.harmfd; eval_pts = (0:(M-1))/(M-1); 
% FPCA_func = eval_fd(FPCA_func_tmp,eval_pts); Y_eval = Y_scores*FPCA_func'+eval_fd(Y_pca.meanfd,eval_pts)'; 
% test = Y_full - Y_eval; Y_sm_vec = reshape(Y_scores',N*M_pc,1);


%convert coefficient function to FD conversion, and convert by FPCA to get
%ture version of x
T_obs = zeros(N,M); 
T_pos = zeros(N,M); 
alpha_obs_cell = cell(1,I0);
T_obs_cell = cell(1,I0);
    for i = 1:I
        T_pos(i,:) = 1:M;
        T_obs(i,:) = T_domain(T_pos(i,:));
        alpha_obs(i,:) = alpha_1(i,T_pos(i,:)); 
        alpha_obs_cell{i} = alpha_obs(i,:);
        T_obs_cell{i} = T_obs(i,:);
    end
    alpha_obs_vec = reshape(alpha_obs',I*M,1);

bspline_basis = create_bspline_basis([0 1], nbasis, 4);
alpha_f = data2fd(T_domain',alpha_obs', bspline_basis,2,0.00001);
alpha_f_eval = eval_fd(eval_pts,alpha_f)';
alpha_f_eps = alpha_1 - alpha_f_eval; %compare the difference between origin coefficient
% alpha_1 and the functional version converted by data2fd
norm(alpha_f_eps,'fro') %the error basically equals to 0 when lambda is set to 0
% subplot(1,2,1),plot(alpha_f_eval')
% subplot(1,2,2),plot(alpha_1')

%use fpca in ocean to convert alpha, just like what we do for y
pca_alpha = fpca(alpha_f);
pc_alpha = proj(alpha_f,pca_alpha);
alpha_f_reco = reco(pca_alpha,pc_alpha,M_pc); %functional object recovered by fpca
alpha_pca_eval = eval_fd(eval_pts,alpha_f_reco)'; %matrix object evaulated by the  
% functional object recovered by fpca
feps_a = alpha_1-alpha_pca_eval; %pca result vs origin data
eeps_a = alpha_f_eval-alpha_pca_eval; %pca result vs functional evaulated data
norm(feps_a,'fro')
norm(eeps_a,'fro')
% subplot(1,2,1),plot(alpha_pca_eval')
% subplot(1,2,2),plot(alpha_1')



A_origin = X;
A = kron(X,diag(ones(M_pc,1)));
Y = Y_pc;
Y_full = Y_full;
temp = zeros(I-I0,50);
real_func = [alpha_1(1:I0,:);temp];
Y_pca = Y_pca;
% real_func_fpca = pc_alpha(:,1:M_pc)';

end