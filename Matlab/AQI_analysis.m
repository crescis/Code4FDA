%---------------------------------------------------------------------------
%  spADMM Analysis of Air‑Quality‑Index (AQI) Data
%
%  This script applies the semi-proximal ADMM (spADMM) algorithm to AQI data
%  under multiple configurations of the model parameters.  For each setting,
%  the script exports the key numerical results to *.csv files, making them
%  easy to import into R for subsequent statistical analysis and visualization.
%--------------------------------------------------------------------------- 
%surfate2021 = readtable("sulfate2021_24.csv");
surfate2021_3col = readtable("spec_sulfate_24_3col.csv",VariableNamingRule="preserve");
G = groupsummary(surfate2021_3col,'Latitude');
G = table2array(G);
nbasis = 50;
Y_f_299 = zeros(299,50);
for i = 1:size(G,1)
    M = G(i,2);
    T_domain = (0:(M-1))/(M-1); 
    bspline_basis = create_bspline_basis([0 1], nbasis, 5);
    if i == 1
        Y_obs = table2array(surfate2021_3col(1:G(i,2),4))';
    else
        Y_obs = table2array(surfate2021_3col((G(i-1,2)+1):(G(i-1,2)+G(i,2)),4))';
    end
    Y_f = data2fd(T_domain',Y_obs', bspline_basis,2,0.0000001);
    eval_pts = (0:(50-1))/(50-1);
    Y_f_eval = eval_fd(eval_pts,Y_f)';
    Y_f_299(i,:) = Y_f_eval;
end
%plot(Y_f_299')
%since the toolbox 'fadM' also defined a 'plot' function to plot basis
%function, the builtin function should be used as following, this rule can
%also applied later if an error occured.
builtin('plot',Y_f_299')
csvwrite('Y_f_299.csv',Y_f_299)
fulldata = readtable('fulldata.csv');
Y_full = table2array(fulldata(:,3:52));
X1 = table2array(fulldata(:,57:80));
X2 = table2array(fulldata(:,82:87));
X = [X1 X2];
geo = table2array(fulldata(:,53:54));
A_origin = X;
M_pc = 5;
A = kron(X,diag(ones(M_pc,1)));
N = size(X,1); M = 50; I = size(X,2);
nbasis = min([50 M]);
T_domain = (0:(M-1))/(M-1); 
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
Y_f = data2fd(T_domain',Y_obs', bspline_basis,2,0.00001);
Y_pca = fpca(Y_f);
Y_pc = proj(Y_f,Y_pca);
Y_pc = Y_pc(:,1:M_pc);
Y = Y_pc;
n = N; p = I; K = M_pc;
alpha1 = 0.1; alpha2 = 1;
Z = [];
for l = 1:n
    Z = blkdiag(Z,A_origin(l,:));
end
Q_origin = buildQ(geo);
Q_origin = sparse(Q_origin);
[m,~] = size(Q_origin);
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
%subplot(1,2,1); plot(Y_full'); subplot(1,2,2); plot(Y_pre')
norm(Y_full-Y_pre,'fro')/(size(Y_full,1)*size(Y_full,2))


aqidata = readtable("aqidata.csv");
G = groupsummary(aqidata,'fullname');
G = table2cell(G);
G = cell2table(G);
nbasis = 100;
Y_f_1002 = zeros(1002,100);
for i = 1:size(G,1)
    M = table2array(G(i,2));
    T_domain = (0:(M-1))/(M-1); 
    bspline_basis = create_bspline_basis([0 1], nbasis, 5);
    if i == 1
        Y_obs = table2array(aqidata(1:table2array(G(i,2)),7))';
    else
        Y_obs = table2array(aqidata((table2array(G(i-1,2))+1):(table2array(G(i-1,2))+table2array(G(i,2))),7))';
    end
    Y_f = data2fd(T_domain',Y_obs', bspline_basis,2,0.00001);
    eval_pts = (0:(100-1))/(100-1);
    Y_f_eval = eval_fd(eval_pts,Y_f)';
    Y_f_1002(i,:) = Y_f_eval;
end
%plot(Y_f_1002')
builtin('plot',Y_f_1002')
csvwrite('Y_f_1002.csv',Y_f_1002)


fulldata = readtable('fulldata_aqi.csv');
Y_full = table2array(fulldata(:,4:103));
X1 = table2array(fulldata(:,115:120));
X2 = table2array(fulldata(:,122:145));
X = [X1 X2];
%Y_full = zscore(Y_full);
%X = zscore(X);
geo = table2array(fulldata(:,110:111));
A_origin = X;
M_pc = 5;
A = kron(X,diag(ones(M_pc,1)));
N = size(X,1); M = 100; I = size(X,2);
nbasis = min([50 M]);
T_domain = (0:(M-1))/(M-1); 
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
Y_f = data2fd(T_domain',Y_obs', bspline_basis,2,0.00001);
Y_pca = fpca(Y_f);
Y_pc = proj(Y_f,Y_pca);
Y_pc = Y_pc(:,1:M_pc);
Y = Y_pc;
n = N; p = I; K = M_pc;
alpha1 = 0.1; alpha2 = 5;
Z = [];
for l = 1:n
    Z = blkdiag(Z,A_origin(l,:));
end
Z = sparse(Z);
Q_origin = buildQ(geo);
Q_origin = sparse(Q_origin);
[m,~] = size(Q_origin);
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
%subplot(1,2,1); plot(Y_full'); subplot(1,2,2); plot(Y_pre')

%turn Y_full from 100 time points to 50
N = 960; M = 100;
nbasis = min([50 M]);
T_domain = (0:(M-1))/(M-1); 
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
Y_f = data2fd(T_domain',Y_obs', bspline_basis,2,0.00001);
eval_pts = (0:(50-1))/(50-1);
Y_f_eval = eval_fd(eval_pts,Y_f)';
norm(Y_f_eval-Y_pre,'fro')/(size(Y_full,1)*size(Y_full,2))
for i = 1:size(Y_pre,1)
    Y_pre(i,:) = Y_pre(i,:) + mean(Y_f_eval,1);
end
%Y_pre = Y_pre + mean(mean(Y_f_eval));
%add the mean value of Y_f_eval
Y_error = Y_f_eval-Y_pre;
Y_error_norm = sum(abs(Y_error).^2,2).^(1/2);
norm(Y_error_norm)/(251*30)
%subplot(1,2,1); plot(Y_pre'); subplot(1,2,2); plot(Y_f_eval')
[Y_error_sort,Y_error_sort_index] = sort(Y_error_norm);
rng(125);
random_num = Y_error_sort_index(randperm(numel(Y_error_sort_index(1:600)),15));
Y_pre_10 = Y_pre(random_num,:);
Y_f_eval_10 = Y_f_eval(random_num,:);
%subplot(1,2,1); plot(Y_pre_10'); subplot(1,2,2); plot(Y_f_eval_10'); ylim([23 50])

% plot(Y_f_eval(Y_error_sort_index(1),:),'LineWidth',2);
% hold on;
% plot(Y_pre(Y_error_sort_index(1),:),'LineWidth',2);
% legend('Real','Prediction','Location','NorthEast', 'FontName','Times New Roman','FontSize',8,'FontWeight','normal');

B_coeff = B*Y_pca.vectors(1:50,1:M_pc)';
%deal with the coefficients
B_vec = reshape(B', size(B,1)*size(B,2),1);
B_clus = reshape(B_vec,M_pc*size(A_origin,2),size(A_origin,1))';
idx = kmeans(B_clus,4);
B1idx = 26:30:28800;
B1 = B_coeff(B1idx,:);
% plot(B1')
B1_rowsum = sum(B1,2);
[B1_rank,B1_rank_index] = sort(B1_rowsum);
random_num = randperm(960,4)
%random_num = [B1_rank_index(1),B1_rank_index(10),B1_rank_index(20),B1_rank_index(30)];
%random_num = B1_rank_index(1:4)
%random_num = [24,406,778,97]
random_num = [B1_rank_index(10),B1_rank_index(80),B1_rank_index(150),B1_rank_index(960)];
B1_sub = B1(random_num,:);
% plot(B1_sub')


% [B1_nonzero_rank,B1_nonzero_rank_index] = sort(sum(B1>0,2));
% random_num = [B1_nonzero_rank_index(930:960)];
% B1_sub = B1(random_num,:);
% plot(B1_sub')
N = 4; M = 45;
Y_obs = zeros(N,M); 
T_obs = zeros(N,M); 
T_pos = zeros(N,M); 
Y_obs_cell = cell(1,N);
T_obs_cell = cell(1,N);
    for i = 1:N
        T_pos(i,:) = 1:M;
        T_obs(i,:) = T_domain(T_pos(i,:));
        Y_obs(i,:) = B1_sub(i,T_pos(i,:)); 
        Y_obs_cell{i} = Y_obs(i,:);
        T_obs_cell{i} = T_obs(i,:);
    end
    Y_obs_vec = reshape(Y_obs',N*M,1);

bspline_basis = create_bspline_basis([0 1], nbasis, 4);
Y_f = data2fd(T_domain',Y_obs', bspline_basis,2,0.0001);
plot(Y_f)
legend('M1','M2','M3','M4','Location','NorthWest', 'FontName','Times New Roman','FontSize',8,'FontWeight','normal'); 


center = readtable("center.csv");
center = table2array(center);
plot(center')
% csvwrite('B1.csv',B1)
% csvwrite('idx.csv',idx)
% idx_2 = kmeans(Y_full(:,1),4);
% csvwrite('idx_2.csv',idx_2)


%choose 30 samples to plot
Y_f_1002 = zeros(30,100);
for i = 1:30
    M = table2array(G(i,2));
    T_domain = (0:(M-1))/(M-1); 
    bspline_basis = create_bspline_basis([0 1], nbasis, 5);
    if i == 1
        Y_obs = table2array(aqidata(1:table2array(G(i,2)),7))';
    else
        Y_obs = table2array(aqidata((table2array(G(i-1,2))+1):(table2array(G(i-1,2))+table2array(G(i,2))),7))';
    end
    Y_f = data2fd(T_domain',Y_obs', bspline_basis,2,0.00001);
    eval_pts = (0:(100-1))/(100-1);
    Y_f_eval = eval_fd(eval_pts,Y_f)';
    Y_f_1002(i,:) = Y_f_eval;
end
% plot(Y_f_1002')
% axis( [0,100,15,60])  
% set(gca,'XTick',[0:100/11:100]) 
% xlabel('Time domain', 'FontName','Times New Roman','FontSize',11,'FontWeight','normal')  
% xticklabels({'2021-1','2021-2','2021-3','2021-4','2021-5','2021-6','2021-7','2021-8','2021-9','2011-10','2021-11','2021-12'})
% ylabel('AQI Index', 'FontName','Times New Roman','FontSize',11,'FontWeight','normal')