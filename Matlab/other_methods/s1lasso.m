function[RMSP,RMSC,Y_full,Y_s1lasso,real_func,alpha_s1lasso] = s1lasso(N,I,M_pc,I0)
%[A_origin,A,b,Y_full,real_func,Y_pca] = datagenerate(N,I,M_pc,I0);
[A_origin,A,b,Y_full,real_func,Y_pca] = datagenerate_s2(N,I,M_pc,I0);
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
base1_alpha = lasso(base1_X,base1_b);
toc
base1_alpha_f = reshape(base1_alpha(:,50),M,I)';
rng(427)
[A_origin,A,b,Y_full,real_func,Y_pca] = datagenerate_s2(N,I,M_pc,I0);
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
base1_y = X*base1_alpha_f;
Y_s1lasso = base1_y;
alpha_s1lasso = base1_alpha_f;
RMSP = (norm(base1_y - Y_full,'fro')/(size(Y_full,1)*size(Y_full,2)))^0.5
%RMSC = (norm(base1_alpha_f-real_func,'fro')/(size(real_func,1)*size(real_func,2)))^0.5 
temp = [];
for i = 1:3
    RMSC = (norm(base1_alpha_f-real_func(:,((i-1)*50+1):i*50),'fro')/(size(real_func,1)*size(real_func,2)))^0.5; 
    temp = [temp,RMSC];
end
RMSC = mean(temp)