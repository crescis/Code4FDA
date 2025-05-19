function[A_final,Q,Q_new,d,algpara,T_inv,x0,y0,z0] = datagenerate_2(A_origin,M_pc,I,N)
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
algpara.lambda = 1;
algpara.lambda1 = 0;
algpara.lambda2 = 3;
[algpara.N_alg,algpara.d_alg] = size(Q_new);
%initial value
T = (1/algpara.alpha)*(A_final'*A_final) + speye(d); %sparse form of matrix helps a lot
T_inv = inv(T);
x0 = rand(algpara.d_alg,1);
z0 = rand(algpara.N_alg,1);
y0 = rand(algpara.N_alg,1);
end