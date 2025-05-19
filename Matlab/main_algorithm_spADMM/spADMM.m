function[B,iter] = spADMM(Q_origin,Z,Y,sppara,n,p,m,K)
tau = sppara.tau; sigma = sppara.sigma; tol = sppara.tol;
kmax = sppara.kmax; k = sppara.k;
alpha1 = sppara.alpha1; alpha2 = sppara.alpha2;
% Q1 = [];
% for l = 1:(n*p-1)
%     Q1 = [sparse(Q1);l*ones(n*p-l,1)];
% end
% Q1 = Q1';
% Q2 = [];
% for l = 2:n*p
%     Q2 = [sparse(Q2) l:n*p];
% end
% Q = sparse(n*p*(n*p-1)/2,n*p);
% for l = 1:size(Q2,2)
%     Q(l,Q1(l)) = 1;
%     Q(l,Q2(l)) = -1;
% end
% Q = sparse(Q);
% only need the value of QtQ
% change the Q calculation after Cinv to ensure memory
Z_temp = Z'*Z;
Z_temp = single(full(Z_temp));
Z_temp = sparse(double(Z_temp));
C1_temp = (n*p*sigma + sigma)*sparse(eye(n*p))-sigma*ones(n*p,1)*ones(n*p,1)';
C2_temp = sigma*ones(n*p,1)*ones(n*p,1)';
C3_temp = (C1_temp-C2_temp);
clear("C1_temp",'C2_temp');
C = Z_temp + C3_temp;
clear("Z_temp",'C3_temp');
C = Z'*Z + (n*p*sigma + sigma)*sparse(eye(n*p)) -sigma*ones(n*p,1)*ones(n*p,1)';
Cinv = (inv(C));
clear('C');
QtQ = (n*p)*eye(n*p)+(-1)*ones(n*p,n*p);
Qs_w = kron(Q_origin,eye(p));
Qs_w = sparse(Qs_w);
% Qs = sparse(Qs_w'*Qs_w);
Qs = logical(Qs_w);
% Qs = Qs_w./Qs_w;
% Qs(isnan(Qs))=0;
QtQ = double(QtQ);
S = (QtQ) - Qs_w'*Qs_w;
clear("QtQ");
B0 = rand(n*p,K); U0 = rand(n*p,K); lam0 = rand(n*p,K);
V0 = rand(m*p,K); mu0 = rand(m*p,K);
B_old = B0; U_old = U0; lam_old = lam0; V_old = V0; mu_old = mu0;
iter = 0;
for k = 1:kmax
    iter = iter + 1;
    B_new = Cinv*(Z'*Y + sigma*(U_old+Qs_w'*V_old+S*B_old) - lam_old - Qs_w'*mu_old); 
    B_new = double(B_new);
    for j = 1:size(U_old,1)
        U_new(j,:) = max(1-alpha1/(sigma*norm(B_new(j,:)+lam_old(j,:)/sigma,2)),0)*(B_new(j,:)+lam_old(j,:)/sigma);
    end
    for j = 1:size(V_old,1)
        omega = Qs_w(j,find(Qs_w(j,:),1));
        V_new(j,:) = max(1-alpha2*omega/(sigma*norm(Qs(j,:)*B_new+mu_old(j,:)/sigma,2)),0)*(Qs(j,:)*B_new+mu_old(j,:)/sigma);
    end
    lam_new = lam_old + tau*sigma*(B_new-U_new);
    mu_new = mu_old + tau*sigma*(Qs*B_new-V_new);
    vareps1 = (norm(Qs*B_new-V_new,'fro')+norm(B_new-U_new,'fro'))/(1+norm(V_new,'fro')+norm(U_new,'fro'));
    vareps2 = norm(Z'*(Z*B_new-Y)+Qs'*mu_new+lam_new,"fro")/(1+norm(Z'*Y,'fro'));
    if (vareps1 < tol) | (vareps2<tol)
        break
    end
    B_old = B_new; U_old = U_new; V_old = V_new; lam_old = lam_new; mu_old = mu_new;
end
B = B_old; iter = iter;