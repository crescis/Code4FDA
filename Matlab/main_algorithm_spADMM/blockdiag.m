function[A_final] = blockdiag(A_origin,M_pc)
A_final = [];
[iter,~] = size(A_origin);
for i = 1:iter
    A_temp = kron(A_origin(i,:),diag(ones(M_pc,1)));
    A_final = sparse(blkdiag(A_final,A_temp));
end
end

