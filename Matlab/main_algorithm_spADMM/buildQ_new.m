function[Q_new] = buildQ_new(Q,m,I)
for i = 1:m
    for j = 1:I
        [k1,k2] = find(Q(i,:));
        Q_new(I*(i-1)+j,(k1-1)*I+j) = Q(i,k1);
        Q_new(I*(i-1)+j,(k2-1)*I+j) = Q(i,k2);
    end
end