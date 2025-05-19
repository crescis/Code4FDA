function[Q] = buildQ(A_origin)
A = A_origin;
[n,~] = size(A);
edge_matrix = zeros(n,n);
weight_matrix = zeros(n,n);

for i=1:n-1
    diff = sum((repmat(A(i,:),n,1) - A) .* (repmat(A(i,:),n,1) - A),2);
    %diff = sum((A(i,:) - A) .* (A(i,:) - A),2);
    %diff = sum(abs(A(i,:) - A),2);
    [~,distance_id] = sort(diff,'ascend');
        for l=1:10                                   %%%%%%%%% number of neighbours
            if distance_id(l) == i || edge_matrix(i,distance_id(l)) == 1
                continue;
            end
            edge_matrix(i,distance_id(l)) = 1;
            
            %weight_matrix(i,distance_id(l)) = exp(-1*5*diff(distance_id(l),1));%for iris dataset
            %weight_matrix(distance_id(l),i) = exp(-1*5*diff(distance_id(l),1));%for iris dataset
            %weight_matrix(i,distance_id(l)) = exp(-2*diff(distance_id(l),1));%for moon dataset
            %weight_matrix(distance_id(l),i) = exp(-2*diff(distance_id(l),1));%for moon dataset
            weight_matrix(i,distance_id(l)) = 1/(sqrt(diff(distance_id(l),1))+0.1);%for moon dataset
            weight_matrix(distance_id(l),i) = 1/(sqrt(diff(distance_id(l),1))+0.1);%for moon dataset
            edge_matrix(distance_id(l),i) = 1;
        end
end
 
%build Q
m = sum(sum(edge_matrix))/2;
Q= zeros(m,n);
counter = 0;
for i=1:n-1
    for j=i+1:n
        if(edge_matrix(i,j) == 1)
            counter = counter + 1;
            Q(counter,i) = weight_matrix(i,j);
            Q(counter,j) = -weight_matrix(i,j);
        end
    end
end

[~,d] = size(A);
alpha = 0.2*d;
Q = alpha*Q;
end