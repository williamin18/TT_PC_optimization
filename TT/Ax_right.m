function [yrk,yr] = Ax_right(A,x,k)
%multiply rank 1 samples A with TT x for cores indices greater than k, 
% yr{k} is a r_k*n matrix: n is number of samples, r_k is the rank of k's TT-core

[~,~,n_samples] = size(A);
[d,m,r] = TTsizes(x);

yl = cell(d-k+1,1);
yr{d-k+1} = ones(1,n_samples);

for i = d:-1:k+1
    xi = reshape(x{i},[r(i), m(i), r(i+1)]);
    xi = reshape(permute(xi, [2 1 3]),m(i),[]);
    Axi = reshape(A(i,:,:),[m(i),n_samples])'*xi;
    Axi = reshape(Axi,n_samples,r(i),r(i+1));
    temp = zeros(r(i),n_samples);
    for j = 1:n_samples
        temp(:,j) = reshape(Axi(j,:,:),r(i),r(i+1))*yr{i-k+1}(:,j);
    end
    yr{i-k} = temp;
end
yrk = yr{1};
end

