function b = multi_r1_times_vec_to_TT(A,x)
d = length(A);

[n_d,N,n_samples] = size(A);
b = cell(n_d,1);
b{1} = reshape(A(1,:,:),N,n_samples);
for j = 1:n_samples
    b{1}(:,j) = b{1}(:,j)*x(j);
end
for i = 2:n_d-1
    b{i} = zeros(n_samples*N,n_samples);
    for j = 1:N
        b{i}((j-1)*n_samples+1:j*n_samples,:) = diag(reshape(A(i,j,:),n_samples,1));
    end
end
b{n_d} =  reshape(permute(A(n_d,:,:),[3 2 1]),n_samples*N,1);
end