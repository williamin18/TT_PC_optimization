function [x] = TT_GD(A,b,x,batch_size,rank,tol,max_epoches,A_test,b_test)
%TT_RIEMANNIAN_GD Summary of this function goes here
%   Detailed explanation goes here
% 1. A contains the totl order sample polynomials
% for m = 4, d = 10, 600 samples, the size of A is 600 * (4^10),
% each row of A is the kronecker product of ten 4*1 vectors
% we dont want to compute A explicitly, so we store A as a 10*4*600 tensor
% 2. b is the vector of corresponding output samples 
% 3. x is the initial guess in TT-format


[d,m,n_samples] = size(A);

x = TTorthogonalizeLR(x);
for epoch = 1:max_epoches
    %shuffle
    new_order = randperm(n_samples);
    A = A(:,:,new_order);
    b = b(new_order);
    for j = 1:batch_size:n_samples-batch_size+1
        A_j = A(:,:,j:j+batch_size-1);
        b_j = b(j:j+batch_size-1);
        
        r = b_j - multi_r1_times_TT(A_j,x) 
        [V,dU] = TT_Riemannian_Gradient(A_j,x,r);
        
        g_norm = 0;
        for i = 1:d
            g_norm = g_norm + norm(dU{i});
        end
        step_size = norm(r)/g_norm;
        x = TT_Riemannian_update(x,V,dU,step_size,rank);
    end
    r_test = multi_r1_times_TT(A_test,x) - b_test;
    test_err = norm(r_test)/norm(b_test)
    if test_err < tol
        break
    end

end


end

