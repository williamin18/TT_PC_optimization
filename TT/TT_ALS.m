function [x] = TT_ALS(A,b,x,batch_size,rank,tol,max_epoches,A_test,b_test)
%TT_ALS Summary of this function goes here
%   Detailed explanation goes here

[d,~,n_samples] = size(A);
[~,m,r] = TTsizes(x);


x = TTorthogonalizeRL(x); %orthogonalize

core_idx = 1; %index for TT-core updated at current iteration
dir = 1; %direction to update next TT-core

%test codes
x0 = x;
r_test = multi_r1_times_TT(A_test,x) - b_test;
test_err = norm(r_test)/norm(b_test)

for epoch = 1:max_epoches
    %shuffle
    new_order = randperm(n_samples);
    A = A(:,:,new_order);
    b = b(new_order);
    for j = 1:batch_size:n_samples-batch_size+1
        A_j = A(:,:,j:j+batch_size-1);
        b_j = b(j:j+batch_size-1);

        residual = b_j - multi_r1_times_TT(A_j,x)
        
        
        %compute the product of A with the fixed TT-cores to compute the
        %partial derivative to update the current TT-core

        [yl,~] = Ax_left(A_j,x,core_idx);    %product of cores of Ax with index larger than the current core index
        [yr,~] = Ax_right(A_j,x,core_idx);   %product of cores of Ax with index smaller than the current core index
        
        
        g = zeros(r(core_idx)*m(core_idx),r(core_idx+1));%partial derivative 
        for i = 1:m(core_idx)
            g((i-1)*r(core_idx)+1:i*r(core_idx),:) = (yr* diag(residual.*reshape(A_j(core_idx,i,:),batch_size,1))*yl)';
        end

        %update 
        
        alpha = norm(residual)^2/norm(g,'fro')^2;
        x{core_idx} = x{core_idx} + alpha*g;

        %orthogonalize to update the next TT-core
        if dir
            if core_idx == d-1
                dir = 0;
            end
            

            [x{core_idx}, R_k] = qr(x{core_idx},0);
            x{core_idx + 1} = h2v(R_k * v2h(x{core_idx + 1}, m(core_idx + 1)), m(core_idx + 1));
            core_idx = core_idx + 1;
        else
            if core_idx == 2
                dir = 1;
            end
            
            [Q_k, R_k] = qr(v2h(x{core_idx}, m(core_idx))', 0);
            x{core_idx} = h2v(Q_k', m(core_idx));
            x{core_idx-1} = x{core_idx-1} * R_k';

            core_idx = core_idx -1;

        end
        
    end
    r_test = multi_r1_times_TT(A_test,x) - b_test;
    test_err = norm(r_test)/norm(b_test);
    if test_err < tol
        break
    end
end
end

