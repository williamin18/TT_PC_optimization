function [x] = TT_ALS_gd_old(A,b,x,batch_size,rank,tol,max_epoches,A_test,b_test)
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

g_old_project = 0;
norm_g_old = 1;

for epoch = 1:max_epoches
    %shuffle
    new_order = randperm(n_samples);
    A = A(:,:,new_order);
    b = b(new_order);
    for j = 1:batch_size:n_samples-batch_size+1
        A_j = A(:,:,j:j+batch_size-1);
        b_j = b(j:j+batch_size-1);

        residual = b_j - multi_r1_times_TT(A_j,x);
        
        %compute the product of A with the fixed TT-cores to compute the
        %partial derivative to update the current TT-core

        [yl,~] = Ax_left(A_j,x,core_idx);    %product of cores of Ax with index larger than the current core index
        [yr,~] = Ax_right(A_j,x,core_idx);   %product of cores of Ax with index smaller than the current core index
        
        
        g = zeros(r(core_idx)*m(core_idx),r(core_idx+1));%partial derivative 
        for i = 1:m(core_idx)
            g((i-1)*r(core_idx)+1:i*r(core_idx),:) = (yr* diag(residual.*reshape(A_j(core_idx,i,:),batch_size,1))*yl)';
        end

        norm_g = norm(g,'fro');
        % g = 0.2*g + 0.8*g_old_project;

        g = g + norm_g/norm_g_old*g_old_project;
        norm_g_old = norm_g;
      

        %compute step size that minizmize Ax_(k+1)-b = r - alpha*(yl* (Ak *g)*yr )
        
        g2 = reshape(g, r(core_idx),m(core_idx),r(core_idx+1));
        g2 = reshape(permute(g2,[2 1 3]),m(core_idx),[]);
        Ag2 = reshape(A_j(core_idx,:,:),[m(core_idx) batch_size])'*g2;
        Ag2 = reshape(Ag2, [batch_size, r(core_idx), r(core_idx+1)] );
        yg = zeros(batch_size,r(core_idx+1));
        for i = 1:r(core_idx+1)
            yg(:,i) = sum(yl.*Ag2(:,:,i),2);
        end
        yg = sum(yg.*yr',2);
        alpha = yg'*residual/(yg'*yg);

        x{core_idx} = x{core_idx} + alpha*g;

        %orthogonalize to update the next TT-core and project gradient to
        %the next TT-core
        if dir
            if core_idx == d-1
                dir = 0;
            end
            

            [x{core_idx}, R_k] = qr(x{core_idx},0);
            x{core_idx + 1} = h2v(R_k * v2h(x{core_idx + 1}, m(core_idx + 1)), m(core_idx + 1));
            core_idx = core_idx + 1;

            g_old_project = h2v( x{core_idx -1}'*g* v2h(x{core_idx},m(core_idx)) ,m(core_idx));
            

        else
            if core_idx == 2
                dir = 1;
            end
            
            [Q_k, R_k] = qr(v2h(x{core_idx}, m(core_idx))', 0);
            x{core_idx} = h2v(Q_k', m(core_idx));
            x{core_idx-1} = x{core_idx-1} * R_k';

            core_idx = core_idx -1;

            g_old_project = x{core_idx}*(v2h(g,m(core_idx+1)) * v2h(x{core_idx + 1},m(core_idx+1))');
        end
        
    end
    training_err = norm(b - multi_r1_times_TT(A,x))/norm(b)
    r_test = multi_r1_times_TT(A_test,x) - b_test;
    test_err = norm(r_test)/norm(b_test)
    if test_err < tol
        break
    end
end
end

