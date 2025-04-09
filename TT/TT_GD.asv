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

%%%%%%%%
r_old = 2 * (b - multi_r1_times_TT(A,x)) ;
dx_TT = 0;
first_iteration = 1;
%%%%%%%%

for epoch = 1:max_epoches
    %shuffle
%     new_order = randperm(n_samples);
%     A = A(:,:,new_order);
%     b = b(new_order);
    for j = 1:batch_size:n_samples-batch_size+1
        A_j = A(:,:,j:j+batch_size-1);
        b_j = b(j:j+batch_size-1);
        
        r = b_j - multi_r1_times_TT(A_j,x) ;
%         [V,dUx] = TT_Newton_Gradient(A_j,x,r);
        
        %%%%%%%%%%%%%%%
        if first_iteration
            beta = -1;
            first_iteration = 0;
        else
            beta = (r'*r)/(r_old'*r_old);
        end
        [V,dUx] = TT_Newton_Gradient_Momentum(A_j,x,r,beta,dx_TT);
        dx_TT = TT_Riemannian_fromGTensor(x,V,dUx);
        r_old = r;
        %%%%%%%%%%%%%%
        % if first_iteration
        %     dUx = TT_Riemannian_Gauge_update(x,V,dUx);
        %     %first_iteration = false;
        % else
        %     dUx = TT_Riemannian_Gauge_update(x,V,dUx);
        %     % [Q_l,Q_r] = TT_Riemannian_projection(x,V,U_old,V_old);
        %     dUx = TT_Riemannian_search_direction_update(x,V,dUx,g,g_norm,g_norm_old);
        %     %dUx = TT_Riemannian_search_direction_update_old(x,dUx,dUx_old,Q_l,Q_r,g_norm,g_norm_old);
        %     %dUx = TT_Riemannian_Gauge_update(x,V,dUx);
        % 
        % 
        % end





        %compute step size
        %min ||b-A(x+alpha*g)||  => min ||r - alpha*A*g||
        % g = TT_Riemannian_fromGTensor(x,V,dUx);
        % Ag = multi_r1_times_TT(A_j,g);
        % alpha = Ag'*r/(Ag'*Ag);


        %update
        x = TT_Riemannian_update(x,V,dUx,1,rank);

    end

    r_test = multi_r1_times_TT(A_test,x) - b_test;
    train_err = norm(multi_r1_times_TT(A_j,x)- b_j)/norm(b_j)
    test_err = norm(r_test)/norm(b_test)
    if test_err < tol
        break
    end

end


end

