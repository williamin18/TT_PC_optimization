function [x,training_err,test_err,epoch] = TT_Newton_GD(A,b,x,rank,tol,max_epoches,A_test,b_test,lambda)
%TT_NEWTON_GD Summary of this function goes here
%   Detailed explanation goes here
d = length(A);
[n_samples,~] = size(A{1});
[~,m,~] = TTsizes(x);

x = TTorthogonalizeLR(x);

r = b - multi_r1_times_TT(A,x);
dx_TT = 0;
beta = 0;

break_counter = 0;
break_limit = 5;
err_old = 100;

for epoch = 1:max_epoches

    [V,dUx] = TT_Newton_Gradient3(A,x,r,beta,dx_TT,lambda);

    %compute Momentum
    dx_TT = TT_Riemannian_fromGTensor(x,V,dUx);

    x = TT_Riemannian_update(x,V,dUx,1,rank);
    
    % r_old = r;
    r = b - multi_r1_times_TT(A,x);
    beta = 1;
    % beta = (r_old'*r_old)/(r'*r);

    
    

    training_err = norm(r)/norm(b);
    r_test = multi_r1_times_TT(A_test,x) - b_test;
    test_err = norm(r_test)/norm(b_test);

    if test_err < tol || test_err/training_err>4
        break
    end

    if   err_old-training_err < tol/1000*d
        break_counter = break_counter+1;
        if break_counter > break_limit
            break
        end
    else
        break_counter = 0;
    end
    err_old = training_err;
end
end

