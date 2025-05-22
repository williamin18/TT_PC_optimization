function [x,training_err,test_err,epoch] = TT_SGD_linear(A,b,x,r_round,tol,max_epoches,A_test,b_test,lambda,batch_size)

% (A,x,b,preconidtioner,r_round,batch_size,n_epochs,err_max)
d = length(A);
[n_samples,~] = size(A{1});
[~,m,~] = TTsizes(x);


break_counter = 0;
break_limit = 5;
err_old = 100;

for epoch = 1:n_epochs
    new_order = randperm(n_samples);
    for i = 1:d
        A{i} = A{i}(new_order,:);
    end
    b = b(new_order);
    for j = 1:batch_size:n_samples-batch_size+1
        A_j = cell(d,1);
        for i = 1:d
            A_j{i} = A{i}(j:j+batch_size-1,:);
        end
        b_j = b(j:j+batch_size-1);

        % preconditioner = ones(batch_size,1);
        % for i = 1:n_d
        %     preconditioner = preconditioner./reshape(dot(A_j(i,:,:),A_j(i,:,:),2),batch_size,1);
        % end
        % preconditioner = preconditioner.^0.5;
        % preconditioner = diag(preconditioner);
        % A_j(1,:,:) = reshape(A_j(1,:,:),m,batch_size)*preconditioner;
        % b_j = preconditioner*b_j;

        r =  b_j - multi_r1_times_TT(A_j,x);

        df = multi_r1_times_vec_to_TT(A_j,r);
        step_size = TTnorm(df);

        training_err = norm(r)/norm(b);
        r_test = multi_r1_times_TT(A_test,x) - b_test;
        test_err = norm(r_test)/norm(b_test);

        if test_err < tol || test_err/training_err>5
            break
        end
        if   err_old-training_err < tol/1000
            break_counter = break_counter+1;
            if break_counter > break_limit
                break
            end
        else
            break_counter = 0;
        end
        err_old = training_err;




        step_size = r'*r/TTdot(df,df);
        %step_size = 1;
        x = TTaxby(1,x,step_size,df);

        % x = TTrounding(x,1e-4,r_round);
        x = TTrounding_Randomize_then_Orthogonalize(x,[1 r_round*ones(1,n_d-1) 1]);

        err = norm(r)/norm(b_j)
    end

end
end

