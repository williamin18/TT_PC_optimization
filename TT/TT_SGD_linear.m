function [x,training_err,test_err,n_iterations] = TT_SGD_linear(A,b,x,r_round,tol,max_epoches,A_test,b_test,lambda,batch_size)

% (A,x,b,preconidtioner,r_round,batch_size,n_epochs,err_max)
d = length(A);
[n_samples,~] = size(A{1});
[~,m,~] = TTsizes(x);


break_counter = 0;
break_limit = 5;
err_old = 100;
beta = 0;
n_iterations = 0;

for epoch = 1:max_epoches
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
        % for i = 1:d
        %     preconditioner = preconditioner./reshape(dot(A_j{i},A_j{i},2),batch_size,1);
        % end
        % preconditioner = preconditioner.^0.5;
        % preconditioner = diag(preconditioner);
        % A_j{1} = preconditioner*A_j{1};
        % b_j = preconditioner*b_j;

        r =  b_j - multi_r1_times_TT(A_j,x);

        if(norm(r)/norm(b_j) < tol)
            break_counter = break_counter+1;
        else
            break_counter = 0;
        end
        if break_counter >= break_limit
            break
        end

        df = multi_r1_times_vec_to_TT(A_j,r);
        df = TTaxby(1,df,-lambda,x);


        Adf = multi_r1_times_TT(A_j,df);
        % step_size = (Adf'*r)/(Adf'*Adf);
        step_size = r'*r/TTdot(df,df);
        x = TTaxby(1,x,step_size,df);
        x = TTrounding_Randomize_then_Orthogonalize(x,[1 r_round*ones(1,d-1) 1]);

        % if beta <= 0 
        %     Adf = multi_r1_times_TT(A_j,df);
        %     step_size = (Adf'*r-lambda^2*TTdot(df,x))/((Adf'*Adf)+lambda^2*TTnorm(df)^2);
        %     x = TTaxby(1,x,step_size,df);
        %     x = TTrounding_Randomize_then_Orthogonalize(x,[1 r_round*ones(1,d-1) 1]);
        %     g = TTrounding_Randomize_then_Orthogonalize(df,[1 r_round*ones(1,d-1) 1]);
        %     beta = 1;
        % else
        %     Adf = multi_r1_times_TT(A_j,df);
        %     Ag_old = multi_r1_times_TT(A_j,g);
        %     Ag2 = [Adf Ag_old];
        %     AtA = Ag2'*Ag2;
        %     reg_matrix = zeros(2,2);
        %     reg_matrix(1,1) = TTdot(df,df);
        %     reg_matrix(2,2) = TTdot(g,g);
        %     reg_matrix(1,2) = TTdot(df,g);
        %     reg_matrix(2,1) = conj(reg_matrix(1,2));
        %     AtA = AtA + lambda^2*reg_matrix;
        %     step_sizes = AtA\(Ag2'*r-lambda^2*[TTdot(df,x);TTdot(g,x)]);
        %     x = TTsum_Randomize_then_Orthogonalize({x,df,g}, [1;step_sizes], tol, r_round);
        %     g = TTsum_Randomize_then_Orthogonalize({df,g}, step_sizes, tol, r_round);
        % end
        training_err = norm(r)/norm(b);
        r_test = multi_r1_times_TT(A_test,x) - b_test;
        n_iterations = n_iterations+1;

    end
    test_err = norm(r_test)/norm(b_test);
    if   err_old-test_err < tol/100
        break_counter = break_counter+1;
        if break_counter > break_limit
            break
        end
    else
        break_counter = 0;
    end
    err_old = test_err;
    if test_err < tol || break_counter > break_limit
        break
    end
end
end

