function [outputArg1,outputArg2] = TT_SGD_linear(A,x,b,preconidtioner,r_round,batch_size,n_epochs,err_max)
%TT_SGD_VEC Summary of this function goes here
%   Detailed explanation goes here
[n_d,m,n_samples] = size(A);

A = reshape(permute(A,[2 1 3]),m,[]);
A = r1_preconditioner*A;
A = reshape(A, [m,n_d,n_samples] );
A = permute(A,[2 1 3]);

break_counter = 0;
break_condition = 3;

for epoch = 1:n_epochs
    new_order = randperm(n_samples);
    A = A(:,:,new_order);
    b = b(new_order);
    for j = 1:batch_size:n_samples-batch_size+1
        A_j = A(:,:,j:j+batch_size-1);
        b_j = b(j:j+batch_size-1);

        preconditioner = ones(batch_size,1);
        for i = 1:n_d
            preconditioner = preconditioner./reshape(dot(A_j(i,:,:),A_j(i,:,:),2),batch_size,1);
        end
        preconditioner = preconditioner.^0.5;
        preconditioner = diag(preconditioner);
        A_j(1,:,:) = reshape(A_j(1,:,:),m,batch_size)*preconditioner;
        b_j = preconditioner*b_j;

        %for k = 1:2
        r =  b_j - multi_r1_times_TT(A_j,x);
        if(norm(r)/norm(b_j) < err_max)
            break_counter = break_counter+1;
        else
            break_counter = 0;
        end
        if break_counter >= break_condition
            break
        end
        df = multi_r1_times_vector_to_TT(A_j,r);
        % df = TTaxby(1,df,-1*lambda,x);
        step_size = r'*r/TTdot(df,df);
        %step_size = 1;
        x = TTaxby(1,x,step_size,df);

        % x = TTrounding(x,1e-4,r_round);
        x = TTrounding_Randomize_then_Orthogonalize(x,[1 r_round*ones(1,n_d-1) 1]);

        err = norm(r)/norm(b_j)
    end

end
end

