
function [x,err] = TT_sgd_optimization(obj,A,x,b,n_epochs)
%A contains the totl order sample polynomials,
%b is the vector of sample voltages
%x is the initial guess of PC-coefficients in TT-format
%for m = 4, d = 10, 600 samples, the size of A is 600 * (4^10),
%each row of A is the kronecker product of ten 4*1 vectors
%we dont want to compute A explicitly, so we store A as a 10*4*600 tensor
[n_d,m,n_samples] = size(A);
alpha = 0.001;
beta1 = 0.9;
beta2 = 0.999;
epsilon = 10^-8;
m_tm1 = cell(n_d,1);
v_tm1 = cell(n_d,1);
for i = 1:n_d
    m_tm1{i} = 0;
    v_tm1{i} = 0;
end
beta1_t = beta1;
beta2_t = beta2;
err = zeros(n_epochs,1);
for epoch = 1:n_epochs
    err_e = zeros(n_samples,1);
    for i = 1:n_samples

        %Compute the gradient of each TT-core seperately
        %y(i) = sum_j1(a1(j1)G1(j1))...sum_jd(ad(jd)Gd(jd))
        %objective function f = 1/2 (y(i)-b(i))^2
        %gradient for TT-core Gk

        [~,N,R] = TTsizes(x);
        Axk_list = cell(n_d,1);
        for k = 1:n_d
            xk = reshape(x{k},[R(k), N(k), R(k+1)]);
            xk = reshape(permute(xk, [2 1 3]),N(k),[]);
            Axk = A(k,:,i)*xk;
            Axk_list{k} = reshape(Axk,R(k),R(k+1));
        end


        yl = cell(n_d,1);
        yr = cell(n_d,1);
        yl{1} = 1;
        yr{n_d} = 1;
        for k = 1:n_d
            yl{k+1} = yl{k}*Axk_list{k};
        end

        for k = n_d:-1:2
            yr{k-1} = Axk_list{k}*yr{k};
        end



        y_m_b =  Axk_list{1}*yr{1} - b(i);
        err_e(i) = abs(y_m_b);

        x2 = x;
        %dfdX = cell{n_d};
        for k = 1:n_d
            g_t = zeros(R(k)*N(k),R(k+1));
            for j = 1:m
                %change to new TT-form
                g_t((j-1)*R(k)+1:j*R(k),:) = y_m_b*A(k,j,i)*(yr{k}*yl{k})';
            end
            m_t = beta1*m_tm1{k} + (1-beta1)*g_t;
            v_t = beta2*v_tm1{k} + (1-beta2)*(g_t.*g_t);
            m_tm1{k} = m_t;
            v_tm1{k} = v_t;

            m_h = m_t/(1-beta1_t);
            v_h = v_t/(1-beta2_t);
            beta1_t = beta1_t*beta1;
            beta2_t = beta2_t*beta2;

            x2{k} = x2{k} - alpha*m_h./(v_h.^0.5 + epsilon);
        end
        x = x2;
    end
    err(epoch) = mean(err_e);
end

end




function [x,err] = TT_sgd_optimization_simple2(obj,A,x,b,n_epochs)
%A contains the totl order sample polynomials,
%b is the vector of sample voltages
%x is the initial guess of PC-coefficients in TT-format
%for m = 4, d = 10, 600 samples, the size of A is 600 * (4^10),
%each row of A is the kronecker product of ten 4*1 vectors
%we dont want to compute A explicitly, so we store A as a 10*4*600 tensor
[n_d,m,n_samples] = size(A);
x = TTorthogonalizeRL(x);
alpha = 0.001;
err = zeros(n_epochs,1);
for epoch = 1:n_epochs
    err_e = zeros(n_samples,1);
    for i = 1:n_samples

        %Compute the gradient of each TT-core seperately
        %y(i) = sum_j1(a1(j1)G1(j1))...sum_jd(ad(jd)Gd(jd))
        %objective function f = 1/2 (y(i)-b(i))^2
        %gradient for TT-core Gk

        [~,N,r] = TTsizes(x);



        yr = cell(n_d,1);
        yr{n_d} = 1;
        for k = n_d:-1:2
            xk = reshape(x{k},[r(k), N(k), r(k+1)]);
            xk = reshape(permute(xk, [2 1 3]),N(k),[]);
            Axk = A(k,:,i)*xk;
            Axk = reshape(Axk,r(k),r(k+1));
            yr{k-1} = Axk*yr{k};
        end


        yl = cell(n_d,1);
        yl{1} = 1;



        for k = 1:n_d
            [~,N,r] = TTsizes(x);

            xk = reshape(x{k},[r(k), N(k), r(k+1)]);
            xk = reshape(permute(xk, [2 1 3]),N(k),[]);
            Axk = A(k,:,i)*xk;
            Axk = reshape(Axk,r(k),r(k+1));
            y =  yl{k}*Axk*yr{k};
            y_m_b = y - b(i);
            g_t = zeros(r(k)*N(k),r(k+1));
            for j = 1:m
                g_t((j-1)*r(k)+1:j*r(k),:) = y_m_b*A(k,j,i)*(yr{k}*yl{k})';
            end
            x{k} = x{k}-alpha*g_t;
            [Q_k, R_k] = qr(x{k}, 0);
            if k<n_d
                x{k+1} = h2v(R_k*v2h(x{k+1},N(k)),N(k));
                x{k} = Q_k;
                [~,N,r] = TTsizes(x);
                xk = reshape(x{k},[r(k), N(k), r(k+1)]);
                xk = reshape(permute(xk, [2 1 3]),N(k),[]);
                Axk = A(k,:,i)*xk;
                Axk = reshape(Axk,r(k),r(k+1));
                yl{k+1} = yl{k}*Axk;
            end

        end
        err_e(i) = abs(y_m_b);

    end
    err(epoch) = mean(err_e);
end

end



function [x,x_regularized,err] = TT_gd_linear(obj,A,x,b,r_round,batch_size,n_epochs,err_max)


[n_d,m,n_samples] = size(A);
lambda = 1/3;
if obj.M ==3
    r1_preconditioner = diag([1; lambda; lambda^2]);
else
    r1_preconditioner = diag([1; lambda; lambda^2; lambda^3]);
end
A = reshape(permute(A,[2 1 3]),m,[]);
A = r1_preconditioner*A;
A = reshape(A, [m,n_d,n_samples] );
A = permute(A,[2 1 3]);



break_counter = 0;
break_condition = 3;
% lambda = 0.01;
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
        r =  b_j - obj.multi_r1_times_TT(A_j,x);
        if(norm(r)/norm(b_j) < err_max)
            break_counter = break_counter+1;
        else
            break_counter = 0;
        end
        if break_counter >= break_condition
            break
        end
        df = obj.multi_r1_times_vector_to_TT(A_j,r);
        % df = TTaxby(1,df,-1*lambda,x);
        step_size = r'*r/TTdot(df,df);
        %step_size = 1;
        x = TTaxby(1,x,step_size,df);

        % x = TTrounding(x,1e-4,r_round);
        x = TTrounding_Randomize_then_Orthogonalize(x,[1 r_round*ones(1,obj.d-1) 1]);

        err = norm(r)/norm(b_j)
        %end
    end
    if break_counter >= break_condition
        break
    end
end
x_regularized = x;
r1_preconditioner = [1; lambda; lambda^2; lambda^3];
[~,~,R] = TTsizes(x);
for i = 1:n_d
    for j = 1:m
        x{i}((j-1)*R(i)+1:j*R(i),:) = x{i}((j-1)*R(i)+1:j*R(i),:)*r1_preconditioner(j);
    end
end
end