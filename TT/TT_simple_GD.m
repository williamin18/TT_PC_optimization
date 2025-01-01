function [x,err] = TT_simple_GD(A,b,x,batch_size,max_rank,tol,n_epochs,A_test,b_test)
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