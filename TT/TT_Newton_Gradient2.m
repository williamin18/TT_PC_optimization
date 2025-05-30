function [V,dUx] = TT_Newton_Gradient2(A,U,residual,beta,dx_old,lambda)
%TT_NEWTON_GRADIENT Summary of this function goes here
%   Detailed explanation goes here
[n_samples,~] = size(A{1});
[d,m,r] = TTsizes(U);
V = U;
Ux = cell(d,1);
Ux{d} = U{d};


%Right Orthogonalize, find every Ux to compute the partial derivative
V{d} = U{d};
for i = d:-1:2
    [Q, R] = qr(v2h(V{i}, m(i))', 'econ');
    V{i} = h2v(Q', m(i));
    V{i-1} = U{i-1} * R';
    Ux{i-1} = V{i-1};
end


[~,yl] = Ax_left(A,U,d);
[~,yr] = Ax_right(A,V,1);

dUx = cell(d,1);

for i = 1:d
    dUx{i} = TTcore_Newton(yl{i},yr{i},A{i},Ux{i},residual,lambda);
end

for i = 1:d-1
    L = U{i}'*dUx{i};
    dUx{i} = dUx{i} - U{i}*L;
    % dUx{i+1} = dUx{i+1} + h2v(L*v2h(V{i+1},m(i+1)),m(i+1));
end

Adx = zeros(n_samples,d);

for i = 1:d
    dUi = reshape(dUx{i},[r(i), m(i), r(i+1)]);
    dUi = reshape(permute(dUi, [2 1 3]),m(i),[]);
    AdU = A{i}*dUi;
    AdU = reshape(AdU,n_samples,r(i),r(i+1));

    Adxi = zeros(n_samples,r(i+1));
    for j = 1:r(i+1)
        Adxi(:,j) = sum(yl{i}.*AdU(:,:,j),2);
    end
    Adx(:,i) = sum(Adxi.*yr{i},2);
end


if beta <= 0
    reg_matrix = zeros(d,d);
    for i = 1:d
        reg_matrix(i,i) = lambda*norm(dUx{i},'fro');
    end
    Adx = [Adx; reg_matrix];
    residual = [residual; zeros(d,1)];
    alpha = Adx\residual;
    for i = 1:d
        dUx{i} = alpha(i)*dUx{i};
    end
else
    dU_old = TT_Riemannian_projection_orthogonal(U,V,dx_old);
    Adx_old = zeros(n_samples,d);

    for i = 1:d
        dUi = reshape(dU_old{i},[r(i), m(i), r(i+1)]);
        dUi = reshape(permute(dUi, [2 1 3]),m(i),[]);
        AdU = A{i}*dUi;
        AdU = reshape(AdU,n_samples,r(i),r(i+1));

        Adxi = zeros(n_samples,r(i+1));
        for j = 1:r(i+1)
            Adxi(:,j) = sum(yl{i}.*AdU(:,:,j),2);
        end
        Adx_old(:,i) = sum(Adxi.*yr{i},2);
    end
    Adx = [Adx Adx_old];
    reg_matrix = zeros(d*2,d*2);
    for i = 1:d
        reg_matrix(i,i) = (lambda*norm(dUx{i},'fro'))^2;
        reg_matrix(i+d,i+d) = (lambda*norm(dU_old{i},'fro'))^2;
        reg_matrix(i,i+d) = lambda^2*sum(conj(dUx{i}).*dU_old{i},"all");
        reg_matrix(i+d,i) = conj(reg_matrix(i,i+d));
    end
    reg_matrix = chol(reg_matrix);
    Adx = [Adx; reg_matrix];
    residual = [residual; zeros(2*d,1)];
    alpha = Adx\residual;
    for i = 1:d
        dUx{i} = alpha(i)*dUx{i}+alpha(i+d)*dU_old{i};
    end
end

end

