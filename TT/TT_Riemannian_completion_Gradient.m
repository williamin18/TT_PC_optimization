function [V,dUx,dx_TT] = TT_Riemannian_completion_Gradient(A,U,residual,beta,dx_old,lambda)
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
    for j = 1:m(i)
        temp = residual.*A{i}(:,j);
        temp = yr{i}.*repelem(temp,1,r(i+1));
        dUx{i}((j-1)*r(i)+1:j*r(i),:) = (yl{i}'*temp);
    end
end
for i = 1:d-1
    dUx{i} = dUx{i} - U{i}*U{i}'*dUx{i};
end

if beta>0
    dU_old = TT_Riemannian_projection(U,V,dx_old);
    for i = 1:d
        dUx{i}  = dUx{i} + beta*dU_old{i};
    end
end

dx_TT = TT_Riemannian_fromGTensor(U,V,dUx);

Adx = multi_r1_times_TT(A,dx_TT);
alpha = Adx\residual;
%alpha = lsqnonneg(Adx,residual);

for i = 1:d
    dUx{i} = alpha*dUx{i};
end
end

