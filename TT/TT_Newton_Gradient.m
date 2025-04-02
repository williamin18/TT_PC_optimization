function [V,dUx] = TT_Newton_Gradient(A,U,residual)
%TT_NEWTON_GRADIENT Summary of this function goes here
%   Detailed explanation goes here
[~,~,n_samples] = size(A);
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
    dUx{i} = TT_Newton(yl{i},yr{i},reshape(A(i,:,:),[m(i),n_samples])',U{i},residual);

    % dUx{i} = zeros(r(i)*m(i),r(i+1));
    % for j = 1:m(i)
    %     temp = residual.*reshape(A(i,j,:),n_samples,1);
    %     temp = yr{i}.*repelem(temp,1,r(i+1));
    %     dUx{i}((j-1)*r(i)+1:j*r(i),:) = (temp'*yl{i})';  
    % end
    % %keep orthogonal direction
    % dUx{i} = dUx{i} - U{i}*U{i}'*dUx{i};
end

Adx = zeros(n_samples,d);
for i = 1:d
    dUi = reshape(dUx{i},[r(i), m(i), r(i+1)]);
    dUi = reshape(permute(dUi, [2 1 3]),m(i),[]);
    AdU = reshape(A(i,:,:), m(i),n_samples)'*dUi;
    AdU = reshape(AdU,n_samples,r(i),r(i+1));

    Adxi = zeros(n_samples,r(i+1));
    for j = 1:r(i+1)
        Adxi(:,j) = sum(yl{i}.*AdU(:,:,j),2);
    end
    Adx(:,i) = sum(Adxi.*yr{i},2);
end

alpha = Adx\residual;
for i = 1:d
    dUx{i} = alpha(i)*dUx{i};
end
end

