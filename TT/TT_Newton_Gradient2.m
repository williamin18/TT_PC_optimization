function [U,dUx] = TT_Newton_Gradient2(A,U,residual)
%TT_NEWTON_GRADIENT Summary of this function goes here
%   Detailed explanation goes here
[~,~,n_samples] = size(A);
[d,m,r] = TTsizes(U);
V = U;


[~,yl] = Ax_left(A,U,d);
[~,yr] = Ax_right(A,U,1);

dUx = cell(d,1);

for i = 1:d
    dUx{i} = TT_Newton(yl{i},yr{i},reshape(A(i,:,:),[m(i),n_samples])',U{i},residual);
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
% alpha = lsqnonneg(Adx,residual);

for i = 1:d
    dUx{i} = alpha(i)*dUx{i};
end
end

