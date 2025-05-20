function [dXk] = TT_Newton_Riemannian(yl,yr,Ak,Xk,residual,lambda)
%optimize yl*(Ak*dXk)*yr = r s.t. dXk'Xk = 0

[n_samples,m] = size(Ak);
[~,r_k] = size(yl);
[~,r_k1] = size(yr);

Y_yl = kron( kron(ones(1,r_k1), ones(1,m)), yl );
Y_Ai = kron( kron(ones(1,r_k1), Ak ),       ones(1,r_k));
Y_yr = kron( kron(yr,           ones(1,m)), ones(1,r_k));

Y = Y_yl.*Y_Ai.*Y_yr;

n = r_k*m*r_k1;
%regularization
if lambda ~= 0
    Y = [Y;lambda*eye(n)];
    residual = [residual; -lambda*reshape(Xk,[],1)];
    % residual = [residual; zeros(n,1)];
end

X = kron(eye(r_k1),Xk');
YtY = Y'*Y;
lambda =  (X/YtY*X') \ (X/YtY*Y'*residual);
dXk = (YtY)\(Y'*residual-X'*lambda);


dXk = reshape(dXk,[r_k*m, r_k1]);



end
