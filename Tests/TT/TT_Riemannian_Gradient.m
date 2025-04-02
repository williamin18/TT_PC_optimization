function [V,dUx,g_norm] = TT_Riemannian_Gradient(A,U,residual)
%TT_RIEMANNIAN_GD computes the Riemannian gradient in implicit format
%This function is based on https://sma.epfl.ch/~anchpcommon/publications/ttcompletion.pdf
%   Inputs :
% 1. A contains the totl order sample polynomials
% for m = 4, d = 10, 600 samples, the size of A is 600 * (4^10),
% each row of A is the kronecker product of ten 4*1 vectors
% we dont want to compute A explicitly, so we store A as a 10*4*600 tensor
% 2. U is the initial guess of unknown in left orthogonal
% TT-format(except the last TT-core)
% 3. residual is the residial vector Ax-b
%   Outputs: 
% 1. V is the right orthogonalized U (except the first TT-core)
% 2. dUx is the direvative of the non-orthogonal part 
% There is another variable Ux used in this function, the unknown can be
% written as U_{k-}Ux_{k}V_{k+} for any index k. The TT-cores in U are left
% orthogonal, TT-cores in V are right orthogonal, the only non-orthogonal
% TT-core is Ux_{k}. The Riemannian gradient is computed based on the partial
% derivative with respect to Ux_{k}
% 
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



%yl and yr are tempory vectors for efficiency, for one sample set i, 
%yl{k}*(A(i,k,:)*Ux{k})*yr{k} = b(i) for the ideal solution 

[~,yl] = Ax_left(A,U,d);
[~,yr] = Ax_right(A,V,1);




%Compute Parial derivatives, dUx_k_j = sum(r(i)A(k,j,i)[yr{k}(:,i)'*yl{k}(i,:)]')
dUx = cell(d,1);
for i = 1:d
    dUx{i} = zeros(r(i)*m(i),r(i+1));
    for j = i:m(i)
        dUx{i}((j-1)*r(i)+1:j*r(i),:) = (yr{i}'* diag(residual.*reshape(A(i,j,:),n_samples,1))*yl{i})'/n_samples;
    end
end

g_norm = 0;
for i = 1:d
    g_norm = g_norm + norm(dUx{i},'fro');
end

%Orthogonalize
% for i = 1:d-1
%     % wrong codes
%     % [Q,~] = qr(Ux{i},'econ');
%     % dUx{i} = dUx{i}-Q*Q'*dUx{i};
%     L = U{i}'*dUx{i};
%     dUx{i} = dUx{i} - U{i}*L;
%     dUx{i+1} = dUx{i+1} + h2v(L*v2h(V{i+1},m(i+1)),m(i+1));
% end
% [Q,~] = qr(Ux{d},'econ');
% dUx{d} = dUx{d} - Q*Q'*dUx{d};






end

