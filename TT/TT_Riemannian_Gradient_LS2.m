function [V,dUx,g_norm] = TT_Riemannian_Gradient_LS2(A,U,residual)
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


%Compute direction that minimizes the residual for each TT-core
dUx = cell(d,1);

g_succeed = 0;
r_newtion = residual;
for i = 1:d

    %minimize yl*(a_i*U_i)*yr-r
    dUx{i} = zeros(r(i)*m(i),r(i+1));

    Y_yl = kron( kron(ones(1,r(i+1)), ones(1,m(i))) , yl{i} );
    Y_Ai = kron( kron(ones(1,r(i+1)), reshape(A(i,:,:),m(i),n_samples)') , ones(1,r(i)));
    Y_yr = kron( kron(yr{i} ,ones(1,m(i)) ) , ones(1,r(i)));

    Y = Y_yl.*Y_Ai.*Y_yr;
    
    gi =  Y\residual;

    
    gi =  reshape(gi,[r(i)*m(i) r(i+1)])+g_succeed;
    if i~=d
        L = U{i}'*gi;
        gi = gi - U{i}*L;
        g_succeed = h2v(L*v2h(V{i+1},m(i+1)),m(i+1));
    end

    dUx{i} = gi;

    % 
    % dUxi = reshape(dUx{i}, r(i),m(i),r(i+1));
    % dUxi = reshape(permute(dUxi,[2 1 3]),m(i),[]);
    % AdU = reshape(A(i,:,:), m(i),n_samples)'*dUxi;
    % AdU = reshape(AdU,n_samples,r(i),r(i+1));
    % 
    % Yi = zeros(n_samples,r(i+1));
    % for j = 1:r(i+1)
    %     Yi(:,j) = sum(yl{i}.*AdU(:,:,j),2);
    % end
    % Yi = sum(Yi.*yr{i},2);
    % alpha = Yi'*residual/(Yi'*Yi);
    % dUx{i} = alpha*dUx{i};
    % 
    % residual = residual-alpha*Yi;
end


g_norm = 0;
for i = 1:d
    g_norm = g_norm + norm(dUx{i},'fro');
end


end

