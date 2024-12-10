function [V,dUx] = TT_Riemannian_Gradient(A,U,b,step_size)
%TT_RIEMANNIAN_GD computes the Riemannian gradient in implicit format
%This function is based on https://sma.epfl.ch/~anchpcommon/publications/ttcompletion.pdf
%   Inputs :
% 1. A contains the totl order sample polynomials
% for m = 4, d = 10, 600 samples, the size of A is 600 * (4^10),
% each row of A is the kronecker product of ten 4*1 vectors
% we dont want to compute A explicitly, so we store A as a 10*4*600 tensor
% 2. U is the initial guess of unknown in left orthogonal
% TT-format(except the last TT-core)
% 3. b is the vector of sample outputs 
% 4. step_size 
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
[d,m,r] = TTsizes(x);
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
yl = cell(d,1);
yr = cell(1,d);
yl{1} = ones(n_samples,1);
yr{d} = ones(n_samples,1);


for i = 1:d-1
    xi = reshape(U{i},[r(i), m(i), r(i+1)]);
    xi = reshape(permute(xi, [2 1 3]),m(i),[]);
    Axi = A(i,:,:)'*xi;
    Axi = reshape(Axi,n_samples,r(i),r(i+1));

    temp = zeros(n_samples,r(i+1));
    for j = 1:r(i+1)
        temp(:,j) = sum(yl{i}.*Axi(:,:,j),2);
    end
    yl{i+1} = temp;

end
for i = d:-1:2
    xi = reshape(V{i},[r(i), m(i), r(i+1)]);
    xi = reshape(permute(xi, [2 1 3]),m(i),[]);
    Axi = A(i,:,:)'*xi;
    Axi = reshape(Axi,n_samples,r(i),r(i+1));
    temp = zeros(n_samples,r(i-1));
    for j = 1:n_samples
        temp(j,:) = Axi(j,:,:)*yl{i}(:,j);
    end
    yr{k-1} = temp;
end



%Compute Parial derivatives, dUx_k_j = sum(r(i)A(k,j,i)[yr{k}(:,i)'*yl{k}(i,:)]')
dUx = cell(d,1);


end

