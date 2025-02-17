function [dUx] = TT_Riemannian_Gauge_update(U,V,dUx)

[d,m,~] = TTsizes(U);

%project first d-1 gradient cores to orthogonal complements of the orthogonal cores
for i = 1:d-1
    L = U{i}'*dUx{i};
    dUx{i} = dUx{i} - U{i}*L;
    dUx{i+1} = dUx{i+1} + h2v(L*v2h(V{i+1},m(i+1)),m(i+1));
end

% [Q,~] = qr(U{d},'econ');
% dUx{d} = dUx{d} - Q*Q'*dUx{d};

end