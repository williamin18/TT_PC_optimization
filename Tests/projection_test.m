n = 4;
d = 10;
N = n*ones(d,1);
r = 3;

x = TTrand(N,r);
U = TTorthogonalizeLR(x);

V{d} = U{d};
for i = d:-1:2
    [Q, R] = qr(v2h(V{i}, n)', 'econ');
    V{i} = h2v(Q', n);
    V{i-1} = U{i-1} * R';
    Ux{i-1} = V{i-1};
end


Ux2 = TT_Riemannian_projection(U,V,x);
