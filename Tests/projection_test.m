n = 4;
d = 10;
N = n*ones(d,1);
r = 3;

x = TTrand(N,r);
U = TTorthogonalizeLR(x);
x = TTrand(N,r);

V{d} = U{d};
Ux{d} = U{d};
for i = d:-1:2
    [Q, R] = qr(v2h(V{i}, n)', 'econ');
    V{i} = h2v(Q', n);
    V{i-1} = U{i-1} * R';
    %Ux{i-1} = V{i-1};
    Ux{i-1} = rand(size(V{i-1}));
end


x_test = TT_Riemannian_fromGTensor(U,V,Ux);

Ux2 = TT_Riemannian_projection(U,V,x_test);
x_test2 = TT_Riemannian_fromGTensor(U,V,Ux2);
err = TTnorm(TTaxby(1,x_test,-1,x_test2))/TTnorm(x_test2);