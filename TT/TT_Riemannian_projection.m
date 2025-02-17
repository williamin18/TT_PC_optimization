function [Q_l,Q_r] = TT_Riemannian_projection(U,V,U2,V2)
%Projection of TT tensors  {U2_1 U2_2 ... U2_k-1 X V2_k+1 ...
%V2_d-1 V2_d} to TT tensors {U_1 U_2 ... U_k-1 Q_lk*X*Q_rk V_k+1 ... V_d-1
%V_d} for every k

[d,m,~] = TTsizes(U);

Q_l = cell(d,1);
Q_l{1} = 1;
Q_r = cell(d,1);
Q_r{d} = 1;

U2k2 = U2{1};
for i = 1:d-1
    U2k2 = h2v( Q_l{i}*v2h(U2k2,m(i)) , m(i) );
    U2k = U2k2;
    U2k2 = U2{i+1};
    Q_l{i+1} = U{i}'*U2k;
end

V2k2 = V2{d};
for i = d:-1:2
    V2k2 = V2k2*Q_r{i};
    V2k = V2k2;
    V2k2 = V2{i-1};
    Q_r{i-1} = v2h(V2k,m(i)) * v2h(V{i},m(i))';
end


end

