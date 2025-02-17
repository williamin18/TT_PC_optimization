function [dUx] = TT_Riemannian_search_direction_update(U,V,dUx,g_old,g_norm,old_g_norm)
%TT_RIEMANNIAN_UPDATE_DIRECTION Summary of this function goes here
%   Detailed explanation goes here



[d,m,~] = TTsizes(U);



Y_l = cell(d,1);
Y_l{1} = 1;
Y_r = cell(d,1);
Y_r{d} = 1;

% gi2 = g_old{1};
% for i = 1:d-1
%     Y_l{i+1} = U{i}'*gi2;
%     gi2 = h2v( Y_l{i+1}* v2h(g_old{i+1},m(i+1)) , m(i+1) );
% end

for i = 1:d-1
    gi = h2v( Y_l{i}* v2h(g_old{i},m(i)) , m(i) );
    Y_l{i+1} = U{i}'*gi;
end

for i = d:-1:2
    gi = g_old{i}*Y_r{i};
    Y_r{i-1} = v2h(gi,m(i)) * v2h(V{i},m(i))';
end


dUx2 = cell(d,1);
for i = 1:d
    dUx2{i} =  h2v( Y_l{i}* v2h(g_old{i},m(i)) , m(i))* (Y_r{i});
end

beta = g_norm/old_g_norm;
for i = 1:d
    dUx{i} = dUx{i} + beta * h2v( Y_l{i}* v2h(g_old{i},m(i)) , m(i))* (Y_r{i});
end



end

