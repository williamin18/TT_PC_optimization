function [dUx] = TT_Riemannian_search_direction_update_old(U,dUx,dUx_old,Q_l,Q_r,g_norm,old_g_norm)
%TT_RIEMANNIAN_UPDATE_DIRECTION Summary of this function goes here
%   Detailed explanation goes here
[d,m,~] = TTsizes(U);
beta = g_norm/old_g_norm;
for i = 1:d
    dUx{i} = dUx{i} + beta * h2v( Q_l{i}* v2h(dUx_old{i},m(i)) , m(i))* (Q_r{i});
end



end

