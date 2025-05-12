function [y] = formRank1Tensor(r1_idx,r1_samples,m,d)
%FORMRANK1TENSOR Summary of this function goes here
%   Detailed explanation goes here
y = cell(d,1);
ref_points = r1_idx(1,:);
ref_out = r1_samples(1);
counter = 2;
for i = 1:d
    temp = r1_samples(counter:(counter+m-1));
    y{i} = [temp(1:(ref_points(i)-1)); ref_out ;temp(ref_points(i):m)]/ref_out;
end
y = TTorthogonalizeLR(y);
end

