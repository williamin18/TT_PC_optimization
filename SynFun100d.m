function [y] = SynFun100d(xi)
%UNTITLED Summary of this function goes here
%   xi is 100*1 uniform distributed between [-1,1]
x = xi/2+1.5;
x(20) = xi(20)+2;
k = 1:100;
y = 3-5/100*k*x+1/100*k*x.^3 + x(1)*x(2)^2+x(2)*x(4)-x(3)*x(5)+x(51)+x(50)*x(54)^2+log(1/300*k*(x.^2+x.^4));
%y = 3-5/100*k*x+1/100*k*x.^3 + x(1)*x(2)^2+x(2)*x(4)-x(3)*x(5)+x(51)+x(50)*x(54)^2+1/300*sum(k'.*log((x.^2+x.^4)));

end

