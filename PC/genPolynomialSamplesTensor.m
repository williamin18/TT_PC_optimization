function [y] = genPolynomialSamplesTensor(samples,order,polynomial)
%GENHERMITE Summary of this function goes here
%   Detailed explanation goes here


switch polynomial
    case "Hermite"
        f = @genHermite;
    case "Legendre"
        f = @genLegendre;
    otherwise
        err('Unsupported polynomial type')
end
[n,d]=size(samples);
samples = samples';
y = zeros(d,order,n);
for i = 0:order
    y(:,i+1,:) = f(samples,i);
end

end

