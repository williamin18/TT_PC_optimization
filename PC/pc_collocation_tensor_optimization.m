function [y_test,PC_coefficients] = pc_collocation_tensor_optimization(x_train,y_train,x_test,order,polynomial,method)
%PC_COLLOCATION_TENSOR Summary of this function goes here
%   Detailed explanation goes here


sample_polynomial_tensor = genPolynomialSamplesTensor(x_train,order,polynomial);

PC_coefficients = TT_ALS%(sample_polynomial_tensor,y_train)
end

