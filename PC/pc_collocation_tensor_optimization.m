function [y_test,PC_coefficients] = pc_collocation_tensor_optimization(x_train,y_train,x_test,order,polynomial,method)
%PC_COLLOCATION_TENSOR Summary of this function goes here
%   Detailed explanation goes here
switch method
    case "ALS"
        f = @TT_ALS;
    case "TT-Newton"
        f = @TT_GD;
    otherwise
        err('Unsupported optimization type')
end

sample_polynomial_tensor = genPolynomialSamplesTensor(x_train,order,polynomial);



PC_coefficients = TT_ALS%(sample_polynomial_tensor,y_train)
test_sample_polynomial_tensor = genPolynomialSamplesTensor(x_test,order,polynomial);
y_test = multi_r1_times_TT(test_sample_polynomial_tensor,PC_coefficients);
end

