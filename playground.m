
% Set random seed for reproducibility
y_value = zeros(10000,1);
for i = 1:10000
% Generate xi values
xi = unifrnd(1, 2, 1, 100); % All ξ_k ~ U([1, 2])
xi(20) = unifrnd(1, 3);     % ξ_20 ~ U([1, 3])

% Compute y(ξ)

y_value(i) = compute_y(xi);
end
mean(y_value)

function y = compute_y(xi)
    % compute_y computes the value of y(ξ) as defined in the given equation.
    % 
    % Inputs:
    % xi - A vector of length d containing the ξ_k values
    %
    % Outputs:
    % y - The computed value of y(ξ)

    % Dimension
    d = 100;

    % Validate the input
    if length(xi) ~= d
        error('Input xi must be a vector of length d = 100.');
    end

    % Compute individual terms
    k = 1:d; % Index vector
    term1 = -(5 / d) * sum(k .* xi);
    term2 = (1 / d) * sum(k .* xi.^3);
    term3 = (xi(1)*xi(2)^2 + xi(2)*xi(4) -xi(3) * xi(5) + xi(51) + xi(50)*xi(54)^2); % Adjust indices for MATLAB (1-based indexing)
    term4 = log(1 / (3 * d) * sum(k .* (xi.^2 + xi.^4)));

    % Combine terms to compute y
    y = 3 + term1 + term2 + term3 + term4;
end