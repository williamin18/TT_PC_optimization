clear variables

m = 3;
d = 8;
xi_roots = genPolynomialRoots(m+1,'Hermite');

n_samples = 1000;
sample_indices = randi(m+1,[n_samples d]);
sample_xi = xi_roots(sample_indices);
