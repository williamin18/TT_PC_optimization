clear variables

m = 3;
d = 8;
n_samples = 1000;


[sample_idx,sample_xi,n_samples] = genTensorCompletionSamples(m,d,n_samples,'Hermite');

[r1_sample_idx,r1_sample_xi] = genRank1Samples(m,d,'Hermite');