clear variables
freq = linspace(0.5e9,14e9,41);
training_data = load('Tests/Filter/data_training.mat');
mc_data = load('Tests/Filter/data_mc.mat');
S11_train = 20*log10(abs(training_data.S11));
training_samples = normcdf(training_data.samples)*2-1;
S11_mc = 20*log10(abs(mc_data.S11));
mc_samples = normcdf(mc_data.samples)*2-1;

freq_low = zeros(1000,1);
freq_low_test = zeros(1000,1);
for i = 41:-1:1
    for j = 1:1000
        if S11_train(j,i) < -10
            freq_low(j) = freq(i);
        end
        if S11_mc(j,i) < -10
            freq_low_test(j) = freq(i);
        end
    end
end
% pdf_freq_idx = 28;
% f = figure(5);
% Hmc = histogram( abs(S21_train(:,pdf_freq_idx)) ,50,'Normalization','pdf', 'DisplayStyle','bar', 'FaceColor',[0.7 0.7 0.7]);

[n_samples,d] = size(training_samples);
m = 5;

% vouts_total = pc_collocation_total(training_samples,vouts_train,samples,m ,'Hermite');
% norm(vouts_total-vouts)/norm(vouts);

N = (m+1)*ones(d,1);
r = 3;
x = TTrand(N,r);
x = TTorthogonalizeLR(x);
x{d} = x{d}/norm( x{d},'fro');
x{1}(1) = freq_low(1,1);
for i = 1:d
    x{i}(1)=1;
end

tic
% [y_out,PC_coefficients] = pc_collocation_total(training_samples,S21_train,mc_data.samples,m,'Legendre');
[y_out,PC_coefficients,training_err,test_err,n_iterations] = ...
    pc_collocation_tensor_optimization(training_samples,freq_low,x,mc_samples,m,'Hermite','TT-ALS',0.3,0.2,r);
toc
