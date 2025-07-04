clear variables


training_data = load('Tests/Filter/data_training.mat');
training_samples = training_data.samples;
[n_samples,~] = size(training_samples);
S21= ((training_data.S21));
[~,n_freq] = size(S21);
freq = linspace(0.5e9,14e9,n_freq);

S21_abs = abs(S21);
S21_dB = 20*log10(abs(S21));



S21_sort = sort(S21_dB);
S21_98 = S21_sort(round(0.02*n_samples),:);
S21_02 = S21_sort(round(0.98*n_samples),:);
S21_50 = S21_sort(round(0.5*n_samples),:);
S21_84 = S21_sort(round(0.16*n_samples),:);
S21_16 = S21_sort(round(0.84*n_samples),:);

figure(1)
hold on
grid on
fill([freq,fliplr(freq)], [S21_02,fliplr(S21_98)],'k','FaceAlpha',.2,EdgeColor='none');
plot(freq, [S21_50],'k');
plot(freq,[S21_16;S21_84],'r--')

% plot(freq,[S21_02;S21_50; S21_98]);

% s_max = max(S21_abs,[],'all');
% s_min = min(S21_abs,[],'all');
% s_scale = (1-2/n_samples)/(s_max-s_min);
% s_bias = s_max*s_scale-(1-1/n_samples);
% s = S21_abs*s_scale-s_bias;
% s =  log(s./(1-s));
low_freq = zeros(n_samples,1);
high_freq = zeros(n_samples,1);
for i = 1:n_samples
    is_good_performance = S21_dB(i,:) > -3;
    good_indices = find(is_good_performance);
    j1 = good_indices(1);
    j2 = good_indices(end);
    low_freq(i) = ( freq(j1)*(-3-S21_dB(i,j1-1)) + freq(j1-1)*(S21_dB(i,j1)+3) )/(S21_dB(i,j1)-S21_dB(i,j1-1));
    high_freq(i) = ( freq(j2)*(-3-S21_dB(i,j2+1)) + freq(j2+1)*(S21_dB(i,j2)+3) )/(S21_dB(i,j2)-S21_dB(i,j2+1));
end
pdf_freq_idx = 20;
f = figure(5);
Hmc = histogram( high_freq ,50,'Normalization','pdf', 'DisplayStyle','bar', 'FaceColor',[0.7 0.7 0.7]);


[n_samples,d] = size(training_samples);
m = 4;
N = (m+1)*ones(d,1);
r = 3;
x = TTrand(N,r);
x = TTorthogonalizeLR(x);
x{d} = x{d}/norm( x{d},'fro');





[y_out,PC_coefficients,training_err,test_err,n_iterations] = ...
    pc_collocation_tensor_optimization(training_samples,low_freq,x,training_samples,m,'Hermite','TT-Newton',0.3,0.2,r);