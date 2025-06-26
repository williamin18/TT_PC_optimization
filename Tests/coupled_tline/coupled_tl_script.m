clear variables

load('Tests/coupled_tline/coupled_tl_28par_2.mat')

[n_samples,d] = size(training_samples);
training_samples = training_samples(1:870,:);
vouts_train = vouts_train(1:870,:);
m = 3;

% tic
% vouts_total = pc_collocation_total(training_samples,vouts_train,samples,m ,'Hermite');
% norm(vouts_total-vouts,"fro")/norm(vouts,"fro")
% toc

N = (m+1)*ones(d,1);
r = 3;
x = TTrand(N,r);
x = TTorthogonalizeLR(x);
x{d} = x{d}/norm( x{d},'fro');
x{1}(1) = vouts_train(1,1);
for i = 1:d
    x{i}(1)=1;
end

% x = cell(d,1);
% for i = 1:d
%     x{i} = zeros(m+1,1);
% end
tic
 [vouts_TT,PC_coefficients,training_err,test_err,n_iterations] = ...
     pc_collocation_tensor_optimization(training_samples,vouts_train,x,samples,m,'Hermite','TT-ALS',0.3,0.2,3);
toc
norm(vouts_TT-vouts,'fro')/norm(vouts,'fro')

vout_mean0 = mean(abs(vouts));
vout_sigma0 = std((vouts));

vout_mean1 = mean(abs(vouts_TT));
vout_sigma1 = std((vouts_TT));


fpoints = linspace(0,5e8,100);
scale = 1e9;
f = figure(2);
subplot(2,1,1);
hold on
plot(fpoints/scale,vout_mean0,'k-','LineWidth',3);
plot(fpoints/scale,vout_mean1,'g--','LineWidth',2);
legend('Monte Carlo simulation','TT-Newton')
grid on
xlabel('Frequency (GHz)','interpreter','LaTex')
ylabel('Mean of $ \mid V_{out} \mid $ (V)','interpreter','LaTex')
set(gca,'GridLineStyle','--')
set(gca, 'FontName', 'Times New Roman')
set(gca,'FontSize',12)
box on;


subplot(2,1,2);
hold on
plot(fpoints/scale,vout_sigma0,'k-','LineWidth',3);
plot(fpoints/scale,vout_sigma1,'g--','LineWidth',2);
ylim([0.0 0.04])
grid on
ylabel('Standard devation of $ V_{out} $ (V)','interpreter','LaTex')
xlabel('Frequency (GHz)','interpreter','LaTex')
set(gca,'GridLineStyle','--')
set(gca, 'FontName', 'Times New Roman')
set(gca,'FontSize',12)
box on;

f.Position = [100 100 675 500];
