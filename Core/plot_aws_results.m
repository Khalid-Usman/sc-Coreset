k = 100;
clc
figure(k), clf, hold on

%%
disp(repmat('-',1,40))
coreset_sizes = [50 100 150 200 500];
D_approx_errs = [0.60 0.35 0.23 0.20 0.12];
%Ru_approx_errs = [7.02e-05 3.88e-05 3.85e-05];
%Rw_approx_errs = [7.02e-05 3.88e-05 3.85e-05];
SC_approx_errs = [0.06 0.07 0.03 0.06 0.02];

%(D_approx_errs < Ru_approx_errs) & (D_approx_errs < Rw_approx_errs)
%disp(repmat('-',1,40))

%%
% starting_idx = 1:length(coreset_sizes);
% starting_idx = [1 2 6 7 11 12 13] % k=10
% starting_idx = [1:3 5 8 9 11 13] % k=20
starting_idx = [1 2 3 4 5]; %k = 50

S = coreset_sizes(starting_idx);
Dx = D_approx_errs(starting_idx);
%Ux = Ru_approx_errs(starting_idx);
%Wx = Rw_approx_errs(starting_idx);
SCx = SC_approx_errs(starting_idx);
%Ux = Ru_approx_errs(starting_idx);
%Ue = Ru_approx_stds(starting_idx);
%Wx = Rw_approx_errs(starting_idx);
%We = Rw_approx_stds(starting_idx);

error_bar = [0.005 0.005 0.005 0.005 0.005];
%error_bar = error_bar(starting_idx);

[~,idx] = sort(S);

%%
plot(S(idx),SCx(idx),'og-','LineWidth',2)
errorbar(S(idx),Dx(idx),error_bar,'or-','LineWidth',2)
%errorbar(S(idx),Ux(idx),error_bar,'ob-','LineWidth',2)
%errorbar(S(idx),Wx(idx),error_bar,'om-','LineWidth',2)
%ylims = get(gca,'ylim');
%line([k k],[0 1],'LineStyle','--','Color',[0.5 0.5 0.5])
%set(gca,'ylim',ylims)
xlabel('Coreset size (number of points)','FontSize',12)
ylabel('Error','FontSize',12)
%title(sprintf('k = %d',k),'FontSize',16)
title('Zeisel, k=25','FontSize',12)
legend({'SC-Coreset', 'Uniform'},'FontSize',10)
%legend({'SC-Coreset', 'SVD-Coreset', 'Uniform Sampling', 'SVD-Coreset'},'FontSize',10)
axis([50 500 0 0.7])