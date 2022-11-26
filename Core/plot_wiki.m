close all, clc
load wiki_unified_k1
i = argmin(cumulative_errors)
cumulative_errors(i) = [];
num_points_streamed(i) = [];
i = argmin(cumulative_errors)
cumulative_errors(i) = [];
num_points_streamed(i) = [];
i = argmin(cumulative_errors)
cumulative_errors(i) = [];
num_points_streamed(i) = [];
i = argmin(cumulative_errors)
cumulative_errors(i) = [];
num_points_streamed(i) = [];
figure, hold on
plot(num_points_streamed,log10(cumulative_errors),'LineWidth',2)
load wiki_unified_k10
plot(num_points_streamed,log10(cumulative_errors),'LineWidth',2)
load('wiki_unified_k100_copy')
plot(num_points_streamed,log10(cumulative_errors),'LineWidth',2)
title('Wikipedia approximation log error','FontSize',18)
xlabel('Number of million points streamed','FontSize',16)
ylabel('log_{10} eps','FontSize',16)
legend({'k = 1','k = 10','k = 100'},'Location','SE','FontSize',16)
set(gca,'XTickLabel',{'0','0.5','1','1.5','2','2.5','3','3.5','4'})
axis([0 3.8e6 -5 0.4])