%#ok<*SUSENS>

figure(2), clf, hold on

plot(-9999,-9999,'or-','LineWidth',3)
plot(-9999,-9999,'ok-','LineWidth',3)
legend({'MATLAB svds','SVD Coreset'},'Location','N','FontSize',24)
axis([0 5000 0 800])

%%
load coreset_runtimes_0204x
plot(K,runtimes/60,'ok-','LineWidth',3)
xlabel('Approximation rank k','FontSize',16)
ylabel('Running time (min)','FontSize',16)
title(sprintf('A[%d x %d], sparsity=%1.3f',10000,100000,0.001),'FontSize',18)

%%
load matlab_runtimes_0127_3x
plot(K,runtimes/60,'r-','LineWidth',3)
plot(K(1:end-1),runtimes(1:end-1)/60,'sr','LineWidth',3,'MarkerSize',5)
plot(K(end),runtimes(end)/60,'xr','LineWidth',3,'MarkerSize',20)
annotation('textarrow',[0.285,0.285],[0.285,0.225],'String',sprintf(' MATLAB\ncrashed'),'FontSize',16)
