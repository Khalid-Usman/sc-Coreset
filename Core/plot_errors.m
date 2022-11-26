%#ok<*SUSENS>

close all, clc
load coreset_errors_0203_1
X = 1:m;
errors = errors{1};

figure(1), hold on

%%
X = X(2:end);
Y = errors(X);
Y = Y./max(Y);
plot(X,Y,'k-','LineWidth',2)

%%
X = X(2:end);
Y = X.*errors(X);
Y = Y./max(Y);
plot(X,Y,'--','LineWidth',2)

%%
X = X(20:end);
Y = X.*log(X).*errors(X);
Y = Y./max(Y);
plot(X,Y,'--','LineWidth',3)

%%
X = X(2:end);
Y = X.^2.*errors(X);
Y = Y./max(Y);
plot(X,Y,':','LineWidth',3)

%%

line([0 2000],[0.365 0.365],'LineStyle','-.','Color',[0.6 0.6 0.6])

axis([0 2000 0 1])
xlabel('Number of iterations N','FontSize',16)
ylabel('f(N)','FontSize',16)
title(sprintf('Coreset error results\nA[%d x %d], k=%d, sparsity=%1.4f',n,d,k,s),'FontSize',16)
legend({'f(N) = eps','f(N) = N eps','f(N) = N logN eps','f(N) = N^2 eps','f(N) = f*(N)+C'},'Location','N','FontSize',20)


% close all, clc
% load coreset_errors_0203_2
% 
% %%
% figure(1), hold on
% 
% for i = 1:length(errors)
%     X = 2:iters(i);
%     Y = errors{i}(X);
%     plot(X,Y,'-','LineWidth',2)
% end
% 
% axis([0 100 0 1])
% title(sprintf('Coreset error results\nA[nx%d], k=%d, sparsity=%1.4f',d,k,s),'FontSize',14)
% xlabel('Number of iterations N','FontSize',16)
% ylabel('eps','FontSize',16)
% legend({'n = 100','n = 200','n = 500','n = 1000','n = 2000','n = 5000'},'Location','NE','FontSize',16)

