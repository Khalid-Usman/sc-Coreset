clf, hold on, grid on
plot(0,0,'og')
ang = 0:0.01:2*pi;
xp = cos(ang);
yp = sin(ang);
plot(xp,yp,'k-','LineWidth',1);
axis equal