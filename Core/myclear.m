diary off
s = dbstatus;
save('dbstatus.mat','s');
clearvars -except mycluster
load('dbstatus.mat')
dbstop(s)
delete('dbstatus.mat')
clear s