myclear, clc
profile off; profile on

mycluster = 'mikhail_c3x32_1';

AWS = true;

% load A1000
% load A4685

fprintf('Starting job ...\n')
batch(parcluster(mycluster),'test_wiki');

fprintf('Done!\n')

