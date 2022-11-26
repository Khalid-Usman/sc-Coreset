clear;

data = dlmread('/home/khalid/Datasets/neurons_fbpca100.txt', ' ');
%data = dlmread('/home/khalid/Datasets/moca_fbpca100.txt', ' ');
%data = dlmread('/home/khalid/Datasets/mca_fbpca75.txt', ',');
%data = dlmread('/home/khalid/Datasets/pbmc68k_fbpca50.txt', ' ');
size(data)
data = sparse(data);
n = size(data,1);
data = data(1:floor(n/20),:);
size(data)
nn = size(data,1);

%%%%    Variables    %%%%
k = 100;
max_error = 0;
max_iter = nn;

%%%%%    Compute Coreset    %%%%%%%
M = [500,1000, 1500, 2000, 2500];
for i = 1:length(M)
    m = M(i)
    start_time = tic;
    
    C = SVDCoreset;
    C.max_iter = max_iter;
    C.max_error = max_error;
    C.max_size = m;
    C.compute(data,k);
    size(C.C)

    running_time = toc(start_time);
    fprintf('Done!\n')
    fprintf('Total running time = %f sec\n',running_time')
end
fprintf('Done')
