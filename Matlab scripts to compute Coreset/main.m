%%%    Data    %%%%
%data = dlmread('/home/khalid/Datasets/pbmc3k_pca50.txt', ',');
A = dlmread('/home/khalid/Datasets/moca_fbpca100.txt', ',');
size(A)
n = size(A,1)
A = A(floor((9*n)/10)+1:floor((10*n)/10),:);
A = sparse(A);
nn = size(A,1);

%%%%    Variables    %%%%
k = 100;
max_error = 0;
max_iter = n;

%%%%%    Compute Coreset    %%%%%%%
M = [500];
for i = 1:length(M)
    m = M(i)
    C = SVDCoreset;
    C.max_iter = max_iter;
    C.max_error = max_error;
    C.max_size = m;
    C.compute(A,k);
    size(C.C)
    X = full(C.C);
    Idx = full(C.coreset_idx);
    fileName = strcat('/home/khalid/Datasets/neurons_coreset_', num2str(m), '_10.txt');
    dlmwrite(fileName, X, 'delimiter', ',');
    idxName = strcat('/home/khalid/Datasets/neurons_coreset_', num2str(m), '_idx_10.txt')
    dlmwrite(idxName, Idx);
end
fprintf('Done')
