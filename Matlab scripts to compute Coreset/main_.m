%%%    Data    %%%%
%data = dlmread('/home/khalid/Datasets/pbmc3k_pca50.txt', ',');
%data = dlmread('/home/khalid/Datasets/moca_fbpca100.txt', ',');
%size(data)
%n = size(data,1)

for j=1:10
    A = dlmread('/home/khalid/Datasets/neurons_fbpca100.txt', ',');
    size(A)
    n = size(A,1) 
    A = A(floor(((j-1)*n)/10)+1:floor((j*n)/10),:);
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
        fileName = strcat('/home/khalid/Datasets/neurons_coreset_', num2str(m), '_', num2str(j), '.txt');
        dlmwrite(fileName, X, 'delimiter', ',');
        idxName = strcat('/home/khalid/Datasets/neurons_coreset_', num2str(m), '_idx_', num2str(j), '.txt')
        dlmwrite(idxName, Idx);
    end
end
fprintf('Done')
