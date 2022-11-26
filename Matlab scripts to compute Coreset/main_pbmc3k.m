%%%%    Data    %%%%
data = dlmread('/home/khalid/Datasets/pbmc3k_pca40.txt', ',');
%data = dlmread('/home/khalid/Datasets/pbmc3k_seurat_preprocessed_withoutHeaders.txt', ',');
%data = transpose(data);
size(data)
data = sparse(data);
n = size(data,1);
%%%%    Variables    %%%%
k = 40;
max_error = 0;
max_iter = n;

%%%%%    Compute Coreset    %%%%%%%
M = [500, 1000];
for i = 1:length(M)
    m = M(i)
    C = SVDCoreset;
    C.max_iter = max_iter;
    C.max_error = max_error;
    C.max_size = m;
    C.compute(data,k);
    size(C.C)
    X = full(C.C);
    Idx = full(C.coreset_idx);
    fileName = strcat('/home/khalid/Datasets/pbmc3k_coreset_40_', num2str(m), '.txt');
    dlmwrite(fileName, X, 'delimiter', ',');
    idxName = strcat('/home/khalid/Datasets/pbmc3k_coreset_40_', num2str(m), '_idx.txt')
    dlmwrite(idxName, Idx);
end
fprintf('Done')
