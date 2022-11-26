%%%%    Data    %%%%
data = dlmread('/home/khalid/Datasets/pbmc68k_fbpca50.txt', ',');
%data = dlmread('/home/khalid/Datasets/pbmc68k_1000_withoutHeaders.txt', ',');
%size(data)
data = sparse(data);
n = size(data,1);
%data = data(1:floor(n/20),:);
size(data)
nn = size(data,1);

%%%%    Variables    %%%%
k = 50; %1000;
max_error = 0;
max_iter = nn;

%%%%%    Compute Coreset    %%%%%%%
M = [5000];
%M = [10, 20, 30, 40, 50];
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
    fileName = strcat('/home/khalid/Datasets/pbmc68k_coreset_', num2str(m), '.txt');
    dlmwrite(fileName, X, 'delimiter', ',');
    idxName = strcat('/home/khalid/Datasets/pbmc68k_coreset_', num2str(m), '_idx.txt')
    dlmwrite(idxName, Idx);
end
fprintf('Done')
