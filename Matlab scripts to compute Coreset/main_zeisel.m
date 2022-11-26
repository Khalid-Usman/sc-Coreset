global CORESET_OUTPUT STREAM_OUTPUT
CORESET_OUTPUT = 1
STREAM_OUTPUT = 0

clear;
clc;

disp(repmat('=',1,100))

if exist('AWS','var') && AWS
    % put in some code to load mat files
    % from S3 or wherever you saved them
else
%     myclear, clc
%     load A1000
end

A = dlmread('/home/khalid/Datasets/zeisel_fbpca40.txt', ',');
%A = dlmread('/home/khalid/Datasets/zeisel_seurat_preprocessed_withoutHeader.txt', ',');
%A = transpose(A);
A = sparse(A);
n = size(A,1);
size(A)

k = 40;

fprintf('Computing SVD ...\n')
svds_opts = struct;
svds_opts.tol = 1e-4;
svds_opts.maxit = n;

%[~,~,V] = svds(A,k,'L',svds_opts);
M = [500, 1000];
%M = [1371 2743 4114 5486 6857];
disp(repmat('_',1,80))

for i = 1:length(M)

    m = M(i)
    max_error = 0;
    epsilon = 1;
    max_iter = n;               %ceil(m^2/epsilon);
    leaf_size = max(1,n)     %max(100,m);     %ceil(max_iter/2);

    D = SVDCoreset;
    D.max_iter = max_iter;
    D.max_size = m;
    D.max_error = max_error;
    D.compute(A,m);

    size(D.C)
    fprintf('preparing to save \n')
    X = full(D.C);
    Idx = full(D.coreset_idx);
    fileName = strcat('/home/khalid/Datasets/zeisel_coreset_40_', num2str(m), '.txt')
    idxName = strcat('/home/khalid/Datasets/zeisel_coreset_40_idx_', num2str(m), '.txt')
    dlmwrite(fileName, X, 'delimiter', ',');
    dlmwrite(idxName, Idx);
end
