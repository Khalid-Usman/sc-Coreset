%%%%    Data    %%%%
A = dlmread('/home/khalid/Datasets/mca_fbpca75.txt', ',');
size(A);
A = sparse(A);
n = size(A,1);
n = floor(n/10)
A = A(1:n,:);
size(A)
n = size(A,1);
%%%%    Variables    %%%%
k = 75;
max_error = 0;
max_iter = n;

%%%%%    Compute Coreset    %%%%%%%
M = [1000 2000 3000 4000 5000];
for i = 1:length(M)
    m = M(i)
    start_time = tic;
    C = SVDCoreset;
    C.max_iter = max_iter;
    C.max_error = max_error;
    C.max_size = m;
    C.compute(A,k);
    size(C.C)
    running_time = toc(start_time)
end
fprintf('Done')
