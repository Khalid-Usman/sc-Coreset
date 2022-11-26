global CORESET_OUTPUT STREAM_OUTPUT
CORESET_OUTPUT = 1
STREAM_OUTPUT = 0

disp(repmat('=',1,100))

k = 1000;
cores = 10;

A = dlmread('/home/khalid/Datasets/neurons_HVG1000.txt', ',');
A = sparse(A);
n = size(A,1);
nn = ceil(n/cores)
size(A)

fprintf('Computing SVD ...\n')
svds_opts = struct;
svds_opts.tol = 1e-4;
svds_opts.maxit = nn;
%[~,~,V] = svds(A,k,'L',svds_opts);

%M = [1000 2000 3000 4000 5000];
M = [500]
%M = [26639 53279 79919 106558 133198];
disp(repmat('_',1,80))

for i = 1:length(M)
    coresetSize = M(i)
    %coresetSize = coresetSize/cores;
    maxIter = ceil(n/cores);

    B{1} = sparse(A(1:floor(n/10),:));
    B{2} = sparse(A(floor(n/10)+1:floor(n/5),:));
    B{3} = sparse(A(floor(n/5)+1:floor((3*n)/10),:));
    B{4} = sparse(A(floor((3*n)/10)+1:floor(2*n/5),:));
    B{5} = sparse(A(floor(2*n/5)+1:floor(n/2),:));
    %B{6} = sparse(A(floor(n/2)+1:floor((3*n)/5),:));
    %B{7} = sparse(A(floor((3*n)/5)+1:floor((7*n)/10),:));
    %B{8} = sparse(A(floor((7*n)/10)+1:floor((4*n)/5),:));
    %B{9} = sparse(A(floor((4*n)/5)+1:floor((9*n)/10),:));
    %B{10} = sparse(A(floor((9*n)/10)+1:n,:));
    clear A;

    delete(gcp('nocreate'))
    parpool(cores,'IdleTimeout', 300)
    parfor i=1:5
        i
        %coreset = computeCoreset(B{1,i}, k, coresetSize, maxIter)
        fileName = strcat('/home/khalid/Datasets/neurons_HVG1000_coreset_parallel_beforePCA_', int2str(coresetSize), '_', int2str(i), '.txt');
        computeCoreset(B{1,i}, k, coresetSize, maxIter, fileName)
        %parsave(fileName, coreset)
    end
end
