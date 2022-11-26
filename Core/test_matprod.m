myclear, close all, clc
global CORESET_OUTPUT STREAM_OUTPUT
CORESET_OUTPUT = 0
STREAM_OUTPUT = 0

disp(repmat('=',1,100))

% if not(exist('AWS','var'))
%     myclear, close all, clc
%     block_no = 1;
%     prefix = 'synth_bow_';
%     block = load([prefix num2str(block_no)],'A');
% end

N = [100 200 500 1000 2000 5000 10000 20000 50000 100000]
d = 10000
rd = 1e-6

j = 1000
rj = 1

% block_no = 1;
% A = SVDCoresetTest.gen_noisy_subspace(n,d,j,rd,rj);
% x1 = 1;
% x2 = n;
%
% rA = SVDCoresetTest.CalcDensity(A)
% sA = SVDCoresetTest.CalcSparsity(A)
%
%
% A = block.A;
% clear block
%
% x1 = find(sum(A,2)>0,1,'first');
% x2 = find(sum(A,2)>0,1,'last')-1;
% A = A(x1:x2,:);
% n = size(A,1);
% d = size(A,2);
% % s = length(A(A>0))/length(A(:));
% s = numel(nonzeros(A))/(n*d);

k = 100
m = 100
disp(repmat('_',1,80))

for i = 1:length(N)
    n = N(i)
    
    block_no = 1;
    A = SVDCoresetTest.gen_noisy_subspace(n,d,j,rd,rj);
    x1 = 1;
    x2 = n;
    
    rA = SVDCoresetTest.CalcDensity(A)
    sA = SVDCoresetTest.CalcSparsity(A)
    
    %     svds_opts.tol = 1e-4;
    %     svds_opts.maxit = 100;
    %     [~,~,V] = svds(A,k,'L',svds_opts);
    
    %%
    fprintf('Computing SVD coreset ...\n')
    
    max_error = 0.01;
    epsilon = 1;
    max_iter = m/epsilon
    leaf_size = ceil(max_iter/2);
    
    fprintf('block=%d: A[%dx%d], rows %d:%d, r=%1.6f, k=%d\n',block_no,n,d,x1,x2,rd,k)
    fprintf('eps=%1.4f, N=%d, L=%d\n',max_error,max_iter,leaf_size)
    
    stream = SVDBufferStream;
    stream.k = k;
    stream.max_iter = max_iter;
    stream.max_error = max_error;
    stream.leaf_size = leaf_size;
    
    %%
    fprintf('Streaming points ...\n')
    
    start_time = tic;
    
    stream.addPoints(A)
    D = stream.getUnifiedCoreset();
    
    num_iter = D.num_iter;
    coreset_size = D.coreset_size
    coreset_errors = D.coreset_errors;
    
    %     fprintf('num iter = %d, coreset size = %d, error = %1.6f\n',num_iter,coreset_size,coreset_errors(end))
    %     coreset_sizes(i) = coreset_size;
    %     D_approx_err = SVDCoreset.ComputeApproxError(A,V,D.C,k)
    %     D_approx_errs(i) = D_approx_err;
    
    running_time = toc(start_time);
    
    fprintf('Done!\n')
    fprintf('Total running time = %f min\n',running_time/60)
    
    disp(repmat('.',1,60))
    
    %%
    fprintf('Computing matrix product approx...\n')
    
    start_time = tic;
    
    M = SVDCoreset;
    M.max_iter = n;
    M.max_error = 0.01;
    M.compute_new(A,k);
    
    %     fprintf('num iter = %d, coreset size = %d, error = %1.6f\n',num_iter,coreset_size,coreset_errors(end))
    %     coreset_sizes(i) = coreset_size;
    %     M_approx_err = SVDCoreset.ComputeApproxError(A,V,M.C,k)
    %     M_approx_errs(i) = M_approx_err;
    
    running_time = toc(start_time);
    
    fprintf('Done!\n')
    fprintf('Total running time = %f min\n',running_time/60)
    
    disp(repmat('-',1,80))
    
end
