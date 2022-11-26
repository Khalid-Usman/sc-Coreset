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

%A = dlmread('/home/khalid/Datasets/mca_fbpca75.txt', ',');
A = dlmread('/home/khalid/Datasets/MCA_HVG_2k.txt', ',');
size(A)
A = transpose(A);
A = sparse(A);
n = size(A,1);
nn = ceil(n/8)
size(A)

k = 1000;

fprintf('Computing SVD ...\n')
svds_opts = struct;
svds_opts.tol = 1e-4;
svds_opts.maxit = n;

%[~,~,V] = svds(A,k,'L',svds_opts);
M = [5000];
%M = [4850 9701 14551 19402 24253];
disp(repmat('_',1,80))

for i = 1:length(M)

   %fprintf('Test #%d/%d:\n',i,length(M))

    m = M(i)
    max_error = 0;
    epsilon = 1;
    max_iter = nn;               %ceil(m^2/epsilon);
    leaf_size = max(1,nn)     %max(100,m);     %ceil(max_iter/2);

    %fprintf('block=%d: A[%dx%d], rows %d:%d, r=%1.6f, s=%1.6f\n',block_no,n,d,x1,x2,r,s)
    %fprintf('k=%d, eps=%1.4f, N=%d, L=%d\n',k,max_error,max_iter,leaf_size)

    stream = SVDBufferStream;
    stream.k = k;
    stream.max_iter = max_iter;
    stream.max_size = m;
    stream.max_error = max_error;
    stream.leaf_size = leaf_size;
    stream 

    start_time = tic;
%
    fprintf('Computing SVD coreset ...\n')
    stream.addPoints(A)
    D = stream.getUnifiedCoreset();

    %size(D.C)
    fprintf('preparing to save \n')
    X = full(D.C);
    Idx = full(D.coreset_idx);
    fileName = strcat('/home/khalid/Datasets/mca_beforePCA_coreset_', num2str(m), '.txt')
    idxName = strcat('/home/khalid/Datasets/mca_beforePCA_coreset_idx_', num2str(m), '.txt')
    dlmwrite(fileName, X, 'delimiter', ',');
    dlmwrite(idxName, Idx);
    %save('/home/khalid/Datasets/clustering/mca_svdCoreset_w_10k.mat', 'w');

%    num_iter = D.num_iter;
%    coreset_size = D.coreset_size;
%    coreset_error = D.coreset_errors(end);
%    coreset_errors(i) = coreset_error;
    %fprintf('num iter = %d, coreset size = %d, error = %1.6f\n',num_iter,coreset_size,coreset_errors(end))
%    coreset_sizes(i) = coreset_size;
 %   D_approx_err = SVDCoreset.ComputeApproxError(A,D.C,k);
 %   D_approx_errs(i) = D_approx_err;

%%
%        fprintf('Computing uniform random sampling ...\n')
%        Ru = UniformRandomSampling(A,V,k,coreset_size);
%        Ru_approx_err = Ru.approx_error
%        Ru_approx_errs(i) = Ru_approx_err;        
%        for ti = 1:num_random_sampling_trials
%            fprintf('.')
%            Ru = UniformRandomSampling(A,V,k,coreset_size);
%            errs(ti) = Ru.approx_error;
%        end
%        fprintf('\n');
%        Ru_approx_err = mean(errs);
%        Ru_approx_std = sqrt(var(errs));
%        Ru_approx_errs(i) = Ru_approx_err;
%        Ru_approx_stds(i) = Ru_approx_std;
%     
%        fprintf('Computing weighted random sampling ...\n')
%        Rw = WeightedRandomSampling(A,V,k,coreset_size);
%        Rw_approx_err = Rw.approx_error
%        Rw_approx_errs(i) = Rw_approx_err;
%        for ti = 1:num_random_sampling_trials
%            fprintf('.')
%            Rw = WeightedRandomSampling(A,V,k,coreset_size);
%            errs(ti) = Rw.approx_error;
%        end
%        fprintf('\n');
%        Rw_approx_err = mean(errs);
%        Rw_approx_std = sqrt(var(errs));
%        Rw_approx_errs(i) = Rw_approx_err;
%        Rw_approx_stds(i) = Rw_approx_std;
        
%
        fprintf('Computing SC coreset ...\n')

%         stream.addPoints(B)
%         SC = stream.getUnifiedCoreset();
%         SC_approx_err = SVDCoreset.ComputeApproxError(B,V1,SC.C,k)
%         SC_approx_errs(i) = SC_approx_err;

        %%

        running_time = toc(start_time);

        fprintf('Done!\n')
        fprintf('Total running time = %f min\n',running_time/60)

        disp(repmat('-',1,60))

end
clear V1
