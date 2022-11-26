%#ok<*NOPTS,*SNASGU>
global CORESET_OUTPUT STREAM_OUTPUT
CORESET_OUTPUT = 1
STREAM_OUTPUT = 0

if not(exist('AWS','var'))
    block_no = 1;
    prefix = 'small_synth_bow_';
    block = load([prefix num2str(block_no)],'A');
end

%%
A = block.A;
clear block

r1 = find(sum(A,2)>0,1,'first')
r2 = find(sum(A,2)>0,1,'last')-1
A = A(r1:r2,:);
n = size(A,1);
d = size(A,2);
% s = length(A(A>0))/length(A(:));
s = numel(nonzeros(A))/(n*d);
k = 100;
max_error = 1e-6;
max_iter = 2*k/0.5;
leaf_size = max_iter/2;

fprintf([repmat('-',1,80) '\n'])
fprintf('block=%d: A[%dx%d], rows %d:%d, s=%1.6f, k=%d\n',block_no,n,d,r1,r2,s,k)
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
coreset_size = D.coreset_size;
coreset_errors = D.coreset_errors;

fprintf('num iter = %d, coreset size = %d, error = %1.6f\n',num_iter,coreset_size,coreset_errors(end))

%%

approx_err = SVDCoreset.compute_error(A,D.C,k);

%%

running_time = toc(start_time);

fprintf('Done!\n')
fprintf('Total running time = %f min\n',running_time/60)
