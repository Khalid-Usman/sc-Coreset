n = [100 200 500 1000 2000 5000 10000 20000 50000];
d = 1000;
s = 1e-3;
c = 1;
k = d-1;
e = 1e-3;

runtimes = zeros(1,length(n));
fprintf('A[nx%d], s=%1.4f, c=%1.4f, k = %d:\n',d,s,c,k)

for i = 1:length(n)
    
    fprintf('n = %d:\n',n(i))
    
    start_time = tic;
    
    % generate points
    A = SVDCoresetTest.generate_points(n(i),d,s,c);
    
    svds_opts = struct;
    svds_opts.tol = e;
    svds_opts.maxit = n(i);
    [~,~,~] = svds(A,k,'L',svds_opts);

    runtimes(i) = toc(start_time);
    
    fprintf('running time = %fs\n',runtimes(i))
    
end
