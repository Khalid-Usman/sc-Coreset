classdef SVDCoresetTest
    
    methods(Static)
        
        % n = number of points
        % d = dimension R^d
        % r = density of A, i.e. P(A_{i,j} > 0)
        % s = sparsity of A, i.e. P(sum_{i=1}^n A_{i,j} > 0) 
        %   = (max nnz(A))/d
        % the kernel density will be r/s <= 1, so r < s
        function A = gen_points(n,d,r,s)
            
            % generate column kernel with density r0 = r/s
            d0 = ceil(d*s);
            r0 = r/s;
            H = sprand(n,d0,r0);
            
            % randomly select columns with sparsity s = d0/d
            idx = randperm(d);
            nzx = sort(idx(1:d0));
            
            % generate sparse matrix A with kernel entries
            A = sparse(n,d);
            A(:,nzx) = H;
            
            % ensure that each point has at least one entry
            zrx = find(sum(H,2)==0);
            zcx = nzx(randi(length(nzx),1,length(zrx)))';
            zix = sub2ind(size(A),zrx,zcx);
            A(zix) = rand(1,length(zix));
            
        end
        
        % n = number of points
        % d = dimension R^d
        % j = subspace dimension (s = k/d)
        % rd = density of A, i.e. P(A_{i,j} > 0)
        % rj = density of H, i.e. H(A_{i,j} > 0)
        % output: A, a random matrix of size n*d
        % where all the points are lying on a j-dimensional subspace,
        % with additional Gaussian noise of variance 1
        function A = gen_noisy_subspace(n,d,j,rd,rj)
            
            % generate the k-dimensional subspace 
            H = sprand(n,j,rj);
            
            % generate sparse matrix A with density rd
            % with additional Gaussian noise of variance 1
            G = sprand(n,d,rd);
            A = sparse(n,d)+G;
            
            % randomly select columns with sparsity s = d0/d
            idx = randperm(d);
            nzx = sort(idx(1:j));
            
            % generate sparse matrix A with kernel entries
            A(:,nzx) = H;
            
            % ensure that each point has at least one entry
            zrx = find(sum(H,2)==0);
            zcx = nzx(randi(length(nzx),1,length(zrx)))';
            zix = sub2ind(size(A),zrx,zcx);
            A(zix) = rand(1,length(zix));
            
            % ensure that (max nnz(A)) = s
            nnz = sum(A~=0,2);
            zgx = (nnz>j);
            zgn = nnz(zgx)-j;

            % for every row of A, remove zgn(i) elements from columns in A(i,zx)
            zx = setdiff(1:d,nzx);
            for i = 1:n
               if zgx(i)
                  % non-zero non-kernel cols
                  zni = zx(A(i,zx) > 0);
                  % number of elements to remove
                  nri = zgn(find(zgx==1)==i);
                  xxi = randperm(length(zni));
                  % indices of elements to remove
                  zxi = zni(xxi(1:nri));
                  A(i,zxi) = 0; 
               end
            end
            
        end
        
        %%
        % A = nxd input points 
        % k = number of largest singular values
        % optional: 'MaxIter' [integer] -- max number of iterations
        % optional: 'MaxError' [float] -- max error
        % optional: 'Stream' [integer] -- use streaming with this leaf size
        % optional: 'Debug' [boolean] -- enable keyboard
        function [runtime,num_iter,coreset_size,coreset_errors] = test(A,k,varargin)
            
            profile off, profile on
            
            n = size(A,1);
            d = size(A,2);

            parser = inputParser;
            parser.addParameter('MaxIter',0,@(x)(isnumeric(x)&&(mod(x,1)==0)&&(x>=0)&&(x<=n)))
            parser.addParameter('MaxError',0,@(x)(isfloat(x)&&(x>=0)&&(x<1)))
            parser.addParameter('Stream',0,@(x)(isnumeric(x)&&(x>=0)&&(mod(x,1)==0)))
            parser.addParameter('Debug',0,@boolean)
            parser.parse(varargin{:})
            max_iter = parser.Results.MaxIter;
            max_error = parser.Results.MaxError;
            leaf_size = parser.Results.Stream;
            is_debug = parser.Results.Debug;
            
            fprintf('A[%dx%d], k=%d, m=%d, eps=%1.6f:\n',n,d,k,max_iter,eps)
            
            start_time = tic;
            
            if leaf_size == 0
                
                % direct computation
                D = SVDCoreset;
                D.DEBUG = is_debug;
                D.max_iter = max_iter;
                D.max_error = max_error;
                D.compute(A,k);
                
            else
                
                stream = SVDBufferStream;
                stream.k = k;
                stream.max_iter = max_iter;
                stream.max_error = max_error;
                stream.leaf_size = leaf_size;
                stream.addPoints(A)
                D = stream.getUnifiedCoreset();
                
            end
            
            num_iter = D.num_iter;
            coreset_size = D.coreset_size;
            coreset_errors = D.coreset_errors;
            
            fprintf('num iter = %d, coreset size = %d, error = %1.6f\n',num_iter,coreset_size,coreset_errors(end))
            
            % compute svds on the coreset for timing
            fprintf('Computing k-rank SVD on coreset ... ')
            svds_tic = tic;
            svds_opts = struct;
            svds_opts.tol = D.svds_convergence_tol;
            svds_opts.maxit = D.svds_max_iterations;
            [~,~,~] = svds(D.C,min(k,D.coreset_size),'L',svds_opts);
            fprintf('completed in %1.2fs\n',toc(svds_tic))
            
            runtime = toc(start_time);
            
            fprintf('Test running time = %.2fs\n',runtime)
            fprintf('Done!\n')
                        
        end
        
        %%
        function r = CalcDensity(A)
            n = size(A,1);
            d = size(A,2);
            r = numel(nonzeros(A))/(n*d);
        end
        
        function s = CalcSparsity(A)
            d = size(A,2);
            s = full(max(sum(A~=0,2))/d);
        end
        
        %% test for number of iterations
        function [iters,sizes,errors,runtimes] = pivot_test(n,d,s,c,k,m,e,l)
            
            fprintf('Pivot test:\n')
            
            arg_keys = {'n','d','s','c','k','m','e','l'};
            args = cell(size(arg_keys));
            for i = 1:length(args)
                args{i} = eval(arg_keys{i});
            end
            pivot_index = argmax(cellfun(@length,args));
            pivot_length = length(args{pivot_index});
            for i = 1:length(args)
                if numel(args{i}) == 1
                    args{i} = repmat(args{i},1,pivot_length);
                elseif numel(args{i}) ~= pivot_length
                    error([],'Invalid size of arg %d!',i)
                end
                eval([arg_keys{i} '=args{i};'])
                fprintf('%s = %s\n',arg_keys{i},mat2str(args{i}))
            end
            
            runtimes = zeros(1,pivot_length);            
            iters = zeros(1,pivot_length);
            sizes = zeros(1,pivot_length);
            errors = cell(1,pivot_length);
            
            for i = 1:pivot_length
                
                fprintf('\n%s\n',repmat('-',1,80))
                fprintf('#%d: ',i);
                
                % generate points
                A = SVDCoresetTest.generate_points(n(i),d(i),s(i),c(i));
                [runtimes(i),iters(i),sizes(i),errors{i}] = SVDCoresetTest.test(A,k(i),'MaxIter',m(i),'MaxError',e(i),'Stream',l(i));
                
            end
            
        end
        
    end
    
end
