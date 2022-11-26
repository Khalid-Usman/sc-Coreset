classdef SVDCoreset < handle
    
    properties
        
        DEBUG = false
        
        max_iter
        max_error
        max_size       
 
        w       % weights
        C       % approx
        
        num_iter
        coreset_size
        coreset_errors
        coreset_idx
        
        approx_error
        
        svds_max_iterations     = 100
        svds_convergence_tol    = 1e-4
        alpha_convergence_tol   = 1e-10
        small_number            = 1e-10
        
        xa,xb   % starting indices of original points (for streaming)
        
    end
    
    %%
    methods
        
        %%
        % input: A -- nxd matrix
        % input: k -- integer 0...n
        function compute(this,A,k)
            
            n = size(A,1);
            d = size(A,2);
            
            assert(isempty(this.max_iter)||isnumeric(this.max_iter)&&(mod(this.max_iter,1)==0)&&(this.max_iter>=0))
            assert(isempty(this.max_error)||(isfloat(this.max_error)&&(this.max_error>=0)&&(this.max_error<1)))
            
            if isempty(this.max_iter) || this.max_iter == 0 || this.max_iter > n
                this.max_iter = n;
            end
            if isempty(this.max_error)
                this.max_error = 0;
            end
            if isempty(this.max_size)
                this.max_size = n;
            end
            
            if k == 0
                
                error('not implemented')
                %this.compute_k0(A)
                
            else
                
                % adjust k for leaf sizes smaller than k
                k = min(k,n);
                
                % compute column-dense kernel
                fprintf('Computing column-dense kernel: d=%d->',d)
                nzx = sum(A,1)>0;
                H = A(:,nzx);
                %H = A;
                fprintf('%d\n',size(H,2));
                
                % compute sparse SVD
                svds_tic = tic;
                svds_opts.tol = this.svds_convergence_tol;
                svds_opts.maxit = this.svds_max_iterations;
                d0 = size(H,2);
                k0 = min(n,d0);
                fprintf('Computing kernel SVD (rank=%d) ... ',k0)
                [U,D,~] = svds(H,k0,'L',svds_opts);
                % [U,D,~] = redsvds(H,k0);
                
                %fprintf('Computing k-rank SVD ... ')
                %[U,S,~] = svds(H,k,'L',svds_opts);
                fprintf('completed in %1.2f min\n',toc(svds_tic)/60)
                
                % compute the residual = sigma_k+1 ... sigma_d
                % sum of squared distances from rows of A to the
                % approximated subspace that is returned from svds
                %sigma = sqrt(abs(norm(A,'fro')^2-norm(D,'fro')^2));
                
                % find z such that the sum of the smallest z squared
                % singular values is less than eps*sigma
                %z = find([0 cumsum(fliplr(diag(D.^2)'))] <= this.max_error*sigma,1,'last')-1;
                
                % TODO: fix definition of sigma
                %R = A(:,k+1:d-z)/sigma;
                %R = bsxfun(@rdivide,R,sqrt(row_inner_products(R)*n));
                %assert(all(full(isnan(R(:))==0)))
                du = size(U,2);
                kx = min(du,k);
                B = D(kx+1:du,kx+1:du);
                P = [U(:,1:kx) U(:,kx+1:du)*B/norm(B,'fro')];
                
                % compute input
                fprintf('Computing input ...\n')
                
                if this.DEBUG
                    fprintf('Type commands to run:\n')
                    keyboard
                else
                    %this.compute_slow(P)
                    this.compute_fast(P)
                end
                
            end
            
            % scale weights
            w = this.w;
            q = norm(diag(sqrt(w))*A,'fro');
            w = norm(A,'fro')*sqrt(w)/q;
            
            % TODO: fix stuff
            w = w./sum(w(:))*numel(w(w>0));
            this.w = w;
            
            % update approximation
            fprintf('Computing approximation ...\n')
            W = diag(this.w);
            C = W*A;
            C(sum(C,2)==0,:) = [];
            this.C = C;
            
            this.approx_error = SVDCoreset.ComputeApproxError(A,C,k);
            
        end
        
        %%
        % input: A -- nxd matrix
        function compute_k0(this,A)
            fprintf('Computing k0 .... \n')
            n = size(X,1);
            d = size(X,2);
            
            % shift to origin
            mu = mean(A,1);
            A = bsxfun(@minus,A,mu);
            assert(all(sum(A,1)<this.small_number))
            
            % normalize
            z = sum(bsxfun(@times,A,A),2);
            A = bsxfun(@rdivide,A,sqrt(z));
            assert(all(sum(sum(bsxfun(@times,A,A),2)))==1)
            
            % setup
            j = 1;
            J = [j];
            c = A(j,:);
            w = sparse(n,1);
            w(j) = 1;
            num_iter = 1;
            m = 1;
            
            fprintf('Computing coreset ... ')
            bstr = '';
            
            % add m-1 additional points to coreset
            for iter = 2:this.max_iter
                
                num_iter = num_iter+1;
                
                % Find the farthest row in X from r.
                D = A*c'; % inner product of points with center
                D(D<this.small_number) = 0;
                j = argmin(D);
                p = A(j,:);
                
                % norm(c-p):
                norm_c = norm(c);
                norm_p = norm(p);
                cp = c*p';
                norm_c_p = sqrt(norm_p^2 + norm_c^2 - 2*cp);
                assert(abs(norm_c_p-norm(c-p))<this.small_number)
                
                % Find the point u on the line from p1 to p2
                % that is closest to the origin.
                v = p - c;
                c1 = p - (v/norm(v))*(v/norm(v)*p');
                
                % norm(c-c1):
                norm_c1 = norm(c1);
                assert(abs(norm_c1-norm(c1))<this.small_number)
                cc1 = c*c1';
                assert(abs(cc1-c*c1')<this.small_number)
                norm_c_c1 = sqrt(norm_c^2 + norm_c1^2 - 2*cc1);
                assert(abs(norm_c_c1-norm(c-c1))<this.small_number)
                
                % update alpha
                alpha = norm(c-c1)/norm(c-p);
                assert(alpha<1)
                
                % update coreset weights
                w = w * (1-abs(alpha));
                w(j) = w(j) + alpha;
                w = w./sum(w);
                
                % check if point was already added
                if not(ismember(j,J))
                    
                    m = m+1;
                    fprintf(' m = %d', m)     
                    % indices of coreset points 1...i
                    J = [J,j];
                    
                end

                % update center
                c = c1;
                
                msg = sprintf('[%d/%d]',iter,this.max_iter);
                msg = [msg sprintf(', alpha = %1.6f\n',alpha)];
                fprintf(msg)
                bstr = repmat(sprintf('\b'),1,length(msg));
                
                % \sum_i z_i ||p_i-c||^2=1+||c||^2
                % so ||c||^2 is the approximation error!
                % if ||c||^2=eps you get (1+eps) approximation
                coreset_error = norm_c^2;
                
                % check alpha convergence
                if abs(alpha) < this.alpha_convergence_tol
                    break
                end
                
                % check error convergence
                if abs(coreset_error-this.max_error)<this.small_number
                    break
                end
                
            end
            
            % update coreset
            this.w = w;
            this.coreset_idx = J;
            this.num_iter = num_iter;
            this.coreset_size = m;
            this.coreset_error = coreset_error;
            
        end
        
        %%
        % input: A -- nxd matrix
        function compute_slow(this,X)
            
            % inefficient reduction to k=0
            XX = row_outer_products(X);
            
            % compute k0 algorithm on XX
            this.compute_k0(XX);
            
            if this.DEBUG
                fprintf('\nidx =\n')
                this.coreset_idx
                dbquit
            end
            
        end
        
        %%
        % input: U -- nxd matrix (non-normalized)
        function compute_fast(this,U)
            fprintf('Computing Fast .... \n')
            size(U)
            n = size(U,1);
            d = size(U,2);
            %assert(abs(norm(U,'fro')^2-size(U,2))<this.small_number)
            
            fprintf('Computing input...\n')
            
            % inner products of non-normalized points
            uu = row_inner_products(U);
            
            % compute the distance from the origin
            % s = ||U||^2 = \sum_j u_i u_i'
            norm_U2 = full(sum(uu));
            norm_U4 = norm_U2^2;
            %z = UU/s;
            
            [~,R] = qr(U);
            norm_UU2 = norm(R*R','fro')^2;
            
            % normalize
            X = bsxfun(@rdivide,U,sqrt(uu));
            %assert(abs(sum(X(:).^2)-n)<this.small_number)
            
            % setup
            j = 1;
            J = [j];
            w = sparse(1,n);
            w(j) = 1;
            Mx = [full((X*X(j,:)').^2)'];
            Mu = [norm(U*X(j,:)')^2];
            Du = w(J)*Mu;
            alpha = 1;
            norm_c2 = 1;
            num_iter = 1;
            coreset_error = norm_c2 + norm_UU2/norm_U4 - 2/norm_U2*Du;
            coreset_errors = [coreset_error];
            
            % compute coreset
            fprintf('Computing coreset ... ')
            bstr = '';
            alphas = [];
            % add N-1 additional points to coreset
            for iter = 2:this.max_iter
                
                num_iter = num_iter+1;
                
                % compute next point from input to coreset:
                % find the farthest row in X from q
                Dx = w(J)*Mx;
                Dx(Dx<this.small_number) = 0;
                j = argmin(Dx);
                
                % check if j is already a coreset point
                % if it is, and there is more than 1 argmin(D)
                % then select a different point
                %     Jmin = find(Dx==Dx(j));
                %     if ismember(j,J) && not(all(ismember(Jmin,J)))
                %        j = Jmin(find(ismember(Jmin,J)==0,1,'first'));
                %     end
                
                % check if point was already added
                if not(ismember(j,J))
                    % maintain inner products of every row in X to the center
                    m = full((X*X(j,:)').^2)';
                    Mx = cat(1,Mx,m);
                    
                    norm_Ux2 = norm(U*X(j,:)')^2;
                    Mu = cat(1,Mu,norm_Ux2);
                    
                    % indices of coreset points 1...m
                    J = cat(2,J,j);
                else
                    % already computed in previous iteration
                    m = Mx((J==j),:);
                    norm_Ux2 = Mu((J==j),:);
                end

                %fprintf('this max size = \n')
                %fprintf(this.max_size)
                %XJ = X(J,:);
                %WJ = diag(w(J));
                
                % calc norm (c-p):
                %norm_c = norm(XJ'*WJ*XJ,'fro')
                %G = sqrt(WJ)*XJ;
                %norm_c = norm(G'*G,'fro');
                %H = G*G';
                %H = sqrt((w(J)*w(J)').*M(1:m,J));
                %norm_c2 = trace(H'*H);
                %norm_c2 = norm_c^2
                %norm_p2 = 1; % norm(x)^2;
                %cp = sum(WJ'*(XJ*x').^2)
                %cp = full(sum((G*x').^2));
                %cp = sum(w(J)'*M(1:m,j));
                %norm_c_p = sqrt(norm_p2 + norm_c2 - 2*cp);
                
                % calc norm(c-c1):
                %comp_pv = (norm_p2/norm_c_p - cp/norm_c_p);
                %norm_c_c1 = norm_c_p - comp_pv;
                
                % calc alpha
                %alpha = norm_c_c1/norm_c_p;
                %assert(alpha<1)
                
                % weighted inner products of coreset points with x
                wm = full(dot(w,m)); % \sum_i w_i(x_i x')^2
                
                % calc alpha in O(1)
                %alpha = (1-wm)/(1+norm_c2-2*wm);
                alpha = (1 - wm - norm_Ux2/norm_U2 + Du/norm_U2) / (1 + norm_c2 - 2*wm);
                alpha(alpha>1) = 1;
                
                % update new norm(c)
                norm_c1_2 = (1-alpha)^2 + alpha^2*norm_c2 + 2*(1-alpha)*alpha*wm;
                assert(norm_c1_2>=0)
                
                % update coreset weights
                %w = w * alpha;
                %w(j) = w(j) + (1-alpha);
                %w = w./sum(w);
                
                % update coreset weights
                ej = sparse(1,j,1,1,n);
                w1 = (1-alpha)*ej + alpha*w;
                assert(abs(1-sum(w1))<this.small_number)
                
                Du1 = alpha*Du + (w1(j)-w(j)*alpha)*norm_Ux2;
                
                % \sum_i z_i ||p_i-c||^2=1+||c||^2
                % so ||c||^2 is the approximation error!
                % if ||c||^2=eps you get (1+eps) approximation
                % run algorithm around center I/d instead of 0
                % d = \sum_j ||u_i u_i'||
                % approximation error is ||c-I/d||^2 = ||c||^2-1/s
                %coreset_error = norm_c1_2 - 1/norm_U2;
                coreset_error = norm_c1_2 + norm_UU2/norm_U4 - 2/norm_U2*Du;
                coreset_error(coreset_error<0) = 0;
                
                coreset_errors(end+1) = coreset_error;
                
                norm_c2 = norm_c1_2;
                w = w1;
                Du = Du1;
                
                msg = sprintf('[%d/%d]',iter,this.max_iter);
                msg = [msg sprintf(', alpha = %1.2f, error = %1.6f',alpha,coreset_error)];
                %fprintf([bstr msg])
                bstr = repmat(sprintf('\b'),1,length(msg));
                
                % check error convergence
                %if coreset_error <= this.max_error
                %    break
                %end
                %this.max_size
                %length(J)              
                if this.max_size == length(J)
                    break
                end

                %                 % check error convergence
                %                 if abs(1-alpha) <= this.alpha_convergence_tol
                %                     break
                %                 end
                
            end
            
            coreset_size = length(J);
            %assert(coreset_size>1)
            
            msg = sprintf('[%d/%d]',iter,this.max_iter);
            msg = [msg sprintf(', alpha = %1.6f, error = %1.6f\n',alpha,coreset_error)];
            fprintf([bstr msg])
            
            % update coreset
            this.w = w;
            this.coreset_idx = J;
            this.num_iter = num_iter;
            this.coreset_size = coreset_size;
            this.coreset_errors = coreset_errors;
            
            if this.DEBUG
                fprintf('\nidx =\n')
                this.coreset_idx
                dbquit
            end
            
        end

        %%
        % input: A -- nxd matrix
        % input: k -- integer 0...n
        function compute_new(this,A,k)

            fprintf('Computing New ..... \n')
            n = size(A,1);
            d = size(A,2);
 
            err = this.max_error;

            Xu = speye(d) * k;
            Xl = speye(d) * -k;
            delta_u = err + 2*err^2;
            delta_l = err - 2*err^2;

            w = sparse(1,n);
            Z = sparse(d,d);

            for iter = 1:this.max_iter

                Mu = inv((Xu + delta_u*A'*A) - Z);
                Ml = inv(Z - (Xl + delta_l*A'*A));

                beta_diff = Inf;
                beta_uj_min = 0;
                j = 0;
                for i = 1:n

                    ai = A(i,:)';
                    beta_li = (ai'*Ml*A'*A*Ml*ai)/(delta_l*trace(A*Ml*A'*A*Ml*A')) - ai'*Ml*ai;
                    beta_ui = (ai'*Mu*A'*A*Mu*ai)/(delta_u*trace(A*Mu*A'*A*Mu*A')) - ai'*Mu*ai;

                    if (beta_li - beta_ui) < beta_diff
                        beta_diff = (beta_li - beta_ui);
                        beta_uj_min = beta_ui;
                        j = i;
                    end

                end

                w(j) = 1/beta_uj_min;
                aj = A(j,:)';
                Z = Z + w(j)^2*aj*aj';

            end 

            this.w = w; 

        end
        
    end
    
    %%
    methods (Static)
        
        % calculate error
        % input: A -- nxd matrix
        % input: V -- svd of A (A = UDV')
        % input: C -- nxd matrix (coreset C = W*A)
        % input: k -- subspace dimension
        function err = ComputeApproxError(A,C,k)
            
            svds_opts = struct;
            svds_opts.tol = 1e-4;
            svds_opts.maxit = 100;
            
            % [~,~,Q] = redsvds(C,k);
            [~,~,Q] = svds(C,k,'L',svds_opts);
            [~,~,V] = svds(A,k,'L',svds_opts);
            
            %approx_dist = full(sum(row_inner_products(X'*X-C'*C)));
            %best_dist = full(sum(row_inner_products(X'*X)));
            %err = abs(approx_dist-best_dist)/n;
            %err(err<0) = 0;
            
            %norm_A_AQ = norm(A,'fro').^2 - norm(A*Q,'fro').^2; % approx dist
            %norm_A_AV = norm(A,'fro').^2 - norm(A*V,'fro').^2; % opt dist
            
            % using norm(X,'fro').^2 = sum(row_inner_products(X))
            approx_dist = full(sum(row_inner_products(A)) - sum(row_inner_products(A*Q)));
            opt_dist = full(sum(row_inner_products(A)) - sum(row_inner_products(A*V)));
            fprintf('printing distances...')
            approx_dist
            opt_dist 
            %err = (approx_dist-opt_dist)/opt_dist;
            err = (approx_dist/opt_dist) - 1 % same as above
            
            % scale by n to get average error per point
            n = size(A,1);
            err = err/n;
            
        end
        
    end
    
    
end

function fprintf(varargin)
global CORESET_OUTPUT
if isempty(CORESET_OUTPUT)
    CORESET_OUTPUT = 1;
end
if CORESET_OUTPUT
    builtin('fprintf',varargin{:});
end
end

% ------------------------------------------------
% reformatted with stylefix.py on 2014/11/21 14:44
