classdef UniformRandomSampling < handle
    
    properties
        
        w       % weights
        C       % approx
        
        coreset_idx
        coreset_size
        approx_error
        
    end
    
    %%
    methods
        
        % input: A -- nxd matrix
        % input: V -- svd of A (A = UDV')
        % input: k -- subspace dimension
        % input: m -- coreset size (integer 1...n)
        function this = UniformRandomSampling(A,V,k,m)
            
            n = size(A,1);
            d = size(A,2);
            
            sample_idx = this.sample(n,m);
            
            w = sparse(1,n);
            w(sample_idx) = 1;
            
            % scale weights
            q = norm(diag(sqrt(w))*A,'fro');
            w = norm(A,'fro')*sqrt(w)/q;

            this.w = w;
            this.coreset_idx = sample_idx';
            
            % update approximation
            W = diag(this.w);
            C = W*A;
            C(sum(C,2)==0,:) = [];
            this.C = C;
            
            this.coreset_size = size(C,1);
            this.approx_error = SVDCoreset.ComputeApproxError(A,V,C,k);
            
        end
        
        % Sample indices uniformly at random
        function sample_idx = sample(this,n,sample_size)
            
            rand_idx = randperm(n);
            sample_idx = rand_idx(1:sample_size);
            
        end
        
    end
    
end

