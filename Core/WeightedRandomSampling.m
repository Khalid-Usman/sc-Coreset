classdef WeightedRandomSampling < handle
    
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
        function this = WeightedRandomSampling(A,V,k,m)
            
            n = size(A,1);
            d = size(A,2);
            
            % find the squared row norms of U
            % and use as the sampling distribution 
            [U,~,~] = svds(A);
            dists = row_inner_products(U);
            [sample_weights,sample_idx] = this.sample(dists,m);
            
            % index back into full set of points
            w = sparse(1,n);
            w(sample_idx) = sample_weights;
            
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

        % Sample indexes with probability proportional to their dist values,
        %  and assign weights proportional to the inverse of this prob.
        % Input:
        %       dists - a vector of n numbers, usually represent distance of points from
        %               a shape (usually opt center)
        %       sampleSize -
        % Output:
        %     indexes - a vector of size sampleSize, contains the index of the chosen
        %               points (corresponding to dists)
        %     weights - a vector of sampleSize weights, where weight[i] is the
        %               weight of the i'th point
        function [sample_weights,sample_idx] = sample(this,dists,sample_size)

            sumDists=sum(dists);
            dists(dists < 0) = 0;
            
            % Represent the probabilities as consecutive sub-intervals on
            % the interval [0,1], such that the length of the i'th sub-interval
            % is proportional to the i'th probability.
            % Then sample a point uniforly at random from [0,1], and choose
            % the i'th index that corresponds to the sub-interval that contains
            % this point.
            % Repeat sampleSize times, and denote by bins[i] the number of times that the
            % i'th sub-interval was selected, for every 1<i<sampleSize.
            ind = randsample(length(dists),sample_size,true,full(dists));
            
            % Remove the duplicated indexes
            ind = sort(ind); % can use count sort in O(n) time.
            h = histc(ind, ind);
            sample_idx = ind(h > 0);
            copies = h(h>0); % how many times each index appeared
            
            sample_weights = (sumDists/sample_size)./dists(sample_idx);
            sample_weights = copies.*sample_weights;
            if sample_size ~= sum(copies)
                error('wrong number of samples');
            end
            
        end
        
    end
    
end

