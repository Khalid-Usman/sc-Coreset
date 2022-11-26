classdef SVDBufferStream < handle
    % Class for maintaining the SVDCoreset on-line.
    % The stream collects the original points, untill we have a set a last
    % leaf of size leaf_size. Then, we compute coreset in the first level from 
    % this set. If there is already a coreset in the same level, we continue 
    % to merge coresets recursively. An item in the streaming tree 
    % (implemented using Stack) consists of two properties: coreset, which is
    % the coreset itself, and level, which is the level of the coreset in the 
    % tree. There are no two coresets in the same level. The coreset of 
    % minimum level is on the top of the stack.
    
    properties
        
        % number of points that have been fed to the stream
        % (useful for figuring out the next "time" value
        num_points_streamed = 0
        
        % Coreset that support .mergedCoreset and .computeCoreset
        % (usually of type coreset_alg)
        coreset_alg = SVDCoresetAlg
        
        % Size of original points to collect before constructing coreset.
        leaf_size
        
        % the current accumulated coreset
        C
        D
        
        k
        max_iter
        max_error
        max_size
 
        is_init = false
        
    end
    
    %%
    methods
        
        function this = SVDBufferStream()
           this.coreset_alg = SVDCoresetAlg;
        end
        
        %% Add a set of points to the stream
        function addPoints(this,P)
            
            if not(this.is_init)
                this.init();
            end
                
            n = size(P,1);
            
            for i = 1:ceil(n/this.leaf_size)
                
                xa = this.leaf_size*(i-1)+1;
                xb = this.leaf_size*i;
                
                % boundary condition fix
                xb = min(xb,n);
                
                %fprintf('Adding points: %d:%d\n',xa,xb)
                Q = P(xa:xb,:);
                fprintf('Khalid')
                this.leaf_size
                size(Q)
                this.addLeaf(Q,xa,xb);
                this.num_points_streamed = this.num_points_streamed + size(Q,1);
                
            end
            
            fprintf(['\n%s\n',repmat('- ',1,30)])
            
        end
        
        %% Unite all the coresets in the tree to a single coreset
        function D = getUnifiedCoreset(this)
            
            D = this.D;
            
        end
        
    end
    
    methods(Access = protected)
        
        % initialize SVD coreset
        function init(this)

            this.coreset_alg.k = this.k;
            if not(isempty(this.max_iter))
                this.coreset_alg.max_iter = this.max_iter;
            end
            if not(isempty(this.max_error))
                this.coreset_alg.max_error = this.max_error;
            end
            if not(isempty(this.max_size))
                this.coreset_alg.max_size = this.max_size;
            end
        end
        
        %% Construct a corest and add it to the tree
        function addLeaf(this,Q,xa,xb)
            
            start_time = tic;
            
            if isempty(this.D)
                
                fprintf('\n%s\n',repmat('- ',1,30))
            
                % compute initial coreset
                fprintf('Computing coreset: [%d:%d]\n',xa,xb)
                
                D = this.coreset_alg.computeCoreset(Q);
                D.xa = xa;
                D.xb = xb;
            
                fprintf('num iter = %d/%d, coreset size = %d, error = %f\n', ...
                    D.num_iter,D.max_iter,D.coreset_size,D.coreset_errors(end));
                
                % set first coreset
                this.D = D;
                
            else
                
                %fprintf('\n%s\n',repmat('. ',1,30))
                
                xxa = min(xa,this.D.xa);
                xxb = max(xb,this.D.xb);
                
                % merge coreset with incoming points
                %fprintf('Merging coresets: [%d:%d--%d:%d] -> [%d:%d]\n',this.D.xa,this.D.xb,xa,xb,xxa,xxb)    
                
                %D = this.coreset_alg.computeCoreset(Q);
                %this.D = this.coreset_alg.mergedCoreset(this.D.C,D.C);
                
                prev_num_iter = this.D.num_iter;
                
                D = this.coreset_alg.mergedCoreset(this.D.C,Q);
                D.xa = xxa;
                D.xb = xxb;
   
                this_num_iter = D.num_iter;
                D.num_iter = prev_num_iter + this_num_iter;
                
                %fprintf('total iter = %d (+%d/%d), coreset size = %d, error = %f\n', ...
                % D.num_iter,this_num_iter,D.max_iter,D.coreset_size,D.coreset_errors(end));
                
                % save accumulated coreset
                this.D = D;
            
            end
            
            %fprintf('Running time = %1.2f min\n',toc(start_time)/60)
            
        end
        
    end
    
end

function fprintf(varargin)
    global STREAM_OUTPUT
    if isempty(STREAM_OUTPUT)
        STREAM_OUTPUT = 1;
    end
    if STREAM_OUTPUT
       builtin('fprintf',varargin{:}); 
    end
end
    
