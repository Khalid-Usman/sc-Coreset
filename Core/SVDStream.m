classdef SVDStream < handle
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
 
        % don't delete data while constructing tree
        save_tree = false % boolean
        
        % saved list of coresets constructured over time.
        % Only if save_tree==true
        coreset_list
        
        % number of points that have been fed to the stream
        % (useful for figuring out the next "time" value
        num_points_streamed = 0
        
        % Coreset that support .mergedCoreset and .computeCoreset
        % (usually of type coreset_alg)
        coreset_alg = SVDCoresetAlg
        
        % Size of original points to collect before constructing coreset.
        leaf_size
        
        % stack to implement the coresets tree
        stack
        
        % The set of last original points collected, since the last coreset
        % construction. Size at most leaf_size.
        % Type = AbstractFunctionSet
        last_leaf
        
        k
        max_iter
        max_error
        max_size 
        
        is_init = false
        
    end
    
    %%
    methods
        
        function this = SVDCoresetStream()
           this.stack = Stack;
           this.coreset_alg = SVDCoresetAlg;
        end
        
        %% Add a set of points to the stream
        function addPoints(this,P)
            
            if not(this.is_init)
                this.init();
            end
                
            n = size(P,1);
            
            for i = 1:ceil(n/this.leaf_size)
                
                fprintf('\n%s\n',repmat('- ',1,30))
                
                xa = this.leaf_size*(i-1)+1;
                xb = this.leaf_size*i;
                
                % boundary condition fix
                xb = min(xb,n);
                
                fprintf('Adding points: [%d:%d]\n',xa,xb)
                Q = P(xa:xb,:);
                this.addLeaf(Q,xa,xb);
                this.num_points_streamed = this.num_points_streamed + size(Q,1);
                
            end
            
            fprintf(['\n%s\n',repmat('- ',1,30)])
            
        end
        
        %% Unite all the coresets in the tree to a single coreset
        function D = getUnifiedCoreset(this)
            
            % copy stack so we can retrieve again
            s = Stack(this.stack);
            
            isFirst = true;
            while not(s.isEmpty())
                C = s.pop().coreset;
                if isFirst
                    isFirst = false;
                    D = C;
                else
                    D = this.coreset_alg.mergedCoreset(D,C);
                end
            end
            if isempty(this.last_leaf)
                nLeaf = 0;
            else
                nLeaf = this.last_leaf.M.n;
            end
            if nLeaf > 0
                D = this.coreset_alg.mergedCoreset(D,this.last_leaf);
            end
            
        end
        
    end
    
    methods(Access = protected)
        
        % initialize SVD coreset
        function init(this)
            
            this.last_leaf = [];
            this.coreset_list = [];
            this.stack.clear();
            
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
            
            % the line that computes coreset
            coreset = this.coreset_alg.computeCoreset(Q);
            coreset.xa = xa;
            coreset.xb = xb;
            
            % --- ADD CORESET ---
            this.addCoreset(coreset);
            
        end
                
        %%
        % Get a coreset and insert it to the streaming tree.
        % Merge with existing coresets untill there are no two coresets in
        % the same level of the tree.
        % 'coreset' is the output of coreset_alg, that i
        % accepted by coreset_alg.mergedCoreset
        function addCoreset(this,coreset)
            
            % addCoreset2List
            if this.save_tree
                i = length(this.coreset_list)+1;
                this.coreset_list{i} = coreset;
            end
            
            level = 1;
            
            while not(this.isCorrectLevel(level))
                
                stackItem = this.stack.pop();
                
                % TODO: there is a horrible bug in the original Stream.m
                % coresets are randomly merged with themselves for no
                % apparent reason
                if coreset.xa == stackItem.coreset.xa || coreset.xb == stackItem.coreset.xb
                    % this coreset already exists, so don't add it
                    fprintf('Warning: duplicate coreset detected, dropping computation')
                    break
                end
                
                % --- MERGE CORESETS ---
                coreset = this.coreset_alg.mergedCoreset(stackItem.coreset,coreset);
                
                % addCoreset2List
                if this.save_tree
                    i = length(this.coreset_list)+1;
                    this.coreset_list{i} = coreset;
                end
                
                level = level+1;
                
            end
            
            newStackItem.level = level;
            newStackItem.coreset = coreset;
            this.stack.push(newStackItem);
            
        end
        
        %%
        % Returns false if integer 'level' is the same as the level of the
        % last coreset that was inserted to the stack. Used in the recursion 
        % to decide whether to merge the current coreset with a previous one, 
        % and go up another level.
        function result = isCorrectLevel(this,level)
            
            if this.stack.isEmpty()
                result = true;
            elseif this.stack.top().level>level
                result = true;
            elseif  this.stack.top().level==level
                result = false;
            else % this.stack.top().level > level
                error('Stack error: this.stack.top().level > level');
            end
            
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
    
