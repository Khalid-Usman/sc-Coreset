classdef SVDCoresetAlg < handle
    
    properties
        
        k
        max_iter
        max_error
        max_size
        
    end
    
    %%
    methods
        
        %% Compute coreset from signal points
        % Input:  SignalPointSet P
        % Output: SVDCoreset D
        function D = computeCoreset(this,A)
            
            D = SVDCoreset;
            D.max_iter = this.max_iter;
            D.max_error = this.max_error;
            D.max_size = this.max_size;

            D.compute(A,this.k);
            
        end
        
        %% Merge two SVDCoresets C1,C2 into a new one D.C.
        %   Input:    SVDCoreset D1,D2
        %   Output:   SVDCoreset D
        function D = mergedCoreset(this,C1,C2)
            
            A = spcat(C1,C2);
            
            D = SVDCoreset;
            D.max_iter = this.max_iter;
            D.max_error = this.max_error;
            D.max_size = this.max_size;

            D.compute(A,this.k);
            
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
