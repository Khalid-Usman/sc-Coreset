classdef Stack < handle
    % implements stack data structure. The itmes are stored in a cell
    % arrray, and can be of any type.
    properties
        stack={};
    end
    
    methods
        function this=Stack(inputStack)
            if nargin>0
                this.stack = inputStack.stack;
            else
                this.stack = {};
            end
        end
        function clear(this)
            this.stack={};
        end
       function result=isEmpty(this)
           if isempty(this.stack)
                result=true;
           else result=false;
           end
        end % function isEmpty
        function item=top(this)
          if this.isEmpty() 
              error('top from empty stack');
          end
          item= this.stack{end};
        end
        function [item]=pop(this)
          if this.isEmpty() 
              error('pop from empty stack');
          end
          item= this.stack{end};
          this.stack(end)=[];
        end % function pop       
        function push(this, item)
            this.stack{end+1}=item;
        end
    end % methods
    
end

