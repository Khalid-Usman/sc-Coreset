    %%
    methods (Static)
        
        % calculate error
        function err = compute_error(X,W)
            
            n = size(X,1);
            d = size(X,2);
            
            C = W*X;
            approx_dist = full(sum(row_inner_products(X'*X-C'*C)));
            best_dist = full(sum(row_inner_products(X'*X)));
            err = abs(approx_dist-best_dist)/n;
            err(err<0) = 0;
            
%             n = size(X,1);
%             d = size(X,2);
%             xx = sparse(n,d^2);
%             wxx = sparse(n,d^2);
%             for i = 1:size(X,1)
%                 xi = X(1,:)';
%                 xx(i,:) = reshape((xi*xi'),[1 d^2]);
%                 wxx(i,:) = reshape(W(i,i)*(xi*xi'),[1 d^2]);
%             end
            
%             n = size(XX,1);
%             WXX = W*XX;
%             err = norm(sum(XX,2)-sum(WXX,2),'fro');
%             
%             %approx_dist = norm(sum(XX,2)-sum(WXX,2),'fro');
%             %best_dist = norm(sum(XX,2),'fro');
%             %err = abs(approx_dist-best_dist);
%             err(err<0) = 0;

        end
        
    end