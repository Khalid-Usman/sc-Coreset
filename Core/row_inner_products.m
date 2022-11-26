function XX = row_inner_products(X)
XX = sum(bsxfun(@times,X,X),2);

