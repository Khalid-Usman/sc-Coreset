function XX = row_outer_products(X)
n = size(X,1);
d = size(X,2);
XX = spalloc(n,d^2,n*d^2);
fprintf('Computing outer products {xx''}: d=%d->%d ... ',d,d^2)
bstr = '';
for i = 1:n
    x = X(i,:);
    xx = x'*x;
    %XX = sprep(XX,i,[1:d^2],xx(:)');
    XX(i,:) = xx(:)';
    msg = sprintf('[%d/%d]\n',i,n);
    fprintf([bstr msg])
    bstr = repmat(sprintf('\b'),1,length(msg));
end
end
