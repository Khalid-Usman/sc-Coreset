load wiki_en1_bow_1
n = size(A,1);
d = size(A,2);
X = sparse(n,d);
for i = 1
    fprintf('adding block #%d\n',i)
    load(['wiki_en1_bow_' num2str(i)])
    x1 = find(sum(A,2)>0,1,'first');
    x2 = find(sum(A,2)>0,1,'last');
    z1 = find(A(x1,:)>0,1,'first');
    X(x1,z1:end) = A(x1,z1:end);
    X(x1+1:x2,:) = A(x1+1:x2,:);
end
A = X;
% A = A(1:1000,:);
% save A1000 A