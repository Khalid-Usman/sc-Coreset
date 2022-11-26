
X = rand(20,2);
mu = mean(X,1);
X = bsxfun(@minus,X,mu);
squared_norm_rows = @(X)(sum(bsxfun(@times,X,X),2))
z = squared_norm_rows(X);
X = bsxfun(@rdivide,X,sqrt(z));
draw_circle
c = X(1,:);
p = X(2,:);
plot([0 c(1)],[0 c(2)],'go-')
plot([0 p(1)],[0 p(2)],'bo-')
v = p - c;
c1 = p - (v/norm(v))*(v/norm(v)*p')                
plot([c(1) c1(1)],[c(2) c1(2)],'ro:')
c = c1
p = X(3,:);
plot([0 p(1)],[0 p(2)],'ko-')
v = p - c;
c1 = p - (v/norm(v))*(v/norm(v)*p')                
plot([c(1) c1(1)],[c(2) c1(2)],'ro:')
