% load new_center_example
X = rand(20,2);
mu = mean(X,1);
X = bsxfun(@minus,X,mu);
squared_norm_rows = @(X)(sum(bsxfun(@times,X,X),2));
z = squared_norm_rows(X);
X = bsxfun(@rdivide,X,sqrt(z));

figure(1), clf, hold on
draw_circle

% c = X(1,:);
c = [-0.5332 -0.846];

p = X(2,:);
plot([0 c(1)],[0 c(2)],'go-','LineWidth',1.5)
plot([0 p(1)],[0 p(2)],'ko-','LineWidth',1.5)
v = p - c;
c1 = p - (v/norm(v))*(v/norm(v)*p');      
plot([p(1) c(1)],[p(2) c(2)],'bo:','LineWidth',1)
plot([0 c1(1)],[0 c1(2)],'bo--','LineWidth',1)
plot([c(1) c1(1)],[c(2) c1(2)],'ro-','LineWidth',2)

%%
c = c1;
p = X(10,:);
plot([0 p(1)],[0 p(2)],'ko-')
v = p - c;
c1 = p - (v/norm(v))*(v/norm(v)*p');               
plot([p(1) c(1)],[p(2) c(2)],'bo:','LineWidth',1)
plot([0 c1(1)],[0 c1(2)],'bo--','LineWidth',1)
plot([c(1) c1(1)],[c(2) c1(2)],'ro-','LineWidth',2)

%%
c = c1;
p = X(3,:);
plot([0 p(1)],[0 p(2)],'ko-')
v = p - c;
c1 = p - (v/norm(v))*(v/norm(v)*p');               
plot([p(1) c(1)],[p(2) c(2)],'bo:','LineWidth',1)
plot([0 c1(1)],[0 c1(2)],'bo--','LineWidth',1)
plot([c(1) c1(1)],[c(2) c1(2)],'ro-','LineWidth',2)
