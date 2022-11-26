%%
function [j,r1,alpha,is_changed] = find_next_point(this,A,r)

if this.k == 0
    
    % Find the farthest row in X from r.
    D = A*r'; % inner product of points with center
    [~,j] = min(D);
    
    p = A(j,:);
    
    % Find the point u on the line from p1 to p2 that is closest to the
    % origin.
    v = p - r;
    v2 = v * v';
    r1 = p - v * ( v * p' ) / v2;
    alpha = norm(r-r1) / norm(r-p);  % -v*r'/sqrrt(v2);
    alpha = min(alpha,1);
    
else % k > 0
    
    % general linearized algorithm:
    
    l = r / norm(r);
    
    distFromCurrentLine = 1 - squared_norm_rows(A*(l'*l));
    [~,j] = max(distFromCurrentLine ./ this.opt_dist);
    
    q = A(j,:);
    
    % u = q^
    % v = q'
    u = q*(this.opt'*this.opt);
    v = u*(l'*l);
    
    % Find closest point on the line spanning v-q to the
    % opt_lineimal line
    d1 = v - q;
    d1 = d1 / norm(d1);
    
    r1 = d1*(u-v)'*d1 + v;
    
    % alternative implementation:
    %     d2 = this.opt;
    %     a = dot(d1,d1);
    %     b = dot(d1,d2);
    %     e = dot(d2,d2);
    %     d = a*e - b*b;
    %     c = dot(d1,q);
    %     f = dot(d2,q);
    %     s = (b*f - c*e) / d;
    %     r = q + d1*s;
    
    A = r1 / [r; q];
    A = A / sum(abs(A));
    alpha = A(2);
    
    if A(1) < 0
        alpha = -alpha;
    end
    
    r1 = (r*(1-abs(alpha)) + q*alpha);
    
end

% check if coreset changed
if abs(alpha) > 1e-10 && abs(alpha) <= 1
    is_changed = true;
else
    % % there is no row in P that is closer to the origin.
    % if P'*r > r'*r
    %     error('alpha is bad from wrong reasons');
    % else
    is_changed = false;
    r1 = r;
    j = [];
    alpha = 0;
    % end
end

end

%%
% input: A -- n-by-d matrix
% N -- number of iterations
% output x -- an n-by-1 vector x of at most N positive entries that sum to 1.
% The other entries are zeroes.
% r -- the residual ||Ax-b|| = ||r|| where b is the mean of A
% Compute the SVDCoreset Algorithm for the function f(x) = -||Ax-b||
% By translating the mean, we can approximate any vector in the convex hull of A.
% By Caratheodory Theorem, zero error (A'x = b) will be for N<=d+1
% If the mean of A is the origin then A'x|| = r.
% Otherwise, A'x-sum(x)*b = r
function [w,coreset_size] = compute(this,A,m,error_threshold)

    this.n = size(A,1);
    this.d = size(A,2);

    if this.k == 0
        this.opt = [];
        this.opt_dist = 0;
    else % k > 0
        try
            [~,~,this.opt] = svds(A,this.k);
            this.opt = this.opt'/norm(this.opt);
            this.opt_dist = 1 - squared_norm_rows((A*this.opt')*this.opt);
        catch e
            disp('svds failed!')
            warning(e.identifier,e.message)
            w = [];
            coreset_size = 0;
            return
        end
    end

    % select starting point
    j0 = 1;
    c0 = A(j0,:);

    w = sparse(this.n,1); % create sparse vector of length n
    w(j0) = 1; % update weight of coreset point j

    if this.PLOT && this.d == 2

        figure(99), clf, hold on, grid on

        %plot(P(:,1),P(:,2),'x','LineWidth',2)
        %for i = 1:this.n
        %  line([P(i,1) A(i,1)],[P(i,2) A(i,2)],'LineStyle',':')
        %end

        plot(0,0,'og')
        ang = 0:0.01:2*pi;
        x = 0;
        y = 0;
        c = 1;
        xp = c*cos(ang);
        yp = c*sin(ang);
        plot(x+xp,y+yp);
        axis equal

        plot(A(:,1),A(:,2),'ok','LineWidth',2)
        plot(A(j0,1),A(j0,2),'og','MarkerSize',20,'LineWidth',4)

    end

    coreset_size = 1;

    %fprintf('.');
    bstr = '';
    for i = 2:m

        % add m additional points to coreset
        % compute next point from input to coreset
        [j,c,alpha,is_changed] = this.find_next_point(A,c0);

        if not(is_changed)
            break;
        end

        msg = sprintf('%d/%d',i,m);

        %fprintf('.');

        % update coreset weights: x = alpha e(j)+(1-alpha)x
        w = w * (1 - abs(alpha));
        w(j) = w(j) + alpha;

        % assert that the sum of x is one and that x is positive
        if this.DEBUG
            s = full(sum(w));
            if abs(s-1) > 1e-10
                warning('sum of x: %d', s)
            end % end if s~=1
            if w ~= abs(w)
                warning('x is negative')
                difsp(w);
                w = abs(w);
            end
        end

        if this.k == 0 % k~=1
            w = w ./ sum(w);
        end

        if this.d == 2
            line([c0(1) A(j,1)],[c0(2) A(j,2)],'LineWidth',1)
            line([0 c(1)],[0 c(2)],'Color','r','LineStyle',':')
            plot(c(1),c(2),'xr','MarkerSize',10,'LineWidth',2)
        end

        % update
        coreset_size = coreset_size+1;
        c0 = c;

        %                 if mod(i,10) == 0
        %                     err = D.calc_error(P)
        %                     if err < error_threshold
        %                         break
        %                     end
        %                     %fprintf('error = %d\n',err)
        %                     msg = [msg ', error = ' num2str(err)];
        %                 end

        fprintf([bstr msg]);
        bstr = repmat(sprintf('\b'),1,length(msg));

    end
    fprintf('\n');

    this.w = w;
    this.C = this.w' * A;

end

%%
function err = calc_error(this,A)

approx = this.w' * A;

if this.k == 0
    approx_dist = full(sum(squared_norm_rows(bsxfun(@minus,A,approx))));
    best_dist = full(sum(squared_norm_rows(A)));
else
    approx = approx / norm(approx);
    approx_dist = full(sum(squared_norm_rows(bsxfun(@minus,A,A*(approx'*approx)))));
    best_dist = full(sum(squared_norm_rows(bsxfun(@minus,A,A*(this.opt'*this.opt)))));
end

err = (approx_dist-best_dist) / this.n;
err(err<0) = 0;

end




% Returns the norm of each row of X
function S = squared_norm_rows(X)
S = sum(X.^2,2);
end

% The distance function used in the SVDCoreset implementation.
% a and b must be vectors of size (n x 1)
function res = dist(a,b)
assert(size(a,1) == size(b,1));
assert(size(a,2) == 1);
assert(size(a,2) == size(b,2));
A = a * transpose(a);
B = b * transpose(b);
res = svds(B-A,1);
end

% Return the dist of column vector a from the origin
function res = dist_from_origin(a)
n = size(a,1);
A = a * transpose(a);
I = eye(n);
res = norm(I-A,2);
end






% check if coreset changed (alpha convergence)
if isempty(error_threshold)
    convergence_threshold = 1e-10;
else
    convergence_threshold = 0;
end
if abs(alpha)<convergence_threshold
    done = true;
else
    done = false;
    msg = sprintf('[%d/%d]',i,this.m);
    msg = [msg sprintf(', alpha = %1.6f',alpha)];
    fprintf([bstr msg]);
    bstr = repmat(sprintf('\b'),1,length(msg));
end

% check if coreset changed (error convergence)
if mod(i,10)==0 && not(isempty(error_threshold))
    err = SVDCoreset.calc_error(A,w);
    if err < error_threshold
        done = true;
    end
end

