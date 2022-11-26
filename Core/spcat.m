% Equivalent to
% >> AB = [A;B]
function C = spcat(A,B)
    if size(A,1)==1 || size(B,1)==1
        % find doesn't work the same way on vectors
        C = [A;B];
    else
        [iA,jA,sA] = find(A);
        [iB,jB,sB] = find(B);
        C = sparse(...
            [iA(:);iB(:)+size(A,1)],...
            [jA(:);jB(:)],...
            [sA(:);sB(:)],...
            size(A,1)+size(B,1),max(size(A,2),size(B,2)));
    end
end