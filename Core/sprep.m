% Equivalent to
% >> A(I,J) = B
% But does not support "end" indexing in I and J
function A = sprep(A,I,J,B)
    I = I(:);
    J = J(:);
    [iA,jA,sA] = find(A);
    [iB,jB,sB] = find(B);
    trimA = ~(ismember(iA,I) & ismember(jA,J));
    tiA = iA(trimA);
    tjA = jA(trimA);
    tsA = sA(trimA);
    iiB = I(iB);
    jjB = J(jB);
    A = sparse(...
        [tiA(:);iiB(:)],...
        [tjA(:);jjB(:)],...
        [tsA(:);sB(:)],...
        size(A,1),size(A,2));
end