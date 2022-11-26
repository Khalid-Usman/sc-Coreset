function coreset = computeCoreset(A, k, coresetSize, maxIter, fName)
    D = SVDCoreset;
    D.max_iter = maxIter;
    D.max_error = 0;
    D.max_size = coresetSize;
    size(A)   
    D.compute(A, k);
    %coreset_size = D.coreset_size, D_approx_error = D.approx_error;
    coreset = full(D.C);
    save(fName, 'coreset')
end
