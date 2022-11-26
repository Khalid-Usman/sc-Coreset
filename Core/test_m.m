%#ok<*NOPTS,*SNASGU>
global CORESET_OUTPUT STREAM_OUTPUT
CORESET_OUTPUT = 0 
STREAM_OUTPUT = 1

n = 5000;
d = 1000;
s = 1e-3;
c = sqrt(s);
k = d-1;
m = n;
e = 0;
l = 0;

[iters,sizes,errors,runtimes] = SVDCoresetTest.pivot_test(n,d,s,c,k,m,e,l);
