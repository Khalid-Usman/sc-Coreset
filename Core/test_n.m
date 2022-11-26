%#ok<*NOPTS,*SNASGU>
global CORESET_OUTPUT STREAM_OUTPUT
CORESET_OUTPUT = 0 
STREAM_OUTPUT = 1

n = 10000;
d = 100000;
s = 0.001;
c = 0.01;
k = [500 1000 2000];
m = 0;
e = 1e-3;
l = 2500;

[iters,sizes,errors,runtimes] = SVDCoresetTest.pivot_test(n,d,s,c,k,m,e,l);
