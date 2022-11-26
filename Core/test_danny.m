k = 10;
A = dlmread('D:\Research\SC_Coreset\Datasets\Basic\Pollen_Original.txt', ',', 1, 1);
A = A';
D = SVDCoreset; 
D.max_iter = 100; 
D.max_error = 0; 
D.compute(A,k); 
coreset_size = D.coreset_size, D_approx_error = D.approx_error
Ru = UniformRandomSampling(A,k,coreset_size); 
Ru_approx_error = Ru.approx_error
Rw = WeightedRandomSampling(A,k,coreset_size); 
Rw_approx_error = Rw.approx_error
