% cluster1 = parcluster('mikhail_c3x64_1');
% cluster2 = parcluster('mikhail_c3x64_2');
cluster = cluster1;

%%
for i = 1:64
    
    while not(strcmp(cluster.Jobs(i).State,'finished'))
        pause(10);
    end
        
    fprintf('Loading job %d ... ',i)
    
    S = load(cluster.Jobs(i),'block_no','stream','D','r1','r2','running_time','num_iter','coreset_size','coreset_errors');
    
    save(sprintf('data/k%d/block%d',S.stream.k,S.block_no),'S')
    
    fprintf('Done!\n')
  
    
end
